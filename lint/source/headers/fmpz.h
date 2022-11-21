#pragma once

#include<vector>
#include "types.h"
#include "asm_inlines.h"
#include "random_state.h"
#include "nmod.h"
#include "my_mpn.h"

struct fmpn_struct {
    slimb length;
    slimb alloc;
    ulimb limbs[];
};

inline void _fmpn_free(fmpn_struct* x)
{
    flint_free(x);
}

inline fmpn_struct* _fmpn_alloc(slimb fit)
{
    fmpn_struct* X = (fmpn_struct*) flint_malloc(sizeof(fmpn_struct) + fit*sizeof(ulimb));
    X->alloc = fit;
    return X;
}

// assuming x is already big, ensure at least "fit" limbs allocated
inline fmpn_struct* _fmpn_fit_alloc(fmpn_struct* X, slimb fit)
{
    if (X->alloc < fit)
    {
        X->alloc = std::max(2*X->alloc, fit);
        X = (fmpn_struct*) flint_realloc(X, sizeof(fmpn_struct) + X->alloc*sizeof(ulimb));
    }
    return X;
}

inline fmpn_struct* _fmpn_fit_alloc_destroy(fmpn_struct* X, slimb fit)
{
    if (X->alloc < fit)
    {
        ulimb n = std::max(2*X->alloc, fit);
        flint_free(X);
        X = (fmpn_struct*) flint_malloc(sizeof(fmpn_struct) + n*sizeof(ulimb));
        X->alloc = n;
    }
    return X;
}




#define COEFF_MAX ((slimb(1) << (FLINT_BITS-2)) - 1)
#define COEFF_MIN (-COEFF_MAX)

/*
        the data
        11xxxxxx    small negative immediate
        11000000    large negative fmpn_struct* (this one is nullptr and not used)
-2^63:  10xxxxxx    large negative fmpn_struct*
2^63-1: 01xxxxxx    large positive fmpn_struct*
        00111111    small positive immediate
        00xxxxxx    small positive immediate

large negative <= -2^62, -(2^62-1) <= small <= 2^62-1, 2^62 <= large positive
*/


inline slimb _fmpz_from_ptr(fmpn_struct* x, ulimb s) {
    FLINT_ASSERT(s == 0 || s == 1);
    return (slimb)rotate_right(((ulimb)(x)) | (s + 1), 2);
}

inline fmpn_struct* _fmpz_to_ptr(slimb a) {
    return (fmpn_struct*)(ulimb(a) << 2);
}

// return a non-ptr when a is not a ptr
inline fmpn_struct* _fmpz_to_ptr_no_check(slimb a) {
    ulimb A = a;
    A = rotate_left(A, 2);
//11xxxxxx  xxxxxx11  xxxxxx10
//10xxxxxx  xxxxxx10  xxxxxx01
//01xxxxxx  xxxxxx01  xxxxxx00
//00xxxxxx  xxxxxx00  ??????11
    return (fmpn_struct*)((A - 1)&(~ulimb(1)));
}


// nullptr is not used so 11000000 is an illegal fmpz
inline bool _fmpz_is_small(slimb x) {
    ulimb X = x;
    return !has_top_bit(X^(X+X));
}

inline bool _fmpz_is_ptr(slimb x) {
    ulimb X = x;
    return has_top_bit(X^(X+X));
}


// |x| < 2^(FLINT_BITS-2)
inline bool _si_is_small(slimb x)
{
    ulimb l = ulimb(1) << (FLINT_BITS - 2);
    ulimb X = (ulimb)x + l;
    return 0 < X && X < 2*l;
}

inline bool _ui_is_small(ulimb x)
{
    return x < (ulimb(1) << (FLINT_BITS - 2));
}

inline slimb _fmpz_neg(slimb x) {
    return _fmpz_is_small(x) ? -x : x ^ (ulimb(3) << (FLINT_BITS - 2));
}

inline slimb _fmpz_abs(slimb x) {
    return x < 0 ? _fmpz_neg(x) : x;
}

inline slimb _fmpz_make_si(slimb a) {
    if (_si_is_small(a))
        return a;
    
    fmpn_struct* x = _fmpn_alloc(1);
    x->length = 1;
    x->limbs[0] = my_abs(a);
    return _fmpz_from_ptr(x, a < 0);
}

inline slimb _fmpz_make_ui(ulimb a) {
    if (_si_is_small(a))
        return a;
    
    fmpn_struct* x = _fmpn_alloc(1);
    x->length = 1;
    x->limbs[0] = a;
    return _fmpz_from_ptr(x, false);
}


slimb _fmpz_copy(slimb a);
bool _fmpz_is_canonical(slimb);

struct fmpzc;

struct fmpz {
    typedef fmpzc source_t;
    typedef fmpzc modulus_t;

    slimb data;

    bool is_small() const {return _fmpz_is_small(data);}
    slimb small() const {FLINT_ASSERT(_fmpz_is_small(data)); return data;}
    bool is_ptr() const {return _fmpz_is_ptr(data);}
    fmpn_struct * ptr() {return _fmpz_to_ptr(data);}
    fmpn_struct * _ptr() {return _fmpz_to_ptr_no_check(data);}
    const ulimb* limbs() const {return _fmpz_to_ptr(data)->limbs;}
    ulimb* limbs() {return _fmpz_to_ptr(data)->limbs;}
    slimb length() const {return _fmpz_to_ptr(data)->length;}
    slimb alloc() {return _fmpz_to_ptr(data)->alloc;}
    slimb size() {return is_ptr() ? length() : !is_zero();}

    int sign() {return my_cmp(data, slimb(0));}
    int cmp_small(slimb a) const {return my_cmp(data, a);}
    bool is_zero() const {return data == 0;}
    bool is_one() const {return data == 1;}
    bool is_positive() const {return data > 0;}
    bool is_negative() const {return data < 0;}
    bool is_even() const {return is_small() ? ::is_even(small()) : ::is_even(limbs()[0]);}
    bool is_odd() const {return is_small() ? ::is_odd(small()) : ::is_odd(limbs()[0]);}

    slimb get_data() const {
        return data;
    }

    slimb steal_data() {
        slimb x = data;
        data = 0;
        return x;
    }

    ~fmpz() {
        if (_fmpz_is_ptr(data))
            _fmpn_free(_fmpz_to_ptr(data));
    }

    fmpz() : data(0) {}
    fmpz(int a) {data = a;}
    fmpz(slimb a) : data(_fmpz_make_si(a)) {}
    fmpz(ulimb a) : data(_fmpz_make_ui(a)) {}
    fmpz(const fmpz& other) : data(_fmpz_copy(other.data)) {}
    fmpz(fmpz&& other) noexcept {
        data = other.data;
        other.data = 0;
    }

	fmpz& operator=(const fmpz& other) {
        if (_fmpz_is_ptr(data))
            _fmpn_free(_fmpz_to_ptr(data));
        data = 0;
        data = _fmpz_copy(other.data);
        return *this;
    }

    fmpz& operator=(fmpz&& other) noexcept {
        std::swap(data, other.data);
        return *this;
    }

    void set_strn(const char* c, size_t n);
    void set_c_str(const char* c) {set_strn(c, strlen(c));}
    void set_limbs(const ulimb* p, ulimb n, int sign);

    fmpz(const ulimb* p, ulimb n, int sign) : data(0) {
        set_limbs(p, n, sign);
    }

    fmpz(const char* c) : data(0) {
        set_strn(c, strlen(c));
    }

    bool operator<(const fmpz& other) const noexcept;
    bool operator<=(const fmpz& other) const noexcept;
};

/*
    bool operator<(const fmpz& other) const noexcept
    {
        return fmpz_cmp(get(), other.get()) < 0;
    }

    bool operator<=(const fmpz& other) const noexcept
    {
        return fmpz_cmp(get(), other.get()) <= 0;
    }
*/


struct fmpzc {
    slimb data;

    fmpzc() : data(0) {}
    fmpzc(slimb d) : data(d) {}
    fmpzc(const fmpz& d) : data(d.get_data()) {}

    bool is_small() const {return _fmpz_is_small(data);}
    slimb small() const {FLINT_ASSERT(_fmpz_is_small(data)); return data;}
    bool is_ptr() const {return _fmpz_is_ptr(data);}
    fmpn_struct * ptr() {return _fmpz_to_ptr(data);}
    fmpn_struct * _ptr() {return _fmpz_to_ptr_no_check(data);}
    const ulimb* limbs() const {return _fmpz_to_ptr(data)->limbs;}
    ulimb* limbs() {return _fmpz_to_ptr(data)->limbs;}
    slimb length() const {return _fmpz_to_ptr(data)->length;}
    slimb alloc() {return _fmpz_to_ptr(data)->alloc;}
    slimb size() {return is_ptr() ? length() : !is_zero();}

    fmpzc abs() {return fmpzc(_fmpz_abs(data));}
    int sign() {return my_cmp(data, slimb(0));}
    int cmp_small(slimb a) const {return my_cmp(data, a);}
    bool is_zero() const {return data == 0;}
    bool is_one() const {return data == 1;}
    bool is_positive() const {return data > 0;}
    bool is_negative() const {return data < 0;}
    bool is_even() const {return is_small() ? ::is_even(small()) : ::is_even(limbs()[0]);}
    bool is_odd() const {return is_small() ? ::is_odd(small()) : ::is_odd(limbs()[0]);}

    void get_limbs(ulimb* p, ulimb n) const;

    // cmp(this, a)
    int cmp_small(slimb a) {
        FLINT_ASSERT(_fmpz_is_small(a));
        return my_cmp(data, a);
    }

    fmpz copy() const {fmpz z; z.data = _fmpz_copy(data); return z;}

/*
    slimb copy_data() const {return _fmpz_copy(data);}

    void reset(slimb d) {
        if (is_ptr())
            flint_free(ptr());
        data = d;
    }

    void clear() {reset(0);}
*/
};

bool fmpz_get_ui_vec(ulimb* x, ulimb xn, fmpzc a);

int fmpz_cmp(fmpzc a, fmpzc b);

inline bool fmpz_is_canonical(fmpzc a) {return _fmpz_is_canonical(a.data);}
inline bool fmpz_is_canonical(fmpz& a) {return _fmpz_is_canonical(a.data);}

std::ostream& fmpz_write(std::ostream& o, fmpzc a);
inline std::ostream& operator<<(std::ostream& o, fmpzc a) {return fmpz_write(o, a);}
inline std::ostream& operator<<(std::ostream& o, const fmpz& a) {return fmpz_write(o, a);}


struct fmpz_ring {
    typedef fmpz elem_t;
};

fmpn_struct* _fmpz_fit_alloc(fmpz& x, slimb fit, int s);
fmpn_struct* _fmpz_fit_destroy(fmpz& x, slimb fit, int s);
slimb _fmpz_normalize_value(fmpn_struct* X, slimb len, int s);

struct fmpz_ulimb_realizer : ulimb_realizer {
    fmpz& _x;
    ulimb _fit;
    int _s;
    fmpz_ulimb_realizer(fmpz& x, ulimb fit, int s) : _x(x), _fit(fit), _s(s) {}
    ulimb* get_limbs() {return _fmpz_fit_alloc(_x, _fit, _s)->limbs;}
};

inline void fmpz_set_small(fmpz& x, slimb a)
{
    FLINT_ASSERT(_si_is_small(a));
    if (UNLIKELY(x.is_ptr()))
        _fmpn_free(x.ptr());
    x.data = a;
}

inline void fmpz_set_small(fmpz& x, ulimb a) {fmpz_set_small(x, (slimb) a);}
inline void fmpz_zero(fmpz& x) {fmpz_set_small(x, slimb(0));}
inline void fmpz_one(fmpz& x) {fmpz_set_small(x, slimb(1));}

inline bool fmpz_is_zero(fmpzc a) {return a.is_zero();}
inline bool fmpz_is_one(fmpzc a) {return a.is_one();}
inline int fmpz_sgn(fmpzc a) {return a.sign();}
inline int fmpz_sign(fmpzc a) {return a.sign();}
inline void fmpz_swap(fmpz& a, fmpz& b) {std::swap(a.data, b.data);}

inline bool is_zero(fmpz_ring& ZZ, const fmpz& a) {return a.is_zero();}
inline bool is_one(fmpz_ring& ZZ, const fmpz& a) {return a.is_one();}
inline void zero(fmpz_ring& ZZ, fmpz& x) {fmpz_zero(x);}
inline void one(fmpz_ring& ZZ, fmpz& x) {fmpz_one(x);}

bool fmpz_equal(fmpzc f, fmpzc g);
bool fmpz_equal_ui(fmpzc a, ulimb b);
bool fmpz_equal_si(fmpzc a, slimb g);
int fmpz_cmp(fmpzc a, fmpzc b);
int fmpz_cmpabs(fmpzc a, fmpzc b);
int fmpz_cmp2abs(fmpzc a, fmpzc b);
int fmpz_cmp_si(fmpzc a, slimb b);
int fmpz_cmp_ui(fmpzc a, ulimb b);

inline bool equal(fmpz_ring& ZZ, const fmpz& a, const fmpz& b) {return fmpz_equal(a, b);}



bool fmpz_is_probable_prime(fmpzc a);
bool fmpz_is_prime(fmpzc a);
void fmpz_next_prime(fmpz& x, fmpzc a, bool proved);
void fmpz_prev_prime(fmpz& x, fmpzc a, bool proved);

void fmpz_set(fmpz& x, fmpzc a);
void fmpz_abs(fmpz& x);
void fmpz_abs(fmpz& x, fmpzc a);
void fmpz_neg(fmpz& x);
void fmpz_neg(fmpz& x, fmpzc a);
void fmpz_add(fmpz& x, fmpzc a, fmpzc b);
void fmpz_sub(fmpz& x, fmpzc a, fmpzc b);
void fmpz_mul(fmpz& x, fmpzc a, fmpzc b);
void fmpz_sqr(fmpz& x, fmpzc a);
void fmpz_add_ui(fmpz& x, fmpzc a, ulimb b);
void fmpz_sub_ui(fmpz& x, fmpzc a, ulimb b);
void fmpz_mul_ui(fmpz& x, fmpzc a, ulimb b);
void fmpz_add_si(fmpz& x, fmpzc a, slimb b);
void fmpz_sub_si(fmpz& x, fmpzc a, slimb b);
void fmpz_mul_si(fmpz& x, fmpzc a, slimb b);

void fmpz_submul_ui(fmpz& x, fmpzc a, ulimb b);



bool fmpz_divides(fmpz& q, fmpzc a, fmpzc b);
bool fmpz_divisible(fmpzc a, fmpzc b);
int fmpz_divisible_ui(fmpzc a, ulimb b);
int fmpz_divisible_si(fmpzc a, slimb b);
void fmpz_divexact(fmpz& q, fmpzc a, fmpzc b);
void fmpz_divexact_ui(fmpz& q, fmpzc a, ulimb b);
void fmpz_divexact_si(fmpz& q, fmpzc a, slimb b);

void fmpz_tdiv_qr(fmpz& q, fmpz& r, fmpzc a, fmpzc b);
void fmpz_fdiv_qr(fmpz& q, fmpz& r, fmpzc a, fmpzc b);
void fmpz_cdiv_qr(fmpz& q, fmpz& r, fmpzc a, fmpzc b);

ulimb fmpz_fdiv_ui(fmpzc a, ulimb b);
void fmpz_tdiv_q_ui(fmpz& q, fmpzc a, ulimb b);
void fmpz_fdiv_q_ui(fmpz& q, fmpzc a, ulimb b);


void fmpz_mod(fmpz& r, fmpzc a, fmpzc b);
void fmpz_sqrmod(fmpz& x, fmpzc a, fmpzc m);
void fmpz_mulmod(fmpz& x, fmpzc a, fmpzc b, fmpzc m);
void fmpz_powmod_ui(fmpz& x, fmpzc a, ulimb e, fmpzc m);
void fmpz_powmod_fmpz(fmpz& x, fmpzc a, fmpzc e, fmpzc m);

int fmpz_kronecker(fmpzc a, fmpzc b);
void fmpz_pow_ui(fmpz& a, fmpzc b, ulimb e);

void fmpz_fib_ui(fmpz& fn, ulimb n);

void fmpz_random_of_bits_unsigned(random_state& state, fmpz& x, ulimb bits);

inline void random_of_bits(random_state& state, fmpz_ring& ZZ, fmpz& x, ulimb bits) {
    fmpz_random_of_bits_unsigned(state, x, bits);
    if (state.get_bit())
        fmpz_neg(x); 
}

void fmpz_randtest(fmpz& f, random_state& state, flint_bitcnt_t bits);
void fmpz_randtest_not_zero(fmpz& f, random_state& state, flint_bitcnt_t bits);
void fmpz_randtest_mod(fmpz& f, random_state& state, fmpzc m);

bool fmpz_test_bit(fmpzc a, ulimb i);
ulimb fmpz_trailing_zeros(fmpzc a);
ulimb fmpz_bits(fmpzc a);
void fmpz_tdiv_q_2exp(fmpz& q, fmpzc a, ulimb e);
void fmpz_mul_2exp(fmpz& x, fmpzc a, ulimb e);

inline void mul_2exp(fmpz_ring& ZZ, fmpz& x, const fmpz& a, ulimb e) {fmpz_mul_2exp(x, a, e);}
inline ulimb nbits(fmpzc a) {return fmpz_bits(a);}
inline ulimb nbits(fmpz_ring& R, fmpzc a) {return fmpz_bits(a);}

double fmpz_get_d(fmpzc a);
double fmpz_get_d_2exp(slimb* exp, fmpzc a);

inline bool fmpz_abs_fits_ui(fmpzc a) {
    return a.is_small() || a.length() <= 1;
}

// abs(a) mod 2^64
inline ulimb fmpz_abs_get_ui(fmpzc a) {
    if (a.is_small())
        return my_abs(a.small());
    else
        return a.limbs()[0];
}

// return a mod 2^64 assuming a >= 0
inline ulimb fmpz_get_ui(fmpzc a) {
    FLINT_ASSERT(!a.is_negative());
    if (a.is_small())
        return a.small();
    else
        return a.limbs()[0];
}

bool fmpz_fits_si(fmpzc x);
slimb fmpz_get_si(fmpzc a);

void fmpz_set_ui(fmpz& x, ulimb a);
void fmpz_set_uiui(fmpz& x, ulimb a1, ulimb a0);
void fmpz_set_uiuiui(fmpz& x, ulimb a2, ulimb a1, ulimb a0);

void fmpz_neg_ui(fmpz& x, ulimb a);
void fmpz_neg_uiui(fmpz& x, ulimb a1, ulimb a0);

void fmpz_set_si(fmpz& x, slimb a);
void fmpz_set_signed_uiui(fmpz& x, ulimb a1, ulimb a0);
void fmpz_set_signed_uiuiui(fmpz& x, ulimb a2, ulimb a1, ulimb a0);

void fmpz_gcd(fmpz& g, fmpzc a, fmpzc b);
void fmpz_gcd(fmpz& g, fmpzc a, fmpzc b, fmpzc c);
void fmpz_gcdx(fmpz& g, fmpz& s, fmpz& t, fmpzc a, fmpzc b);
void fmpz_gcdinv(fmpz& g, fmpz& s, fmpzc a, fmpzc b);
void fmpz_invmod(fmpz& x, fmpzc a, fmpzc b);

inline void randtest(fmpz_ring& ZZ, fmpz& x, random_state& state, flint_bitcnt_t bits) {fmpz_randtest(x, state, bits);}
inline void randtest_not_zero(fmpz_ring& ZZ, fmpz& x, random_state& state, flint_bitcnt_t bits) {fmpz_randtest_not_zero(x, state, bits);}
inline void swap(fmpz_ring& ZZ, fmpz& x, fmpz& y) {std::swap(x.data, y.data);}
inline void set(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_set(x, a);}
inline void neg(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_neg(x, a);}
inline void neg(fmpz_ring& ZZ, fmpz& x) {fmpz_neg(x);}
inline void add(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_add(x, a, b);}
inline void add(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_add(x, x, a);}
inline void sub(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_sub(x, a, b);}
inline void sub(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_sub(x, x, a);}
inline void mul(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_mul(x, a, b);}
inline void sqr(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_sqr(x, a);}
inline bool divisible(fmpz_ring& ZZ, fmpzc a, fmpzc b) {return fmpz_divisible(a, b);}
inline void divexact(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_divexact(x, a, b);}
inline void mod(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_mod(x, a, b);}
inline void tdiv_qr(fmpz_ring& ZZ, fmpz& q, fmpz& r, fmpzc a, fmpzc b) {fmpz_tdiv_qr(q, r, a, b);}
inline void fdiv_qr(fmpz_ring& ZZ, fmpz& q, fmpz& r, fmpzc a, fmpzc b) {fmpz_fdiv_qr(q, r, a, b);}
inline void cdiv_qr(fmpz_ring& ZZ, fmpz& q, fmpz& r, fmpzc a, fmpzc b) {fmpz_cdiv_qr(q, r, a, b);}
inline void gcd(fmpz_ring& ZZ, fmpz& g, fmpzc a, fmpzc b) {fmpz_gcd(g, a, b);}
inline void gcd(fmpz_ring& ZZ, fmpz& g, fmpzc a, fmpzc b, fmpzc c) {fmpz_gcd(g, a, b, c);}
inline void gcdx(fmpz_ring& ZZ, fmpz& g, fmpz& s, fmpz& t, fmpzc a, fmpzc b) {fmpz_gcdx(g, s, t, a, b);}
inline void mulmod(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b, fmpzc m) {fmpz_mulmod(x, a, b, m);}
inline void sqrmod(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc m) {fmpz_sqrmod(x, a, m);}
inline void powmod(fmpz_ring& ZZ, fmpz& x, fmpzc a, ulimb e, fmpzc m) {fmpz_powmod_ui(x, a, e, m);}
inline void powmod(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc e, fmpzc m) {fmpz_powmod_fmpz(x, a, e, m);}

inline fmpz operator -(fmpzc a) {fmpz x; fmpz_neg(x, a); return x;}
inline fmpz operator -(fmpz&& x) {fmpz_neg(x, x); return x;}
inline fmpz operator +(fmpzc a, fmpzc b) {fmpz x; fmpz_add(x, a, b); return x;}
inline fmpz operator -(fmpzc a, fmpzc b) {fmpz x; fmpz_sub(x, a, b); return x;}
inline fmpz operator *(fmpzc a, fmpzc b) {fmpz x; fmpz_mul(x, a, b); return x;}

