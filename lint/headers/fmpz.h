#pragma once

#include<vector>
#include "types.h"
#include "asm_inlines.h"
#include "rand_state.h"
#include "nmod.h"
#include "my_mpn.h"
#include "generics.h"


struct fmpn_struct {
    slimb length;
    slimb alloc;
    ulimb limbs[];
};

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


inline void _fmpn_free(fmpn_struct* x)
{
    flint_free(x);
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
    assert(s == 0 || s == 1);
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

bool _fmpz_is_canonical(slimb);
slimb _fmpz_copy(slimb a);

struct fmpzc;

struct fmpz {
    typedef fmpzc source_t;

    slimb data;

    int sign() {return my_cmp(data, slimb(0));}
    bool is_zero() const {return data == 0;}
    bool is_one() const {return data == 1;}
    bool is_positive() const {return data > 0;}
    bool is_negative() const {return data < 0;}
    bool is_small() const {return _fmpz_is_small(data);}
    slimb small() {return data;}
    bool is_ptr() const {return _fmpz_is_ptr(data);}
    fmpn_struct * ptr() {return _fmpz_to_ptr(data);}
    fmpn_struct * _ptr() {return _fmpz_to_ptr_no_check(data);}
    ulimb* limbs() {return _fmpz_to_ptr(data)->limbs;}
    slimb length() {return _fmpz_to_ptr(data)->length;}
    slimb alloc() {return _fmpz_to_ptr(data)->alloc;}
    slimb size() {return is_ptr() ? length() : !is_zero();}

    slimb get_data() const
    {
        return data;
    }

    slimb steal_data()
    {
        slimb x = data;
        data = 0;
        return x;
    }

    ~fmpz()
    {
        if (_fmpz_is_ptr(data))
            flint_free(_fmpz_to_ptr(data));
    }

    fmpz() : data(0) {}
    fmpz(int a) {data = a;}
    fmpz(slimb a) {data = a;} // TODO
    fmpz(ulimb a) {data = a;} // TODO

    fmpz(const fmpz& other)
    {
        data = _fmpz_copy(other.data);
    }

    fmpz(fmpz&& other)
    {
        data = other.data;
        other.data = 0;
    }

	fmpz& operator=(const fmpz& other)
    {
        if (_fmpz_is_ptr(data))
            flint_free(_fmpz_to_ptr(data));
        data = 0;
        data = _fmpz_copy(other.data);
        return *this;
    }

    fmpz& operator=(fmpz&& other) noexcept
    {
        data = other.data;
        other.data = 0;
        return *this;
    }

    void set_strn(const char* c, size_t n);
    void set_limbs(ulimb* p, ulimb n, int sign);

    fmpz(ulimb* p, ulimb n, int sign) : data(0) {
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
    fmpzc(slimb d) : data(d) {}
    fmpzc(const fmpz& d) : data(d.get_data()) {}

    int sign() const {return my_cmp(data, slimb(0));}
    bool is_zero() const {return data == 0;}
    bool is_one() const {return data == 1;}
    bool is_positive() const {return data > 0;}
    bool is_negative() const {return data < 0;}
    bool is_small() const {return _fmpz_is_small(data);}
    slimb small() const {return data;}
    bool is_ptr() const {return _fmpz_is_ptr(data);}
    fmpn_struct * ptr() const {return _fmpz_to_ptr(data);}
    fmpn_struct * _ptr() const {return _fmpz_to_ptr_no_check(data);}
    const ulimb* limbs() const {return _fmpz_to_ptr(data)->limbs;}
    slimb length() const {return _fmpz_to_ptr(data)->length;}
    slimb alloc() const {return _fmpz_to_ptr(data)->alloc;}
    slimb size() const {return is_ptr() ? length() : !is_zero();}
    void get_limbs(ulimb* p, ulimb n) const;

    // cmp(this, a)
    int cmp_small(slimb a) {
        FLINT_ASSERT(_fmpz_is_small(a));
        return my_cmp(data, a);
    }

    void clear() {
        if (is_ptr())
            flint_free(ptr());
        data = 0;
    }

    slimb copy_abs() const {
        if (is_small())
            return my_abs(small());
        fmpn_struct* A = ptr();
        auto Alen = A->length;
        fmpn_struct* X = _fmpn_alloc(Alen);
        MPN_COPY(X->limbs, A->limbs, Alen);
        X->length = Alen;
        return _fmpz_from_ptr(X, 0);
    }
};

int fmpz_cmp(fmpzc a, fmpzc b);

inline bool fmpz_is_canonical(fmpzc a) {return _fmpz_is_canonical(a.data);}
inline bool fmpz_is_canonical(fmpz& a) {return _fmpz_is_canonical(a.data);}

std::ostream& fmpz_write(std::ostream& o, fmpzc a);
inline std::ostream& operator<<(std::ostream& o, fmpzc a) {return fmpz_write(o, a);}
inline std::ostream& operator<<(std::ostream& o, const fmpz& a) {return fmpz_write(o, a);}


struct fmpz_ring {
    typedef fmpz elem_t;
};

fmpn_struct* _fmpz_promote(fmpz& x, slimb fit, int s);
fmpn_struct* _fmpz_fit_destroy(fmpz& x, slimb fit, int s);
slimb _fmpz_normalize_value(fmpn_struct* X, slimb len, int s);
void _fmpz_demote(fmpz& x, slimb a);


struct fmpz_ulimb_realizer : ulimb_realizer {
    fmpz& _x;
    ulimb _fit;
    int _s;
    fmpz_ulimb_realizer(fmpz& x, ulimb fit, int s) : _x(x), _fit(fit), _s(s) {};
    ulimb* get_limbs() {return _fmpz_promote(_x, _fit, _s)->limbs;}
};


inline void _fmpz_demote(fmpz& x, ulimb a) {_fmpz_demote(x, (slimb) a);}
inline void fmpz_zero(fmpz& x) {_fmpz_demote(x, slimb(0));}
inline void fmpz_one(fmpz& x) {_fmpz_demote(x, slimb(1));}

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


void fmpz_mul_2exp(fmpz& x, fmpzc a, ulimb e);
inline void mul_2exp(fmpz_ring& ZZ, fmpz& x, const fmpz& a, ulimb e) {fmpz_mul_2exp(x, a, e);}

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

void fmpz_pow_ui(fmpz& a, fmpzc b, ulimb e);

void fmpz_fib_ui(fmpz& fn, ulimb n);


void fmpz_randtest(fmpz& f, rand_state& state, flint_bitcnt_t bits);
void fmpz_randtest_not_zero(fmpz& f, rand_state& state, flint_bitcnt_t bits);
void fmpz_randtest_mod(fmpz& f, rand_state& state, fmpzc m);

ulimb fmpz_bits(fmpzc a);

double fmpz_get_d(fmpzc a);
double fmpz_get_d_2exp(slimb* exp, fmpzc a);

inline bool fmpz_abs_fits_ui(fmpzc x) {return x.small() || x.length() < 2;}

ulimb fmpz_get_ui(fmpzc x);

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

inline void randtest(fmpz_ring& ZZ, fmpz& x, rand_state& state, flint_bitcnt_t bits) {fmpz_randtest(x, state, bits);}
inline void randtest_not_zero(fmpz_ring& ZZ, fmpz& x, rand_state& state, flint_bitcnt_t bits) {fmpz_randtest_not_zero(x, state, bits);}
inline void swap(fmpz_ring& ZZ, fmpz& x, fmpz& y) {std::swap(x.data, y.data);}
inline void set(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_set(x, a);}
inline void neg(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_neg(x, a);}
inline void neg(fmpz_ring& ZZ, fmpz& x) {fmpz_neg(x);}
inline void add(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_add(x, a, b);}
inline void sub(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_sub(x, a, b);}
inline void mul(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_mul(x, a, b);}
inline void sqr(fmpz_ring& ZZ, fmpz& x, fmpzc a) {fmpz_sqr(x, a);}
inline bool divisible(fmpz_ring& ZZ, fmpzc a, fmpzc b) {return fmpz_divisible(a, b);}
inline void divexact(fmpz_ring& ZZ, fmpz& x, fmpzc a, fmpzc b) {fmpz_divexact(x, a, b);}
inline void tdiv_qr(fmpz_ring& ZZ, fmpz& q, fmpz& r, fmpzc a, fmpzc b) {fmpz_tdiv_qr(q, r, a, b);}
inline void fdiv_qr(fmpz_ring& ZZ, fmpz& q, fmpz& r, fmpzc a, fmpzc b) {fmpz_fdiv_qr(q, r, a, b);}
inline void cdiv_qr(fmpz_ring& ZZ, fmpz& q, fmpz& r, fmpzc a, fmpzc b) {fmpz_cdiv_qr(q, r, a, b);}
inline void gcd(fmpz_ring& ZZ, fmpz& g, fmpzc a, fmpzc b) {fmpz_gcd(g, a, b);}
inline void gcd(fmpz_ring& ZZ, fmpz& g, fmpzc a, fmpzc b, fmpzc c) {fmpz_gcd(g, a, b, c);}
inline void gcdx(fmpz_ring& ZZ, fmpz& g, fmpz& s, fmpz& t, fmpzc a, fmpzc b) {fmpz_gcdx(g, s, t, a, b);}

inline fmpz operator -(fmpzc a) {fmpz x; fmpz_neg(x, a); return x;}
inline fmpz operator -(fmpz&& x) {fmpz_neg(x, x); return x;}
inline fmpz operator +(fmpzc a, fmpzc b) {fmpz x; fmpz_add(x, a, b); return x;}
inline fmpz operator -(fmpzc a, fmpzc b) {fmpz x; fmpz_sub(x, a, b); return x;}
inline fmpz operator *(fmpzc a, fmpzc b) {fmpz x; fmpz_mul(x, a, b); return x;}



#if 0
/*
        the fmpz
        11xxxxxx    small negative immediate
        11000000    large negative fmpn_struct* (this one is nullptr and not used)
-2^63:  10xxxxxx    large negative fmpn_struct*
2^63-1: 01xxxxxx    large positive fmpn_struct*
        00111111    small positive immediate
        00xxxxxx    small positive immediate

large negative <= -2^62, -(2^62-1) <= small <= 2^62-1, 2^62 <= large positive
*/

struct fmpn_struct {
    slimb length;
    slimb alloc;
    ulimb d[];
};

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

//s = 0 -> 01b  (positive)
//s = 1 -> 10b  (negative)
#define PTR_TO_COEFF(x, s) (((ulimb) (x) >> 2) | ((((ulimb)(1+(s))) << (FLINT_BITS-2))))
#define COEFF_TO_PTR(x) ((fmpn_struct*) ((x) << 2))
#define COEFF_MAX ((SWORD(1) << (FLINT_BITS-2)) - 1)
#define COEFF_MIN (-COEFF_MAX)
// these also could be !COEFF_IS_SMALL
#define COEFF_IS_PTR(x) ((x) < COEFF_MIN || COEFF_MAX < (x))
#define COEFF_IS_MPZ(x) ((x) < COEFF_MIN || COEFF_MAX < (x))
// nullptr is not used so 11000000 is an illegal fmpz
//#define COEFF_IS_SMALL(x) (COEFF_MIN-1 <= (x) && (x) <= COEFF_MAX)
//#define COEFF_IS_SMALL(x) ( ( ((ulimb)(x))^(2*((ulimb)(x))) ) < UWORD(1) << (FLINT_BITS - 1) )

// nullptr is not used so 11000000 is an illegal fmpz
inline bool _fmpz_is_small(slimb x)
{
    ulimb X = x;
    return (X^(X+X)) < UWORD(1) << (FLINT_BITS - 1);
}

// |x| < 2^(FLINT_BITS-2)
inline bool _word_is_small(slimb x)
{
    ulong X = (ulong)x + (UWORD(1) << (FLINT_BITS - 2));
    return 0 < X && X < UWORD(1) << (FLINT_BITS - 1);
}


typedef slimb fmpz;
typedef fmpz fmpz_t[1];

fmpn_struct* _fmpz_promote(fmpz_t x, slimb fit, int s);
fmpz _fmpz_normalize_value(fmpn_struct* X, slimb len, int s);

/* init **********************************************************************/
void fmpz_init_set(fmpz_t x, const fmpz_t a);
bool fmpz_is_canonical(const fmpz_t a);

/* io ************************************************************************/
std::ostream& operator<<(std::ostream& o, const fmpz_t a);
size_t fmpz_sizeinbase(const fmpz_t f, int b);
char * fmpz_get_str(char * str, int b, const fmpz_t f);

/* rand **********************************************************************/
void fmpz_randbits(fmpz_t f, rand_state& state, flint_bitcnt_t bits);
void fmpz_randm(fmpz_t f, rand_state& state, const fmpz_t m);
void fmpz_randtest(fmpz_t f, rand_state& state, flint_bitcnt_t bits);
void fmpz_randtest_unsigned(fmpz_t f, rand_state& state, flint_bitcnt_t bits);
void fmpz_randtest_not_zero(fmpz_t f, rand_state& state, flint_bitcnt_t bits);
void fmpz_randtest_mod(fmpz_t f, rand_state& state, const fmpz_t m);
void fmpz_randtest_mod_signed(fmpz_t f, rand_state& state, const fmpz_t m);
void fmpz_randprime(fmpz_t f, rand_state& state, flint_bitcnt_t bits, int proved);

/* get_set *******************************************************************/
slimb fmpz_get_si(const fmpz_t f);
ulimb fmpz_get_ui(const fmpz_t f);
ulimb fmpz_get_nmod(const fmpz_t f, nmod_t mod);
void fmpz_get_signed_uiui(ulimb * hi, ulimb * lo, const fmpz_t x);
void fmpz_set_signed_uiuiui(fmpz_t r, ulimb hi, ulimb mid, ulimb lo);
void fmpz_get_ui_array(ulimb * out, slimb n, const fmpz_t in);
void fmpz_set_ui_array(fmpz_t out, const ulimb * in, slimb n);
void fmpz_get_signed_ui_array(ulimb * out, slimb n, const fmpz_t in);
void fmpz_set_signed_ui_array(fmpz_t out, const ulimb * in, slimb n);
double fmpz_get_d(const fmpz_t f);
double fmpz_get_d_2exp(slimb * exp, const fmpz_t f);
void fmpz_set_d(fmpz_t f, double c);
void fmpz_set_d_2exp(fmpz_t f, double m, slimb exp);
int fmpz_set_str(fmpz_t f, const char * str, int b);
void fmpz_set(fmpz_t f, const fmpz_t g);

/* bits **********************************************************************/
int fmpz_abs_fits_ui(const fmpz_t f);
int fmpz_fits_si(const fmpz_t f);
slimb fmpz_size(const fmpz_t f);
flint_bitcnt_t fmpz_bits(const fmpz_t f);
flint_bitcnt_t fmpz_val2(const fmpz_t x);
int fmpz_bit_pack(mp_ptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits, const fmpz_t coeff, int negate, int borrow);
int fmpz_bit_unpack(fmpz_t coeff, mp_srcptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits, int negate, int borrow);
void fmpz_bit_unpack_unsigned(fmpz_t coeff, mp_srcptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits);
void fmpz_setbit(fmpz_t f, ulimb i);
int fmpz_tstbit(const fmpz_t f, ulimb i);
void fmpz_clrbit(fmpz_t f, ulimb i);
void fmpz_complement(fmpz_t r, const fmpz_t f);
void fmpz_combit(fmpz_t f, ulimb i);
void fmpz_and(fmpz_t r, const fmpz_t a, const fmpz_t b);
void fmpz_or(fmpz_t r, const fmpz_t a, const fmpz_t b);
void fmpz_xor(fmpz_t r, const fmpz_t a, const fmpz_t b);
flint_bitcnt_t fmpz_popcnt(const fmpz_t c);

/* cmp ***********************************************************************/
int fmpz_equal(const fmpz_t f, const fmpz_t g);
int fmpz_equal_si(const fmpz_t f, slimb g);
int fmpz_equal_ui(const fmpz_t f, ulimb g);
int fmpz_cmp(const fmpz_t f, const fmpz_t g);
int fmpz_cmp_ui(const fmpz_t f, ulimb g);
int fmpz_cmp_si(const fmpz_t f, slimb g);
int fmpz_cmpabs(const fmpz_t f, const fmpz_t g);
int fmpz_cmp2abs(const fmpz_t f, const fmpz_t g);

/* add_sub_mul ***************************************************************/
void fmpz_neg(fmpz_t f1, const fmpz_t f2);
void fmpz_abs(fmpz_t f1, const fmpz_t f2);
void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulimb x);
void fmpz_mul_si(fmpz_t f, const fmpz_t g, slimb x);
void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulimb exp);
void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulimb x);
void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulimb x);
void fmpz_addmul_si(fmpz_t f, const fmpz_t g, slimb x);
void fmpz_submul_si(fmpz_t f, const fmpz_t g, slimb x);
void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulimb x);
void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulimb x);
void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_fmma(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d);
void fmpz_fmms(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d);

/* pow_root ******************************************************************/
void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulimb e);
int fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e);
void fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulimb exp, const fmpz_t m);
void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m);
int fmpz_is_square(const fmpz_t f);
int fmpz_root(fmpz_t r, const fmpz_t f, slimb n);
int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f);
void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g);

double fmpz_dlog(const fmpz_t x);
slimb fmpz_flog(const fmpz_t x, const fmpz_t b);
slimb fmpz_flog_ui(const fmpz_t x, ulimb b);
slimb fmpz_clog(const fmpz_t x, const fmpz_t b);
slimb fmpz_clog_ui(const fmpz_t x, ulimb b);

int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p);
void fmpz_sqrt(fmpz_t f, const fmpz_t g);

slimb _fmpz_remove(fmpz_t x, const fmpz_t f, double finv);
slimb fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f);

/* divexact ******************************************************************/
void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slimb h);
void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulimb h);
int fmpz_divisible(const fmpz_t f, const fmpz_t g);
int fmpz_divides(fmpz_t q, const fmpz_t g, const fmpz_t h);
int fmpz_divisible_si(const fmpz_t f, slimb g);

/* divrem ********************************************************************/
ulimb fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulimb h);
void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_smod(fmpz_t f, const fmpz_t g, const fmpz_t h);
void _fmpz_smod(fmpz_t r, const fmpz_t a, const fmpz_t m, int sign, fmpz_t t);

void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);
void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_fdiv_r(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulimb h);
void fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, slimb h);
void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulimb exp);
void fmpz_fdiv_r_2exp(fmpz_t f, const fmpz_t g, ulimb exp);
ulimb fmpz_fdiv_ui(const fmpz_t g, ulimb h);

void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);
void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, slimb h);
void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulimb h);
void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulimb exp);
void fmpz_cdiv_r_2exp(fmpz_t f, const fmpz_t g, ulimb exp);
ulimb fmpz_cdiv_ui(const fmpz_t g, ulimb h);

void fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h);
void fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, ulimb h);
void fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, slimb h);
void fmpz_tdiv_r_2exp(fmpz_t f, const fmpz_t g, ulimb exp);
ulimb fmpz_tdiv_ui(const fmpz_t g, ulimb h);
void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulimb exp);

void fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulimb exp);
void fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, slimb x, ulimb exp);

void fmpz_ndiv_qr(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b);

/* fac ***********************************************************************/
void fmpz_fac_ui(fmpz_t f, ulimb n);
void fmpz_fib_ui(fmpz_t f, ulimb n);
void fmpz_bin_uiui(fmpz_t res, ulimb n, ulimb k);
void _fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulimb a, ulimb b);
void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulimb n);
void fmpz_rfac_uiui(fmpz_t r, ulimb x, ulimb n);

/* gcd ***********************************************************************/
void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_gcd3(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c);
void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h);
void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g);
void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g);
void fmpz_xgcd_canonical_bezout(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g);
void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, fmpz_t r2, fmpz_t r1, const fmpz_t L);
int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h);
int fmpz_jacobi(const fmpz_t a, const fmpz_t p);
int fmpz_kronecker(const fmpz_t a, const fmpz_t n);
void fmpz_divides_mod_list(fmpz_t xstart, fmpz_t xstride, fmpz_t xlength, const fmpz_t a, const fmpz_t b, const fmpz_t n);

/* crt ***********************************************************************/
void _fmpz_CRT_ui_precomp(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulimb r2, ulimb m2, ulimb m2inv, const fmpz_t m1m2, ulimb c, int sign);
void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulimb r2, ulimb m2, int sign);
void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, fmpz_t r2, fmpz_t m2, int sign);

/* multi_CRT *****************************************************************/
typedef struct
{
    slimb a_idx; /* index of A */
    slimb b_idx; /* index of B */
    slimb c_idx; /* index of C */
    fmpz_t b_modulus;
    fmpz_t c_modulus;
} _fmpz_multi_CRT_instr;

typedef struct
{
    _fmpz_multi_CRT_instr * prog; /* straight line program */
    fmpz * moduli, * fracmoduli;
    fmpz_t final_modulus;
    slimb moduli_count;
    flint_bitcnt_t min_modulus_bits;
    slimb length; /* length of prog */
    slimb alloc;  /* alloc of prog */
    slimb localsize; /* length of outputs required in fmpz_multi_CRT_run */
    slimb temp1loc, temp2loc, temp3loc, temp4loc;
    int good;   /* the moduli are good for CRT, essentially relatively prime */
} fmpz_multi_CRT_struct;

typedef fmpz_multi_CRT_struct fmpz_multi_CRT_t[1];

void fmpz_multi_CRT_init(fmpz_multi_CRT_t CRT);
int fmpz_multi_CRT_precompute(fmpz_multi_CRT_t CRT, const fmpz * moduli, slimb len);
void fmpz_multi_CRT_precomp(fmpz_t output, const fmpz_multi_CRT_t P, const fmpz * inputs, int sign);
int fmpz_multi_CRT(fmpz_t output, const fmpz * moduli, const fmpz * values, slimb len, int sign);
void fmpz_multi_CRT_clear(fmpz_multi_CRT_t P);
void _fmpz_multi_CRT_precomp(fmpz * outputs, const fmpz_multi_CRT_t P, const fmpz * inputs, int sign);

#if 0
void fmpz_multi_crt_init(fmpz_multi_crt_t CRT);
int fmpz_multi_crt_precompute(fmpz_multi_crt_t CRT, const fmpz * moduli, slimb len);
int fmpz_multi_crt_precompute_p(fmpz_multi_crt_t CRT, const fmpz * const * moduli, slimb len);
void fmpz_multi_crt_precomp(fmpz_t output, const fmpz_multi_crt_t P, const fmpz * inputs);
void fmpz_multi_crt_precomp_p(fmpz_t output, const fmpz_multi_crt_t P, const fmpz * const * inputs);
int fmpz_multi_crt(fmpz_t output, const fmpz * moduli, const fmpz * values, slimb len);
void fmpz_multi_crt_clear(fmpz_multi_crt_t P);
inline slimb _fmpz_multi_crt_local_size(const fmpz_multi_crt_t CRT) {return CRT->localsize;}
void _fmpz_multi_crt_run(fmpz * outputs, const fmpz_multi_crt_t CRT, const fmpz * inputs);
void _fmpz_multi_crt_run_p(fmpz * outputs, const fmpz_multi_crt_t CRT, const fmpz * const * inputs);
#endif

/* multi_mod *****************************************************************/
typedef struct
{
    slimb in_idx;
    slimb out_idx;
    fmpz_t modulus;
} _fmpz_multi_mod_instr;

typedef struct
{
    _fmpz_multi_mod_instr * prog; /* straight line program */
    fmpz * moduli;
    slimb moduli_count;
    flint_bitcnt_t min_modulus_bits;
    slimb length; /* length of prog */
    slimb alloc;  /* alloc of prog */
    slimb localsize; /* length of tmp required in _fmpz_multi_mod_precomp */
    slimb temp1loc;
    int good;   /* the moduli are good for MOD, none are zero */
} fmpz_multi_mod_struct;

typedef fmpz_multi_mod_struct fmpz_multi_mod_t[1];

void fmpz_multi_mod_init(fmpz_multi_mod_t P);
void fmpz_multi_mod_clear(fmpz_multi_mod_t P);
int fmpz_multi_mod_precompute(fmpz_multi_mod_t P, const fmpz * f, slimb r);
void fmpz_multi_mod_precomp(fmpz * outputs, const fmpz_multi_mod_t P, const fmpz_t input, int sign);
void _fmpz_multi_mod_precomp(fmpz * outputs, const fmpz_multi_mod_t P, const fmpz_t input, int sign, fmpz * tmp);

/* multi_CRT_ui **************************************************************/
typedef struct {
    nmod_t mod;
    ulimb i0, i1, i2;
} crt_lut_entry;

typedef struct {
    nmod_t mod;
    nmod_t mod0, mod1, mod2;
} mod_lut_entry;

typedef struct {
    fmpz_multi_CRT_t crt_P;
    fmpz_multi_mod_t mod_P;
    ulimb * packed_multipliers;
    slimb * step;
    slimb * crt_offsets;
    slimb crt_offsets_alloc;
    slimb * mod_offsets;
    slimb mod_offsets_alloc;
    crt_lut_entry * crt_lu;
    slimb crt_lu_alloc;
    slimb crt_klen;
    mod_lut_entry * mod_lu;
    slimb mod_lu_alloc;
    slimb mod_klen;
    slimb num_primes;
} fmpz_comb_struct;

typedef fmpz_comb_struct fmpz_comb_t[1];

typedef struct {
    slimb Alen, Tlen;
    fmpz * A, * T;
} fmpz_comb_temp_struct;

typedef fmpz_comb_temp_struct fmpz_comb_temp_t[1];

void fmpz_comb_temp_init(fmpz_comb_temp_t CT, const fmpz_comb_t C);
void fmpz_comb_temp_clear(fmpz_comb_temp_t CT);
void fmpz_comb_init(fmpz_comb_t C, mp_srcptr primes, slimb num_primes);
void fmpz_comb_clear(fmpz_comb_t C);
void fmpz_multi_mod_ui(ulimb * out, const fmpz_t in, const fmpz_comb_t C, fmpz_comb_temp_t CT);
void fmpz_multi_CRT_ui(fmpz_t output, mp_srcptr residues, const fmpz_comb_t comb, fmpz_comb_temp_t temp, int sign);

/* prime *********************************************************************/
void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t m, const fmpz_t n);
void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t m, const fmpz_t n);
void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t n);
void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t Un, const fmpz_t Un1, const fmpz_t A, const fmpz_t B, const fmpz_t n);
void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t k, const fmpz_t n);
void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1, const fmpz_t Vm, const fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t Dinv, const fmpz_t n);

int fmpz_is_probabprime_lucas(const fmpz_t n);
int fmpz_is_probabprime_BPSW(const fmpz_t n);
int fmpz_is_strong_probabprime(const fmpz_t n, const fmpz_t a);
int fmpz_is_probabprime(const fmpz_t p);
int fmpz_is_prime_pseudosquare(const fmpz_t n);
void _fmpz_nm1_trial_factors(const fmpz_t n, mp_ptr pm1, slimb * num_pm1, ulimb limit);

int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pm1, slimb num_pm1);
void _fmpz_np1_trial_factors(const fmpz_t n, mp_ptr pp1, slimb * num_pp1, ulimb limit);
int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pm1, slimb num_pm1);
int fmpz_is_prime(const fmpz_t p);
int fmpz_divisor_in_residue_class_lenstra(fmpz_t fac, const fmpz_t n, const fmpz_t r, const fmpz_t s);
void fmpz_nextprime(fmpz_t res, const fmpz_t n, int proved);
void fmpz_primorial(fmpz_t res, ulimb n);
void fmpz_euler_phi(fmpz_t res, const fmpz_t n);
int fmpz_moebius_mu(const fmpz_t n);
void fmpz_divisor_sigma(fmpz_t res, const fmpz_t n, ulimb k);

/* misc **********************************************************************/
ulimb fmpz_abs_ubound_ui_2exp(slimb * exp, const fmpz_t x, int bits);
ulimb fmpz_abs_lbound_ui_2exp(slimb * exp, const fmpz_t x, int bits);
ulimb nmod_pow_fmpz(ulimb a, const fmpz_t exp, nmod_t mod);
ulimb n_powmod2_fmpz_preinv(ulimb a, const fmpz_t exp, ulimb n, ulimb ninv);


/**************************** inlines ****************************************/

inline void fmpz_clear(fmpz_t x)
{
    if (COEFF_IS_PTR(*x))
        flint_free(COEFF_TO_PTR(*x));
}

inline void fmpz_init(fmpz_t x)
{
    *x = 0;
}

inline void fmpz_init_set_si(fmpz_t x, slimb a)
{
    if (COEFF_MIN <= a && a <= COEFF_MAX)
    {
        *x = a;
    }
    else
    {
        fmpn_struct* X = _fmpn_alloc(1);
        X->length = 1;
        X->d[0] = (a > 0) ? a : -a;
        *x = PTR_TO_COEFF(X, a < 0 ? 1 : 0);
    }
}

// ensure x is small
inline void _fmpz_demote(fmpz_t x)
{
    if (!_fmpz_is_small(*x))
    {
        free(COEFF_TO_PTR(*x));
        *x = 0;
    }
}

inline int fmpz_is_even(const fmpz_t f)
{
    if (!COEFF_IS_PTR(*f))
        return !((*f) & WORD(1));
    else
        return !((COEFF_TO_PTR(*f)->d[0]) & UWORD(1));
}

inline int fmpz_is_odd(const fmpz_t f)
{
    if (!COEFF_IS_PTR(*f))
        return ((*f) & WORD(1));
    else
        return ((COEFF_TO_PTR(*f)->d[0]) & UWORD(1));
}

inline void fmpz_swap(fmpz_t x, fmpz_t y)
{
    fmpz t = *y;
    *y = *x;
    *x = t;
}

inline int fmpz_is_zero(const fmpz_t x)
{
    return *x == 0;
}

inline int fmpz_is_one(const fmpz_t x)
{
    return *x == 1;
}

inline int fmpz_sgn(const fmpz_t a)
{
    return *a > 0 ? 1 : *a < 0 ? -1 : 0;
}

inline int fmpz_is_pm1(const fmpz_t a)
{
    return *a == 1 || *a == -1;
}

inline void fmpz_zero(fmpz_t x)
{
    _fmpz_demote(x);
    *x = 0;
}

inline void fmpz_one(fmpz_t x)
{
    _fmpz_demote(x);
    *x = 1;
}

inline void fmpz_abs(fmpz_t x, const fmpz_t a)
{
    if (*a < 0)
        fmpz_neg(x, a);
    else
        fmpz_set(x, a);
}

inline void fmpz_add_si(fmpz_t z, const fmpz_t x, slimb y)
{
    if (y >= 0)
        fmpz_add_ui(z, x, y);
    else
        fmpz_sub_ui(z, x, -y);
}

inline void fmpz_sub_si(fmpz_t z, const fmpz_t x, slimb y)
{
    if (y >= 0)
        fmpz_sub_ui(z, x, y);
    else
        fmpz_add_ui(z, x, -y);
}

inline void fmpz_set_si(fmpz_t x, slimb a)
{
    if (LIKELY(_word_is_small(a)))
    {
        _fmpz_demote(x);
        *x = a;
    }
    else
    {
        fmpn_struct* X = _fmpz_promote(x, 2, a < 0);
        X->d[0] = FLINT_ABS(a);
        X->length = 1;
    }
}

inline void fmpz_set_ui(fmpz_t x, ulimb a)
{
    if (LIKELY(a <= COEFF_MAX))
    {
        _fmpz_demote(x);
        *x = (slimb) a;
    }
    else
    {
        fmpn_struct* X = _fmpz_promote(x, 1, 0);
        X->d[0] = a;
        X->length = 1;
    }
}

inline void fmpz_set_pm_ui(fmpz_t x, ulimb a, int s)
{
    FLINT_ASSERT(s == 0 || s == 1);

    if (LIKELY(a <= COEFF_MAX))
    {
        _fmpz_demote(x);
        *x = s ? -(slimb) a : (slimb) a;
    }
    else
    {
        fmpn_struct* X = _fmpz_promote(x, 1, s);
        X->d[0] = a;
        X->length = 1;
    }
}

inline void fmpz_neg_ui(fmpz_t x, ulimb a)
{
    if (a > COEFF_MAX)
    {
        fmpn_struct* X = _fmpz_promote(x, 1, 1);
        X->d[0] = a;
        X->length = 1;
    }
    else
    {
        _fmpz_demote(x);
        *x = -(slimb) a;
    }
}

inline void fmpz_set_uiui(fmpz_t x, ulimb hi, ulimb lo)
{
    if (hi == 0)
    {
        fmpz_set_ui(x, lo);
    }
    else
    {
        fmpn_struct* X = _fmpz_promote(x, 2, 0);
        X->d[0] = lo;
        X->d[1] = hi;
        X->length = 2;
    }
}

inline void fmpz_neg_uiui(fmpz_t x, ulimb hi, ulimb lo)
{
    if (hi == 0)
    {
        fmpz_neg_ui(x, lo);
    }
    else
    {
        fmpn_struct* X = _fmpz_promote(x, 2, 1);
        X->d[0] = lo;
        X->d[1] = hi;
        X->length = 2;
    }
}

inline void fmpz_set_signed_uiui(fmpz_t x, ulimb hi, ulimb lo)
{
    if (LIKELY(SIGN_EXT(lo) == hi))
    {
        fmpz_set_si(x, (slimb)lo);
        return;
    }

    ulimb s = SIGN_EXT(hi);
    fmpn_struct* X = _fmpz_promote(x, 2, -s);
    SUB_DDMMSS(X->d[1], X->d[0], hi^s, lo^s, s, s);
    X->length = 1 + (X->d[1] != 0);
}

inline int fmpz_abs_fits_ui(const fmpz_t a)
{
    return (!COEFF_IS_PTR(*a)) || COEFF_TO_PTR(*a)->length <= 1;
}

inline int fmpz_fits_si(const fmpz_t a)
{
    return (!COEFF_IS_PTR(*a)) || (COEFF_TO_PTR(*a)->length == 1 && COEFF_TO_PTR(*a)->d[0] <= SWORD_MAX);
}


inline void fmpz_mul2_uiui(fmpz_t f, const fmpz_t g, ulimb h1, ulimb h2)
{
    ulimb hi, lo;
    umul_ppmm(hi, lo, h1, h2);
    if (!hi)
    {
        fmpz_mul_ui(f, g, lo);
    }
    else
    {
        fmpz_mul_ui(f, g, h1);
        fmpz_mul_ui(f, f, h2);
    }
}

inline void fmpz_divexact2_uiui(fmpz_t f, const fmpz_t g, ulimb h1, ulimb h2)
{
    ulimb hi, lo;
    umul_ppmm(hi, lo, h1, h2);
    if (hi == 0)
    {
        fmpz_divexact_ui(f, g, lo);
    }
    else
    {
        fmpz_divexact_ui(f, g, h1);
        fmpz_divexact_ui(f, f, h2);
    }
}


inline void fmpz_get_uiui(ulimb * hi, ulimb * low, const fmpz_t a)
{
    FLINT_ASSERT(*a >= 0);
    if (!COEFF_IS_MPZ(*a))
    {
        *low = *a;
        *hi  = 0;
    }
    else
    {
        fmpn_struct * A = COEFF_TO_PTR(*a);
        FLINT_ASSERT(A->length <= 2);
        *low = A->d[0];
        *hi  = A->length < 2 ? 0 : A->d[1];
    }
}

inline void fmpz_negmod(fmpz_t r, const fmpz_t a, const fmpz_t mod)
{
   if (fmpz_is_zero(a))
      fmpz_zero(r);
   else
      fmpz_sub(r, mod, a);
}

inline void fmpz_set_ui_smod(fmpz_t f, ulimb x, ulimb m)
{
    FLINT_ASSERT(x < m);
    if (x <= m / 2)
        fmpz_set_ui(f, x);
    else
        fmpz_neg_ui(f, m - x);
}




#if 0
FMPZ_INLINE
void flint_mpz_add_uiui(mpz_ptr a, mpz_srcptr b, ulimb c1, ulimb c0)
{
    ulimb d[2];
    mpz_t c;
    d[0] = c0;
    d[1] = c1;
    c->_mp_d = d;
    c->_mp_alloc = 2;
    c->_mp_size = d[1] != 0 ? 2 : d[0] != 0;
    mpz_add(a, b, c);
}

FMPZ_INLINE
void flint_mpz_add_signed_uiui(mpz_ptr a, mpz_srcptr b, ulimb c1, ulimb c0)
{
    ulimb d[2];
    ulimb c2 = FLINT_SIGN_EXT(c1);
    mpz_t c;
    sub_ddmmss(d[1], d[0], c2^c1, c2^c0, c2, c2);
    c->_mp_d = d;
    c->_mp_alloc = 2;
    c->_mp_size = d[1] != 0 ? 2 : d[0] != 0;
    if (c2 != 0)
        c->_mp_size = -c->_mp_size;
    mpz_add(a, b, c);
}

FMPZ_INLINE
void flint_mpz_add_uiuiui(mpz_ptr a, mpz_srcptr b, ulimb c2, ulimb c1, ulimb c0)
{
    ulimb d[3];
    mpz_t c;
    d[0] = c0;
    d[1] = c1;
    d[2] = c2;
    c->_mp_d = d;
    c->_mp_alloc = 3;
    c->_mp_size = d[2] != 0 ? 3 : d[1] != 0 ? 2 : d[0] != 0;
    mpz_add(a, b, c);
}
#endif

#endif
