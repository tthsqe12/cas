#pragma once

#include "types.h"
#include "buffers.h"
#include "fmpz.h"
#include "generic/generics.h"
#include "fp.h"

////////////////////////////////////////////////////////////////////////////////

struct fp_poly_poly;
struct fp_poly_product;

struct fp_poly : my_vec<ulimb> {
    typedef const fp_poly& source_t;
    typedef fp_poly_product product_t;
    typedef fp_poly modulus_t;
    typedef fp_poly argument_t;

    slimb degree() const {return slimb(length())-1;}
    ulimb* coeff(ulimb n, ulimb N) {return data() + N*n;}
    const ulimb* coeff(ulimb n, ulimb N) const {return data() + N*n;}
    void normalize_length(ulimb len, ulimb N);
};

struct fp_poly_poly : my_vec<fp_poly> {

    fp_poly& coeff(ulimb i) {return (*this)[i];}
    const fp_poly& coeff(ulimb i) const {return (*this)[i];}
};

struct fp_poly_product {
    fp_poly_poly bases;
    chunk<slimb> exps;

    ulimb length() const {return bases.length();}
    void set_length(ulimb n) {bases.set_length(n);}

    fp_poly& base(ulimb i) {return bases[i];}
    const fp_poly& base(ulimb i) const {return bases[i];}

    slimb& exp(slimb i) {return exps[i];}
    const slimb& exp(slimb i) const {return exps[i];}

    void fit_alloc(ulimb n) {exps.fit_alloc(n); bases.fit_alloc(n);}
    void fit_alloc_destroy(ulimb n) {exps.fit_alloc_destroy(n); bases.fit_alloc_destroy(n);}

    void append(const fp_poly& a, slimb e) {
        ulimb n = length();
        fit_alloc(n + 1);
        exp(n) = e;
        fp_poly acopy(a);
        std::swap(base(n), acopy);
        set_length(n+1);
    }
};

/*
    n = B->npoints is the number of points a_1, ..., a_n that have been added
    to the sequence. The polynomials A and S are then defined as

        A = x^n
        S = a_1*x^(n-1) + a_2*x^(n-2) + ... + a_n

    We maintain polynomials U0, V0, U1, V1 such that

        U0*A + V0*S = R0   deg(R0) >= n/2
        U1*A + V1*S = R1   deg(R1) < n/2

    where R0 and R1 are consecutive euclidean remainders and U0, V0, U1, V1 are
    the corresponding Bezout coefficients. Note that
        deg(U1) < deg(V1) = deg(A) - deg(R0) <= n/2

    The U0 and U1 are not stored explicitly. The points a_1, ..., a_n are stored
    in B->points, which is used merely as a resizable array.

    The main usage of this function is the rational reconstruction of a series

     a1    a2    a3         -U1
    --- + --- + --- + ... = ---- maybe
     x    x^2   x^3          V1

    It can be seen that

     a1    a2        an   -U1     R1
    --- + --- + ... --- = --- + -------
     x    x^2       x^n    V1   V1*x^n

    Thus the error is O(1/x^(n+1)) iff deg(R1) < deg(V1).
*/
struct fp_berlekamp_massey {
    slimb npoints;
    fp_poly R0, R1, V0, V1, qt;

    ulimb* points_data;
    ulimb points_alloc;
    ulimb points_length;
    void points_fit_alloc(ulimb l, ulimb N);
    // location of the i^th point
    ulimb* point(ulimb i, ulimb N) {return points_data + points_alloc - N*(i+1);}
    const ulimb* point(ulimb i, ulimb N) const {return points_data + points_alloc - N*(i+1);}

    fp_berlekamp_massey();
    fp_berlekamp_massey(fp_ring& ctx);
    ~fp_berlekamp_massey();

    void start_over(fp_ring& ctx);
    void add_zeros(fp_ring& ctx, ulimb count);
    void add_points(const fp_ring_base& ctx, const ulimb* a, ulimb count);
    void add_point(fp_ring& ctx, const ulimb* a);
    void add_point_ui(fp_ring& ctx, ulong a);
    bool reduce(fp_ring& ctx);

    const ulimb point_count() {return points_length;}
    const fp_poly& v_poly() {return V1;}
    const fp_poly& r_poly() {return R1;}

    std::ostream& write(std::ostream& o, const fp_ring_base& ctx, const char* var) const;
};


typedef fp_poly fp_poly_modulus;
typedef fp_poly fp_poly_argument;

struct _fp_poly {
    ulimb* coeffs;
    ulimb len;

    ulimb length() {return len;}
    ulimb* data() {return coeffs;}
    _fp_poly() : coeffs(nullptr), len(0) {}
    _fp_poly(const fp_poly& a) : coeffs((ulimb*)a.data()), len(a.length()) {}
};

struct _fp_poly_mat22 {
    _fp_poly _11;
    _fp_poly _12;
    _fp_poly _21;
    _fp_poly _22;
    int det;
};

////////////////////////////////////////////////////////////////////////////////

//// poly

void fp_poly_make_monic_with_lcinv(fp_ring& ctx, fp_poly& x, const fp_poly& a, const ulimb* alcinv);
void fp_poly_make_monic(fp_ring& ctx, fp_poly& x, const fp_poly& a);
bool fp_poly_is_monic(fp_ring& ctx, const fp_poly& a);

bool _fp_poly_is_canonical(fp_ring& ctx, _fp_poly& a);
bool fp_poly_is_canonical(const fp_ring_base& R, const fp_poly& a);

// input_output.cpp

std::ostream& _fp_poly_write(
    std::ostream& o,
    const fp_ring_base& ctx,
    _fp_poly a,
    const char* var);

inline std::ostream& fp_poly_write(
    std::ostream& o,
    const fp_ring_base& ctx,
    const fp_poly& a,
    const char* var)
{
    return _fp_poly_write(o, ctx, _fp_poly(a), var);
}

void fp_poly_set_strn(fp_ring& ctx, fp_poly& x, const char* s, ulimb sn);

inline void set(fp_ring& ctx, fp_poly& x, const char* s) {fp_poly_set_strn(ctx, x, s, strlen(s));}

inline std::ostream& operator<<(std::ostream& o, const with<fp_ring, fp_poly>& e) {
    return fp_poly_write(o, e.parent(), e.elem(), "x");
}

inline std::ostream& operator<<(std::ostream& o, const with<fp_ring, fmpz>& e) {
    return fmpz_write(o, e.elem());
}

inline std::ostream& operator<<(std::ostream& o, const fp_ring_base& R) {
    o << "ZZ/" << R.characteristic() << "ZZ";
    return o;
}

inline std::ostream& operator<<(std::ostream& o, const with<fp_ring, fp_berlekamp_massey>& e) {
    return e.elem().write(o, e.parent(), "x");
}



void fp_poly_random_of_length(random_state& state, fp_ring& R, fp_poly& x, ulimb len);
void fp_poly_random_up_to_length(random_state& state, fp_ring& R, fp_poly& x, ulimb len);
void fp_poly_random_monic_of_length(random_state& state, fp_ring& R, fp_poly& x, ulimb len);

inline void random_of_length(random_state& state, fp_ring& R, fp_poly& x, ulimb len) {fp_poly_random_of_length(state, R, x, len);}
inline void random_up_to_length(random_state& state, fp_ring& R, fp_poly& x, ulimb len) {fp_poly_random_up_to_length(state, R, x, len);}
inline void random_monic_of_length(random_state& state, fp_ring& R, fp_poly& x, ulimb len) {fp_poly_random_monic_of_length(state, R, x, len);}

inline bool fp_poly_is_zero(const fp_ring_base& R, const fp_poly& a) {return a.length() == 0;}
bool fp_poly_is_one(const fp_ring_base& R, const fp_poly& a);
bool fp_poly_equal(const fp_ring_base& R, const fp_poly& a, const fp_poly& b);
bool fp_poly_equal_si(fp_ring& ctx, const fp_poly& a, slimb b);
bool fp_poly_equal_ui(fp_ring& ctx, const fp_poly& a, ulimb b);

inline bool is_zero(fp_ring& R, const fp_poly& a) {return fp_poly_is_zero(R, a);}
inline bool is_one(fp_ring& R, const fp_poly& a) {return fp_poly_is_one(R, a);}
inline bool equal(fp_ring& R, const fp_poly& a, const fp_poly& b) {return fp_poly_equal(R, a, b);}
inline bool equal(fp_ring& R, const fp_poly& a, int b) {return fp_poly_equal_si(R, a, b);}
inline bool equal(fp_ring& R, const fp_poly& a, slimb b) {return fp_poly_equal_si(R, a, b);}
inline bool equal(fp_ring& R, const fp_poly& a, ulimb b) {return fp_poly_equal_ui(R, a, b);}

// get_set.cpp
void fp_poly_set(const fp_ring_base& R, fp_poly& x, const fp_poly& a);
void fp_poly_set_fmpz(fp_ring& R, fp_poly& x, fmpzc a);

inline void set(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_set(R, x, a);}
inline void set(fp_ring& R, fp_poly& x, fmpzc a) {fp_poly_set_fmpz(R, x, a);}

inline void swap(fp_ring& R, fp_poly& x, fp_poly& y) {std::swap(x, y);}

inline void fp_poly_zero(const fp_ring_base& R, fp_poly& x) {x.set_length(0);}
void fp_poly_one(const fp_ring_base& R, fp_poly& x);

inline void zero(fp_ring& R, fp_poly& x) {fp_poly_zero(R, x);}
inline void one(fp_ring& R, fp_poly& x) {fp_poly_one(R, x);}

void fp_poly_gen(const fp_ring_base& R, fp_poly& x);

void fp_poly_deflate_inplace(const fp_ring_base& ctx, fp_poly& x, ulimb e);
void fp_poly_make_monic(fp_ring& ctx, fp_poly& x, const fp_poly& a);
void fp_poly_derivative(fp_ring& ctx, fp_poly& x, const fp_poly& a);


// add_sub_mul_div.cpp
void fp_poly_neg(const fp_ring_base& R, fp_poly& x, const fp_poly& a);
void fp_poly_add(const fp_ring_base& R, fp_poly& x, const fp_poly& a, const fp_poly& b);
void fp_poly_sub(const fp_ring_base& R, fp_poly& x, const fp_poly& a, const fp_poly& b);
void fp_poly_mul(fp_ring& R, fp_poly& x, _fp_poly a, _fp_poly b);
void fp_poly_sqr(fp_ring& R, fp_poly& x, const fp_poly& a);
void fp_poly_pow_ui(fp_ring& R, fp_poly& x, const fp_poly& a, ulimb e);
void fp_poly_addmul(fp_ring& ctx, fp_poly& x, const fp_poly& a, const fp_poly& b);
void fp_poly_submul(fp_ring& ctx, fp_poly& x, const fp_poly& a, const fp_poly& b);

void fp_poly_add_inplace_ui(const fp_ring_base& R, fp_poly& x, ulimb a);
void fp_poly_sub_inplace_ui(const fp_ring_base& R, fp_poly& x, ulimb a);


void fp_poly_divexact(fp_ring& R, fp_poly& q, const fp_poly& a, const fp_poly& b);
bool fp_poly_divides(fp_ring& R, fp_poly& q, const fp_poly& a, const fp_poly& b);
void fp_poly_divrem(fp_ring& R, fp_poly& q, fp_poly& r, const fp_poly& a, const fp_poly& b);
void fp_poly_divrem_inplace(fp_ring& R, fp_poly& q, fp_poly& a, const fp_poly& b);
void fp_poly_mod(fp_ring& R, fp_poly& r, const fp_poly& a, const fp_poly& b);
void fp_poly_mulmod(fp_ring& ctx, fp_poly& x, const fp_poly& a, const fp_poly& b, const fp_poly& m);
void fp_poly_sqrmod(fp_ring& ctx, fp_poly& x, const fp_poly& a, const fp_poly& m);
void fp_poly_powmod_fmpz(fp_ring& ctx, fp_poly& x, const fp_poly& a, fmpzc e, const fp_poly& m);

// gcd
void fp_poly_gcd(fp_ring& R, fp_poly& g, const fp_poly& a, const fp_poly& b);
void fp_poly_gcdc(fp_ring& R, fp_poly& g, fp_poly& abar, fp_poly& bbar, const fp_poly& a, const fp_poly& b);
void fp_poly_gcdx(fp_ring& R, fp_poly& g, fp_poly& s, fp_poly& t, const fp_poly& a, const fp_poly& b);
void fp_poly_gcdinv(fp_ring& R, fp_poly& g, fp_poly& s, const fp_poly& a, const fp_poly& b);
int fp_poly_hgcd(fp_ring& ctx, fp_poly& m11, fp_poly& m12, fp_poly& m21, fp_poly& m22, fp_poly& A, fp_poly& B, const fp_poly& a, const fp_poly& b);


/*****************************************************************************/

inline void deflate(fp_ring& R, fp_poly& x, ulimb e) {fp_poly_deflate_inplace(R, x, e);}
inline void derivative(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_derivative(R, x, a);}
inline void unit_normalize(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_make_monic(R, x, a);}
inline void make_monic(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_make_monic(R, x, a);}

inline void neg(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_neg(R, x, a);}
inline void neg(fp_ring& R, fp_poly& x) {fp_poly_neg(R, x, x);}

inline void add(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b) {fp_poly_add(R, x, a, b);}
inline void add(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_add(R, x, x, a);}
inline void sub(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b) {fp_poly_sub(R, x, a, b);}
inline void sub(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_sub(R, x, x, a);}
inline void mul(fp_ring& R, fp_poly& x, _fp_poly a, _fp_poly b) {fp_poly_mul(R, x, a, b);}
inline void mul(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b) {fp_poly_mul(R, x, a, b);}
inline void mul(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly t(std::move(x)); fp_poly_mul(R, x, t, a);}
inline void sqr(fp_ring& R, fp_poly& x, const fp_poly& a) {fp_poly_sqr(R, x, a);}
inline void addmul(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b) {fp_poly_addmul(R, x, a, b);}
inline void submul(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b) {fp_poly_submul(R, x, a, b);}
inline void add(fp_ring& R, fp_poly& x, ulimb a) {fp_poly_add_inplace_ui(R, x, a);}
inline void sub(fp_ring& R, fp_poly& x, ulimb a) {fp_poly_sub_inplace_ui(R, x, a);}


void fp_poly_shift_left_inplace(fp_ring& ctx, fp_poly& x, ulimb e);
void fp_poly_shift_left(fp_ring& ctx, fp_poly& x, const fp_poly& a, ulimb e);
void fp_poly_shift_right_inplace(fp_ring& ctx, fp_poly& x, ulimb e);
void fp_poly_shift_right(fp_ring& ctx, fp_poly& x, const fp_poly& a, ulimb e);

inline void shift_left(fp_ring& R, fp_poly& x, ulimb e) {fp_poly_shift_left_inplace(R, x, e);}
inline void shift_left(fp_ring& R, fp_poly& x, const fp_poly& a, ulimb e) {fp_poly_shift_left(R, x, a, e);}
inline void shift_right(fp_ring& R, fp_poly& x, ulimb e) {fp_poly_shift_right_inplace(R, x, e);}
inline void shift_right(fp_ring& R, fp_poly& x, const fp_poly& a, ulimb e) {fp_poly_shift_right(R, x, a, e);}


inline void pow(fp_ring& R, fp_poly& x, const fp_poly& a, fmpzc e) {fp_poly_pow_ui(R, x, a, fmpz_abs_get_ui(e));}

inline void divexact(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b) {fp_poly_divexact(R, x, a, b);}
inline bool divides(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b) {return fp_poly_divides(R, x, a, b);}
inline void divrem(fp_ring& R, fp_poly& q, fp_poly& r, const fp_poly& a, const fp_poly& b) {fp_poly_divrem(R, q, r, a, b);}
inline void divrem(fp_ring& R, fp_poly& q, fp_poly& a, const fp_poly& b) {fp_poly_divrem_inplace(R, q, a, b);}
inline void mod(fp_ring& R, fp_poly& r, const fp_poly& a, const fp_poly& b) {fp_poly_mod(R, r, a, b);}

inline void mulmod(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& b, const fp_poly& m) {fp_poly_mulmod(R, x, a, b, m);}
inline void sqrmod(fp_ring& R, fp_poly& x, const fp_poly& a, const fp_poly& m) {fp_poly_sqrmod(R, x, a, m);}
inline void powmod(fp_ring& R, fp_poly& x, const fp_poly& a, fmpzc e, const fp_poly& m) {fp_poly_powmod_fmpz(R, x, a, e, m);}

inline void gcd(fp_ring& R, fp_poly& g, const fp_poly& a, const fp_poly& b) {fp_poly_gcd(R, g, a, b);}
inline void gcdc(fp_ring& R, fp_poly& g, fp_poly& abar, fp_poly& bbar, const fp_poly& a, const fp_poly& b) {fp_poly_gcdc(R, g, abar, bbar, a, b);}
inline void gcdx(fp_ring& R, fp_poly& g, fp_poly& s, fp_poly& t, const fp_poly& a, const fp_poly& b) {fp_poly_gcdx(R, g, s, t, a, b);}
inline void gcdinv(fp_ring& R, fp_poly& g, fp_poly& s, const fp_poly& a, const fp_poly& b) {fp_poly_gcdinv(R, g, s, a, b);}
inline int hgcd(fp_ring& ctx, fp_poly& m11, fp_poly& m12, fp_poly& m21, fp_poly& m22, fp_poly& A, fp_poly& B, const fp_poly& a, const fp_poly& b) {return fp_poly_hgcd(ctx, m11, m12, m21, m22, A, B, a, b);}

// factor

std::ostream& fp_poly_product_write(std::ostream& o, const fp_ring_base& ctx, const fp_poly_product& a, const char* var);
std::ostream& fp_poly_poly_write(std::ostream& o, const fp_ring_base& ctx, const fp_poly_poly& a, const char* var);

inline std::ostream& operator<<(std::ostream& o, const with<fp_ring, fp_poly_product>& e) {
    return fp_poly_product_write(o, e.parent(), e.elem(), "x");
}

inline std::ostream& operator<<(std::ostream& o, const with<fp_ring, fp_poly_poly>& e) {
    return fp_poly_poly_write(o, e.parent(), e.elem(), "x");
}

void fp_poly_factor_squarefree(fp_ring& ctx, fp_poly_product& f, const fp_poly& a);
void fp_poly_factor(fp_ring& ctx, fp_poly_product& f, const fp_poly& a);
bool fp_poly_is_squarefree(fp_ring& ctx, const fp_poly& f);
bool fp_poly_is_irreducible(fp_ring& ctx, const fp_poly& f);


inline bool is_monic(fp_ring& ctx, const fp_poly& a) {return fp_poly_is_monic(ctx, a);}
inline void factor_squarefree(fp_ring& ctx, fp_poly_product& f, const fp_poly& a) {fp_poly_factor_squarefree(ctx, f, a);}
inline void factor(fp_ring& ctx, fp_poly_product& f, const fp_poly& a) {fp_poly_factor(ctx, f, a);}
inline bool is_squarefree(fp_ring& ctx, const fp_poly& a) {return fp_poly_is_squarefree(ctx, a);}
inline bool is_irreducible(fp_ring& ctx, const fp_poly& a) {return fp_poly_is_irreducible(ctx, a);}


