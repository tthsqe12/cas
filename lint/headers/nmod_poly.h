#pragma once

#include "types.h"
#include "ulimb_extras.h"
#include "fmpz.h"
#include "generics.h"

// length-less vector
template <class T>
struct chunk {
    T* _data;
    size_t _alloc;

    T& operator[](ulimb i) const {return _data[i];}
    T* data() const {return _data;}
    size_t alloc() const {return _alloc;}

    chunk() : _data(nullptr), _alloc(0) {};

    ~chunk() {
        for (ulimb i = _alloc; i > 0; i--)
            _data[i-1].~T();
        std::free(_data);
    }

    chunk(const chunk<T>& other) {
        _alloc = other._alloc;
        _data = reinterpret_cast<T*>(std::malloc(_alloc*sizeof(T)));
        for (ulimb i = 0; i < _alloc; i++)
            new (&_data[i]) T(other._data[i]);
    }

    chunk(chunk<T>&& other) {
        _data = other._data;
        _alloc = other._alloc;
        other._data = nullptr;
        other._alloc = 0;
    }

	chunk<T>& operator=(chunk<T>&& other)
    {
        std::swap(_data, other._data);
        std::swap(_alloc, other._alloc);
        return *this;
    }

    T* fit_alloc(size_t n) {
        if (n <= _alloc)
            return _data;
        size_t newalloc = std::max(n, _alloc + _alloc/2);
        T* newdata = reinterpret_cast<T*>(std::realloc(_data, newalloc*sizeof(T)));
        for (ulimb i = _alloc; i < newalloc; i++)
            new (&newdata[i]) T();
        _alloc = newalloc;
        _data = newdata;
        return _data;
    }

    // TODO free + malloc for trivial types
    T* fit_alloc_destroy(size_t n) {
        if (n <= _alloc)
            return _data;
        size_t newalloc = std::max(n, _alloc + _alloc/2);
        T* newdata = reinterpret_cast<T*>(std::realloc(_data, newalloc*sizeof(T)));
        for (ulimb i = _alloc; i < newalloc; i++)
            new (&newdata[i]) T();
        _alloc = newalloc;
        _data = newdata;
        return _data;
    }
};

struct nmod_poly_product;

struct nmod_poly {
    typedef const nmod_poly& source_t;
    typedef nmod_poly_product product_t;

    ulimb len = 0;
    chunk<ulimb> coeffs;

    bool is_constant() const {return len <= 1;}
    slimb degree() const {return slimb(len)-1;}
    ulimb length() const {return len;}
    ulimb* data() const {return coeffs.data();}
    ulimb* fit_alloc(ulimb n, ulimb N) {return coeffs.fit_alloc(n*N);}
    ulimb* fit_alloc_destroy(ulimb n, ulimb N) {return coeffs.fit_alloc_destroy(n*N);}
    void set_length(ulimb n) {len = n;}
    void normalize_length(ulimb len, ulimb N);
};

struct nmod_bpoly {
    ulimb len = 0;
    chunk<nmod_poly> coeffs;

    ulimb length() const {return len;}
    void set_length(ulimb n) {len = n;}

    nmod_poly* fit_alloc(ulimb n) {return coeffs.fit_alloc(n);}
    nmod_poly* fit_alloc_destroy(ulimb n) {return coeffs.fit_alloc(n);}
};

struct nmod_poly_product {
    nmod_bpoly bases;
    chunk<ulimb> exps;

    ulimb length() const {return bases.length();}
    void set_length(ulimb n) {bases.set_length(n);}
    nmod_poly& base(ulimb i) const {return bases.coeffs.data()[i];}
    ulimb& exp(ulimb i) const {return exps[i];}

    void fit_alloc(ulimb n) {exps.fit_alloc(n); bases.fit_alloc(n);}
    void fit_alloc_destroy(ulimb n) {exps.fit_alloc_destroy(n); bases.fit_alloc_destroy(n);}
};

// immutable thread safe base
struct nmod_ring_base {
    fmpzc modulus;
    ulimb stride;
    const ulimb* modulus_limbs;

    fmpzc characteristic() {return modulus;}

    void complete_modulus_data() {
        FLINT_ASSERT_ALWAYS(modulus.cmp_small(1) > 0);
        stride = modulus.size();
        modulus_limbs = modulus.is_ptr() ? modulus.limbs() :
                                 reinterpret_cast<const ulimb*>(&modulus.data); // TODO dirty and bad
    }

    nmod_ring_base(fmpz&& m) : modulus(m.steal_data()) {
        complete_modulus_data();
    }

    nmod_ring_base(fmpzc m) : modulus(m.copy_abs()) {
        complete_modulus_data();
    }

    nmod_ring_base(const nmod_ring_base& R) : modulus(R.modulus.copy_abs()) {
        complete_modulus_data();
    }

    ~nmod_ring_base() {
        modulus.clear();
    }
};

// mutable non thread safe
struct nmod_ring : nmod_ring_base {
    typedef fmpz elem_t;
    typedef nmod_poly poly_t;
    typedef nmod_poly_product poly_product_t;

    ulimb* tmp8N;
    ulimb* mul_tmp2N;
    ulimb* mul_tmp2Np1;
    ulimb* mul_tmpNp1;

    void complete_tmp_data() {
        ulimb N = stride;
        ulimb* p = my_alloc<ulimb>(8*N + 2*N + 2*N+1 + N+1);
        tmp8N       = p; p += 8*N;
        mul_tmp2N   = p; p += 2*N;
        mul_tmp2Np1 = p; p += 2*N+1;
        mul_tmpNp1  = p; p += N+1;
    }

    nmod_ring(fmpz&& m) : nmod_ring_base(std::move(m)) {
        complete_tmp_data();
    }

    nmod_ring(fmpzc m) : nmod_ring_base(m) {
        complete_tmp_data();
    }

    ~nmod_ring() {
        my_free(tmp8N);
    }
};

inline void ui_vec_zero(ulimb* x, ulimb n)
{
    for (ulimb i = 0; i < n; i++)
        x[i] = 0;
}

inline void ui_vec_set_nz(ulimb* x, const ulimb* a, ulimb n)
{
    FLINT_ASSERT(n > 0);
    ulimb i = 0;
    do {
        x[i] = a[i];
    } while (++i < n);
}

inline void ui_vec_set(ulimb* x, const ulimb* a, ulimb n)
{
    if (n > 0)
        ui_vec_set_nz(x, a, n);
}

inline bool ui_vec_is_zero(const ulimb* a, ulimb N)
{
    for (ulimb i = 0; i < N; i++)
        if (a[i] != 0)
            return false;
    return true;
}

inline bool nmod_is_canonical(const nmod_ring_base& ctx, const ulimb* a)
{
    return mpn_cmp(a, ctx.modulus_limbs, ctx.stride) < 0;
}

inline bool nmod_is_zero(const nmod_ring_base& ctx, const ulimb* a)
{
    return ui_vec_is_zero(a, ctx.stride);
}

inline bool nmod_is_one(const nmod_ring_base& ctx, const ulimb* a)
{
    return a[0] == 1 && ui_vec_is_zero(a + 1, ctx.stride - 1);
}

inline void nmod_one(const nmod_ring_base& ctx, ulimb* x)
{
    x[0] = 1;
    ui_vec_zero(x + 1, ctx.stride - 1);
}

inline void nmod_set(const nmod_ring_base& ctx, ulimb* x, const ulimb* a)
{
    ui_vec_set(x, a, ctx.stride);
}

void nmod_set_ui(nmod_ring& ctx, ulimb* x, ulimb a);
void nmod_set_ui_array(nmod_ring& ctx, ulimb* x, const ulimb* a, ulimb an);
void nmod_neg(const nmod_ring_base& ctx, ulimb* x, const ulimb* a);
void nmod_set_fmpz(nmod_ring& ctx, ulimb* x, fmpzc a);
void nmod_add(const nmod_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b);
void nmod_sub(const nmod_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b);
void nmod_mul(nmod_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b);
void nmod_mul_ui(nmod_ring& ctx, ulimb* x, const ulimb* a, ulimb b);


bool nmod_poly_is_canonical(const nmod_ring_base& R, const nmod_poly& a);
std::ostream& nmod_poly_write(std::ostream& o, const nmod_ring_base& R, const nmod_poly& a, const char* var);

// input_output.cpp

void nmod_poly_set_strn(nmod_ring& ctx, nmod_poly& x, const char* s, ulimb sn);

inline void set(nmod_ring& ctx, nmod_poly& x, const char* s) {nmod_poly_set_strn(ctx, x, s, strlen(s));}

inline std::ostream& operator<<(std::ostream& o, const with<nmod_ring, nmod_poly>& e) {
    return nmod_poly_write(o, e.parent(), e.elem_const(), "x");
}

inline std::ostream& operator<<(std::ostream& o, const with<nmod_ring, fmpz>& e) {
    return fmpz_write(o, e.elem());
}

inline std::ostream& operator<<(std::ostream& o, const nmod_ring_base& R) {
    o << "ZZ/" << R.modulus << "ZZ";
    return o;
}


void nmod_poly_randtest(nmod_ring& R, nmod_poly& x, rand_state& state, ulimb len);

inline void randtest(nmod_ring& R, nmod_poly& x, rand_state& state, ulimb len) {nmod_poly_randtest(R, x, state, len);}

inline bool nmod_poly_is_zero(const nmod_ring_base& R, const nmod_poly& a) {return a.length() == 0;}
bool nmod_poly_equal(const nmod_ring_base& R, const nmod_poly& a, const nmod_poly& b);

inline bool is_zero(nmod_ring& R, const nmod_poly& a) {return nmod_poly_is_zero(R, a);}
inline bool equal(nmod_ring& R, const nmod_poly& a, const nmod_poly& b) {return nmod_poly_equal(R, a, b);}

// get_set.cpp
void nmod_poly_set(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a);
void nmod_poly_set_fmpz(nmod_ring& R, nmod_poly& x, fmpzc a);

inline void set(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_set(R, x, a);}
inline void set(nmod_ring& R, nmod_poly& x, fmpzc a) {nmod_poly_set_fmpz(R, x, a);}

inline void swap(nmod_ring& R, nmod_poly& x, nmod_poly& y) {std::swap(x, y);}

inline void nmod_poly_zero(const nmod_ring_base& R, nmod_poly& x) {x.set_length(0);}
void nmod_poly_one(const nmod_ring_base& R, nmod_poly& x);

inline void zero(nmod_ring& R, nmod_poly& x) {nmod_poly_zero(R, x);}
inline void one(nmod_ring& R, nmod_poly& x) {nmod_poly_one(R, x);}

void nmod_poly_gen(const nmod_ring_base& R, nmod_poly& x);

void nmod_poly_deflate_inplace(const nmod_ring_base& ctx, nmod_poly& x, ulimb e);
void nmod_poly_make_monic(nmod_ring& ctx, nmod_poly& x, const nmod_poly& a);
void nmod_poly_derivative(nmod_ring& ctx, nmod_poly& x, const nmod_poly& a);


// add_sub_mul_div.cpp
void nmod_poly_neg(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a);
void nmod_poly_add(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_sub(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_mul(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_sqr(nmod_ring& R, nmod_poly& x, const nmod_poly& a);
void nmod_poly_pow_ui(nmod_ring& R, nmod_poly& x, const nmod_poly& a, ulimb e);

void nmod_poly_divexact(nmod_ring& R, nmod_poly& q, const nmod_poly& a, const nmod_poly& b);
bool nmod_poly_divides(nmod_ring& R, nmod_poly& q, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_divrem(nmod_ring& R, nmod_poly& q, nmod_poly& r, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_gcd(nmod_ring& R, nmod_poly& g, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_gcdc(nmod_ring& R, nmod_poly& g, nmod_poly& abar, nmod_poly& bbar, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_gcdx(nmod_ring& R, nmod_poly& g, nmod_poly& s, nmod_poly& t, const nmod_poly& a, const nmod_poly& b);

/*****************************************************************************/

inline void deflate(nmod_ring& R, nmod_poly& x, ulimb e) {nmod_poly_deflate_inplace(R, x, e);}
inline void derivative(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_derivative(R, x, a);}
inline void unit_normalize(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_make_monic(R, x, a);}

inline void neg(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_neg(R, x, a);}
inline void neg(nmod_ring& R, nmod_poly& x) {nmod_poly_neg(R, x, x);}

inline void add(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_add(R, x, a, b);}
inline void add(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_add(R, x, x, a);}
inline void sub(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_sub(R, x, a, b);}
inline void sub(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_sub(R, x, x, a);}
inline void mul(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_mul(R, x, a, b);}
inline void mul(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly t(std::move(x)); nmod_poly_mul(R, x, t, a);}
inline void sqr(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_sqr(R, x, a);}

inline void pow(nmod_ring& R, nmod_poly& x, const nmod_poly& a, fmpzc e) {nmod_poly_pow_ui(R, x, a, fmpz_get_ui(e));}

inline void divexact(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_divexact(R, x, a, b);}
inline bool divides(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {return nmod_poly_divides(R, x, a, b);}
inline void divrem(nmod_ring& R, nmod_poly& q, nmod_poly& r, const nmod_poly& a, const nmod_poly& b) {nmod_poly_divrem(R, q, r, a, b);}
inline void gcd(nmod_ring& R, nmod_poly& g, const nmod_poly& a, const nmod_poly& b) {nmod_poly_gcd(R, g, a, b);}
inline void gcdc(nmod_ring& R, nmod_poly& g, nmod_poly& abar, nmod_poly& bbar, const nmod_poly& a, const nmod_poly& b) {nmod_poly_gcdc(R, g, abar, bbar, a, b);}
inline void gcdx(nmod_ring& R, nmod_poly& g, nmod_poly& s, nmod_poly& t, const nmod_poly& a, const nmod_poly& b) {nmod_poly_gcdx(R, g, s, t, a, b);}

// factor
std::ostream& nmod_poly_product_write(std::ostream& o, const nmod_ring_base& ctx, const nmod_poly_product& a, const char* var);

inline std::ostream& operator<<(std::ostream& o, const with<nmod_ring, nmod_poly_product>& e) {
    return nmod_poly_product_write(o, e.parent(), e.elem_const(), "x");
}

void nmod_poly_factor_squarefree(nmod_ring& ctx, nmod_poly_product& f, const nmod_poly& a);

inline void factor_squarefree(nmod_ring& ctx, nmod_poly_product& f, const nmod_poly& a) {nmod_poly_factor_squarefree(ctx, f, a);}
