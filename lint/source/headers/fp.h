#pragma once

#include "types.h"
#include "buffers.h"
#include "fmpz.h"

struct fp_poly;
struct fp_poly_poly;
struct fp_poly_product;
struct fp_mpoly;
struct fp_mpoly_ring;

// immutable thread safe base
struct fp_ring_base {
    ulimb _stride;
    fmpz _modulus;
    fmpz _modulus_half;
    const ulimb* modulus_limbs;
    ulimb _modulus_nbits;

    ulimb stride() const {return _stride;}

    fmpzc characteristic() const {return _modulus;}
    fmpzc characteristic_half() const {return _modulus_half;}
    ulimb characteristic_nbits() const {return _modulus_nbits;}
    fmpzc cardinality() const {return _modulus;}
    fmpzc cardinality_half() const {return _modulus_half;}
    ulimb cardinality_nbits() const {return _modulus_nbits;}

    void complete_modulus_data() {
        fmpz_abs(_modulus);
        FLINT_ASSERT_ALWAYS(_modulus.cmp_small(1) > 0);
        _stride = _modulus.size();
        modulus_limbs = _modulus.is_ptr() ? _modulus.limbs() :
                   reinterpret_cast<const ulimb*>(&_modulus.data); // TODO dirty, bad, and not trivially movable

        fmpz_tdiv_q_2exp(_modulus_half, _modulus, 1);
        _modulus_nbits = nbits(_modulus);
    }

    fp_ring_base(fmpz&& m) : _modulus(std::move(m)) {
        complete_modulus_data();
    }

    fp_ring_base(fmpzc m) {
        fmpz_abs(_modulus, m);
        complete_modulus_data();
    }

    void reset_modulus(fmpz&& m) {
        _modulus = std::move(m);
        complete_modulus_data();        
    }

    void reset_modulus(fmpzc m) {
        _modulus = m.copy();
        complete_modulus_data();        
    }
};

struct fp_tmp_elem;
struct fp_dotter;
struct fp_divider;
struct _fp_poly;
struct _fp_poly_mat22;

// mutable non thread safe
struct fp_ring : fp_ring_base {
    typedef fmpz elem_t;
    typedef ulimb coeff_t;
    typedef fp_tmp_elem tmp_elem_t;

    typedef fp_dotter dotter_t;
    typedef fp_divider divider_t;

    typedef fp_poly poly_t;
    typedef fp_poly_product poly_product_t;
    typedef fp_mpoly mpoly_t;
    typedef fp_mpoly_ring mpoly_ring_t;

    typedef _fp_poly _poly_t;
    typedef _fp_poly_mat22 _poly_mat22_t;

    static constexpr int tmp_elem_count = 13;

    ulimb* tmp_ptr;
    ulimb* tmp_base;
    ulimb* mul_tmp2N;
    ulimb* mul_tmp2Np1;
    ulimb* mul_tmpNp2;

    ~fp_ring() {
        FLINT_ASSERT_ALWAYS(tmp_ptr == tmp_base);
        my_free(tmp_base);
    }

    void complete_tmp_data() {
        ulimb N = stride();
        ulimb* p = my_alloc<ulimb>(tmp_elem_count*N + 2*N + 2*N+1 + N+2);
        tmp_ptr = p;
        tmp_base    = p; p += tmp_elem_count*N;
        mul_tmp2N   = p; p += 2*N;
        mul_tmp2Np1 = p; p += 2*N+1;
        mul_tmpNp2  = p; p += N+2;
    }

    fp_ring(fmpz&& m) : fp_ring_base(std::move(m)) {
        complete_tmp_data();
    }

    fp_ring(fmpzc m) : fp_ring_base(m) {
        complete_tmp_data();
    }

    ulimb* alloc_tmp_space(ulimb n) {
        ulimb* res = tmp_ptr;
        tmp_ptr += n;
        FLINT_ASSERT(tmp_base < tmp_ptr && tmp_ptr <= tmp_base + tmp_elem_count*stride());
        return res;
    }

    void free_tmp_space(ulimb n) {
        tmp_ptr -= n;
        FLINT_ASSERT(tmp_base <= tmp_ptr && tmp_ptr < tmp_base + tmp_elem_count*stride());
    }

    void set_modulus(fmpzc m) {
        FLINT_ASSERT_ALWAYS(tmp_base == tmp_ptr);
        reset_modulus(m);
        my_free(tmp_base);
        complete_tmp_data();
    }
};

struct fp_tmp_elem {
    fp_ring* _ctx;
    ulimb* _data;
    ulimb* data() {return _data;}
    fp_tmp_elem(fp_ring& R) : _ctx(&R) {_data = _ctx->alloc_tmp_space(_ctx->stride());}
    ~fp_tmp_elem() {_ctx->free_tmp_space(_ctx->stride());}
};

//// coeffs

std::string monomial_format(
    const ulimb* a,
    ulimb M);

inline void ui_vec_zero(ulimb* x, ulimb n)
{
    for (ulimb i = 0; i < n; i++)
        x[i] = 0;
}

inline void ui_vec_swap(ulimb* x, ulimb* y, ulimb n)
{
    for (ulimb i = 0; i < n; i++)
        std::swap(x[i], y[i]);
}

inline bool ui_vec_equal(const ulimb* a, const ulimb* b, ulimb n)
{
    for (ulimb i = 0; i < n; i++) {
        if (a[i] != b[i])
            return false;
    }
    return true;
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

bool fp_is_canonical(const fp_ring_base& ctx, const ulimb* a);

inline void fp_zero(const fp_ring_base& ctx, ulimb* x)
{
    return ui_vec_zero(x, ctx.stride());
}

inline bool fp_is_zero(const fp_ring_base& ctx, const ulimb* a)
{
    return ui_vec_is_zero(a, ctx.stride());
}

inline bool fp_is_one(const fp_ring_base& ctx, const ulimb* a)
{
    return a[0] == 1 && ui_vec_is_zero(a + 1, ctx.stride() - 1);
}

inline void fp_one(const fp_ring_base& ctx, ulimb* x)
{
    x[0] = 1;
    ui_vec_zero(x + 1, ctx.stride() - 1);
}

inline bool fp_equal(const fp_ring_base& ctx, const ulimb* a, const ulimb* b)
{
    return ui_vec_equal(a, b, ctx.stride());
}

inline void fp_set(const fp_ring_base& ctx, ulimb* x, const ulimb* a)
{
    ui_vec_set(x, a, ctx.stride());
}

bool fp_equal_si(fp_ring& ctx, const ulimb* a, slimb b);
bool fp_equal_ui(fp_ring& ctx, const ulimb* a, ulimb b);

void fp_set_ui(fp_ring& ctx, ulimb* x, ulimb a);
void fp_set_ui_array(fp_ring& ctx, ulimb* x, const ulimb* a, ulimb an);
void fp_neg(const fp_ring_base& ctx, ulimb* x, const ulimb* a);
void fp_set_fmpz(fp_ring& ctx, ulimb* x, fmpzc a);
void fp_add(const fp_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b);
void fp_sub(const fp_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b);
void fp_add_ui(const fp_ring_base& ctx, ulimb* x, const ulimb* a, ulimb b);
void fp_sub_ui(const fp_ring_base& ctx, ulimb* x, const ulimb* a, ulimb b);
void fp_mul(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b);
void fp_mul_ui(fp_ring& ctx, ulimb* x, const ulimb* a, ulimb b);
void fp_inv(fp_ring& ctx, ulimb* x, const ulimb* a);

void fp_dot_mul(fp_ring& ctx, const ulimb* a, const ulimb* b);
void fp_dot_madd(fp_ring& ctx, const ulimb* a, const ulimb* b);
void fp_dot_reduce(fp_ring& ctx, ulimb* x);

struct fp_dotter {
    fp_ring* _ctx;
    ulimb* _data;

    fp_dotter(fp_ring& R) : _ctx(&R) {
        ulimb N = _ctx->stride();
        _data = _ctx->alloc_tmp_space(2*N+1);
    }

    ~fp_dotter() {
        ulimb N = _ctx->stride();
        _ctx->free_tmp_space(2*N+1);
    }

    void zero(fp_ring& ctx) {
        ulimb N = ctx.stride();
        ui_vec_zero(_data, 2*N + 1);
    }

    void mul(fp_ring& ctx, const ulimb* a, const ulimb* b) {
        ulimb N = ctx.stride();
        _data[2*N] = 0;
        my_mpn_mul_n(_data, a, b, N);
    }

    void madd(fp_ring& ctx, const ulimb* a, const ulimb* b) {
        ulimb N = ctx.stride();
        my_mpn_mul_n(ctx.mul_tmp2N, a, b, N);
        _data[2*N] += mpn_add_n(_data, _data, ctx.mul_tmp2N, 2*N);
    }

    void add(fp_ring& ctx, const ulimb* a) {
        ulimb N = ctx.stride();
        _data[2*N] += mpn_add(_data, _data, 2*N, a, N);
    }

    void sub(fp_ring& ctx, const ulimb* a) {
        ulimb N = ctx.stride();
        ulimb* s = ctx.mul_tmpNp2;
        mpn_sub_n(s, ctx.modulus_limbs, a, N);
        _data[2*N] += mpn_add(_data, _data, 2*N, s, N);
    }

    void reduce(fp_ring& ctx, ulimb* x) {
        ulimb N = ctx.stride();
        my_mpn_tdiv_qr(ctx.mul_tmpNp2, x, _data, 2*N+1, ctx.modulus_limbs, N);
    }
};

struct fp_divider {
    fp_ring* _ctx;
    ulimb* _data;

    fp_divider(fp_ring& ctx) : _ctx(&ctx) {
        _data = _ctx->alloc_tmp_space(_ctx->stride());
    }

    // inv can throw, so not in a constructor that can leak tmp space
    void set_divisor(fp_ring& ctx, const ulimb* a, bool neg = false)  {
        fp_inv(ctx, _data, a);
        if (neg)
            fp_neg(ctx, _data, _data);
    }

    ~fp_divider() {
        _ctx->free_tmp_space(_ctx->stride());
    }

    bool divides(fp_ring& ctx, ulimb* x, const ulimb* a) {
        fp_mul(ctx, x, a, _data);
        return true;
    }
};

void fp_inv(fp_ring& ctx, ulimb* x, const ulimb* a);

void fp_random(random_state& state, fp_ring& ctx, ulimb* x);

void fp_vec_scalar_mul(fp_ring& ctx, ulimb* X, const ulimb* A, ulimb n, const ulimb* b);

inline void fp_vec_zero(fp_ring& ctx, ulimb* x, ulimb n) {
    ui_vec_zero(x, n*ctx.stride());
}

inline void fp_vec_set(fp_ring& ctx, ulimb* x, const ulimb* a, ulimb n) {
    ui_vec_set(x, a, n*ctx.stride());
}

inline ulimb fp_vec_normalized_length(fp_ring& ctx, const ulimb* a, ulimb n) {
    ulimb N = ctx.stride();
    while (n > 0 && ui_vec_is_zero(a + N*(n - 1), N))
        n--;
    return n;
}

void fp_vec_add(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b, ulimb n);
void fp_vec_sub(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b, ulimb n);
void fp_vec_neg(fp_ring& ctx, ulimb* x, const ulimb* a, ulimb n);

inline std::ostream& fp_write(std::ostream& o, const fp_ring_base& ctx, const ulimb* a) {return mpn_write(o, a, ctx.stride());}


inline bool is_zero(fp_ring& ctx, const ulimb* a) {return fp_is_zero(ctx, a);}
inline bool is_one(fp_ring& ctx, const ulimb* a) {return fp_is_one(ctx, a);}
inline void zero(fp_ring& ctx, ulimb* x) {fp_zero(ctx, x);}
inline void one(fp_ring& ctx, ulimb* x) {fp_one(ctx, x);}
inline void set(fp_ring& ctx, ulimb* x, const ulimb* a) {fp_set(ctx, x, a);}
inline void neg(fp_ring& ctx, ulimb* x, const ulimb* a) {fp_neg(ctx, x, a);}
inline void add(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b) {fp_add(ctx, x, a, b);}
inline void sub(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b) {fp_sub(ctx, x, a, b);}
inline void mul(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b) {fp_mul(ctx, x, a, b);}
inline void inv(fp_ring& ctx, ulimb* x, const ulimb* a) {fp_inv(ctx, x, a);}

inline void vec_set(fp_ring& ctx, ulimb*x, const ulimb* a, ulimb n) {fp_vec_set(ctx, x, a, n);}
inline void vec_neg(fp_ring& ctx, ulimb*x, const ulimb* a, ulimb n) {fp_vec_neg(ctx, x, a, n);}
inline void vec_add(fp_ring& ctx, ulimb*x, const ulimb* a, const ulimb* b, ulimb n) {fp_vec_add(ctx, x, a, b, n);}
inline void vec_sub(fp_ring& ctx, ulimb*x, const ulimb* a, const ulimb* b, ulimb n) {fp_vec_sub(ctx, x, a, b, n);}
inline void vec_scalar_mul(fp_ring& ctx, ulimb*x, const ulimb* a, ulimb n, const ulimb* b) {fp_vec_scalar_mul(ctx, x, a, n, b);}
inline ulimb vec_normalized_length(fp_ring& ctx, const ulimb* a, ulimb n) {return fp_vec_normalized_length(ctx, a, n);}

