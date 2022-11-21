#include "fp_poly.h"
#include "fmpz.h"

bool fp_is_canonical(const fp_ring_base& ctx, const ulimb* a)
{
    if (mpn_cmp(a, ctx.modulus_limbs, ctx.stride()) < 0)
        return true;

    std::cerr << "fp coeff is out of range" << std::endl;
    return false;
}

bool fp_equal_ui(fp_ring& ctx, const ulimb* a, ulimb b)
{
    ulimb N = ctx.stride();
    if (N > 1)
        return a[0] == b && ui_vec_is_zero(a + 1, N - 1);
    else
        return a[0] == b % ctx.modulus_limbs[0];
}

bool fp_equal_si(fp_ring& ctx, const ulimb* a, slimb b)
{
    ulimb N = ctx.stride();
    if (b >= 0)
        return fp_equal_ui(ctx, a, b);
    ulimb* t = ctx.mul_tmpNp2;
    fp_add_ui(ctx, t, a, -b);
    return fp_is_zero(ctx, t);
}

void fp_set_ui(fp_ring& ctx, ulimb* x, ulimb a)
{
    ulimb N = ctx.stride();

    if (N > 1)
    {
        x[0] = a;
        ui_vec_zero(x + 1, N - 1);
    }
    else
    {
        x[0] = a % ctx.modulus_limbs[0];
    }
}

void fp_set_ui_array(fp_ring& ctx, ulimb* x, const ulimb* a, ulimb an)
{
    ulimb N = ctx.stride();

    if (UNLIKELY(an >= 2*N))
    {
        tmp_allocator push;
        ulimb* t = push.recursive_alloc<ulimb>(an - N + 1);
        my_mpn_tdiv_qr(t, x, a, an, ctx.modulus_limbs, N);
    }
    else if (an >= N)
    {
        ulimb* t = ctx.mul_tmp2Np1;
        my_mpn_tdiv_qr(t, x, a, an, ctx.modulus_limbs, N);        
    }
    else
    {
        ui_vec_set(x, a, an);
        ui_vec_zero(x + an, N - an);
    }
}

void fp_neg(const fp_ring_base& ctx, ulimb* x, const ulimb* a)
{
    if (ui_vec_is_zero(a, ctx.stride()))
        ui_vec_zero(x, ctx.stride());
    else
        mpn_sub_n(x, ctx.modulus_limbs, a, ctx.stride());
}

void fp_set_fmpz(fp_ring& ctx, ulimb* x, fmpzc a)
{
    if (a.is_small())
        fp_set_ui(ctx, x, my_abs(a.small()));
    else
        fp_set_ui_array(ctx, x, a.limbs(), a.length());

    if (a.is_negative())
        fp_neg(ctx, x, x);
}

void fp_add(const fp_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    auto c = mpn_add_n(x, a, b, ctx.stride());
    if (c != 0 || mpn_cmp(x, ctx.modulus_limbs, ctx.stride()) >= 0)
        mpn_sub_n(x, x, ctx.modulus_limbs, ctx.stride());
}

void fp_add_ui(const fp_ring_base& ctx, ulimb* x, const ulimb* a, ulimb b)
{
    ulimb N = ctx.stride();
    if (N == 1 && b >= ctx.modulus_limbs[0])
        b = b%ctx.modulus_limbs[0];
    auto c = mpn_add_1(x, a, N, b);
    if (c != 0 || mpn_cmp(x, ctx.modulus_limbs, N) >= 0)
        mpn_sub_n(x, x, ctx.modulus_limbs, N);
}

void fp_sub(const fp_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    auto c = mpn_sub_n(x, a, b, ctx.stride());
    if (c != 0)
        mpn_add_n(x, x, ctx.modulus_limbs, ctx.stride());
}

void fp_sub_ui(const fp_ring_base& ctx, ulimb* x, const ulimb* a, ulimb b)
{
    ulimb N = ctx.stride();
    if (N == 1 && b >= ctx.modulus_limbs[0])
        b = b%ctx.modulus_limbs[0];
    auto c = mpn_sub_1(x, a, ctx.stride(), b);
    if (c != 0)
        mpn_add_n(x, x, ctx.modulus_limbs, ctx.stride());
}

void fp_mul(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    ulimb N = ctx.stride();
    my_mpn_mul_n(ctx.mul_tmp2N, a, b, N);
    my_mpn_tdiv_qr(ctx.mul_tmpNp2, x, ctx.mul_tmp2N, 2*N, ctx.modulus_limbs, N);
}

void fp_mul_ui(fp_ring& ctx, ulimb* x, const ulimb* a, ulimb b)
{
    ulimb N = ctx.stride();
    ulimb* t = ctx.mul_tmp2N;
    t[N] = mpn_mul_1(t, a, N, b);
    my_mpn_tdiv_qr(ctx.mul_tmpNp2, x, ctx.mul_tmp2N, N+1, ctx.modulus_limbs, N);
}

void fp_dot_mul(fp_ring& ctx, const ulimb* a, const ulimb* b)
{
    ulimb N = ctx.stride();
    ctx.mul_tmp2Np1[2*N] = 0;
    my_mpn_mul_n(ctx.mul_tmp2Np1, a, b, N);
}

void fp_dot_madd(fp_ring& ctx, const ulimb* a, const ulimb* b)
{
    ulimb N = ctx.stride();
    my_mpn_mul_n(ctx.mul_tmp2N, a, b, N);
    ulimb* t = ctx.mul_tmp2Np1;
    t[2*N] += mpn_add_n(t, t, ctx.mul_tmp2N, 2*N);
}

void fp_dot_reduce(fp_ring& ctx, ulimb* x)
{
    ulimb N = ctx.stride();
    my_mpn_tdiv_qr(ctx.mul_tmpNp2, x, ctx.mul_tmp2Np1, 2*N+1, ctx.modulus_limbs, N);
}

void fp_inv(fp_ring& ctx, ulimb* x, const ulimb* a)
{
    FLINT_ASSERT(fp_is_canonical(ctx, a));
    FLINT_ASSERT(!fp_is_zero(ctx, a));

    ulimb N = ctx.stride();
    ulimb* g = ctx.mul_tmpNp2;
    ulimb* A = ctx.mul_tmp2N + N*0;
    ulimb* M = ctx.mul_tmp2N + N*1;
    ulimb* u = ctx.mul_tmp2Np1;
    mp_size_t usize;

    MPN_COPY(A, a, N);
    MPN_COPY(M, ctx.modulus_limbs, N);

    auto gn = mpn_gcdext(g, u, &usize, A, N, M, N);
    if (gn != 1 || g[0] != 1)
        throw with<fp_ring, fmpz>(ctx, fmpz(g, gn, 0));

    if (usize < 0)
    {
        usize = -usize;
        FLINT_ASSERT(usize <= N);
        mpn_sub(x, ctx.modulus_limbs, N, u, usize);
    }
    else
    {
        FLINT_ASSERT(usize <= N);
        MPN_COPY(x, u, usize);
        MPN_ZERO(x + usize, N - usize);
    }
}

void fp_divexact(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    fp_tmp_elem t(ctx);
    fp_inv(ctx, t.data(), b);
    fp_mul(ctx, x, a, t.data());
}

void fp_random(random_state& state, fp_ring& ctx, ulimb* x)
{
    ulimb N = ctx.stride();
    ulimb i = 0;
    ulimb* t = ctx.mul_tmpNp2;
    do {
        t[i] = state.get_limb();
    } while (++i < N+1);
    ulimb q[3];
    my_mpn_tdiv_qr(q, x, t, N+1, ctx.modulus_limbs, N);
}

void fp_vec_scalar_mul(fp_ring& ctx, ulimb* X, const ulimb* A, ulimb n, const ulimb* b)
{
    ulimb N = ctx.stride();
    for (ulimb i = 0; i < n; i++)
        fp_mul(ctx, X + N*i, A + N*i, b);
}


void fp_vec_add(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b, ulimb n) {
    ulimb N = ctx.stride();
    for (ulimb i = 0; i < n; i++)
        fp_add(ctx, x + N*i, a + N*i, b + N*i);
}

void fp_vec_sub(fp_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b, ulimb n) {
    ulimb N = ctx.stride();
    for (ulimb i = 0; i < n; i++)
        fp_sub(ctx, x + N*i, a + N*i, b + N*i);
}

void fp_vec_neg(fp_ring& ctx, ulimb* x, const ulimb* a, ulimb n) {
    ulimb N = ctx.stride();
    for (ulimb i = 0; i < n; i++)
        fp_neg(ctx, x + N*i, a + N*i);
}

