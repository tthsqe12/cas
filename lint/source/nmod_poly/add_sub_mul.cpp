#include "nmod_poly.h"
#include "generics.h"

void nmod_set_ui(nmod_ring& ctx, ulimb* x, ulimb a)
{
    ulimb N = ctx.stride;

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

void nmod_set_ui_array(nmod_ring& ctx, ulimb* x, const ulimb* a, ulimb an)
{
    ulimb N = ctx.stride;

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

void nmod_neg(const nmod_ring_base& ctx, ulimb* x, const ulimb* a)
{
    if (ui_vec_is_zero(a, ctx.stride))
        ui_vec_zero(x, ctx.stride);
    else
        mpn_sub_n(x, ctx.modulus_limbs, a, ctx.stride);
}

void nmod_set_fmpz(nmod_ring& ctx, ulimb* x, fmpzc a)
{
    if (a.is_small())
        nmod_set_ui(ctx, x, my_abs(a.small()));
    else
        nmod_set_ui_array(ctx, x, a.limbs(), a.length());

    if (a.is_negative())
        nmod_neg(ctx, x, x);
}

void nmod_add(const nmod_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    auto c = mpn_add_n(x, a, b, ctx.stride);
    if (c != 0 || mpn_cmp(x, ctx.modulus_limbs, ctx.stride) >= 0)
        mpn_sub_n(x, x, ctx.modulus_limbs, ctx.stride);
}

void nmod_sub(const nmod_ring_base& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    auto c = mpn_sub_n(x, a, b, ctx.stride);
    if (c != 0)
        mpn_add_n(x, x, ctx.modulus_limbs, ctx.stride);
}

void nmod_mul(nmod_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    ulimb N = ctx.stride;
    my_mpn_mul_n(ctx.mul_tmp2N, a, b, N);
    my_mpn_tdiv_qr(ctx.mul_tmpNp1, x, ctx.mul_tmp2N, 2*N, ctx.modulus_limbs, N);
}

void nmod_mul_ui(nmod_ring& ctx, ulimb* x, const ulimb* a, ulimb b)
{
    ulimb N = ctx.stride;
    ulimb* t = ctx.mul_tmp2N;
    t[N] = mpn_mul_1(t, a, N, b);
    my_mpn_tdiv_qr(ctx.mul_tmpNp1, x, ctx.mul_tmp2N, N+1, ctx.modulus_limbs, N);
}

void nmod_dot_mul(nmod_ring& ctx, const ulimb* a, const ulimb* b)
{
    ulimb N = ctx.stride;
    ctx.mul_tmp2Np1[2*N] = 0;
    my_mpn_mul_n(ctx.mul_tmp2Np1, a, b, N);
}

void nmod_dot_madd(nmod_ring& ctx, const ulimb* a, const ulimb* b)
{
    ulimb N = ctx.stride;
    my_mpn_mul_n(ctx.mul_tmp2N, a, b, N);
    ulimb* t = ctx.mul_tmp2Np1;
    t[2*N] += mpn_add_n(t, t, ctx.mul_tmp2N, 2*N);
}

void nmod_dot_reduce(nmod_ring& ctx, ulimb* x)
{
    ulimb N = ctx.stride;
    my_mpn_tdiv_qr(ctx.mul_tmpNp1, x, ctx.mul_tmp2Np1, 2*N, ctx.modulus_limbs, N);
}

void nmod_inv(nmod_ring& ctx, ulimb* x, const ulimb* a)
{
    ulimb N = ctx.stride;
    ulimb* g = ctx.mul_tmpNp1;
    ulimb* A = ctx.mul_tmp2N + N*0;
    ulimb* M = ctx.mul_tmp2N + N*1;
    ulimb* u = ctx.mul_tmp2Np1;
    mp_size_t usize;

    MPN_COPY(A, a, N);
    MPN_COPY(M, ctx.modulus_limbs, N);

    auto gn = mpn_gcdext(g, u, &usize, A, N, M, N);
    if (gn != 1 || g[0] != 1)
        throw with<nmod_ring, fmpz>(ctx, fmpz(g, gn, 0));

    if (usize < 0)
    {
        usize = -usize;
        FLINT_ASSERT_ALWAYS(usize <= N);
        mpn_sub(x, ctx.modulus_limbs, N, u, usize);
    }
    else
    {
        FLINT_ASSERT_ALWAYS(usize <= N);
        MPN_COPY(x, u, usize);
        MPN_ZERO(x + usize, N - usize);
    }
}

void nmod_divexact(nmod_ring& ctx, ulimb* x, const ulimb* a, const ulimb* b)
{
    ulimb* t = ctx.tmp8N;
    nmod_inv(ctx, t, b);
    nmod_mul(ctx, x, a, t);
}

void nmod_randtest(nmod_ring& ctx, ulimb* x, rand_state& state)
{
    ulimb N = ctx.stride;
    ulimb i = 0;
    ulimb* t = ctx.mul_tmpNp1;
    do {
        t[i] = state.get_limb();
    } while (++i < N+1);
    ulimb q[3];
    my_mpn_tdiv_qr(q, x, t, N+1, ctx.modulus_limbs, N);
}

void nmod_vec_scalar_mul(nmod_ring& ctx, ulimb* X, ulimb* A, ulimb* b, ulimb n)
{
    ulimb N = ctx.stride;
    for (ulimb i = 0; i < n; i++)
        nmod_mul(ctx, X + N*i, A + N*i, b);
}


void nmod_poly::normalize_length(ulimb n, ulimb N)
{
    ulimb* d = coeffs.data();

    while (n > 0 && ui_vec_is_zero(d + N*(n-1), N))
        n--;

    set_length(n);
}

bool nmod_poly_is_canonical(const nmod_ring_base& ctx, const nmod_poly& a)
{
    ulimb N = ctx.stride;
    const ulimb* ad = a.data();
    ulimb an = a.length();

    FLINT_ASSERT(an*N <= a.coeffs._alloc);

    for (ulimb i = an; i > 0; i--)
    {
        if (!nmod_is_canonical(ctx, ad + N*(i-1)))
            return false;
        if (i == an && nmod_is_zero(ctx, ad + N*(i-1)))
            return false;
    }
    return true;
}

void nmod_poly_randtest(nmod_ring& ctx, nmod_poly& x, rand_state& state, ulimb len)
{
    ulimb N = ctx.stride;

    if (len < 1)
    {
        x.set_length(0);
        return;
    }

    ulimb* xd = x.fit_alloc(len, N);
    ulimb i = 0;
    do {
        nmod_randtest(ctx, xd + i*N, state);
    } while (++i < len);

    if (ui_vec_is_zero(xd + (len-1)*N, N))
        (xd + (len-1)*N)[0] = 1;

    x.set_length(len);
    return;
}

void nmod_poly_one(const nmod_ring_base& ctx, nmod_poly& x)
{
    ulimb N = ctx.stride;
    ulimb* xd = x.fit_alloc_destroy(1, N);
    x.set_length(1);
    nmod_one(ctx, xd + N*0);
}

bool nmod_poly_equal(const nmod_ring_base& ctx, const nmod_poly& a, const nmod_poly& b)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();
    if (an != bn)
        return false;
    for (ulimb i = 0; i < N*an; i++)
        if (ad[i] != bd[i])
            return false;
    return true;
}

int nmod_poly_cmp(const nmod_ring_base& ctx, const nmod_poly& a, const nmod_poly& b)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();
    if (an != bn)
        return an > bn ? 1 : -1;
    for (ulimb i = 0; i < N*an; i++)
        if (ad[i] != bd[i])
            return ad[i] > bd[i] ? 1 : -1;
    return 0;
}

void nmod_poly_set(const nmod_ring_base& ctx, nmod_poly& x, const nmod_poly& a)
{
    if (&x == &a)
        return;
    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb* xd = x.fit_alloc(an, N);
    x.set_length(an);
    ui_vec_set(xd, a.data(), N*an);    
}

void nmod_poly_set_fmpz(nmod_ring& ctx, nmod_poly& x, fmpzc a)
{
    ulimb N = ctx.stride;
    ulimb* xd = x.fit_alloc(1, N);
    nmod_set_fmpz(ctx, xd, a);
    x.normalize_length(1, N);    
}

void nmod_poly_neg(const nmod_ring_base& ctx, nmod_poly& x, const nmod_poly& a)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb* xd = x.fit_alloc(an, N);
    x.set_length(an);
    ulimb* ad = a.data();
    for (ulimb i = 0; i < an; i++)
        nmod_neg(ctx, xd + N*i, ad + N*i);
}

void nmod_poly_gen(const nmod_ring_base& ctx, nmod_poly& x)
{
    ulimb N = ctx.stride;
    ulimb* xd = x.fit_alloc_destroy(2, N);
    x.set_length(2);
    ui_vec_zero(xd + N*0, N);
    nmod_one(ctx, xd + N*1);
}

void nmod_poly_add(const nmod_ring_base& ctx, nmod_poly& x, const nmod_poly& a, const nmod_poly& b)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    ulimb xn = std::max(an, bn);
    ulimb* xd = x.fit_alloc(xn, N);
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    xn = std::min(an, bn);
    for (ulimb i = 0; i < xn; i++)
        nmod_add(ctx, xd + N*i, ad + N*i, bd + N*i);

    if (an == bn)
    {
        x.normalize_length(an, N);
    }
    else if (an >= bn)
    {
        ui_vec_set_nz(xd + N*bn, ad + N*bn, N*(an - bn));
        x.set_length(an);
    }
    else
    {
        ui_vec_set(xd + N*an, bd + N*an, N*(bn - an));
        x.set_length(bn);
    }
}


void nmod_poly_sub(const nmod_ring_base& ctx, nmod_poly& x, const nmod_poly& a, const nmod_poly& b)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    ulimb xn = std::max(an, bn);
    ulimb* xd = x.fit_alloc(xn, N);
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    xn = std::min(an, bn);
    for (ulimb i = 0; i < xn; i++)
        nmod_sub(ctx, xd + N*i, ad + N*i, bd + N*i);

    if (an == bn)
    {
        x.normalize_length(an, N);
    }
    else if (an > bn)
    {
        ui_vec_set(xd + N*bn, ad + N*bn, N*(an - bn));
        x.set_length(an);
    }
    else
    {
        ulimb i = an;
        do {
            nmod_neg(ctx, xd + N*i, bd + N*i);
        } while (++i < bn);

        x.set_length(bn);
    }
}

void nmod_poly_make_monic_with_lcinv(nmod_ring& ctx,
    nmod_poly& x,
    const nmod_poly& a, const ulimb* alcinv)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    const ulimb* ad = a.data();

    FLINT_ASSERT(an > 0);

    if (nmod_is_one(ctx, alcinv)) {
        nmod_poly_set(ctx, x, a);
        return;
    }

    ulimb* xd = x.fit_alloc(an, N);
    x.set_length(an);
    nmod_one(ctx, xd + N*(an-1));
    for (ulimb i = 0; i < an-1; i++)
        nmod_mul(ctx, xd + N*i, ad + N*i, alcinv);
}


void nmod_poly_make_monic(nmod_ring& ctx, nmod_poly& x, const nmod_poly& a)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    const ulimb* ad = a.data();

    if (an < 1) {
        x.set_length(0);
        return;
    }

    if (nmod_is_one(ctx, ad + N*(an-1))) {
        nmod_poly_set(ctx, x, a);
        return;
    }

    if (an == 1) {
        nmod_poly_one(ctx, x);
        return;
    }

    nmod_inv(ctx, ctx.tmp8N + N*0, ad + N*(an-1));
    nmod_poly_make_monic_with_lcinv(ctx, x, a, ctx.tmp8N + N*0);
}

void nmod_poly_mul(nmod_ring& ctx, nmod_poly& x, const nmod_poly& a, const nmod_poly& b)
{
    FLINT_ASSERT(&x != &a); FLINT_ASSERT(&x != &b);

    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    if (an < 1 || bn < 1)
    {
        x.set_length(0);
        return;
    }

    ulimb xn = an + bn - 1;
    ulimb* xd = x.fit_alloc(xn, N);

    for (ulimb i = 0; i < xn; i++)
    {
        ulong jstart = i > bn - 1 ? i - (bn - 1) : 0;
        ulong jstop = std::min(i + 1, an);
        FLINT_ASSERT(jstart < jstop);

        ulimb j = jstart;
        nmod_dot_mul(ctx, ad + N*j, bd + (i - j)*N);
        for (j++; j < jstop; j++)
            nmod_dot_madd(ctx, ad + N*j, bd + (i - j)*N);

        nmod_dot_reduce(ctx, xd + N*i);
    }

    x.normalize_length(xn, N);
}

void nmod_poly_sqr(nmod_ring& ctx, nmod_poly& x, const nmod_poly& a)
{
    nmod_poly_mul(ctx, x, a, a);
}

void nmod_poly_pow_ui(nmod_ring& R, nmod_poly& x, const nmod_poly& a, ulimb n)
{
    generic_pow_ui_binexp<nmod_ring, nmod_poly>(R, x, a, n);
}


void nmod_poly_divexact(nmod_ring& ctx,
    nmod_poly& q,
    const nmod_poly& a,
    const nmod_poly& b)
{
    nmod_poly r;
    nmod_poly_divrem(ctx, q, r, a, b);
    FLINT_ASSERT_ALWAYS(nmod_poly_is_zero(ctx, r));
}

bool nmod_poly_divides(nmod_ring& ctx,
    nmod_poly& q,
    const nmod_poly& a,
    const nmod_poly& b)
{
    nmod_poly r;
    nmod_poly_divrem(ctx, q, r, a, b);
    return nmod_poly_is_zero(ctx, r);
}

void nmod_poly_divrem(nmod_ring& ctx,
    nmod_poly& q,
    nmod_poly& r,
    const nmod_poly& a,
    const nmod_poly& b)
{
    FLINT_ASSERT(&q != &a); FLINT_ASSERT(&q != &b);
    FLINT_ASSERT(&r != &a); FLINT_ASSERT(&r != &b);
    FLINT_ASSERT(nmod_poly_is_canonical(ctx, a));
    FLINT_ASSERT(nmod_poly_is_canonical(ctx, b));

    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* bd = b.data();

    if (UNLIKELY(an < bn || bn < 1))
    {
        nmod_poly_set(ctx, r, a);
        q.set_length(0);
        return;
    }

    ulimb* blcinv = ctx.tmp8N + N*0;
    nmod_inv(ctx, blcinv, bd + N*(bn - 1));

    q.set_length(an - (bn - 1));
    ulimb* qd = q.fit_alloc(an - (bn - 1), N);

    if (UNLIKELY(bn == 1))
    {
        nmod_inv(ctx, blcinv, bd + N*(bn - 1));
        for (ulimb i = 0; i < an; i++)
            nmod_mul(ctx, qd + i*N, a.data() + i*N, blcinv);
        r.set_length(0);
        return;
    }

    ulimb* rd = r.fit_alloc(an, N);
    ui_vec_set(rd, a.data(), an*N);

    ulimb* t = ctx.tmp8N + N*1;
    for (ulimb i = an - 1; i >= bn - 1; i--)
    {
        nmod_mul(ctx, qd + (i - (bn - 1))*N, rd + (i)*N, blcinv);

        for (ulimb j = 0; j < bn - 1; j++)
        {
            nmod_mul(ctx, t, bd + j*N, qd + (i - (bn - 1))*N);
            nmod_sub(ctx, rd + (j+i-(bn-1))*N,
                          rd + (j+i-(bn-1))*N, t);
        }
    }

    r.normalize_length(bn - 1, N);
}

static slimb _nmod_poly_gcd_euclidean_inplace(nmod_ring& ctx,
    ulimb* A, ulimb Alen,
    ulimb* B, ulimb Blen)
{
    ulimb N = ctx.stride;
    ulimb i;
    ulimb* u  = ctx.tmp8N + N*0;
    ulimb* q0 = ctx.tmp8N + N*1;
    ulimb* q1 = ctx.tmp8N + N*2;
    ulimb s = 0;

again:

    if (Alen < Blen)
    {
        std::swap(A, B);
        std::swap(Alen, Blen);
        s =~ s;
    }

    FLINT_ASSERT(Alen >= Blen);

    if (UNLIKELY(Blen < 2))
    {
        if (Blen < 1)
        {
            nmod_inv(ctx, u, A + N*(Alen - 1));
            for (i = 0; i < Alen - 1; i++)
                nmod_mul(ctx, A + N*i, A + N*i, u);
            nmod_one(ctx, A + N*(Alen - 1));
            return slimb(Alen)^s;
        }
        else
        {
            nmod_one(ctx, B + N*0);
            return slimb(Blen)^~s;
        }
    }

    if (Alen == Blen)
    {
        // generate one quotient
        nmod_inv(ctx, u, B + N*(Blen - 1));
        nmod_mul(ctx, q0, A + N*(Alen - 1), u);

        for (i = 0; i < Blen - 1; i++)
        {
            nmod_mul(ctx, u, q0, B + N*i);
            nmod_sub(ctx, A + N*i, A + N*i, u);
        }

        Alen -= 1;
    }
    else
    {
        // generate two quotients
        nmod_inv(ctx, u, B + N*(Blen - 1));
        nmod_mul(ctx, q1, A + N*(Alen - 1), u);
        nmod_mul(ctx, q0, q1, B + N*(Blen - 2));
        nmod_sub(ctx, q0, q0, A + N*(Alen - 2));
        nmod_mul(ctx, q0, q0, u);

        nmod_neg(ctx, q1, q1);

        nmod_mul(ctx, u, q0, B + N*0);
        nmod_add(ctx, A + N*(-1 + Alen - Blen),
                      A + N*(-1 + Alen - Blen), u);

        for (i = 0; i < Blen - 1; i++)
        {
            nmod_dot_mul(ctx, q1, B + N*i);
            nmod_dot_madd(ctx, q0, B + N*(i + 1));
            nmod_dot_reduce(ctx, u);
            nmod_add(ctx, A + N*(i + Alen - Blen),
                          A + N*(i + Alen - Blen), u);
        }

        Alen -= 2;
    }

    while (Alen > 0 && nmod_is_zero(ctx, A + N*(Alen - 1)))
        Alen--;

    goto again;
}

void nmod_poly_gcd(nmod_ring& ctx,
    nmod_poly& g,
    const nmod_poly& a,
    const nmod_poly& b)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    tmp_allocator push;

    if (an < 1) {
        nmod_poly_make_monic(ctx, g, b);
        return;
    }
    if (bn < 1) {
        nmod_poly_make_monic(ctx, g, a);
        return;
    }

    ulimb* A = push.recursive_alloc<ulimb>(N*(an + bn));
    ulimb* B = A + N*an;

    ui_vec_set(A, a.data(), an*N);
    ui_vec_set(B, b.data(), bn*N);

    slimb gn = _nmod_poly_gcd_euclidean_inplace(ctx, A, an, B, bn);

    if (gn < 0) {
        gn = ~gn;
        std::swap(A, B);
    }

    ulimb* gd = g.fit_alloc_destroy(gn, N);
    g.set_length(gn);    
    ui_vec_set(gd, A, gn);
}

void nmod_poly_gcdc(
    nmod_ring& ctx,
    nmod_poly& g, nmod_poly& abar, nmod_poly& bbar,
    const nmod_poly& a, const nmod_poly& b)
{
    nmod_poly t0, t1, t2;
    nmod_poly_gcd(ctx, t0, a, b);
    nmod_poly_divexact(ctx, t1, a, t0);
    nmod_poly_divexact(ctx, t2, b, t0);
    std::swap(g, t0);
    std::swap(abar, t1);
    std::swap(bbar, t2);
}

/*
Want to find m:    [m11 m12] [a] = [g]
                   [m21 m22] [b]   [0]
m starts out as [1 0; 0 1].
If [A; B] is multiplied by M on the left, then so is m
*/
void nmod_poly_gcdx(nmod_ring& ctx,
    nmod_poly& g, nmod_poly& m11, nmod_poly& m12,
    const nmod_poly& a_,
    const nmod_poly& b_)
{
    ulimb N = ctx.stride;
    nmod_poly q, r, a, b, m21, m22;

    nmod_poly_set(ctx, a, a_);
    nmod_poly_set(ctx, b, b_);

    nmod_poly_one(ctx, m11);
    nmod_poly_zero(ctx, m12);
    nmod_poly_zero(ctx, m21);
    nmod_poly_one(ctx, m22);

    if (b.length() == 0)
    {
        if (a.length() == 0)
        {
            nmod_poly_zero(ctx, g);
            return;
        }

        nmod_set(ctx, m22.data() + N*0, a.data() + N*(a.length()-1));
        nmod_inv(ctx, m11.data() + N*0, a.data() + N*(a.length()-1));
        nmod_poly_make_monic_with_lcinv(ctx, g, a, m11.data() + N*0);
        return;
    }

    if (a.length() < b.length())
    {
        std::swap(a, b);
        std::swap(m11, m21);
        std::swap(m12, m22);
    }

    while (b.length() > 0)
    {
        nmod_poly_divrem(ctx, q, r, a, b);
        std::swap(a, r);
        std::swap(a, b);

        nmod_poly_mul(ctx, r, m21, q);
        nmod_poly_sub(ctx, m11, m11, r);
        std::swap(m11, m21);
        FLINT_ASSERT(m11.length() <= m21.length());

        nmod_poly_mul(ctx, r, m22, q);
        nmod_poly_sub(ctx, m12, m12, r);
        std::swap(m12, m22);
        FLINT_ASSERT(m12.length() <= m22.length());
    }

    FLINT_ASSERT(a.length() > 0);

    ulimb* lcinv = ctx.tmp8N + N*0;
    nmod_inv(ctx, lcinv, a.data() + N*(a.length()-1));
    nmod_vec_scalar_mul(ctx, m11.data(), m11.data(), lcinv, m11.length());
    nmod_vec_scalar_mul(ctx, m12.data(), m12.data(), lcinv, m12.length());
    nmod_vec_scalar_mul(ctx, m21.data(), m21.data(), a.data() + N*(a.length()-1), m21.length());
    nmod_vec_scalar_mul(ctx, m22.data(), m22.data(), a.data() + N*(a.length()-1), m22.length());
    nmod_poly_make_monic_with_lcinv(ctx, g, a, lcinv);
    return;
}

void nmod_poly_deflate_inplace(
    const nmod_ring_base& ctx,
    nmod_poly& x,
    ulimb e)
{
    FLINT_ASSERT(e > 1);
    ulimb N = ctx.stride;
    ulimb xn = x.length();
    ulimb* xd = x.data();

    if (xn == 0)
        return;

    ulimb i;
    for (i = 1; e*i < xn; i++)
        ui_vec_set(xd + N*i, xd + N*e*i, N);

    x.set_length(i);
}

void nmod_poly_derivative(
    nmod_ring& ctx,
    nmod_poly& x,
    const nmod_poly& a)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();

    if (an <= 1)
    {
        x.set_length(0);
        return;
    }

    ulimb* xd = x.fit_alloc(an - 1, N);
    const ulimb* ad = a.data();

    for (ulimb i = 1; i < an; i++)
        nmod_mul_ui(ctx, xd + N*(i-1), ad + N*i, i);

    x.normalize_length(an - 1, N);
}


/*
static void _append_factor(
    nmod_ring& ctx,
    nmod_poly_product& f,
    const nmod_poly& a, ulimb e)
{
    if (a.is_constant())
        return;

    ulimb n = f.length();
    f.fit_alloc(n+1);
    f.set_length(n+1);
    f.exp(n) = e;
    nmod_poly_make_monic(ctx, f.base(n), a);
}

static void _append_factor_pow_normalize(
    nmod_ring& ctx,
    nmod_poly_product& f,
    nmod_poly& a, ulimb e,
    ulimb i)
{
    FLINT_ASSERT(e > 0);
    FLINT_ASSERT(a.length() > 0);

    nmod_poly g, lbar, abar;

    while (i < f.length() && a.length() > 0)
    {
        // (g*lbar)^l[i].exp * (g*abar)^e
        // (lbar)^l[i].exp * (g)^(l[i].exp+e) * (abar)^e
        nmod_poly_gcd_cofactors(ctx, g, lbar, abar, f.base(i), a);
        if (g.is_constant())
        {
            i++;
        }
        else if (lbar.is_constant())
        {
            std::swap(a, abar);
            std::swap(f.base(i), g);
            f.exp(i) += e;
        }
        else if (abar.is_constant())
        {
            std::swap(a, g);
            e += f.exp(i);
            nmod_poly_make_monic(ctx, f.base(i), lbar);
        }
        else
        {
            std::swap(a, abar);
            _append_factor_pow_normalize(ctx, f, g, e, i);
        }   
    }

    _append_factor(ctx, f, a, e);
}

void nmod_poly_factor_squarefree(
    nmod_ring& ctx,
    nmod_poly_product& f,
    const nmod_poly& A)
{
    f.set_length(0);

    if (A.length() < 2)
        return;

    if (A.length() == 2)
    {
        f.fit_alloc_destroy(1);
        nmod_poly_make_monic(ctx, f.base(0), A);
        f.exp(0) = 1;
        f.set_length(1);
        return;
    }

    nmod_poly B, C, U, V, W, G;

    nmod_poly_set(ctx, C, A);
    nmod_poly_derivative(ctx, G, C);
    nmod_poly_gcd_cofactors(ctx, C, W, V, C, G);

    ulimb k, p = fmpz_abs_fits_ui(ctx.modulus) ? fmpz_get_ui(ctx.modulus) : 1;
    for (k = 1; k <= p-2 && !(nmod_poly_derivative(ctx, G, W),
                             nmod_poly_sub(ctx, U, V, G),
                             nmod_poly_is_zero(ctx, U)); k++)
    {
        nmod_poly_gcd_cofactors(ctx, G, W, V, W, U);
        _append_factor(ctx, f, G, k);
        nmod_poly_divexact(ctx, U, C, W);
        std::swap(C, U);
    }
    _append_factor(ctx, f, W, k);

    FLINT_ASSERT(!nmod_poly_is_zero(ctx, C));
    FLINT_ASSERT((nmod_poly_derivative(ctx, U, C), nmod_poly_is_zero(ctx, U)));
    
    if (C.length() > 1)
    {
        nmod_poly_product Tf;
        nmod_poly_deflate_inplace(ctx, C, p);
        nmod_poly_factor_squarefree(ctx, Tf, C);
        for (ulimb i = Tf.length(); i > 0; i--)
            _append_factor_pow_normalize(ctx, f, Tf.base(i-1), p*Tf.exp(i-1), 0);
    }
}
*/

void nmod_poly_factor_squarefree(
    nmod_ring& ctx,
    nmod_poly_product& f,
    const nmod_poly& a)
{
    f.set_length(0);
    nmod_poly acopy(a);
    generic_finite_field_poly_factor_squarefree<nmod_ring, nmod_poly>(ctx, f, acopy, 1, -ulimb(1));
}
