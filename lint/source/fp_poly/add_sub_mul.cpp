#include "fmpz.h"
#include "fp_poly.h"
#include "generic/generics.h"


void fp_poly::normalize_length(ulimb n, ulimb N)
{
    ulimb* d = data();
    while (n > 0 && ui_vec_is_zero(d + N*(n-1), N))
        n--;
    set_length(n);
}

bool fp_poly_is_canonical(const fp_ring_base& ctx, const fp_poly& a)
{
    ulimb N = ctx.stride();
    const ulimb* ad = a.data();
    ulimb an = a.length();

    if (an*N > a._alloc)
    {
        std::cerr << "poly length exceeds alloc" << std::endl;
        return false;
    }

    for (ulimb i = an; i > 0; i--)
    {
        if (!fp_is_canonical(ctx, ad + N*(i-1)))
        {
            std::cerr << "poly coefficient is not canonical" << std::endl;
            return false;
        }
        if (i == an && fp_is_zero(ctx, ad + N*(i-1)))
        {
            std::cerr << "leading coefficient is zero" << std::endl;
            return false;
        }
    }
    return true;
}


void fp_poly_one(const fp_ring_base& ctx, fp_poly& x)
{
    ulimb N = ctx.stride();
    ulimb* xd = x.fit_alloc_destroy(1, N);
    x.set_length(1);
    fp_one(ctx, xd + N*0);
}

bool fp_poly_is_one(const fp_ring_base& ctx, const fp_poly& x)
{
    return x.length() == 1 && fp_is_one(ctx, x.data() + ctx.stride()*0);
}

bool fp_poly_equal_si(fp_ring& ctx, const fp_poly& a, slimb b)
{
    if (a.length() > 1)
        return false;

    if (a.length() == 1)
        return fp_equal_si(ctx, a.data(), b);

    fmpz bb(b);
    return fmpz_divisible(bb, ctx.characteristic());
}

bool fp_poly_equal(const fp_ring_base& ctx, const fp_poly& a, const fp_poly& b)
{
    ulimb N = ctx.stride();
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

int fp_poly_cmp(const fp_ring_base& ctx, const fp_poly& a, const fp_poly& b)
{
    ulimb N = ctx.stride();
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

void fp_poly_set(const fp_ring_base& ctx, fp_poly& x, const fp_poly& a)
{
    if (&x == &a)
        return;
    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb* xd = x.fit_alloc(an, N);
    x.set_length(an);
    ui_vec_set(xd, a.data(), N*an);    
}

void fp_poly_set_fmpz(fp_ring& ctx, fp_poly& x, fmpzc a)
{
    ulimb N = ctx.stride();
    ulimb* xd = x.fit_alloc(1, N);
    fp_set_fmpz(ctx, xd, a);
    x.normalize_length(1, N);    
}

void fp_poly_neg(const fp_ring_base& ctx, fp_poly& x, const fp_poly& a)
{
    FLINT_ASSERT(fp_poly_is_canonical(ctx, a));

    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb* xd = x.fit_alloc(an, N);
    x.set_length(an);
    const ulimb* ad = a.data();
    for (ulimb i = 0; i < an; i++)
        fp_neg(ctx, xd + N*i, ad + N*i);
}

void fp_poly_gen(const fp_ring_base& ctx, fp_poly& x)
{
    ulimb N = ctx.stride();
    ulimb* xd = x.fit_alloc_destroy(2, N);
    x.set_length(2);
    ui_vec_zero(xd + N*0, N);
    fp_one(ctx, xd + N*1);
}

void fp_poly_add(const fp_ring_base& ctx, fp_poly& x, const fp_poly& a, const fp_poly& b)
{
    FLINT_ASSERT(fp_poly_is_canonical(ctx, a));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, b));

    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();
    ulimb xn = std::max(an, bn);
    ulimb* xd = x.fit_alloc(xn, N);
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    xn = std::min(an, bn);
    for (ulimb i = 0; i < xn; i++)
        fp_add(ctx, xd + N*i, ad + N*i, bd + N*i);

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


void fp_poly_sub(const fp_ring_base& ctx, fp_poly& x, const fp_poly& a, const fp_poly& b)
{
    FLINT_ASSERT(fp_poly_is_canonical(ctx, a));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, b));

    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();
    ulimb xn = std::max(an, bn);
    ulimb* xd = x.fit_alloc(xn, N);
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    xn = std::min(an, bn);
    for (ulimb i = 0; i < xn; i++)
        fp_sub(ctx, xd + N*i, ad + N*i, bd + N*i);

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
            fp_neg(ctx, xd + N*i, bd + N*i);
        } while (++i < bn);

        x.set_length(bn);
    }
}

void fp_poly_make_monic_with_lcinv(fp_ring& ctx,
    fp_poly& x,
    const fp_poly& a, const ulimb* alcinv)
{
    ulimb N = ctx.stride();
    ulimb an = a.length();
    const ulimb* ad = a.data();

    FLINT_ASSERT(an > 0);

    if (fp_is_one(ctx, alcinv)) {
        fp_poly_set(ctx, x, a);
        return;
    }

    ulimb* xd = x.fit_alloc(an, N);
    x.set_length(an);
    fp_one(ctx, xd + N*(an-1));
    for (ulimb i = 0; i < an-1; i++)
        fp_mul(ctx, xd + N*i, ad + N*i, alcinv);
}


void fp_poly_make_monic(fp_ring& ctx, fp_poly& x, const fp_poly& a)
{
    ulimb N = ctx.stride();
    ulimb an = a.length();
    const ulimb* ad = a.data();

    if (an < 1) {
        x.set_length(0);
        return;
    }

    if (fp_is_one(ctx, ad + N*(an-1))) {
        fp_poly_set(ctx, x, a);
        return;
    }

    if (an == 1) {
        fp_poly_one(ctx, x);
        return;
    }

    fp_tmp_elem t(ctx);
    fp_inv(ctx, t.data(), ad + N*(an-1));
    fp_poly_make_monic_with_lcinv(ctx, x, a, t.data());
}

bool fp_poly_is_monic(fp_ring& ctx, const fp_poly& a)
{
    return a.length() > 0 && fp_is_one(ctx, a.data() + ctx.stride()*(a.length()-1));
}

void fp_poly_mul(fp_ring& ctx, fp_poly& x, _fp_poly a, _fp_poly b)
{
//    FLINT_ASSERT(&x != &a); FLINT_ASSERT(&x != &b);
    FLINT_ASSERT(_fp_poly_is_canonical(ctx, a));
    FLINT_ASSERT(_fp_poly_is_canonical(ctx, b));

    ulimb N = ctx.stride();
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
        fp_dot_mul(ctx, ad + N*j, bd + (i - j)*N);
        for (j++; j < jstop; j++)
            fp_dot_madd(ctx, ad + N*j, bd + (i - j)*N);

        fp_dot_reduce(ctx, xd + N*i);
    }

    x.normalize_length(xn, N);
    FLINT_ASSERT(fp_poly_is_canonical(ctx, x));
}

void fp_poly_submul(fp_ring& ctx, fp_poly& x, const fp_poly& a, const fp_poly& b)
{
    FLINT_ASSERT(&x != &a); FLINT_ASSERT(&x != &b);
    FLINT_ASSERT(fp_poly_is_canonical(ctx, x));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, a));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, b));

    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();
    ulimb xn = x.length();
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    if (an < 1 || bn < 1)
        return;

    ulimb abn = an + bn - 1;
    ulimb* xd = x.fit_alloc(abn /*std::max(xn, abn)*/, N);

    fp_tmp_elem t(ctx);
    fp_dotter dotter(ctx);
    for (ulimb i = 0; i < abn; i++)
    {
        ulong jstart = i > bn - 1 ? i - (bn - 1) : 0;
        ulong jstop = std::min(i + 1, an);
        FLINT_ASSERT(jstart < jstop);

        ulimb j = jstart;
        dotter.mul(ctx, ad + N*j, bd + (i - j)*N);
        for (j++; j < jstop; j++)
            dotter.madd(ctx, ad + N*j, bd + (i - j)*N);

        dotter.reduce(ctx, t.data());
        if (i < xn)
            fp_sub(ctx, xd + N*i, xd + N*i, t.data());
        else
            fp_neg(ctx, xd + N*i, t.data());
    }

    x.normalize_length(std::max(xn, abn), N);

    FLINT_ASSERT(fp_poly_is_canonical(ctx, x));
}


void fp_poly_sqr(fp_ring& ctx, fp_poly& x, const fp_poly& a)
{
    fp_poly_mul(ctx, x, a, a);
}

void fp_poly_pow_ui(fp_ring& R, fp_poly& x, const fp_poly& a, ulimb n)
{
    generic_pow_binexp<fp_ring, fp_poly>(R, x, a, n);
}


void fp_poly_divexact(fp_ring& ctx,
    fp_poly& q,
    const fp_poly& a,
    const fp_poly& b)
{
    fp_poly r;
    if (&q == &a || &q == &b)
    {
        fp_poly t;
        fp_poly_divrem(ctx, t, r, a, b);
        std::swap(q, t);
    }
    else
        fp_poly_divrem(ctx, q, r, a, b);
    FLINT_ASSERT_ALWAYS(fp_poly_is_zero(ctx, r));
}

bool fp_poly_divides(fp_ring& ctx,
    fp_poly& q,
    const fp_poly& a,
    const fp_poly& b)
{
    if (a.length() < b.length())
    {
        if (a.length() == 0)
            return true;

        fp_tmp_elem t(ctx);
        /* only can say no if leading coeff is invertible */
        fp_inv(ctx, t.data(), b.data() + ctx.stride()*b.degree());
        return false;
    }

    fp_poly r;
    fp_poly_divrem(ctx, q, r, a, b);
    return fp_poly_is_zero(ctx, r);
}

void fp_poly_divrem(fp_ring& ctx,
    fp_poly& q,
    fp_poly& r,
    const fp_poly& a,
    const fp_poly& b)
{
    FLINT_ASSERT(&q != &r);
    FLINT_ASSERT(&q != &a); FLINT_ASSERT(&q != &b);
    FLINT_ASSERT(&r != &a); FLINT_ASSERT(&r != &b);
    FLINT_ASSERT(fp_poly_is_canonical(ctx, a));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, b));

    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* bd = b.data();

    if (UNLIKELY(an < bn || bn < 1))
    {
        fp_poly_set(ctx, r, a);
        q.set_length(0);
        return;
    }

    fp_tmp_elem blcinv(ctx);
    fp_inv(ctx, blcinv.data(), bd + N*(bn - 1));

    q.set_length(an - (bn - 1));
    ulimb* qd = q.fit_alloc(an - (bn - 1), N);

    if (UNLIKELY(bn == 1))
    {
        fp_inv(ctx, blcinv.data(), bd + N*(bn - 1));
        for (ulimb i = 0; i < an; i++)
            fp_mul(ctx, qd + i*N, a.data() + i*N, blcinv.data());
        r.set_length(0);
        return;
    }

    ulimb* rd = r.fit_alloc(an, N);
    ui_vec_set(rd, a.data(), an*N);

    fp_tmp_elem t(ctx);
    for (ulimb i = an - 1; i >= bn - 1; i--)
    {
        fp_mul(ctx, qd + (i - (bn - 1))*N, rd + (i)*N, blcinv.data());

        for (ulimb j = 0; j < bn - 1; j++)
        {
            fp_mul(ctx, t.data(), bd + j*N, qd + (i - (bn - 1))*N);
            fp_sub(ctx, rd + (j+i-(bn-1))*N,
                        rd + (j+i-(bn-1))*N, t.data());
        }
    }

    r.normalize_length(bn - 1, N);

    FLINT_ASSERT(fp_poly_is_canonical(ctx, q));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, r));
}

void fp_poly_divrem_inplace(fp_ring& ctx,
    fp_poly& q,
    fp_poly& a,
    const fp_poly& b)
{
    FLINT_ASSERT(&q != &a); FLINT_ASSERT(&q != &b); FLINT_ASSERT(&a != &b);
    FLINT_ASSERT(fp_poly_is_canonical(ctx, a));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, b));

    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* bd = b.data();

    if (UNLIKELY(an < bn || bn < 1))
    {
        q.set_length(0);
        return;
    }

    fp_tmp_elem blcinv(ctx);
    fp_inv(ctx, blcinv.data(), bd + N*(bn - 1));

    q.set_length(an - (bn - 1));
    ulimb* qd = q.fit_alloc(an - (bn - 1), N);

    if (UNLIKELY(bn == 1))
    {
        fp_inv(ctx, blcinv.data(), bd + N*(bn - 1));
        for (ulimb i = 0; i < an; i++)
            fp_mul(ctx, qd + i*N, a.data() + i*N, blcinv.data());
        a.set_length(0);
        return;
    }

    ulimb* rd = a.data();

    fp_tmp_elem t(ctx);
    for (ulimb i = an - 1; i >= bn - 1; i--)
    {
        fp_mul(ctx, qd + (i - (bn - 1))*N, rd + (i)*N, blcinv.data());

        for (ulimb j = 0; j < bn - 1; j++)
        {
            fp_mul(ctx, t.data(), bd + j*N, qd + (i - (bn - 1))*N);
            fp_sub(ctx, rd + (j+i-(bn-1))*N,
                        rd + (j+i-(bn-1))*N, t.data());
        }
    }

    a.normalize_length(bn - 1, N);

    FLINT_ASSERT(fp_poly_is_canonical(ctx, q));
    FLINT_ASSERT(fp_poly_is_canonical(ctx, a));
}


void fp_poly_mod(fp_ring& ctx,
    fp_poly& r,
    const fp_poly& a,
    const fp_poly& b)
{
    fp_poly q, t;
    fp_poly_divrem(ctx, q, t, a, b);
    std::swap(r, t);
}

void fp_poly_rem(fp_ring& ctx,
    fp_poly& r,
    const fp_poly& a,
    const fp_poly& b)
{
    fp_poly q;
    if (&r == &a || &r == &b)
    {
        fp_poly t;
        fp_poly_divrem(ctx, q, t, a, b);
        std::swap(r, t);
    }
    else
        fp_poly_divrem(ctx, q, r, a, b);
    FLINT_ASSERT_ALWAYS(fp_poly_is_zero(ctx, r));
}

void fp_poly_mulmod(fp_ring& ctx,
    fp_poly& x,
    const fp_poly& a,
    const fp_poly& b,
    const fp_poly& m)
{
    fp_poly q, t;
    fp_poly_mul(ctx, t, a, b);
    fp_poly_divrem(ctx, q, x, t, m);
}

void fp_poly_sqrmod(fp_ring& ctx,
    fp_poly& x,
    const fp_poly& a,
    const fp_poly& m)
{
    fp_poly q, t;
    fp_poly_mul(ctx, t, a, a);
    fp_poly_divrem(ctx, q, x, t, m);
}

void fp_poly_powmod_fmpz(fp_ring& ctx,
    fp_poly& x,
    const fp_poly& a,
    fmpzc e,
    const fp_poly& m)
{
    generic_powmod_binexp<fp_ring, fp_poly>(ctx, x, a, e, m);
}



void fp_poly_deflate_inplace(
    const fp_ring_base& ctx,
    fp_poly& x,
    ulimb e)
{
    FLINT_ASSERT(e > 1);
    ulimb N = ctx.stride();
    ulimb xn = x.length();
    ulimb* xd = x.data();

    if (xn == 0)
        return;

    ulimb i;
    for (i = 1; e*i < xn; i++)
        ui_vec_set(xd + N*i, xd + N*e*i, N);

    x.set_length(i);
}

void fp_poly_derivative(
    fp_ring& ctx,
    fp_poly& x,
    const fp_poly& a)
{
    ulimb N = ctx.stride();
    ulimb an = a.length();

    if (an <= 1)
    {
        x.set_length(0);
        return;
    }

    ulimb* xd = x.fit_alloc(an - 1, N);
    const ulimb* ad = a.data();

    for (ulimb i = 1; i < an; i++)
        fp_mul_ui(ctx, xd + N*(i-1), ad + N*i, i);

    x.normalize_length(an - 1, N);
}

