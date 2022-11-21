#include "fp_poly.h"
#include "generic/generics.h"
#include "generic/generic_poly_gcd.h"

bool _fp_poly_is_canonical(fp_ring& ctx, _fp_poly& a)
{
    return _generic_poly_is_canonical<fp_ring>(ctx, a);
}

void fp_poly_gcd(fp_ring& ctx, fp_poly& g, const fp_poly& a, const fp_poly& b)
{
    return generic_poly_gcd<fp_ring>(ctx, g, a, b);
}

int fp_poly_hgcd(fp_ring& ctx,
    fp_poly& m11, fp_poly& m12, fp_poly& m21, fp_poly& m22,
    fp_poly& A, fp_poly& B,
    const fp_poly& a, const fp_poly& b)
{
    return generic_poly_hgcd<fp_ring>(ctx, m11, m12, m21, m22, A, B, a, b);
}

void fp_poly_gcdc(
    fp_ring& ctx,
    fp_poly& g, fp_poly& abar, fp_poly& bbar,
    const fp_poly& a, const fp_poly& b)
{
    fp_poly t0, t1, t2;
    fp_poly_gcd(ctx, t0, a, b);
    fp_poly_divexact(ctx, t1, a, t0);
    fp_poly_divexact(ctx, t2, b, t0);
    std::swap(g, t0);
    std::swap(abar, t1);
    std::swap(bbar, t2);
}

void fp_poly_gcdx(fp_ring& ctx,
    fp_poly& g, fp_poly& m11, fp_poly& m12,
    const fp_poly& a,
    const fp_poly& b)
{
    return generic_poly_gcdx<fp_ring>(ctx, g, m11, m12, a, b);
#if 0
    ulimb N = ctx.stride();
    fp_poly q, r, a, b, m21, m22;

    fp_poly_set(ctx, a, a_);
    fp_poly_set(ctx, b, b_);

    fp_poly_one(ctx, m11);
    fp_poly_zero(ctx, m12);
    fp_poly_zero(ctx, m21);
    fp_poly_one(ctx, m22);

    if (b.length() == 0)
    {
        if (a.length() == 0)
        {
            zero(ctx, g);
            return;
        }

        set(ctx, m22.data() + N*0, a.data() + N*(a.length()-1));
        inv(ctx, m11.data() + N*0, a.data() + N*(a.length()-1));
        fp_poly_make_monic_with_lcinv(ctx, g, a, m11.data() + N*0);
        return;
    }

    if (a.length() < b.length())
    {
        std::swap(a, b);
        std::swap(m11, m21);
        std::swap(m12, m22);
    }
    else if (a.length() == b.length())
    {
            // [b      ] = [0  1] [a]
            // [a - q b]   [1 -q] [b]
            fp_poly_divrem_inplace(ctx, q, a, b); std::swap(a, b);
            fp_poly_submul(ctx, m11, q, m21);     std::swap(m11, m21);
            fp_poly_submul(ctx, m12, q, m22);     std::swap(m12, m22);
    }

    while (a.length() >= POLY_HGCD_CUTOFF)
    {
FLINT_ASSERT(a.length() > b.length());
        fp_poly A, B, n11, n12, n21, n22, t11, t12, t21, t22;
        int ndet = fp_poly_hgcd(ctx, n11, n12, n21, n22, A, B, a, b);

        // multiply m by n^-1 on the left
        // [n22 -n12 ] [m11 m12]
        // [-n21  n11] [m21 m22]

        mul(ctx, t11, n22, m11); submul(ctx, t11, n12, m21);
        mul(ctx, t12, n22, m12); submul(ctx, t12, n12, m22);
        mul(ctx, t21, n11, m21); submul(ctx, t21, n21, m11);
        mul(ctx, t22, n11, m22); submul(ctx, t22, n21, m12);
        if (ndet < 0)
        {
            neg(ctx, m11, t11);
            neg(ctx, m12, t12);
            neg(ctx, m21, t21);
            neg(ctx, m22, t22);
        }
        else
        {
            set(ctx, m11, t11);
            set(ctx, m12, t12);
            set(ctx, m21, t21);
            set(ctx, m22, t22);
        }

        fp_poly_set(ctx, a, A);
        fp_poly_set(ctx, b, B);
FLINT_ASSERT(a.length() > b.length());
        if (b.length() < 1)
            break;

        fp_poly_divrem_inplace(ctx, q, a, b); std::swap(a, b);
        fp_poly_submul(ctx, m11, q, m21);     std::swap(m11, m21);
        fp_poly_submul(ctx, m12, q, m22);     std::swap(m12, m22);
FLINT_ASSERT(a.length() > b.length());
    }


    while (b.length() > 0)
    {
        fp_poly_divrem_inplace(ctx, q, a, b); std::swap(a, b);
        fp_poly_submul(ctx, m11, q, m21);     std::swap(m11, m21);
        fp_poly_submul(ctx, m12, q, m22);     std::swap(m12, m22);
    }

    FLINT_ASSERT(a.length() > 0);

    fp_tmp_elem lcinv(ctx);
    fp_inv(ctx, lcinv.data(), a.data() + N*(a.length()-1));
    fp_vec_scalar_mul(ctx, m11.data(), m11.data(), m11.length(), lcinv.data());
    fp_vec_scalar_mul(ctx, m12.data(), m12.data(), m12.length(), lcinv.data());
    fp_poly_make_monic_with_lcinv(ctx, g, a, lcinv.data());
    return;
#endif
}

void fp_poly_gcdinv(fp_ring& ctx,
    fp_poly& g, fp_poly& m11,
    const fp_poly& a,
    const fp_poly& b)
{
    return generic_poly_gcdinv<fp_ring>(ctx, g, m11, a, b);
}

