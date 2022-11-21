#pragma once
#include "generic_poly_scalar.h"

template <typename RingT>
void _generic_poly_mul(
    RingT& ctx,
    typename RingT::_poly_t& X,
    typename RingT::_poly_t A,
    typename RingT::_poly_t B)
{
    ulimb N = ctx.stride();
    ulimb an = A.len;
    ulimb bn = B.len;
    const ulimb* ad = A.coeffs;
    const ulimb* bd = B.coeffs;
    ulimb* xd = X.coeffs;

    if (UNLIKELY(an < 1 || bn < 1))
    {
        X.len = 0;
        return;
    }

    typename RingT::dotter_t dotter(ctx);
    ulimb xn = an + bn - 1;

    for (ulimb i = 0; i < xn; i++)
    {
        ulong jstart = i > bn - 1 ? i - (bn - 1) : 0;
        ulong jstop = std::min(i + 1, an);
        FLINT_ASSERT(jstart < jstop);

        ulimb j = jstart;
        dotter.mul(ctx, ad + N*j, bd + (i - j)*N);
        for (j++; j < jstop; j++)
            dotter.madd(ctx, ad + N*j, bd + (i - j)*N);

        dotter.reduce(ctx, xd + N*i);
    }

    X.len = vec_normalized_length(ctx, xd, xn);
}

template <typename RingT>
void _generic_poly_addmul(
    RingT& ctx,
    typename RingT::_poly_t& X,
    typename RingT::_poly_t A,
    typename RingT::_poly_t B)
{
    ulimb N = ctx.stride();
    ulimb an = A.len;
    ulimb bn = B.len;
    const ulimb* ad = A.coeffs;
    const ulimb* bd = B.coeffs;
    ulimb* xd = X.coeffs;
    ulimb xn = X.len;

    if (UNLIKELY(an < 1 || bn < 1))
        return;

    fp_dotter dotter(ctx);
    ulimb abn = an + bn - 1;

    for (ulimb i = 0; i < abn; i++)
    {
        ulong jstart = i > bn - 1 ? i - (bn - 1) : 0;
        ulong jstop = std::min(i + 1, an);
        FLINT_ASSERT(jstart < jstop);

        ulimb j = jstart;
        dotter.mul(ctx, ad + N*j, bd + (i - j)*N);
        for (j++; j < jstop; j++)
            dotter.madd(ctx, ad + N*j, bd + (i - j)*N);
        if (i < xn)
            dotter.add(ctx, xd + N*i);
        dotter.reduce(ctx, xd + N*i);
    }

    if (abn >= xn)
        X.len = fp_vec_normalized_length(ctx, xd, abn);

    FLINT_ASSERT(_generic_poly_is_canonical<RingT>(ctx, X));
}

template <typename RingT>
void _generic_poly_submul(
    RingT& ctx,
    typename RingT::_poly_t& X,
    typename RingT::_poly_t A,
    typename RingT::_poly_t B)
{
    ulimb N = ctx.stride();
    ulimb an = A.len;
    ulimb bn = B.len;
    const ulimb* ad = A.coeffs;
    const ulimb* bd = B.coeffs;
    ulimb* xd = X.coeffs;
    ulimb xn = X.len;

    if (UNLIKELY(an < 1 || bn < 1))
        return;

    fp_dotter dotter(ctx);
    ulimb abn = an + bn - 1;

    for (ulimb i = 0; i < abn; i++)
    {
        ulong jstart = i > bn - 1 ? i - (bn - 1) : 0;
        ulong jstop = std::min(i + 1, an);
        FLINT_ASSERT(jstart < jstop);

        ulimb j = jstart;
        dotter.mul(ctx, ad + N*j, bd + (i - j)*N);
        for (j++; j < jstop; j++)
            dotter.madd(ctx, ad + N*j, bd + (i - j)*N);
        if (i < xn)
            dotter.sub(ctx, xd + N*i);
        dotter.reduce(ctx, xd + N*i);
        neg(ctx, xd + N*i, xd + N*i);
    }

    if (abn >= xn)
        X.len = fp_vec_normalized_length(ctx, xd, abn);

    FLINT_ASSERT(_generic_poly_is_canonical<RingT>(ctx, X));
}

