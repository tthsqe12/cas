#pragma once

template <typename RingT>
bool _generic_poly_is_canonical(
    RingT& ctx,
    typename RingT::_poly_t a)
{
    if (a.len == 0)
        return true;

    if (is_zero(ctx, a.coeffs + ctx.stride()*(a.len-1)))
    {
        std::cerr << "_poly leading coefficient is zero" << std::endl;
        return false;
    }

    return true;
}

template <typename RingT>
bool generic_poly_is_canonical(
    RingT& ctx,
    const typename RingT::poly_t& a)
{
    return _generic_poly_is_canonical<RingT>(ctx, (typename RingT::_poly_t)(a));
}

// B = div(A, x^m)
template <typename RingT>
inline void _generic_poly_shift_right_readonly(
    RingT& ctx,
    typename RingT::_poly_t& B,
    typename RingT::_poly_t& A,
    ulimb m)
{
    B.coeffs = A.coeffs + m*ctx.stride();
    B.len = A.len >= m ? A.len - m : 0;
}

// B = rem(A, x^m)
template <typename RingT>
inline void _generic_poly_truncate_readonly(
    RingT& ctx,
    typename RingT::_poly_t& B,
    typename RingT::_poly_t& A,
    ulimb m)
{
    B.coeffs = A.coeffs;
    B.len = A.len < m ? A.len : m;
}

// B = A
template <typename RingT>
inline void _generic_poly_set(
    RingT& ctx,
    typename RingT::_poly_t& B,
    typename RingT::_poly_t A)
{
    B.len = A.len;
    vec_set(ctx, B.coeffs, A.coeffs, A.len);
}

template <typename RingT>
inline void _generic_poly_zero(
    RingT& ctx,
    typename RingT::_poly_t& A)
{
    A.len = 0;
}

template <typename RingT>
inline void _generic_poly_one(
    RingT& ctx,
    typename RingT::_poly_t& A)
{
    A.len = 1;
    one(ctx, A.coeffs);
}

template <typename RingT>
void _generic_poly_add(
    RingT& ctx,
    typename RingT::_poly_t& X,
    typename RingT::_poly_t& A,
    typename RingT::_poly_t& B)
{
    ulimb an = A.len;
    ulimb bn = B.len;
    const ulimb* ad = A.coeffs;
    const ulimb* bd = B.coeffs;
    ulimb* xd = X.coeffs;
    ulimb N = ctx.stride();

    if (an > bn)
    {
        fp_vec_add(ctx, xd, ad, bd, bn);
        fp_vec_set(ctx, xd + N*bn, ad + N*bn, an - bn);
        X.len = an;
    }
    else if (an < bn)
    {
        fp_vec_add(ctx, xd, ad, bd, an);
        fp_vec_set(ctx, xd + N*an, bd + N*an, bn - an);
        X.len = bn;
    }
    else
    {
        fp_vec_add(ctx, xd, ad, bd, bn);
        X.len = fp_vec_normalized_length(ctx, xd, bn);
    }
}

// a += b*X^e
template <typename RingT>
void _generic_poly_add_inplace_shift_left(
    RingT& ctx,
    typename RingT::_poly_t& A,
    typename RingT::_poly_t& B,
    ulimb e)
{
    ulimb an = A.len;
    ulimb bn = B.len;
    ulimb* ad = A.coeffs;
    ulimb* bd = B.coeffs;
    ulimb N = ctx.stride();

    FLINT_ASSERT(_generic_poly_is_canonical<RingT>(ctx, A));
    FLINT_ASSERT(_generic_poly_is_canonical<RingT>(ctx, B));

    if (an >= bn + e)
    {
        //    |   bn   |  e  |
        //  |       an       |
        fp_vec_add(ctx, ad + N*e, ad + N*e, bd, bn);
        A.len = fp_vec_normalized_length(ctx, ad, an);
    }
    else
    {
        if (bn < 1)
            return;

        if (an >= e)
        {
            //  |   bn    |   e   |
            //        |    an     |
            fp_vec_add(ctx, ad + N*e, ad + N*e, bd, an - e);
            fp_vec_set(ctx, ad + N*an, bd + N*(an - e), bn + e - an);
        }
        else
        {
            //  |   bn  |    e    |
            //            |  an   |
            fp_vec_zero(ctx, ad + N*an, e - an);
            fp_vec_set(ctx, ad + N*e, bd, bn);
        }
        A.len = bn + e;
    }
    FLINT_ASSERT(_generic_poly_is_canonical<RingT>(ctx, A));
}

template <typename RingT>
void _generic_poly_sub(
    RingT& ctx,
    typename RingT::_poly_t& X,
    typename RingT::_poly_t A,
    typename RingT::_poly_t B)
{
    ulimb an = A.len;
    ulimb bn = B.len;
    const ulimb* ad = A.coeffs;
    const ulimb* bd = B.coeffs;
    ulimb* xd = X.coeffs;
    ulimb N = ctx.stride();

    if (an > bn)
    {
        fp_vec_sub(ctx, xd, ad, bd, bn);
        fp_vec_set(ctx, xd + N*bn, ad + N*bn, an - bn);
        X.len = an;
    }
    else if (an < bn)
    {
        fp_vec_sub(ctx, xd, ad, bd, an);
        fp_vec_neg(ctx, xd + N*an, bd + N*an, bn - an);
        X.len = bn;
    }
    else
    {
        fp_vec_sub(ctx, xd, ad, bd, bn);
        X.len = fp_vec_normalized_length(ctx, xd, bn);
    }
}
