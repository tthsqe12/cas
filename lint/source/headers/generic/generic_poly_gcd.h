#pragma once
#include "generic/generic_poly_divrem.h"

template <typename RingT>
void _generic_poly_mat22_one(
    RingT& ctx,
    typename RingT::_poly_mat22_t& M)
{
    M.det = 1;
    M._11.len = 1;
    M._12.len = 0;
    M._21.len = 0;
    M._22.len = 1;
    one(ctx, M._11.coeffs);
    one(ctx, M._22.coeffs);
}

#define POLY_HGCD_CUTOFF 5

template <typename RingT>
void _generic_poly_mat22_mul(
    RingT& ctx,
    typename RingT::_poly_mat22_t& C,
    typename RingT::_poly_mat22_t& A,
    typename RingT::_poly_mat22_t& B,
    typename RingT::_poly_t& T0,
    typename RingT::_poly_t& T1)
{
    ulimb l = 5;
    C.det = A.det*B.det;
    if (A._11.len < l || A._12.len < l || A._21.len < l || A._22.len < l ||
        B._11.len < l || B._12.len < l || B._21.len < l || B._22.len < l)
    {
        _generic_poly_mul<RingT>(ctx, C._11, A._11, B._11);
        _generic_poly_mul<RingT>(ctx, T0, A._12, B._21);
        _generic_poly_add<RingT>(ctx, C._11, C._11, T0);

        _generic_poly_mul<RingT>(ctx, C._12, A._11, B._12);
        _generic_poly_mul<RingT>(ctx, T0, A._12, B._22);
        _generic_poly_add<RingT>(ctx, C._12, C._12, T0);

        _generic_poly_mul<RingT>(ctx, C._21, A._21, B._11);
        _generic_poly_mul<RingT>(ctx, T0, A._22, B._21);
        _generic_poly_add<RingT>(ctx, C._21, C._21, T0);

        _generic_poly_mul<RingT>(ctx, C._22, A._21, B._12);
        _generic_poly_mul<RingT>(ctx, T0, A._22, B._22);
        _generic_poly_add<RingT>(ctx, C._22, C._22, T0);
    }
    else
    {
        _generic_poly_sub<RingT>(ctx, T0, A._11, A._21);
        _generic_poly_sub<RingT>(ctx, T1, B._22, B._12);
        _generic_poly_mul<RingT>(ctx, C._21, T0, T1);

        _generic_poly_add<RingT>(ctx, T0, A._21, A._22);
        _generic_poly_sub<RingT>(ctx, T1, B._12, B._11);
        _generic_poly_mul<RingT>(ctx, C._22, T0, T1);

        _generic_poly_sub<RingT>(ctx, T0, T0, A._11);
        _generic_poly_sub<RingT>(ctx, T1, B._22, T1);
        _generic_poly_mul<RingT>(ctx, C._12, T0, T1);

        _generic_poly_sub<RingT>(ctx, T0, A._12, T0);
        _generic_poly_mul<RingT>(ctx, C._11, T0, B._22);

        _generic_poly_mul<RingT>(ctx, T0, A._11, B._11);

        _generic_poly_add<RingT>(ctx, C._12, T0, C._12);
        _generic_poly_add<RingT>(ctx, C._21, C._12, C._21);
        _generic_poly_add<RingT>(ctx, C._12, C._12, C._22);
        _generic_poly_add<RingT>(ctx, C._22, C._21, C._22);
        _generic_poly_add<RingT>(ctx, C._12, C._12, C._11);
        _generic_poly_sub<RingT>(ctx, T1, T1, B._21);
        _generic_poly_mul<RingT>(ctx, C._11, A._22, T1);

        _generic_poly_sub<RingT>(ctx, C._21, C._21, C._11);
        _generic_poly_mul<RingT>(ctx, C._11, A._12, B._21);

        _generic_poly_add<RingT>(ctx, C._11, C._11, T0);
    }
}
    
// {A, B} = M^-1*{a, b}
template <typename RingT>
void _generic_poly_mat22_apply_inverse(
    RingT& ctx,
    typename RingT::_poly_t& A,
    typename RingT::_poly_t& B,
    typename RingT::_poly_mat22_t& M,
    typename RingT::_poly_t& a,
    typename RingT::_poly_t& b,
    typename RingT::_poly_t& T) // temp
{
    FLINT_ASSERT(M.det == 1 || M.det == -1);

    _generic_poly_mul<RingT>(ctx, B, M._21, a);
    _generic_poly_mul<RingT>(ctx, T, M._11, b);
    M.det < 0 ? _generic_poly_sub<RingT>(ctx, B, B, T) : _generic_poly_sub<RingT>(ctx, B, T, B);

    _generic_poly_mul<RingT>(ctx, A, M._22, a);
    _generic_poly_mul<RingT>(ctx, T, M._12, b);
    M.det < 0 ? _generic_poly_sub<RingT>(ctx, A, T, A) : _generic_poly_sub<RingT>(ctx, A, A, T);
}

template <typename RingT>
void _generic_poly_mat22_mul_elementary(
    RingT& ctx,
    typename RingT::_poly_mat22_t& R,
    typename RingT::_poly_t& q)
{
    _generic_poly_addmul<RingT>(ctx, R._12, q, R._11);
    std::swap(R._11, R._12);
    _generic_poly_addmul<RingT>(ctx, R._22, q, R._21);
    std::swap(R._21, R._22);
    R.det = -R.det;
}


template <typename RingT>
void _generic_poly_mat22_set(
    RingT& ctx,
    typename RingT::_poly_mat22_t& M,
    typename RingT::_poly_mat22_t& N)
{
    M.det = N.det;
    _generic_poly_set<RingT>(ctx, M._11, N._11);
    _generic_poly_set<RingT>(ctx, M._12, N._12);
    _generic_poly_set<RingT>(ctx, M._21, N._21);
    _generic_poly_set<RingT>(ctx, M._22, N._22);
}


/*
    hgcd: M.{A,B} = {a,b} with B.len <= a.len/2 < A.len
          The inplace version has aliased {A,B} and {a,b}
*/
template <typename RingT>
bool _generic_poly_hgcd_inplace_basecase(
    RingT& ctx,
    typename RingT::_poly_mat22_t& M,
    typename RingT::_poly_t& A,
    typename RingT::_poly_t& B,
    typename RingT::_poly_t& Q) // temp space
{
    FLINT_ASSERT(A.len > B.len);
    _generic_poly_mat22_one(ctx, M);
    bool swapped = false;
    ulimb m = A.len/2;
    while (B.len > m)
    {
        _generic_poly_divrem_inplace<RingT, 2>(ctx, Q, A, B);
        std::swap(A, B);
        _generic_poly_mat22_mul_elementary<RingT>(ctx, M, Q);
        swapped = !swapped;
    }
    return swapped;
}

template <typename RingT>
bool _generic_poly_hgcd_basecase(
    RingT& ctx,
    typename RingT::_poly_mat22_t& M,
    typename RingT::_poly_t& A,
    typename RingT::_poly_t& B,
    typename RingT::_poly_t& a,
    typename RingT::_poly_t& b,
    typename RingT::_poly_t& Q) // temp space
{
    _generic_poly_set<RingT>(ctx, A, a);
    _generic_poly_set<RingT>(ctx, B, b);
    return _generic_poly_hgcd_inplace_basecase<RingT>(ctx, M, A, B, Q);
}

template <typename RingT>
void _generic_poly_hgcd_inplace_recursive(
    RingT& ctx,
    typename RingT::_poly_mat22_t* M, // could be nullptr
    typename RingT::_poly_t& a,
    typename RingT::_poly_t& b)
{
    ulimb N = ctx.stride();
    const slimb m = a.len/2;

    if (b.len <= m)
    {
        if (M)
            _generic_poly_mat22_one<RingT>(ctx, *M);
        return;
    }

    // read only
    typename RingT::_poly_t a0, b0, s, t, c0, d0;

    // mutable
    typename RingT::_poly_t a2, b2, q, T0;
    typename RingT::_poly_mat22_t R, S;

    tmp_allocator push;
    a2.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
    b2.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
     q.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
    T0.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);

    R._11.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    R._12.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    R._21.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    R._22.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._11.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._12.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._21.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._22.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);

    _generic_poly_shift_right_readonly<RingT>(ctx, a0, a, m);
    _generic_poly_shift_right_readonly<RingT>(ctx, b0, b, m);
    _generic_poly_truncate_readonly<RingT>(ctx, s, a, m);
    _generic_poly_truncate_readonly<RingT>(ctx, t, b, m);

    if (a0.len < POLY_HGCD_CUTOFF)
        _generic_poly_hgcd_inplace_basecase<RingT>(ctx, R, a0, b0, q);
    else
        _generic_poly_hgcd_inplace_recursive<RingT>(ctx, &R, a0, b0);

    _generic_poly_mat22_apply_inverse<RingT>(ctx, a2, b2, R, s, t, T0);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, a2, a0, m);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, b2, b0, m);

    if (b2.len <= m)
    {
        _generic_poly_set<RingT>(ctx, a, a2);
        _generic_poly_set<RingT>(ctx, b, b2);
        if (M)
            _generic_poly_mat22_set<RingT>(ctx, *M, R);
        return;
    }

    slimb k = 2*m - b2.len + 1;

    _generic_poly_divrem_inplace<RingT>(ctx, q, a2, b2);

    _generic_poly_shift_right_readonly<RingT>(ctx, c0, b2, k);
    _generic_poly_shift_right_readonly<RingT>(ctx, d0, a2, k);
    _generic_poly_truncate_readonly<RingT>(ctx, s, b2, k);
    _generic_poly_truncate_readonly<RingT>(ctx, t, a2, k);

    if (c0.len < POLY_HGCD_CUTOFF)
        _generic_poly_hgcd_inplace_basecase<RingT>(ctx, S, c0, d0, T0);
    else 
        _generic_poly_hgcd_inplace_recursive<RingT>(ctx, &S, c0, d0);

    _generic_poly_mat22_apply_inverse<RingT>(ctx, a, b, S, s, t, T0);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, a, c0, k);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, b, d0, k);

    if (M)
    {
        _generic_poly_mat22_mul_elementary<RingT>(ctx, R, q);
        _generic_poly_mat22_mul<RingT>(ctx, *M, R, S, a2, b2);
    }
}

template <typename RingT>
void _generic_poly_hgcd_recursive(
    RingT& ctx,
    typename RingT::_poly_mat22_t* M, // could be nullptr
    typename RingT::_poly_t& A,
    typename RingT::_poly_t& B,
    typename RingT::_poly_t& a,
    typename RingT::_poly_t& b)
{
    ulimb N = ctx.stride();
    const slimb m = a.len/2;

    if (b.len <= m)
    {
        _generic_poly_set<RingT>(ctx, A, a);
        _generic_poly_set<RingT>(ctx, B, b);
        if (M)
            _generic_poly_mat22_one<RingT>(ctx, *M);
        return;
    }

    // read only
    typename RingT::_poly_t a0, b0, s, t;

    // mutable
    typename RingT::_poly_t a2, b2, a3, b3, q, T0;
    typename RingT::_poly_mat22_t R, S;

    tmp_allocator push;
    a2.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
    b2.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
    a3.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
    b3.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
     q.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);
    T0.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1), N);

    R._11.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    R._12.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    R._21.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    R._22.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._11.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._12.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._21.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);
    S._22.coeffs = push.recursive_alloc<typename RingT::coeff_t>((a.len+1) ,N);

    _generic_poly_shift_right_readonly<RingT>(ctx, a0, a, m);
    _generic_poly_shift_right_readonly<RingT>(ctx, b0, b, m);
    _generic_poly_truncate_readonly<RingT>(ctx, s, a, m);
    _generic_poly_truncate_readonly<RingT>(ctx, t, b, m);

    if (a0.len < POLY_HGCD_CUTOFF)
        _generic_poly_hgcd_basecase<RingT>(ctx, R, a3, b3, a0, b0, q);
    else
        _generic_poly_hgcd_recursive<RingT>(ctx, &R, a3, b3, a0, b0);

    _generic_poly_mat22_apply_inverse<RingT>(ctx, a2, b2, R, s, t, T0);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, a2, a3, m);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, b2, b3, m);

    if (b2.len <= m)
    {
        _generic_poly_set<RingT>(ctx, A, a2);
        _generic_poly_set<RingT>(ctx, B, b2);
        if (M)
            _generic_poly_mat22_set(ctx, *M, R);
        return;
    }

    slimb k = 2*m - b2.len + 1;

    _generic_poly_divrem_inplace<RingT>(ctx, q, a2, b2);

    _generic_poly_shift_right_readonly<RingT>(ctx, a0, b2, k);
    _generic_poly_shift_right_readonly<RingT>(ctx, b0, a2, k);
    _generic_poly_truncate_readonly<RingT>(ctx, s, b2, k);
    _generic_poly_truncate_readonly<RingT>(ctx, t, a2, k);

    if (a0.len < POLY_HGCD_CUTOFF)
        _generic_poly_hgcd_inplace_basecase<RingT>(ctx, S, a0, b0, T0);
    else 
        _generic_poly_hgcd_inplace_recursive<RingT>(ctx, &S, a0, b0);

    _generic_poly_mat22_apply_inverse<RingT>(ctx, A, B, S, s, t, T0);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, A, a0, k);
    _generic_poly_add_inplace_shift_left<RingT>(ctx, B, b0, k);

    if (M)
    {
        _generic_poly_mat22_mul_elementary<RingT>(ctx, R, q);
        _generic_poly_mat22_mul<RingT>(ctx, *M, R, S, a2, b2);
    }
}


template <typename RingT>
int generic_poly_hgcd(
    RingT& ctx,
    typename RingT::poly_t& m11,
    typename RingT::poly_t& m12,
    typename RingT::poly_t& m21,
    typename RingT::poly_t& m22,
    typename RingT::poly_t& x,
    typename RingT::poly_t& y,
    const typename RingT::poly_t& a,
    const typename RingT::poly_t& b)
{
    FLINT_ASSERT(a.length() > b.length());

    if (b.length() < 1)
    {
        set(ctx, x, a);
        set(ctx, y, b);
        one(ctx, m11); zero(ctx, m12);
        zero(ctx, m21); one(ctx, m22);
        return 1;
    }

    typename RingT::_poly_mat22_t M;
    typename RingT::_poly_t A, B, _a(a), _b(b);

    ulimb n = (a.length())/2;
    ulimb N = ctx.stride();
    M._11.coeffs = m11.fit_alloc(n+1, N);
    M._12.coeffs = m12.fit_alloc(n+1, N);
    M._21.coeffs = m21.fit_alloc(n+1, N);
    M._22.coeffs = m22.fit_alloc(n+1, N);
    A.coeffs = x.fit_alloc(a.length(), N);
    B.coeffs = y.fit_alloc(a.length(), N);

    if (n < 5)
    {
        tmp_allocator push;
        typename RingT::_poly_t Q;
        Q.coeffs = push.recursive_alloc<typename RingT::coeff_t>(a.length(), N);
        if (_generic_poly_hgcd_basecase<RingT>(ctx, M, A, B, _a, _b, Q))
        {
            // basecase swapped some pointers that it didn't own
            std::swap(x, y);
            std::swap(m11, m12);
            std::swap(m21, m22);
        }
    }
    else
    {
        _generic_poly_hgcd_recursive<RingT>(ctx, &M, A, B, _a, _b);
        // this version didn't swap anything it shouldn't have
    }

    x.set_length(A.len);
    y.set_length(B.len);
    m11.set_length(M._11.len);
    m12.set_length(M._12.len);
    m21.set_length(M._21.len);
    m22.set_length(M._22.len);

    FLINT_ASSERT(generic_poly_is_canonical<RingT>(ctx, m11));
    FLINT_ASSERT(generic_poly_is_canonical<RingT>(ctx, m12));
    FLINT_ASSERT(generic_poly_is_canonical<RingT>(ctx, m21));
    FLINT_ASSERT(generic_poly_is_canonical<RingT>(ctx, m22));
    FLINT_ASSERT(generic_poly_is_canonical<RingT>(ctx, x));
    FLINT_ASSERT(generic_poly_is_canonical<RingT>(ctx, y));

    return M.det;
}

template <typename RingT>
void generic_poly_gcd(
    RingT& ctx,
    typename RingT::poly_t& g,
    const typename RingT::poly_t& a,
    const typename RingT::poly_t& b)
{
    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();

    if (an < 1) {
        make_monic(ctx, g, b);
        return;
    }
    else if (bn < 1) {
        make_monic(ctx, g, a);
        return;
    }
    else if (an == 1 || bn == 1) {
        // might try inverting some lcs
        one(ctx, g);
        return;
    }

    typename RingT::_poly_t A, B, Q;

    tmp_allocator push;
    ulimb abn = std::max(an, bn);
    A.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    B.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    Q.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);

    _generic_poly_set<RingT>(ctx, A, (typename RingT::_poly_t)(a));
    _generic_poly_set<RingT>(ctx, B, (typename RingT::_poly_t)(b));

    if (A.len <= B.len)
    {
        if (A.len == B.len)
            _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
        std::swap(A, B);
    }

    while (A.len >= POLY_HGCD_CUTOFF)
    {
        _generic_poly_hgcd_inplace_recursive<RingT>(ctx, nullptr, A, B);
        if (B.len < 1)
            break;
        _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
        std::swap(A, B);        
    }

    while (B.len > 1)
    {
        _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
        std::swap(A, B);
    }

    FLINT_ASSERT(A.len >= 1);

    if (B.len < 1)
    {
        typename RingT::coeff_t* gd = g.fit_alloc_destroy(A.len, N);
        g.set_length(A.len);
        inv(ctx, Q.coeffs + N*0, A.coeffs + N*(A.len - 1));
        vec_scalar_mul(ctx, gd, A.coeffs, A.len - 1, Q.coeffs + N*0);
        one(ctx, gd + N*(A.len - 1));
    }
    else
    {
        one(ctx, g);
    }
}


/*
    Want to find m:    [m11 m12] [a] = [g]
                       [m21 m22] [b]   [0]
    m starts out as [1 0; 0 1].
    If [A; B] is multiplied by M on the left, then so is m
    second column of m is not needed
*/
template<typename RingT>
void generic_poly_gcdinv(
    RingT& ctx,
    typename RingT::poly_t& g,
    typename RingT::poly_t& m11,
    const typename RingT::poly_t& a,
    const typename RingT::poly_t& b)
{
    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();

    if (bn < 1)
    {
        typename RingT::coeff_t* m11d = m11.fit_alloc_destroy(1, N);
        m11.set_length(1);

        if (an < 1)
        {
            zero(ctx, g);
            one(ctx, m11d + N*0);
            return;
        }

        inv(ctx, m11d + N*0, a.data() + N*(an-1));

        typename RingT::coeff_t* gd = g.fit_alloc(an, N);
        g.set_length(an);
        vec_scalar_mul(ctx, gd, a.data(), an-1, m11d + N*0);
        one(ctx, gd + N*(an-1));
        return;
    }

    tmp_allocator push;
    typename RingT::_poly_t A, B, Q, M11, M21;

    ulimb abn = std::max(an, bn);
    A.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    B.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    Q.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    M11.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    M21.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);

    _generic_poly_set<RingT>(ctx, A, (typename RingT::_poly_t)(a));
    _generic_poly_set<RingT>(ctx, B, (typename RingT::_poly_t)(b));
    _generic_poly_one<RingT>(ctx, M11);
    _generic_poly_zero<RingT>(ctx, M21);

    if (A.len <= B.len)
    {
        if (A.len == B.len)
            _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B); // no submul because m21 is zero
        std::swap(A, B);
        std::swap(M11, M21);
    }

    if (A.len >= POLY_HGCD_CUTOFF)
    {
        typename RingT::_poly_mat22_t R;
        typename RingT::_poly_t T0, T1;

        R._11.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        R._12.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        R._21.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        R._22.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        T0.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        T1.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);

        while (A.len >= POLY_HGCD_CUTOFF)
        {
            _generic_poly_hgcd_inplace_recursive<RingT>(ctx, &R, A, B);
            // multiply m by r^-1 on the left   [r22  -r12] [m11 m12]
            //                                  [-r21  r11] [m21 m22]
            _generic_poly_mul<RingT>(ctx, T0, R._22, M11);
            _generic_poly_mul<RingT>(ctx, T1, R._12, M21);
            _generic_poly_mul<RingT>(ctx, Q, R._21, M11);
            R.det < 0 ? _generic_poly_sub<RingT>(ctx, M11, T1, T0) : _generic_poly_sub<RingT>(ctx, M11, T0, T1);
            _generic_poly_mul<RingT>(ctx, T0, R._11, M21);
            R.det < 0 ? _generic_poly_sub<RingT>(ctx, M21, Q, T0) : _generic_poly_sub<RingT>(ctx, M21, T0, Q);
            if (B.len < 1)
                break;
            _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
            std::swap(A, B);
            _generic_poly_submul<RingT>(ctx, M11, Q, M21);
            std::swap(M11, M21);
        }
    }

    while (B.len > 0)
    {
        _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
        std::swap(A, B);
        _generic_poly_submul<RingT>(ctx, M11, Q, M21);
        std::swap(M11, M21);
    }

    FLINT_ASSERT(A.len > 0);

    typename RingT::coeff_t* lcinv = Q.coeffs;
    inv(ctx, lcinv, A.coeffs + N*(A.len-1));

    typename RingT::coeff_t* m11d = m11.fit_alloc_destroy(M11.len, N);
    m11.set_length(M11.len);
    vec_scalar_mul(ctx, m11d, M11.coeffs, M11.len, lcinv);

    typename RingT::coeff_t* gd = g.fit_alloc_destroy(A.len, N);
    g.set_length(A.len);
    vec_scalar_mul(ctx, gd, A.coeffs, A.len-1, lcinv);
    one(ctx, gd + N*(A.len-1));
}


/*
    Supposedly it is faster to only calculate m11 with gcdinv and find m12 with 
    m12 = (g - m11*a)/b. This is not done here.
*/
template<typename RingT>
void generic_poly_gcdx(
    RingT& ctx,
    typename RingT::poly_t& g,
    typename RingT::poly_t& m11,
    typename RingT::poly_t& m12,
    const typename RingT::poly_t& a,
    const typename RingT::poly_t& b)
{
    ulimb N = ctx.stride();
    ulimb an = a.length();
    ulimb bn = b.length();

    if (bn < 1)
    {
        typename RingT::coeff_t* m11d = m11.fit_alloc_destroy(1, N);
        m11.set_length(1);
        m12.set_length(0);

        if (an < 1)
        {
            zero(ctx, g);
            one(ctx, m11d + N*0);
            return;
        }

        inv(ctx, m11d + N*0, a.data() + N*(an-1));

        typename RingT::coeff_t* gd = g.fit_alloc(an, N);
        g.set_length(an);
        vec_scalar_mul(ctx, gd, a.data(), an-1, m11d + N*0);
        one(ctx, gd + N*(an-1));
        return;
    }

    tmp_allocator push;
    typename RingT::_poly_t A, B, Q, M11, M21, M12, M22;

    ulimb abn = std::max(an, bn);
    A.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    B.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    Q.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    M11.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    M21.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    M12.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
    M22.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);

    _generic_poly_set<RingT>(ctx, A, (typename RingT::_poly_t)(a));
    _generic_poly_set<RingT>(ctx, B, (typename RingT::_poly_t)(b));
    _generic_poly_one<RingT>(ctx, M11);
    _generic_poly_zero<RingT>(ctx, M21);
    _generic_poly_zero<RingT>(ctx, M12);
    _generic_poly_one<RingT>(ctx, M22);

    if (A.len <= B.len)
    {
        if (A.len == B.len)
        {
            _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
            _generic_poly_submul<RingT>(ctx, M12, Q, M22);            
        }
        std::swap(A, B);
        std::swap(M11, M21);
        std::swap(M12, M22);
    }

    if (A.len >= POLY_HGCD_CUTOFF)
    {
        typename RingT::_poly_mat22_t R;
        typename RingT::_poly_t T0, T1;

        T0.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        T1.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        R._11.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        R._12.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        R._21.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);
        R._22.coeffs = push.recursive_alloc<typename RingT::coeff_t>(abn, N);

        while (A.len >= POLY_HGCD_CUTOFF)
        {
            _generic_poly_hgcd_inplace_recursive<RingT>(ctx, &R, A, B);
            // multiply m by r^-1 on the left   [r22  -r12] [m11 m12]
            //                                  [-r21  r11] [m21 m22]
            _generic_poly_mul<RingT>(ctx, T0, R._22, M11);
            _generic_poly_mul<RingT>(ctx, T1, R._12, M21);
            _generic_poly_mul<RingT>(ctx, Q, R._21, M11);
            R.det < 0 ? _generic_poly_sub<RingT>(ctx, M11, T1, T0) : _generic_poly_sub<RingT>(ctx, M11, T0, T1);
            _generic_poly_mul<RingT>(ctx, T0, R._11, M21);
            R.det < 0 ? _generic_poly_sub<RingT>(ctx, M21, Q, T0) : _generic_poly_sub<RingT>(ctx, M21, T0, Q);

            _generic_poly_mul<RingT>(ctx, T0, R._22, M12);
            _generic_poly_mul<RingT>(ctx, T1, R._12, M22);
            _generic_poly_mul<RingT>(ctx, Q, R._21, M12);
            R.det < 0 ? _generic_poly_sub<RingT>(ctx, M12, T1, T0) : _generic_poly_sub<RingT>(ctx, M12, T0, T1);
            _generic_poly_mul<RingT>(ctx, T0, R._11, M22);
            R.det < 0 ? _generic_poly_sub<RingT>(ctx, M22, Q, T0) : _generic_poly_sub<RingT>(ctx, M22, T0, Q);

            if (B.len < 1)
                break;

            _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
            std::swap(A, B);
            _generic_poly_submul<RingT>(ctx, M11, Q, M21);
            std::swap(M11, M21);
            _generic_poly_submul<RingT>(ctx, M12, Q, M22);
            std::swap(M12, M22);
        }
    }

    while (B.len > 0)
    {
        _generic_poly_divrem_inplace<RingT>(ctx, Q, A, B);
        std::swap(A, B);
        _generic_poly_submul<RingT>(ctx, M11, Q, M21);
        std::swap(M11, M21);
        _generic_poly_submul<RingT>(ctx, M12, Q, M22);
        std::swap(M12, M22);
    }

    FLINT_ASSERT(A.len > 0);

    typename RingT::coeff_t* lcinv = Q.coeffs;
    inv(ctx, lcinv, A.coeffs + N*(A.len-1));

    typename RingT::coeff_t* m11d = m11.fit_alloc_destroy(M11.len, N);
    m11.set_length(M11.len);
    vec_scalar_mul(ctx, m11d, M11.coeffs, M11.len, lcinv);

    typename RingT::coeff_t* m12d = m12.fit_alloc_destroy(M12.len, N);
    m12.set_length(M12.len);
    vec_scalar_mul(ctx, m12d, M12.coeffs, M12.len, lcinv);

    typename RingT::coeff_t* gd = g.fit_alloc_destroy(A.len, N);
    g.set_length(A.len);
    vec_scalar_mul(ctx, gd, A.coeffs, A.len-1, lcinv);
    one(ctx, gd + N*(A.len-1));
}
