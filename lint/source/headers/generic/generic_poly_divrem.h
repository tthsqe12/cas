#pragma once
#include "generic_poly_mul.h"

// A is updated with the remainder
// Q needs space an - (bn - 1)
// parameter L is the max size of a dot product over ctx, O(L^2) code gen :(
template <typename RingT, int L = 4>
void _generic_poly_divrem_inplace(
    RingT& ctx,
    typename RingT::_poly_t& Q,
    typename RingT::_poly_t& A,
    typename RingT::_poly_t& B)
{
    ulimb N = ctx.stride();
    ulimb an = A.len;
    ulimb bn = B.len;
    typename RingT::coeff_t* ad = A.coeffs;
    const typename RingT::coeff_t* bd = B.coeffs;
    typename RingT::coeff_t* qd = Q.coeffs;
    ulimb i = an - (bn - 1);    // A.len = i + bn - 1;
    Q.len = i;
    typename RingT::tmp_elem_t t(ctx);

    if (UNLIKELY(an < bn || bn < 1))
    {
        Q.len = 0;
        return;
    }
    if (UNLIKELY(bn == 1))
    {
        fp_inv(ctx, t.data(), bd + N*(bn - 1));
        fp_vec_scalar_mul(ctx, Q.coeffs, ad, an, t.data());
        A.len = 0;
        return;
    }

    typename RingT::divider_t divider(ctx);
    typename RingT::dotter_t dotter(ctx);

    divider.set_divisor(ctx, bd + N*(bn - 1)); // can throw

    while (i >= L && bn >= L && L >= 3) // bn not changing
    {
        i -= L;
/*      generate L (= 3) quotient terms:
                       q2 = a4/b2   q1 = (a3-q2*b1)/b2    q0 = (a2-q2*b0-q1*b1)/b2
                       q2*x^2        + q1*x              + q0
                    ---------------------------------
b2*x^2 + b1*x + b0  |  a4*x^4 + a3*x^3 + a2*x^2 + ...
*/
        divider.divides(ctx, qd + N*(i+L-1), ad + N*(i+L-1+bn-1));
        for (int k = L-2; k >= 0; k--)
        {
            int l = L-1;
            dotter.mul(ctx, qd + N*(i+l), bd + N*(bn-1-(l-k)));
            for (l--; l > k; l--)
                dotter.madd(ctx, qd + N*(i+l), bd + N*(bn-1-(l-k)));
            dotter.reduce(ctx, qd + N*(i+k));
            sub(ctx, qd + N*(i+k), ad + N*(i+k+bn-1), qd + N*(i+k));
            divider.divides(ctx, qd + N*(i+k), qd + N*(i+k));
        }

        for (int j = -(L-1); j < 0; j++)
        {
            int k = -j;
            dotter.mul(ctx, qd + N*(i+L-1-k), bd + N*(j+k));
            for (k++; k < L; k++)
                dotter.madd(ctx, qd + N*(i+L-1-k), bd + N*(j+k));
            dotter.reduce(ctx, t.data());
            sub(ctx, ad + N*(j+L-1+i), ad + N*(j+L-1+i), t.data());
        }

        for (ulimb j = 0; j < bn-(L-1); j++)
        {
            int k = 0;
            dotter.mul(ctx, qd + N*(i+L-1-k), bd + N*(j+k));
            for (k++; k < L; k++)
                dotter.madd(ctx, qd + N*(i+L-1-k), bd + N*(j+k));
            dotter.reduce(ctx, t.data());
            sub(ctx, ad + N*(j+L-1+i), ad + N*(j+L-1+i), t.data());
        }
    }

    FLINT_ASSERT(bn >= 2);
    while (i >= 2)
    {
        i -= 2;
        // generate two quotient terms
        ulimb* q0 = qd + (i+1)*N;
        ulimb* q1 = qd + (i+0)*N;
        divider.divides(ctx, q0, ad + N*(i+1+bn-1));
        mul(ctx, q1, q0, bd + N*(bn-2));
        sub(ctx, q1, ad + N*(i+0+bn-1), q1);
        divider.divides(ctx, q1, q1);
        mul(ctx, t.data(), q1, bd + N*0);
        sub(ctx, ad + N*i, ad + N*i, t.data());
        for (ulimb j = 0; j < bn-1; j++)
        {
            dotter.mul(ctx, q0, bd + N*(j+0));
            dotter.madd(ctx, q1, bd + N*(j+1));
            dotter.reduce(ctx, t.data());
            sub(ctx, ad + N*(j+1+i), ad + N*(j+1+i), t.data());
        }
    }

    while (i >= 1)
    {
        i -= 1;
        //  generate one quotient term
        ulimb* q0 = qd + N*i;
        divider.divides(ctx, q0, ad + N*(i+bn-1));
        for (ulimb j = 0; j < bn-1; j++)
        {
            mul(ctx, t.data(), q0, bd + N*j);
            sub(ctx, ad + N*(j+i), ad + N*(j+i), t.data());
        }
    }

    A.len = vec_normalized_length(ctx, ad, bn - 1);
}
