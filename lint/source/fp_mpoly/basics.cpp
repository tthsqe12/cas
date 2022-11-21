#include "fp_mpoly.h"
#include "misc.h"


void fmpz_vec_add(fmpz* x, const fmpz* a, const fmpz* b, ulimb len)
{
    for (ulimb i = 0; i < len; i++)
        fmpz_add(x[i], a[i], b[i]);
}

slimb fmpz_vec_max_bits(const fmpz* a, ulimb len) {
    ulimb bits = 0;
    slimb have_neg = 0;
    for (ulimb i = 0; i < len; i++)
    {
        have_neg |= a[i].data;
        bits = std::max(bits, nbits(a[i]));
    }
    return copy_sign(bits, have_neg);
}

ulimb ui_vec_max_bits(const ulimb* a, ulimb an) {
    ulimb m = 0;
    for (ulimb i = 0; i < an; i++)
        m |= a[i];
    return nbits(m);
}



void fp_mpoly_combine_like_terms(
    fp_mpoly_ring& ctx,
    fp_mpoly& X)
{
    ulimb M = ctx.m.stride(X.bits());
    ulimb N = ctx.c.stride();
    ulimb* Xexps = X.exps();
    ulimb* Xcoeffs = X.coeffs();
    slimb out = -1;

    for (slimb in = 0; in < X.length(); in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && mpoly_monomial_equal(Xexps + M*out, Xexps + M*in, M))
        {
            fp_add(ctx.c, Xcoeffs + N*out, Xcoeffs + N*out, Xcoeffs + N*in);
        }
        else
        {
            if (out < 0 || !fp_is_zero(ctx.c, Xcoeffs + N*out))
                out++;

            if (out != in)
            {
                mpoly_monomial_set(Xexps + M*out, Xexps + M*in, M);
                ui_vec_set(Xcoeffs + N*out, Xcoeffs + N*in, N);
            }
        }
    }

    if (out < 0 || !fp_is_zero(ctx.c, Xcoeffs + N*out))
        out++;

    X.set_length(out);
}


/*
    sort terms in [left, right) by exponent
    assuming that bits in position > pos are already sorted
*/
static void _mpoly_radix_sort(
    ulimb* Xexps, ulimb* Xcoeffs,
    ulimb left, ulimb right,
    ulimb pos, ulimb M, ulimb N)
{
    ulimb off = pos/FLINT_BITS;
    ulimb bit = pos%FLINT_BITS;
    ulimb mask = UWORD(1) << bit;
    ulimb mid, check;

    FLINT_ASSERT(left <= right);
    FLINT_ASSERT(pos < M*FLINT_BITS);

    /* do nothing on lists of 0 or 1 elements */
    if (left + 1 >= right)
        return;

    /* find first 'zero' */
    mid = left;
    while (mid < right && ((Xexps + M*mid)[off] & mask) != 0)
        mid++;

    /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                 [mid,right)    does match cmpmask in position pos 'zero' */
    check = mid;
    while (++check < right)
    {
        if (((Xexps + M*check)[off] & mask) != 0)
        {
            ui_vec_swap(Xcoeffs + N*mid, Xcoeffs + N*check, N);
            mpoly_monomial_swap(Xexps + M*check, Xexps + M*mid, M);
            mid++;
        }
    }

    if (pos > 0)
    {
        --pos;
        _mpoly_radix_sort(Xexps, Xcoeffs, left,  mid, pos, M, N);
        _mpoly_radix_sort(Xexps, Xcoeffs, mid, right, pos, M, N);
    }
}

/*
    sort the terms in A by exponent
    assuming that the exponents are valid (other than being in order)
*/
void fp_mpoly_sort_terms(
    fp_mpoly_ring& ctx,
    fp_mpoly& X)
{
    ulimb M = ctx.m.stride(X.bits());
    ulimb N = ctx.c.stride();
    _mpoly_radix_sort(X.exps(), X.coeffs(), 0, X.length(), M*FLINT_BITS - 1, M, N);
}

