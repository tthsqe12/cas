#include "fmpz.h"
#include "fp_poly.h"

void fp_poly_random_of_length(random_state& state, fp_ring& ctx, fp_poly& x, ulimb len)
{
    ulimb N = ctx.stride();

    if (len < 1)
    {
        x.set_length(0);
        return;
    }

    ulimb* xd = x.fit_alloc(len, N);
    ulimb i = 0;
    do {
        fp_random(state, ctx, xd + i*N);
    } while (++i < len);

    if (ui_vec_is_zero(xd + (len-1)*N, N))
        (xd + (len-1)*N)[0] = 1;

    x.set_length(len);
    return;
}

void fp_poly_random_monic_of_length(random_state& state, fp_ring& ctx, fp_poly& x, ulimb len)
{
    ulimb N = ctx.stride();

    if (len < 1)
    {
        x.set_length(0);
        return;
    }

    ulimb* xd = x.fit_alloc(len, N);
    ulimb i = len;
    fp_one(ctx, xd + (i-1)*N);
    for (i--; i > 0; i--)
        fp_random(state, ctx, xd + (i-1)*N);

    x.set_length(len);
    return;
}

void fp_poly_random_up_to_length(random_state& state, fp_ring& ctx, fp_poly& x, ulimb len)
{
    ulimb N = ctx.stride();
    ulimb* xd = x.fit_alloc(len, N);
    for (ulimb i = 0; i < len; i++)
        fp_random(state, ctx, xd + i*N);
    x.normalize_length(len, N);
    return;
}


