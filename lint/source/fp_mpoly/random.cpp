#include "fp_mpoly.h"

void fp_mpoly_random_up_to_length_and_degree_bound(
    random_state& state,
    fp_mpoly_ring& ctx,
    fp_mpoly& X,
    ulimb l,
    ulimb d)
{
    ulimb N = ctx.c.stride();
    ulimb n = ctx.m.nvars();
    ulimb* e = ctx.m.leaf_alloc_ui(n);

    X.m.zero(ctx.m);
    ulimb* Xcoeffs = X.c.fit_alloc(l, N);
    for (ulimb i = 0; i < l; i++)
    {
        for (ulimb v = 0; v < n; v++)
            e[v] = state.get_mod(d+1);
        X.m.push_monomial(ctx.m, e);
        fp_random(state, ctx.c, Xcoeffs + N*i);
    }

    fp_mpoly_sort_terms(ctx, X);
    fp_mpoly_combine_like_terms(ctx, X);
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, X));
}

void fp_mpoly_random_up_to_length_and_degree_bits(
    random_state& state,
    fp_mpoly_ring& ctx,
    fp_mpoly& X,
    ulimb l,
    ulimb dbits)
{
    ulimb N = ctx.c.stride();
    ulimb n = ctx.m.nvars();
    fmpz* e = ctx.m.leaf_alloc_fmpz(n);

    X.m.zero(ctx.m);
    ulimb* Xcoeffs = X.c.fit_alloc(l, N);
    for (ulimb i = 0; i < l; i++)
    {
        for (ulimb v = 0; v < n; v++)
            fmpz_random_of_bits_unsigned(state, e[v], state.get_mod(dbits + 1));
        X.m.push_monomial(ctx.m, e);
        fp_random(state, ctx.c, Xcoeffs + N*i);
    }

    fp_mpoly_sort_terms(ctx, X);
    fp_mpoly_combine_like_terms(ctx, X);
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, X));
}
