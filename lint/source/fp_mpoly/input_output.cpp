#include "fp_mpoly.h"

bool mpoly_ctx::write_monomial(
    std::ostream& o,
    bool first,
    const ulimb* a,
    ulimb bits,
    const char* x) const
{
    if (UNLIKELY(bits > FLINT_BITS))
    {
        ulimb v = 0; do {
            const ulimb* e = a + bits/FLINT_BITS*v;
            ulimb en = mpn_normalize(a, bits/FLINT_BITS);
            if (en == 0)
                continue;
            if (!first)
                o << "*";
            if (en == 1 && e[0] == 1) {
                o << x << v;
                first = false;
            } else {
                o << x << v << "^";
                mpn_write(o, e, en);
                first = false;
            }
        } while (++v < nvars());
    }
    else
    {
        ulimb mask = var_mask_sp(bits);
        auto shifts = var_shifts_sp(bits);
        auto offsets = var_offsets_sp(bits);
        ulimb v = 0; do {
            ulimb e = (a[offsets[v]]>>shifts[v])&mask;
            if (e == 0)
                continue;
            if (!first)
                o << "*";
            if (e == 1) {
                o << x << v;
                first = false;
            } else {
                o << x << v << "^" << e;
                first = false;
            }
        } while (++v < nvars());
    }
    return first;
}
    

std::ostream& fp_mpoly_write(
    std::ostream& o,
    const fp_mpoly_ring& ctx,
    const fp_mpoly& a,
    const char* x)
{
    ulimb M = ctx.m.stride(a.bits());
    ulimb N = ctx.c.stride();
    bool sum_first = true;
    for (ulimb i = 0; i < a.length(); i++)
    {
        if (!sum_first)
            o << " + ";
        sum_first = false;
        bool prod_first = fp_is_one(ctx.c, a.coeffs() + N*i);
        if (!prod_first)
            fp_write(o, ctx.c, a.coeffs() + N*i);
        prod_first = ctx.m.write_monomial(o, prod_first, a.exps() + M*i, a.bits(), x);
        if (prod_first)
            o << "1";
    }
    if (sum_first)
        o << "0";
    return o;
}

