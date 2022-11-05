#include "nmod_poly.h"
#include "generic_parser.h"

std::ostream& nmod_poly_write(std::ostream& o, const nmod_ring_base& ctx, const nmod_poly& a, const char* var)
{
    ulimb N = ctx.stride;
    ulimb an = a.length();
    const ulimb* ad = a.data();
    bool sum_first = true;
    for (ulimb i = an; i > 0; i--)
    {
        if (i < an && ui_vec_is_zero(ad + N*(i-1), N))
            continue;

        if (!sum_first)
            o << " + ";

        sum_first = false;

        bool prod_first = false;
        if (nmod_is_one(ctx, ad + N*(i-1)))
            prod_first = true;
        else
            mpn_write(o, ad + N*(i-1), N);
        if (i-1 > 0)
        {
            if (!prod_first)
                o << "*";
            o << var;
            if (i-1 > 1)
                o << "^" << i-1;
            prod_first = false;
        }
        if (prod_first)
            o << "1";
    }
    if (sum_first)
        o << "0";
    return o;
}

// faster parser of only the nmod_poly_write format
ulimb nmod_poly_read(nmod_ring& ctx, nmod_poly& x, const char* s, ulimb sn)
{
    nmod_poly_zero(ctx, x);
    return -1;
}

void nmod_poly_set_strn(nmod_ring& ctx, nmod_poly& x, const char* s, ulimb sn)
{
    if (sn == nmod_poly_read(ctx, x, s, sn))
        return;
    elem_parser<nmod_ring, nmod_poly> E(ctx);
    nmod_poly_gen(ctx, x);
    E.add_terminal(x, "x");
    E.parse(x, s, sn);
}

std::ostream& nmod_poly_product_write(
    std::ostream& o,
    const nmod_ring_base& ctx,
    const nmod_poly_product& f,
    const char* var)
{
    o << "1";
    for (ulimb i = 0; i < f.length(); i++)
        nmod_poly_write(o << " * (", ctx, f.base(i), var) << ")^" << f.exp(i);
    return o;
}

