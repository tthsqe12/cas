#include "fp_poly.h"
#include "fmpz.h"
#include "generic/generic_parser.h"

std::ostream& _fp_poly_write(
    std::ostream& o,
    const fp_ring_base& ctx,
    _fp_poly a,
    const char* var)
{
    ulimb N = ctx.stride();
    ulimb an = a.len;
    const ulimb* ad = a.coeffs;
    bool sum_first = true;
    for (ulimb i = an; i > 0; i--)
    {
        if (i < an && fp_is_zero(ctx, ad + N*(i-1)))
            continue;
        if (!sum_first)
            o << " + ";
        sum_first = false;
        bool prod_first = false;
        if (fp_is_one(ctx, ad + N*(i-1)))
            prod_first = true;
        else
            fp_write(o, ctx, ad + N*(i-1));
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



// faster parser of only the fp_poly_write format
ulimb fp_poly_read(fp_ring& ctx, fp_poly& x, const char* s, ulimb sn)
{
    fp_poly_zero(ctx, x);
    return -1;
}

void fp_poly_set_strn(fp_ring& ctx, fp_poly& x, const char* s, ulimb sn)
{
    if (sn == fp_poly_read(ctx, x, s, sn))
        return;
    elem_parser<fp_ring, fp_poly> E(ctx);
    fp_poly_gen(ctx, x);
    E.add_terminal(x, "x");
    E.parse(x, s, sn);
}

std::ostream& fp_poly_poly_write(
    std::ostream& o,
    const fp_ring_base& ctx,
    const fp_poly_poly& f,
    const char* var)
{
    o << "{";
    for (ulimb i = 0; i < f.length(); i++)
    {
        if (i > 0)
            o << ", ";
        fp_poly_write(o, ctx, f.coeff(i), var);
    }
    o << "}";
    return o;
}

std::ostream& fp_poly_product_write(
    std::ostream& o,
    const fp_ring_base& ctx,
    const fp_poly_product& f,
    const char* var)
{
    o << "1";
    for (ulimb i = 0; i < f.length(); i++)
        fp_poly_write(o << " * (", ctx, f.base(i), var) << ")^" << f.exp(i);
    return o;
}



