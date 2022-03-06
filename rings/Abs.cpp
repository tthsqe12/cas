#include "globalstate.h"
#include "arithmetic.h"
#include "ex_print.h"

ex num_AbsInt(er X)
{
    assert(eis_int(X));

    if (fmpz_sgn(eint_data(X)) >= 0)
        return ecopy(X);

    ex Z = emake_int();
    fmpz_neg(eint_data(Z), eint_data(X));
    return efix_int(Z);
}

ex num_AbsRat(er X)
{
    assert(eis_rat(X));

    if (fmpq_sgn(erat_data(X)) >= 0)
        return ecopy(X);

    ex Z = emake_rat();
    fmpq_neg(erat_data(Z), erat_data(X));
    return efix_rat(Z);
}

ex num_AbsDouble(er X)
{
    assert(eis_double(X));

    double x = edouble_number(X);

    if (x >= 0)
        return ecopy(X);

    return emake_double(-x);
}

ex num_AbsReal(er X)
{
    assert(eis_real(X));

    if (arb_is_nonnegative(ereal_data(X)))
        return ecopy(X);

    ex Z = emake_real();
    arb_abs(ereal_data(Z), ereal_data(X));
    return efix_real(Z);
}

ex num_Abs(er X)
{
//std::cout << "num_Abs: " << ex_tostring_full(X) << std::endl;
    assert(eis_number(X));

    uint32_t tx = etype(X);
    switch (tx)
    {
        case ETYPE_INT:
            return num_AbsInt(X);
        case ETYPE_RAT:
            return num_AbsRat(X);
        case ETYPE_DOUBLE:
            return num_AbsDouble(X);
        case ETYPE_REAL:
            return num_AbsReal(X);
        default:
            return ecopy(X);
    }
}