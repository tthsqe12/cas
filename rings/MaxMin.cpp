#include "globalstate.h"
#include "arithmetic.h"
#include "ex_print.h"
#include "timing.h"




ex num_MinIntInt(er X, er Y)
{
    assert(eis_int(X));
    assert(eis_int(Y));

    if (fmpz_cmp(eint_data(X), eint_data(Y)) < 0)
        return ecopy(X);
    else
        return ecopy(Y);
}

ex num_Min(er X, er Y)
{
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_MinIntInt(X, Y);

        default:
            assert(false);
            return nullptr;
    }
}

ex num_MaxIntInt(er X, er Y)
{
    assert(eis_int(X));
    assert(eis_int(Y));

    if (fmpz_cmp(eint_data(X), eint_data(Y)) < 0)
        return ecopy(Y);
    else
        return ecopy(X);
}

ex num_Max(er X, er Y)
{
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_MaxIntInt(X, Y);

        default:
            assert(false);
            return nullptr;
    }
}
