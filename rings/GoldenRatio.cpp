#include "globalstate.h"
#include "ex.h"
#include <arb.h>

ex ncode_sGoldenRatio(er e, slong prec)
{
    if (!eis_sym(e, gs.sym_sGoldenRatio.get()))
        return ecopy(e);

    ex z = emake_real();
    arb_sqrt_ui(ereal_data(z), 5, prec);
    arb_add_ui(ereal_data(z), ereal_data(z), 1, prec);
    arb_mul_2exp_si(ereal_data(z), ereal_data(z), -1);
    return z;
}
