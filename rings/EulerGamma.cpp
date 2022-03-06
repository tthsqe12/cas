#include "globalstate.h"
#include "ex.h"
#include <arb.h>

ex ncode_sEulerGamma(er e, slong prec)
{
    if (!eis_sym(e, gs.sym_sEulerGamma.get()))
        return ecopy(e);

    ex z = emake_real();
    arb_const_euler(ereal_data(z), prec);
    return z;
}
