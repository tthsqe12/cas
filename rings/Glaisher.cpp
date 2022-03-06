#include "globalstate.h"
#include "ex.h"
#include <arb.h>

ex ncode_sGlaisher(er e, slong prec)
{
    if (!eis_sym(e, gs.sym_sGlaisher.get()))
        return ecopy(e);

    ex z = emake_real();
    arb_const_glaisher(ereal_data(z), prec);
    return z;
}

