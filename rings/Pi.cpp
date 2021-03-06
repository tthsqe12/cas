#include "globalstate.h"
#include "ex.h"
#include <arb.h>

ex ncode_sPi(er e, slong prec)
{
    if (!eis_sym(e, gs.sym_sPi.get()))
        return ecopy(e);

    ex z = emake_real();
    arb_const_pi(ereal_data(z), prec);
    return z;
}
