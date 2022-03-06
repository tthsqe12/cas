#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"

ex ncode_sE(er e, slong prec)
{
    if (!eis_sym(e, gs.sym_sE.get()))
        return ecopy(e);

    ex z = emake_real();
    arb_const_e(ereal_data(z), prec);
    return z;
}
