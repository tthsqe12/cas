#include "globalstate.h"
#include "ex.h"
#include <arb.h>

ex ncode_sKhinchin(er e, slong prec)
{
    if (!eis_sym(e, gs.sym_sKhinchin.get()))
        return ecopy(e);

    ex z = emake_real();
    arb_const_khinchin(ereal_data(z), prec);
    return z;
}
