#include "globalstate.h"
#include "eval.h"
#include "code.h"
#include "arithmetic.h"


ex dcode_sHypergeometricPFQ(er e)
{
//std::cout << "dcode_sHypergeometricPFQ: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sHypergeometricPFQ.get()));

    return ecopy(e);
}

ex dcode_sMeijerG(er e)
{
//std::cout << "dcode_sMeijerG: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMeijerG.get()));

    return ecopy(e);
}

ex dcode_iNumericHypergeometricPFQ(er e)
{
//std::cout << "dcode_iNumericHypergeometricPFQ: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_iNumericHypergeometricPFQ.get()));

    return ecopy(e);
}
