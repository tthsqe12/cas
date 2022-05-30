#include "eval.h"
#include "code.h"
#include "ex_cont.h"
#include "ex_print.h"

ex dcode_sRe(er e)
{
//std::cout << "dcode_sRe: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sRe.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er e1 = echild(e,1);
    switch (etype(e1))
    {
        case ETYPE_INT:
        case ETYPE_RAT:
        case ETYPE_DOUBLE:
        case ETYPE_REAL:
        {
            return emake_cint(0);
        }

        case ETYPE_CMPLX:
        {
            return ecopy(ecmplx_real(e1));
        }

        default:
        {
            return ecopy(e);
        }
    }
}

ex dcode_sIm(er e)
{
//std::cout << "dcode_sIm: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sIm.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er e1 = echild(e,1);
    switch (etype(e1))
    {
        case ETYPE_INT:
        case ETYPE_RAT:
        case ETYPE_DOUBLE:
        case ETYPE_REAL:
        {
            return emake_cint(0);
        }

        case ETYPE_CMPLX:
        {
            return ecopy(ecmplx_imag(e1));
        }

        default:
        {
            return ecopy(e);
        }
    }
}
