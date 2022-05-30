#include <cmath>
#include <cfloat>

#include "timing.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "flint/arith.h"

ex dcode_sHarmonicNumber(er e)
{
//std::cout << "dcode_sHarmonicNumber: " << ex_tostring_full(e) << std::endl;

    if (elength(e) != 1)
        return ecopy(e);

    if (eis_int(echild(e,1)))
    {
        er x = echild(e,1);
        if (fmpz_sgn(eint_data(x)) < 0)
            return gs.const_complexinfinity.copy();

        if (!eint_is_intsm(x))
            return ecopy(x);

        ulong a = eintsm_get(x);
        if (a < 4)
        {
            if (a < 2)
                return ecopy(x);
            else if (a == 2)
                return emake_crat(3, 2);
            else
                return emake_crat(11, 6);
        }
        else
        {
            ex z = emake_rat();
            fmpq_harmonic_ui(erat_data(z), a);
            return efix_rat(z);
        }
    }
    else
    {
        return ecopy(e);
    }
}


ex dcode_sFibonacci(er e)
{
//std::cout << "dcode_sFibonacci: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFibonacci.get()));

    if (elength(e) == 1)
    {
        er x = echild(e,1);
        if (eis_int(x) && fmpz_fits_si(eint_data(x)))
        {
            slong n = fmpz_get_si(eint_data(x));
            ex z = emake_int();
            if (n >= 0)
            {
                fmpz_fib_ui(eint_data(z), n);
            }
            else
            {
                fmpz_fib_ui(eint_data(z), -n);
                if (n % 2 == 0)
                    fmpz_neg(eint_data(z), eint_data(z));
            }
            return efix_int(z);
        }
        else
        {
            return ecopy(e);
        }
    }
    else if (elength(e) == 2)
	{
        return ecopy(e);
    }
	else
	{
		return _handle_message_argt(e, (1 << 0) + (2 << 8));
	}
}
