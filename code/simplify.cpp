/*
SimplifyCount[p_] :=
 Which[Head[p] === Symbol, 1,
  IntegerQ[p], 
  If[p == 0, 1, Floor[N[Log[2, Abs[p]]/Log[2, 10]]] + If[p > 0, 1, 2]],
  Head[p] === Rational, 
  SimplifyCount[Numerator[p]] + SimplifyCount[Denominator[p]] + 1,
  Head[p] === Complex, 
  SimplifyCount[Re[p]] + SimplifyCount[Im[p]] + 1, NumberQ[p], 2,
  True, SimplifyCount[Head[p]] + 
   If[Length[p] == 0, 0, Plus @@ (SimplifyCount /@ (List @@ p))]]
*/

#include "globalstate.h"
#include "code.h"

static ulong complexity(er e)
{
    if (eis_leaf(e))
    {
        if (eis_number(e))
        {
            if (eis_int(e))
            {
                return fmpz_sizeinbase(eint_data(e), 10) +
                       (fmpz_sgn(eint_data(e)) < 0);
            }
            else if (eis_rat(e))
            {
                return fmpz_sizeinbase(fmpq_numref(erat_data(e)), 10) +
                       fmpz_sizeinbase(fmpq_denref(erat_data(e)), 10) +
                       (fmpq_sgn(erat_data(e)) < 0) + 1;
            }
            else if (eis_cmplx(e))
            {
                return complexity(ecmplx_real(e)) + complexity(ecmplx_imag(e)) + 1;
            }
            else
            {
                return 2;
            }
        }
        else
        {
            return 1;
        }
    }
    else if (eis_node(e))
    {
        ulong s = complexity(echild(e,0));
        for (ulong i = 0; i < elength(e); i++)
            s += complexity(echild(e,i+1));
        return s;
    }
    else
    {
        return 1;
    }
}

ex dcode_sSimplify(er e)
{
//std::cout << "dcode_sSimplify: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSimplify.get()));

    if (elength(e) == 1)
    {
        er f = echild(e,1);
        return emake_list(ecopy(f), emake_int_ui(complexity(f)));
    }
    else
    {
        return ecopy(e);
    }
}


ex dcode_sFullSimplify(er e)
{
//std::cout << "dcode_sFullSimplify: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFullSimplify.get()));

    return ecopy(e);
}
