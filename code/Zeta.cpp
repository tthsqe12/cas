#include "globalstate.h"
#include "eval.h"
#include "code.h"
#include "flintarb_wrappers.h"
#include "flint/arith.h"


ex dcode_sZeta(er e)
{
//std::cout << "dcode_sZeta: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sZeta.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er e1 = echild(e,1);
    if (eis_intsm(e1))
    {
        slong a = eintsm_get(e1);
        if (a < -10000 || a > 1000)
            return ecopy(e);

        if (a <= 0)
        {
            if (a & 1)
            {
                xfmpq_t r;
                arith_bernoulli_number(r.data, 1 - a);
                fmpq_div_si(r.data, r.data, a - 1);
                return emake_rat_move(r.data);
            }
            else
            {
                return (a == 0) ? emake_crat(-1,2) : emake_cint(0);
            }
        }
        else
        {
            if (a & 1)
            {
                return (a == 1) ? gs.const_complexinfinity.copy() : ecopy(e);
            }
            else
            {
                // Zeta[a] = (-1)^(a/2-1)*2^(a-1)*BernoulliB[a]/a!*Pi^a
                xfmpq_t u, v;
                arith_bernoulli_number(u.data, a);
                fmpz_fac_ui(fmpq_numref(v.data), a);
                fmpq_div_2exp(v.data, v.data, a - 1);
                fmpq_div(u.data, u.data, v.data);
                if ((a & 2) == 0)
                    fmpz_neg(fmpq_numref(u.data), fmpq_numref(u.data));
                uex p(emake_node(gs.sym_sPower.copy(), gs.sym_sPi.copy(), ecopy(e1)));
                ex t = emake_rat_move(u.data);
                return emake_node(gs.sym_sTimes.copy(), t, p.release());
            }
        }
    }
    else if (eis_real(e1))
    {
        ex z = emake_real();
        slong p = ereal_number(e1).wprec();
        arb_zeta(ereal_data(z), ereal_data(e1), p + EXTRA_PRECISION_BASIC);
        return efix_real(z);
    }
    else if (eis_double(e1))
    {
        slong p = 53;
        arb_set_d(gs.tmpreal[7].data, edouble_number(e1));
        arb_zeta(gs.tmpreal[6].data, gs.tmpreal[7].data, p + EXTRA_PRECISION_BASIC);
        return emake_double(arf_get_d(arb_midref(gs.tmpreal[6].data), ARF_RND_NEAR));
    }
    else
    {
        return ecopy(e);
    }
}
