#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"

class cmplx_fxn_Log : public cmplx_fxn {
public:

    slong extra_prec()
    {
        return EXTRA_PRECISION_BASIC;
    }

    ex eval_exact(er z, er def)
    {
        if (!eis_zero(ecmplx_real(z)))
            return ecopy(def);

        if (num_Less2(eget_cint(0), ecmplx_imag(z)))
        {
            uex r(emake_node(gs.sym_sLog.copy(), ecopy(ecmplx_imag(z))));
            ex t = emake_cmplx(emake_cint(0), emake_crat(1,2));
            t = emake_node(gs.sym_sTimes.copy(), t, gs.sym_sPi.copy());
            return emake_node(gs.sym_sPlus.copy(), t, r.release());
        }
        else if (num_Less2(ecmplx_imag(z), eget_cint(0)))
        {
            uex r(emake_node(gs.sym_sLog.copy(), num_Minus(ecmplx_imag(z))));
            ex t = emake_cmplx(emake_cint(0), emake_crat(-1,2));
            t = emake_node(gs.sym_sTimes.copy(), t, gs.sym_sPi.copy());
            return emake_node(gs.sym_sPlus.copy(), t, r.release());
        }
        else
        {
            return ecopy(def);
        }
    }

    ex eval_double(std::complex<double> x, er def)
    {
        std::complex<double> r = log(x);
        return emake_cmplx(emake_double(r.real()), emake_double(r.imag()));
    }

    ex eval_double_imag(double y, er def)
    {
        if (y > 0)
            return emake_cmplx(emake_double(log(y)), gs.const_double_halfpi.copy());
        else if (y < 0)
            return emake_cmplx(emake_double(log(-y)), gs.const_double_mhalfpi.copy());
        else
            return gs.const_indeterminate.copy();
    }

    ex eval_arb_imag(const xarb_t&y, er def)
    {
        slong p = y.wprec() + EXTRA_PRECISION_BASIC;
        ex Z = emake_real();
        uex z(Z);
        uex pi(emake_real());
        arb_const_pi(ereal_data(pi.get()), p);
        arb_mul_2exp_si(ereal_data(pi.get()), ereal_data(pi.get()), -1);
        if (arb_is_positive(y.data))
        {
            arb_log(ereal_data(Z), y.data, p);
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        else if (arb_is_negative(y.data))
        {
            arb_neg(ereal_data(pi.get()), ereal_data(pi.get()));
            arb_neg(ereal_data(Z), y.data);
            arb_log(ereal_data(Z), ereal_data(Z), p);
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        else
        {
            return gs.const_indeterminate.copy();
        }
    }

    void acb_eval(acb_t w, const acb_t z, slong prec)
    {
        acb_log(w, z, prec);
    }
};



ex num_Log(er X, er def)
{
    switch (etype(X))
    {
    case ETYPE_INT:
    {
        if (fmpz_is_zero(eint_data(X)))
            return gs.const_minfinity.copy();
        else if (fmpz_is_one(eint_data(X)))
            return emake_cint(0);
        else
            return ecopy(def);
    }

    case ETYPE_RAT:
    {
        return ecopy(def);
    }

    case ETYPE_DOUBLE:
    {
        ex z = emake_double();
        double d = edouble_number(X);
        if (d > 0)
        {
            edouble_number(z) = log(d);
            return z;
        }
        else if (d < 0)
        {
            edouble_number(z) = log(-d);
            return emake_cmplx(z, gs.const_double_pi.copy());
        }
        else
        {
            eclear(z);
            return gs.const_indeterminate.copy();
        }
    }

    case ETYPE_REAL:
    {
        slong p = ereal_number(X).wprec() + EXTRA_PRECISION_BASIC;
        ex Z = emake_real();
        uex z(Z);
        if (arb_is_positive(ereal_data(X)))
        {
            arb_log(ereal_data(Z), ereal_data(X), p);
            return efix_real(z.release());
        }
        else if (arb_is_negative(ereal_data(X)))
        {
            arb_neg(ereal_data(Z), ereal_data(X));
            arb_log(ereal_data(Z), ereal_data(Z), p);
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), ereal_number(z.get()).wprec());
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        else
        {
            return emake_nan_Indeterminate();
        }
    }

    case ETYPE_CMPLX:
    {
        cmplx_fxn_Log f;
        return f.eval(X, def);
    }

    default:
    {
        return emake_nan_Indeterminate();
    }
    }
}