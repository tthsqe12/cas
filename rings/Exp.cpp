#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"


class cmplx_fxn_Exp : public cmplx_fxn {
public:

    slong extra_prec()
    {
        return EXTRA_PRECISION_BASIC;
    }

    ex eval_exact(er z, er def)
    {
        return ecopy(def);
    }

    ex eval_double(std::complex<double> x, er def)
    {
        std::complex<double> r = exp(x);
        return emake_cmplx(r.real(), r.imag());
    }

    ex eval_double_imag(double y, er def)
    {
        return emake_cmplx(cos(y), sin(y));
    }

    ex eval_arb_imag(const xarb_t&y, er def)
    {
        uex wx(emake_real());
        uex wy(emake_real());
        slong p = y.wprec() + EXTRA_PRECISION_BASIC;
        arb_sin_cos(ereal_data(wy.get()), ereal_data(wx.get()), y.data, p);
        return emake_cmplx(wx.release(), wy.release());
    }

    void acb_eval(acb_t w, const acb_t z, slong prec)
    {
        acb_exp(w, z, prec);
    }
};


ex num_Exp(er X, er def)
{
    switch (etype(X))
    {
    case ETYPE_INT:
    {
        if (fmpz_is_zero(eint_data(X)))
            return emake_cint(1);
        else if (fmpz_is_one(eint_data(X)))
            return gs.sym_sE.copy();
    }
    case ETYPE_RAT:
    {
        return emake_node(gs.sym_sPower.copy(), gs.sym_sE.copy(), ecopy(X));
    }

    case ETYPE_DOUBLE:
    {
        return emake_double(exp(edouble_number(X)));
    }

    case ETYPE_REAL:
    {
        ex z = emake_real();
        slong p = ereal_number(X).wprec();
        arb_exp(ereal_data(z), ereal_data(X), p + EXTRA_PRECISION_BASIC);
        return efix_real(z);
    }

    case ETYPE_CMPLX:
    {
        cmplx_fxn_Exp f;
        return f.eval(X, def);
    }

    default:
    {
        return emake_nan_Indeterminate();
    }
    }
}