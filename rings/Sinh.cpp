#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"


class cmplx_fxn_Sinh : public cmplx_fxn {
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
        std::complex<double> r = sinh(x);
        return emake_cmplx(r.real(), r.imag());
    }

    ex eval_double_imag(double y, er def)
    {
        ex t = emake_double(sin(y));
        return emake_cmplx(emake_cint(0), t);
    }

    ex eval_arb_imag(const xarb_t&y, er def)
    {
        uex wy(emake_real());
        slong p = y.wprec() + EXTRA_PRECISION_BASIC;
        arb_sin(ereal_data(wy.get()), y.data, p);
        ex t = efix_real(wy.release());
        return emake_cmplx(emake_cint(0), t);
    }

    void acb_eval(acb_t w, const acb_t z, slong prec)
    {
        acb_sinh(w, z, prec);
    }
};


ex num_Sinh(er X, er def)
{
    switch (etype(X))
    {
    case ETYPE_INT:
    {
        if (fmpz_is_zero(eint_data(X)))
            return emake_cint(0);
    }
    case ETYPE_RAT:
    {
        return ecopy(def);
    }

    case ETYPE_DOUBLE:
    {
        return emake_double(sinh(edouble_number(X)));
    }

    case ETYPE_REAL:
    {
        ex z = emake_real();
        slong p = ereal_number(X).wprec();
        arb_sinh(ereal_data(z), ereal_data(X), p + EXTRA_PRECISION_BASIC);
        return efix_real(z);
    }

    case ETYPE_CMPLX:
    {
        cmplx_fxn_Sinh f;
        return f.eval(X, def);
    }

    default:
    {
        return emake_nan_Indeterminate();
    }
    }
}