#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"


class cmplx_fxn_ArcCosh : public cmplx_fxn {
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
        std::complex<double> z = std::acosh(x);
        return emake_cmplx(z.real(), z.imag());
    }

    ex eval_double_imag(double y, er def)
    {
        er s;
        if (y >= 0)
        {
            s = gs.const_double_halfpi.get();
        }
        else
        {
            s = gs.const_double_mhalfpi.get();
            y = -y;
        } 
        ex t = emake_double(asinh(y));
        return emake_cmplx(t, ecopy(s));
    }

    ex eval_arb_imag(const xarb_t&y, er def)
    {
        uex wx(emake_real());
        uex wy(emake_real());
        slong p = y.wprec() + EXTRA_PRECISION_BASIC;
        arb_asinh(ereal_data(wx.get()), y.data, p);
        arb_const_pi(ereal_data(wy.get()), p);
        arb_mul_2exp_si(ereal_data(wy.get()), ereal_data(wy.get()), -WORD(1));    

        if (arb_is_positive(y.data))
        {
        }
        else if (arb_is_negative(y.data))
        {
            arb_neg(ereal_data(wx.get()), ereal_data(wx.get()));
            arb_neg(ereal_data(wy.get()), ereal_data(wy.get()));
        }
        else
        {
            arb_abs(ereal_data(wx.get()), ereal_data(wx.get()));
        }
        
        return emake_cmplx(wx.release(), wy.release());
    }

    void acb_eval(acb_t w, const acb_t z, slong prec)
    {
        acb_acosh(w, z, prec);
    }
};


ex num_ArcCosh(er X, er def)
{
    switch (etype(X))
    {
    case ETYPE_INT:
    {
        if (fmpz_is_zero(eint_data(X)))
            return emake_cint(1);
    }
    case ETYPE_RAT:
    {
        return ecopy(def);
    }

    case ETYPE_DOUBLE:
    {
        double d = edouble_number(X);
        if (d >= 1)
        {
            return emake_double(acosh(d));
        }
        else if (d <= -1)
        {
            ex t = emake_double(std::acosh(-d));
            return emake_cmplx(t, gs.const_double_pi.copy());
        }
        else
        {
            ex t = emake_double(std::acos(d));
            return emake_cmplx(emake_cint(0), t);
        }
    }

    case ETYPE_REAL:
    {
        slong p = ereal_number(X).wprec() + EXTRA_PRECISION_BASIC;
        ex Z = emake_real();
        uex z(Z);
        xarb_t one;
        arb_one(one.data);
        if (arb_ge(ereal_data(X), one.data))
        {
            arb_acosh(ereal_data(Z), ereal_data(X), p);
            return efix_real(z.release());
        }
        arb_neg(one.data, one.data);
        if (arb_le(ereal_data(X), one.data))
        {
            arb_neg(ereal_data(Z), ereal_data(X));
            arb_acosh(ereal_data(Z), ereal_data(Z), p);
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), ereal_number(z.get()).wprec());
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        arb_acos(ereal_data(Z), ereal_data(X), p);
        ex t = efix_real(z.release());
        return emake_cmplx(emake_cint(0), t);
    }

    case ETYPE_CMPLX:
    {
        cmplx_fxn_ArcCosh f;
        return f.eval(X, def);
    }

    default:
    {
        return emake_nan_Indeterminate();
    }
    }
}