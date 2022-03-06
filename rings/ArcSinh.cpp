#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"


class cmplx_fxn_ArcSinh : public cmplx_fxn {
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
        std::complex<double> z = std::asinh(x);
        return emake_cmplx(z.real(), z.imag());
    }

    ex eval_double_imag(double y, er def)
    {
        ex wx, wy;
        if (y >= 1)
        {
            wx = emake_double(acosh(y));
            wy = gs.const_double_halfpi.copy();
        }
        else if (y <= -1)
        {
            wx = emake_double(-acosh(-y));
            wy = gs.const_double_mhalfpi.copy();
        }
        else 
        {
            wy = emake_double(asin(y));
            wx = emake_cint(0);
        }
        return emake_cmplx(wx, wy);
    }

    ex eval_arb_imag(const xarb_t&y, er def)
    {
        slong p = y.wprec() + EXTRA_PRECISION_BASIC;
        ex Z = emake_real();
        uex z(Z);
        xarb_t one;
        arb_one(one.data);
        if (arb_gt(y.data, one.data))
        {
            arb_acosh(ereal_data(Z), y.data, p);
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), p);
            arb_mul_2exp_si(ereal_data(pi.get()), ereal_data(pi.get()), -WORD(1));
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        arb_neg(one.data, one.data);
        if (arb_lt(y.data, one.data))
        {
            arb_neg(ereal_data(Z), y.data);
            arb_acosh(ereal_data(Z), ereal_data(Z), p);
            arb_neg(ereal_data(Z), ereal_data(Z));
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), p);
            arb_neg(ereal_data(pi.get()), ereal_data(pi.get()));
            arb_mul_2exp_si(ereal_data(pi.get()), ereal_data(pi.get()), -WORD(1));
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        arb_asin(ereal_data(Z), y.data, p);
        ex t = efix_real(z.release());
        return emake_cmplx(emake_cint(0), t);
    }

    void acb_eval(acb_t w, const acb_t z, slong prec)
    {
        acb_asinh(w, z, prec);
    }
};


ex num_ArcSinh(er X, er def)
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
        return emake_double(asinh(edouble_number(X)));
    }
    case ETYPE_REAL:
    {
        slong p = ereal_number(X).wprec() + EXTRA_PRECISION_BASIC;
        ex Z = emake_real();
        arb_asinh(ereal_data(Z), ereal_data(X), p);
        return efix_real(Z);
    }
    case ETYPE_CMPLX:
    {
        cmplx_fxn_ArcSinh f;
        return f.eval(X, def);
    }

    default:
    {
        return ecopy(def);
    }
    }
}