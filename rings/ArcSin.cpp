#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"


class cmplx_fxn_ArcSin : public cmplx_fxn {
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
        std::complex<double> z = std::asin(x);
        return emake_cmplx(z.real(), z.imag());
    }

    ex eval_double_imag(double y, er def)
    {
        ex t = emake_double(asinh(y));
        return emake_cmplx(emake_cint(0), t);
    }

    ex eval_arb_imag(const xarb_t&y, er def)
    {
        uex wy(emake_real());
        slong p = y.wprec() + EXTRA_PRECISION_BASIC;
        arb_asinh(ereal_data(wy.get()), y.data, p);
        ex t = efix_real(wy.release());
        return emake_cmplx(emake_cint(0), t);
    }

    void acb_eval(acb_t w, const acb_t z, slong prec)
    {
        acb_asin(w, z, prec);
    }
};


ex num_ArcSin(er X, er def)
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
        double d = edouble_number(X);
        if (d > 1)
        {
            ex t = emake_double(-acosh(d));
            return emake_cmplx(gs.const_double_halfpi.copy(), t);
        }
        if (d < -1)
        {
            ex t = emake_double(acosh(-d));
            return emake_cmplx(gs.const_double_mhalfpi.copy(), t);
        }
        return emake_double(asin(d));
    }
    case ETYPE_REAL:
    {
        slong p = ereal_number(X).wprec() + EXTRA_PRECISION_BASIC;
        ex Z = emake_real();
        uex z(Z);
        xarb_t one;
        arb_one(one.data);
        if (arb_gt(ereal_data(X), one.data))
        {
            arb_acosh(ereal_data(Z), ereal_data(X), p);
            arb_neg(ereal_data(Z), ereal_data(Z));                    
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), ereal_number(z.get()).wprec());
            arb_mul_2exp_si(ereal_data(pi.get()), ereal_data(pi.get()), -WORD(1));
            ex t = efix_real(z.release());
            return emake_cmplx(pi.release(), t);
        }
        arb_neg(one.data, one.data);
        if (arb_lt(ereal_data(X), one.data))
        {
            arb_neg(ereal_data(Z), ereal_data(X));
            arb_acosh(ereal_data(Z), ereal_data(Z), p);
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), ereal_number(z.get()).wprec());
            arb_neg(ereal_data(pi.get()), ereal_data(pi.get()));
            arb_mul_2exp_si(ereal_data(pi.get()), ereal_data(pi.get()), -WORD(1));
            ex t = efix_real(z.release());
            return emake_cmplx(pi.release(), t);
        }
        arb_asin(ereal_data(Z), ereal_data(X), p);
        return efix_real(z.release());
    }
    case ETYPE_CMPLX:
    {
        cmplx_fxn_ArcSin f;
        return f.eval(X, def);
    }

    default:
    {
        return ecopy(def);
    }
    }
}