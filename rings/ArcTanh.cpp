#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"


class cmplx_fxn_ArcTanh : public cmplx_fxn {
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
        std::complex<double> z = std::atanh(x);
        return emake_cmplx(z.real(), z.imag());
    }

    ex eval_double_imag(double y, er def)
    {
        return emake_cmplx(emake_cint(0), emake_double(atan(y)));
    }

    ex eval_arb_imag(const xarb_t&y, er def)
    {
        uex wy(emake_real());
        slong p = y.wprec() + EXTRA_PRECISION_BASIC;
        arb_atan(ereal_data(wy.get()), y.data, p);
        ex t = efix_real(wy.release());
        return emake_cmplx(emake_cint(0), t);
    }

    void acb_eval(acb_t w, const acb_t z, slong prec)
    {
        acb_atanh(w, z, prec);
    }
};


ex num_ArcTanh(er X, er def)
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
        double y = edouble_number(X);
        ex wx, wy;
        if (y > 1)
        {
            wx = emake_double(0.5*log1p(2/(y-1)));
            wy = gs.const_double_mhalfpi.copy();
            return emake_cmplx(wx, wy);
        }
        else if (y < -1)
        {
            wx = emake_double(-0.5*log1p(2/(-y-1)));
            wy = gs.const_double_halfpi.copy();
            return emake_cmplx(wx, wy);
        }
        else
        {
            return emake_double(atanh(y));
        }
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
            /* 1/2 Log[1 + 2/(y-1)] */
            arb_sub_ui(ereal_data(Z), ereal_data(X), 1, p);
            arb_ui_div(ereal_data(Z), 2, ereal_data(Z), p);
            arb_log1p(ereal_data(Z), ereal_data(Z), p);
            arb_mul_2exp_si(ereal_data(Z), ereal_data(Z), -WORD(1));
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), p);
            arb_mul_2exp_si(ereal_data(pi.get()), ereal_data(pi.get()), -WORD(1));
            arb_neg(ereal_data(pi.get()), ereal_data(pi.get()));
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        arb_neg(one.data, one.data);
        if (arb_lt(ereal_data(X), one.data))
        {
            /* -1/2 Log[1 + 2/(-y - 1)] */
            arb_add_ui(ereal_data(Z), ereal_data(X), UWORD(1), p);
            arb_neg(ereal_data(Z), ereal_data(Z));
            arb_ui_div(ereal_data(Z), UWORD(2), ereal_data(Z), p);
            arb_log1p(ereal_data(Z), ereal_data(Z), p);
            arb_mul_2exp_si(ereal_data(Z), ereal_data(Z), -WORD(1));
            arb_neg(ereal_data(Z), ereal_data(Z));
            uex pi(emake_real());
            arb_const_pi(ereal_data(pi.get()), p);
            arb_mul_2exp_si(ereal_data(pi.get()), ereal_data(pi.get()), -WORD(1));
            ex t = efix_real(z.release());
            return emake_cmplx(t, pi.release());
        }
        bool negate = arb_is_negative(ereal_data(X));
            /* -1/2 Log[1 + 2y/(-y - 1)] */
            /* 1/2 Log[1 + 2y/(-y + 1)] */
        if (negate)
            arb_add_ui(ereal_data(Z), ereal_data(X), UWORD(1), p);
        else
            arb_sub_ui(ereal_data(Z), ereal_data(X), UWORD(1), p);
        arb_neg(ereal_data(Z), ereal_data(Z));
        arb_div(ereal_data(Z), ereal_data(X), ereal_data(Z), p);
        arb_mul_2exp_si(ereal_data(Z), ereal_data(Z), WORD(1));
        arb_log1p(ereal_data(Z), ereal_data(Z), p);
        if (negate)
            arb_neg(ereal_data(Z), ereal_data(Z));
        arb_mul_2exp_si(ereal_data(Z), ereal_data(Z), -WORD(1));
        return efix_real(z.release());
    }
    case ETYPE_CMPLX:
    {
        cmplx_fxn_ArcTanh f;
        return f.eval(X, def);
    }
    default:
    {
        return ecopy(def);
    }
    }
}