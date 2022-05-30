#include "globalstate.h"
#include "eval.h"
#include "code.h"
#include "arithmetic.h"
#include "flintarb_wrappers.h"


// (n/2)! for n odd
ex _factorial_half_int(slong n)
{
    uex sqrtpi(emake_node(gs.sym_sPower.copy(), gs.sym_sPi.copy(), emake_crat(1,2)));

    if (n == -1)
        return sqrtpi.release();

    xfmpq_t z(1,1);        
    if (n > 0)
    {
        fmpz_mul_2exp(fmpq_denref(z.data), fmpq_denref(z.data), (n + 1)/2);
        for (slong i = 3; i <= n; i += 2)
            fmpz_mul_ui(fmpq_numref(z.data), fmpq_numref(z.data), i);
    }
    else
    {
        n = -n;
        fmpz_mul_2exp(fmpq_numref(z.data), fmpq_numref(z.data), (n - 1)/2);
        if (n & 2)
            fmpz_neg(fmpq_numref(z.data), fmpq_numref(z.data));
        for (slong i = 3; i < n; i += 2)
            fmpz_mul_ui(fmpq_denref(z.data), fmpq_denref(z.data), i);
    }
    ex t = emake_rat_move(z.data);
    return emake_node(gs.sym_sTimes.copy(), t, sqrtpi.release());
}


ex dcode_sGamma(er e)
{
//std::cout << "dcode_sGamma: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sGamma.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_int(X))
    {
        if (fmpz_sgn(eint_data(X)) <= 0)
        {
            return emake_nan_ComplexInfinity();
        }
        else if (fmpz_fits_si(eint_data(X)))
        {
            ex z = emake_int();
            fmpz_fac_ui(eint_data(z), fmpz_get_ui(eint_data(X)) - 1);
            return efix_int(z);
        }
    }
    else if (eis_rat(X))
    {
        if (fmpz_equal_int(fmpq_denref(erat_data(X)), 2))
        {
            slong n = *fmpq_numref(erat_data(X));
            if (-10000 < n && n < 10000 && (n & 1))
                return _factorial_half_int(n - 2);
        }
    }
    else if (eis_real(X))
    {
        xarb_t z;
        slong p = ereal_number(X).wprec();
        arb_gamma(z.data, ereal_number(X).data, p + EXTRA_PRECISION_BASIC);
        return emake_real_move(z);
    }
    else if (eis_double(X))
    {
        xarb_t z;
        slong p = 53;
        arb_set_d(z.data, edouble_number(X));
        arb_gamma(z.data, z.data, p + EXTRA_PRECISION_BASIC);
        return emake_double(arf_get_d(arb_midref(z.data), ARF_RND_NEAR));
    }

    return ecopy(e);
}


ex dcode_sFactorial(er e)
{
//std::cout << "dcode_sFactorial: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFactorial.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_int(X))
    {
        if (fmpz_sgn(eint_data(X)) < 0)
        {
            return emake_nan_ComplexInfinity();
        }
        else if (fmpz_fits_si(eint_data(X)))
        {
            ex z = emake_int();
            fmpz_fac_ui(eint_data(z), fmpz_get_ui(eint_data(X)));
            return efix_int(z);
        }
    }
    else if (eis_rat(X))
    {
        if (fmpz_equal_int(fmpq_denref(erat_data(X)), 2))
        {
            slong n = *fmpq_numref(erat_data(X));
            if (-10000 < n && n < 10000 && (n & 1))
                return _factorial_half_int(n);
        }
    }
    else if (eis_real(X))
    {
        xarb_t z;
        slong p = ereal_number(X).wprec();
        arb_add_ui(z.data, ereal_number(X).data, 1, p + EXTRA_PRECISION_BASIC);
        arb_gamma(z.data, z.data, p + EXTRA_PRECISION_BASIC);
        return emake_real_move(z);
    }
    else if (eis_double(X))
    {
        xarb_t z;
        slong p = 53;
        arb_set_d(z.data, edouble_number(X));
        arb_add_ui(z.data, z.data, 1, p + EXTRA_PRECISION_BASIC);
        arb_gamma(z.data, z.data, p + EXTRA_PRECISION_BASIC);
        return emake_double(arf_get_d(arb_midref(z.data), ARF_RND_NEAR));
    }

    return ecopy(e);
}


ex ncode_sGamma(er e, slong prec)
{
//std::cout << "ncode_sGamma(" << prec << "): " << e << std::endl;
    if (!ehas_head_sym_length(e, gs.sym_sGamma.get(), 1))
        return ecopy(e);

    er e1 = echild(e,1);

    if (eis_int(e1) || eis_rat(e1))
    {
        uex z(emake_real());
        if (eis_int(e1))
            arb_gamma_fmpz(ereal_data(z.get()), eint_data(e1), prec);
        else
            arb_gamma_fmpq(ereal_data(z.get()), erat_data(e1), prec);
        return efix_real(z.release());
    }

    return ecopy(e);
}

