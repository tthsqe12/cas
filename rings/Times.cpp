#include "globalstate.h"
#include "eval.h"
#include "code.h"
#include "arithmetic.h"

ex num_TimesIntInt(er X, er Y)
{
//std::cout << "num_TimesIntInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_int(Y));

    ex z = emake_int();
    fmpz_mul(eint_data(z), eint_data(X), eint_data(Y));
    return efix_int(z);
}

ex num_TimesIntRat(er X, er Y)
{
//std::cout << "num_TimesIntRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_rat(Y));

    ex z = emake_rat();
    fmpq_mul_fmpz(erat_data(z), erat_data(Y), eint_data(X));
    return efix_rat(z);
}

ex num_TimesRatInt(ex Y, ex X)
{
//std::cout << "num_TimesRatInt: " << ex_tostring_full(Y) << ", " << ex_tostring_full(X) << std::endl;
    assert(eis_int(X));
    assert(eis_rat(Y));

    ex z = emake_rat();
    fmpq_mul_fmpz(erat_data(z), erat_data(Y), eint_data(X));
    return efix_rat(z);
}

ex num_TimesRatRat(er X, er Y)
{
//std::cout << "num_TimesRatRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_rat(Y));

    ex z = emake_rat();
    fmpq_mul(erat_data(z), erat_data(Y), erat_data(X));
    return efix_rat(z);
}

ex num_TimesIntDouble(er X, er Y)
{
//std::cout << "num_TimesIntDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_double(Y));

	if (eis_zero(X))
		return ecopy(X);

	return emake_double(fmpz_get_d(eint_data(X)) * edouble_number(Y));
}


/*
double fmpq_get_d(const fmpq_t x)
{
    mpq_t z;
    double d;
    flint_mpq_init_set_readonly(z, x);
    d = mpq_get_d(z);
    flint_mpq_clear_readonly(z);
    return d;
}
*/

ex num_TimesRatDouble(er X, er Y)
{
//std::cout << "num_TimesRatDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_double(Y));

	return emake_double(fmpq_get_d(erat_data(X)) * edouble_number(Y));
}

ex num_TimesDoubleDouble(er X, er Y)
{
//std::cout << "num_TimesDoubleDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_double(Y));

	return emake_double(edouble_number(X) * edouble_number(Y));
}


ex num_TimesIntReal(er X, er Y)
{
//std::cout << "num_TimesIntFlt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_real(Y));

	if (eis_zero(X))
		return ecopy(X);

	ex z = emake_real();
    slong p = arb_rel_accuracy_bits(ereal_data(Y));
    arb_mul_fmpz(ereal_data(z), ereal_data(Y), eint_data(X), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_TimesRatReal(er X, er Y)
{
//std::cout << "num_TimesRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_real(Y));

	ex z = emake_real();
    slong p = arb_rel_accuracy_bits(ereal_data(Y));
	arb_div_fmpz(ereal_data(z), ereal_data(Y), fmpq_denref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    arb_mul_fmpz(ereal_data(z), ereal_data(z), fmpq_numref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_TimesDoubleReal(er X, er Y)
{
//std::cout << "num_TimesDoubleReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_real(Y));

	return emake_double(edouble_number(X) * num_todouble(Y));    
}   


ex num_TimesRealReal(er X, er Y)
{
//std::cout << "num_TimesRealReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_real(Y));

	ex z = emake_real();
    slong p = arb_rel_accuracy_bits(ereal_data(Y));
    p = std::max(p, arb_rel_accuracy_bits(ereal_data(X)));
    arb_mul(ereal_data(z), ereal_data(X), ereal_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}


ex num_TimesIntCmplx(er X, er Y)
{
//std::cout << "num_TimesIntComplex: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_cmplx(Y));
   
	uex re(num_Times(X, ecmplx_real(Y)));
	ex im = num_Times(X, ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_TimesRatCmplx(er X, er Y)
{
//std::cout << "num_TimesRatCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_cmplx(Y));

	uex re(num_Times(X, ecmplx_real(Y)));
	ex im = num_Times(X, ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_TimesDoubleCmplx(er X, er Y)
{
//std::cout << "num_TimesDoubleCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_cmplx(Y));

	uex re(num_Times(X, ecmplx_real(Y)));
	ex im = num_Times(X, ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_TimesRealCmplx(er X, er Y)
{
//std::cout << "num_TimesRealCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_cmplx(Y));

	uex re(num_Times(X, ecmplx_real(Y)));
	ex im = num_Times(X, ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_TimesCmplxCmplx(er X, er Y)
{
//std::cout << "num_TimesCmplxCmplx: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_cmplx(Y));

	uex ad(num_Times(ecmplx_real(X), ecmplx_imag(Y)));
	uex bc(num_Times(ecmplx_imag(X), ecmplx_real(Y)));
	uex ac(num_Times(ecmplx_real(X), ecmplx_real(Y)));
	uex bd(num_Times(ecmplx_imag(X), ecmplx_imag(Y)));
	uex re(num_Minus(ac.get(), bd.get()));
	uex im(num_Plus(ad.get(), bc.get()));
    return emake_cmplx(re.release(), im.release());
}

ex num_TimesIntNan(er X, er Y)
{
//std::cout << "num_TimesIntNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_TimesRatNan(er X, er Y)
{
//std::cout << "num_TimesRatNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_TimesDoubleNan(er X, er Y)
{
//std::cout << "num_TimesDoubleNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_TimesRealNan(er X, er Y)
{
//std::cout << "num_TimesRealNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_TimesCmplxNan(er X, er Y)
{
//std::cout << "num_TimesCmplxNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_TimesNanNan(er X, er Y)
{
//std::cout << "num_TimesNanNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}




ex num_Times(er X, er Y)
{
//std::cout << "num_Times2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_TimesIntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_TimesIntRat(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_TimesIntRat(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_TimesRatRat(Y, X);

        case ETYPE_number * ETYPE_INT + ETYPE_DOUBLE:
            return num_TimesIntDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
            return num_TimesIntDouble(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_DOUBLE:
            return num_TimesRatDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return num_TimesRatDouble(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_DOUBLE:
            return num_TimesDoubleDouble(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return num_TimesIntReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return num_TimesIntReal(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_REAL:
            return num_TimesRatReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return num_TimesRatReal(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_REAL:
            return num_TimesDoubleReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_DOUBLE:
            return num_TimesDoubleReal(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_TimesRealReal(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_CMPLX:
            return num_TimesIntCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_INT:
            return num_TimesIntCmplx(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_CMPLX:
            return num_TimesRatCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_RAT:
            return num_TimesRatCmplx(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_CMPLX:
            return num_TimesDoubleCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_DOUBLE:
            return num_TimesDoubleCmplx(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_CMPLX:
            return num_TimesRealCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_REAL:
            return num_TimesRealCmplx(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_CMPLX:
            return num_TimesCmplxCmplx(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_NAN:
            return num_TimesIntNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_INT:
            return num_TimesIntNan(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_NAN:
            return num_TimesRatNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_RAT:
            return num_TimesRatNan(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_NAN:
            return num_TimesDoubleNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_DOUBLE:
            return num_TimesDoubleNan(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_NAN:
            return num_TimesRealNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_REAL:
            return num_TimesRealNan(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_NAN:
            return num_TimesCmplxNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_CMPLX:
            return num_TimesCmplxNan(Y, X);
        case ETYPE_number * ETYPE_NAN + ETYPE_NAN:
            return num_TimesNanNan(X, Y);

        default:
            assert(false);
            return nullptr;
    }
}


ex ncode_sTimes(er e, slong prec)
{
    if (!ehas_head_sym(e, gs.sym_sTimes.get()))
        return ecopy(e);

    size_t n = elength(e);
    prec += FLINT_BIT_COUNT(n);
    uex f; f.init_push_backr(gs.sym_sTimes.get(), n);
    for (size_t i = 0; i < n; i++)
        f.push_back(eval_num(echild(e,i+1), prec));
    return ex_canonicalize_times(f.release());
}
