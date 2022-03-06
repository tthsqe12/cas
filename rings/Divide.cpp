#include "globalstate.h"
#include "arithmetic.h"
#include "ex.h"

ex num_DivideInt(er X)
{
//std::cout << "num_DivideInt: " << ex_tostring_full(X) << std::endl;
    assert(eis_int(X));

    if (fmpz_is_pm1(eint_data(X)))
    {
        return ecopy(X);
    }
    else if (fmpz_is_zero(eint_data(X)))
    {
        return gs.const_complexinfinity.copy();
    }
    else
    {
        ex z = emake_rat();
        if (fmpz_sgn(eint_data(X)) < 0)
        {
            fmpz_set_si(fmpq_numref(erat_data(z)), -1);
            fmpz_neg(fmpq_denref(erat_data(z)), eint_data(X));
        }
        else
        {
            fmpz_one(fmpq_numref(erat_data(z)));
            fmpz_set(fmpq_denref(erat_data(z)), eint_data(X));
        }
        return efix_rat(z);
    }
}

ex num_DivideRat(er X)
{
//std::cout << "num_DivideRat: " << ex_tostring_full(X) << std::endl;
    assert(eis_rat(X));

    xfmpq_t z;
    if (fmpz_sgn(fmpq_numref(erat_data(X))) < 0)
    {
        fmpz_neg(fmpq_denref(z.data), fmpq_numref(erat_data(X)));
        fmpz_neg(fmpq_numref(z.data), fmpq_denref(erat_data(X)));
    }
    else
    {
        fmpz_set(fmpq_denref(z.data), fmpq_numref(erat_data(X)));
        fmpz_set(fmpq_numref(z.data), fmpq_denref(erat_data(X)));
    }
    return emake_rat_move(z.data);      
}


ex num_DivideReal(er X)
{
//std::cout << "num_DivideReal: " << ex_tostring_full(X) << std::endl;
    assert(eis_real(X));

	ex z = emake_real();
    slong p = ereal_number(X).wprec();
    arb_ui_div(ereal_data(z), 1, ereal_data(X), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_DivideCmplx(er Y)
{
//std::cout << "num_DivideCmplx: " << ex_tostring(Y) << std::endl;
    assert(eis_cmplx(Y));

	uex yre2(num_Times(ecmplx_real(Y), ecmplx_real(Y)));
	uex yim2(num_Times(ecmplx_imag(Y), ecmplx_imag(Y)));
	uex yabs2(num_Plus(yre2.get(), yim2.get()));
	uex re(num_Divide(ecmplx_real(Y), yabs2.get()));
	uex im(num_Divide(ecmplx_imag(Y), yabs2.get()));
	ex t = num_Minus(im.get());
    return emake_cmplx(re.release(), t);
}

ex num_Divide(er X)
{
//std::cout << "num_Divide: " << ex_tostring(X) << std::endl;
    assert(eis_number(X));

    uint32_t tx = etype(X);
    switch (tx)
    {
        case ETYPE_INT:
            return num_DivideInt(X);
        case ETYPE_RAT:
            return num_DivideRat(X);
        case ETYPE_REAL:
            return num_DivideReal(X);
        case ETYPE_CMPLX:
            return num_DivideCmplx(X);
        default:
            assert(false);
            return nullptr;
    }
}



ex num_DivideIntInt(er X, er Y)
{
//std::cout << "num_DivideIntInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_int(Y));

    if (fmpz_is_zero(eint_data(Y)))
        return fmpz_is_zero(eint_data(X)) ? gs.const_indeterminate.copy()
										  : gs.const_complexinfinity.copy();

    if (fmpz_is_one(eint_data(Y)) || fmpz_is_zero(eint_data(X)))
        return ecopy(X);

	ex z = emake_rat();

	fmpz_t u;

    fmpz_init(u);
    fmpz_gcd(u, eint_data(X), eint_data(Y));
    if (fmpz_sgn(eint_data(Y)) < 0)
        fmpz_neg(u, u);

    if (!fmpz_is_one(u))
    {
        fmpz_divexact(fmpq_numref(erat_data(z)), eint_data(X), u);
        fmpz_divexact(fmpq_denref(erat_data(z)), eint_data(Y), u);
    }
	else
	{
        fmpz_set(fmpq_numref(erat_data(z)), eint_data(X));
        fmpz_set(fmpq_denref(erat_data(z)), eint_data(Y));
	}
    fmpz_clear(u);

    return efix_rat(z);
}

ex num_DivideIntRat(er X, er Y)
{
//std::cout << "num_DivideIntRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_rat(Y));

    xfmpq_t z;
    if (fmpz_sgn(fmpq_denref(erat_data(Y))) < 0)
    {
        fmpz_neg(fmpq_denref(z.data), fmpq_numref(erat_data(Y)));
        fmpz_neg(fmpq_numref(z.data), fmpq_denref(erat_data(Y)));
    }
    else
    {
        fmpz_set(fmpq_denref(z.data), fmpq_numref(erat_data(Y)));
        fmpz_set(fmpq_numref(z.data), fmpq_denref(erat_data(Y)));
    }
    fmpq_mul_fmpz(z.data, z.data, eint_data(X));
    return emake_rat_move(z.data);
}

ex num_DivideRatInt(er X, er Y)
{
//std::cout << "num_DivideRatInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_int(Y));

    if (fmpz_is_zero(eint_data(Y)))
        return gs.const_complexinfinity.copy();

    ex z = emake_rat();
    fmpq_div_fmpz(erat_data(z), erat_data(X), eint_data(Y));
    return efix_rat(z);
}

ex num_DivideRatRat(er X, er Y)
{
//std::cout << "num_DivideRatRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_rat(Y));

    ex z = emake_rat();
    fmpq_div(erat_data(z), erat_data(X), erat_data(Y));
    return efix_rat(z);
}

ex num_DivideIntDouble(er X, er Y)
{
//std::cout << "num_DivideIntDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_double(Y));

	return emake_double(fmpz_get_d(eint_data(X))/edouble_number(Y));    
}

ex num_DivideDoubleInt(er X, er Y)
{
//std::cout << "num_DivideDoubleInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(Y));
    assert(eis_double(X));

	return emake_double(edouble_number(X)/fmpz_get_d(eint_data(Y)));    
}

ex num_DivideRatDouble(er X, er Y)
{
//std::cout << "num_DivideRatDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_double(Y));

	return emake_double(num_todouble(X)/edouble_number(Y));
}

ex num_DivideDoubleRat(er X, er Y)
{
//std::cout << "num_DivideDoubleRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_rat(Y));

	return emake_double(edouble_number(X)/fmpq_get_d(erat_data(Y)));
}

ex num_DivideDoubleDouble(er X, er Y)
{
//std::cout << "num_DivideDoubleDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_double(Y));

	return emake_double(edouble_number(X)/edouble_number(Y));    
}

ex num_DivideIntReal(er X, er Y)
{
//std::cout << "num_MinusIntFlt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_real(Y));

	ex z = emake_real();
    slong p = ereal_number(Y).wprec();
	arb_set_round_fmpz(gs.tmpreal[0].data, eint_data(X), p + EXTRA_PRECISION_BASIC);
    arb_div(ereal_data(z), gs.tmpreal[0].data, ereal_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_DivideRealInt(er X, er Y)
{
//std::cout << "num_MinusIntFlt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_int(Y));

	ex z = emake_real();
    slong p = ereal_number(X).wprec();
    arb_div_fmpz(ereal_data(z), ereal_data(X), eint_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_DivideRatReal(er X, er Y)
{
//std::cout << "num_MinusRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_real(Y));

	ex z = emake_real();
    slong p = ereal_number(Y).wprec();
	arb_mul_fmpz(gs.tmpreal[0].data, ereal_data(Y), fmpq_denref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
	arb_set_round_fmpz(gs.tmpreal[1].data, fmpq_numref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    arb_div(ereal_data(z), gs.tmpreal[1].data, gs.tmpreal[0].data, p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_DivideRealRat(er X, er Y)
{
//std::cout << "num_MinusRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_rat(Y));

	ex z = emake_real();
    slong p = ereal_number(X).wprec();
    arb_mul_fmpz(gs.tmpreal[0].data, ereal_data(X), fmpq_denref(erat_data(Y)), p + EXTRA_PRECISION_BASIC);
    arb_div_fmpz(ereal_data(z), gs.tmpreal[0].data, fmpq_numref(erat_data(Y)), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_DivideDoubleReal(er X, er Y)
{
//std::cout << "num_MinusDoubleReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_real(Y));

	return emake_double(edouble_number(X)/num_todouble(Y));    
}   

ex num_DivideRealDouble(er X, er Y)
{
//std::cout << "num_DivideRealDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(Y));
    assert(eis_real(X));

	return emake_double(num_todouble(X)/edouble_number(Y));    
}

ex num_DivideRealReal(er X, er Y)
{
//std::cout << "num_DivideRealReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_real(Y));

    slong p = ereal_number(Y).wprec();
    p = FLINT_MAX(p, ereal_number(X).wprec());
    xarb_t z;
    arb_div(z.data, ereal_data(X), ereal_data(Y), p + EXTRA_PRECISION_BASIC);
    return emake_real_move(z);
}


ex num_DivideIntCmplx(er X, er Y)
{
//std::cout << "num_DivideIntCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_cmplx(Y));

	uex yre2(num_Times(ecmplx_real(Y), ecmplx_real(Y)));
	uex yim2(num_Times(ecmplx_imag(Y), ecmplx_imag(Y)));
	uex yabs2(num_Plus(yre2.get(), yim2.get()));
	uex s(num_Divide(X, yabs2.get()));
	uex re(num_Times(s.get(), ecmplx_real(Y)));
	uex im(num_Times(s.get(), ecmplx_imag(Y)));
	ex t = num_Minus(im.get());
    return emake_cmplx(re.release(), t);
}

ex num_DivideCmplxInt(er X, er Y)
{
//std::cout << "num_DivideCmplxInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_int(Y));

	uex re(num_Divide(ecmplx_real(X), Y));
	uex im(num_Divide(ecmplx_imag(X), Y));
    return emake_cmplx(re.release(), im.release());
}

ex num_DivideRatCmplx(er X, er Y)
{
//std::cout << "num_MinusRatCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_cmplx(Y));

	uex yre2(num_Times(ecmplx_real(Y), ecmplx_real(Y)));
	uex yim2(num_Times(ecmplx_imag(Y), ecmplx_imag(Y)));
	uex yabs2(num_Plus(yre2.get(), yim2.get()));
	uex s(num_Divide(X, yabs2.get()));
	uex re(num_Times(s.get(), ecmplx_real(Y)));
	uex im(num_Times(s.get(), ecmplx_imag(Y)));
	ex t = num_Minus(im.get());
    return emake_cmplx(re.release(), t);
}

ex num_DivideCmplxRat(er X, er Y)
{
//std::cout << "num_MinusIntComplex: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_rat(Y));

	uex re(num_Divide(ecmplx_real(X), Y));
	ex im = num_Divide(ecmplx_imag(X), Y);
    return emake_cmplx(re.release(), im);
}

ex num_DivideDoubleCmplx(er X, er Y)
{
//std::cout << "num_MinusDoubleCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_cmplx(Y));

	uex yre2(num_Times(ecmplx_real(Y), ecmplx_real(Y)));
	uex yim2(num_Times(ecmplx_imag(Y), ecmplx_imag(Y)));
	uex yabs2(num_Plus(yre2.get(), yim2.get()));
	uex s(num_Divide(X, yabs2.get()));
	uex re(num_Times(s.get(), ecmplx_real(Y)));
	uex im(num_Times(s.get(), ecmplx_imag(Y)));
	ex t = num_Minus(im.get());
    return emake_cmplx(re.release(), t);
}

ex num_DivideCmplxDouble(er X, er Y)
{
//std::cout << "num_MinusDoubleCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_double(Y));

    uex re(num_Divide(ecmplx_real(X), Y));
    ex im = num_Divide(ecmplx_imag(X), Y);
    return emake_cmplx(re.release(), im);
}

ex num_DivideRealCmplx(er X, er Y)
{
//std::cout << "num_MinusRealCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_cmplx(Y));

	uex re(num_Minus(X, ecmplx_real(Y)));
	ex im = num_Minus(ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_DivideCmplxReal(er X, er Y)
{
//std::cout << "num_MinusCmplxReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_real(Y));

	uex re(num_Divide(ecmplx_real(X), Y));
	ex im = num_Divide(ecmplx_imag(X), Y);
    return emake_cmplx(re.release(), im);
}

ex num_DivideCmplxCmplx(er X, er Y)
{
//std::cout << "num_MinusCmplxCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_cmplx(Y));

	uex z(num_Divide(Y));
	return num_Times(X, z.get());
}

ex num_DivideIntNan(er X, er Y)
{
//std::cout << "num_DivideIntNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_DivideNanInt(er X, er Y)
{
//std::cout << "num_DivideNanInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_int(Y));

    return emake_nan_Indeterminate();
}

ex num_DivideRatNan(er X, er Y)
{
//std::cout << "num_DivideRatNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_DivideNanRat(er X, er Y)
{
//std::cout << "num_DivideNanRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_rat(Y));

    return emake_nan_Indeterminate();
}

ex num_DivideDoubleNan(er X, er Y)
{
//std::cout << "num_DivideDoubleNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_DivideNanDouble(er X, er Y)
{
//std::cout << "num_DivideNanDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_double(Y));

	return emake_nan_Indeterminate();
}

ex num_DivideRealNan(er X, er Y)
{
//std::cout << "num_DivideRealNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_DivideNanReal(er X, er Y)
{
//std::cout << "num_DivideNanReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_real(Y));

	return emake_nan_Indeterminate();
}

ex num_DivideCmplxNan(er X, er Y)
{
//std::cout << "num_DivideCmplxNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_DivideNanCmplx(er X, er Y)
{
//std::cout << "num_DivideNanCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_cmplx(Y));

    return emake_nan_Indeterminate();
}

ex num_DivideNanNan(er X, er Y)
{
//std::cout << "num_DivideNanNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_Divide(er X, er Y)
{
//std::cout << "num_Minus2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_DivideIntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_DivideIntRat(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_DivideRatInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_DivideRatRat(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_DOUBLE:
            return num_DivideIntDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
            return num_DivideDoubleInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_DOUBLE:
            return num_DivideRatDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return num_DivideDoubleRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_DOUBLE:
            return num_DivideDoubleDouble(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return num_DivideIntReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return num_DivideRealInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_REAL:
            return num_DivideRatReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return num_DivideRealRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_REAL:
            return num_DivideDoubleReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_DOUBLE:
            return num_DivideRealDouble(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_DivideRealReal(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_CMPLX:
            return num_DivideIntCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_INT:
            return num_DivideCmplxInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_CMPLX:
            return num_DivideRatCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_RAT:
            return num_DivideCmplxRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_CMPLX:
            return num_DivideDoubleCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_DOUBLE:
            return num_DivideCmplxDouble(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_CMPLX:
            return num_DivideRealCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_REAL:
            return num_DivideCmplxReal(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_CMPLX:
            return num_DivideCmplxCmplx(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_NAN:
            return num_DivideIntNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_INT:
            return num_DivideNanInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_NAN:
            return num_DivideRatNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_RAT:
            return num_DivideNanRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_NAN:
            return num_DivideDoubleNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_DOUBLE:
            return num_DivideNanDouble(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_NAN:
            return num_DivideRealNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_REAL:
            return num_DivideNanReal(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_NAN:
            return num_DivideCmplxNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_CMPLX:
            return num_DivideNanCmplx(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_NAN:
            return num_DivideNanNan(X, Y);

        default:
            assert(false);
            return nullptr;
    }
}
