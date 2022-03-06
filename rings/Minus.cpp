#include "globalstate.h"
#include "arithmetic.h"
#include "ex_print.h"
#include "timing.h"



ex num_MinusInt(er X)
{
//std::cout << "num_MinusInt: " << ex_tostring_full(X) << std::endl;
    assert(eis_int(X));

    ex Z = emake_int();
    fmpz_neg(eint_data(Z), eint_data(X));
    return efix_int(Z);
}

ex num_MinusRat(er X)
{
//std::cout << "num_MinusRat: " << ex_tostring_full(X) << std::endl;
    assert(eis_rat(X));

    ex Z = emake_rat();
    fmpq_neg(erat_data(Z), erat_data(X));
    return efix_rat(Z);
}

ex num_MinusDouble(er X)
{
//std::cout << "num_MinusDouble: " << ex_tostring_full(X) << std::endl;
    assert(eis_double(X));

    return emake_double(-eto_double(X)->number);
}

ex num_MinusReal(er X)
{
//std::cout << "num_MinusReal: " << ex_tostring_full(X) << std::endl;
    assert(eis_real(X));

    xarb_t z;
    arb_neg(z.data, eto_real(X)->number.data);
    return emake_real_move(z);        
}

ex num_MinusCmplx(er X)
{
//std::cout << "num_MinusCmplx: " << ex_tostring_full(X) << std::endl;
    assert(eis_cmplx(X));

	uex re(num_Minus(ecmplx_real(X)));
	ex im = num_Minus(ecmplx_imag(X));
    return emake_cmplx(re.release(), im);
}

ex num_Minus(er X)
{
//std::cout << "num_Minus1: " << ex_tostring_full(X) << std::endl;
    assert(eis_number(X));

    uint32_t tx = etype(X);
    switch (tx)
    {
        case ETYPE_INT:
            return num_MinusInt(X);
        case ETYPE_RAT:
            return num_MinusRat(X);
        case ETYPE_DOUBLE:
            return num_MinusDouble(X);
        case ETYPE_REAL:
            return num_MinusReal(X);
        case ETYPE_CMPLX:
            return num_MinusCmplx(X);
        default:
            assert(false);
            return nullptr;
    }
}


ex num_MinusIntInt(er X, er Y)
{
//std::cout << "num_MinusIntInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_int(Y));

    ex z = emake_int();
    fmpz_sub(eint_data(z), eint_data(X), eint_data(Y));
    return efix_int(z);
}

ex num_MinusRatInt(er X, er Y)
{
//std::cout << "num_MinusIntRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(Y));
    assert(eis_rat(X));

    ex z = emake_rat();
    fmpq_sub_fmpz(erat_data(z), erat_data(X), eint_data(Y));
    return efix_rat(z);
}

ex num_MinusIntRat(er X, er Y)
{
//std::cout << "num_MinusIntRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_rat(Y));

    ex z = emake_rat();
    fmpq_sub_fmpz(erat_data(z), erat_data(Y), eint_data(X));
	fmpz_neg(fmpq_numref(erat_data(z)), fmpq_numref(erat_data(z)));
    return efix_rat(z);
}


ex num_MinusRatRat(er X, er Y)
{
//std::cout << "num_MinusRatRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_rat(Y));

    ex z = emake_rat();
    fmpq_sub(erat_data(z), erat_data(Y), erat_data(X));
    return efix_rat(z);
}

ex num_MinusIntDouble(er X, er Y)
{
//std::cout << "num_MinusIntDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_double(Y));

	return emake_double(fmpz_get_d(eint_data(X)) - edouble_number(Y));    
}

ex num_MinusDoubleInt(er X, er Y)
{
//std::cout << "num_MinusIntDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(Y));
    assert(eis_double(X));

	return emake_double(edouble_number(X) - fmpz_get_d(eint_data(Y)));    
}

ex num_MinusRatDouble(er X, er Y)
{
//std::cout << "num_MinusRatDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_double(Y));

	return emake_double(num_todouble(X) - edouble_number(Y));
}

ex num_MinusDoubleRat(er X, er Y)
{
//std::cout << "num_MinusRatDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(Y));
    assert(eis_double(X));

	return emake_double(edouble_number(X) - num_todouble(Y));
}

ex num_MinusDoubleDouble(er X, er Y)
{
//std::cout << "num_MinusDoubleDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_double(Y));

	return emake_double(edouble_number(X) - edouble_number(Y));    
}


ex num_MinusIntReal(er X, er Y)
{
//std::cout << "num_MinusIntReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_real(Y));

    ex z = emake_real();
    slong p = arb_rel_accuracy_bits(ereal_data(Y));
    arb_sub_fmpz(ereal_data(z), ereal_data(Y), eint_data(X), p + EXTRA_PRECISION_BASIC);
    arb_neg(ereal_data(z), ereal_data(z));
    return efix_real(z);
}

ex num_MinusRealInt(er X, er Y)
{
//std::cout << "num_MinusIntFlt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_int(Y));

    ex z = emake_real();
    slong p = arb_rel_accuracy_bits(ereal_data(X));
    arb_sub_fmpz(ereal_data(z), ereal_data(X), eint_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_MinusRatReal(er X, er Y)
{
//std::cout << "num_MinusRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_real(Y));

    ex z = emake_real();
    slong p = arb_rel_accuracy_bits(ereal_data(Y));
    arb_mul_fmpz(ereal_data(z), eto_real(Y)->number.data, fmpq_denref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    arb_sub_fmpz(gs.tmpreal[0].data, ereal_data(z), fmpq_numref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    arb_div_fmpz(ereal_data(z), gs.tmpreal[0].data, fmpq_denref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
	arb_neg(ereal_data(z), ereal_data(z));
    return efix_real(z);
}

ex num_MinusRealRat(er X, er Y)
{
//std::cout << "num_MinusRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_rat(Y));

    ex z = emake_real();
    slong p = arb_rel_accuracy_bits(eto_real(X)->number.data);
    arb_mul_fmpz(ereal_data(z), eto_real(X)->number.data, fmpq_denref(erat_data(Y)), p + EXTRA_PRECISION_BASIC);
    arb_sub_fmpz(gs.tmpreal[0].data, ereal_data(z), fmpq_numref(erat_data(Y)), p + EXTRA_PRECISION_BASIC);
    arb_div_fmpz(ereal_data(z), gs.tmpreal[0].data, fmpq_denref(erat_data(Y)), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_MinusDoubleReal(er X, er Y)
{
//std::cout << "num_MinusDoubleReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_real(Y));

	return emake_double(edouble_number(X) - num_todouble(Y));
}   

ex num_MinusRealDouble(er X, er Y)
{
//std::cout << "num_MinusDoubleReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(Y));
    assert(eis_real(X));

	return emake_double(num_todouble(X) - edouble_number(Y));
}

ex num_MinusRealReal(er X, er Y)
{
//std::cout << "num_MinusRealReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_real(Y));

	ex z = emake_real();
    slong p = arb_rel_accuracy_bits(ereal_data(Y));
    p = FLINT_MAX(p, arb_rel_accuracy_bits(ereal_data(X)));
    arb_sub(ereal_data(z), ereal_data(X), ereal_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}


ex num_MinusIntCmplx(er X, er Y)
{
//std::cout << "num_MinusIntComplex: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_cmplx(Y));

	uex re(num_Minus(X, ecmplx_real(Y)));
	ex im = num_Minus(ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_MinusCmplxInt(er X, er Y)
{
//std::cout << "num_MinusIntComplex: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_int(Y));
 
	ex re = num_Minus(ecmplx_real(X), Y);
    return emake_cmplx(re, ecopy(ecmplx_imag(X)));
}

ex num_MinusRatCmplx(er X, er Y)
{
//std::cout << "num_MinusRatCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(Y));
    assert(eis_rat(X));
  
	uex re(num_Minus(X, ecmplx_real(Y)));
	ex im = num_Minus(ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);

}

ex num_MinusCmplxRat(er X, er Y)
{
//std::cout << "num_MinusIntComplex: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_rat(Y));
   
	ex re = num_Minus(ecmplx_real(X), Y);
    return emake_cmplx(re, ecopy(ecmplx_imag(X)));
}

ex num_MinusDoubleCmplx(er X, er Y)
{
//std::cout << "num_MinusDoubleCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_cmplx(Y));
   
	uex re(num_Minus(X, ecmplx_real(Y)));
	ex im = num_Minus(ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_MinusCmplxDouble(er X, er Y)
{
//std::cout << "num_MinusDoubleCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_double(Y));

	ex re = num_Minus(ecmplx_real(X), Y);
    return emake_cmplx(re, ecopy(ecmplx_imag(X)));
}

ex num_MinusRealCmplx(er X, er Y)
{
//std::cout << "num_MinusRealCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_cmplx(Y));
   
	uex re(num_Minus(X, ecmplx_real(Y)));
	ex im = num_Minus(ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_MinusCmplxReal(er X, er Y)
{
//std::cout << "num_MinusCmplxReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_real(Y));

	ex re = num_Minus(ecmplx_real(X), Y);
    return emake_cmplx(re, ecopy(ecmplx_imag(X)));
}

ex num_MinusCmplxCmplx(er X, er Y)
{
//std::cout << "num_MinusCmplxCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_cmplx(Y));

	uex re(num_Minus(ecmplx_real(X), ecmplx_real(Y)));
	ex im = num_Minus(ecmplx_imag(X), ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_MinusIntNan(er X, er Y)
{
//std::cout << "num_MinusIntNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_MinusNanInt(er X, er Y)
{
//std::cout << "num_MinusIntNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_int(Y));

    return emake_nan_Indeterminate();
}

ex num_MinusRatNan(er X, er Y)
{
//std::cout << "num_MinusRatNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_MinusNanRat(er X, er Y)
{
//std::cout << "num_MinusRatNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_rat(Y));

    return emake_nan_Indeterminate();
}

ex num_MinusDoubleNan(er X, er Y)
{
//std::cout << "num_MinusDoubleNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_MinusNanDouble(er X, er Y)
{
//std::cout << "num_MinusDoubleNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_double(Y));

	return emake_nan_Indeterminate();
}

ex num_MinusRealNan(er X, er Y)
{
//std::cout << "num_MinusRealNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_MinusNanReal(er X, er Y)
{
//std::cout << "num_MinusRealNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_real(Y));

	return emake_nan_Indeterminate();
}

ex num_MinusCmplxNan(er X, er Y)
{
//std::cout << "num_MinusCmplxNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_MinusNanCmplx(er X, er Y)
{
//std::cout << "num_MinusCmplxNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_cmplx(Y));

    return emake_nan_Indeterminate();
}

ex num_MinusNanNan(er X, er Y)
{
//std::cout << "num_MinusNanNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}


ex num_Minus(er X, er Y)
{
//std::cout << "num_Minus2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_MinusIntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_MinusIntRat(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_MinusRatInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_MinusRatRat(Y, X);

        case ETYPE_number * ETYPE_INT + ETYPE_DOUBLE:
            return num_MinusIntDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
            return num_MinusDoubleInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_DOUBLE:
            return num_MinusRatDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return num_MinusDoubleRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_DOUBLE:
            return num_MinusDoubleDouble(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return num_MinusIntReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return num_MinusRealInt(X, Y);
		case ETYPE_number * ETYPE_RAT + ETYPE_REAL:
            return num_MinusRatReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return num_MinusRealRat(X, Y);
		case ETYPE_number * ETYPE_DOUBLE + ETYPE_REAL:
            return num_MinusDoubleReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_DOUBLE:
            return num_MinusRealDouble(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_MinusRealReal(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_CMPLX:
            return num_MinusIntCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_INT:
            return num_MinusCmplxInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_CMPLX:
            return num_MinusRatCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_RAT:
            return num_MinusCmplxRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_CMPLX:
            return num_MinusDoubleCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_DOUBLE:
            return num_MinusCmplxDouble(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_CMPLX:
            return num_MinusRealCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_REAL:
            return num_MinusCmplxReal(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_CMPLX:
            return num_MinusCmplxCmplx(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_NAN:
            return num_MinusIntNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_INT:
            return num_MinusNanInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_NAN:
            return num_MinusRatNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_RAT:
            return num_MinusNanRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_NAN:
            return num_MinusDoubleNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_DOUBLE:
            return num_MinusNanDouble(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_NAN:
            return num_MinusRealNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_REAL:
            return num_MinusNanReal(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_NAN:
            return num_MinusCmplxNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_CMPLX:
            return num_MinusNanCmplx(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_NAN:
            return num_MinusNanNan(X, Y);

        default:
            assert(false);
            return nullptr;
    }
}
