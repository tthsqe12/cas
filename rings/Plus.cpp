#include "globalstate.h"
#include "eval.h"
#include "code.h"
#include "arithmetic.h"

ex enum_add_int_si(er x, slong y)
{
    ex z = emake_int();
    fmpz_add_si(eint_data(z), eint_data(x), y);
    return efix_int(z);
}

ex enum_add_rat_si(er x, slong y)
{
    ex z = emake_rat();
    fmpq_add_si(erat_data(z), erat_data(x), y);
    return efix_rat(z);
}

ex enum_add_double_si(er x, slong y)
{
	return emake_double(edouble_number(x) + y);
}

ex enum_add_real_si(er x, slong y)
{
    ex z = emake_real();
    slong p = ereal_number(x).wprec();
    arb_add_si(ereal_data(z), ereal_data(x), y, p + EXTRA_PRECISION_BASIC);
	return efix_real(z);
}

ex enum_add_cmplx_si(er x, slong y)
{
    ex z = ex_add_er_si(ecmplx_real(x), y);
    return emake_cmplx(z, ecopy(ecmplx_imag(x)));
}


ex ex_add_er_si(er x, slong y)
{
    uint32_t tx = etype(x);
    switch (tx)
    {
        case ETYPE_INT:
            return enum_add_int_si(x, y);
        case ETYPE_RAT:
            return enum_add_int_si(x, y);
        case ETYPE_DOUBLE:
            return enum_add_double_si(x, y);
        case ETYPE_REAL:
            return enum_add_real_si(x, y);
        case ETYPE_CMPLX:
            return enum_add_cmplx_si(x, y);
        default:
            ex t = emake_int_si(y);
            return ex_addx(ecopy(x), t);
    }
}

ex ex_addx(ex a, ex b)
{
    return ex_canonicalize_plus(emake_node(gs.sym_sPlus.copy(), a, b));
}

ex ex_addr(er a, er b)
{
    return ex_canonicalize_plus(emake_node(gs.sym_sPlus.copy(), ecopy(a), ecopy(b)));
}


ex num_PlusIntInt(er X, er Y)
{
//std::cout << "num_PlusIntInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_int(Y));

    ex Z = emake_int();
    fmpz_add(eint_data(Z), eint_data(X), eint_data(Y));
    return efix_int(Z);
}

ex num_PlusIntRat(er X, er Y)
{
//std::cout << "num_PlusIntRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_rat(Y));

    ex Z = emake_rat();
    fmpq_add_fmpz(erat_data(Z), erat_data(Y), eint_data(X));
    return efix_rat(Z);
}

ex num_PlusRatInt(er Y, er X)
{
//std::cout << "num_PlusRatInt: " << ex_tostring_full(Y) << ", " << ex_tostring_full(X) << std::endl;
    assert(eis_int(X));
    assert(eis_rat(Y));

    ex Z = emake_rat();
    fmpq_add_fmpz(erat_data(Z), erat_data(Y), eint_data(X));
    return efix_rat(Z);
}

ex num_PlusRatRat(er X, er Y)
{
//std::cout << "num_PlusRatRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_rat(Y));

    ex Z = emake_rat();
    fmpq_add(erat_data(Z), erat_data(Y), erat_data(X));
    return efix_rat(Z);
}

ex num_PlusIntDouble(er X, er Y)
{
//std::cout << "num_PlusIntDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_double(Y));

	return emake_double(edouble_number(Y) + fmpz_get_d(eint_data(X)));
}

ex num_PlusRatDouble(er X, er Y)
{
//std::cout << "num_PlusRatDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_double(Y));

	return emake_double(edouble_number(Y) + num_todouble(X));
}

ex num_PlusDoubleDouble(er X, er Y)
{
//std::cout << "num_PlusDoubleDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_double(Y));

	return emake_double(edouble_number(X) + edouble_number(Y));    
}


ex num_PlusIntReal(er X, er Y)
{
//std::cout << "num_PlusIntReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_real(Y));

    slong p = arb_rel_accuracy_bits(eto_real(Y)->number.data);
    xarb_t z;
    arb_add_fmpz(z.data, eto_real(Y)->number.data, eint_data(X), p + EXTRA_PRECISION_BASIC);
    return emake_real_move(z);
}

ex num_PlusRatReal(er X, er Y)
{
//std::cout << "num_PlusRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_real(Y));

    slong p = arb_rel_accuracy_bits(eto_real(Y)->number.data);
    xarb_t z;
    arb_mul_fmpz(z.data, eto_real(Y)->number.data, fmpq_denref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    arb_add_fmpz(z.data, z.data, fmpq_numref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    arb_div_fmpz(z.data, z.data, fmpq_denref(erat_data(X)), p + EXTRA_PRECISION_BASIC);
    return emake_real_move(z);
}

ex num_PlusDoubleReal(er X, er Y)
{
//std::cout << "num_PlusDoubleReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_real(Y));

	return emake_double(edouble_number(X) + num_todouble(Y));    
}   


ex num_PlusRealReal(er X, er Y)
{
//std::cout << "num_PlusRealReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_real(Y));

    slong p = arb_rel_accuracy_bits(eto_real(Y)->number.data);
    p = FLINT_MAX(p, arb_rel_accuracy_bits(eto_real(X)->number.data));
    xarb_t z;
    arb_add(z.data, eto_real(X)->number.data, eto_real(Y)->number.data, p + EXTRA_PRECISION_BASIC);
    return emake_real_move(z);
}


ex num_PlusIntCmplx(er X, er Y)
{
//std::cout << "num_PlusIntComplex: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_cmplx(Y));

	ex re = num_Plus(X, ecmplx_real(Y));
    return emake_cmplx(re, ecopy(ecmplx_imag(Y)));
}

ex num_PlusRatCmplx(er X, er Y)
{
//std::cout << "num_PlusRatCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_cmplx(Y));

	ex re = num_Plus(X, ecmplx_real(Y));
    return emake_cmplx(re, ecopy(ecmplx_imag(Y)));
}

ex num_PlusDoubleCmplx(er X, er Y)
{
//std::cout << "num_PlusDoubleCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_cmplx(Y));

	ex re = num_Plus(X, ecmplx_real(Y));
    return emake_cmplx(re, ecopy(ecmplx_imag(Y)));
}

ex num_PlusRealCmplx(er X, er Y)
{
//std::cout << "num_PlusRealCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_cmplx(Y));

	ex re = num_Plus(X, ecmplx_real(Y));
    return emake_cmplx(re, ecopy(ecmplx_imag(Y)));
}

ex num_PlusCmplxCmplx(er X, er Y)
{
//std::cout << "num_PlusCmplxCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_cmplx(Y));
 
	uex re(num_Plus(ecmplx_real(X), ecmplx_real(Y)));
	ex im = num_Plus(ecmplx_imag(X), ecmplx_imag(Y));
    return emake_cmplx(re.release(), im);
}

ex num_PlusIntNan(er X, er Y)
{
//std::cout << "num_PlusIntNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_nan(Y));
 
    return emake_nan_Indeterminate();
}

ex num_PlusRatNan(er X, er Y)
{
//std::cout << "num_PlusRatNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_PlusDoubleNan(er X, er Y)
{
//std::cout << "num_PlusDoubleNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_PlusRealNan(er X, er Y)
{
//std::cout << "num_PlusRealNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_nan(Y));

	return emake_nan_Indeterminate();
}

ex num_PlusCmplxNan(er X, er Y)
{
//std::cout << "num_PlusCmplxNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}

ex num_PlusNanNan(er X, er Y)
{
//std::cout << "num_PlusNanNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_nan(Y));

    return emake_nan_Indeterminate();
}


ex num_Plus(er X, er Y)
{
//std::cout << "num_Plus2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_PlusIntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_PlusIntRat(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_PlusIntRat(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_PlusRatRat(Y, X);

        case ETYPE_number * ETYPE_INT + ETYPE_DOUBLE:
            return num_PlusIntDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
            return num_PlusIntDouble(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_DOUBLE:
            return num_PlusRatDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return num_PlusRatDouble(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_DOUBLE:
            return num_PlusDoubleDouble(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return num_PlusIntReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return num_PlusIntReal(Y, X);
		case ETYPE_number * ETYPE_RAT + ETYPE_REAL:
            return num_PlusRatReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return num_PlusRatReal(Y, X);
		case ETYPE_number * ETYPE_DOUBLE + ETYPE_REAL:
            return num_PlusDoubleReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_DOUBLE:
            return num_PlusDoubleReal(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_PlusRealReal(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_CMPLX:
            return num_PlusIntCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_INT:
            return num_PlusIntCmplx(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_CMPLX:
            return num_PlusRatCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_RAT:
            return num_PlusRatCmplx(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_CMPLX:
            return num_PlusDoubleCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_DOUBLE:
            return num_PlusDoubleCmplx(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_CMPLX:
            return num_PlusRealCmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_REAL:
            return num_PlusRealCmplx(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_CMPLX:
            return num_PlusCmplxCmplx(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_NAN:
            return num_PlusIntNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_INT:
            return num_PlusIntNan(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_NAN:
            return num_PlusRatNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_RAT:
            return num_PlusRatNan(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_NAN:
            return num_PlusDoubleNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_DOUBLE:
            return num_PlusDoubleNan(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_NAN:
            return num_PlusRealNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_REAL:
            return num_PlusRealNan(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_NAN:
            return num_PlusCmplxNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_CMPLX:
            return num_PlusCmplxNan(Y, X);
        case ETYPE_number * ETYPE_NAN + ETYPE_NAN:
            return num_PlusNanNan(X, Y);

        default:
            assert(false);
            return nullptr;
    }
}

ex ncode_sPlus(er e, slong prec)
{
    if (!ehas_head_sym(e, gs.sym_sPlus.get()))
        return ecopy(e);

    uex f; f.init_push_backr(gs.sym_sPlus.get(), elength(e));
    for (ulong i = 0; i < elength(e); i++)
        f.push_back(eval_num(echild(e,i+1), prec + 1));
    return ex_canonicalize_plus(f.release());
}
