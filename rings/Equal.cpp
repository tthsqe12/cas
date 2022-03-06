#include "globalstate.h"
#include "arithmetic.h"
#include "ex_print.h"
#include "timing.h"


//1.0 == 1.0 + (2^6)*$MachineEpsilon
//1.0 != 1.0 + (2^6+1)*$MachineEpsilon
static bool Equal_to_7bits(double x, double y)
{
    double tol = 0.9999999999999858;
    if (x > 0)
        return y > 0 && (x < y ? x >= tol*y : y >= tol*x);
    else if (x < 0)
        return y < 0 && (x > y ? x <= tol*y : y <= tol*x);
    else
        return x == 0 && y == 0;
}

ex num_EqualIntInt(er X, er Y)
{
//std::cout << "num_EqualIntInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_int(Y));

    return emake_boole(fmpz_equal(eint_data(X), eint_data(Y)));
}

ex num_EqualIntRat(er X, er Y)
{
//std::cout << "num_EqualIntRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_rat(Y));

    return gs.sym_sFalse.copy();
}


ex num_EqualRatRat(er X, er Y)
{
//std::cout << "num_EqualRatRat: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_rat(Y));

    return emake_boole(fmpq_equal(erat_data(X), erat_data(Y)));
}

ex num_EqualIntDouble(er X, er Y)
{
//std::cout << "num_EqualIntDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_double(Y));

	return emake_boole(Equal_to_7bits(fmpz_get_d(eint_data(X)), edouble_number(Y)));
}


ex num_EqualRatDouble(er X, er Y)
{
//std::cout << "num_EqualRatDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_double(Y));

	return emake_boole(Equal_to_7bits(fmpq_get_d(erat_data(X)), edouble_number(Y)));
}


ex num_EqualDoubleDouble(er X, er Y)
{
//std::cout << "num_EqualDoubleDouble: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_double(Y));
    return emake_boole(Equal_to_7bits(edouble_number(X), edouble_number(Y)));
}


ex num_EqualIntReal(er X, er Y)
{
//std::cout << "num_EqualIntReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_real(Y));
    return emake_boole(arb_contains_fmpz(ereal_data(Y), eint_data(X)));
}


ex num_EqualRatReal(er X, er Y)
{
//std::cout << "num_EqualRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_real(Y));
    return emake_boole(arb_contains_fmpq(ereal_data(Y), erat_data(X)));
}

ex num_EqualDoubleReal(er X, er Y)
{
//std::cout << "num_EqualDoubleReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_real(Y));
    double y = arf_get_d(arb_midref(ereal_data(Y)), ARF_RND_DOWN);
    return emake_boole(Equal_to_7bits(edouble_number(X), y));
}

ex num_EqualRealReal(er X, er Y)
{
//std::cout << "num_EqualRealReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_real(Y));
    return emake_boole(arb_overlaps(ereal_data(X), ereal_data(Y)));
}


ex num_Equal_Cmplx(er X, er Y)
{
//std::cout << "num_Equal_Cmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(Y));
    ex c = num_Equal(eget_cint(0), ecmplx_imag(Y));
    if (c == nullptr || etor(c) == gs.sym_sFalse.get())
        return c;
    eclear(c);
    return num_Equal(X, ecmplx_real(Y));
}

ex num_EqualCmplxCmplx(er X, er Y)
{
//std::cout << "num_EqualCmplxCmplx: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_cmplx(Y));
    ex c = num_Equal(ecmplx_imag(X), ecmplx_imag(Y));
    if (c == nullptr || etor(c) == gs.sym_sFalse.get())
        return c;
    eclear(c);
    return num_Equal(ecmplx_real(X), ecmplx_real(Y));
}

ex num_EqualIntNan(er X, er Y)
{
//std::cout << "num_EqualIntNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_nan(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualNanInt(er X, er Y)
{
//std::cout << "num_EqualIntNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_int(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualRatNan(er X, er Y)
{
//std::cout << "num_EqualRatNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_nan(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualNanRat(er X, er Y)
{
//std::cout << "num_EqualRatNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_rat(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualDoubleNan(er X, er Y)
{
//std::cout << "num_EqualDoubleNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_double(X));
    assert(eis_nan(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualNanDouble(er X, er Y)
{
//std::cout << "num_EqualDoubleNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_double(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualRealNan(er X, er Y)
{
//std::cout << "num_EqualRealNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_nan(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualNanReal(er X, er Y)
{
//std::cout << "num_EqualRealNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_real(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualCmplxNan(er X, er Y)
{
//std::cout << "num_EqualCmplxNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_cmplx(X));
    assert(eis_nan(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualNanCmplx(er X, er Y)
{
//std::cout << "num_EqualCmplxNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_cmplx(Y));

    return gs.sym_sFalse.copy();
}

ex num_EqualNanNan(er X, er Y)
{
//std::cout << "num_EqualNanNan: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_nan(X));
    assert(eis_nan(Y));

    return gs.sym_sFalse.copy();
}


// could return nullptr for unevaluated
ex num_Equal(er X, er Y)
{
//std::cout << "num_Equal2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_EqualIntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_EqualIntRat(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_EqualIntRat(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_EqualRatRat(Y, X);

        case ETYPE_number * ETYPE_INT + ETYPE_DOUBLE:
            return num_EqualIntDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
            return num_EqualIntDouble(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_DOUBLE:
            return num_EqualRatDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return num_EqualRatDouble(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_DOUBLE:
            return num_EqualDoubleDouble(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return num_EqualIntReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return num_EqualIntReal(Y, X);
		case ETYPE_number * ETYPE_RAT + ETYPE_REAL:
            return num_EqualRatReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return num_EqualRatReal(Y, X);
		case ETYPE_number * ETYPE_DOUBLE + ETYPE_REAL:
            return num_EqualDoubleReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_DOUBLE:
            return num_EqualDoubleReal(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_EqualRealReal(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_CMPLX:
        case ETYPE_number * ETYPE_RAT + ETYPE_CMPLX:
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_CMPLX:
        case ETYPE_number * ETYPE_REAL + ETYPE_CMPLX:
            return num_Equal_Cmplx(X, Y);
            return num_Equal_Cmplx(X, Y);
            return num_Equal_Cmplx(X, Y);
            return num_Equal_Cmplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_INT:
        case ETYPE_number * ETYPE_CMPLX + ETYPE_RAT:
        case ETYPE_number * ETYPE_CMPLX + ETYPE_DOUBLE:
        case ETYPE_number * ETYPE_CMPLX + ETYPE_REAL:
            return num_Equal_Cmplx(Y, X);
            return num_Equal_Cmplx(Y, X);
            return num_Equal_Cmplx(Y, X);
            return num_Equal_Cmplx(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_CMPLX:
            return num_EqualCmplxCmplx(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_NAN:
            return num_EqualIntNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_INT:
            return num_EqualNanInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_NAN:
            return num_EqualRatNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_RAT:
            return num_EqualNanRat(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_NAN:
            return num_EqualDoubleNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_DOUBLE:
            return num_EqualNanDouble(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_NAN:
            return num_EqualRealNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_REAL:
            return num_EqualNanReal(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_NAN:
            return num_EqualCmplxNan(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_CMPLX:
            return num_EqualNanCmplx(X, Y);
        case ETYPE_number * ETYPE_NAN + ETYPE_NAN:
            return num_EqualNanNan(X, Y);

        default:
            assert(false);
            return nullptr;
    }
}
