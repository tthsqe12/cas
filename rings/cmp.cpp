#include "globalstate.h"
#include "arithmetic.h"
#include "ex_print.h"
#include "timing.h"

int num_CmpIntInt(er x, er y)
{
//std::cout << "num_CmpIntInt: " << ex_tostring_full(x) << ", " << ex_tostring_full(y) << std::endl;
    assert(eis_int(x));
    assert(eis_int(y));
    return fmpz_cmp(eint_data(x), eint_data(y));
}

int num_CmpRatInt(er x, er y)
{
    assert(eis_rat(x));
    assert(eis_int(y));
    return fmpq_cmp_fmpz(erat_data(x), eint_data(y));
}

int num_CmpRatRat(er x, er y)
{
//std::cout << "num_CmpRatRat: " << ex_tostring_full(x) << ", " << ex_tostring_full(y) << std::endl;
    assert(eis_rat(x));
    assert(eis_rat(y));
    return fmpq_cmp(erat_data(x), erat_data(y));
}

int num_CmpIntDouble(er X, er Y)
{
    double x = fmpz_get_d(eint_data(X));
    double y = edouble_number(Y);
    return (x > y) ? 1 : -1;
}

int num_CmpRatDouble(er X, er Y)
{
    return 0;
}

int num_CmpDoubleDouble(er X, er Y)
{
    double x = edouble_number(X);
    double y = edouble_number(Y);
    if (x < y)
        return -1;
    else if (x > y)
        return 1;
    else
        return 0;
}


int num_CmpIntReal(er X, er Y)
{
    arf_t f;
    arf_init(f);
    arf_set_fmpz(f, eint_data(X));
    int cmp = arf_cmp(f, arb_midref(ereal_data(Y)));
    arf_clear(f);
    return cmp > 0 ? 1 : -1;
}

int num_CmpRatReal(er X, er Y)
{
    arf_t f;
    arf_init(f);
    slong p = arf_bits(arb_midref(ereal_data(Y)));
    arf_fmpz_div_fmpz(f, fmpq_numref(erat_data(X)), fmpq_denref(erat_data(X)), p + 10, ARF_RND_NEAR);
    int cmp = arf_cmp(f, arb_midref(ereal_data(Y)));
    arf_clear(f);
    return cmp > 0 ? 1 : -1;
}

int num_CmpDoubleReal(er X, er Y)
{
    return 0;
}

int num_CmpRealReal(er X, er Y)
{
    int cmp = arf_cmp(arb_midref(ereal_data(X)), arb_midref(ereal_data(Y)));
    if (cmp != 0)
        return cmp > 0 ? 1 : -1;
    return mag_cmp(arb_radref(ereal_data(X)), arb_radref(ereal_data(Y)));
}

int num_CmpIntCplx(er X, er Y)
{
    int cmp = num_Cmp2(X, ecmplx_real(Y));
    if (cmp != 0)
        return cmp > 0 ? 1 : -1;
    return num_Cmp2(eget_cint(0), ecmplx_imag(Y));
}

int num_CmpRatCplx(er X, er Y)
{
    int cmp = num_Cmp2(X, ecmplx_real(Y));
    if (cmp != 0)
        return cmp > 0 ? 1 : -1;
    return num_Cmp2(eget_cint(0), ecmplx_imag(Y));
}

int num_CmpRealCplx(er X, er Y)
{
    int cmp = num_Cmp2(X, ecmplx_real(Y));
    if (cmp != 0)
        return cmp > 0 ? 1 : -1;
    return num_Cmp2(eget_cint(0), ecmplx_imag(Y));
}

int num_CmpCplxCplx(er X, er Y)
{
    int cmp = num_Cmp2(ecmplx_real(X), ecmplx_real(Y));
    if (cmp != 0)
        return cmp > 0 ? 1 : -1;
    return num_Cmp2(ecmplx_imag(X), ecmplx_imag(Y));
}

int num_Cmp2(er X, er Y)
{
//std::cout << "num_Cmp2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_CmpIntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return -num_CmpRatInt(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_CmpRatInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_CmpRatRat(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_DOUBLE:
            return num_CmpIntDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
            return -num_CmpIntDouble(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_DOUBLE:
            return num_CmpRatDouble(X, Y);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return -num_CmpRatDouble(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_DOUBLE:
            return num_CmpDoubleDouble(X, Y);


        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return num_CmpIntReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return -num_CmpIntReal(Y, X);
        case ETYPE_number * ETYPE_RAT + ETYPE_REAL:
            return num_CmpRatReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return -num_CmpRatReal(Y, X);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_REAL:
            return num_CmpDoubleReal(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_DOUBLE:
            return -num_CmpDoubleReal(Y, X);
        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_CmpRealReal(X, Y);



        case ETYPE_number * ETYPE_INT + ETYPE_CMPLX:
            return num_CmpIntCplx(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_CMPLX:
            return num_CmpRatCplx(X, Y);
        case ETYPE_number * ETYPE_REAL + ETYPE_CMPLX:
            return num_CmpRealCplx(X, Y);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_INT:
            return -num_CmpIntCplx(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_RAT:
            return -num_CmpRatCplx(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_REAL:
            return -num_CmpRealCplx(Y, X);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_CMPLX:
            return num_CmpCplxCplx(X, Y);

        default:
            return 0;
    }
}


bool num_Less_IntInt(er x, er y)
{
//std::cout << "num_CmpIntInt: " << ex_tostring_full(x) << ", " << ex_tostring_full(y) << std::endl;
    assert(eis_int(x));
    assert(eis_int(y));
    return fmpz_cmp(eint_data(x), eint_data(y)) < 0;
}

bool num_Less_IntRat(er x, er y)
{
    assert(eis_int(x));
    assert(eis_rat(y));
    fmpq_t t;
    *fmpq_numref(t) = *(eint_data(x));
    *fmpq_denref(t) = WORD(1);
    return fmpq_cmp(t, erat_data(y)) < 0;
}

bool num_Less_RatInt(er y, er x)
{
    assert(eis_int(x));
    assert(eis_rat(y));
    fmpq_t t;
    *fmpq_numref(t) = *(eint_data(x));
    *fmpq_denref(t) = WORD(1);
    return fmpq_cmp(erat_data(y), t) < 0;
}

bool num_Less_RatRat(er x, er y)
{
//std::cout << "num_CmpRatRat: " << ex_tostring_full(x) << ", " << ex_tostring_full(y) << std::endl;
    assert(eis_rat(x));
    assert(eis_rat(y));
    return fmpq_cmp(erat_data(x), erat_data(y)) < 0;
}

bool num_Less_RealReal(er X, er Y)
{
    return arb_lt(ereal_data(X), ereal_data(Y));
}

bool num_Less2(er X, er Y)
{
//std::cout << "num_Less2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_Less_IntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_Less_IntRat(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_Less_RatInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_Less_RatRat(X, Y);

        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_Less_RealReal(X, Y);

        default:
            assert(false);
            return false;
    }
}


bool num_LessEqual_IntInt(er x, er y)
{
//std::cout << "num_CmpIntInt: " << ex_tostring_full(x) << ", " << ex_tostring_full(y) << std::endl;
    assert(eis_int(x));
    assert(eis_int(y));
    return fmpz_cmp(eint_data(x), eint_data(y)) <= 0;
}

bool num_LessEqual_IntRat(er x, er y)
{
    assert(eis_int(x));
    assert(eis_rat(y));
    fmpq_t t;
    *fmpq_numref(t) = *(eint_data(x));
    *fmpq_denref(t) = WORD(1);
    return fmpq_cmp(t, erat_data(y)) <= 0;
}

bool num_LessEqual_RatInt(er y, er x)
{
    assert(eis_int(x));
    assert(eis_rat(y));
    fmpq_t t;
    *fmpq_numref(t) = *(eint_data(x));
    *fmpq_denref(t) = WORD(1);
    return fmpq_cmp(erat_data(y), t) <= 0;
}

bool num_LessEqual_RatRat(er x, er y)
{
//std::cout << "num_CmpRatRat: " << ex_tostring_full(x) << ", " << ex_tostring_full(y) << std::endl;
    assert(eis_rat(x));
    assert(eis_rat(y));
    return fmpq_cmp(erat_data(x), erat_data(y)) <= 0;
}

bool num_LessEqual_RealReal(er X, er Y)
{
    return arb_le(ereal_data(X), ereal_data(Y));
}

bool num_LessEqual2(er X, er Y)
{
//std::cout << "num_Cmp2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_LessEqual_IntInt(X, Y);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_LessEqual_IntRat(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_LessEqual_RatInt(X, Y);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_LessEqual_RatRat(X, Y);

        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_LessEqual_RealReal(X, Y);

        default:
            assert(false);
            return false;
    }
}


