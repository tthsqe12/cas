#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"

double ereal_get_double(er e)
{
    double d = arf_get_d(arb_midref(ereal_data(e)), ARF_RND_NEAR);
    switch(std::fpclassify(d))
    {
        case FP_NORMAL:
        case FP_ZERO:
        default:
            return d;
        case FP_INFINITE:
        case FP_NAN:
        case FP_SUBNORMAL:
            return NAN;
    }
}

slong cmplx_fxn::extra_prec()
{
    return 0;
}

ex cmplx_fxn::eval_IntDouble(er X, er Y, er def)
{
    assert(eis_int(X));
    assert(eis_double(Y));
    double y = edouble_number(Y);
    if (fmpz_is_zero(eint_data(X)))
        return eval_double_imag(y, def);
    double x = eint_get_double(X);
    return eval_double(std::complex<double>(x, y), def);
}

ex cmplx_fxn::eval_DoubleInt(er X, er Y, er def)
{
    assert(eis_double(X));
    assert(eis_int(Y));
    double x = edouble_number(X);
    double y = eint_get_double(Y);
    return eval_double(std::complex<double>(x, y), def);
}

ex cmplx_fxn::eval_RatDouble(er X, er Y, er def)
{
    assert(eis_rat(X));
    assert(eis_double(Y));
    double x = erat_get_double(X);
    double y = edouble_number(Y);
    return eval_double(std::complex<double>(x, y), def);
}

ex cmplx_fxn::eval_DoubleRat(er X, er Y, er def)
{
    assert(eis_double(X));
    assert(eis_rat(Y));
    double x = edouble_number(X);
    double y = num_todouble(Y);
    return eval_double(std::complex<double>(x, y), def);
}

ex cmplx_fxn::eval_DoubleDouble(er X, er Y, er def)
{
    assert(eis_double(X));
    assert(eis_double(Y));
    double x = edouble_number(X);
    double y = edouble_number(Y);
    return eval_double(std::complex<double>(x, y), def);
}

ex cmplx_fxn::eval_IntReal(er X, er Y, er def)
{
    assert(eis_int(X));
    assert(eis_real(Y));
    if (fmpz_is_zero(eint_data(X)))
        return eval_arb_imag(ereal_number(Y), def);

    uex wre(emake_real());
    uex wim(emake_real());
    slong p = ereal_number(Y).wprec() + extra_prec();
    acb_t z, w;
    arb_init(acb_realref(z));
    arb_set_round_fmpz(acb_realref(z), eint_data(X), p);
    *acb_imagref(z) = *ereal_data(Y);
    *acb_realref(w) = *ereal_data(wre.get());
    *acb_imagref(w) = *ereal_data(wim.get());
    acb_eval(w, z, p);
    *ereal_data(wre.get()) = *acb_realref(w);
    *ereal_data(wim.get()) = *acb_imagref(w);
    arb_clear(acb_realref(z));
    return emake_cmplx(wre.release(), wim.release());
}

ex cmplx_fxn::eval_RealInt(er X, er Y, er def)
{
    assert(eis_real(X));
    assert(eis_int(Y));
    uex wre(emake_real());
    uex wim(emake_real());
    slong p = ereal_number(X).wprec() + extra_prec();
    acb_t z, w;
    *acb_realref(z) = *ereal_data(X);
    arb_init(acb_imagref(z));
    arb_set_round_fmpz(acb_imagref(z), eint_data(Y), p);
    *acb_realref(w) = *ereal_data(wre.get());
    *acb_imagref(w) = *ereal_data(wim.get());
    acb_eval(w, z, p);
    *ereal_data(wre.get()) = *acb_realref(w);
    *ereal_data(wim.get()) = *acb_imagref(w);
    arb_clear(acb_imagref(z));
    return emake_cmplx(wre.release(), wim.release());
}

ex cmplx_fxn::eval_RatReal(er X, er Y, er def)
{
    assert(eis_rat(X));
    assert(eis_real(Y));
    uex wre(emake_real());
    uex wim(emake_real());
    slong p = ereal_number(Y).wprec() + extra_prec();
    acb_t z, w;
    arb_init(acb_realref(z));
    arb_set_fmpq(acb_realref(z), erat_data(X), p);
    *acb_imagref(z) = *ereal_data(Y);
    *acb_realref(w) = *ereal_data(wre.get());
    *acb_imagref(w) = *ereal_data(wim.get());
    acb_eval(w, z, p);
    *ereal_data(wre.get()) = *acb_realref(w);
    *ereal_data(wim.get()) = *acb_imagref(w);
    arb_clear(acb_realref(z));
    return emake_cmplx(wre.release(), wim.release());
}

ex cmplx_fxn::eval_RealRat(er X, er Y, er def)
{
    assert(eis_real(X));
    assert(eis_rat(Y));
    uex wre(emake_real());
    uex wim(emake_real());
    slong p = ereal_number(X).wprec() + extra_prec();
    acb_t z, w;
    *acb_realref(z) = *ereal_data(X);
    arb_init(acb_imagref(z));
    arb_set_fmpq(acb_imagref(z), erat_data(Y), p);
    *acb_realref(w) = *ereal_data(wre.get());
    *acb_imagref(w) = *ereal_data(wim.get());
    acb_eval(w, z, p);
    *ereal_data(wre.get()) = *acb_realref(w);
    *ereal_data(wim.get()) = *acb_imagref(w);
    arb_clear(acb_imagref(z));
    return emake_cmplx(wre.release(), wim.release());
}

ex cmplx_fxn::eval_DoubleReal(er X, er Y, er def)
{
    assert(eis_double(X));
    assert(eis_real(Y));
    double x = edouble_number(X);
    double y = ereal_get_double(Y);
    return eval_double(std::complex<double>(x, y), def);
}

ex cmplx_fxn::eval_RealDouble(er X, er Y, er def)
{
    assert(eis_real(X));
    assert(eis_double(Y));
    double x = num_todouble(X);
    double y = edouble_number(Y);
    return eval_double(std::complex<double>(x, y), def);
}

ex cmplx_fxn::eval_RealReal(er X, er Y, er def)
{
    assert(eis_real(X));
    assert(eis_real(Y));
    acb_t z, w;
    *acb_realref(z) = *ereal_data(X);
    *acb_imagref(z) = *ereal_data(Y);
    uex wre(emake_real());
    uex wim(emake_real());
    *acb_realref(w) = *ereal_data(wre.get());
    *acb_imagref(w) = *ereal_data(wim.get());
    slong p = FLINT_MIN(ereal_number(X).wprec(),
                        ereal_number(Y).wprec());
    acb_eval(w, z, p + extra_prec());
    *ereal_data(wre.get()) = *acb_realref(w);
    *ereal_data(wim.get()) = *acb_imagref(w);
    return emake_cmplx(wre.release(), wim.release());
}

ex cmplx_fxn::eval(er a, er def)
{
    assert(eis_cmplx(a));
    er X = ecmplx_real(a);
    er Y = ecmplx_imag(a);
    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
        default:
            return eval_exact(a, def);
        case ETYPE_number * ETYPE_INT + ETYPE_DOUBLE:
            return eval_IntDouble(X, Y, def);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
           return eval_DoubleInt(X, Y, def);
        case ETYPE_number * ETYPE_RAT + ETYPE_DOUBLE:
            return eval_RatDouble(X, Y, def);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return eval_DoubleRat(X, Y, def);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_DOUBLE:
            return eval_DoubleDouble(X, Y, def);
        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return eval_IntReal(X, Y, def);
        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return eval_RealInt(X, Y, def);
        case ETYPE_number * ETYPE_RAT + ETYPE_REAL:
            return eval_RatReal(X, Y, def);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return eval_RealRat(X, Y, def);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_REAL:
            return eval_DoubleReal(X, Y, def);
        case ETYPE_number * ETYPE_REAL + ETYPE_DOUBLE:
            return eval_RealDouble(X, Y, def);
        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return eval_RealReal(X, Y, def);
    }
}




double num_todouble(er x)
{
//std::cout << "num_todouble: " << ex_tostring_full(x) << std::endl;
    assert(eis_number(x));

    uint32_t tx = etype(x);
    switch (tx)
    {
        case ETYPE_INT:
            return fmpz_get_d(eint_data(x));
        case ETYPE_RAT:
            return fmpq_get_d(erat_data(x));
        case ETYPE_DOUBLE:
            return eto_double(x)->number;
        case ETYPE_REAL:
        {
            double d = arf_get_d(arb_midref(ereal_data(x)), ARF_RND_NEAR);
            switch(std::fpclassify(d))
            {
                case FP_NORMAL:
                case FP_ZERO:
                default:
                    return d;
                case FP_INFINITE:
                case FP_NAN:
                case FP_SUBNORMAL:
                    return NAN;
            }
        }
        default:
            return NAN;
    }
}


double econvert_todouble(er e)
{
    if (eis_number(e))
    {
        return num_todouble(e);
    }
    else
    {
        return NAN;
    }
}

