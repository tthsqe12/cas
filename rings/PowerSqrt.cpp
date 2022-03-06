#include "globalstate.h"
#include "eval.h"
#include "code.h"
#include "arithmetic.h"
#include "flintarb_wrappers.h"
#include "xacb_t.h"


void fmpz_factor_combine_like_exps(fmpz_factor_t f)
{
    slong i, j, in, out;

    for (i = 1; i < f->num; i++)
    for (j = i; j > 0 && f->exp[j - 1] > f->exp[j]; j--)
    {
        fmpz_swap(f->p + j - 1, f->p + j);
        ULONG_SWAP(f->exp[j - 1], f->exp[j]);
    }

    out = -1;
    for (in = 0; in < f->num; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && f->exp[out] == f->exp[in])
        {
            fmpz_mul(f->p + out, f->p + out, f->p + in);
        }
        else
        {
            out++;
            f->exp[out] = f->exp[in];
            fmpz_swap(f->p + out, f->p + in);
        }
    }

    out++;

    f->num = out;
}

/* (x/y)^(p/q) */
ex _num_PowerRatRat(const fmpz_t x, const fmpz_t y,
                    const fmpz_t p, const fmpz_t q)
{
    if (fmpz_sgn(p) == 0)
        return emake_cint(1);

    if (fmpz_sgn(p) < 0)
    {
        xfmpz_t mp;
        fmpz_neg(mp.data, p);
        return _num_PowerRatRat(y, x, mp.data, q);
    }

    xfmpz_factor_t fx, fy;

    assert(!fmpz_is_zero(x));
    assert(!fmpz_is_zero(y));

    fmpz_factor_smooth(fx.data, x, 12, 0);
    fmpz_factor_smooth(fy.data, y, 12, 0);

    assert(fmpz_sgn(p) > 0);
    assert(fmpz_sgn(q) >= 0);

/*
    (+- b1^e1 * ... *bn^en)^(p/q)

    = (+-)^(p/q) b1^(e1*p/q) * ... *bn^(en*p/q)

    then factor out integers to make all exponents mod 1
*/
    xfmpq_t a(1,1), pow;
    xfmpz_t P, Q, R;
    for (slong i = 0; i < fx.data->num; i++)
    {
        fmpz_one(R.data);
        fmpz_mul_ui(P.data, p, fx.data->exp[i]);
        fmpz_fdiv_qr(Q.data, R.data, P.data, q);
        if (fmpz_is_zero(R.data))
        {
            assert(fmpz_fits_si(Q.data));
            fmpz_pow_ui(R.data, fx.data->p + i, fmpz_get_ui(Q.data));
            fmpq_mul_fmpz(a.data, a.data, R.data);
            fx.data->exp[i] = fx.data->exp[fx.data->num - 1];
            fmpz_swap(fx.data->p + i, fx.data->p + fx.data->num - 1);
            fx.data->num--;
        }
        else
        {
            fx.data->exp[i] = fmpz_get_ui(P.data);
        }
    }
    fmpz_factor_combine_like_exps(fx.data);

    for (slong i = 0; i < fy.data->num; i++)
    {
        fmpz_one(R.data);
        fmpz_mul_ui(P.data, p, fy.data->exp[i]);
        fmpz_fdiv_qr(Q.data, R.data, P.data, q);
        if (fmpz_is_zero(R.data))
        {
            assert(fmpz_fits_si(Q.data));
            fmpz_pow_ui(R.data, fy.data->p + i, fmpz_get_ui(Q.data));
            fmpq_div_fmpz(a.data, a.data, R.data);
            fy.data->exp[i] = fy.data->exp[fy.data->num - 1];
            fmpz_swap(fy.data->p + i, fy.data->p + fy.data->num - 1);
            fy.data->num--;
        }
        else
        {
            fy.data->exp[i] = fmpz_get_ui(P.data);
        }
    }
    fmpz_factor_combine_like_exps(fy.data);

    std::vector<wex> v;

    if (fx.data->sign*fy.data->sign < 0)
    {
        fmpz_fdiv_qr(Q.data, R.data, p, q);
        fmpq_set_fmpz_frac(pow.data, R.data, q);
        if (fmpz_is_odd(Q.data))
            fmpq_neg(a.data, a.data);
    }

    if (!fmpq_is_one(a.data))
        v.push_back(wex(emake_rat_move(a.data)));

    if (fx.data->sign*fy.data->sign < 0)
    {
        if (fmpz_is_one(fmpq_numref(pow.data)) &&
            fmpz_cmp_ui(fmpq_denref(pow.data), 2) == 0)
        {
            v.push_back(wex(gs.const_i.copy()));
        }
        else
        {
            ex e = emake_rat_move(pow.data);
            v.push_back(wex(emake_node(gs.sym_sPower.copy(), emake_cint(-1), e)));            
        }
    }

    for (slong i = 0, j = 0; i < fx.data->num || j < fy.data->num; )
    {
        if (i < fx.data->num && (j >= fy.data->num || fx.data->exp[i] < fy.data->exp[j]))
        {
            fmpz_set_ui(R.data, fx.data->exp[i]);
            fmpq_set_fmpz_frac(pow.data, R.data, q);
            if (!fmpq_is_zero(pow.data))
            {
                uex e2(emake_rat_move(pow.data));
                ex e1 = emake_int_move(fx.data->p + i);
                v.push_back(wex(emake_node(gs.sym_sPower.copy(), e1, e2.release())));
            }
            i++;
        }
        else if (j < fy.data->num && (i >= fx.data->num || fy.data->exp[j] < fx.data->exp[i]))
        {
            fmpz_set_ui(R.data, fy.data->exp[j]);
            fmpq_set_fmpz_frac(pow.data, R.data, q);
            if (!fmpq_is_zero(pow.data))
            {
                fmpq_neg(pow.data, pow.data);
                uex e2(emake_rat_move(pow.data));
                ex e1 = emake_int_move(fy.data->p + j);
                v.push_back(wex(emake_node(gs.sym_sPower.copy(), e1, e2.release())));
            }
            j++;
        }
        else
        {
            fmpz_set_ui(R.data, fx.data->exp[i]);
            fmpq_set_fmpz_frac(pow.data, R.data, q);
            if (!fmpq_is_zero(pow.data))
            {
                uex e2(emake_rat_move(pow.data));
                fmpq_set_fmpz_frac(a.data, fx.data->p + i, fy.data->p + j);
                ex e1 = emake_rat_move(a.data);
                v.push_back(wex(emake_node(gs.sym_sPower.copy(), e1, e2.release())));
            }
            i++;
            j++;
        }
    }

    ex r = emake_node_times(v);
    return r;
}



ex num_PowerIntInt(er X, er Y, er def)
{
//std::cout << "num_PowerIntInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_int(Y));

    if (fmpz_is_zero(eint_data(X)))
    {
        int s = fmpz_sgn(eint_data(Y));
        if (s > 0)
            return ecopy(X);
        else if (s < 0)
            return gs.const_complexinfinity.copy();
        else
            return gs.const_indeterminate.copy();
    }
    else if (fmpz_is_pm1(eint_data(X)))
    {
        if (fmpz_is_one(eint_data(X)) || fmpz_is_even(eint_data(Y)))
            return emake_cint(1);
        else
            return emake_cint(-1);
    }
    else if (fmpz_fits_si(eint_data(Y)))
    {
        slong n = fmpz_get_si(eint_data(Y));
        if (n > 0)
        {
            ex z = emake_int();
            fmpz_pow_ui(eint_data(z), eint_data(X), n);
            return efix_int(z);
        }
        else if (n < 0)
        {
            ex Z = emake_rat();
            fmpz_pow_ui(fmpq_denref(erat_data(Z)), eint_data(X), -n);
            if (fmpz_sgn(fmpq_denref(erat_data(Z))) < 0)
            {
                fmpz_set_si(fmpq_numref(erat_data(Z)), -1);
                fmpz_neg(fmpq_denref(erat_data(Z)), fmpq_denref(erat_data(Z)));
            }
            else
            {
                fmpz_set_ui(fmpq_numref(erat_data(Z)), 1);
            }
            return efix_rat(Z);
        }
        else
        {
            return emake_cint(1);
        }
    }
    else
    {
        return fmpz_sgn(eint_data(Y)) > 0 ? gs.const_overflow.copy() : gs.const_underflow.copy();
    }
}

ex num_PowerIntRat(er X, er Y, er def)
{
    fmpz * x = eint_data(X);
    fmpq * y = erat_data(Y);

    if (fmpz_is_zero(x))
    {
        int s = fmpq_sgn(y);
        if (s > 0)
            return ecopy(X);
        else if (s < 0)
            return gs.const_complexinfinity.copy();
        else
            return gs.const_indeterminate.copy();
    }

    uex r(_num_PowerRatRat(x, eget_cint_data(1), fmpq_numref(y), fmpq_denref(y)));
    if (ex_same(r.get(), def))
        return ecopy(def);
    else
        return r.release();
}

ex num_PowerRatInt(er X, er Y, er def)
{
//std::cout << "num_PowerRatInt: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_rat(X));
    assert(eis_int(Y));

    if (fmpz_fits_si(eint_data(Y)))
    {
        slong n = fmpz_get_si(eint_data(Y));
        ex z = emake_rat();
        fmpq_pow_si(erat_data(z), erat_data(X), n);
        return efix_rat(z);
    }
    else
    {
        return emake_nan_Overflow();
    }
}

ex num_PowerRatRat(er X, er Y, er def)
{
    fmpq * x = erat_data(X);
    fmpq * y = erat_data(Y);

    uex r(_num_PowerRatRat(fmpq_numref(x), fmpq_denref(x),
                           fmpq_numref(y), fmpq_denref(y)));

    if (ex_same(r.get(), def))
        return ecopy(def);
    else
        return r.release();
}

ex num_PowerDoubleInt(er X, er Y, er def)
{
    assert(eis_double(X));
    assert(eis_int(Y));

    double rr = 1;
    double x = edouble_number(X);
    double xx = x;
    const fmpz * y = eint_data(Y);

    if (fmpz_fits_si(y))
    {
        slong yy = fmpz_get_si(y);

        for (ulong yu = FLINT_ABS(yy); yu > 0; yu /= 2)
        {
            if (yu & 1)
                rr *= xx;
            xx *= xx;
        }

        if (yy < 0)
            rr = 1/rr;
    }
    else
    {
        rr = pow(xx, fmpz_get_d(y));
    }

    return emake_double(rr);
}


ex num_PowerDoubleRat(er X, er Y, er def)
{
    assert(eis_double(X));
    assert(eis_rat(Y));

    double rr = 1;
    double x = edouble_number(X);
    double xx = x;
    const fmpq * y = erat_data(Y);
    double dy = fmpq_get_d(y);

    if (xx < 0)
    {
        rr = pow(-xx, dy);
        double pi = edouble_number(gs.const_double_pi.get());
        return emake_cmplx(emake_double(rr*cos(pi*dy)),
                           emake_double(rr*sin(pi*dy)));
    }
    else
    {
        rr = pow(xx, dy);
        return emake_double(rr);
    }
}


ex num_PowerRealInt(er X, er Y, er def)
{
//std::cout << "num_TimesRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_int(Y));

	ex z = emake_real();
    slong p = ereal_number(X).wprec();
    arb_pow_fmpz(ereal_data(z), ereal_data(X), eint_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_PowerRealReal(er X, er Y, er def)
{
//std::cout << "num_TimesRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_real(Y));

    ex z = emake_real();
    slong p = ereal_number(X).wprec();
    arb_pow(ereal_data(z), ereal_data(X), ereal_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_PowerIntReal(er X, er Y, er def)
{
//std::cout << "num_PowerIntReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_int(X));
    assert(eis_real(Y));

    slong p = ereal_number(Y).wprec();
    ex z = emake_real();
	arb_set_round_fmpz(ereal_data(z), eint_data(X), p + EXTRA_PRECISION_BASIC);
    arb_pow(ereal_data(z), ereal_data(z), ereal_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(z);
}

ex num_PowerRealRat(er X, er Y, er def)
{
//std::cout << "num_TimesRatReal: " << ex_tostring_full(X) << ", " << ex_tostring_full(Y) << std::endl;
    assert(eis_real(X));
    assert(eis_rat(Y));

    slong p = ereal_number(X).wprec();
    ex Z = emake_real();
    arb_pow_fmpq(ereal_data(Z), ereal_data(X), erat_data(Y), p + EXTRA_PRECISION_BASIC);
    return efix_real(Z);
}

ex num_PowerIntCplx(er X, er Y, er def)
{
    assert(false);
    return nullptr;
}

/***********************/
class cmplx_fxn_PowerInt : public cmplx_fxn {
public:
    const fmpz * P;
    cmplx_fxn_PowerInt(er PP, er ddef) : P(eint_data(PP)) {}
    slong extra_prec();
    ex eval_exact(er z, er def);
    ex eval_double(std::complex<double> x, er def);
    ex eval_double_imag(double y, er def);
    ex eval_arb_imag(const xarb_t& y, er def);
    void acb_eval(acb_t w, const acb_t z, slong prec);
};

slong cmplx_fxn_PowerInt::extra_prec()
{
    return EXTRA_PRECISION_BASIC;
}

ex cmplx_fxn_PowerInt::eval_exact(er z, er def)
{
    return ecopy(def);
}

ex cmplx_fxn_PowerInt::eval_double(std::complex<double> x, er def)
{
    std::complex<double> r(1.0, 0.0);

    if (fmpz_fits_si(P))
    {
        slong p = fmpz_get_si(P);

        if (p < 0)
            x = std::complex<double>(1.0, 0.0)/x;

        for (ulong pu = FLINT_ABS(p); pu > 0; pu /= 2)
        {
            if (pu & 1)
                r *= x;
            x *= x;
        }
    }
    else
    {
        r = pow(x, fmpz_get_d(P));
    }

    return emake_cmplx(emake_double(r.real()), emake_double(r.imag()));
}

ex cmplx_fxn_PowerInt::eval_double_imag(double y, er def)
{
    return gs.sym_sNull.copy();
}

ex cmplx_fxn_PowerInt::eval_arb_imag(const xarb_t& y, er def)
{
    return gs.sym_sNull.copy();
}

void cmplx_fxn_PowerInt::acb_eval(acb_t w, const acb_t z, slong prec)
{
    acb_mul(w, z, z, prec);
}

ex num_PowerCplxInt(er X, er Y, er def)
{
    cmplx_fxn_PowerInt f(Y, def);
    return f.eval(X, def);
}

/***********************/
class cmplx_fxn_PowerRat : public cmplx_fxn {
public:
    const fmpq * P;
    cmplx_fxn_PowerRat(er PP, er ddef) : P(erat_data(PP)) {}
    slong extra_prec();
    ex eval_exact(er z, er def);
    ex eval_double(std::complex<double> x, er def);
    ex eval_double_imag(double y, er def);
    ex eval_arb_imag(const xarb_t& y, er def);
    void acb_eval(acb_t w, const acb_t z, slong prec);
};

slong cmplx_fxn_PowerRat::extra_prec()
{
    return EXTRA_PRECISION_BASIC;
}

ex cmplx_fxn_PowerRat::eval_exact(er z, er def)
{
    return ecopy(def);
}

ex cmplx_fxn_PowerRat::eval_double(std::complex<double> x, er def)
{
    std::complex<double> r(1.0, 0.0);

    if (fmpq_equal_int(P, 1, 2))
    {
        r = sqrt(x);
    }
    else if (fmpq_equal_int(P, -1, 2))
    {
        r = std::complex<double>(1.0, 0.0)/sqrt(x);
    }
    else
    {
        r = pow(x, fmpq_get_d(P));
    }

    return emake_cmplx(emake_double(r.real()), emake_double(r.imag()));
}

ex cmplx_fxn_PowerRat::eval_double_imag(double y, er def)
{
    return gs.sym_sNull.copy();
}

ex cmplx_fxn_PowerRat::eval_arb_imag(const xarb_t& y, er def)
{
    return gs.sym_sNull.copy();
}

void cmplx_fxn_PowerRat::acb_eval(acb_t w, const acb_t z, slong prec)
{
    if (!fmpz_fits_si(fmpq_numref(P)) || !fmpz_abs_fits_ui(fmpq_denref(P)))
    {
        arb_t p;
        arb_init(p);
        arb_set_fmpq(p, P, prec);
        acb_pow_arb(w, z, p, prec);
        arb_clear(p);
        return;
    }

    slong n = fmpz_get_si(fmpq_numref(P));
    ulong m = fmpz_get_ui(fmpq_denref(P));
    if (m == 2)
    {
        if (n < 0)
        {
            acb_rsqrt(w, z, prec);
            n = -n;
            if (n != 1)
            {
                acb_pow_si(w, w, n, prec);
                return;
            }
        }

        acb_sqrt(w, z, prec);
    }
    else
    {
        acb_root_ui(w, z, m, prec);
    }

    if (n != 1)
        acb_pow_fmpz(w, w, fmpq_numref(P), prec);
}

ex num_PowerCplxRat(er X, er Y, er def)
{
    cmplx_fxn_PowerRat f(Y, def);
    return f.eval(X, def);
}

/*******************************/
ex num_PowerRatCplx(er X, er Y, er def)
{
    assert(false);
    return nullptr;
}

ex num_PowerCplxCplx(er X, er Y, er def)
{
    assert(false);
    return nullptr;
}


ex num_Power(er X, er Y, er def)
{
//std::cout << "num_Power2: " << ex_tostring(X) << ", " << ex_tostring(Y) << std::endl;

    assert(eis_number(X));
    assert(eis_number(Y));
    uint32_t tx = etype(X);
    uint32_t ty = etype(Y);
    switch (ETYPE_number * tx + ty)
    {
        case ETYPE_number * ETYPE_INT + ETYPE_INT:
            return num_PowerIntInt(X, Y, def);

        case ETYPE_number * ETYPE_INT + ETYPE_RAT:
            return num_PowerIntRat(X, Y, def);
        case ETYPE_number * ETYPE_RAT + ETYPE_INT:
            return num_PowerRatInt(X, Y, def);
        case ETYPE_number * ETYPE_RAT + ETYPE_RAT:
            return num_PowerRatRat(X, Y, def);

        case ETYPE_number * ETYPE_DOUBLE + ETYPE_INT:
            return num_PowerDoubleInt(X, Y, def);
        case ETYPE_number * ETYPE_DOUBLE + ETYPE_RAT:
            return num_PowerDoubleRat(X, Y, def);

        case ETYPE_number * ETYPE_REAL + ETYPE_INT:
            return num_PowerRealInt(X, Y, def);
        case ETYPE_number * ETYPE_REAL + ETYPE_RAT:
            return num_PowerRealRat(X, Y, def);

        case ETYPE_number * ETYPE_REAL + ETYPE_REAL:
            return num_PowerRealReal(X, Y, def);
        case ETYPE_number * ETYPE_INT + ETYPE_REAL:
            return num_PowerIntReal(X, Y, def);


        case ETYPE_number * ETYPE_INT + ETYPE_CMPLX:
            return num_PowerIntCplx(X, Y, def);
        case ETYPE_number * ETYPE_RAT + ETYPE_CMPLX:
            return num_PowerRatCplx(X, Y, def);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_INT:
            return num_PowerCplxInt(X, Y, def);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_RAT:
            return num_PowerCplxRat(X, Y, def);
        case ETYPE_number * ETYPE_CMPLX + ETYPE_CMPLX:
            return num_PowerCplxCplx(X, Y, def);

        default:
            assert(false);
            return nullptr;
    }
}

/***********************/
double cospi(double x)
{
    return cos(x*3.141592653589793238463);
}

double sinpi(double x)
{
    return sin(x*3.141592653589793238463);
}

ex emake_cmplx(std::complex<double> z)
{
    return emake_cmplx(z.real(), z.imag());
}

class cmplx_fxn_CisPi : public cmplx_fxn {
public:
    cmplx_fxn_CisPi() {}
    slong extra_prec();
    ex eval_exact(er z, er def);
    ex eval_double(std::complex<double> x, er def);
    ex eval_double_imag(double y, er def);
    ex eval_arb_imag(const xarb_t& y, er def);
    void acb_eval(acb_t w, const acb_t z, slong prec);
};

slong cmplx_fxn_CisPi::extra_prec()
{
    return EXTRA_PRECISION_BASIC;
}

ex cmplx_fxn_CisPi::eval_exact(er z, er def)
{
    return ecopy(def);
}

ex cmplx_fxn_CisPi::eval_double(std::complex<double> x, er def)
{
    return emake_cmplx(exp(std::complex<double>(0.0, 3.141592653589793238463)*x));
}

ex cmplx_fxn_CisPi::eval_double_imag(double y, er def)
{
    return emake_double(exp(-y));
}

ex cmplx_fxn_CisPi::eval_arb_imag(const xarb_t& y, er def)
{
    ex z = emake_real();
    arb_neg(ereal_data(z), y.data);
    arb_exp(ereal_data(z), ereal_data(z), y.wprec());
    return efix_real(z);
}

void cmplx_fxn_CisPi::acb_eval(acb_t w, const acb_t z, slong prec)
{
    acb_exp_pi_i(w, z, prec);
}

ex ncode_sPower(er e, slong prec)
{
    if (!ehas_head_sym_length(e, gs.sym_sPower.get(), 2))
        return ecopy(e);

    er e1 = echild(e,1);
    er e2 = echild(e,2);
    if (eis_int(e2))
    {
        return ex_powx(eval_num(e1, prec + 1), ecopy(e2));
    }
    else if (eis_sym(e1, gs.sym_sE.get()))
    {
        return ex_expx(eval_num(e2, prec + 1));
    }
    else if (eis_int(e1, -1))
    {
        if (eis_rat(e2))
        {
            xacb_t w;
            arb_sin_cos_pi_fmpq(acb_imagref(w.data), acb_realref(w.data), erat_data(e2), prec + 1);
            return emake_cmplx_move(w.data);
        }
        else
        {
            ex a = eval_num(e1, prec + 1);
            return ex_powx(ecopy(e1), a);
        }
    }
    else
    {
        uex a(eval_num(e1, prec + 1));
        ex b = eis_int_or_rat(e2) ? ecopy(e2) : eval_num(e2, prec + 1);
        return ex_powx(a.release(), b);
    }
}




/*********************************************/
ex num_Sqrt(er x)
{
    switch (etype(x))
    {
        case ETYPE_INT:
        {
            return num_PowerIntRat(x, eget_crat(1,2), x);
        }
        case ETYPE_RAT:
        {
            return num_PowerRatRat(x, eget_crat(1,2), x);
        }
        case ETYPE_DOUBLE:
        {
            double d = edouble_number(x);
            bool negative = d < 0;
            ex r = emake_double(sqrt(d < 0 ? -d : d));
            if (!negative)
                return r;
            else
                return emake_cmplx(emake_cint(0), r);
        }
        case ETYPE_REAL:
        {
            ex z = emake_real();
            slong p = ereal_number(x).wprec();
            arb_sqrt(ereal_data(z), ereal_data(x), p + EXTRA_PRECISION_BASIC);
            if (my_arb_is_ok(ereal_data(z)))
                return z;

            arb_neg(ereal_data(z), ereal_data(x));
            arb_sqrt(ereal_data(z), ereal_data(z), p + EXTRA_PRECISION_BASIC);
            if (my_arb_is_ok(ereal_data(z)))
                return emake_cmplx(emake_cint(0), z);

            eclear(z);

            xacb_t w;
            acb_set_arb(w.data, ereal_data(x));
            acb_sqrt(w.data, w.data, p + EXTRA_PRECISION_BASIC);
            return emake_cmplx_move(w.data);
        }
        default:
        {
            return ecopy(x);
        }
    }
}
