#include <cmath>
#include <cfloat>

#include "timing.h"
#include "ex_print.h"
#include "ex_cont.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "flint/arith.h"

slong max_extra_precision = 2000;

ex eval_num_limit_prec(er e)
{
    if (eis_parray(e))
    {
        uex r(emake_parray_double());
        return parray_set_ex(eto_parray_double(r.get()), e) ? r.release() : ecopy(e);
    }
    else if (eis_node(e))
    {
        size_t n = elength(e);
        uex r; r.init_push_backr(echild(e,0), n);
        bool changed = false;
        for (size_t i = 0; i < n; i++)
        {
            er ei = echild(e,i+1);
            ex ri = eval_num_limit_prec(ei);
            changed = changed || ei != etor(ri);
            r.push_back(ri);
        }
        return changed ? r.release() : ecopy(e);
    }

    switch (etype(e))
    {
        case ETYPE_INT:
        {
            return emake_double(fmpz_get_d(eint_data(e)));
        }
        case ETYPE_RAT:
        {
            return emake_double(fmpq_get_d(erat_data(e)));
        }
        case ETYPE_REAL:
        {
            return emake_double(arf_get_d(arb_midref(ereal_data(e)), ARF_RND_NEAR));
        }
        case ETYPE_CMPLX:
        {
            uex re(eval_num_limit_prec(ecmplx_real(e)));
            ex im = eval_num_limit_prec(ecmplx_imag(e));
            return emake_cmplx(re.release(), im);
        }
        default:
        {
            return ecopy(e);
        }
    }
}

ex eval_num_limit_prec(er e, double p)
{
    if (!eis_real(e))
        return ecopy(e);

    double t = 3.3219280948873623479*p + 1;
    slong n = t;
    t -= n;
    mag_t z1, z2;
    mag_init(z1);
    mag_init(z2);
    arb_get_mag_lower(z1, ereal_data(e));
    mag_set_d_lower(z2, pow(2.0, -t));
    mag_mul_2exp_si(z2, z2, -n);
    mag_mul_lower(z1, z1, z2);
    xarb_t z;
    arb_set(z.data, ereal_data(e));
    if (mag_cmp(z1, arb_radref(z.data)) > 0)
        mag_swap(z1, arb_radref(z.data));

    return emake_real_move(z);
}


ex eval_num(er e, slong prec)
{
    if (eis_sym(e))
    {
        if (esym_ncode(e) != nullptr)
            return esym_ncode(e)(e, prec);
        else
            return ecopy(e);
    }

    switch (etype(e))
    {
        case ETYPE_INT:
        {
            xarb_t z;
            arb_set_round_fmpz(z.data, eint_data(e), prec + EXTRA_PRECISION_BASIC);
            z.limit_prec(prec + 4);
            return emake_real_move(z);
        }
        case ETYPE_RAT:
        {
            xarb_t z;
            arb_set_fmpq(z.data, erat_data(e), prec + EXTRA_PRECISION_BASIC);
            z.limit_prec(prec + 4);
            return emake_real_move(z);
        }
        case ETYPE_REAL:
        {
            xarb_t z;
            arb_set(z.data, ereal_data(e));
            z.limit_prec(prec + 4);
            return emake_real_move(z);
        }
        default:
        {
            return ecopy(e);
        }

        case ETYPE_NODE:;
    }

    ulong len = elength(e);
    ex h = eval_num(echild(e,0), prec);

    uex f; f.init_push_backx(h, len);

    er hh = esymbolic_head(etor(h));
    uint32_t attr = enormalattr(hh);
    uint32_t test = ATTR_NHoldFirst;
    for (ulong i = 0; i < len; i++)
    {
        er ei = echild(e,i+1);
        f.push_back((attr & test) ? ecopy(ei) : eval_num(ei, prec));
        test = ATTR_NHoldRest;
    }

    if (esym_ncode(hh) != nullptr)
        f.setnz(esym_ncode(hh)(f.get(), prec));
    else
        f.setz(eval(f.release()));

    return f.release();
}

ex eval_num_extra(er e, slong prec, slong extra)
{
    slong cur_prec = prec;
    uex f(eval_num(e, cur_prec));
    while (eis_real(f.get()))
    {
        slong fprec = ereal_number(f.get()).prec_bits();
        if (fprec >= prec || cur_prec - prec >= extra)
            return f.release();
        cur_prec = prec + std::min(1 + (cur_prec - fprec + 4)/4*5, extra);
        f.setnz(eval_num(e, cur_prec));
    }
    return f.release();
}

ex topeval_num(er e, double p) // prec in decimal
{
    uex f(eval_num_extra(e, 3.3219280948873623479*p + 4, max_extra_precision));
    f.setnz(eval_num_limit_prec(f.get(), p));
    return f.release();
}

ex topeval_num(er e) // machine precision
{
    uex f(eval_num_extra(e, 53, max_extra_precision));
    f.setnz(eval_num_limit_prec(f.get()));
    return f.release();
}

ex dcode_sN(er e)
{
    if (elength(e) == 1)
    {
        return topeval_num(echild(e,1));
    }
    else if (elength(e) == 2)
    {
        uex p(eval_num(echild(e,2), 53));
        p.setnz(eval_num_limit_prec(p.get()));
        if (!eis_double(p.get()))
        {
            _gen_message(echild(e,0), "precbd", NULL, ecopychild(e,2));
            return ecopy(e);
        }
        return topeval_num(echild(e,1), edouble_number(p.get()));
    }
    else
    {
        return ecopy(e);
    }
}


ex dcode_sPrecision(er e)
{
//std::cout << "dcode_sPrecision: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sPrecision.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    if (!eis_number(x))
        return ecopy(e);

    if (!eis_real(x))
        return gs.const_infinity.copy();

    if (arb_contains_zero(ereal_data(x)))
        return emake_cint(0);

    ex z = emake_real();
    arb_get_rad_arb(ereal_data(z), ereal_data(x));
    arb_div(ereal_data(z), ereal_data(z), ereal_data(x), 32);
    arb_mul_2exp_si(ereal_data(z), ereal_data(z), 1);
    arb_abs(ereal_data(z), ereal_data(z));
    arb_log_base_ui(ereal_data(z), ereal_data(z), 10, 32 + fmpz_bits(ARF_EXPREF(arb_midref(ereal_data(z)))));
    arb_add_error_2exp_si(ereal_data(z), -30);
    arb_neg(ereal_data(z), ereal_data(z));
    return efix_real(z);
}

ex dcode_sAccuracy(er e)
{
//std::cout << "dcode_sAccuracy: " << ex_tostring_full(e) << std::endl; 
    assert(ehas_head_sym(e, gs.sym_sAccuracy.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    if (!eis_number(x))
        return ecopy(e);

    if (!eis_real(x))
        return gs.const_infinity.copy();

    ex z = emake_real();
    arb_get_rad_arb(ereal_data(z), ereal_data(x));
    arb_mul_2exp_si(ereal_data(z), ereal_data(z), 1);
    arb_abs(ereal_data(z), ereal_data(z));
    arb_log_base_ui(ereal_data(z), ereal_data(z), 10, 32 + fmpz_bits(ARF_EXPREF(arb_midref(ereal_data(z)))));
    arb_add_error_2exp_si(ereal_data(z), -30);
    arb_neg(ereal_data(z), ereal_data(z));
    return efix_real(z);
}

ex ex_logr(er a)
{
    return ex_canonicalize_log(emake_node(gs.sym_sLog.copy(), ecopy(a)));
}

ex ex_logx(ex a)
{
    return ex_canonicalize_log(emake_node(gs.sym_sLog.copy(), a));
}

ex ex_expr(er a)
{
    return ex_canonicalize_exp(emake_node(gs.sym_sLog.copy(), ecopy(a)));
}

ex ex_expx(ex a)
{
    return ex_canonicalize_exp(emake_node(gs.sym_sLog.copy(), a));
}

ex ex_max(ex a, ex b)
{
    return ex_canonicalize_max(emake_node(gs.sym_sMax.copy(), a, b));
}

ex ex_mulx(ex a, ex b)
{
    return ex_canonicalize_times(emake_node(gs.sym_sTimes.copy(), a, b));
}

ex ex_mulx(ex a, ex b, ex c)
{
    return ex_canonicalize_times(emake_node(gs.sym_sTimes.copy(), a, b, c));
}

ex ex_mulx(ex a, ex b, ex c, ex d)
{
    return ex_canonicalize_times(emake_node(gs.sym_sTimes.copy(), a, b, c, d));
}

ex ex_mulr(er a, er b)
{
    return ex_canonicalize_times(emake_node(gs.sym_sTimes.copy(), ecopy(a), ecopy(b)));
}

ex ex_mul_six(ex a, slong b)
{
    return ex_canonicalize_times(emake_node(gs.sym_sTimes.copy(), a, emake_int_si(b)));
}

ex ex_subx(ex a, ex b)
{
    return ex_canonicalize_plus(emake_node(gs.sym_sPlus.copy(), a, ex_negate(b)));
}

ex ex_subr(er a, er b)
{
    return ex_canonicalize_plus(emake_node(gs.sym_sPlus.copy(), ecopy(a), ex_negate(ecopy(b))));
}

ex ex_divx(ex a, ex b)
{
    return ex_canonicalize_times(emake_node(gs.sym_sTimes.copy(), a, ex_reciprocal(b)));
}

ex ex_divr(er a, er b)
{
    return ex_canonicalize_times(emake_node(gs.sym_sTimes.copy(), ecopy(a), ex_reciprocal(ecopy(b))));
}


ex ex_powx(ex a, ex b)
{
    return ex_canonicalize_power(emake_node(gs.sym_sPower.copy(), a, b));
}

ex ex_powr(er a, er b)
{
    return ex_canonicalize_power(emake_node(gs.sym_sPower.copy(), ecopy(a), ecopy(b)));
}


ex ex_floor(ex A)
{
    uex a(A);
    if (eis_int(A))
    {
        return a.release();
    }
    else if (eis_rat(A))
    {
        ex z = emake_int();
        fmpz_fdiv_q(eint_data(z), fmpq_numref(erat_data(A)), fmpq_denref(erat_data(A)));
        return z;
    }
    return gs.sym_sNull.copy();
}
