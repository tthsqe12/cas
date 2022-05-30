#include <cmath>
#include <cfloat>
#include <complex>

#include "timing.h"
#include "uex.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "flint/arith.h"


static ex _exp_to_power(ex E)
{
    uex e(E);
    return emake_node(gs.sym_sPower.copy(), gs.sym_sE.copy(), e.copychild(1));
}

static ex _can_simplify_exp_arg(slong & sign, er x)
{
    if (!ehas_head_sym_length(x, gs.sym_sTimes.get(), 2))
        return nullptr;

    if (!eis_sym(echild(x,2), gs.sym_sPi.get()))
        return nullptr;

    er y = echild(x, 1);

    if (!eis_cmplx(y))
        return nullptr;

    if (!eis_zero(ecmplx_real(y)))
        return nullptr;

    y = ecmplx_imag(y);

    if (eis_int(y))
    {
        sign = fmpz_is_even(eint_data(y)) ? 1 : -1;
        return emake_cint(0);
    }
    else if (eis_rat(y))
    {
        xfmpq_t q;
        fmpz_fdiv_qr(fmpq_denref(q.data), fmpq_numref(q.data),
                     fmpq_numref(erat_data(y)), fmpq_denref(erat_data(y)));

        if (fmpz_is_zero(fmpq_denref(q.data)))
            return nullptr;

        sign = fmpz_is_even(fmpq_denref(q.data)) ? 1 : -1;
        fmpz_set(fmpq_denref(q.data), fmpq_denref(erat_data(y)));
        ex t = emake_rat_move(q.data);
        t = emake_cmplx(emake_cint(0), t);
        return emake_node(gs.sym_sTimes.copy(), t, gs.sym_sPi.copy());
    }
    else
    {
        return nullptr;
    }
}

ex ex_canonicalize_exp(ex e)
{
    uex E(e);
    er X = echild(e,1);

    if (ehas_head_sym_length(X, gs.sym_sLog.get(),1))
    {
        return ecopychild(X,1);
    }
    else if (ehas_head_sym(X, gs.sym_sTimes.get()))
    {
        slong sign;
        ex r = _can_simplify_exp_arg(sign, X);
        if (r != nullptr)
        {
            r = emake_node(gs.sym_sPower.copy(), gs.sym_sE.copy(), r);
            if (sign != 1)
                r = emake_node(gs.sym_sTimes.copy(), emake_cint(sign), r);
            return r;
        }

        ulong n = elength(X);
        if (n > 1)
        {
            ulong c = 0;
            for (ulong i = 0; i < n; i++)
            {
                if (!ehas_head_sym_length(echild(X,i+1), gs.sym_sLog.get(),1))
                    continue;

                if (c == 0)
                    c = i+1;
                else
                    return _exp_to_power(ecopy(e));
            }
            if (c == 0)
                return _exp_to_power(ecopy(e));
            uex f(ecopy(X));
            f.removechild(c);
            return emake_node(gs.sym_sPower.copy(), ecopychild(X,c,1), f.release());
        }
    }
    else if (ehas_head_sym(X, gs.sym_sPlus.get()))
    {
        slong sign;
        ulong n = elength(X);
        for (ulong i = 0; i < n; i++)
        {
            er ei = echild(X, i+1);
            ex r = _can_simplify_exp_arg(sign, ei);
            if (r != nullptr)
            {
                r = ereplacechild(X, i+1, r);
                r = emake_node(gs.sym_sPower.copy(), gs.sym_sE.copy(), r);
                if (sign != 1)
                    r = emake_node(gs.sym_sTimes.copy(), emake_cint(sign), r);
                return r;
            }
        }
    }
    else if (ehas_head_sym_length(X, gs.sym_sArcSinh.get(), 1))
    {
        ex t = _arctrig_helper(ecopychild(X,1), 1, 1, 2);
        t = ex_powx(t, emake_crat(1,2));
        return ex_addx(ecopychild(X,1), t);
    }
    else if (ehas_head_sym_length(X, gs.sym_sArcCosh.get(), 1))
    {
        uex s(ex_addx(emake_cint(1), ecopychild(X,1)));
        s.setz(ex_powx(s.release(), emake_crat(1,2)));
        ex t = ex_addx(emake_cint(-1), ecopychild(X,1));
        t = ex_powx(t, emake_crat(1,2));
        t = ex_mulx(t, s.release());
        return ex_addx(ecopychild(X,1), t);
    }

    return _exp_to_power(ecopy(e));
}

ex ex_canonicalize_log(ex E)
{
    uex e(E);

    er X = echild(E,1);
    if (eis_sym(X, gs.sym_sE.get()))
        return emake_cint(1);

    return e.release();
}

ex ex_canonicalize_trig(ex E)
{
    return E;
#if 0
    if (eis_node(X) && elength(X) == 1)
    {
        if (ehas_head_sym(X, gs.sym_sArcTan.get()))
        {
            return ecopychild(X,1);
        }
        else if (ehas_head_sym(X, gs.sym_sArcSin.get()) ||
                 ehas_head_sym(X, gs.sym_sArcCos.get()))
        {
            ex a = _arctrig_helper(ecopychild(X,1), 1, -1, 2);
            a = ex_powx(a, emake_crat(1,2));
            return ehas_head_sym(X, gs.sym_sArcSin.get()) ?
                    ex_divx(ecopychild(X,1), a) :
                    ex_divx(a, ecopychild(X,1));
        }
        else if (ehas_head_sym(X, gs.sym_sArcCot.get()))
        {
            return ex_reciprocal(ecopychild(X,1));
        }
        else if (ehas_head_sym(X, gs.sym_sArcSec.get()) ||
                 ehas_head_sym(X, gs.sym_sArcCsc.get()))
        {
            ex a = _arctrig_helper(ecopychild(X,1), 1, -1, -2);
            a = ex_powx(a, emake_crat(1,2));
            a = ex_mulx(a, ecopychild(X,1));
            if (ehas_head_sym(X, gs.sym_sArcCsc.get()))
                a = ex_reciprocal(a);
            return a;   
        }
        return ecopy(e);
    }
    else if (ehas_head_sym(X, gs.sym_sTimes.get()) && elength(X) > 0)
    {
        bool switched;
        ex f;
        if (eis_cmplx(echild(X,1)) && eis_zero(ecmplx_real(echild(X,1))))
        {
            uex a(ecopy(X));
            a.replacechild(1, ecopy(ecmplx_imag(echild(X,1))));
            a.setz(eval_fakeabs(switched, a.release()));
            a.setz(emake_node(gs.sym_sTanh.copy(), a.release()));
            f = emake_node(gs.sym_sTimes.copy(), gs.const_i.copy(), a.release());
            return switched ? ex_negate(f) : f;
        }
        else
        {
            f = eval_fakeabs(switched, ecopy(X));
            if (!switched)
            {
                eclear(f);
                return ecopy(e);
            }
            return ex_negate(emake_node(gs.sym_sTan.copy(), f));
        }
    }
    else
    {
        return ecopy(e);
    }
#endif
}

ex ex_canonicalize_trigh(ex E)
{
    return E;
}

ex ex_canonicalize_arctrig(ex E)
{
    return E;
}

ex ex_canonicalize_arctrigh(ex E)
{
    return E;
}


ex dcode_sExp(er e)
{
//std::cout << "dcode_sExp: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sExp.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Exp(X, e);

    return ex_canonicalize_exp(ecopy(e));
}

ex dcode_sLog(er e)
{
//std::cout << "dcode_sLog: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sLog.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Log(X, e);

    return ex_canonicalize_log(ecopy(e));
}


ex dcode_sCos(er e)
{
//std::cout << "dcode_sCos: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCos.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Cos(X, e);

    return ex_canonicalize_trig(ecopy(e));
}


ex dcode_sSin(er e)
{
//std::cout << "dcode_sSin: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSin.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Sin(X, e);

    return ex_canonicalize_trig(ecopy(e));
}


ex dcode_sCosh(er e)
{
//std::cout << "dcode_sCosh: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCosh.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Cosh(X, e);

    return ex_canonicalize_trigh(ecopy(e));
}

ex dcode_sSinh(er e)
{
//std::cout << "dcode_sSinh: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSinh.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Sinh(X, e);

    return ex_canonicalize_trigh(ecopy(e));
}

ex dcode_sTan(er e)
{
//std::cout << "dcode_sTan: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sTan.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Tan(X, e);

    return ex_canonicalize_trig(ecopy(e));
}

ex dcode_sTanh(er e)
{
//std::cout << "dcode_sTanh: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sTanh.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_Tanh(X, e);

    return ex_canonicalize_trigh(ecopy(e));
}

ex dcode_sLog10(er e)
{
    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}


ex dcode_sArcCos(er e)
{
//std::cout << "dcode_sArcCos: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcCos.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_ArcCos(X, e);

    return ex_canonicalize_arctrig(ecopy(e));
}

ex dcode_sArcCosh(er e)
{
//std::cout << "dcode_sArcCosh: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcCosh.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_ArcCosh(X, e);

    return ex_canonicalize_arctrigh(ecopy(e));
}


ex dcode_sArcSin(er e)
{
//std::cout << "dcode_sArcSin: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcSin.get()));

    if (elength(e) != 1)
    {
        return _handle_message_argx1(e);
    }

    er X = echild(e,1);
    if (eis_number(X))
        return num_ArcSin(X, e);

    return ex_canonicalize_arctrig(ecopy(e));
}

ex dcode_sArcSinh(er e)
{
//std::cout << "dcode_sArcSinh: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcSinh.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_ArcSinh(X, e);

    return ex_canonicalize_arctrigh(ecopy(e));
}


ex ucode_ArcTan(er e)
{
    assert(eis_node(e));
    return ecopy(e);
}

ex dcode_sArcTan(er e)
{
//std::cout << "dcode_sArcTan: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcTan.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_ArcTan(X, e);

    return ex_canonicalize_arctrig(ecopy(e));
}


ex dcode_sArcTanh(er e)
{
//std::cout << "dcode_sArcTanh: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcTanh.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er X = echild(e,1);
    if (eis_number(X))
        return num_ArcTanh(X, e);

    return ex_canonicalize_arctrigh(ecopy(e));
}

ex dcode_sArcCot(er e)
{
//std::cout << "dcode_sArcCot: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcCot.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sArcCoth(er e)
{
//std::cout << "dcode_sArcCoth: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcCoth.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sArcCsc(er e)
{
//std::cout << "dcode_sArcCsc: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcCsc.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sArcCsch(er e)
{
//std::cout << "dcode_sArcCoth: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcCsch.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sArcSec(er e)
{
//std::cout << "dcode_sArcSec: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcSec.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sArcSech(er e)
{
//std::cout << "dcode_sArcSech: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sArcSech.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}


ex dcode_sCot(er e)
{
//std::cout << "dcode_sCot: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCot.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sCoth(er e)
{
//std::cout << "dcode_sCoth: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCoth.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}


ex dcode_sCsc(er e)
{
//std::cout << "dcode_sCot: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCsc.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sCsch(er e)
{
//std::cout << "dcode_sCoth: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCsch.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}


ex dcode_sSec(er e)
{
//std::cout << "dcode_sSec: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSec.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}

ex dcode_sSech(er e)
{
//std::cout << "dcode_sSech: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSech.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}
