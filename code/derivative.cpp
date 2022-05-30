#include "timing.h"
#include "uex.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"


ex eval_diff_Log(ex F, er e)
{
    ex d;
    uex f(F);
    d = ex_reciprocal(ecopychild(e,1));
    return ex_mulx(f.release(), d);
}

ex eval_diff_Exp(ex F, er e)
{
    return ex_mulx(F, ecopy(e));
}

ex eval_diff_Sin(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sCos.copy(), ecopychild(e,1));
    return ex_mulx(f.release(), d);
}

ex eval_diff_Cos(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sSin.copy(), ecopychild(e,1));
    return ex_mulx(emake_cint(-1), f.release(), d);
}

ex eval_diff_Tan(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sSec.copy(), ecopychild(e,1));
    d = emake_node(gs.sym_sPower.copy(), d, emake_cint(2));
    return ex_mulx(f.release(), d);
}

ex eval_diff_Csc(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sCot.copy(), ecopychild(e,1));
    return ex_mulx(emake_cint(-1), f.release(), d, ecopy(e));
}

ex eval_diff_Sec(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sTan.copy(), ecopychild(e,1));
    return ex_mulx(f.release(), ecopy(e), d);
}

ex eval_diff_Cot(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sCsc.copy(), ecopychild(e,1));
    d = emake_node(gs.sym_sPower.copy(), d, emake_cint(2));
    return ex_mulx(emake_cint(-1), f.release(), d);
}

ex eval_diff_Sinh(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sCosh.copy(), ecopychild(e,1));
    return ex_mulx(f.release(), d);
}

ex eval_diff_Cosh(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sSinh.copy(), ecopychild(e,1));
    return ex_mulx(f.release(), d);
}

ex eval_diff_Tanh(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sSech.copy(), ecopychild(e,1));
    d = emake_node(gs.sym_sPower.copy(), d, emake_cint(2));
    return ex_mulx(f.release(), d);
}

ex eval_diff_Csch(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sCoth.copy(), ecopychild(e,1));
    return ex_mulx(emake_cint(-1), f.release(), d, ecopy(e));
}

ex eval_diff_Sech(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sTanh.copy(), ecopychild(e,1));
    return ex_mulx(emake_cint(-1), f.release(), ecopy(e), d);
}

ex eval_diff_Coth(ex F, er e)
{
    ex d;
    uex f(F);
    d = emake_node(gs.sym_sCsch.copy(), ecopychild(e,1));
    d = emake_node(gs.sym_sPower.copy(), d, emake_cint(2));
    return ex_mulx(emake_cint(-1), f.release(), d);
}


/* return a + b*E^i */
ex _arctrig_helper(ex E, slong a, slong b, slong i)
{
    if (i != 1)
        E = ex_powx(E, emake_cint(i));
    if (b != 1)
        E = ex_mulx(emake_cint(b), E);
    if (a != 0)
        E = ex_addx(emake_cint(a), E);
    return E;
}

ex eval_diff_ArcSin(ex F, er e)
{
    ex d;
    uex f(F);
    d = _arctrig_helper(ecopychild(e,1), 1, -1, 2);
    d = emake_node(gs.sym_sPower.copy(), d, emake_crat(-1,2));
    return ex_mulx(f.release(), d);
}

ex eval_diff_ArcCos(ex F, er e)
{
    ex d;
    uex f(F);
    d = _arctrig_helper(ecopychild(e,1), 1, -1, 2);
    d = emake_node(ecopy(gs.sym_sPower.get()), d, emake_crat(-1,2));
    return ex_mulx(emake_cint(-1), f.release(), d);
}

ex eval_diff_ArcTan(ex F, er e)
{
    ex d;
    uex f(F);
    d = _arctrig_helper(ecopychild(e,1), 1, 1, 2);
    d = ex_reciprocal(d);
    return ex_mulx(f.release(), d);
}

ex eval_diff_ArcCsc(ex F, er e)
{
    ex d;
    uex f(F);
    uex g(_arctrig_helper(ecopychild(e,1), 0, 1, -2));
    d = _arctrig_helper(ecopychild(e,1), 1, -1, -2);
    d = emake_node(gs.sym_sPower.copy(), d, emake_crat(-1,2));
    return ex_mulx(emake_cint(-1), f.release(), g.release(), d);
}

ex eval_diff_ArcSec(ex F, er e)
{
    ex d;
    uex f(F);
    uex g(_arctrig_helper(ecopychild(e,1), 0, 1, -2));
    d = _arctrig_helper(ecopychild(e,1), 1, -1, -2);
    d = emake_node(gs.sym_sPower.copy(), d, emake_crat(-1,2));
    return ex_mulx(f.release(), g.release(), d);
}

ex eval_diff_ArcCot(ex F, er e)
{
    ex d;
    uex f(F);
    d = _arctrig_helper(ecopychild(e,1), 1, 1, 2);
    d = ex_reciprocal(d);
    return ex_mulx(emake_cint(-1), f.release(), d);
}



ex eval_diff_ArcSinh(ex F, er e)
{
    ex d;
    uex f(F);
    d = _arctrig_helper(ecopychild(e,1), 1, 1, 2);
    d = emake_node(gs.sym_sPower.copy(), d, emake_crat(-1,2));
    return ex_mulx(f.release(), d);
}

ex eval_diff_ArcCosh(ex F, er e)
{
    uex f(F);
    uex d1(_arctrig_helper(ecopychild(e,1), -1, 1, 1));
    uex d2(_arctrig_helper(ecopychild(e,1), 1, 1, 1));
    d1.setz(ex_powx(d1.release(), emake_crat(-1,2)));
    d2.setz(ex_powx(d2.release(), emake_crat(-1,2)));
    return ex_mulx(f.release(), d1.release(), d2.release());
}

ex eval_diff_ArcTanh(ex F, er e)
{
    ex d;
    uex f(F);
    d = _arctrig_helper(ecopychild(e,1), 1, -1, 2);
    d = ex_reciprocal(d);
    return ex_mulx(f.release(), d);
}

ex eval_diff_ArcCsch(ex F, er e)
{
    ex d;
    uex f(F);
    uex g(_arctrig_helper(ecopychild(e,1), 0, 1, -2));
    d = _arctrig_helper(ecopychild(e,1), 1, 1, -2);
    d = ex_powx(d, emake_crat(-1,2));
    return ex_mulx(emake_cint(-1), f.release(), g.release(), d);
}

ex eval_diff_ArcSech(ex F, er e)
{
    ex d;
    uex f(F);
    uex d1(_arctrig_helper(ecopychild(e,1), 1, +1, 1));
    uex d2(_arctrig_helper(ecopychild(e,1), 1, -1, 1));
    d2.setz(ex_divx(d2.release(), d1.copy()));
    d2.setz(ex_powx(d2.release(), emake_crat(-1,2)));
    d1.setz(ex_reciprocal(d1.release()));
    ex g = ex_reciprocal(ecopychild(e,1));
    return ex_mulx(f.release(), g, d2.release(), d1.release());
}

ex eval_diff_ArcCoth(ex F, er e)
{
    ex d;
    uex f(F);
    d = _arctrig_helper(ecopychild(e,1), 1, -1, 2);
    d = ex_reciprocal(d);
    return ex_mulx(emake_cint(-1), f.release(), d);
}

ex eval_diff(er e, er var);

void _diff_arg(std::vector<wex> &v, er e, er child, std::vector<size_t> &idx, er var)
{
	if (ehas_head_sym(child, gs.sym_sList.get()))
	{
		for (size_t i = 1; i <= elength(child); i++)
		{
			idx.push_back(i);
			idx.push_back(elength(child));
			_diff_arg(v, e, echild(child,i), idx, var);
			idx.pop_back();
			idx.pop_back();
		}
	}
	else
	{
		uex f(eval_diff(child, var));
		if (!eis_zero(f.get()))
		{
			assert(idx.size() >= 2);
			assert((idx.size() % 2) == 0);
			uex cur(emake_cint(1));
			for (size_t j = idx.size(); j >= 2; j -= 2)
			{
				uex new_cur; new_cur.init_push_backr(j > 2 ? gs.sym_sList.get() : gs.sym_sDerivative.get(), idx[j - 1]);
				for (size_t k = 1; k <= idx[j - 1]; k++)
					new_cur.push_back(k == idx[j - 2] ? cur.copy() : emake_cint(0));
				cur.setnz(new_cur.release());
			}
			uex ne(ecopy(e));
			er e0 = echild(e,0);
			if (eis_node(e0) && elength(e0) == 1 &&
                ehas_head_sym_length(echild(e0,0), gs.sym_sDerivative.get(), elength(cur.get())))
			{
				uex new_cur; new_cur.init_push_backr(gs.sym_sDerivative.get(), elength(cur.get()));
				for (size_t k = 1; k <= elength(cur.get()); k++)
					new_cur.push_back(ex_addr(cur.child(k), echild(e0,0,k)));
				ne.replacechild(0, emake_node(new_cur.release(), ecopychild(e0,1)));
			}
			else
			{
				ne.replacechild(0, emake_node(cur.release(), ecopy(e0)));
			}
			v.push_back(ex_mulx(f.release(), ne.release()));
		}
	}
}

ex eval_diff(er e, er var)
{
//std::cout << "eval_diff: e: " << ex_tostring_full(e) << " var: " << ex_tostring_full(var) << std::endl;

    if (!eis_node(e))
    {
        return emake_cint(eis_sym(e, var) ? 1 : 0);
    }
	if (ex_same(e, var))
	{
        return emake_cint(1);
	}
    size_t n = elength(e);
    er h = echild(e,0);
    if (h == gs.sym_sPlus.get())
    {
        uex d; d.init_push_backr(gs.sym_sPlus.get(), n);
        for (size_t i = 1; i <= n; i++)
            d.push_back(eval_diff(echild(e,i),var));
        assert(evalid_bucketflags(d.get()));
        return ex_canonicalize_plus(d.release());
    }
    else if (h == gs.sym_sTimes.get())
    {
        uex d; d.init_push_backr(gs.sym_sPlus.get(), n);
        assert(evalid_bucketflags(d.get()));
        for (size_t i = 1; i <= n; i++)
        {
            uex f(ecopy(e));
            f.replacechild(i, eval_diff(echild(e,i), var));
            assert(evalid_bucketflags(echild(f.get(),i)));
            assert(evalid_bucketflags(f.get()));
            d.push_back(ex_canonicalize_times(f.release()));
        }
        assert(evalid_bucketflags(d.get()));
        return ex_canonicalize_plus(d.release());
    }
    else if (n == 2 && h == gs.sym_sPower.get())
    {
        er a = echild(e,1);
        er b = echild(e,2);
        uex y(eval_diff(b, var));
        if (eis_zero(y.get()))
        {
            y.reset(_arctrig_helper(ecopy(b), -1, 1, 1));
            y.reset(ex_powx(ecopy(a), y.release()));
            ex x = eval_diff(a,var);
            return ex_mulx(ecopy(b), x, y.release());
        }
        else
        {
            ex x = ex_logr(a);
            y.setz(ex_mulx(y.release(), x));
            uex t1(ex_reciprocal(ecopy(a)));
            ex t2 = eval_diff(a,var);
            x = ex_mulx(ecopy(b), t2, t1.release());
            return ex_mulx(ecopy(e), ex_addx(x, y.release()));
        }
    }
    else if (h == gs.sym_sList.get())
    {
        uex d; d.init_push_backr(gs.sym_sList.get(), n);
        for (ulong i = 0; i < n; i++)
            d.push_back(eval_diff(echild(e,i+1), var));
        return d.release();
    }
    else if (h == gs.sym_sMinus.get())
    {
        uex d; d.init_push_backr(gs.sym_sPlus.get(), n);
        for (ulong i = 0; i < n; i++)
        {
            ex f = eval_diff(echild(e,i+1),var);
            d.push_back(i == 0 ? f : ex_negate(f));
        }
        return ex_canonicalize_plus(d.release());
    }
    else if (n == 1 && ehas_head_sym(h, gs.sym_sDerivative.get()))
    {
        ex f = eval_diff(echild(e,1),var);
        if (eis_zero(f))
            return f;
        else
            return emake_node(ecopy(h), f);
    }
    else if (n == 1)
    {
        ex (*df)(ex, er) = nullptr;
             if (false) df = nullptr;
#define CHECK(name) else if (h == gs.sym_s##name.get()) df = eval_diff_##name;
        CHECK(Log)
        CHECK(Exp)
        CHECK(Sin)
        CHECK(Cos)
        CHECK(Tan)
        CHECK(Csc)
        CHECK(Sec)
        CHECK(Cot)
        CHECK(Sinh)
        CHECK(Cosh)
        CHECK(Tanh)
        CHECK(Csch)
        CHECK(Sech)
        CHECK(Coth)
        CHECK(ArcSin)
        CHECK(ArcCos)
        CHECK(ArcTan)
        CHECK(ArcCsc)
        CHECK(ArcSec)
        CHECK(ArcCot)
        CHECK(ArcSinh)
        CHECK(ArcCosh)
        CHECK(ArcTanh)
        CHECK(ArcCsch)
        CHECK(ArcSech)
        CHECK(ArcCoth)
#undef CHECK

        if (df != nullptr)
        {
            ex F = eval_diff(echild(e,1), var);
            if (eis_zero(F))
                return F;
            else
                return df(F, e);
        }
    }

	std::vector<wex> v;
	ex hd = eval_diff(echild(e,0), var);


	if (eis_zero(hd))
	{
		eclear(hd);
	}
	else
	{
		wex t(ecopy(e));
		t.replacechild(0, hd);
		v.push_back(t);

std::cout << "initial push back: " << ex_tostring(v.back().get()) << std::endl;

	}

	std::vector<size_t> idx;
	for (size_t i = 1; i <= elength(e); i++)
	{
		idx.clear();
		idx.push_back(i);
		idx.push_back(elength(e));
		_diff_arg(v, e, echild(e,i), idx, var);
	}

	return ex_canonicalize_plus(emake_node(gs.sym_sPlus.copy(), v));
}


ex dcode_sD(er e)
{
//std::cout << "dcode_sD: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sD.get()));

	if (elength(e) == 0)
		return _handle_message_argm(e, 1);

	uex c(ecopychild(e,1));
	for (ulong i = 1; i < elength(e); i++)
	{
		if (ehas_head_sym(echild(e,i+1), gs.sym_sList.get()))
		{
			return ecopy(e);
		}
		else
		{
			c.setnz(eval_diff(c.get(), echild(e,i+1)));
		}
	}
	return c.release();
}
