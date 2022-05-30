#include <cmath>
#include <cfloat>

#include "uex.h"
#include "timing.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "polynomial.h"

#include <flint/fmpq_mpoly_factor.h>
#include <flint/fmpz.h>

#include "djrat.h"

ex ex_numerator(er e);

ex ex_denominator(er e);

ex eval_fakeabs(bool&switched, ex E)
{
    uex e(E);

    switched = false;
    if (eis_number(E))
    {
        if ((eis_int(E) && fmpz_sgn(eint_data(E)) < 0) ||
            (eis_rat(E) && fmpq_sgn(erat_data(E)) < 0))
        {
            switched = true;
            e.setnz(num_Minus(e.get()));
        }
        return e.release();
    }

    if (ehas_head_sym(E, gs.sym_sTimes.get()) &&
        elength(E) > 1 &&
        eis_number(echild(E,1)))
    {
        er F = echild(E,1);
        if ((eis_int(F) && fmpz_sgn(eint_data(F)) < 0) ||
            (eis_rat(F) && fmpq_sgn(erat_data(F)) < 0))
        {
            switched = true;
            if (eis_int(F, -1))
            {
                if (elength(E) > 2)
                    e.removechild(1);
                else
                    e.setnz(e.copychild(2));
            }
            else
            {
                e.replacechild(1, num_Minus(F));
            }
        }
    }

    return e.release();
}


static void split_numerator_denominator(std::vector<wex> & a, std::vector<wex> & b, er e)
{
    if (ehas_head_sym(e, gs.sym_sTimes.get()))
    {
        for (size_t i = 1; i <= elength(e); i++)
            split_numerator_denominator(a, b, echild(e,i));
    }
    else if (ehas_head_sym_length(e, gs.sym_sPower.get(), 2))
    {
        bool negative = false;
        if (ehas_head_sym(echild(e,2), gs.sym_sPlus.get()))
        {
            for (size_t i = 1; i <= elength(echild(e,2)); i++)
            {
                ex newpow = eval_fakeabs(negative, ecopychild(e,2,i));
                (negative ? b : a).push_back(emake_node(gs.sym_sPower.copy(), ecopychild(e,1), newpow));
            }
        }
        else
        {
            ex newpow = eval_fakeabs(negative, ecopychild(e,2));
            if (negative)
            {
                b.push_back(emake_node(gs.sym_sPower.copy(), ecopychild(e,1), newpow));
            }
            else
            {
                eclear(newpow);
                a.push_back(ecopy(e));
            }
        }
    }
    else if (eis_int(e))
    {
        a.push_back(ecopy(e));
    }
    else if (eis_rat(e))
    {
        b.push_back(emake_int_copy(fmpq_denref(erat_data(e))));
        if (!fmpz_is_one(fmpq_numref(erat_data(e))))
            a.push_back(emake_int_copy(fmpq_numref(erat_data(e))));
    }
    else if (eis_cmplx(e))
    {
        uex rd(ex_denominator(ecmplx_real(e)));
        uex id(ex_denominator(ecmplx_imag(e)));

        if (eis_int(rd.get()) && eis_int(id.get()))
        {
            ex z = emake_int();
            fmpz_lcm(eint_data(z), rd.int_data(), id.int_data());
            uex cd(efix_int(z));
            if (!eis_one(cd.get()))
            {
                b.push_back(cd.copy());
                a.push_back(ex_mulr(cd.get(), e));
                return;
            }
        }

        a.push_back(ecopy(e));
    }
    else
    {
        a.push_back(ecopy(e));
    }
}


ex ex_numerator(er e)
{
    if (eis_int(e))
    {
        return ecopy(e);
    }
    else if (eis_rat(e))
    {
        return emake_int_copy(fmpq_numref(erat_data(e)));
    }
    else
    {
        std::vector<wex> num, den;
        split_numerator_denominator(num, den, e);
        return emake_node_times(num);
    }
}

ex ex_denominator(er e)
{
    if (eis_int(e))
    {
        return ecopy(e);
    }
    else if (eis_rat(e))
    {
        return emake_int_copy(fmpq_denref(erat_data(e)));
    }
    else
    {
        std::vector<wex> num, den;
        split_numerator_denominator(num, den, e);
        return emake_node_times(den);
    }
}


ex dcode_sNumerator(er e)
{
//std::cout << "dcode_sNumerator: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sNumerator.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ex_numerator(echild(e,1));
}

ex dcode_sDenominator(er e)
{
//std::cout << "dcode_sDenominator: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDenominator.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ex_denominator(echild(e,1));
}



ex fmpz_mpoly_to_ex(
    const fmpz_mpoly_t a,
    const fmpz_mpoly_ctx_t ctx,
    const std::vector<wex> vars)
{
    slong n = ctx->minfo->nvars;
    assert(n == vars.size());

    uex e; e.init_push_backr(gs.sym_sPlus.get(), fmpz_mpoly_length(a, ctx));

    xfmpz_t c;
    std::vector<xfmpz_t> exps;
    std::vector<fmpz *> expps;

    for (slong i = 0; i < n; i++)
        exps.push_back(xfmpz_t());
    for (slong i = 0; i < n; i++)
        expps.push_back(exps[i].data);

    for (slong i = fmpz_mpoly_length(a, ctx) - 1; i >= 0; i--)
    {
        fmpz_mpoly_get_term_coeff_fmpz(c.data, a, i, ctx);
        fmpz_mpoly_get_term_exp_fmpz(expps.data(), a, i, ctx);

        slong l = 1;
        for (slong j = 0; j < n; j++)
            l += !fmpz_is_zero(expps[j]);

        uex t; t.init_push_backr(gs.sym_sTimes.get(), l);

        t.push_back(emake_int_move(c));
        for (slong j = 0; j < n; j++)
        {
            if (fmpz_is_zero(expps[j]))
                continue;
            else if (fmpz_is_one(expps[j]))
                t.push_back(vars[j].copy());
            else
                t.push_back(ex_powx(vars[j].copy(), emake_int_move(expps[j])));
        }
        e.push_back(ex_canonicalize_times(t.release()));
    }

    return eupdate_timestamp(ex_canonicalize_plus(e.release()));
}

ex nmod_mpoly_to_ex(const nmod_mpoly_t a, const nmod_mpoly_ctx_t ctx, const std::vector<wex> vars)
{
    slong n = ctx->minfo->nvars;
    assert(n == vars.size());

    uex e; e.init_push_backr(gs.sym_sPlus.get(), a->length);

    std::vector<xfmpz_t> exps;
    std::vector<fmpz *> expps;
    for (slong i = 0; i < n; i++)
        exps.push_back(xfmpz_t());
    for (slong i = 0; i < n; i++)
        expps.push_back(exps[i].data);

    for (slong i = nmod_mpoly_length(a, ctx) - 1; i >= 0; i--)
    {
        nmod_mpoly_get_term_exp_fmpz(expps.data(), a, i, ctx);

        slong l = 1;
        for (slong j = 0; j < n; j++)
            l += !fmpz_is_zero(expps[j]);

        uex t; t.init_push_backr(gs.sym_sTimes.get(), l);

        ulong c = nmod_mpoly_get_term_coeff_ui(a, i, ctx);
        t.push_back(emake_int_ui(c));
        for (slong j = 0; j < n; j++)
        {
            if (fmpz_is_zero(expps[j]))
                continue;
            else if (fmpz_is_one(expps[j]))
                t.push_back(vars[j].copy());
            else
                t.push_back(ex_powx(vars[j].copy(), emake_int_move(expps[j])));
        }
        e.push_back(ex_canonicalize_times(t.release()));
    }

    return eupdate_timestamp(ex_canonicalize_plus(e.release()));
}



ex fmpz_mod_mpoly_to_ex(
    const fmpz_mod_mpoly_t a,
    const fmpz_mod_mpoly_ctx_t ctx,
    const std::vector<wex> vars)
{
    slong n = ctx->minfo->nvars;
    assert(n == vars.size());

    uex e; e.init_push_backr(gs.sym_sPlus.get(), fmpz_mod_mpoly_length(a, ctx));

    xfmpz_t c;
    std::vector<xfmpz_t> exps;
    std::vector<fmpz *> expps;

    for (slong i = 0; i < n; i++)
        exps.push_back(xfmpz_t());
    for (slong i = 0; i < n; i++)
        expps.push_back(exps[i].data);

    for (slong i = fmpz_mod_mpoly_length(a, ctx) - 1; i >= 0; i--)
    {
        fmpz_mod_mpoly_get_term_coeff_fmpz(c.data, a, i, ctx);
        fmpz_mod_mpoly_get_term_exp_fmpz(expps.data(), a, i, ctx);

        slong l = 1;
        for (slong j = 0; j < n; j++)
            l += !fmpz_is_zero(expps[j]);

        uex t; t.init_push_backr(gs.sym_sTimes.get(), l);

        t.push_back(emake_int_move(c));
        for (slong j = 0; j < n; j++)
        {
            if (fmpz_is_zero(expps[j]))
                continue;
            else if (fmpz_is_one(expps[j]))
                t.push_back(vars[j].copy());
            else
                t.push_back(ex_powx(vars[j].copy(), emake_int_move(expps[j])));
        }
        e.push_back(ex_canonicalize_times(t.release()));
    }

    return eupdate_timestamp(ex_canonicalize_plus(e.release()));
}




void poly_perm_vars(poly & a, poly & b, std::vector<size_t> p)
{
    a.coeffs = b.coeffs;
    a.exps.clear();
    std::vector<xfmpz_t> t(a.vars.size());
    for (size_t i = 0; i < b.coeffs.size(); i++)
    {
        for (size_t j = 0; j < a.vars.size(); j++)
        {
            fmpz_zero(t[j].data);
        }

        for (size_t j = 0; j < b.vars.size(); j++)
        {
            fmpz_set(t[p[j]].data, b.exps[i*b.vars.size() + j].data);
        }

        for (size_t j = 0; j < a.vars.size(); j++)
        {
            a.exps.push_back(t[j]);
        }
    }
}


void poly_add(poly & a, poly & b, poly & c)
{
    a.vars.clear();

    std::vector<size_t> bperm;
    for (size_t i = 0; i < b.vars.size(); i++)
    {
        a.vars.push_back(wex(b.vars[i].copy()));
        bperm.push_back(i);
    }

    std::vector<size_t> cperm;
    for (size_t i = 0; i < c.vars.size(); i++)
    {
        for (size_t j = 0; j < a.vars.size(); j++)
        {
            if (ex_same(a.vars[j].get(), c.vars[i].get()))
            {
                cperm.push_back(j);
                break;
            }
        }
        if (cperm.size() <= i)
        {
            cperm.push_back(a.vars.size());
            a.vars.push_back(wex(c.vars[i].copy()));
        }
    }

    poly B(0), C(0);
    B.vars = a.vars;
    poly_perm_vars(B, b, bperm);
    poly_sort_terms(B);
    C.vars = a.vars;
    poly_perm_vars(C, c, cperm);
    poly_sort_terms(C);

    poly_add_vars_match(a, B, C);
}




void poly_mul(poly & a, poly & b, poly & c)
{
    a.vars.clear();

    std::vector<size_t> bperm;
    for (size_t i = 0; i < b.vars.size(); i++)
    {
        a.vars.push_back(wex(b.vars[i].copy()));
        bperm.push_back(i);
    }

    std::vector<size_t> cperm;
    for (size_t i = 0; i < c.vars.size(); i++)
    {
        for (size_t j = 0; j < a.vars.size(); j++)
        {
            if (ex_same(a.vars[j].get(), c.vars[i].get()))
            {
                cperm.push_back(j);
                break;
            }
        }
        if (cperm.size() <= i)
        {
            cperm.push_back(a.vars.size());
            a.vars.push_back(wex(c.vars[i].copy()));
        }
    }

    poly B(0), C(0);
    B.vars = a.vars;
    poly_perm_vars(B, b, bperm);
    poly_sort_terms(B);
    C.vars = a.vars;
    poly_perm_vars(C, c, cperm);
    poly_sort_terms(C);

    poly_mul_vars_match(a, B, C);
}

void poly_append(poly & b, poly & c)
{
    std::vector<wex> vars;

    std::vector<size_t> bperm;
    for (size_t i = 0; i < b.vars.size(); i++)
    {
        vars.push_back(wex(b.vars[i].copy()));
        bperm.push_back(i);
    }

    bool must_change = false;

    std::vector<size_t> cperm;
    for (size_t i = 0; i < c.vars.size(); i++)
    {
        for (size_t j = 0; j < vars.size(); j++)
        {
            if (ex_same(vars[j].get(), c.vars[i].get()))
            {
                cperm.push_back(j);
                break;
            }
        }
        if (cperm.size() <= i)
        {
            cperm.push_back(vars.size());
            vars.push_back(wex(c.vars[i].copy()));
            must_change = true;
        }
    }

    poly C(0);
    C.vars = vars;
    poly_perm_vars(C, c, cperm);

    if (must_change)
    {
//std::cout << "must_change" << std::endl;
        poly B(0);
        B.vars = vars;
        poly_perm_vars(B, b, bperm);
        b.swap(B);
    }

    assert(b.vars.size() == vars.size());
    assert(C.vars.size() == vars.size());

    for (size_t i = 0; i < C.coeffs.size(); i++)
    {
        b.coeffs.push_back(C.coeffs[i]);
        for (size_t j = 0; j < vars.size(); j++)
        {
            b.exps.push_back(C.exps[i*vars.size() + j]);
        }
    }
}


ex enumber_reduce_mod(er e, const xfmpz_t & m)
{
    assert(fmpz_sgn(m.data) >= 0);

    if (fmpz_is_zero(m.data))
        return ecopy(e);

    if (eis_number(e))
    {
        switch (etype(e))
        {
            case ETYPE_INT:
            {
                ex x = emake_int();
                fmpz_mod(eint_data(x), eint_data(e), m.data);
                return efix_int(x);
            }
            case ETYPE_RAT:
            {
                ex x = emake_int();
                if (fmpq_mod_fmpz(eint_data(x), erat_data(e), m.data))
                    return efix_int(x);
                eclear(x);
                return nullptr;
            }
            default:
            {
                return nullptr;
            }
        }
    }
    else
    {
        return ecopy(e);
    }
}

void monomial_swap(xfmpz_t * A, xfmpz_t * B, size_t nvars);

int poly_reduce_mod(poly & A, const xfmpz_t & m)
{
    size_t nvars = A.vars.size();
    xfmpz_t * Aexps = A.exps.data();
    size_t Alen = 0;
    for (size_t i = 0; i < A.coeffs.size(); i++)
    {
        ex x = enumber_reduce_mod(A.coeffs[i].get(), m);
        if (x == nullptr)
            return 1;
        if (eis_zero(x))
        {
            eclear(x);
        }
        else
        {
            A.coeffs[Alen].reset(x);
            if (Alen != i)
                monomial_swap(Aexps + Alen*nvars, Aexps + i*nvars, nvars);
            Alen++;
        }
    }

    A.coeffs.resize(Alen);
    A.exps.resize(Alen*nvars);
    return 0;
}

/* return 0 for success, nonzero for failure */
int poly_set_any_ex(poly & p, er e, const xfmpz_t & m)
{
    int r;

    if (!eis_node(e))
    {
        if (eis_number(e))
        {
            p.vars.clear();
            p.exps.clear();
            p.coeffs.clear();
            ex c = enumber_reduce_mod(e, m);
            if (c == nullptr)
                return 1;
            if (eis_zero(c))
                eclear(c);
            else
                p.coeffs.push_back(wex(c));
            return 0;
        }
        else
        {
            p.vars.clear();
            p.exps.clear();
            p.coeffs.clear();
            p.vars.push_back(wex(ecopy(e)));
            p.coeffs.push_back(wex(emake_cint(1)));
            p.exps.push_back(xfmpz_t(ulong(1)));
        }
//std::cout << " from " << ex_tostring(e) << " got " << p.tostring() << std::endl;
        return 0;
    }
    else if (ehas_head_sym(e, gs.sym_sPlus.get()))
    {
        p.vars.clear();
        p.exps.clear();
        p.coeffs.clear();
        poly q(0);
        for (size_t i = 1; i <= elength(e); i++)
        {
            r = poly_set_any_ex(q, echild(e,i), m);
            if (r != 0)
            {
//std::cout << " from " << ex_tostring(e) << " failed 5" << std::endl;
                return r;
            }
            poly_append(p, q);
        }

//timeit_t timer;
//timeit_start(timer);
        poly_sort_terms(p);
//timeit_stop(timer);
//std::cout << "sort time: " << timer->wall << std::endl;


//timeit_start(timer);
        poly_combine_like_terms(p);
//timeit_stop(timer);
//std::cout << "combine time: " << timer->wall << std::endl;

        r = poly_reduce_mod(p, m);
        if (r != 0)
        {
//std::cout << " from " << ex_tostring(e) << " failed 6" << std::endl;
            return r;
        }

//std::cout << " from " << ex_tostring(e) << " got " << p.tostring() << std::endl;
        return 0;
    }
    else if (ehas_head_sym(e, gs.sym_sTimes.get()))
    {
        p.vars.clear();
        p.exps.clear();
        p.coeffs.clear();
        p.coeffs.push_back(wex(emake_cint(1)));
        poly t(0), q(0);
        for (size_t i = 1; i <= elength(e); i++)
        {
            r = poly_set_any_ex(q, echild(e,i), m);
            if (r != 0)
            {
//std::cout << " from " << ex_tostring(e) << " failed 4" << std::endl;
                return r;
            }
            poly_mul(t, p, q);
            p.swap(t);
        }

        r = poly_reduce_mod(p, m);
        if (r != 0)
        {
//std::cout << " from " << ex_tostring(e) << " failed 7" << std::endl;
            return r;
        }


//std::cout << " from " << ex_tostring(e) << " got " << p.tostring() << std::endl;
        return 0;
    }
    else if (ehas_head_sym_length(e, gs.sym_sPower.get(), 2)
                && eis_int(echild(e,2))
                && fmpz_sgn(eint_data(echild(e,2))) >= 0)
    {
        poly q(0);
        r = poly_set_any_ex(q, echild(e,1), m);
        if (r != 0)
        {
//std::cout << " from " << ex_tostring(e) << " failed 3" << std::endl;
            return r;
        }
        r = poly_pow(p, q, fmpz_get_ui(eint_data(echild(e,2))));
        if (r != 0)
        {
//std::cout << " from " << ex_tostring(e) << " failed 2" << std::endl;
            return r;
        }


        r = poly_reduce_mod(p, m);
        if (r != 0)
        {
//std::cout << " from " << ex_tostring(e) << " failed 8" << std::endl;
            return r;
        }


//std::cout << " from " << ex_tostring(e) << " got " << p.tostring() << std::endl;
        return 0;
    }
    else
    {
//std::cout << " from " << ex_tostring(e) << " failed 1" << std::endl;
        return 1;
    }
}


/* return 0 for successs, nonzero for failure */
int poly_to_fmpq_mpoly(fmpq_mpoly_t a, const fmpq_mpoly_ctx_t ctx, const poly & b)
{
    size_t n = b.vars.size();
    assert(n == ctx->zctx->minfo->nvars);

    fmpz ** exps = new fmpz*[n];
    for (slong i = 0; i < n; i++)
    {
        exps[i] = new fmpz;
        fmpz_init(exps[i]);
    }

    fmpq_mpoly_resize(a, 0, ctx);
    for (size_t i = 0; i < b.coeffs.size(); i++)
    {
        for (size_t j = 0; j < n; j++)
            fmpz_set(exps[j], b.exps[n*i + j].data);

        er bi = b.coeffs[i].get();
        if (eis_int(bi))
        {
            fmpq_mpoly_push_term_fmpz_fmpz(a, eint_data(bi), exps, ctx);
        }
        else if (eis_rat(bi))
        {
            fmpq_mpoly_push_term_fmpq_fmpz(a, erat_data(bi), exps, ctx);
        }
        else
        {
            assert(false);
            return 1;
        }
    }

    for (slong i = 0; i < n; i++)
    {
        fmpz_clear(exps[i]);
        delete exps[i];
    }
    delete[] exps;

    fmpq_mpoly_reduce(a, ctx);
    assert(fmpq_mpoly_is_canonical(a, ctx));

    return 0;
}

/* return 0 for successs, nonzero for failure */
int poly_to_nmod_mpoly(nmod_mpoly_t a, const nmod_mpoly_ctx_t ctx, const poly & b)
{
    size_t n = b.vars.size();
    assert(n == ctx->minfo->nvars);

    fmpz ** exps = new fmpz*[n];
    for (slong i = 0; i < n; i++)
    {
        exps[i] = new fmpz;
        fmpz_init(exps[i]);
    }

    nmod_mpoly_resize(a, 0, ctx);
    for (size_t i = 0; i < b.coeffs.size(); i++)
    {
        for (size_t j = 0; j < n; j++)
            fmpz_set(exps[j], b.exps[n*i + j].data);

        er bi = b.coeffs[i].get();
        ulong c;
        if (eis_int(bi))
        {
            c = fmpz_fdiv_ui(eint_data(bi), ctx->mod.n);
        }
        else if (eis_rat(bi))
        {
            c = nmod_div(fmpz_fdiv_ui(fmpq_numref(erat_data(bi)), ctx->mod.n),
                         fmpz_fdiv_ui(fmpq_denref(erat_data(bi)), ctx->mod.n), ctx->mod);
        }
        else
        {
//std::cout << "poly_to_nmod_mpoly failing" << std::endl;
            return 1;
        }

        nmod_mpoly_push_term_ui_fmpz(a, c, exps, ctx);
    }

    for (slong i = 0; i < n; i++)
    {
        fmpz_clear(exps[i]);
        delete exps[i];
    }
    delete[] exps;

    assert(nmod_mpoly_is_canonical(a, ctx));

    return 0;
}









/*****************************************************************************/

static void split_base_power(wex & base, xfmpq_t & power)
{
    er b = base.get();
    if (!ehas_head_sym_length(b, gs.sym_sPower.get(), 2))
        return;

    er e = echild(b, 2);

    if (eis_int(e) && fmpz_sgn(eint_data(e)) > 0)
    {
        base.reset(ecopychild(b,1));
        fmpq_mul_fmpz(power.data, power.data, eint_data(e));
        return;
    }
    else if (eis_rat(e) && fmpq_sgn(erat_data(e)) > 0)
    {
        base.reset(ecopychild(b,1));
        fmpq_mul(power.data, power.data, erat_data(e));
        return;
    }
}


/* power can be negative */
static ex split_base_intpower(er b, fmpz_t power)
{
    if (ehas_head_sym_length(b, gs.sym_sPower.get(), 2))
    {
        er e = echild(b,2);

        if (eis_int(e))
        {
            fmpz_set(power, eint_data(e));
            return ecopychild(b,1);
        }
        else if (eis_rat(e))
        {
            xfmpq_t q(1,1);
            fmpz_set(power, fmpq_numref(erat_data(e)));
            fmpz_set(fmpq_denref(q.data), fmpq_denref(erat_data(e)));
            return emake_node(gs.sym_sPower.copy(), ecopychild(b,1), emake_rat_move(q.data));
        }
        else if (ehas_head_sym(e, gs.sym_sTimes.get())
                     && elength(e) >= 2
                     && (eis_rat(echild(e,1)) || eis_int(echild(e,1))))
        {
            uex f(ecopy(e));
            if (eis_int(echild(e,1)))
            {
                fmpz_set(power, eint_data(echild(e,1)));
                f.removechild(1);
                if (elength(f.get()) == 1)
                    f.reset(f.copychild(1));
            }
            else
            {
                xfmpq_t q(1,1);
                fmpz_set(power, fmpq_numref(erat_data(echild(e,1))));
                fmpz_set(fmpq_denref(q.data), fmpq_denref(erat_data(echild(e,1))));
                f.replacechild(1, emake_rat_move(q.data));
            }
            return emake_node(gs.sym_sPower.copy(), ecopychild(b,1), f.release());
        }
    }
    fmpz_one(power);
    return ecopy(b);
}

/* power must be >= 0 */
static ex split_base_posintpower(er b, fmpz_t power)
{
    if (ehas_head_sym_length(b, gs.sym_sPower.get(), 2))
    {
        er e = echild(b,2);

        if (eis_int(e))
        {
            if (fmpz_sgn(eint_data(e)) >= 0)
            {
                fmpz_set(power, eint_data(e));
                return ecopychild(b,1);
            }
            else
            {
                fmpz_neg(power, eint_data(e));
                return emake_node(gs.sym_sPower.copy(), ecopychild(b,1), emake_cint(-1));
            }
        }
        else if (eis_rat(e))
        {
            xfmpq_t q(1,1);
            fmpz_set(power, fmpq_numref(erat_data(e)));
            fmpz_set(fmpq_denref(q.data), fmpq_denref(erat_data(e)));
            if (fmpz_sgn(power) < 0)
            {
                fmpz_neg(power, power);
                fmpz_set_si(fmpq_numref(q.data), -1);
            }
            return emake_node(gs.sym_sPower.copy(), ecopychild(b,1), emake_rat_move(q.data));
        }
        else if (ehas_head_sym(e, gs.sym_sTimes.get())
                     && elength(e) >= 2
                     && (eis_rat(echild(e,1)) || eis_int(echild(e,1))))
        {
            uex f(ecopy(e));
            if (eis_int(echild(e,1)))
            {
                fmpz_set(power, eint_data(echild(e,1)));
                if (fmpz_sgn(power) >= 0)
                {
                    f.removechild(1);
                    if (elength(f.get()) == 1)
                        f.reset(f.copychild(1));
                }
                else
                {
                    fmpz_neg(power, power);
                    f.replacechild(1, emake_cint(-1));
                }
            }
            else
            {
                xfmpq_t q(1,1);
                fmpz_set(power, fmpq_numref(erat_data(echild(e,1))));
                fmpz_set(fmpq_denref(q.data), fmpq_denref(erat_data(echild(e,1))));
                if (fmpz_sgn(power) < 0)
                {
                    fmpz_neg(power, power);
                    fmpz_set_si(fmpq_numref(q.data), -1);
                }
                f.replacechild(1, emake_rat_move(q.data));
            }
            return emake_node(gs.sym_sPower.copy(), ecopychild(b,1), f.release());
        }
    }
    fmpz_one(power);
    return ecopy(b);
}

template <class R>
static void set_map(
    std::vector<typename R::xmpoly> & map,
    std::vector<wex> & org,
    std::vector<wex> & target,
    const typename R::mpoly_ctx ctx)
{
    map.clear();
    for (size_t i = 0; i < org.size(); i++)
    {
        wex Obase(org[i].copy());
        xfmpq_t Opow(1, 1);
        split_base_power(Obase, Opow);

        bool found = false;
        for (size_t j = 0; j < target.size(); j++)
        {
            wex Tbase(target[j].copy());
            xfmpq_t Tpow(1, 1);
            split_base_power(Tbase, Tpow);

            if (ex_same(Obase.get(), Tbase.get()))
            {
                xfmpq_t q;
                fmpq_div(q.data, Opow.data, Tpow.data);
                if (fmpz_is_one(fmpq_denref(q.data)))
                {
					map.push_back(typename R::xmpoly(ctx));
                    R::mpoly_gen(map.back().data, j, ctx);
                    R::mpoly_pow_fmpz(map.back().data, map.back().data, fmpq_numref(q.data), ctx);
                    found = true;
                    break;
                }
            }
        }
        assert(found);
    }
    assert(map.size() == org.size());
}

static void merge_vars(
    std::vector<wex> & Avars,
    const std::vector<wex> & Bvars,
    const std::vector<wex> & Cvars)
{
    Avars = Bvars;

    for (size_t i = 0; i < Cvars.size(); i++)
    {
        wex Cbase(Cvars[i].copy());
        xfmpq_t Cpow(1, 1);
        split_base_power(Cbase, Cpow);

        bool found = false;
        for (size_t j = 0; j < Avars.size(); j++)
        {
            wex Abase(Avars[j].copy());
            xfmpq_t Apow(1,1);
            split_base_power(Abase, Apow);

            if (ex_same(Abase.get(), Cbase.get()))
            {
                xfmpq_t g;
                fmpq_gcd(g.data, Cpow.data, Apow.data);
                if (fmpq_is_one(g.data))
                    Avars[j].reset(Abase.copy());
                else
                    Avars[j].reset(emake_node(gs.sym_sPower.copy(), Abase.copy(), emake_rat_move(g.data)));
                found = true;
                break;
            }
        }

        if (!found)
        {
            Avars.push_back(wex(Cvars[i].copy()));
        }
    }
}

static bool vars_match(
    const std::vector<wex> & Bvars,
    const std::vector<wex> & Cvars)
{
	size_t n = Cvars.size();
	if (Bvars.size() != n)
		return false;

	for (size_t i = 0; i < n; i++)
	{
		if (!ex_same(Bvars[i].get(), Cvars[i].get()))
			return false;
	}

	return true;
}


#define RATPOLY_FLAG_EXPAND 1
#define RATPOLY_FLAG_TRIG   2

class fmpq_ratpoly {
public:
    std::vector<wex> vars;
    fmpq_polyfactor data;
    fmpz_mpoly_ctx_t ctx;

    constexpr static ex (*mpoly_to_ex)(const fmpz_mpoly_t, const fmpz_mpoly_ctx_t, const std::vector<wex>) = &fmpz_mpoly_to_ex;

	constexpr static bool (*polyfactor_add_const)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpq_t c, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_add_fmpq;
    constexpr static bool (*polyfactor_add)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpq_polyfactor & c, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_add;
	constexpr static bool (*polyfactor_mul_const)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpq_t c, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_mul_fmpq;
    constexpr static bool (*polyfactor_mul)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpq_polyfactor & c, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_mul;
	constexpr static bool (*polyfactor_gcd_const)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpq_t c, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_gcd_fmpq;
    constexpr static bool (*polyfactor_gcd)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpq_polyfactor & c, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_gcd;
    constexpr static bool (*polyfactor_pow)(fmpq_polyfactor & a, const fmpz_t, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_pow;
    constexpr static bool (*polyfactor_factorize)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_factorize;
    constexpr static bool (*polyfactor_expand_numerator)(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpz_mpoly_ctx_t ctx) = &fmpq_polyfactor_expand_numerator;
    constexpr static void (*polyfactor_partial_fractions)(std::vector<fmpq_polyfactor> &, const fmpq_polyfactor &, slong, const fmpz_mpoly_ctx_t) = &fmpq_polyfactor_partial_fractions;

    constexpr static void (*mpoly_gen)(fmpz_mpoly_t, slong, const fmpz_mpoly_ctx_t) = &fmpz_mpoly_gen;
    constexpr static int (*mpoly_pow_fmpz)(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_t, const fmpz_mpoly_ctx_t) = &fmpz_mpoly_pow_fmpz;

	typedef fmpq_polyfactor polyfactor;
	typedef xfmpz_mpoly_t xmpoly;
	typedef fmpz_mpoly_ctx_t mpoly_ctx;

    fmpq_ratpoly()
    {
        fmpz_mpoly_ctx_init(ctx, 1, ORD_LEX);
    }

    fmpq_ratpoly(int)
    {
        fmpz_mpoly_ctx_init(ctx, 1, ORD_LEX);
    }

    fmpq_ratpoly(const fmpq_ratpoly & other)
    {
        fmpz_mpoly_ctx_init(ctx, std::max(size_t(1), other.vars.size()), ORD_LEX);
    }

    ~fmpq_ratpoly()
    {
        fmpz_mpoly_ctx_clear(ctx);
    }

    void set_ctx()
    {
        if (vars.empty() || ctx->minfo->nvars == vars.size())
            return;
        data.reset(ctx);
        fmpz_mpoly_ctx_clear(ctx);
        fmpz_mpoly_ctx_init(ctx, vars.size(), ORD_LEX);
    }

	void swap(fmpq_ratpoly & b)
	{
	    std::swap(vars, b.vars);
	    fmpq_polyfactor_swap(data, b.data);
	    fmpz_mpoly_ctx_struct t = *ctx;
	    *ctx = *b.ctx;
	    *b.ctx = t;
	}

    void set_fmpz(const fmpz_t a)
    {
        vars.clear();
        fmpz_set(fmpq_numref(data.sign), a);
        fmpz_one(fmpq_denref(data.sign));
        data.length = 0;
    }

    void set_fmpq(const fmpq_t a)
    {
        vars.clear();
        fmpq_set(data.sign, a);
        data.length = 0;
    }

    ex sign_to_ex() const
    {
        return emake_rat_copy(data.sign);
    }

    ex base_to_ex(slong i)
    {
        assert(i < data.length);
        return fmpz_mpoly_to_ex(data.base + i, ctx, vars);
    }

	void const_zero()
	{
        fmpq_zero(data.sign);
	}
	void const_one()
	{
        fmpq_one(data.sign);
	}
	void const_set(const fmpq_ratpoly & b)
	{
        fmpq_set(data.sign, b.data.sign);
	}
	void const_add(const fmpq_ratpoly & b, const fmpq_ratpoly & c)
	{
        fmpq_add(data.sign, b.data.sign, c.data.sign);
	}
	void const_gcd(const fmpq_ratpoly & b, const fmpq_ratpoly & c)
	{
        fmpq_gcd(data.sign, b.data.sign, c.data.sign);
	}
	void const_mul(const fmpq_ratpoly & b, const fmpq_ratpoly & c)
	{
        fmpq_mul(data.sign, b.data.sign, c.data.sign);
	}

	void map(fmpq_polyfactor & bm, fmpq_ratpoly & b)
	{
		std::vector<xfmpz_mpoly_t> Bmap;
		set_map<fmpq_ratpoly>(Bmap, b.vars, vars, ctx);
		fmpq_polyfactor_map(bm, ctx, b.data, b.ctx, Bmap);
	}
};

class nmod_ratpoly {
public:
    std::vector<wex> vars;
    nmod_polyfactor data;
    nmod_mpoly_ctx_t ctx;

    constexpr static ex (*mpoly_to_ex)(const nmod_mpoly_t, const nmod_mpoly_ctx_t, const std::vector<wex>) = &nmod_mpoly_to_ex;

	constexpr static bool (*polyfactor_add_const)(nmod_polyfactor & a, const nmod_polyfactor & b, mp_limb_t c, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_add_nmod;
    constexpr static bool (*polyfactor_add)(nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_polyfactor & c, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_add;
	constexpr static bool (*polyfactor_mul_const)(nmod_polyfactor & a, const nmod_polyfactor & b, mp_limb_t c, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_mul_nmod;
    constexpr static bool (*polyfactor_mul)(nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_polyfactor & c, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_mul;
	constexpr static bool (*polyfactor_gcd_const)(nmod_polyfactor & a, const nmod_polyfactor & b, mp_limb_t c, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_gcd_nmod;
    constexpr static bool (*polyfactor_gcd)(nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_polyfactor & c, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_gcd;
    constexpr static bool (*polyfactor_pow)(nmod_polyfactor & a, const fmpz_t, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_pow;
    constexpr static bool (*polyfactor_factorize)(nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_factorize;
    constexpr static bool (*polyfactor_expand_numerator)(nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_mpoly_ctx_t ctx) = &nmod_polyfactor_expand_numerator;

    constexpr static void (*mpoly_gen)(nmod_mpoly_t, slong, const nmod_mpoly_ctx_t) = &nmod_mpoly_gen;
    constexpr static int (*mpoly_pow_fmpz)(nmod_mpoly_t, const nmod_mpoly_t, const fmpz_t, const nmod_mpoly_ctx_t) = &nmod_mpoly_pow_fmpz;

	typedef nmod_polyfactor polyfactor;
	typedef xnmod_mpoly_t xmpoly;
	typedef nmod_mpoly_ctx_t mpoly_ctx;

    nmod_ratpoly(mp_limb_t modulus)
    {
        nmod_mpoly_ctx_init(ctx, 1, ORD_LEX, modulus);
    }

    nmod_ratpoly(const nmod_ratpoly & other)
    {
        nmod_mpoly_ctx_init(ctx, 1, ORD_LEX, other.ctx->mod.n);
    }

    ~nmod_ratpoly()
    {
        nmod_mpoly_ctx_clear(ctx);
    }

    void set_ctx()
    {
		mp_limb_t modulus = nmod_mpoly_ctx_modulus(ctx);
        if (vars.empty())
            return;
        data.reset(ctx);
        nmod_mpoly_ctx_clear(ctx);
        nmod_mpoly_ctx_init(ctx, vars.size(), ORD_LEX, modulus);
    }

	void swap(nmod_ratpoly & b)
	{
	    std::swap(vars, b.vars);
	    nmod_polyfactor_swap(data, b.data);
	    nmod_mpoly_ctx_struct t = *ctx;
	    *ctx = *b.ctx;
	    *b.ctx = t;
	}

    void set_fmpz(const fmpz_t a)
    {
        vars.clear();
		data.sign = fmpz_fdiv_ui(a, ctx->mod.n);
        data.length = 0;
    }

    void set_fmpq(const fmpq_t a)
    {
        vars.clear();
		data.sign = nmod_div(fmpz_fdiv_ui(fmpq_numref(a), ctx->mod.n),
							 fmpz_fdiv_ui(fmpq_denref(a), ctx->mod.n), ctx->mod);
        data.length = 0;
    }

    ex sign_to_ex() const
    {
        return emake_int_ui(data.sign);
    }

    ex base_to_ex(slong i)
    {
        assert(i < data.length);
        return nmod_mpoly_to_ex(data.base + i, ctx, vars);
    }

	void const_zero()
	{
		data.sign = 0;
	}
	void const_one()
	{
		data.sign = 1;
	}
	void const_set(const nmod_ratpoly & b)
	{
        data.sign = b.data.sign;
	}
	void const_add(const nmod_ratpoly & b, const nmod_ratpoly & c)
	{
        data.sign = nmod_add(b.data.sign, c.data.sign, ctx->mod);
	}
	void const_gcd(const nmod_ratpoly & b, const nmod_ratpoly & c)
	{
		data.sign = (b.data.sign != 0) || (c.data.sign != 0);
	}
	void const_mul(const nmod_ratpoly & b, const nmod_ratpoly & c)
	{
        data.sign = nmod_mul(b.data.sign, c.data.sign, ctx->mod);
	}

	void map(nmod_polyfactor & bm, nmod_ratpoly & b)
	{
		std::vector<xnmod_mpoly_t> Bmap;
		set_map<nmod_ratpoly>(Bmap, b.vars, vars, ctx);
		nmod_polyfactor_map(bm, ctx, b.data, b.ctx, Bmap);
	}
};


class fmpz_mod_ratpoly {
public:
    std::vector<wex> vars;
    fmpz_mod_polyfactor data;
    fmpz_mod_mpoly_ctx_t ctx;

    constexpr static ex (*mpoly_to_ex)(const fmpz_mod_mpoly_t, const fmpz_mod_mpoly_ctx_t, const std::vector<wex>) = &fmpz_mod_mpoly_to_ex;

	constexpr static bool (*polyfactor_add_const)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_add_fmpz_mod;
    constexpr static bool (*polyfactor_add)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_polyfactor & c, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_add;
	constexpr static bool (*polyfactor_mul_const)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_mul_fmpz_mod;
    constexpr static bool (*polyfactor_mul)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_polyfactor & c, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_mul;
	constexpr static bool (*polyfactor_gcd_const)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_gcd_fmpz_mod;
    constexpr static bool (*polyfactor_gcd)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_polyfactor & c, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_gcd;
    constexpr static bool (*polyfactor_pow)(fmpz_mod_polyfactor & a, const fmpz_t, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_pow;
    constexpr static bool (*polyfactor_factorize)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_factorize;
    constexpr static bool (*polyfactor_expand_numerator)(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_mpoly_ctx_t ctx) = &fmpz_mod_polyfactor_expand_numerator;

    constexpr static void (*mpoly_gen)(fmpz_mod_mpoly_t, slong, const fmpz_mod_mpoly_ctx_t) = &fmpz_mod_mpoly_gen;
    constexpr static int (*mpoly_pow_fmpz)(fmpz_mod_mpoly_t, const fmpz_mod_mpoly_t, const fmpz_t, const fmpz_mod_mpoly_ctx_t) = &fmpz_mod_mpoly_pow_fmpz;

	typedef fmpz_mod_polyfactor polyfactor;
	typedef xfmpz_mod_mpoly_t xmpoly;
	typedef fmpz_mod_mpoly_ctx_t mpoly_ctx;

    fmpz_mod_ratpoly(const fmpz_t modulus)
    {
        fmpz_mod_mpoly_ctx_init(ctx, 1, ORD_LEX, modulus);
    }

    fmpz_mod_ratpoly(const fmpz_mod_ratpoly & other)
    {
        fmpz_mod_mpoly_ctx_init(ctx, 1, ORD_LEX, fmpz_mod_mpoly_ctx_modulus(other.ctx));
    }

    ~fmpz_mod_ratpoly()
    {
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    void set_ctx()
    {
        if (vars.empty())
            return;
        data.reset(ctx);
        mpoly_ctx_clear(ctx->minfo);
        mpoly_ctx_init(ctx->minfo, vars.size(), ORD_LEX);
    }

	void swap(fmpz_mod_ratpoly & b)
	{
	    std::swap(vars, b.vars);
	    fmpz_mod_polyfactor_swap(data, b.data);
	    fmpz_mod_mpoly_ctx_struct t = *ctx;
	    *ctx = *b.ctx;
	    *b.ctx = t;
	}

    void set_fmpz(const fmpz_t a)
    {
        vars.clear();
        fmpz_mod_set_fmpz(data.sign, a, ctx->ffinfo);
        data.length = 0;
    }

    void set_fmpq(const fmpq_t a)
    {
        xfmpz_t anum, aden;
        vars.clear();
        fmpz_mod_set_fmpz(anum.data, fmpq_numref(a), ctx->ffinfo);
        fmpz_mod_set_fmpz(aden.data, fmpq_denref(a), ctx->ffinfo);
        fmpz_mod_divides(data.sign, anum.data, aden.data, ctx->ffinfo);
        data.length = 0;
    }

    ex sign_to_ex() const
    {
        return emake_int_copy(data.sign);
    }

    ex base_to_ex(slong i)
    {
        assert(i < data.length);
        return fmpz_mod_mpoly_to_ex(data.base + i, ctx, vars);
    }

	void const_zero()
	{
        fmpz_zero(data.sign);
	}
	void const_one()
	{
        fmpz_one(data.sign);
	}
	void const_set(const fmpz_mod_ratpoly & b)
	{
        fmpz_set(data.sign, b.data.sign);
	}
	void const_add(const fmpz_mod_ratpoly & b, const fmpz_mod_ratpoly & c)
	{
        fmpz_mod_add(data.sign, b.data.sign, c.data.sign, ctx->ffinfo);
	}
	void const_gcd(const fmpz_mod_ratpoly & b, const fmpz_mod_ratpoly & c)
	{
        fmpz_set_ui(data.sign, !fmpz_is_zero(b.data.sign) ||
                               !fmpz_is_zero(c.data.sign));
	}
	void const_mul(const fmpz_mod_ratpoly & b, const fmpz_mod_ratpoly & c)
	{
        fmpz_mod_mul(data.sign, b.data.sign, c.data.sign, ctx->ffinfo);
	}

	void map(fmpz_mod_polyfactor & bm, fmpz_mod_ratpoly & b)
	{
		std::vector<xfmpz_mod_mpoly_t> Bmap;
		set_map<fmpz_mod_ratpoly>(Bmap, b.vars, vars, ctx);
		fmpz_mod_polyfactor_map(bm, ctx, b.data, b.ctx, Bmap);
	}
};

template <class R>
ex ratpoly_get_ex(R & p)
{
    uex r; r.init_push_backr(gs.sym_sTimes.get(), p.data.length + 1);
    r.push_back(p.sign_to_ex());
    if (!p.vars.empty() && p.data.length > 0)
    {
        for (slong i = 0; i < p.data.length; i++)
        {
            uex t1(p.base_to_ex(i));
            ex t2 = emake_int_copy(p.data.exp + i);
            r.push_back(ex_powx(t1.release(), t2));
        }
    }
    return eupdate_timestamp(ex_canonicalize_times(r.release()));
}

template <class R>
ex ratpoly_get_ex_list(R & p)
{
    uex r; r.init_push_backr(gs.sym_sList.get(), p.data.length + 1);
    r.push_back(emake_list(p.sign_to_ex(), emake_cint(1)));
    if (!p.vars.empty() && p.data.length > 0)
    {
        for (slong i = 0; i < p.data.length; i++)
        {
            uex t1(p.base_to_ex(i));
            ex t2 = emake_int_copy(p.data.exp + i);
            r.push_back(emake_list(t1.release(), t2));
        }
    }
    return eupdate_timestamp(r.release());
}


template <class R>
std::string ratpoly_tostring(R & p)
{
    wex e(ratpoly_get_ex<R>(p));
    return ex_tostring(e.get());
}

template<class R>
void ratpoly_add(R & a, R & b, R & c)
{
//std::cout << "----------------------" << std::endl;
//std::cout << "ratpoly_add b: " << ratpoly_tostring<R>(b) << std::endl;
//std::cout << "            c: " << ratpoly_tostring<R>(c) << std::endl;

    if (c.vars.empty())
    {
        if (b.vars.empty())
        {
            a.vars.clear();
            a.const_add(b, c);
        }
        else
        {
            a.vars = b.vars;
            a.set_ctx();
            R::polyfactor_add_const(a.data, b.data, c.data.sign, a.ctx);
        }
    }
    else if (b.vars.empty())
    {
        a.vars = c.vars;
        a.set_ctx();
        R::polyfactor_add_const(a.data, c.data, b.data.sign, a.ctx);
    }
	else if (vars_match(b.vars, c.vars))
	{
		a.vars = b.vars;
	    a.set_ctx();
	    R::polyfactor_add(a.data, b.data, c.data, a.ctx);		
	}
	else
	{
	    merge_vars(a.vars, b.vars, c.vars);
	    a.set_ctx();
	    typename R::polyfactor bm, cm;
		a.map(bm, b);
		a.map(cm, c);
	    R::polyfactor_add(a.data, bm, cm, a.ctx);
	}

//std::cout << "ratpoly_add a: " << ratpoly_tostring<R>(a) << std::endl;
}

template <class R>
void ratpoly_gcd(R & a, R & b, R & c)
{
//std::cout << "ratpoly_gcd b: " << ratpoly_tostring<R>(b) << std::endl;
//std::cout << "            c: " << ratpoly_tostring<R>(c) << std::endl;

    if (c.vars.empty())
    {
        if (b.vars.empty())
        {
            a.vars.clear();
            a.const_gcd(b, c);
        }
        else
        {
            a.vars = b.vars;
            a.set_ctx();
            R::polyfactor_gcd_const(a.data, b.data, c.data.sign, a.ctx);
        }
    }
    else if (b.vars.empty())
    {
        a.vars = c.vars;
        a.set_ctx();
        R::polyfactor_gcd_const(a.data, c.data, b.data.sign, a.ctx);
    }
	else if (vars_match(b.vars, c.vars))
	{
		a.vars = b.vars;
	    a.set_ctx();
	    R::polyfactor_gcd(a.data, b.data, c.data, a.ctx);		
	}
	else
	{
	    merge_vars(a.vars, b.vars, c.vars);
	    a.set_ctx();
	    typename R::polyfactor bm, cm;
		a.map(bm, b);
		a.map(cm, c);
	    R::polyfactor_gcd(a.data, bm, cm, a.ctx);
	}

//std::cout << "ratpoly_gcd a: " << ratpoly_tostring<R>(a) << std::endl;
}

template <class R>
void ratpoly_mul(R & a, R & b, R & c)
{
//std::cout << "----------------------" << std::endl;
//std::cout << "ratpoly_mul b: " << ratpoly_tostring<R>(b) << std::endl;
//std::cout << "            c: " << ratpoly_tostring<R>(c) << std::endl;

    if (c.vars.empty())
    {
        if (b.vars.empty())
        {
            a.vars.clear();
			a.const_mul(b, c);
        }
        else
        {
            a.vars = b.vars;
            a.set_ctx();
            R::polyfactor_mul_const(a.data, b.data, c.data.sign, a.ctx);
        }
    }
    else if (b.vars.empty())
    {
        a.vars = c.vars;
        a.set_ctx();
        R::polyfactor_mul_const(a.data, c.data, b.data.sign, a.ctx);
    }
	else if (vars_match(b.vars, c.vars))
	{
		a.vars = b.vars;
	    a.set_ctx();
	    R::polyfactor_mul(a.data, b.data, c.data, a.ctx);		
	}
	else
	{
	    merge_vars(a.vars, b.vars, c.vars);
	    a.set_ctx();
	    typename R::polyfactor bm, cm;
		a.map(bm, b);
		a.map(cm, c);
	    R::polyfactor_mul(a.data, bm, cm, a.ctx);
	}

//std::cout << "ratpoly_mul a: " << ratpoly_tostring<R>(a) << std::endl;
}

template <class R>
void ratpoly_pow(R & a, const fmpz_t power)
{
//std::cout << "ratpoly_pow a: " << ratpoly_tostring<R>(a) << std::endl;
    R::polyfactor_pow(a.data, power, a.ctx);
//std::cout << "ratpoly_pow a: " << ratpoly_tostring<R>(a) << std::endl;
}

template <class R>
void ratpoly_factor(R & a, R & b)
{
//std::cout << "ratpoly_factor b: " << ratpoly_tostring<R>(b) << std::endl;

    if (b.vars.empty())
	{
		a.vars.clear();
		a.const_set(b);
	}
	else
	{
		a.vars = b.vars;
	    a.set_ctx();
	    R::polyfactor_factorize(a.data, b.data, a.ctx);
	}

//std::cout << "ratpoly_factor a: " << ratpoly_tostring<R>(b) << std::endl;
}


template <class R>
void ratpoly_expand_numerator(R & a, R & b)
{
//std::cout << "ratpoly_expand_numerator b: " << ratpoly_tostring<R>(b) << std::endl;

    if (b.vars.empty())
	{
		a.vars.clear();
		a.const_set(b);
	}
	else
	{
		a.vars = b.vars;
	    a.set_ctx();
	    R::polyfactor_expand_numerator(a.data, b.data, a.ctx);
	}

//std::cout << "ratpoly_expand_numerator a: " << ratpoly_tostring<R>(a) << std::endl;
}

void fmpq_polyfactor_print_pretty(const fmpq_polyfactor & a, const fmpz_mpoly_ctx_t ctx);

template <class R>
ex ratpoly_apart_ex(R & a)
{
    if (a.vars.empty())
        return ratpoly_get_ex<R>(a);

    std::vector<typename R::polyfactor> v;
    R::polyfactor_partial_fractions(v, a.data, 0, a.ctx);
    uex r; r.init_push_backr(gs.sym_sPlus.get(), v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        fmpq_polyfactor_swap(a.data, v[i]);
        r.push_back(ratpoly_get_ex<R>(a));
        fmpq_polyfactor_swap(a.data, v[i]);
    }

std::cout << "r: " << ex_tostring(r.get()) << std::endl;

    return r.release();
}

template <class R> void ratpoly_set_ex(
    R & p,
    er e,
    uint32_t flags,
    std::vector<wex> & gvars);

template <class R> void ratpoly_set_ex(R & p, er e, uint32_t flags)
{
    std::vector<wex> gvars;
    ratpoly_set_ex(p, e, flags, gvars);
}

template <class R> void ratpoly_set_ex_sum(
    R & p,
    er e,
    uint32_t flags,
    std::vector<wex> & gvars)
{
    std::vector<R> r(FLINT_BITS, R(p));

    r[0].const_zero();
    r[0].data.length = 0;
    int n = 1;

    for (ulong i = 0; i < elength(e); i++)
    {
        ratpoly_set_ex<R>(r[n], echild(e,i+1), flags, r[n-1].vars);
        n++;
        for (ulong j = i+2; j % 2 == 0 && n > 1; j /= 2, n--)
        {
            ratpoly_add<R>(r[n], r[n-1], r[n-2]);
            r[n-2].swap(r[n]);
        }
    }

    for (; n > 1; n--)
    {
        ratpoly_add<R>(r[n], r[n-1], r[n-2]);
        r[n-2].swap(r[n]);
    }

    p.swap(r[0]);
}

template <class R> void ratpoly_set_ex(
    R & p,
    er e,
    uint32_t flags,
    std::vector<wex> & gvars)
{
    if (!eis_node(e))
    {
        if (eis_number(e))
        {
            p.vars = gvars;
            p.set_ctx();
            if (eis_int(e))
            {
                p.set_fmpz(eint_data(e));
                return;
            }
            if (eis_rat(e))
            {
                p.set_fmpq(erat_data(e));
                return;
            }

            std::cout << "throw something here" << std::endl;
            abort();
        }
        else
        {
            for (size_t i = 0; i < gvars.size(); i++)
            {
                if (ex_same(e, gvars[i].get()))
                {
                    p.vars = gvars;
                    p.set_ctx();
                    p.const_one();
                    p.data.set_length(1, p.ctx);
                    R::mpoly_gen(p.data.base + 0, i, p.ctx);
                    fmpz_one(p.data.exp + 0);
                    return;
                }
            }

            p.vars.clear();
            p.vars.push_back(wex(ecopy(e)));
            p.set_ctx();
			p.const_one();
	        p.data.set_length(1, p.ctx);
	        R::mpoly_gen(p.data.base + 0, 0, p.ctx);
	        fmpz_one(p.data.exp + 0);
            return;
        }
    }
    else if (ehas_head_sym(e, gs.sym_sPlus.get()))
    {
        ratpoly_set_ex_sum<R>(p, e, flags, gvars);
        return;
    }
    else if (ehas_head_sym(e, gs.sym_sTimes.get()))
    {
        R q(p), t(p);
        p.vars = gvars;
        p.set_ctx();
		p.const_one();
        p.data.length = 0;
        for (size_t i = 1; i <= elength(e); i++)
        {
            ratpoly_set_ex<R>(q, echild(e,i), flags, p.vars);
            ratpoly_mul<R>(t, p, q);
            if ((flags & RATPOLY_FLAG_EXPAND) == 0)
                p.swap(t);
            else
                ratpoly_expand_numerator<R>(p, t);
        }
        return;
    }
    else if (ehas_head_sym_length(e, gs.sym_sPower.get(), 2)
                && eis_int(echild(e,2))
                && ((flags & RATPOLY_FLAG_EXPAND) == 0 ||
                    fmpz_sgn(eint_data(echild(e,2))) > 0) )
    {
        ratpoly_set_ex<R>(p, echild(e,1), flags, gvars);

        ratpoly_pow<R>(p, eint_data(echild(e,2)));

        if ((flags & RATPOLY_FLAG_EXPAND) != 0)
        {
            R t(p);
            ratpoly_expand_numerator<R>(t, p);
            p.swap(t);
        }
        return;
    }
    else
    {
        p.vars.clear();
        p.vars.push_back(wex(gs.sym_sNull.copy()));
        p.set_ctx();
        p.const_one();
        p.data.set_length(1, p.ctx);
        R::mpoly_gen(p.data.base + 0, 0, p.ctx);
        p.vars.back() = wex(((flags & RATPOLY_FLAG_EXPAND) == 0)
                                 ? split_base_intpower(e, p.data.exp + 0)
                                 : split_base_posintpower(e, p.data.exp + 0));
        return;
    }
}


bool getOption(wex & optvalue, er e, size_t idx, er optname)
{
    for (size_t i = idx + 1; i <= elength(e); i++)
    {
        if (ehas_head_sym_length(echild(e,i), gs.sym_sRule.get(), 2))
        {
            if (eis_sym(echild(e,i,1), optname))
            {
                optvalue.reset(ecopychild(e,i,2));
                return true;
            }
        }
        else
        {
            _gen_message(echild(e,0), "nonopt", NULL, ecopychild(e,i), emake_int_ui(idx), ecopy(e));
            return false;
        }
    }
    return true; 
}


ex dcode_sApart(er e)
{
//std::cout << "dcode_sApart: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sApart.get()));
    wex optModulus(emake_cint(0));
    wex optExtension(gs.sym_sNone.copy());
    wex optTrig(gs.sym_sFalse.copy());

    if (elength(e) < 1)
        return _handle_message_argx1(e);

    if (!getOption(optExtension, e, 1, gs.sym_sExtension.get()) ||
        !getOption(optModulus, e, 1, gs.sym_sModulus.get()) ||
        !getOption(optTrig, e, 1, gs.sym_sTrig.get()))
    {
        return ecopy(e);
    }

    fmpq_ratpoly p(0);
    ratpoly_set_ex<fmpq_ratpoly>(p, echild(e,1), 0);
    return ratpoly_apart_ex<fmpq_ratpoly>(p);
}

ex dcode_sTogether(er e)
{
//std::cout << "dcode_sTogether: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sTogether.get()));
    wex optModulus(emake_cint(0));
    wex optExtension(gs.sym_sNone.copy());
    wex optTrig(gs.sym_sFalse.copy());

    if (elength(e) < 1)
        return _handle_message_argx1(e);

    if (!getOption(optExtension, e, 1, gs.sym_sExtension.get()) ||
        !getOption(optModulus, e, 1, gs.sym_sModulus.get()) ||
        !getOption(optTrig, e, 1, gs.sym_sTrig.get()))
    {
        return ecopy(e);
    }

    if (eis_int(optModulus.get()))
    {
        xfmpz_t m(eint_number(optModulus.get()));
        fmpz_abs(m.data, m.data);
        if (!fmpz_is_zero(m.data) && (!fmpz_is_probabprime(m.data) || !fmpz_is_prime(m.data)))
        {
            _gen_message(echild(e,0), "modp", NULL, gs.sym_sModulus.copy(), optModulus.copy());
            return ecopy(e);
        }

        if (fmpz_is_zero(m.data))
        {
            fmpq_ratpoly p(0);
            ratpoly_set_ex<fmpq_ratpoly>(p, echild(e,1), 0);
            return ratpoly_get_ex<fmpq_ratpoly>(p);
        }
        else if (fmpz_abs_fits_ui(m.data))
        {
            nmod_ratpoly p(fmpz_get_ui(m.data));
            ratpoly_set_ex<nmod_ratpoly>(p, echild(e,1), 0);
            return ratpoly_get_ex<nmod_ratpoly>(p);
        }
        else
        {
            fmpz_mod_ratpoly p(m.data);
            ratpoly_set_ex<fmpz_mod_ratpoly>(p, echild(e,1), 0);
            return ratpoly_get_ex<fmpz_mod_ratpoly>(p);
        }
    }
    else
    {
        _gen_message(echild(e,0), "bmod", NULL, ecopychild(e,0));
        return ecopy(e);
    }
}

ex dcode_sCancel(er e)
{
//std::cout << "dcode_sCancel: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCancel.get()));
    wex optModulus(emake_cint(0));
    wex optExtension(gs.sym_sNone.copy());
    wex optTrig(gs.sym_sFalse.copy());

    if (elength(e) == 0)
    {
        return _handle_message_argx1(e);
    }

    if (!getOption(optExtension, e, 1, gs.sym_sExtension.get()) ||
        !getOption(optModulus, e, 1, gs.sym_sModulus.get()) ||
        !getOption(optTrig, e, 1, gs.sym_sTrig.get()))
    {
        return ecopy(e);
    }

    if (eis_int(optModulus.get()))
    {
        xfmpz_t m(eint_number(optModulus.get()));
        fmpz_abs(m.data, m.data);
        if (!fmpz_is_zero(m.data) && (!fmpz_is_probabprime(m.data) || !fmpz_is_prime(m.data)))
        {
            _gen_message(echild(e,0), "modp", NULL, gs.sym_sModulus.copy(), optModulus.copy());
            return ecopy(e);
        }

        if (fmpz_is_zero(m.data))
        {
            fmpq_ratpoly p(0);
            er f = echild(e,1);
            if (ehas_head_sym(f, gs.sym_sPlus.get()))
            {
                uex s; s.init_push_backr(gs.sym_sPlus.get(), elength(f));
                for (size_t i = 1; i <= elength(f); i++)
                {
                    ratpoly_set_ex<fmpq_ratpoly>(p, echild(f,i), 0);
                    s.push_back(ratpoly_get_ex<fmpq_ratpoly>(p));
                }
                return s.release();
            }
            else
            {
                ratpoly_set_ex<fmpq_ratpoly>(p, f, 0);
                return ratpoly_get_ex<fmpq_ratpoly>(p);
            }
        }
        else
        {
            _gen_message(echild(e,0), "priml", "Prime `1` is too large for this implementation.");
            return ecopy(e);
        }
    }
    else
    {
        _gen_message(echild(e,0), "bmod", NULL, ecopychild(e,0));
        return ecopy(e);
    }
}


ex dcode_sExpand(er e)
{
//std::cout << "dcode_sExpand: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sExpand.get()));
    wex optModulus(emake_cint(0));
    wex optExtension(gs.sym_sNone.copy());
    wex optTrig(gs.sym_sFalse.copy());

    if (elength(e) == 0)
    {
        return _handle_message_argx1(e);
    }

    if (!getOption(optExtension, e, 1, gs.sym_sExtension.get()) ||
        !getOption(optModulus, e, 1, gs.sym_sModulus.get()) ||
        !getOption(optTrig, e, 1, gs.sym_sTrig.get()))
    {
        return ecopy(e);
    }

    if (eis_int(optModulus.get()))
    {
        xfmpz_t m(eint_number(optModulus.get()));
        fmpz_abs(m.data, m.data);
        if (!fmpz_is_zero(m.data) && (!fmpz_is_probabprime(m.data) || !fmpz_is_prime(m.data)))
        {
            _gen_message(echild(e,0), "modp", NULL, gs.sym_sModulus.copy(), optModulus.copy());
            return ecopy(e);
        }

        if (fmpz_is_zero(m.data))
        {
            fmpq_ratpoly p(0);
            ratpoly_set_ex<fmpq_ratpoly>(p, echild(e,1), RATPOLY_FLAG_EXPAND);
            return ratpoly_get_ex<fmpq_ratpoly>(p);
        }
        else if (fmpz_abs_fits_ui(m.data))
        {
            nmod_ratpoly p(fmpz_get_ui(m.data));
            ratpoly_set_ex<nmod_ratpoly>(p, echild(e,1), RATPOLY_FLAG_EXPAND);
            return ratpoly_get_ex<nmod_ratpoly>(p);
        }
        else
        {
            fmpz_mod_ratpoly p(m.data);
            ratpoly_set_ex<fmpz_mod_ratpoly>(p, echild(e,1), RATPOLY_FLAG_EXPAND);
            return ratpoly_get_ex<fmpz_mod_ratpoly>(p);
        }
    }
    else
    {
        _gen_message(echild(e,0), "bmod", NULL, ecopychild(e,0));
        return ecopy(e);
    }
}


ex dcode_sFactorTerms(er e)
{
//std::cout << "dcode_sFactorTerms: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFactorTerms.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopy(e);
}


ex dcode_sPolynomialGCD(er e)
{
//std::cout << "dcode_sPolynomialGCD: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sPolynomialGCD.get()));
    wex optModulus(emake_cint(0));
    wex optExtension(gs.sym_sNone.copy());
    wex optTrig(gs.sym_sFalse.copy());

    if (elength(e) == 0)
    {
        return ecopy(e);
    }
/*
    if (!getOption(optExtension, e, 1, gs.sym_sExtension.get()) ||
        !getOption(optModulus, e, 1, gs.sym_sModulus.get()) ||
        !getOption(optTrig, e, 1, gs.sym_sTrig.get()))
    {
        return ecopy(e);
    }
*/
    size_t len = elength(e);

    if (eis_int(optModulus.get()))
    {
        xfmpz_t m(eint_number(optModulus.get()));
        fmpz_abs(m.data, m.data);
        if (!fmpz_is_zero(m.data) && (!fmpz_is_probabprime(m.data) || !fmpz_is_prime(m.data)))
        {
            _gen_message(echild(e,0), "modp", NULL, gs.sym_sModulus.copy(), optModulus.copy());
            return ecopy(e);
        }

        if (fmpz_is_zero(m.data))
        {
            fmpq_ratpoly q(0), p(0), t(0);
            ratpoly_set_ex<fmpq_ratpoly>(p, echild(e,1), 0);
            for (size_t i = 2; i <= len; i++)
            {
                ratpoly_set_ex<fmpq_ratpoly>(q, echild(e,i), 0);
                ratpoly_gcd<fmpq_ratpoly>(t, p, q);
                p.swap(t);
            }
            return ratpoly_get_ex<fmpq_ratpoly>(p);
        }
        else
        {
            _gen_message(echild(e,0), "priml", "Prime `1` is too large for this implementation.");
            return ecopy(e);
        }
    }
    else
    {
        _gen_message(echild(e,0), "bmod", NULL, ecopychild(e,0));
        return ecopy(e);
    }
}

ex dcode_sFactor(er e)
{
//std::cout << "dcode_sFactor: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFactor.get()));
    wex optModulus(emake_cint(0));
    wex optExtension(gs.sym_sNone.copy());
    wex optTrig(gs.sym_sFalse.copy());

    if (elength(e) == 0)
        return _handle_message_argx1(e);

    if (!getOption(optExtension, e, 1, gs.sym_sExtension.get()) ||
        !getOption(optModulus, e, 1, gs.sym_sModulus.get()) ||
        !getOption(optTrig, e, 1, gs.sym_sTrig.get()))
    {
        return ecopy(e);
    }

    if (eis_int(optModulus.get()))
    {
        xfmpz_t m(eint_number(optModulus.get()));
        fmpz_abs(m.data, m.data);
        if (!fmpz_is_zero(m.data) && !fmpz_is_probabprime(m.data))
        {
            _gen_message(echild(e,0), "modp", NULL, gs.sym_sModulus.copy(), optModulus.copy());
            return ecopy(e);
        }

        if (fmpz_is_zero(m.data))
        {
            fmpq_ratpoly q(0), p(0);
            ratpoly_set_ex<fmpq_ratpoly>(p, echild(e,1), 0);
            ratpoly_factor<fmpq_ratpoly>(q, p);
            return ratpoly_get_ex<fmpq_ratpoly>(q);
        }
        else if (fmpz_abs_fits_ui(m.data))
        {
            nmod_ratpoly q(fmpz_get_ui(m.data)), p(fmpz_get_ui(m.data));
            ratpoly_set_ex<nmod_ratpoly>(p, echild(e,1), 0);
            ratpoly_factor<nmod_ratpoly>(q, p);
            return ratpoly_get_ex<nmod_ratpoly>(q);
        }
        else
        {
            fmpz_mod_ratpoly q(m.data), p(m.data);
            ratpoly_set_ex<fmpz_mod_ratpoly>(p, echild(e,1), 0);
            ratpoly_factor<fmpz_mod_ratpoly>(q, p);
            return ratpoly_get_ex<fmpz_mod_ratpoly>(q);
        }
    }
    else
    {
        _gen_message(echild(e,0), "bmod", NULL, ecopychild(e,0));
        return ecopy(e);
    }
}

ex dcode_sFactorList(er e)
{
//std::cout << "dcode_sFactorList: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFactorList.get()));
    wex optModulus(emake_cint(0));
    wex optExtension(gs.sym_sNone.copy());
    wex optTrig(gs.sym_sFalse.copy());

    if (elength(e) == 0)
    {
        return _handle_message_argx1(e);
    }

    if (!getOption(optExtension, e, 1, gs.sym_sExtension.get()) ||
        !getOption(optModulus, e, 1, gs.sym_sModulus.get()) ||
        !getOption(optTrig, e, 1, gs.sym_sTrig.get()))
    {
        return ecopy(e);
    }

    if (eis_int(optModulus.get()))
    {
        xfmpz_t m(eint_number(optModulus.get()));
        fmpz_abs(m.data, m.data);
        if (!fmpz_is_zero(m.data) && (!fmpz_is_probabprime(m.data) || !fmpz_is_prime(m.data)))
        {
            _gen_message(echild(e,0), "modp", NULL, gs.sym_sModulus.copy(), optModulus.copy());
            return ecopy(e);
        }

        if (fmpz_is_zero(m.data))
        {
            fmpq_ratpoly q(0), p(0);
            ratpoly_set_ex<fmpq_ratpoly>(p, echild(e,1), 0);
            ratpoly_factor<fmpq_ratpoly>(q, p);
            return ratpoly_get_ex_list<fmpq_ratpoly>(q);
        }
        else if (fmpz_abs_fits_ui(m.data))
        {
            nmod_ratpoly q(fmpz_get_ui(m.data)), p(fmpz_get_ui(m.data));
            ratpoly_set_ex<nmod_ratpoly>(p, echild(e,1), 0);
            ratpoly_factor<nmod_ratpoly>(q, p);
            return ratpoly_get_ex_list<nmod_ratpoly>(q);
        }
        else
        {
            fmpz_mod_ratpoly q(m.data), p(m.data);
            ratpoly_set_ex<fmpz_mod_ratpoly>(p, echild(e,1), 0);
            ratpoly_factor<fmpz_mod_ratpoly>(q, p);
            return ratpoly_get_ex_list<fmpz_mod_ratpoly>(q);
        }
    }
    else
    {
        _gen_message(echild(e,0), "bmod", NULL, ecopychild(e,0));
        return ecopy(e);
    }
}

static ex _coefficient(er e, er vlist, const xfmpz_t * explist)
{
    if (ehas_head_sym(e, gs.sym_sPlus.get()))
    {
        uex t(emake_cint(0));
        for (size_t i = 0; i < elength(e); i++)
        {
            ex t2 = _coefficient(echild(e,i+1), vlist, explist);
            t.setz(ex_addx(t.release(), t2));
        }
        return t.release();
    }

    poly p(elength(vlist));
    if (!ex_to_polynomial(p, e, vlist))
        return emake_cint(0);

    ex t = p.get_coeff(explist);
    return t;
}

ex dcode_sCoefficient(er e)
{
//std::cout << "dcode_sCoefficient: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCoefficient.get()));

    if (elength(e) == 3 && eis_int(echild(e,3)))
    {
        std::vector<xfmpz_t> explist;
        explist.push_back(xfmpz_t(eint_data(echild(e,3))));
        wex vlist = emake_list(ecopychild(e,2));
        return _coefficient(echild(e,1), vlist.get(), explist.data());
    }

    return emake_cint(0);
}

static ex _exponent(er e, er v)
{
    if (!eis_node(e))
        return emake_cint(ex_same(e, v) ? 1 : 0);

    er h = echild(e,0);
    size_t n = elength(e);

    if (eis_sym(h, gs.sym_sPlus.get()) || eis_sym(h, gs.sym_sMinus.get()))
    {
        uex t(_exponent(echild(e,1), v));
        for (size_t i = 1; i < n; i++)
        {
            ex s = _exponent(echild(e,i+1), v);
            t.setz(ex_max(s, t.release()));
        }
        return t.release();
    }
    else if (n > 0 && eis_sym(h, gs.sym_sTimes.get()))
    {
        uex t(_exponent(echild(e,1), v));
        for (size_t i = 1; i < n; i++)
        {
            ex s = _exponent(echild(e,i+1), v);
            t.setz(ex_addx(s, t.release()));
        }
        return t.release();
    }
    else if (n > 0 && eis_sym(h, gs.sym_sDivide.get()))
    {
        uex t(_exponent(echild(e,1), v));
        for (size_t i = 1; i < n; i++)
        {
            ex s = _exponent(echild(e,i+1), v);
            t.setz(ex_subx(s, t.release()));
        }
        return t.release();
    }
    else if (n == 2 && eis_sym(h, gs.sym_sPower.get()))
    {
        uex t(_exponent(echild(e,1), v));
        return ex_mulx(t.release(), ecopychild(e,2));
    }
    else
    {
        return emake_cint(0);
    }
}

ex dcode_sExponent(er e)
{
//std::cout << "dcode_sExponent: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sExponent.get()));

    if (elength(e) != 2)
        return _handle_message_argx(e, 2);

    return _exponent(echild(e,1), echild(e,2));
}

ex dcode_sCoefficientRules(er e)
{
//std::cout << "dcode_sTogether: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCoefficientRules.get()));

    if (elength(e) == 2 && ehas_head_sym(echild(e,2), gs.sym_sList.get()))
    {
        poly p(elength(echild(e,2)));
        if (!ex_to_polynomial(p, echild(e,1), echild(e,2)))
        {
            return ecopy(e);
        }
//        compiledpoly prog;
//        compile_poly(prog, p);
//        std::cout << "compiled: " << prog.tostring() << std::endl;
/*
        std::vector<wex> stack;
        for (size_t i = 0; i < p.nvars; i++)
        {
            stack.push_back(emake_node(gs.sym_sSlot.copy(), emake_int_ui(i+1)));
        }
        stack.push_back(wex(gs.sym_sNull.copy()));
        stack.push_back(wex(gs.sym_sNull.copy()));
        stack.push_back(wex(gs.sym_sNull.copy()));
        stack.push_back(wex(gs.sym_sNull.copy()));
        eval_poly_ex(stack, prog.prog.data(), prog.prog.size());
std::cout << "output from program: " << ex_tostring(stack[p.nvars].get()) << std::endl;
*/
        uex r; r.init_push_backr(gs.sym_sList.get(), p.size());
        for (size_t i = 0; i < p.size(); i++)
        {
            uex ve; ve.init_push_backr(gs.sym_sList.get(), p.vars.size());
            for (size_t j = 0; j < p.vars.size(); j++)
                ve.push_back(emake_int_copy(p.exps.data()[i*p.vars.size() + j].data));
            r.push_back(emake_node(gs.sym_sRule.copy(), ve.release(), p.coeffs[i].copy()));
        }
        return r.release();
    }
    else
    {
        return ecopy(e);
    }
}
