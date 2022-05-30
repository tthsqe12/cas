#include <cmath>
#include <cfloat>

#include "timing.h"
#include "uex.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"

ex ex_negate(ex E)
{
    uex e(E);
//printf("entering eval_change_sign: "); expr_printout(e);printf("\n");

    if (eis_number(e.get()))
        return num_Minus(e.get());

    if (eis_node(E))
    {
        ulong n = elength(E);
        if (n > 0 && ehas_head_sym(E, gs.sym_sTimes.get()))
        {
            if (eis_number(echild(E,1)))
            {
                ex f = ex_negate(ecopychild(E,1));
                if (eis_one(f))
                {
                    e.removechild(1);
                    if (elength(e.get())==1)
                    {
                        return ecopychild(e.get(), 1);
                    }
                }
                else
                {
                    e.replacechild(1, f);
                }
            }
            else
            {
                e.insertchild(1, emake_cint(-1));
            }
            return e.release();
        }
        else if (ehas_head_sym(E, gs.sym_sPlus.get()))
        {
            uex f; f.init_push_backr(gs.sym_sPlus.get(), n);
            for (ulong i = 0; i < n; i++)
            {
                f.push_back(ex_negate(ecopychild(E,i+1)));
            }
            return f.release();
        }
    }

    return emake_node(gs.sym_sTimes.copy(), emake_cint(-1), e.release());
}


ex ex_reciprocal(ex E)
{
    uex e(E);
//std::cout << "etake_reciprocal: " << ex_tostring_full(e.get()) << std::endl;

    if (eis_number(E))
        return num_Divide(e.get());

    if (eis_node(E))
    {
        ulong n = elength(E);
        if (n == 2 && eis_sym(echild(E,0), gs.sym_sPower.get()))
        {
            ex f = ex_negate(ecopychild(E,2));
            if (eis_one(f))
            {
                return ecopychild(E,1);
            }
            else
            {
                e.replacechild(2, f);
                return e.release();
            }
        }
        else if (eis_sym(echild(E,0), gs.sym_sTimes.get()))
        {
            uex f; f.init_push_backr(gs.sym_sTimes.get(), n);
            for (ulong i = 0; i < n; i++)
                f.push_back(ex_reciprocal(ecopychild(E,i+1)));
            return f.release();
        }
    }

    return emake_node(gs.sym_sPower.copy(), e.release(), emake_cint(-1));
}



ex ex_canonicalize_min(ex EE)
{
    uex e(EE);
//std::cout << " 1 e: " << e.get() << std::endl;
    e.reset(eflatten_sym(e.release(), gs.sym_sMin.get()));
    e.reset(ex_sort(e.release()));
//std::cout << "  sorted e: " << ex_tostring(e.get()) << std::endl;

    std::vector<uex> f;
    bool changed = false;
    er E = e.get();
    ulong en = elength(e.get());
    for (size_t ei=1; ei<=en; ei++)
    {
        if (f.empty())
            goto no_match;

        if (eis_number(f.back().get()) && eis_number(echild(E,ei)))
        {
            uex sum(num_Min(f.back().get(), echild(E,ei)));
            f.pop_back();
            if (!eis_int(sum.get(), 0))
                f.push_back(std::move(sum));
            changed = true;
            continue;
        }
        else if (ex_same(f.back().get(), echild(E,ei)))
        {
            changed = true;
            continue;
        }

    no_match:

        if (!ex_same(echild(E,ei), gs.const_infinity.get()))
        {
            f.push_back(uex(ecopychild(E,ei)));
        }
        else
        {
            changed = true;
        }
    }

    if (changed)
    {
        e.setnz(emake_node(gs.sym_sMax.copy(), f));
    }

    if (elength(e.get())<=1)
    {
        e.setnz(elength(e.get())==1 ? ecopychild(e.get(),1) : gs.const_infinity.copy());
    }

    return e.release();
}

ex ex_canonicalize_max(ex EE)
{
    uex e(EE);
    e.reset(eflatten_sym(e.release(), gs.sym_sMax.get()));
    e.reset(ex_sort(e.release()));

    std::vector<uex> f;
    bool changed = false;
    er E = e.get();
    ulong en = elength(e.get());
    for (size_t ei=1; ei<=en; ei++)
    {
        if (f.empty())
            goto no_match;

        if (eis_number(f.back().get()) && eis_number(echild(E,ei)))
        {
            uex sum(num_Max(f.back().get(), echild(E,ei)));
            f.pop_back();
            if (!eis_int(sum.get(), 0))
                f.push_back(std::move(sum));
            changed = true;
            continue;
        }
        else if (ex_same(f.back().get(), echild(E,ei)))
        {
            changed = true;
            continue;
        }

    no_match:

        if (!ex_same(echild(E,ei), gs.const_minfinity.get()))
        {
            f.push_back(uex(ecopychild(E,ei)));
        }
        else
        {
            changed = true;
        }
    }

    if (changed)
    {
        e.setnz(emake_node(gs.sym_sMax.copy(), f));
    }

    if (elength(e.get())<=1)
    {
        e.setnz(elength(e.get())==1 ? ecopychild(e.get(),1) : gs.const_minfinity.copy());
    }

    return e.release();
}


ex ex_canonicalize_plus(ex EE)
{
    uex e(EE);
//std::cout << " 1 e: " << e.get() << std::endl;
    e.reset(eflatten_sym(e.release(), gs.sym_sPlus.get()));
    e.reset(ex_sort(e.release()));
//std::cout << "  sorted e: " << ex_tostring(e.get()) << std::endl;

    std::vector<uex> f;
    bool changed = false;
    er E = e.get();
    ulong en = elength(e.get());
    for (size_t ei=1; ei<=en; ei++) {
//std::cout << "looking at: " << ex_tostring_full(echild(E,ei)) << std::endl;
//std::cout << "f: " << exvector_tostring_full(f) << std::endl;

        if (f.empty()) {goto no_match;}

    // both nontimes
        if ( !ehas_head_sym(f.back().get(), gs.sym_sTimes.get())
           &&!ehas_head_sym(echild(E,ei), gs.sym_sTimes.get())) {

            if (eis_number(f.back().get()) && eis_number(echild(E,ei)))
            {
                uex sum(num_Plus(f.back().get(), echild(E,ei)));
                f.pop_back();
                if (!eis_int(sum.get(), 0))
                    f.push_back(std::move(sum));
                changed = true;
                continue;
            }
            else if (ex_same(f.back().get(), echild(E,ei)))
            {
                f.back().reset(emake_node(gs.sym_sTimes.copy(), emake_cint(2), ecopychild(E,ei)));
                changed = true;
                continue;
            }

            goto no_match;

    // both times
        }
        else if ( ehas_head_sym(f.back().get(), gs.sym_sTimes.get())
               && ehas_head_sym(echild(E,ei), gs.sym_sTimes.get()))
        {

            er x1 = eget_cint(1), x2 = eget_cint(1);
            size_t n1 = elength(f.back().get()), n2 = elength(echild(E,ei));

            size_t k1=1, k2=1;
            if (k1<=n1 && eis_number(echild(f.back().get(),k1))) {x1=echild(f.back().get(),k1); k1++;}
            if (k2<=n2 && eis_number(echild(E,ei,k2))) {x2=echild(E,ei,k2); k2++;}

            if (n1-k1!=n2-k2) {goto no_match;}
            for (; k1<=n1; k1++,k2++) {
                if (!ex_same(echild(f.back().get(),k1), echild(E,ei,k2))) {goto no_match;};
            }

            uex sum(num_Plus(x1, x2));
            if (eis_int(sum.get(), 0))
            {
                f.pop_back();
            }
            else if (eis_int(sum.get(), 1))
            {
                if (x1 != eget_cint(1))
                {
                    f.back().removechild(1);
                    if (elength(f.back().get())==1)
                    {
                        f.back().reset(ecopychild(f.back().get(),1));
                    }
                }
            } else {
                if (x1 != eget_cint(1)) {
                    f.back().replacechild(1, sum.release());
                } else {
                    f.back().insertchild(1, sum.release());
                }
            }
            changed = true;
            continue;

    // out nontimes, ei times
        }
        else if (   !ehas_head_sym(f.back().get(), gs.sym_sTimes.get())
                    && ehas_head_sym_length(echild(E,ei), gs.sym_sTimes.get(), 2)
                    && eis_number(echild(E,ei,1))
                    && ex_same(f.back().get(), echild(E,ei,2))         )
        {
            uex sum(num_Plus(echild(E,ei,1), eget_cint(1)));
            f.pop_back();
            if (!eis_int(sum.get(), 0))
            {
                f.push_back(uex(emake_node(ecopy(gs.sym_sTimes.get()), sum.release(), ecopychild(E,ei,2))));
            }
            changed = true;
            continue;

    // out times, ei nontimes
        }
        else if (   !ehas_head_sym(echild(E,ei), gs.sym_sTimes.get())
                    && ehas_head_sym_length(f.back().get(), gs.sym_sTimes.get(), 2)
                    && eis_number(echild(f.back().get(),1))
                    && ex_same(echild(E,ei), echild(f.back().get(),2))    )
        {
            uex sum(num_Plus(echild(f.back().get(),1), eget_cint(1)));
            f.pop_back();
            if (!eis_int(sum.get(), 0))
            {
                f.push_back(uex(emake_node(gs.sym_sTimes.copy(), sum.release(), ecopychild(E,ei))));
            }
            changed = true;
            continue;
        }

no_match:

        if (!eis_int(echild(E,ei), 0))
        {
            f.push_back(uex(ecopychild(E,ei)));
        }
        else
        {
            changed = true;
        }
    }

//std::cout << "changed: " << changed << std::endl;

    if (changed)
    {
        e.reset(emake_node(gs.sym_sPlus.copy(), f));
    }

    if (elength(e.get())<=1)
    {
        e.reset(elength(e.get())==1 ? ecopychild(e.get(),1) : emake_cint(0));
    }

//std::cout << "returning dcode_sPlus: " << ex_tostring_full(e.get()) << std::endl;
    return e.release();
}

static int append_pow_normalize(
    uex & m1pow,
    std::vector<std::pair<xfmpz_t, wex>>& l,
    fmpz_t a,
    uex & e,
    slong i)
{
    bool changed = false;

    if (eis_zero(e.get()))
        return changed;

    assert(!fmpz_is_zero(a));

    if (fmpz_sgn(a) < 0)
    {
        fmpz_neg(a, a);
        m1pow.setz(ex_addx(m1pow.release(), e.copy()));
    }

    xfmpz_t g, lbar, abar;

    while (i < l.size() && !fmpz_is_one(a))
    {
        fmpz_gcd(g.data, l[i].first.data, a);
        fmpz_divexact(lbar.data, l[i].first.data, g.data);
        fmpz_divexact(abar.data, a, g.data);

        if (fmpz_is_one(g.data))
        {
            i++;
        }
        else if (fmpz_is_one(lbar.data))
        {
            changed = true;
            fmpz_set(a, abar.data);
            fmpz_set(l[i].first.data, g.data);
            l[i].second.reset(ex_addx(l[i].second.copy(), e.copy()));
            if (eis_zero(l[i].second.get()))
            {
                std::swap(l[i], l.back());
                l.pop_back();
            }
        }
        else if (fmpz_is_one(abar.data))
        {
            changed = true;
            fmpz_set(l[i].first.data, lbar.data);
            fmpz_set(a, g.data);
            e.setz(ex_addx(l[i].second.copy(), e.release()));
            if (eis_zero(e.get()))
                return changed;
        }
        else
        {
            changed = true;
            fmpz_set(a, abar.data);
            append_pow_normalize(m1pow, l, g.data, e, i);
        }
    }

    if (!fmpz_is_one(a))
        l.push_back(std::pair<xfmpz_t, wex>(a, e.copy()));

    return changed;
}


ex ex_canonicalize_rat_powers(ex EE)
{
    bool changed = false;
    uex e(EE);

    if (!ehas_head_sym(e.get(), gs.sym_sTimes.get()))
        return e.release();

    std::vector<std::pair<xfmpz_t, wex>> l;
    std::vector<wex> f;
    xfmpz_t b;
    uex p(emake_cint(1));
    uex m1pow(emake_cint(0));

    for (ulong i = 0; i < elength(e.get()); i++)
    {
        er ei = echild(e.get(),i+1);
        if (eis_int(ei))
        {
            fmpz_set(b.data, eint_data(ei));
            p.setnz(emake_cint(1));
            changed |= append_pow_normalize(m1pow, l, b.data, p, 0);
        }
        else if (eis_rat(ei))
        {
            fmpz_set(b.data, fmpq_numref(erat_data(ei)));
            p.setnz(emake_cint(1));
            changed |= append_pow_normalize(m1pow, l, b.data, p, 0);
            fmpz_set(b.data, fmpq_denref(erat_data(ei)));
            p.setnz(emake_cint(-1));
            changed |= append_pow_normalize(m1pow, l, b.data, p, 0);
        }
        else if (ehas_head_sym_length(ei, gs.sym_sPower.get(), 2) &&
                 eis_int(echild(ei,1)))
        {
            fmpz_set(b.data, eint_data(echild(ei,1)));
            p.setnz(ecopychild(ei,2));
            changed |= append_pow_normalize(m1pow, l, b.data, p, 0);
        }
        else if (ehas_head_sym_length(ei, gs.sym_sPower.get(), 2) &&
                 eis_rat(echild(ei,1)))
        {
            fmpz_set(b.data, fmpq_numref(erat_data(echild(ei,1))));
            p.setnz(ecopychild(ei,2));
            changed |= append_pow_normalize(m1pow, l, b.data, p, 0);
            fmpz_set(b.data, fmpq_denref(erat_data(echild(ei,1))));
            p.setnz(ex_negate(ecopychild(ei,2)));
            changed |= append_pow_normalize(m1pow, l, b.data, p, 0);
        }
        else
        {
            f.push_back(wex(ecopy(ei)));
        }
    }

    if (!changed)
        return e.release();

    for (size_t i = 0; i < l.size(); i++)
        f.push_back(emake_node(gs.sym_sPower.copy(),
                               emake_int_move(l[i].first),
                               l[i].second.copy()));

    f.push_back(emake_node(gs.sym_sPower.copy(),
                           emake_cint(-1),
                           m1pow.release()));

    return emake_node_times(f);
}

ex ex_canonicalize_times(ex EE)
{
    uex e(EE);

    e.reset(eflatten_sym(e.release(), gs.sym_sTimes.get()));
    e.reset(ex_sort(e.release()));

    std::vector<uex> f;
    bool changed = false;
    er E = e.get();
    ulong en = elength(e.get());
    for (size_t ei=1; ei<=en; ei++)
    {
//std::cout << "looking at: " << ex_tostring_full(echild(E,ei)) << std::endl;
//std::cout << "f: " << exvector_tostring_full(f) << std::endl;

        if (f.empty()) {goto no_match;}

    // both nonpower
        if (!ehas_head_sym(f.back().get(), gs.sym_sPower.get()) &&
            !ehas_head_sym(echild(E,ei), gs.sym_sPower.get()))
        {
            if (eis_number(f.back().get()) && eis_number(echild(E,ei)))
            {
                uex t(num_Times(f.back().get(), echild(E,ei)));
                f.pop_back();
                if (!eis_int(t.get(), 1))
                    f.push_back(std::move(t));
                changed = true;
                continue;
            }
            else if (ex_same(f.back().get(), echild(E,ei)))
            {
                f.back().reset(emake_node(gs.sym_sPower.copy(), ecopychild(E,ei), emake_cint(2)));
                changed = true;
                continue;
            }

            goto no_match;

    // both power
        }
        else if (ehas_head_sym_length(f.back().get(), gs.sym_sPower.get(), 2) &&
                 ehas_head_sym_length(echild(E,ei), gs.sym_sPower.get(), 2) &&
                 ex_same(echild(f.back().get(),1), echild(E,ei,1)))
        {
            uex t(ex_addx(ecopychild(f.back().get(),2), ecopychild(E,ei,2)));
            f.pop_back();
            if (eis_zero(t.get()))
            {
            }
            else if (eis_one(t.get()))
            {
                f.push_back(uex(ecopy(echild(E,ei,1))));
            }
            else
            {
                f.push_back(uex(emake_node(gs.sym_sPower.copy(), ecopychild(E,ei,1), t.release())));
            }
            changed = true;
            continue;

    // out nonpower, ei power
        }
        else if (!eis_number(f.back().get()) &&
                 !ehas_head_sym(f.back().get(), gs.sym_sPower.get()) &&
                 ehas_head_sym_length(echild(E,ei), gs.sym_sPower.get(), 2) &&
                 ex_same(f.back().get(), echild(E,ei,1)))
        {
            uex t(ex_addx(ecopychild(E,ei,2), emake_cint(1)));
            f.pop_back();
            if (!eis_int(t.get(), 0))
                f.push_back(uex(emake_node(ecopy(gs.sym_sPower.get()), ecopychild(E,ei,1), t.release())));
            changed = true;
            continue;
        }
    // out times, ei nontimes
        else if (!ehas_head_sym(echild(E,ei), gs.sym_sPower.get()) &&
                 ehas_head_sym_length(f.back().get(), gs.sym_sPower.get(), 2) &&
                 ex_same(echild(E,ei), echild(f.back().get(),1)))
        {
            uex t(ex_addx(ecopychild(f.back().get(),2), emake_cint(1)));
            f.pop_back();
            if (!eis_int(t.get(), 0))
                f.push_back(uex(emake_node(gs.sym_sPower.copy(), ecopychild(E,ei), t.release())));
            changed = true;
            continue;
        }
    // a*I*(-1)^r -> a*(-1)^(r+1/2)
        else if (eis_cmplx(f.back().get()) &&
                 eis_zero(ecmplx_real(f.back().get())) &&
                 ehas_head_sym_length(echild(E,ei), gs.sym_sPower.get(), 2) &&
                 eis_int(echild(E,ei,1), -1) &&
                 eis_rat(echild(E,ei,2)))
        {
            f.back().setnz(ecopy(ecmplx_imag(f.back().get())));
            ex t = ex_addx(ecopychild(E,ei,2), emake_crat(1, 2));
            f.emplace_back(ex_powx(ecopychild(E,ei,1), t));
            changed = true;
            continue;
        }

no_match:

        if (!eis_int(echild(E,ei), 1))
            f.emplace_back(ecopychild(E,ei));
        else
            changed = true;
    }

    if (changed)
        e.reset(emake_node(gs.sym_sTimes.copy(), f));

    if (elength(e.get())<=1)
    {
        e.reset(elength(e.get())==1 ? ecopychild(e.get(),1) : emake_cint(1));
    }
    else if (eis_int(echild(e.get(),1), 0))
    {
        e.reset(emake_cint(0));
    }

//std::cout << "returning dcode_sTimes: " << ex_tostring_full(e.get()) << std::endl;
    return ex_canonicalize_rat_powers(e.release());
}

static ex _denest_power(er A, er B)
{
    ex p = ex_mulx(ecopy(B), ecopychild(A,2));
    if (eis_int(p, 1))
    {
        eclear(p);
        return ecopychild(A,1);
    }
    else
    {
        return ex_powx(ecopychild(A,1), p);
    }
}


ex ex_canonicalize_power(ex EE)
{
    uex e(EE);
    assert(ehas_head_sym(e.get(), gs.sym_sPower.get()));
    assert(elength(e.get()) == 2);

    er A = echild(e.get(),1);
    er B = echild(e.get(),2);

    if (eis_number(B))
    {
        if (eis_number(A))
            return num_Power(A, B, e.get());

        if (eis_int(B, 1))
            return ecopy(A);

        if (eis_int(B, 0))
            return emake_cint(1);

        if (eis_int(B))
        {
            if (ehas_head_sym_length(A, gs.sym_sPower.get(), 2))
            {
                return _denest_power(A, B);
            }
            else if (ehas_head_sym(A, gs.sym_sTimes.get()) && elength(A) > 0)
            {
                ulong n = elength(A);
                uex a; a.init_push_backr(gs.sym_sTimes.get(), n);
                for (ulong i = 0; i < n; i++)
                    a.push_back(ex_powx(ecopychild(A,i+1), ecopy(B)));
                return ex_canonicalize_times(a.release());
            }
            else
            {
                return e.release();
            }
        }
    }

    if (ehas_head_sym_length(A, gs.sym_sPower.get(), 2) && eis_rat(echild(A, 2)))
    {
        fmpq * x = erat_data(echild(A, 2));
        if (0 <= fmpz_sgn(fmpq_numref(x)) &&
            fmpz_cmp(fmpq_numref(x), fmpq_denref(x)) <= 0)
        {
            return _denest_power(A, B);
        }
    }

    return e.release();
}


ex ex_power_expand(er e)
{
    if (!eis_node(e))
        return ecopy(e);

    ulong n = elength(e);
    bool changed = false;

    uex f; f.init_push_backx(ex_power_expand(echild(e,0)), n);
    for (ulong i = 0; i < n; i++)
    {
        er ei = echild(e,i+1);
        ex fi = ex_power_expand(ei);
        if (ei != etor(fi))
            changed = true;
        f.push_back(fi);
    }
    if (!changed)
        f.setnz(ecopy(e));
    e = f.get();

    if (ehas_head_sym(e, gs.sym_sPlus.get()))
    {
        return ex_canonicalize_plus(f.release());
    }
    else if (ehas_head_sym(e, gs.sym_sTimes.get()))
    {
        return ex_canonicalize_times(f.release());
    }
    else if (ehas_head_sym_length(e, gs.sym_sPower.get(), 2))
    {
        er b = echild(e, 1);
        if (!eis_node(b))
            return f.release();

        if (ehas_head_sym_length(b, gs.sym_sPower.get(), 2))
        {
            ex t = ex_mulr(echild(b,2), echild(e,2));
            return ex_powx(ecopychild(b,1), t);
        }
        else if (ehas_head_sym(b, gs.sym_sTimes.get()))
        {
            n = elength(b);
            uex g; g.init_push_backr(echild(b,0), n);
            for (ulong i = 0; i < n; i++)
            {
                er c = echild(b,i+1);
                if (ehas_head_sym_length(c, gs.sym_sPower.get(), 2))
                {
                    ex t = ex_mulr(echild(c,2), echild(e,2));
                    g.push_back(ex_powx(ecopychild(c,1), t));
                }
                else
                {
                    g.push_back(ex_powr(c, echild(e,2)));
                }
            }
            return ex_canonicalize_times(g.release());
        }
        else
        {
            return f.release();
        }
    }
    else if (ehas_head_sym_length(e, gs.sym_sLog.get(), 1))
    {
        er b = echild(e, 1);
        if (!eis_node(b))
            return f.release();

        if (ehas_head_sym_length(b, gs.sym_sPower.get(), 2))
        {
            ex t = emake_node(gs.sym_sLog.copy(), ecopychild(b,1));
            return ex_mulx(ecopychild(b,2), t);
        }
        else if (ehas_head_sym(b, gs.sym_sTimes.get()))
        {
            n = elength(b);
            uex g; g.init_push_backr(gs.sym_sPlus.get(), n);
            for (ulong i = 0; i < n; i++)
            {
                er c = echild(b,i+1);
                if (ehas_head_sym_length(c, gs.sym_sPower.get(), 2))
                {
                    ex t = emake_node(gs.sym_sLog.copy(), ecopychild(c,1));
                    g.push_back(ex_mulx(ecopychild(c,2), t));
                }
                else
                {
                    g.push_back(emake_node(gs.sym_sLog.copy(), ecopy(c)));
                }
            }
            return ex_canonicalize_plus(g.release());
        }
        else
        {
            return f.release();
        }
    }
    else
    {
        return f.release();
    }
}

ex dcode_sPowerExpand(er e)
{
    if (elength(e) != 1)
        return ecopy(e);

    return ex_power_expand(echild(e,1));
}


ex dcode_sMin(er e)
{
//std::cout << "dcode_sMin: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMin.get()));

    ulong n = elength(e);

    if (n < 2)
    {
        return (n == 0) ? gs.const_infinity.copy() : ecopychild(e,1);
    }
    else if (n == 2 && eis_number(echild(e,1)) && eis_number(echild(e,2)))
    {
        return num_Min(echild(e,1), echild(e,2));
    }

    return ex_canonicalize_min(ecopy(e));
}

ex dcode_sMax(er e)
{
//std::cout << "dcode_sMax: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMax.get()));

    ulong n = elength(e);

    if (n < 2)
    {
        return (n == 0) ? gs.const_minfinity.copy() : ecopychild(e,1);
    }
    else if (n == 2 && eis_number(echild(e,1)) && eis_number(echild(e,2)))
    {
        return num_Max(echild(e,1), echild(e,2));
    }

    return ex_canonicalize_max(ecopy(e));
}

ex dcode_sPlus(er e)
{
//std::cout << "dcode_sPlus: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sPlus.get()));

    ulong n = elength(e);

    if (n <= 1)
        return (n == 0) ? emake_cint(0) : ecopychild(e,1);

    if (eis_number(echild(e,1)))
    {
        if (n == 2 && eis_number(echild(e,2)))
            return num_Plus(echild(e,1), echild(e,2));
    }

    return ex_canonicalize_plus(ecopy(e));
}


ex dcode_sTimes(er e)
{
//std::cout << "dcode_sTimes: " << ex_tostring_full(EE) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sTimes.get()));

    ulong n = elength(e);

    if (n < 2)
        return n == 0 ? emake_cint(1) : ecopychild(e,1);    // Times[] -> 1, Times[a] -> a

    if (eis_number(echild(e,1)) && eis_number(echild(e,2)))
    {
        if (n > 2)
        {
            if (eis_number(echild(e,3)))
            {
                if (n >= 4)
                {
                    if (n == 4 && eis_number(echild(e,4)))
                    {
                        uex t1(num_Times(echild(e, 1), echild(e, 2)));
                        uex t2(num_Times(echild(e, 3), echild(e, 4)));
                        return num_Times(t1.get(), t2.get());
                    }
                }
                else
                {
                    uex t1(num_Times(echild(e,1), echild(e,2)));
                    return num_Times(t1.get(), echild(e, 3));
                }
            }
        }
        else
        {
            return num_Times(echild(e,1), echild(e,2));
        }
    }

    return ex_canonicalize_times(ecopy(e));
}



ex dcode_sMinus(er e)
{
//std::cout << "dcode_sMinus: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMinus.get()));

    ulong n = elength(e);
    if (n <= 1)
        return (n == 0) ? ecopy(e) : ex_negate(ecopychild(e,1));

    if (n == 2 && eis_number(echild(e,1)) && eis_number(echild(e,2)))
        return num_Minus(echild(e,1), echild(e,2));

    uex f; f.init_push_backr(gs.sym_sPlus.get(), n);
    f.push_back(ecopychild(e,1));
    for (ulong i=1; i < n; i++)
        f.push_back(ex_negate(ecopychild(e,i+1)));
    return ex_canonicalize_plus(f.release());
}



ex dcode_sDivide(er e)
{
//std::cout << "dcode_sDivide: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDivide.get()));

    ulong n = elength(e);
    if (n <= 1)
        return (n == 0) ? ecopy(e) : ex_reciprocal(ecopychild(e,1));

    if (n == 2 && eis_number(echild(e,1)) && eis_number(echild(e,2)))
        return num_Divide(echild(e,1), echild(e,2));

    uex f; f.init_push_backr(gs.sym_sTimes.get(), n);
    f.push_back(ecopychild(e,1));
    for (ulong i=1; i < n; i++)
        f.push_back(ex_reciprocal(ecopychild(e,i+1)));
    return ex_canonicalize_times(f.release());
}


ex dcode_sPower(er e)
{
//std::cout << "dcode_sPower: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sPower.get()));

    if (elength(e) == 2)
    {
        return ex_canonicalize_power(ecopy(e));
    }
    else
    {
        ex r = emake_cint(1);
        for (size_t i = elength(e); i != 0; i--)
            r = ex_powx(ecopychild(e,i), r);
        return r;
    }
}


ex dcode_sSqrt(er e)
{
//std::cout << "dcode_sSqrt: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSqrt.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);
    if (eis_number(x))
        return num_Sqrt(x);

    if (ehas_head_sym(x, gs.sym_sTimes.get()))
    {
        if (elength(x) == 0 || !eis_number(echild(x,1)))
            return emake_node(gs.sym_sPower.copy(), ecopy(x), emake_crat(1,2));

        uex y(num_Sqrt(echild(x,1)));

        if (eis_int(y.get()) || eis_rat(y.get()) || eis_double(y.get()) || eis_real(y.get()))
        {
            uex z(ecopy(x));
            z.removechild(1);
            z.setz(emake_node(gs.sym_sPower.copy(), z.release(), emake_crat(1,2)));
            return emake_node(gs.sym_sTimes.copy(), y.release(), z.release());
        }
    }

    return emake_node(gs.sym_sPower.copy(), ecopy(x), emake_crat(1,2));
}
