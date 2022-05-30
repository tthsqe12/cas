#include "uex.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "ex_cont.h"


static ex _dcode_map_gen(
    er f,
    er e,
    slong & edepth,
    slong elevel,
    const fmpz_t a,
    const fmpz_t b)
{
    ex ne;

    if (fmpz_sgn(b) >= 0 && fmpz_cmp_si(b, elevel) < 0)
        return ecopy(e);

    if (eis_parray(e))
        e = eparray_get_normal(e);

    if (eis_node(e))
    {
        bool changed = false;
        std::vector<wex> v;
        slong child_depth;

        edepth = 2;
        for (ulong i = 0; i < elength(e); i++)
        {
            child_depth = 0;
            ex ei = _dcode_map_gen(f, echild(e,i+1), child_depth, elevel + 1, a, b);
            changed |= (etor(ei) != echild(e,i+1));
            v.emplace_back(ei);
            edepth = std::max(edepth, 1 + child_depth);
        }
        ne = changed ? emake_node(ecopychild(e,0), v) : ecopy(e);
    }
    else
    {
        edepth = 1;
        ne = ecopy(e);
    }

    if ((fmpz_sgn(a) >= 0 ? fmpz_cmp_si(a, elevel) <= 0 : fmpz_cmp_si(a, -edepth) <= 0) &&
        (fmpz_sgn(b) >= 0 ? fmpz_cmp_si(b, elevel) >= 0 : fmpz_cmp_si(b, -edepth) >= 0))
    {
        return emake_node(ecopy(f), ne);
    }
    else
    {
        return ne;
    }
}

static ex _dcode_apply_gen(
    er f,
    er e,
    slong & edepth,
    slong elevel,
    const fmpz_t a,
    const fmpz_t b)
{
    ex ne;

    if (fmpz_sgn(b) >= 0 && fmpz_cmp_si(b, elevel) < 0)
        return ecopy(e);

    if (eis_parray(e))
        e = eparray_get_normal(e);

    if (eis_node(e))
    {
        bool changed = false;
        std::vector<wex> v;
        slong child_depth;

        edepth = 2;
        for (ulong i = 0; i < elength(e); i++)
        {
            child_depth = 0;
            ex ei = _dcode_apply_gen(f, echild(e,i+1), child_depth, elevel + 1, a, b);
            changed |= (etor(ei) != echild(e,i+1));
            v.emplace_back(ei);
            edepth = std::max(edepth, 1 + child_depth);
        }
        ne = changed ? emake_node(ecopychild(e,0), v) : ecopy(e);
    }
    else
    {
        edepth = 1;
        return ecopy(e);
    }

    if ((fmpz_sgn(a) >= 0 ? fmpz_cmp_si(a, elevel) <= 0 : fmpz_cmp_si(a, -edepth) <= 0) &&
        (fmpz_sgn(b) >= 0 ? fmpz_cmp_si(b, elevel) >= 0 : fmpz_cmp_si(b, -edepth) >= 0))
    {
        uex NE(ne);
        NE.replacechild(0, ecopy(f));
        return NE.release();
    }
    else
    {
        return ne;
    }
}

bool parse_level(const fmpz * & n1, const fmpz * & n2, er l)
{
    if (ehas_head_sym(l, gs.sym_sList.get()))
    {
        if (elength(l) == 1 && eis_int(echild(l,1)))
        {
            n1 = n2 = eint_data(echild(l,1));
            return true;
        }
        else if (elength(l) == 2 && eis_int(echild(l,1)) && eis_int(echild(l,2)))
        {
            n1 = eint_data(echild(l,1));
            n2 = eint_data(echild(l,2));
            return true;
        }
        else if (elength(l) == 2 && eis_int(echild(l,1)) && eis_infinity(echild(l,2)))
        {
            n1 = eint_data(echild(l,1));
            n2 = eint_data(gs.const_max_uword.get());
            return true;
        }
        else
        {
            return false;
        }
    }
    else if (eis_int(l))
    {
        n1 = eget_cint_data(1);
        n2 = eint_data(l);
        return true;
    }
    else if (eis_infinity(l))
    {
        n1 = eget_cint_data(1);
        n2 = eint_data(gs.const_max_uword.get());
        return true;
    }
    else
    {
        return false;
    }
}


ex dcode_sMap(er e)
{
//std::cout << "dcode_sMap: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMap.get()));

    slong edepth;
    const fmpz * n1 = eget_cint_data(1);
    const fmpz * n2 = eget_cint_data(1);

    if (elength(e) == 2)
    {
        return _dcode_map_gen(echild(e,1), echild(e,2), edepth, 0, n1, n2);
    }
    else if (elength(e) == 3)
    {
        if (!parse_level(n1, n2, echild(e,3)))
        {
            _gen_message(echild(e,0), "level", NULL, ecopychild(e,3));            
            return ecopy(e);
        }

        return _dcode_map_gen(echild(e,1), echild(e,2), edepth, 0, n1, n2);
    }
    else
    {
        return ecopy(e);
    }
}


ex dcode_sApply(er e)
{
//std::cout << "dcode_sApply: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sApply.get()));

    slong edepth;
    const fmpz * n1 = eget_cint_data(0);
    const fmpz * n2 = eget_cint_data(0);

    if (elength(e) == 2)
    {
        return _dcode_apply_gen(echild(e,1), echild(e,2), edepth, 0, n1, n2);
    }
    else if (elength(e) == 3)
    {
        if (!parse_level(n1, n2, echild(e,3)))
        {
            _gen_message(echild(e,0), "level", NULL, ecopychild(e,3));            
            return ecopy(e);
        }

        return _dcode_apply_gen(echild(e,1), echild(e,2), edepth, 0, n1, n2);
    }
    else
    {
        return ecopy(e);
    }
}



ex dcode_sIdentity(er e)
{
//std::cout << "dcode_sIdentity: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sIdentity.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return ecopychild(e,1);
}


ex dcode_sThrough(er e)
{
//std::cout << "dcode_sThrough: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sThrough.get()));

    if (elength(e) != 1 && elength(e) != 2)
        return _handle_message_argt(e, (1 << 0) + (2 << 8));

    er f = echild(e,1);

    if (!eis_node(f))
        return ecopy(e);

    er h = echild(f,0);

    if (eis_parray(h))
        h = eparray_get_normal(h);

    if (!eis_node(h))
        return ecopy(e);

    if (elength(e) == 2 && !ex_same(echild(h,0), echild(e,2)))
        return ecopy(f);

    uex r; r.init_push_backr(echild(h,0), elength(h));

    for (ulong i = 0; i < elength(h); i++)
        r.push_back(ereplacechild(f,0, ecopychild(h, i+1)));

    return r.release();
}

ex dcode_sComposition(er e)
{
//std::cout << "dcode_sComposition: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sComposition.get()));

    if (elength(e) == 0)
        return gs.sym_sIdentity.copy();

    return ecopy(e);
}

ex scode_sComposition(er e)
{
//std::cout << "scode_sComposition: " << ex_tostring_full(e) << std::endl;
    assert(eis_node(e));

    er h = echild(e,0);

    if (!ehas_head_sym(h, gs.sym_sComposition.get()))
        return ecopy(e);

    if (elength(h) == 0)
        return ecopy(e);

    ulong i = elength(h);
    uex f(ereplacechild(e,0, ecopychild(h,i)));

    for (i--; i > 0; i--)
        f.setz(emake_node(ecopychild(h,i), f.release()));

    return f.release();
}

ex dcode_sRightComposition(er e)
{
//std::cout << "dcode_sRightComposition: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sRightComposition.get()));

    if (elength(e) == 0)
        return gs.sym_sIdentity.copy();

    return ecopy(e);
}

ex scode_sRightComposition(er e)
{
//std::cout << "scode_sRightComposition: " << ex_tostring_full(e) << std::endl;
    assert(eis_node(e));

    er h = echild(e,0);

    if (!ehas_head_sym(h, gs.sym_sRightComposition.get()))
        return ecopy(e);

    if (elength(h) == 0)
        return ecopy(e);

    ulong i = 0;
    uex f(ereplacechild(e,0, ecopychild(h,i+1)));

    for (i++; i < elength(h); i++)
        f.setz(emake_node(ecopychild(h,i+1), f.release()));

    return f.release();
}
