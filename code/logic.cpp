#include "globalstate.h"
#include "code.h"
#include "eval.h"

ex dcode_sNot(er e)
{
//std::cout << "dcode_sNot: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sNot.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    if (eis_sym(echild(e,1), gs.sym_sTrue.get()))
        return gs.sym_sFalse.copy();
    else if (eis_sym(echild(e,1), gs.sym_sFalse.get()))
        return gs.sym_sTrue.copy();
    else
        return ecopy(e);
}

ex dcode_sAnd(er e)
{
//std::cout << "dcode_sAnd: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sAnd.get()));
    std::vector<wex> u;

    ulong n = elength(e);
    for (ulong i = 0; i < n; i++)
    {
        ex c = eval(ecopychild(e, i + 1));
        if (eis_sym(c, gs.sym_sTrue.get()))
            eclear(c);
        else if (eis_sym(c, gs.sym_sFalse.get()))
            return c;
        else
            u.emplace_back(c);
    }

    return u.empty() ? gs.sym_sTrue.copy() : emake_node(gs.sym_sAnd.copy(), u);
}

ex dcode_sOr(er e)
{
//std::cout << "dcode_sOr: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sOr.get()));
    std::vector<wex> u;

    ulong n = elength(e);
    for (ulong i = 0; i < n; i++)
    {
        ex c = eval(ecopychild(e, i + 1));
        if (eis_sym(c, gs.sym_sFalse.get()))
            eclear(c);
        else if (eis_sym(c, gs.sym_sTrue.get()))
            return c;
        else
            u.emplace_back(c);
    }

    return u.empty() ? gs.sym_sFalse.copy() : emake_node(gs.sym_sOr.copy(), u);
}

ex dcode_sNand(er e)
{
//std::cout << "dcode_sNand: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sNand.get()));
    std::vector<wex> u;

    ulong n = elength(e);
    for (ulong i = 0; i < n; i++)
    {
        ex c = eval(ecopychild(e, i + 1));
        if (eis_sym(c, gs.sym_sTrue.get()))
        {
            eclear(c);
        }
        else if (eis_sym(c, gs.sym_sFalse.get()))
        {
            eclear(c);
            return gs.sym_sTrue.copy();
        }
        else
        {
            u.emplace_back(c);
        }
    }

    return u.empty() ? gs.sym_sFalse.copy() : emake_node(gs.sym_sNand.copy(), u);
}

ex dcode_sNor(er e)
{
//std::cout << "dcode_sNor: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sNor.get()));
    std::vector<wex> u;

    ulong n = elength(e);
    for (ulong i = 0; i < n; i++)
    {
        ex c = eval(ecopychild(e, i + 1));
        if (eis_sym(c, gs.sym_sFalse.get()))
        {
            eclear(c);
        }
        else if (eis_sym(c, gs.sym_sTrue.get()))
        {
            eclear(c);
            return gs.sym_sFalse.copy();
        }
        else
        {
            u.emplace_back(c);
        }
    }

    return u.empty() ? gs.sym_sTrue.copy() : emake_node(gs.sym_sNor.copy(), u);
}

ex dcode_sXor(er e)
{
//std::cout << "dcode_sXor: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sXor.get()));
    std::vector<wex> u;
    int count = 0;

    ulong n = elength(e);
    for (ulong i = 0; i < n; i++)
    {
        er c = echild(e, i + 1);
        if (eis_sym(c, gs.sym_sTrue.get()))
            count++;
        else if (!eis_sym(c, gs.sym_sFalse.get()))
            u.emplace_back(ecopy(c));
    }

    if (u.empty())
    {
        return emake_boole(count & 1);
    }
    else
    {
        ex r = emake_node(gs.sym_sXor.copy(), u);
        return (count & 1) ? emake_node(gs.sym_sNot.copy(), r) : r;
    }
}
