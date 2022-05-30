#include "globalstate.h"
#include "code.h"
#include "ex_cont.h"

ex dcode_sPart(er e)
{
//std::cout << "dcode_sPart: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sPart.get()));

    if (elength(e) == 0)
        return _handle_message_argm(e, 1);

    uex r(ecopychild(e,1));
    for (ulong n = 2; n <= elength(e); n++)
    {
        er idx = echild(e, n);
        if (eis_node(r.get()))
        {
            if (!eis_intm(idx))
            {
                _gen_message(gs.sym_sPart.get(), "pkspec1", NULL, ecopy(idx));
                return ecopy(e);
            }
            ulong i = eintm_get(idx); // TODO
            if (i > elength(r.get()))
            {
                _gen_message(gs.sym_sPart.get(), "partw", NULL, ecopy(idx), r.copy());
                return ecopy(e);
            }
            r.setnz(r.copychild(i));
        }
        else if (eis_leaf(r.get()))
        {
            _gen_message(gs.sym_sPart.get(), "partd", NULL, ecopy(e));
            return ecopy(e);
        }
        else if (eis_hmap(r.get()))
        {
            if (ehas_head_sym_length(idx, gs.sym_sKey.get(), 1))
            {
                wex key(ecopychild(idx,1));
                auto it = eto_hmap(r.get()).map_data.find(key);
                if (it == eto_hmap(r.get()).map_data.end())
                {
                    ex t = emake_str("KeyAbsent");
                    return emake_node(gs.sym_sMissing.copy(), t, key.copy());
                }
                else
                {
                    r.setnz(it->second.second.copy());
                }
            }
            else
            {
                if (!eis_intm(idx))
                {
                    _gen_message(gs.sym_sPart.get(), "pkspec1", NULL, ecopy(idx));
                    return ecopy(e);
                }
                ulong i = eintm_get(idx) - 1;
                if (i >= eto_hmap(r.get()).map_data.size())
                {
                    _gen_message(gs.sym_sPart.get(), "partw", NULL, ecopy(idx), r.copy());
                    return ecopy(e);
                }
                for (std::pair<wex, std::pair<ulong, wex>> it : eto_hmap(r.get()).map_data)
                {
                    if (i == it.second.first)
                    {
                        r.setnz(it.second.second.copy());
                        break;
                    }
                }
            }
        }
        else
        {
            _gen_message(gs.sym_sPart.get(), "partd", NULL, ecopy(e));
            return ecopy(e);                
        }
    } 

    return r.release();
}

static ex _head(er x)
{
    switch (etype(x))
    {
        case ETYPE_INT:
        {
            return gs.sym_sInteger.copy();
        }
        case ETYPE_RAT:
        {
            return gs.sym_sRational.copy();
        }
        case ETYPE_DOUBLE:
        case ETYPE_REAL:
        {
            return gs.sym_sReal.copy();
        }
        case ETYPE_CMPLX:
        {
            return gs.sym_sComplex.copy();
        }
        case ETYPE_NAN:
        {
            return gs.sym_sIndeterminate.copy();
        }
        case ETYPE_SYM:
        {
            return gs.sym_sSymbol.copy();
        }
        case ETYPE_STR:
        {
            return gs.sym_sString.copy();
        }
        case ETYPE_NODE:
        {
            return ecopychild(x,0);
        }
        default:
        {
            assert(false);
            return nullptr;
        }
    }
}

ex dcode_sHead(er e)
{
//std::cout << "dcode_sHead: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sHead.get()));

    if (unlikely(elength(e) < 1 || elength(e) > 2))
        return _handle_message_argt(e, (1 << 0) + (2 << 8));

    ex h = _head(echild(e,1));
    return elength(e) == 1 ? h : emake_node(ecopychild(e,2), h);
}

ex dcode_sFirst(er e)
{
//std::cout << "dcode_sFirst: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFirst.get()));

    if (unlikely(elength(e) < 1 || elength(e) > 2))
        return _handle_message_argt(e, (1 << 0) + (2 << 8));

    er f = elength(e) == 1 ? nullptr : echild(e,2);
    er x = echild(e,1);

    if (eis_node(x))
    {
        if (elength(x) > 0)
            return ecopychild(x,1);

        if (f != nullptr)
            return ecopy(f);

        _gen_message(echild(e,0), "nofirst", "`1` has a length of zero and no first element.", ecopy(e));
        return ecopy(e);
    }
    else if (eis_parray(x))
    {
        return eparray_child(x, 0);
    }
    else
    {
        if (f != nullptr)
            return ecopy(f);

        _gen_message(echild(e,0), "normal", NULL, emake_cint(1), ecopy(e));
        return ecopy(e);
    }
}

ex dcode_sLast(er e)
{
//std::cout << "dcode_sLast: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sLast.get()));

    if (elength(e) != 1 && elength(e) != 2)
        return _handle_message_argt(e, (1 << 0) + (2 << 8));

    er f = elength(e) == 1 ? nullptr : echild(e,2);
    er x = echild(e,1);

    if (eis_node(x))
    {
        if (elength(x) > 0)
            return ecopychild(x,elength(x));

        if (f != nullptr)
            return ecopy(f);

        _gen_message(echild(e,0), "nolast", "`1` has a length of zero and no last element.", ecopy(e));
        return ecopy(e);
    }
    else if (eis_parray(x))
    {
        return eparray_child(x, eparray_dimension(x,0) - 1);
    }
    else
    {
        if (f != nullptr)
            return ecopy(f);

        _gen_message(echild(e,0), "normal", NULL, emake_cint(1), ecopy(e));
        return ecopy(e);
    }
}

ex dcode_sRest(er e)
{
//std::cout << "dcode_sRest: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sRest.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    if (eis_node(x))
    {
        if (elength(x) < 1)
        {
            _gen_message(echild(e,0), "norest", "Cannot take Rest of expression `1` with length zero.", ecopy(e));
            return ecopy(e);
        }

        uex r(ecopy(x));
        eremovepart(r, 1);
        return r.release();
    }
    else if (eis_parray(x))
    {
        return eparray_rest(x);
    }
    else
    {
        _gen_message(echild(e,0), "normal", NULL, emake_cint(1), ecopy(e));
        return ecopy(e);
    }
}

ex dcode_sMost(er e)
{
//std::cout << "dcode_sMost: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMost.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    er x = echild(e,1);

    if (eis_node(x))
    {
        if (elength(x) < 1)
        {
            _gen_message(echild(e,0), "nomost", "Cannot take Most of expression `1` with length zero.", ecopy(e));
            return ecopy(e);
        }

        uex r(ecopy(x));
        eremovepart(r, elength(x));
        return r.release();
    }
    else if (eis_parray(x))
    {
        return eparray_most(x);
    }
    else
    {
        _gen_message(echild(e,0), "normal", NULL, emake_cint(1), ecopy(e));
        return ecopy(e);
    }
}
