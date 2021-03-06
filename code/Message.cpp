#include "globalstate.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "serialize.h"

ex dcode_sMessage(er e);

er esym_find_message_tag(er x, er tag)
{
    assert(eis_sym(x));
    assert(eis_str(tag));

//std::cout << "searching " << ex_tostring_full(x) << " for " << ex_tostring_full(tag) << std::endl;

    er l = etor(esym_messages(x));
    if (l == nullptr)
        return nullptr;

    assert(eis_node(l));
    for (size_t i = 1; i <= elength(l); i++)
    {
        if (ehas_head_sym_length(echild(l,i), gs.sym_sRule.get(), 2)
            && ex_same(echild(l,i,1), tag))
        {
            return echild(l,i,2);
        }
    }
    return nullptr;
}

/*
    status 0: on
           1: off
           2: don't care - keep 
*/
void esym_set_message_tag(er x, er tag, er s, int status)
{
    assert(eis_sym(x));
    assert(eis_str(tag));
    assert(eis_str(s));

    er l = etor(esym_messages(x));
    if (l == nullptr)
    {
        ex t = ecopy(s);
        if (status == 1)
        {
            t = emake_node(gs.sym_sOff.copy(), t);
        }
        t = emake_node(ecopy(gs.sym_sRule.get()), ecopy(tag), t);
        t = emake_node(ecopy(gs.sym_sList.get()), t);
        esym_set_messages(x, t);
        return;
    }

    for (size_t i = 1; i <= elength(l); i++)
    {
        if (ehas_head_sym_length(echild(l,i), gs.sym_sRule.get(), 2)
            && ex_same(echild(l,i,1), tag))
        {
            er li2 = echild(l,i,2);
            if (eis_str(li2))
            {
                if (status != 1 && ex_same(li2, s))
                {
                    return;
                }

                uex k(ecopy(l));
                uex k1(ecopychild(k.get(),i));
                ex t = ecopy(s);
                if (status == 1)
                {
                    t = emake_node(gs.sym_sOff.copy(), t);
                }
                k1.replacechild(2, t);
                k.replacechild(i, k1.release());
                esym_set_messages(x, k.release());
                return;
            }
            else
            {
                assert(eis_node(li2));

                if (status != 0 && ex_same(echild(li2,1), s))
                    return;

                uex k(ecopy(l));
                uex k1(ecopychild(k.get(),i));
                uex k2(ecopychild(k1.get(),2));
                ex t = ecopy(s);
                if (status == 1)
                {
                    t = emake_node(gs.sym_sOff.copy(), t);
                }
                k2.replacechild(1, t);
                k1.replacechild(2, k2.release());
                k.replacechild(i, k1.release());
                esym_set_messages(x, k.release());
                return;
            }
        }
    }

    uex k(ecopy(l));
    ex t = ecopy(s);
    if (status == 1)
    {
        t = emake_node(gs.sym_sOff.copy(), t);
    }
    k.insertchild(elength(l) + 1, emake_node(gs.sym_sRule.copy(), ecopy(tag), t));
    esym_set_messages(x, k.release());
}

void esym_set_message_tag(er x, er tag, const char * str, int status)
{
    assert(eis_sym(x));
    assert(eis_str(tag));
    uex s(emake_str(str));
    esym_set_message_tag(x, tag, s.get(), status);
}

void esym_set_message_tag(er x, const char * tag, const char * str, int status)
{
    assert(eis_sym(x));
    uex t(emake_str(tag));
    uex s(emake_str(str));
    esym_set_message_tag(x, t.get(), s.get(), status);
}

void _gen_message(er x, const char * tag, const char * str)
{
    assert(eis_sym(x));
    wex t(emake_str(tag));
    if (str != nullptr)
    {
        wex s(emake_str(str));
        esym_set_message_tag(x, t.get(), s.get(), 2);
    }
    ex mn = emake_node(gs.sym_sMessageName.copy(), ecopy(x), t.copy());
    t.reset(emake_node(gs.sym_sMessage.copy(), mn));
    eclear(dcode_sMessage(t.get()));
}

void _gen_message(er x, const char * tag, const char * str, ex E1)
{
    uex e1(E1);
    assert(eis_sym(x));
    wex t(emake_str(tag));
    if (str != nullptr)
    {
        wex s(emake_str(str));
        esym_set_message_tag(x, t.get(), s.get(), 2);
    }
    ex mn = emake_node(gs.sym_sMessageName.copy(), ecopy(x), t.copy());
    t.reset(emake_node(gs.sym_sMessage.copy(), mn, e1.release()));
    eclear(dcode_sMessage(t.get()));
}

void _gen_message(er x, const char * tag, const char * str, ex E1, ex E2)
{
    uex e1(E1);
    uex e2(E2);
    assert(eis_sym(x));
    wex t(emake_str(tag));
    if (str != nullptr)
    {
        wex s(emake_str(str));
        esym_set_message_tag(x, t.get(), s.get(), 2);
    }
    ex mn = emake_node(gs.sym_sMessageName.copy(), ecopy(x), t.copy());
    t.reset(emake_node(gs.sym_sMessage.copy(), mn, e1.release(), e2.release()));
    eclear(dcode_sMessage(t.get()));
}

void _gen_message(er x, const char * tag, const char * str, ex E1, ex E2, ex E3)
{
    uex e1(E1);
    uex e2(E2);
    uex e3(E3);
    assert(eis_sym(x));
    wex t(emake_str(tag));
    if (str != nullptr)
    {
        wex s(emake_str(str));
        esym_set_message_tag(x, t.get(), s.get(), 2);
    }
    ex mn = emake_node(gs.sym_sMessageName.copy(), ecopy(x), t.copy());
    t.reset(emake_node(gs.sym_sMessage.copy(), mn, e1.release(), e2.release(), e3.release()));
    eclear(dcode_sMessage(t.get()));
}

void _gen_message(er x, const char * tag, const char * str, ex E1, ex E2, ex E3, ex E4)
{
    uex e1(E1);
    uex e2(E2);
    uex e3(E3);
    uex e4(E4);
    assert(eis_sym(x));
    wex t(emake_str(tag));
    if (str != nullptr)
    {
        wex s(emake_str(str));
        esym_set_message_tag(x, t.get(), s.get(), 2);
    }
    ex mn = emake_node(gs.sym_sMessageName.copy(), ecopy(x), t.copy());
    t.reset(emake_node(gs.sym_sMessage.copy(), mn, e1.release(), e2.release(), e3.release(), e4.release()));
    eclear(dcode_sMessage(t.get()));
}


// wrong number of arguments
ex _handle_message_argct(er e)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    if (elength(e) == 1)
    {
        _gen_message(echild(e,0), "argctu", nullptr, ecopychild(e,0));
    }
    else
    {
        ex t = emake_int_ui(elength(e));
        _gen_message(echild(e,0), "argct", nullptr, ecopychild(e,0), t);
    }
    return ecopy(e);
}


// exact number of arguments
ex _handle_message_argx1(er e)
{
	return _handle_message_argx(e, 1);
}

ex _handle_message_argx2(er e)
{
	return _handle_message_argx(e, 2);
}

ex _handle_message_argx(er e, uint32_t n1)
{
    assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    uex t(emake_int_ui(elength(e)));
    if (n1 == 1)
    {
        _gen_message(echild(e,0), "argx", nullptr, ecopychild(e,0), t.release());        
    }
    else
    {
        ex t1 = emake_cint(n1); // n1 should be small
        if (elength(e) == 1)
            _gen_message(echild(e,0), "argr", nullptr, ecopychild(e,0), t1);
        else
            _gen_message(echild(e,0), "argrx", nullptr, ecopychild(e,0), t.release(), t1);
    }
    return ecopy(e);
}

// more arguments
ex _handle_message_argm(er e, uint32_t n1)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    uex t(emake_int_ui(elength(e)));
    ex t1 = emake_cint(n1);
    if (elength(e) == 1)
        _gen_message(echild(e,0), "argmu", nullptr, ecopychild(e,0), t1);
    else
        _gen_message(echild(e,0), "argm", nullptr, ecopychild(e,0), t.release(), t1);
    return ecopy(e);
}

// this or that arguments
ex _handle_message_argt(er e, uint32_t n1n2)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    uex t(emake_int_ui(elength(e)));
    ex t1 = emake_cint((n1n2>>0)&255);
    ex t2 = emake_cint((n1n2>>8)&255);
    if (elength(e) == 1)
        _gen_message(echild(e,0), "argtu", nullptr, ecopychild(e,0), t1, t2);        
    else
        _gen_message(echild(e,0), "argt", nullptr, ecopychild(e,0), t.release(), t1, t2);
    return ecopy(e);
}

// between arguments
ex _handle_message_argb(er e, uint32_t n1n2)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    uex t(emake_int_ui(elength(e)));
    ex t1 = emake_cint((n1n2>>0)&255);
    ex t2 = emake_cint((n1n2>>8)&255);
    if (elength(e) == 1)
        _gen_message(echild(e,0), "argbu", nullptr, ecopychild(e,0), t1, t2);
    else
        _gen_message(echild(e,0), "argb", nullptr, ecopychild(e,0), t.release(), t1, t2);
    return ecopy(e);
}

ex _handle_message_sym(er e, uint32_t n1)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    _gen_message(echild(e,0), "sym", nullptr, ecopychild(e,n1), emake_cint(n1));
    return ecopy(e);
}

ex _handle_message_intpm(er e, uint32_t n1)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    _gen_message(echild(e,0), "intpm", NULL, ecopy(e), emake_cint(n1));
    return ecopy(e);
}

ex _handle_message_intnm(er e, uint32_t n1)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    _gen_message(echild(e,0), "intnm", NULL, ecopy(e), emake_cint(n1));
    return ecopy(e);
}

ex _handle_message_intm(er e, uint32_t n1)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    _gen_message(echild(e,0), "intm", NULL, ecopy(e), emake_cint(n1));
    return ecopy(e);
}

ex _handle_message_ilsmn(er e, uint32_t n1)
{
	assert(eis_node(e));
	assert(eis_sym(echild(e,0)));
    _gen_message(echild(e,0), "ilsmn", nullptr, ecopy(e), emake_cint(n1));
    return ecopy(e);
}


ex dcode_sMessage(er e)
{
//std::cout << "dcode_sMessage: " << ex_tostring_full(e.get()) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sMessage.get()));

    if (elength(e) == 0)
    {
        _gen_message(gs.sym_sMessage.get(), "argm", NULL, gs.sym_sMessage.copy(), emake_cint(0), emake_cint(1));
        return gs.sym_sNull.copy();
    }

    if (!gs.muffler_stack.empty())
    {
        return gs.sym_sNull.copy();
    }

    er E1 = echild(e,1);
    if (ehas_head_sym_length(E1, gs.sym_sMessageName.get(), 2)
        && eis_sym(echild(E1,1)) && eis_str(echild(E1,2)))
    {
        int count = 0;
        for (size_t i = 0; i < gs.current_message_list.size(); i++)
        {
            if (ex_same(E1, gs.current_message_list[i].get()))
            {
                count++;
                if (count >= 3)
                {
                    return gs.sym_sNull.copy();
                }
            }
        }

        std::string xtag;
        xtag.append(esym_name_data(echild(E1,1)), esym_name_length(echild(E1,1)));
        xtag.push_back(':');
        xtag.push_back(':');
        xtag.append(estr_data(echild(E1,2)), estr_length(echild(E1,2)));
        xtag.push_back(':');
        xtag.push_back(' ');

        gs.current_message_list.push_back(uex(ecopy(E1)));
        er m = esym_find_message_tag(echild(E1,1), echild(E1,2));
        if (m == nullptr)
        {
            m = esym_find_message_tag(gs.sym_sGeneral.get(), echild(E1,2));
            if (m != nullptr)
            {
                esym_set_message_tag(echild(E1,1), echild(E1,2), eis_str(m) ? m : echild(m,1), eis_str(m) ? 0 : 1);
            }
        }
        uex r(ecopy(e));
        if (m == nullptr)
        {
            r.replacechild(1, emake_str("-- Message text not found --"));
        }
        else if (eis_str(m))
        {
            r.replacechild(1, ecopy(m));
        }
        else
        {
            return gs.sym_sNull.copy();
        }
        r.replacechild(0, gs.sym_sStringForm.copy());

        if (fp_out_type == FP_OUT_BINARY)
        {
            r.reset(ex_to_exbox_standard(r.get()));
            ex t = emake_str(xtag);
            r.reset(emake_node(gs.sym_sMessagePacket.copy(), t, r.release()));
            swrite_byte(fp_out, CMD_EXPR);
            swrite_ex(fp_out, r.get());
            fflush(fp_out);
        }
        else if (fp_out_type == FP_OUT_TEXT)
        {
            fwrite(xtag.c_str(), 1, xtag.length(), fp_out);
            fflush(fp_out);
            std::string s = ex_tostring(r.get());
            fwrite(s.c_str(), 1, s.length(), fp_out);
            fputc('\n', fp_out);
            fflush(fp_out);
        }
    }
    else
    {
        _gen_message(gs.sym_sMessage.get(), "name", "Message name `1` is not of the form symbol::name.", ecopy(E1));
    }

    return gs.sym_sNull.copy();
}
