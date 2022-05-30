#include "globalstate.h"
#include "code.h"
#include "ex_parse.h"
#include "ex_cont.h"

static bool _string_join(std::vector<wex> &l, std::vector<size_t> &bad, er e)
{
    size_t n = elength(e);
    bool changed = false;
    for (size_t i = 1; i <= n; i++)
    {
        er ei = echild(e,i);
        if (eis_str(ei))
        {
            if (l.empty() || !eis_str(l.back().get()))
            {
                l.emplace_back(ecopy(ei));
            }
            else
            {
                std::string s = estr_mkstring(l.back().get());
                s.append(estr_data(ei), estr_length(ei));
                l.back().reset(emake_str(s));
                changed = true;
            }
        }
        else if (ehas_head_sym(ei, gs.sym_sList.get()))
        {
            _string_join(l, bad, ei);
            changed = true;
        }
        else
        {
            l.emplace_back(ecopy(ei));
            bad.push_back(l.size());
        }
    }
    return changed;
}

ex dcode_sStringJoin(er e)
{
//std::cout << "dcode_sToCharacterCode: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sStringJoin.get()));

    std::vector<wex> l;
    std::vector<size_t> bad;
    bool changed = _string_join(l, bad, e);
    if (bad.empty())
    {
        if (l.empty())
            return emake_str("", 0);
        assert(l.size() == 1);
        assert(eis_str(l[0].get()));
        return l[0].copy();
    }
    else if (!changed)
    {
        for (size_t i = 0; i < bad.size(); i++)
        {
            ex t = emake_int_ui(bad[i]);
            _gen_message(gs.sym_sToCharacterCode.get(), "string", NULL, t, ecopy(e));
        }
        return ecopy(e);
    }
    else
    {
        return emake_node(gs.sym_sStringJoin.copy(), l);
    }
}

ex dcode_sStringLength(er e)
{
//std::cout << "dcode_sStringLength: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sStringLength.get()));

    return ecopy(e);
}

ex dcode_sToCharacterCode(er e)
{
//std::cout << "dcode_sToCharacterCode: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sToCharacterCode.get()));

    if (elength(e) != 1)
    {
        return _handle_message_argx1(e);
    }

    er X = echild(e,1);
    if (eis_str(X))
    {
        std::vector<wex> v;
        size_t n = estr_length(X);
        const char * a = estr_data(X);
        for (size_t i = 0; i < n;)
        {
            char16_t c;
            i += readonechar16(c, a + i);
            v.emplace_back(emake_int_ui(c));
        }
        return emake_node(gs.sym_sList.copy(), v);
    }
    else if (ehas_head_sym(X, gs.sym_sList.get()))
    {
        uex l; l.init_push_backr(gs.sym_sList.get(), elength(X));
        for (size_t j = 1; j <= elength(X); j++)
        {
            er Y = echild(X, j);
            if (!eis_str(Y))
            {
                _gen_message(gs.sym_sToCharacterCode.get(), "strse", NULL, ecopy(e));
                return ecopy(e);
            }
            std::vector<wex> v;
            size_t n = estr_length(Y);
            const char * a = estr_data(Y);
            for (size_t i = 0; i < n;)
            {
                char16_t c;
                i += readonechar16(c, a + i);
                v.emplace_back(emake_int_ui(c));
            }
            l.push_back(emake_node(gs.sym_sList.copy(), v));
        }
        return l.release();
    }
    else
    {
        _gen_message(gs.sym_sToCharacterCode.get(), "strse", NULL, ecopy(e));
        return ecopy(e);
    }
}


static ex _bad_character_code(er e, ex pos, ex f)
{
    _gen_message(gs.sym_sFromCharacterCode.get(), "notunicode", "A character code, which should be a non-negative integer less than 65536, is expected at position `1` in `2`.", pos, f);
    return ecopy(e);
}

ex dcode_sFromCharacterCode(er e)
{
//std::cout << "dcode_sFromCharacterCode: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sFromCharacterCode.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    std::string s;
    er e1 = echild(e,1);
    if (eis_leaf(e1))
    {
        if (eis_int(e1) && 0 <= *eint_data(e1) && *eint_data(e1) < 65536)
        {
            stdstring_pushback_char16(s, *eint_data(e1));
            return emake_str(s);
        }
        else
        {
            return _bad_character_code(e, emake_cint(1), ecopy(e));
        }
    }
    else if (ehas_head_sym(e1, gs.sym_sList.get()))
    {
        for (ulong i = 0; i < elength(e1); i++)
        {
            er ei = echild(e1, i+1);
            if (eis_int(ei) && 0 <= *eint_data(ei) && *eint_data(ei) < 65536)
            {
                stdstring_pushback_char16(s, *eint_data(ei));
            }
            else if (ehas_head_sym(ei, gs.sym_sList.get()))
            {
                uex r; r.init_push_backr(gs.sym_sList.get(), elength(e1));
                for (ulong j = 0; j < elength(e1); j++)
                {
                    er ej = echild(e1, j+1);
                    s.clear();
                    if (eis_int(ej) && 0 <= *eint_data(ej) && *eint_data(ej) < 65536)
                    {
                        stdstring_pushback_char16(s, *eint_data(ej));                        
                    }
                    else if (ehas_head_sym(ej, gs.sym_sList.get()))
                    {
                        for (ulong k = 0; k < elength(ej); k++)
                        {
                            er ejk = echild(ej, k+1);
                            if (eis_int(ejk) && 0 <= *eint_data(ejk) && *eint_data(ejk) < 65536)
                            {
                                stdstring_pushback_char16(s, *eint_data(ejk));
                            }
                            else
                            {
                                return _bad_character_code(e, emake_node(gs.sym_sList.copy(), emake_int_ui(j+1), emake_int_ui(k+1)), ecopy(e1));
                            }
                        }
                    }
                    else
                    {
                        return _bad_character_code(e, emake_int_ui(j+1), ecopy(e1));
                    }
                    r.push_back(emake_str(s));
                }
                return r.release();
            }
            else
            {
                return _bad_character_code(e, emake_int_ui(i+1), ecopy(e1));
            }
        }
        return emake_str(s);
    }
    else if (eis_parray_fmpz(e1) && eto_parray_fmpz(e1).dimensions.rank == 1)
    {
        fmpz * a = eto_parray_fmpz(e1).array->data;
        slong n = eto_parray_fmpz(e1).dimensions.get_index(0);
        for (slong i = 0; i < n; i++)
        {
            if (0 <= a[i] && a[i] < 65536)
            {
                stdstring_pushback_char16(s, a[i]);
            }
            else
            {
                return _bad_character_code(e, emake_cint(i+1), ecopy(e1));
            }
        }
        return emake_str(s);
    }
    else if (eis_parray_fmpz(e1) && eto_parray_fmpz(e1).dimensions.rank == 1)
    {
        fmpz * a = eto_parray_fmpz(e1).array->data;
        slong n = eto_parray_fmpz(e1).dimensions.get_index(0);
        slong m = eto_parray_fmpz(e1).dimensions.get_index(1);
        uex r; r.init_push_backr(gs.sym_sList.get(), n);
        for (ulong i = 0; i < n; i++)
        {
            s.clear();
            for (ulong j = 0; j < m; j++)
            {
                if (0 <= a[m*i+j] && a[m*i+j] < 65536)
                {
                    stdstring_pushback_char16(s, a[m*i+j]);
                }
                else
                {
                    return _bad_character_code(e, emake_node(gs.sym_sList.copy(), emake_cint(i+1), emake_cint(j+1)), ecopy(e1));
                }
            }
            r.push_back(emake_str(s));
        }
        return r.release();
    }
    else
    {
        return _bad_character_code(e, emake_cint(1), ecopy(e));
    }
}
