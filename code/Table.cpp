#include "globalstate.h"
#include "code.h"
#include "eval.h"
#include "local.h"
#include "ex_cont.h"
#include "assign.h"

/* return
    nullptr for failure (bad input)
*/
static ex dcode_sTable_nested(er e, size_t pos)
{
//std::cout << "dcode_sTable_nested: " << ex_tostring_full(e) << " pos: " << pos << std::endl;
    assert(ehas_head_sym(e, gs.sym_sTable.get()));
    size_t n = elength(e);
    assert(n >= 2);
    assert(pos <= n);

    er it = echild(e, pos);
    if (ehas_head_sym(it, gs.sym_sList.get()) && elength(it) >= 2 && elength(it) <= 4)
    {
        uex var(ecopychild(it,1));
        uex sym;
        int value_type = 0;
        if (eis_sym(var.get()))
        {
            sym.setz(var.copy());
        }
        else if (eis_node(var.get()))
        {
            var.setz(eval(var.release()));
            if (eis_sym(var.get()))
            {
                sym.setz(var.copy());
            }
            else if (eis_node(var.get()))
            {
                if (eis_sym(var.child(0)))
                {
                    value_type = 1;
                    sym.setz(var.copychild(0));
                }
                else
                {
                    value_type = 2;
                    sym.setz(ecopy(esymbolic_head(var.child(0))));
                }
            }
            else
            {
                _gen_message(gs.sym_sTable.get(), "itraw", NULL, var.copy());
                return nullptr;
            }
        }
        else
        {
            _gen_message(gs.sym_sTable.get(), "itraw", NULL, var.copy());
            return nullptr;
        }

        // sym is the symbol to localize
        // variable the value to set
        // value_type = 0 -> ovalue
        // value_type = 1 -> dvalue
        // value_type = 2 -> svalue
        sym_localize L(sym.get());

        uex lower(eval(ecopychild(it,2)));
        if (elength(it) == 2 && ehas_head_sym(lower.get(), gs.sym_sList.get()))
        {
            size_t count = 0;
            uex r; r.init_push_backr(gs.sym_sList.get(), elength(lower.get()));
            while (++count <= elength(lower.get()))
            {
                if (value_type == 0)
                    esym_assign_ovaluex(sym.get(), lower.copychild(count));
                else if (value_type == 1)
                    esym_assign_dvaluex(sym.get(), var.copy(), lower.copychild(count));
                else if (value_type == 2)
                    esym_assign_svaluex(sym.get(), var.copy(), lower.copychild(count));

                if (pos >= n)
                {
                    assert(pos == n);
                    r.push_back(eval(ecopychild(e,1)));
                }
                else
                {
                    ex s = dcode_sTable_nested(e, pos + 1);
                    if (s != nullptr)
                    {
                        r.push_back(s);
                    }
                    else
                    {
                        return nullptr;
                    }
                }
            }
			r.setnz(econvert_parray_shallow(r.get()));
            return r.release();
        }
        else
        {
            uex upper(emake_cint(1));
            uex step(emake_cint(1));
            if (elength(it) == 2)
            {
                upper.swap(lower);
            }
            else
            {
                assert(elength(it) >= 3);
                upper.reset(eval(ecopychild(it,3)));
                if (elength(it) != 3)
                {
                    assert(elength(it) == 4);
                    step.reset(eval(ecopychild(it,4)));
                }
            }
            uex scount(ex_subx(upper.copy(), lower.copy()));
            scount.setz(ex_divx(scount.release(), step.copy()));
            scount.setz(ex_floor(scount.release()));
//std::cout << "scount after floor: " << ex_tostring(scount.get()) << std::endl;
            if (!eis_int(scount.get()) || !fmpz_fits_si(scount.int_data()))
            {
                _gen_message(gs.sym_sTable.get(), "iterb", NULL, ecopychild(e, pos));
                return nullptr;
            }
            slong countlimit = fmpz_get_si(scount.int_data());
            slong count = -2;
            uex r; r.init_push_backr(gs.sym_sList.get(), std::max((slong)(0), countlimit + 1));
            while (++count < countlimit)
            {
                ex itval = ex_addx(lower.copy(), ex_mul_six(step.copy(), count + 1));
                if (value_type == 0)
                    esym_assign_ovaluex(sym.get(), itval);
                else if (value_type == 1)
                    esym_assign_dvaluex(sym.get(), var.copy(), itval);
                else if (value_type == 2)
                    esym_assign_svaluex(sym.get(), var.copy(), itval);

                if (pos >= n)
                {
                    assert(pos == n);
                    r.push_back(eval(ecopychild(e,1)));
                }
                else
                {
                    ex s = dcode_sTable_nested(e, pos + 1);
                    if (s != nullptr)
                        r.push_back(s);
                    else
                        return nullptr;
                }
            }
			r.setnz(econvert_parray_shallow(r.get()));
            return r.release();
        }
    }

//std::cout << "got here!!!!!!!!!" << std::endl;

    uex scount(eval(ecopy(it)));
    if (ehas_head_sym(scount.get(), gs.sym_sList.get()))
    {
        if (elength(scount.get()) != 1)
            return nullptr;
        scount.setnz(scount.copychild(1));
    }

    scount.setz(ex_floor(scount.release()));
//std::cout << "scount after floor: " << ex_tostring(scount.get()) << std::endl;
    if (!eis_intsm(scount.get()))
    {
        _gen_message(gs.sym_sTable.get(), "iterb", NULL, ecopychild(e, pos));
        return nullptr;
    }

    slong countlimit = fmpz_get_si(scount.int_data());
    slong count = -1;
//std::cout << "countlimit: " << countlimit<< std::endl;
    uex r; r.init_push_backr(gs.sym_sList.get(), std::max((slong)(0), countlimit));
    while (++count < countlimit)
    {
//std::cout << "evaluating table index: " << count + 1 << std::endl;
        if (pos >= n)
        {
            assert(pos == n);

            ex t = eval(ecopychild(e,1));

//std::cout << "pushing back t = " << ex_tostring_full(etor(t)) << std::endl;
            r.push_back(t);
//printf("pushed back\n");
        }
        else
        {
            ex s = dcode_sTable_nested(e, pos + 1);
            if (s != nullptr)
            {
                r.push_back(s);
            }
            else
            {
                return nullptr;
            }
        }
    }

	r.setnz(econvert_parray_shallow(r.get()));
    return r.release();
}

/*
    Table[expr, n]
    Table[expr, {i,imax}]
    Table[expr, {i,imin,imax}]
    Table[expr, {i,imin,imax,istep}]
    Table[expr, {i,{i1,i2,...,in}}]
*/
ex dcode_sTable(er e)
{
//std::cout << "dcode_sTable: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sTable.get()));
    size_t n = elength(e);

    if (n < 2)
    {
        if (n == 0)
        {
            return ecopy(e);
        }
        return eval(ecopychild(e,1));
    }

    ex r = dcode_sTable_nested(e, 2);
//std::cout << "table returning " << ex_tostring_full(etor(r)) << std::endl;
    return (r != nullptr) ? r : ecopy(e);
}
