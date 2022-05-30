#include "globalstate.h"
#include "code.h"
#include "eval.h"
#include "assign.h"
#include "local.h"

/* return
    0 for success
    1 for break
    2 for failure (bad input)
*/
static int dcode_sDo_nested(er e, size_t pos)
{
//std::cout << "dcode_sDo_nested: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDo.get()));
    size_t n = elength(e);
    assert(n >= 2);
    assert(pos <= n);

    er it = echild(e, pos);
    if (ehas_head_sym(it, gs.sym_sList.get()) && elength(it) >= 2 && elength(it) <= 4)
    {
        uex var(ecopychild(it,1));
        uex sym;
        ex symv;
        int value_type = 0;
        if (eis_sym(var.get()))
        {
            symv = var.copy();
        }
        else if (eis_node(var.get()))
        {
            var.setz(eval(var.release()));
            if (eis_sym(var.get()))
            {
                symv = var.copy();
            }
            else if (eis_node(var.get()))
            {
                if (eis_sym(var.child(0)))
                {
                    value_type = 1;
                    symv = var.copychild(0);
                }
                else
                {
                    value_type = 2;
                    symv = ecopy(esymbolic_head(var.child(0)));
                }
            }
            else
            {
                _gen_message(gs.sym_sDo.get(), "itraw", NULL, var.copy());
                return 2;
            }
        }
        else
        {
            _gen_message(gs.sym_sDo.get(), "itraw", NULL, var.copy());
            return 2;
        }

        sym.setz(symv);

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
                    try
                    {
                        eclear(eval(ecopychild(e, 1)));
                    }
                    catch (const exception_sym_sBreak &)
                    {
                        return 1;
                    }
                    catch (const exception_sym_sContinue &)
                    {
                    }
                }
                else
                {
                    int retcode = dcode_sDo_nested(e, pos + 1);
                    if (retcode != 0)
                        return retcode;
                }
            }
            return 0;
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
                _gen_message(gs.sym_sDo.get(), "iterb", NULL, emake_node(gs.sym_sList.copy(), var.copy(), lower.copy(), upper.copy(), step.copy()));
                return 2;
            }
            slong countlimit = fmpz_get_si(scount.int_data());
            slong count = -2;
//std::cout << "countlimit: " << countlimit<< std::endl;
            while (++count < countlimit)
            {
                ex itval;
                if (eis_one(step.copy()))
                    itval = ex_add_er_si(lower.get(), count + 1);
                else
                    itval = ex_addx(lower.copy(), ex_mul_six(step.copy(), count + 1));

                if (value_type == 0)
                    esym_assign_ovaluex(sym.get(), itval);
                else if (value_type == 1)
                    esym_assign_dvaluex(sym.get(), var.copy(), itval);
                else if (value_type == 2)
                    esym_assign_svaluex(sym.get(), var.copy(), itval);

                if (pos >= n)
                {
                    try
                    {
                        eclear(eval(ecopychild(e, 1)));
                    }
                    catch (const exception_sym_sBreak &)
                    {
                        return 1;
                    }
                    catch (const exception_sym_sContinue &)
                    {
                    }
                }
                else
                {
                    int retcode = dcode_sDo_nested(e, pos + 1);
//std::cout<<"print pos: "<< pos << " retcode:" << retcode << std::endl;
                    if (retcode != 0)
                        return retcode;
                }
            }
            return 0;
        }
    }

    if (ehas_head_sym(it, gs.sym_sList.get()))
    {
        if (elength(it) != 1)
            return 2;
        it = echild(it, 1);
    }
    uex scount(ex_floor(ecopy(it)));
//std::cout << "scount after floor: " << ex_tostring(scount.get()) << std::endl;
    if (!eis_int(scount.get()) || !fmpz_fits_si(scount.int_data()))
    {
        _gen_message(gs.sym_sDo.get(), "iterb", NULL, ecopychild(e, pos));
        return 2;
    }
    slong countlimit = fmpz_get_si(scount.int_data());
    slong count = -1;
//std::cout << "countlimit: " << countlimit<< std::endl;
    while (++count < countlimit)
    {
//std::cout << "count: " << count<< std::endl;

        if (pos >= n)
        {
            try
            {
                eclear(eval(ecopychild(e, 1)));
            }
            catch (const exception_sym_sBreak &)
            {
                return 1;
            }
            catch (const exception_sym_sContinue &)
            {
            }
        }
        else
        {
            int retcode = dcode_sDo_nested(e, pos + 1);
//std::cout<<"print pos: "<< pos << " retcode:" << retcode << std::endl;
            if (retcode != 0)
                return retcode;
        }
    }

    return 0;
}

/*
    Do[expr, n]
    Do[expr, {i,imax}]
    Do[expr, {i,imin,imax}]
    Do[expr, {i,imin,imax,istep}]
    Do[expr, {i,{i1,i2,...,in}}]
*/
ex dcode_sDo(er e)
{
//std::cout << "dcode_sDo: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDo.get()));
    size_t n = elength(e);

    if (n < 2)
    {
        if (n == 0)
        {
            return ecopy(e);
        }
        return eval(ecopychild(e,1));
    }

    int retcode = dcode_sDo_nested(e, 2);
    return retcode != 2 ? gs.sym_sNull.copy() : ecopy(e);
}
