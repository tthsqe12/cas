#include "globalstate.h"
#include "code.h"
#include "timing.h"
#include "ex_print.h"
#include "ex_cont.h"


bool _coprimeq(er e, er f)
{
    if (eis_int(e) && eis_int(f))
        return fmpz_coprime(eint_data(e), eint_data(f));

    return false;
}

ex dcode_sCoprimeQ(er e)
{
//std::cout << "dcode_sCoprimeQ: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sCoprimeQ.get()));

    size_t n = elength(e);
    if (n < 1)
        return gs.sym_sFalse.copy();

    if (n <= 2)
        return emake_boole(_coprimeq(echild(e,1), echild(e,n)));

    for (size_t i = 0; i + 1 < n; i++)
    for (size_t j = i + 1; j < n; j++)
        if (!_coprimeq(echild(e,1+i), echild(e,1+j)))
            return gs.sym_sFalse.copy();

    return gs.sym_sTrue.copy();
}


ex dcode_sGCD(er e)
{
//std::cout << "dcode_sGCD: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sGCD.get()));

    size_t n = elength(e);
    if (n == 2)
    {
        er x = echild(e,1);
        er y = echild(e,2);
        if (eis_int(x) && eis_int(y))
        {
            fmpz_t z;
            fmpz_init(z);
            fmpz_gcd(z, eint_data(x), eint_data(y));
            return emake_int_clear(z);
        }
        else if (eis_int(x) && eis_rat(y))
        {
            ex z = emake_rat();
            _fmpq_gcd(fmpq_numref(erat_data(z)), fmpq_denref(erat_data(z)),
                    eint_data(x), eint_data(eget_cint(1)),
                    fmpq_numref(erat_data(y)), fmpq_denref(erat_data(y)));
            return efix_rat(z);
        }
        else if (eis_rat(x) && eis_int(y))
        {
            ex z = emake_rat();
            _fmpq_gcd(fmpq_numref(erat_data(z)), fmpq_denref(erat_data(z)),
                    eint_data(y), eint_data(eget_cint(1)),
                    fmpq_numref(erat_data(x)), fmpq_denref(erat_data(x)));
            return efix_rat(z);
        }
        else if (eis_rat(x) && eis_rat(y))
        {
            ex z = emake_rat();
            _fmpq_gcd(fmpq_numref(erat_data(z)), fmpq_denref(erat_data(z)),
                    fmpq_numref(erat_data(x)), fmpq_denref(erat_data(x)),
                    fmpq_numref(erat_data(y)), fmpq_denref(erat_data(y)));
            return efix_rat(z);
        }
        else
        {
            return ecopy(e);
        }
    }
    else
    {
        xfmpq_t t, g(0,1);
        for (size_t i = 1; i <= n; i++)
        {
            er x = echild(e,i);
            if (eis_int(x))
            {
                _fmpq_gcd(fmpq_numref(t.data), fmpq_denref(t.data),
                    fmpq_numref(g.data), fmpq_denref(g.data),
                    eint_data(x), eint_data(eget_cint(1)));
            }
            else if (eis_rat(x))
            {
                _fmpq_gcd(fmpq_numref(t.data), fmpq_denref(t.data),
                    fmpq_numref(g.data), fmpq_denref(g.data),
                    fmpq_numref(erat_data(x)), fmpq_denref(erat_data(x)));
            }
            else
            {
                return ecopy(e);
            }
            fmpq_swap(t.data, g.data);
        }
        return emake_rat_move(g.data);
    }
    return ecopy(e);
}

ex dcode_sExtendedGCD(er e)
{
    ulong n = elength(e);
    if (n == 2)
    {
        if (eis_int(echild(e,1)) && eis_int(echild(e,2)))
        {
            uex g(emake_int()), a(emake_int()), b(emake_int());
            fmpz_zero(eint_data(a.get())); // xgcd may leave output undefined
            fmpz_zero(eint_data(b.get())); //
            fmpz_xgcd(eint_data(g.get()), eint_data(a.get()), eint_data(b.get()), eint_data(echild(e,1)), eint_data(echild(e,2)));
            ex t = emake_node(gs.sym_sList.copy(), a.release(), b.release());
            return emake_node(gs.sym_sList.copy(), g.release(), t);
        }
        else
        {
            return ecopy(e);
        }
    }
    else
    {
        return ecopy(e);
    }
}

ex dcode_sLCM(er e)
{
//std::cout << "dcode_sLCM: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sLCM.get()));

    ulong n = elength(e);
    if (n == 2)
    {
        er x = echild(e,1);
        er y = echild(e,2);
        if (eis_int(x) && eis_int(y))
        {
            ex Z = emake_int();
            fmpz_lcm(eint_data(Z), eint_data(x), eint_data(y));
            return efix_int(Z);
        }
        else
        {
            return ecopy(e);
        }
    }
    else
    {
        xfmpz_t t, g;
        for (ulong i = 1; i <= n; i++)
        {
            er x = echild(e,i);
            if (eis_int(x))
                fmpz_lcm(t.data, g.data, eint_data(x));
            else
                return ecopy(e);

            fmpz_swap(t.data, g.data);
        }
        return emake_int_move(g);
    }
    return ecopy(e);
}
