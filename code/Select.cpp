#include "eval.h"
#include "code.h"
#include "ex_cont.h"

static ex _select(er e, er f, size_t limit, er def)
{
    if (eis_parray(e))
        e = eparray_get_normal(e);

    if (!eis_node(e))
    {
        _gen_message(gs.sym_sSelect.get(), "normal", nullptr, emake_cint(1), ecopy(e));
        return ecopy(def);
    }

    std::vector<wex> r;
    for (size_t i = 0; i < elength(e) && r.size() < limit; i++)
    {
        wex res(eval(emake_node(ecopy(f), ecopychild(e,i+1))));
        if (eis_sym(res.get(), gs.sym_sTrue.get()))
            r.emplace_back(ecopychild(e,i+1));
    }
    return emake_node(ecopychild(e,0), r);
}

ex dcode_sSelect(er e)
{
//std::cout << "dcode_sSelect: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_sSelect.get()));

    if (elength(e) == 1)
    {
        return ecopy(e);
    }
    else if (elength(e) == 2)
    {
        return _select(echild(e,1), echild(e,2), WORD_MAX, e);
    }
    else if (elength(e) == 3)
    {
        er n = echild(e,3);
        ulong limit = UWORD_MAX;

        if (!eis_infinity(n))
        {
            if (!eis_intn(n))
            {
                _gen_message(echild(e,0), "innf", NULL, emake_cint(3), ecopy(e));
                return ecopy(e);
            }

            if (fmpz_sgn(eint_data(n)) <= 0)
                limit = 0;
            if (eis_intsm(n))
                limit = eintsm_get(n);
        }

        return _select(echild(e,1), echild(e,2), limit, e);
    }
    else
    {
        return _handle_message_argb(e, (1 << 0) + (3 << 8));
    }
}
