#include <fstream>
#include "globalstate.h"
#include "code.h"
#include "ex_parse.h"
#include "eval.h"

ex dcode_sGet(er e)
{
//std::cout << "dcode_sGet: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sGet.get()));

    size_t n = elength(e);
    if (n < 1 || n > 3)
        return _handle_message_argb(e, 1 + (3 << 8));

    er e1 = echild(e,1);
    if (n == 1 && eis_str(e1))
    {
        wex_vector v;
        {
            std::ifstream is;
            is.open(estr_data(e1), std::ios::in);
            if ((is.rdstate() & std::ifstream::failbit) != 0)
            {
                _gen_message(echild(e,0), "noopen", NULL, ecopy(e1));
                return gs.sym_s$Failed.copy();
            }

            syntax_report sr;
            ex_parse_file(v, is, true, sr);
            if (sr.have_error())
            {
                _gen_message(gs.sym_sGeneral.get(), "sntx", "`1` near `2`.", sr.translate_error(), sr.near_error());
                return gs.sym_s$Failed.copy();
            }
        }

        assert(!v.empty());
        size_t i = 0;
        while (i + 1 < v.size())
            eclear(eval(v[i++].copy()));
        return eval(v[i].copy());
    }
    else
    {
        return ecopy(e);
    }
}
