#include "globalstate.h"
#include "code.h"
#include "ex_cont.h"

ex dcode_sLength(er e)
{
//std::cout << "dcode_sLength: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sLength.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    size_t r;
    er x = echild(e,1);
    if (eis_node(x))
    {
        return emake_int_ui(elength(x));
    }
    else if (eis_packed(x))
    {
        assert(eis_parray(x));
        return emake_int_si(eparray_dimension(x,0));
    }
    else if (eis_hmap(x))
    {
        return emake_int_ui(eto_hmap(x).map_data.size());
    }
    else
    {
        return emake_cint(0);
    }
}
