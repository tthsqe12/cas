#include "globalstate.h"
#include "code.h"

ex dcode_iComparePatterns(er e)
{
//std::cout << "dcode_iComparePatterns: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_iComparePatterns.get()));

    if (elength(e) != 2)
        return _handle_message_argx2(e);

    assert(false);
    return nullptr;
}
