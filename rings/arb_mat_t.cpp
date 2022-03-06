#include "xarb_t.h"
#include "xarb_mat_t.h"

std::ostream& operator<<(std::ostream& o, const arb_mat_t a)
{
    o << "{";
    for (slong i = 0; i < arb_mat_nrows(a); i++)
    {
        if (i > 0)
            o << ", ";

        o << "{";
        for (slong j = 0; j < arb_mat_ncols(a); j++)
        {
            if (j > 0)
                o << ", ";

            o << arb_mat_entry(a, i, j);
        }
        o << "}";
    }
    o << "}";
    return o;
}
