#include "xacb_t.h"
#include "xacb_mat_t.h"

std::ostream& operator<<(std::ostream& o, const acb_mat_t a)
{
    o << "{";
    for (slong i = 0; i < acb_mat_nrows(a); i++)
    {
        if (i > 0)
            o << ", ";

        o << "{";
        for (slong j = 0; j < acb_mat_ncols(a); j++)
        {
            if (j > 0)
                o << ", ";

            o << acb_mat_entry(a, i, j);
        }
        o << "}";
    }
    o << "}";
    return o;
}
