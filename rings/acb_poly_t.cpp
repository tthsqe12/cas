#include "xacb_poly_t.h"

std::ostream& operator<<(std::ostream& o, const acb_poly_t a)
{
    bool first = true;
    for (slong i = 0; i < a->length; i++)
    {
        if (!first)
            o << " + ";
        first = false;
        o << "(" << a->coeffs + i << ")*#^" << i;
    }
    if (first)
        o << "0";
    return o;
}
