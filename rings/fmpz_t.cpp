#include "flintarb_wrappers.h"

std::ostream& operator<<(std::ostream& o, const fmpz_t a)
{
    xfstr s(fmpz_get_str(NULL, 10, a));
    o << std::string(s.data);
    return o;
}
