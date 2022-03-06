#include "flintarb_wrappers.h"
#include "xfmpz_t.h"

std::ostream& operator<<(std::ostream& o, const fmpq_t a)
{
    o << fmpq_numref(a);
    if (!fmpz_is_one(fmpq_denref(a)))
        o << "/" << fmpq_denref(a);
    return o;
}

