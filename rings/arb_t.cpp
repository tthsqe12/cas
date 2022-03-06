#include "flintarb_wrappers.h"
#include "xfmpq_t.h"
#include "rarb_t.h"
#include <complex>
#include "flint/arith.h"

std::ostream& operator<<(std::ostream& o, const arb_t a)
{
    xfstr s(arb_get_str(a, 40, 0));
    o << s.data;
    return o;
}

template <>
std::complex<double> rarb_t::convert<std::complex<double>>(const xarb_t& a)
{
    return std::complex<double>(arf_get_d(arb_midref(a.data), ARF_RND_DOWN), 0.0);
}


void rarb_t::cscpi_series(xarb_poly_t& x, slong ord)
{
    arb_poly_fit_length(x.data, ord);
    xfmpz_t u;
    xfmpq_t v;
    for (slong n = 0; 2*n < ord; n++)
    {
        arb_ptr t = x.data->coeffs + 2*n;
        arb_const_pi(t, prec);
        arb_pow_ui(t, t, 2*n, prec);
        fmpz_fac_ui(u.data, 2*n);
        arith_bernoulli_number(v.data, 2*n);
        fmpq_div_fmpz(v.data, v.data, u.data);
        fmpz_one(u.data);
        fmpz_mul_2exp(u.data, u.data, 2*n);
        fmpz_sub_ui(u.data, u.data, 2);
        fmpq_mul_fmpz(v.data, v.data, u.data);
        arb_mul_fmpz(t, t, fmpq_numref(v.data), prec);
        arb_div_fmpz(t, t, fmpq_denref(v.data), prec);
        if (!(n&1))
            arb_neg(t, t);
    }
    for (slong n = 0; 2*n + 1 < ord; n++)
        arb_zero(x.data->coeffs + 2*n + 1);
    _arb_poly_set_length(x.data, ord);
    _arb_poly_normalise(x.data);
}
