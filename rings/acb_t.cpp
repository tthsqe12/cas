#include "flintarb_wrappers.h"
#include "xfmpq_t.h"
#include "rarb_t.h"
#include "racb_t.h"
#include <complex>
#include "flint/arith.h"

std::ostream& operator<<(std::ostream& o, const acb_t a)
{
    o << acb_realref(a);
    if (!arb_is_zero(acb_imagref(a)))
    {
        o << " + I*";
        o << acb_imagref(a);
    }
    return o;
}

template <>
std::complex<double> racb_t::convert<std::complex<double>>(const xacb_t& a)
{
    return std::complex<double>(arf_get_d(arb_midref(acb_realref(a.data)), ARF_RND_DOWN),
                                arf_get_d(arb_midref(acb_imagref(a.data)), ARF_RND_DOWN));
}


void racb_t::cscpi_series(xacb_poly_t& x, slong ord)
{
    acb_poly_fit_length(x.data, ord);
    xfmpz_t u;
    xfmpq_t v;
    for (slong n = 0; 2*n < ord; n++)
    {
        arb_zero(acb_imagref(x.data->coeffs + 2*n));
        arb_ptr t = acb_realref(x.data->coeffs + 2*n);
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
        acb_zero(x.data->coeffs + 2*n + 1);
    _acb_poly_set_length(x.data, ord);
    _acb_poly_normalise(x.data);
}
