#include "ex.h"
#include "globalstate.h"
#include "arithmetic.h"

ex ncode_sCatalan(er e, slong prec)
{
    if (!eis_sym(e, gs.sym_sCatalan.get()))
        return ecopy(e);

    ex z = emake_real();
    arb_const_catalan(ereal_data(z), prec);
    return z;
}

ex ncode_sChampernowneNumber(er e, slong prec)
{
//std::cout << "ncode_sChampernowneNumber(" << prec << "): " << e << std::endl;

    if (!ehas_head_sym(e, gs.sym_sChampernowneNumber.get()))
        return ecopy(e);

    const fmpz * B = eget_cint_data(10);

    if (elength(e) == 1)
    {
        if (!eis_int(echild(e,1)))
            return ecopy(e);
        B = eint_data(echild(e,1));

        if (fmpz_cmp_ui(B, 1) <= 0)
            return ecopy(e);
    }
    else if (elength(e) > 1)
    {
        return ecopy(e);
    }

    xarb_t z, s, d, u, v;
    xfmpz_t a, t0, t1, t2, t3, t4, t6;

    if (!fmpz_abs_fits_ui(B))
    {
        // only need one term, TODO add error
        fmpz_sub_ui(t0.data, B, 1);
        fmpz_mul(t2.data, t0.data, t0.data);
        arb_fmpz_div_fmpz(z.data, B, t2.data, prec + 1);
        return emake_real_move(z);
    }

    // target is abserror < 2^-p
    ulong b = fmpz_get_ui(B);
    ulong p = prec + 2 + fmpz_bits(B);

    double mlog2u = (b - 1)*log2(b);

    std::vector<xfmpz_t> ck;
    ck.push_back(xfmpz_t(1));
    ulong k = 1;
    while (true)
    {
        k++;
        fmpz_pow_ui(t0.data, B, k - 1);
        fmpz_mul_ui(t0.data, t0.data, k);
        fmpz_add(t0.data, t0.data, ck[k - 2].data);
        if (fmpz_cmp_ui(t0.data, 1 + (p + k + 2)/mlog2u) > 0)
            break;
        ck.push_back(t0);
    }
    ulong n = k - 1;

    arb_zero_pm_one(s.data);
    arb_mul_2exp_si(s.data, s.data, -slong(p + n + 2)); // tail error
    for (k = n; k != 0; k--)
    {
        ulong w1 = 5 + p + k;
        ulong w2 = mlog2u*fmpz_get_d(ck[k - 1].data);
        ulong w = 5 + (w1 > w2 ? w1 - w2 : 0);
        flint_bitcnt_t ckbits = fmpz_bits(ck[k - 1].data);
        fmpz_pow_ui(a.data, B, k);
        fmpz_sub_ui(t0.data, a.data, 1);
        fmpz_mul(t1.data, t0.data, a.data);
        fmpz_add_ui(t1.data, t1.data, 1);
        fmpz_mul(t2.data, t0.data, t0.data);
        fmpz_divexact_ui(t0.data, t0.data, b - 1);
        fmpz_mul(t3.data, t0.data, t0.data);
        fmpz_mul(t0.data, a.data, B);
        fmpz_sub_ui(t0.data, t0.data, 1);
        fmpz_mul(t4.data, t0.data, a.data);
        fmpz_add_ui(t4.data, t4.data, 1);
        fmpz_divexact_ui(t0.data, t0.data, b - 1);
        fmpz_mul(t6.data, t0.data, t0.data);
        fmpz_mul(t0.data, t2.data, t6.data);
        fmpz_mul(t2.data, t6.data, t1.data);
        fmpz_submul(t2.data, t3.data, t4.data);
        arb_fmpz_div_fmpz(d.data, t2.data, t0.data, w);
        arb_ui_pow_ui(u.data, b, b - 1, w + ckbits);
        arb_pow_fmpz(v.data, u.data, ck[k - 1].data, w + ckbits);
        arb_div(z.data, d.data, v.data, w);
        arb_add(s.data, s.data, z.data, w);
    }
    fmpz_set_ui(t2.data, b - 1);
    fmpz_mul_ui(t2.data, t2.data, b - 1);
    arb_fmpz_div_fmpz(z.data, B, t2.data, prec + 3);
    arb_sub(z.data, z.data, s.data, prec + 3);
    return emake_real_move(z);
}
