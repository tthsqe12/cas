#include <iostream>
#include "xfmpzi_t.h"
#include "xfmpqi_t.h"

std::ostream& operator<<(std::ostream& o, const xfmpzi_t& a)
{
    o << a.re;
    if (!fmpz_is_zero(a.im.data))
        o << " + I*" << a.im;

    return o;
}

std::ostream& operator<<(std::ostream& o, const xfmpqi_t& a)
{
    if (fmpz_is_one(a.den.data))
        return o << a.num;

    o << a.num.re << "/" << a.den;
    if (!fmpz_is_zero(a.num.im.data))
        o << " + I*" << a.num.im << "/" << a.den;

    return o;
}

void fmpzi_divexact(
    fmpz_t xre, fmpz_t xim,
    const fmpz_t are, const fmpz_t aim,
    const fmpz_t bre, const fmpz_t bim,
    fmpz_t t, fmpz_t s)
{
    xfmpz_t r;
    fmpz_mul(s, bre, bre);
    fmpz_addmult(s, bim, bim, t);

    fmpz_mul(xre, are, bre);
    fmpz_addmult(xre, aim, bim, t);
    fmpz_tdiv_qr(xre, r.data, xre, s);
    assert(fmpz_is_zero(r.data));

    fmpz_mul(xim, aim, bre);
    fmpz_submult(xim, are, bim, t);
    fmpz_tdiv_qr(xim, r.data, xim, s);
    assert(fmpz_is_zero(r.data));
}



// return k with canonical_unit(a) = i^-k
ulong fmpzi_canonical_unit_i_pow(const fmpz_t ax, const fmpz_t ay)
{
    int s = fmpz_cmp(ax, ay);
    if (s == 0)
    {
        int t = fmpz_sgn(ax);
        return 2*ulong(t < 0);
    }
    else
    {
        int t = fmpz_cmpabs(ax, ay);
        return ulong(t <= 0) + 2*ulong(s <= 0);
    }
}

void fmpzi_mul_i_pow(fmpz_t zx, fmpz_t zy, ulong k)
{
    if (k & 1)
    {
        fmpz_neg_inplace((k & 2) ? zx : zy);
        fmpz_swap(zx, zy);
    }
    else if (k & 2)
    {
        fmpz_neg_inplace(zx);
        fmpz_neg_inplace(zy);
    }
}

void _shortest_l_infinity(
    fmpz_t u1, fmpz_t u2,
    fmpz_t t1, fmpz_t t2,
    const fmpz_t c, const fmpz_t b, const fmpz_t a)
{
    assert(fmpz_cmp_si(c, 0) > 0);
    assert(fmpz_cmp(a, b) > 0);
    assert(fmpz_cmp_si(b, 0) >= 0);

    fmpz_add(u1, b, c);
    fmpz_sub(u2, b, c);
    if (fmpz_cmp(a, c) <= 0)
    {
        fmpz_zero(u1);
        fmpz_set(u2, a);
        fmpz_zero(t1);
        fmpz_one(t2);
        return;
    }
    else if (fmpz_sgn(u2) <= 0)
    {
        fmpz_set(u1, c);
        fmpz_set(u2, b);
        fmpz_one(t1);
        fmpz_zero(t2);
        return;
    }
    else if (fmpz_cmp(a, u1) <= 0)
    {
        fmpz_set(u1, c);
        fmpz_sub(u2, b, a);
        fmpz_one(t1);
        fmpz_set_si(t2, -1);
        return;
    }

    _fmpq_cfrac_list_t s{NULL, -1, 0, WORD_MAX, 0, 0};
    _fmpz_mat22_t m{1, 0, 0, 1, 1};
    _fmpq_ball_t x;
    fmpz_init_set(x->left_num, a);
    fmpz_init(x->left_den); fmpz_swap(x->left_den, u1);
    fmpz_init_set(x->right_num, a);
    fmpz_init(x->right_den); fmpz_swap(x->right_den, u2);
    x->exact = 0;
    _fmpq_ball_get_cfrac(s, m, 1, x);
    _fmpq_ball_clear(x);

    xfmpz_t v11, v12, v21, v22, Q;
    fmpz_mul(v11.data, m->_11, c);
    fmpz_mul(v21.data, m->_11, c);
    fmpz_fmms(v12.data, m->_11, b, m->_21, a);
    fmpz_fmms(v22.data, m->_12, b, m->_22, a);
    int vcmp = fmpz_cmpabs(v11.data, v12.data);

    fmpz_set(u1, v11.data);
    fmpz_set(u2, v12.data);
    fmpz_set(t1, m->_11);
    fmpz_neg(t2, m->_21);

    int triesleft = 2;
    while (--triesleft >= 0 && vcmp < 0)
    {
        fmpz_tdiv_q(Q.data, v22.data, v12.data);
        assert(fmpz_cmp_si(Q.data, 0) < 0);
        assert(fmpz_cmpabs(u1, u2) < 0);
        fmpz_submul(m->_12, m->_11, Q.data); fmpz_swap(m->_12, m->_11);
        fmpz_submul(m->_22, m->_21, Q.data); fmpz_swap(m->_22, m->_21);
        fmpz_submul(v21.data, v11.data, Q.data); fmpz_swap(v21.data, v11.data);
        fmpz_submul(v22.data, v12.data, Q.data); fmpz_swap(v22.data, v12.data);
        vcmp = fmpz_cmpabs(v11.data, v12.data);
        if (fmpz_cmpabs(vcmp < 0 ? v12.data : v11.data, u2) < 0)
        {
            fmpz_set(u1, v11.data);
            fmpz_set(u2, v12.data);
            fmpz_set(t1, m->_11);
            fmpz_neg(t2, m->_21);
        }
    }

    fmpz_clear(m->_11);
    fmpz_clear(m->_12);
    fmpz_clear(m->_21);
    fmpz_clear(m->_22);
}

void fmpzi_gcd_fmpz(
    fmpz_t gx, fmpz_t gy,
    const fmpz_t Ax, const fmpz_t Ay,
    const fmpz_t b)
{
    if (fmpz_is_zero(b))
    {
        fmpz_set(gx, Ax);
        fmpz_set(gy, Ay);
        fmpzi_unit_canonicalize(gx, gy);
        return;
    }
    xfmpz_t A, B, C, ax, ay, ga, ua, va, g, u, v, axog, ayog, m1, m2, m3, m4;
    // ax ay           g*1 g*B
    // -ay ax  row ops   0  g*A
    //  b  0   ------>   0   0
    //  0  b             0   0
//  ax = cmpabs(a.x, b) < 0 ? a.x : smod(a.x, b)
//  ay = cmpabs(a.y, b) < 0 ? a.y : smod(a.y, b)
    fmpz_smod(ax.data, Ax, b);
    fmpz_smod(ay.data, Ay, b);
//  if iszero(ax)
//    return fmpzi(gcd(ay, b), zero(fmpz))
//  elseif iszero(ay)
//    return fmpzi(gcd(ax, b), zero(fmpz))
//  end



    if (fmpz_is_zero(ax.data) || fmpz_is_zero(ay.data))
    {
        fmpz_gcd(gx, fmpz_is_zero(ax.data) ? ay.data : ax.data, b);
        fmpz_zero(gy);
        return;
    }

//  ga, ua, va = gcdx(ax, ay)
    fmpz_xgcd(ga.data, ua.data, va.data, ax.data, ay.data);
//  g, u, _ = gcdx(ga, b)
    fmpz_xgcd(g.data, u.data, v.data, ga.data, b);
//  axog = _divexact(ax, g)
    fmpz_divexact(axog.data, ax.data, g.data);
//  ayog = _divexact(ay, g)
    fmpz_divexact(ayog.data, ay.data, g.data);
//  m1 = ayog*ua - axog*va
    fmpz_fmms(m1.data, ayog.data, ua.data, axog.data, va.data);
//  m2 = _divexact(ax,ga)*axog + _divexact(ay,ga)*ayog
    fmpz_fmma(m2.data, ax.data, axog.data, ay.data, ayog.data);
    fmpz_divexact(m2.data, m2.data, ga.data);
//  A = gcd(m2, _divexact(b, g))
    fmpz_divexact(B.data, b, g.data);
    fmpz_gcd(A.data, m2.data, B.data);
//  v, _ = _shortest_l∞(fmpz(1), mod(m1*u, A), A)
    fmpz_one(C.data);
    fmpz_mul(B.data, m1.data, u.data);
    fmpz_mod(B.data, B.data, A.data);

    _shortest_l_infinity(gx, gy, u.data, v.data, C.data, B.data, A.data);
//  z = isone(g) ? fmpzi(v[1], v[2]) : fmpzi(v[1]*g, v[2]*g)
    fmpz_mul(gx, gx, g.data);
    fmpz_mul(gy, gy, g.data);
//  return unit_canonicalize!(z)
    fmpzi_unit_canonicalize(gx, gy);
//end
}

void fmpzi_gcd(
    fmpz_t gx, fmpz_t gy,
    const fmpz_t ax, const fmpz_t ay,
    const fmpz_t bx, const fmpz_t by)
{
    if (fmpz_is_zero(ax))
    {
        fmpzi_gcd_fmpz(gx, gy, bx, by, ay);
        return;
    }
    else if (fmpz_is_zero(ay))
    {
        fmpzi_gcd_fmpz(gx, gy, bx, by, ax);
        return;
    }
    else if (fmpz_is_zero(bx))
    {
        fmpzi_gcd_fmpz(gx, gy, ax, ay, by);
        return;
    }
    else if (fmpz_is_zero(by))
    {
        fmpzi_gcd_fmpz(gx, gy, ax, ay, bx);
        return;
    }

    xfmpz_t A, B, C, axog, ayog, bxog, byog, ga, ua, va, gb, ub, vb, g, u, v, t, m1, m2, m3, m4;
    //  ax ay           g*1 g*B
    // -ay ax  row ops   0  g*A
    //  bx by  ------>   0   0
    // -by bx            0   0
//  ga, ua, va = gcdx(a.x, a.y)
    fmpz_xgcd(ga.data, ua.data, va.data, ax, ay);
//  gb, ub, vb = gcdx(b.x, b.y)
    fmpz_xgcd(gb.data, ub.data, vb.data, bx, by);
//  g, u, v = gcdx(ga, gb)
    fmpz_xgcd(g.data, u.data, v.data, ga.data, gb.data);
//  axog = _divexact(a.x, g)
    fmpz_divexact(axog.data, ax, g.data);
//  ayog = _divexact(a.y, g)
    fmpz_divexact(ayog.data, ay, g.data);
//  bxog = _divexact(b.x, g)
    fmpz_divexact(bxog.data, bx, g.data);
//  byog = _divexact(b.y, g)
    fmpz_divexact(byog.data, by, g.data);
//  m1 = ayog*ua - axog*va
    fmpz_fmms(m1.data, ayog.data, ua.data, axog.data, va.data);
//  m2 = _divexact(a.x,ga)*axog + _divexact(a.y,ga)*ayog
    fmpz_fmma(m2.data, ax, axog.data, ay, ayog.data);
    fmpz_divexact(m2.data, m2.data, ga.data);
//  m3 = byog*ub - bxog*vb
    fmpz_fmms(m3.data, byog.data, ub.data, bxog.data, vb.data);
//  m4 = _divexact(b.x,gb)*bxog + _divexact(b.y,gb)*byog
    fmpz_fmma(m4.data, bx, bxog.data, by, byog.data);
    fmpz_divexact(m4.data, m4.data, gb.data);
//  (m1, m3) = (m1*u + m3*v, _divexact(ga,g)*m3 - _divexact(gb,g)*m1)
    fmpz_fmma(t.data, m3.data, ga.data, m1.data, gb.data);
    fmpz_fmma(m1.data, m1.data, u.data, m3.data, v.data);
    fmpz_divexact(m3.data, t.data, g.data);
//  A = gcd(m2, m3, m4)
    fmpz_gcd(A.data, m2.data, m3.data);
    fmpz_gcd(A.data, A.data, m4.data);
    fmpz_fdiv_qr(t.data, B.data, m1.data, A.data);
    fmpz_one(C.data);
//  v, _ = _shortest_l∞(fmpz(1), mod(m1, A), A)
    _shortest_l_infinity(gx, gy, u.data, v.data, C.data, B.data, A.data);
//  z = isone(g) ? fmpzi(v[1], v[2]) : fmpzi(v[1]*g, v[2]*g)
    fmpz_mul(gx, gx, g.data);
    fmpz_mul(gy, gy, g.data);
//  return unit_canonicalize!(z)
    fmpzi_unit_canonicalize(gx, gy);
}
