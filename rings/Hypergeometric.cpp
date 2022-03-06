#include "globalstate.h"
#include "arithmetic.h"
#include "ex_print.h"
#include "code.h"
#include "timing.h"
#include "flintarb_wrappers.h"
#include <flint/fmpq.h>
#include <flint/arith.h>
#include "rfmpq_t.h"
#include "rfmpzi_t.h"
#include "rfmpqi_t.h"
#include "racb_mat_t.h"
#include "racb_poly_t.h"
#include "to_from_ex.h"

#define T_elem_t typename T::elem_t
#define T_int_ring_t typename T::int_ring_t
#define T_int_ring_elem_t typename T::int_ring_t::elem_t

#if 0
/** hyp_maj *******************************************************************/
// α > 0, β >= 0, i, j > 0, k >= 0)
//    [1] z^i/(1 - α*z^j)^k, or
//    [2] log((1 + β*z^j)/(1 - α*z^j))
template <class T>
class hyp_maj {
public:
    // coeff, α, i, j, k
    std::vector<std::tuple<T, xarb_t, slong, slong, slong>> rats;
    // coeff, α, β, j
    std::vector<std::tuple<T, xarb_t, xarb_t, slong>> logs;

    get_series(arb_poly_t res, slong prec);
};

template <class T>
std::ostream& operator<<(std::ostream& o, hyp_maj<T>& a)
{
    first = true;
    for (auto& r in a.rats)
    {
        if (first)
            o << " + ";
        first = false;
        o << std::get<0>(r) << "*z^" << std::get<2>(r) << "/(1 - ";
        o << std::get<1>(r) << "*z^" << std::get<3>(r) << ")^" << std::get<4>(r);
    }

    for (auto& r in a.logs)
    {
        if (first)
            o << " + ";
        first = false;
        o << std::get<0>(r) << "*Log[1 + " << std::get<2>(r)  << "*z^" << std::get<3>(r) << ")"
        o << "/(1 - " << std::get<1>(r)  << "*z^" << std::get<3>(r) << ")]";
    }

    if (first)
        o << "0";

    return o;
}

template <>
hyp_maj<xarb_t>::get_series(arb_poly_t res, slong ord, slong prec)
{
    arb_poly_fit_length(res, ord);
    _arb_vec_zero(res->coeffs, ord);

    xarb_t t, t2;

    // z^i/(1-a*z^j)^k = Sum[Binomial[n+k-1,n]*a^n*z^(j*n+i),{n,0,Infinity}]
    xfmpz_t bin;
    for (auto& r : rats)
    {
        slong i = std::get<2>(r);
        slong j = std::get<3>(r);
        slong k = std::get<4>(r);
        fmpz_one(bin.data);
        for (slong n = 0; j*n + i < ord; n++)
        {
            arb_pow_ui(t.data, std::get<1>(r).data, n, prec);
            arb_mul_fmpz(t.data, t.data, bin.data);
            arb_addmul(res->coeffs + j*n + i, std::get<0>(r).data, t.data, prec);
            fmpz_mul_ui(bin.data, bin.data, n + k);
            fmpz_divexact_ui(bin.data, bin.data, n + 1);
        }
    }

    // Log[(1 + b*z^j)/(1 - a*z^j)] = Sum[(a^n-(-b)^n)/n*z^(j*n), {n,1,Infinity}]
    for (auto& r : rats)
    {
        slong j = std::get<3>(r);
        for (slong n = 1; j*n < ord; n++)
        {
            arb_pow_ui(t.data, std::get<1>(r).data, n, prec);
            arb_pow_ui(t2.data, std::get<2>(r).data, n, prec);
            if (n & 1)
                arb_add(t.data, t.data, t2.data);
            else
                arb_sub(t.data, t.data, t2.data);
            arb_div_ui(t.data, t.data, n);
            arb_addmul(res->coeffs + j*n, std::get<0>(r).data, t.data, prec);
        }
    }

    _arb_poly_set_length(res, ord);
    _arb_poly_normalise(res);
}
#endif

/** bs_product ****************************************************************/

template <class T>
class bs_product {
public:
    slong top;
    T_elem_t tmp;
    std::vector<std::tuple<T_elem_t, slong, T_elem_t>> stack;

    bs_product(T& R) : top(-1), tmp(R.new_elem()) {}

    const T_elem_t& product()
    {
        return std::get<2>(stack.at(top));
    }

    void one(T& R)
    {
        top = -1;
    }

    T_elem_t& push(T& R, T_elem_t x);
};

template <class T>
T_elem_t& bs_product<T>::push(
    T& R,
    T_elem_t x)
{
    slong n = top;
    auto & M = stack;
    n++;
    assert(n <= M.size());
    if (n >= M.size())
        M.emplace_back(R.new_elem(), 1, R.new_elem());
    std::get<1>(M[n]) = 1;
    R.swap(std::get<0>(M[n]), x);

    while (n > 0 && std::get<1>(M[n]) >= std::get<1>(M[n-1]))
    {
        R.mul(tmp, std::get<0>(M[n]), std::get<0>(M[n-1]));
        R.swap(std::get<0>(M[n-1]), tmp);
        std::get<1>(M[n-1]) += std::get<1>(M[n]);
        n--;
    }

    if (n > 0)
        R.mul(std::get<2>(M[n]), std::get<0>(M[n]), std::get<2>(M[n-1]));
    else
        R.set(std::get<2>(M[n]), std::get<0>(M[n]));

    top = n;
    return std::get<2>(M[n]);
}

/** diffequ *******************************************************************/

template <class E>
class upoly {
public:
    std::vector<E> coeffs;
};

template <class E>
class bpoly {
public:
    std::vector<upoly<E>> coeffs;
};

template <class E>
class BiBPoly {
public:
    bpoly<E> xinside;
    bpoly<E> xoutside;
};

template <class E>
class LEquation {
public:
    BiBPoly<E> data;
};

template <class E>
class REquation {
public:
    BiBPoly<E> data;
};

template <class E>
class DEquation {
public:
    BiBPoly<E> data;
};

template <class E>
ex to_ex(const upoly<E>& a, er x)
{
    size_t n = a.coeffs.size();
    uex e; e.init_push_backr(gs.sym_sPlus.get(), n); 
    for (size_t i = n; i > 0; i--)
    {
        ex c = to_ex(a.coeffs[i - 1]);
        ex p = ex_powx(ecopy(x), emake_int_ui(i - 1));
        e.push_back(ex_mulx(c, p));
    }
    return ex_canonicalize_plus(e.release());
}

template <class E>
ex to_ex(const bpoly<E>& a, er x, er y)
{
    size_t n = a.coeffs.size();
    uex e; e.init_push_backr(gs.sym_sPlus.get(), n); 
    for (size_t i = n; i > 0; i--)
    {
        ex c = to_ex(a.coeffs[i - 1], x);
        ex p = ex_powx(ecopy(y), emake_int_ui(i - 1));
        e.push_back(emake_node(gs.sym_sNonCommutativeMultiply.copy(), c, p));
    }
    return ex_canonicalize_plus(e.release());
}

template <class E>
std::ostream& operator<<(std::ostream& o, const upoly<E>& a)
{
    wex x(emake_str("x"));
    wex e(to_ex(a, x.get()));
    o << e;
    return o;
}

template <class E>
std::ostream& operator<<(std::ostream& o, const bpoly<E>& a)
{
    wex x(emake_str("x"));
    wex y(emake_str("y"));
    wex e(to_ex(a, x.get(), y.get()));
    o << e;
    return o;
}

template <class E>
std::ostream& operator<<(std::ostream& o, const BiBPoly<E>& a)
{
    o << "{" << a.xinside << ", " << a.xoutside << "}";
    return o;
}

template <class E>
std::ostream& operator<<(std::ostream& o, LEquation<E>& a)
{
    o << a.data;
    return o;
}

template <class E>
std::ostream& operator<<(std::ostream& o, REquation<E>& a)
{
    o << a.data;
    return o;
}

template <class E>
std::ostream& operator<<(std::ostream& o, DEquation<E>& a)
{
    o << a.data;
    return o;
}

template <class T>
bool upoly_is_zero(
    T& R,
    upoly<T_elem_t>& a)
{
    return a.coeffs.empty();
}

template <class T>
bool bpoly_is_zero(
    T& R,
    bpoly<T_elem_t>& a)
{
    return a.coeffs.empty();
}

template <class E>
slong upoly_degree(upoly<E>& a)
{
    return (slong)a.coeffs.size() - 1;
}

template <class E>
slong bpoly_degree_outside(bpoly<E>& a)
{
    return (slong)a.coeffs.size() - 1;
}

template <class T>
void normalize(
    T& R,
    upoly<T_elem_t>& a)
{
    size_t n = a.coeffs.size();
    if (n > 0 && R.is_zero(a.coeffs[n - 1]))
    {
        do {
            n--;
        } while (n > 0 && R.is_zero(a.coeffs[n - 1]));
        a.coeffs.resize(n);
    }
}

template <class T>
void normalize(
    T& R,
    bpoly<T_elem_t>& a)
{
    size_t n = a.coeffs.size();
    if (n > 0 && upoly_is_zero<T>(R, a.coeffs[n - 1]))
    {
        do {
            n--;
        } while (n > 0 && upoly_is_zero<T>(R, a.coeffs[n - 1]));
        a.coeffs.resize(n);
    }
}

template <class T>
void upoly_coeff(
    T& R,
    T_elem_t& c,
    upoly<T_elem_t>& a,
    slong i)
{
    if (i < 0 || i >= a.coeffs.size())
        R.zero(c);
    else
        R.set(c, a.coeffs[i]);
}

template <class T>
void upoly_add_term(
    T& R,
    upoly<T_elem_t>& x,
    size_t i,
    T_elem_t& c)
{
    if (i >= x.coeffs.size())
        x.coeffs.resize(i + 1);
    R.add(x.coeffs[i], c);
    normalize<T>(R, x);
}

template <class T>
void bpoly_add_term(
    T& R,
    bpoly<T_elem_t>& x,
    size_t i, // outside index
    size_t j, // inside index
    T_elem_t& c)
{
    if (i >= x.coeffs.size())
        x.coeffs.resize(i + 1);
    upoly_add_term<T>(R, x.coeffs[i], j, c);
    normalize<T>(R, x);
}

template <class T>
void bpoly_swap_vars(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& y)
{
    x.coeffs.clear();
    for (size_t i = y.coeffs.size(); i > 0; i--)
    {
        auto & yi = y.coeffs[i - 1];
        for (size_t j = yi.coeffs.size(); j > 0; j--)
            bpoly_add_term(R, x, j - 1, i - 1, yi.coeffs[j - 1]);
    }
}

template <class T>
void complete_xinside(
    T& R,
    BiBPoly<T_elem_t>& a)
{
    if (a.xinside.coeffs.empty())
        bpoly_swap_vars<T>(R, a.xinside, a.xoutside);
}

template <class T>
void complete_xoutside(
    T& R,
    BiBPoly<T_elem_t>& a)
{
    if (a.xoutside.coeffs.empty())
        bpoly_swap_vars<T>(R, a.xoutside, a.xinside);
}

template <class T>
void complete(
    T& R,
    BiBPoly<T_elem_t>& a)
{
    if (a.xinside.coeffs.empty())
        bpoly_swap_vars<T>(R, a.xinside, a.xoutside);
    else if (a.xoutside.coeffs.empty())
        bpoly_swap_vars<T>(R, a.xoutside, a.xinside);
}

template <class T>
bpoly<T_elem_t>& coeffs_xinside(
    T& R,
    LEquation<T_elem_t>& a)
{
    complete_xinside<T>(R, a.data);
    return a.data.xinside.coeffs;
}

template <class T>
bpoly<T_elem_t>& coeffs_xoutside(
    T& R,
    LEquation<T_elem_t>& a)
{
    complete_xoutside<T>(R, a.data);
    return a.data.xoutside;
}

template <class T>
bool upoly_is_equal(
    T& R,
    upoly<T_elem_t>& a,
    upoly<T_elem_t>& b)
{
    size_t n = a.coeffs.size();
    if (n != b.coeffs.size())
        return false;
    for (size_t i = 0; i < n; i++)
        if (!R.is_equal(a.coeffs[i], b.coeffs[i]))
            return false;
    return true;
}

template <class T>
bool bpoly_is_equal(
    T& R,
    bpoly<T_elem_t>& a,
    bpoly<T_elem_t>& b)
{
    size_t n = a.coeffs.size();
    if (n != b.coeffs.size())
        return false;
    for (size_t i = 0; i < n; i++)
        if (!upoly_is_equal<T>(R, a.coeffs[i], b.coeffs[i]))
            return false;
    return true;
}

template <class T>
void upoly_zero(
    T& R,
    upoly<T_elem_t>& x)
{
    x.coeffs.clear();
}

template <class T>
void upoly_one(
    T& R,
    upoly<T_elem_t>& x)
{
    x.coeffs.resize(1);
    R.one(x.coeffs.at(0));
}

template <class T>
void upoly_set(
    T& R,
    upoly<T_elem_t>& x,
    T_elem_t& a)
{
    x.coeffs.resize(1);
    R.set(x.coeffs.at(0), a);
    normalize(R, x);
}

template <class T>
void upoly_set(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        R.set(x.coeffs[i], a.coeffs[i]);
}

template <class T>
void upoly_add(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    upoly<T_elem_t>& b)
{
    assert(&x != &a && &x != &b);
    size_t an = a.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)            
            R.set(x.coeffs[i], b.coeffs[i]);
        else if (i >= bn)
            R.set(x.coeffs[i], a.coeffs[i]);
        else
            R.add(x.coeffs[i], a.coeffs[i], b.coeffs[i]);
    }
    normalize<T>(R, x);
}

template <class T>
void upoly_add(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& b)
{
    size_t an = x.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)            
            R.set(x.coeffs[i], b.coeffs[i]);
        else if (i < bn)
            R.add(x.coeffs[i], b.coeffs[i]);
    }
    normalize<T>(R, x);
}

template <class T>
void bpoly_add(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& b)
{
    size_t an = x.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)            
            upoly_set<T>(R, x.coeffs[i], b.coeffs[i]);
        else if (i < bn)
            upoly_add<T>(R, x.coeffs[i], b.coeffs[i]);
    }
    normalize<T>(R, x);
}

template <class T>
void upoly_madd_scalar(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    upoly<T_elem_t>& b,
    T_elem_t& c)
{
    assert(&x != &a && &x != &b);
    size_t an = a.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)            
            R.mul(x.coeffs[i], b.coeffs[i], c);
        else if (i >= bn)
            R.set(x.coeffs[i], a.coeffs[i]);
        else
            R.madd(x.coeffs[i], a.coeffs[i], b.coeffs[i], c);
    }
    normalize<T>(R, x);
}

template <class T>
void upoly_sub(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    upoly<T_elem_t>& b)
{
    assert(&x != &a && &x != &b);
    size_t an = a.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)            
            R.neg(x.coeffs[i], b.coeffs[i]);
        else if (i >= bn)
            R.set(x.coeffs[i], a.coeffs[i]);
        else
            R.sub(x.coeffs[i], a.coeffs[i], b.coeffs[i]);
    }
    normalize<T>(R, x);
}

template <class T>
void upoly_sub(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& b)
{
    assert(&x != &b);
    size_t an = x.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)            
            R.neg(x.coeffs[i], b.coeffs[i]);
        else if (i < bn)
            R.sub(x.coeffs[i], b.coeffs[i]);
    }
    normalize<T>(R, x);
}

template <class T>
void upoly_neg(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a)
{
    size_t n = a.coeffs.size();
    x.coeffs.resize(n);
    for (size_t i = 0; i < n; i++)
        R.neg(x.coeffs[i], a.coeffs[i]);
}

template <class T>
void upoly_neg(
    T& R,
    upoly<T_elem_t>& x)
{
    for (auto& xi : x.coeffs)
        R.neg(xi);
}


template <class T>
void bpoly_add(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& a,
    bpoly<T_elem_t>& b)
{
    assert(&x != &a && &x != &b);
    size_t an = a.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)
            upoly_set<T>(R, x.coeffs[i], b.coeffs[i]);
        else if (i >= bn)
            upoly_set<T>(R, x.coeffs[i], a.coeffs[i]);
        else
            upoly_add<T>(R, x.coeffs[i], a.coeffs[i], b.coeffs[i]);
    }
    normalize<T>(R, x);
}

template <class T>
void bpoly_sub(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& a,
    bpoly<T_elem_t>& b)
{
    assert(&x != &a && &x != &b);
    size_t an = a.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = std::max(an, bn);
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        if (i >= an)
            upoly_neg<T>(R, x.coeffs[i], b.coeffs[i]);
        else if (i >= bn)
            upoly_set<T>(R, x.coeffs[i], a.coeffs[i]);
        else
            upoly_sub<T>(R, x.coeffs[i], a.coeffs[i], b.coeffs[i]);
    }
    normalize<T>(R, x);
}

template <class T>
void upoly_evaluate(
    T& R,
    T_elem_t& z,
    upoly<T_elem_t>& a,
    T_elem_t& b)
{
    R.zero(z);
    for (size_t i = a.coeffs.size(); i > 0; i--)
    {
        R.mul(z, b);
        R.add(z, a.coeffs[i - 1]);
    }
}

template <class T>
void upoly_evaluate(
    T& R,
    T_elem_t& z,
    upoly<T_elem_t>& a,
    slong b)
{
    size_t i = a.coeffs.size();
    if (i > 1)
    {
        R.mul(z, a.coeffs[i - 1], b);
        R.add(z, a.coeffs[i - 2]);
        for (i -= 2; i > 0; i--)
        {
            R.mul(z, b);
            R.add(z, a.coeffs[i - 1]);
        }
    }
    else if (i == 1)
    {
        R.set(z, a.coeffs[0]);
    }
    else
    {
        R.zero(z);
    }
}

template <class T>
void upoly_mul_linear(
    T& R,
    upoly<T_elem_t>& a,
    T_elem_t& b)
{
    T_elem_t t, s;
    R.zero(s);
    for (auto & ai : a.coeffs)
    {
        R.mul(t, ai, b);
        R.add(s, t);
        R.swap(ai, s);
    }
    if (!R.is_zero(s))
        a.coeffs.push_back(std::move(s));
}

template <class T>
void upoly_mul_linear(
    T& R,
    upoly<T_elem_t>& a,
    T_elem_t& b,
    T_elem_t& c)
{
    T_elem_t s;
    R.zero(s);
    for (auto & ai : a.coeffs)
    {
        R.mul(s, s, c);
        R.madd(s, ai, b);
        R.swap(ai, s);
    }
    R.mul(s, s, c);
    if (!R.is_zero(s))
        a.coeffs.push_back(std::move(s));
}

template <class T>
void upoly_mul(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    upoly<T_elem_t>& b)
{
    assert(&x != &a && &x != &b);
    size_t an = a.coeffs.size();
    size_t bn = b.coeffs.size();
    size_t xn = an + bn - 1;
    if (an < 1 || bn < 1)
    {
        x.coeffs.clear();
        return;
    }
    x.coeffs.resize(xn);
    for (size_t i = 0; i < xn; i++)
    {
        // 0 <= j <= bn - 1
        // 0 <= i - j <= an - 1
        // i >= j >= -an+1 + i
        size_t jstop = std::min(i, bn - 1);
        size_t j = i + 1 >= an ? i + 1 - an : 0;
        R.mul(x.coeffs.at(i), a.coeffs.at(i - j), b.coeffs.at(j));
        for (j++ ; j <= jstop; j++)
            R.madd(x.coeffs.at(i), a.coeffs.at(i - j), b.coeffs.at(j));
    }
    normalize<T>(R, x);
}

template <class T>
void upoly_divexact(
    T& R,
    upoly<T_elem_t>& q,
    upoly<T_elem_t>& A,
    upoly<T_elem_t>& b)
{
    assert(&q != &A && &q != &b);
    upoly<T_elem_t> a;
    upoly_set<T>(R, a, A);
    slong n = a.coeffs.size() - 1;
    slong m = b.coeffs.size() - 1;
    assert(n >= m && m > 0);
    q.coeffs.resize(n - m + 1);
    slong i = 0;
    slong j = n - m;
    slong l = 0;
    while (l < m && R.is_zero(b.coeffs.at(l)))
        l++;
    m -= l;
    n -= l;
    while (i <= j)
    {
        R.divexact(q.coeffs.at(i), a.coeffs.at(l + i), b.coeffs.at(l));
        for (slong k = 0; k <= m; k++)
            R.msub(a.coeffs.at(l + i + k), q.coeffs[i], b.coeffs[l + k]);
        if (++i > j)
            break;
        R.divexact(q.coeffs.at(j), a.coeffs.at(l + j + m), b.coeffs.at(l + m));
        for (slong k = 0; k <= m; k++)
            R.msub(a.coeffs.at(l + j + k), q.coeffs.at(j), b.coeffs.at(l + k));
        if (--j < i)
            break;
    }

    upoly_mul<T>(R, a, q, b);
    assert(upoly_is_equal<T>(R, a, A));
}

// R[x][y] * R[x]
template <class T>
void bpoly_mul_inside(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& a,
    upoly<T_elem_t>& b)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        upoly_mul<T>(R, x.coeffs[i], a.coeffs[i], b);
}

// R[x][y] / R[x]
template <class T>
void bpoly_divexact_inside(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& a,
    upoly<T_elem_t>& b)
{
    assert(&x != &a);
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        upoly_divexact<T>(R, x.coeffs[i], a.coeffs[i], b);
}

// aliasing OK
template <class T>
void upoly_derivative(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 1; i < an; i++)
        R.mul(x.coeffs.at(i - 1), a.coeffs.at(i), i);
    x.coeffs.resize(an - 1);
}

template <class T>
void upoly_mul_scalar(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    T_elem_t& b)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        R.mul(x.coeffs[i], a.coeffs[i], b);
    normalize<T>(R, x);
}

template <class T>
void upoly_mul_scalar(
    T& R,
    upoly<T_elem_t>& x,
    T_elem_t& a)
{
    for (auto& xi : x.coeffs)
        R.mul(xi, a);
    normalize<T>(R, x);
}

template <class T>
void upoly_mul_scalar(
    T& R,
    upoly<T_elem_t>& x,
    slong a)
{
    for (auto& xi : x.coeffs)
        R.mul(xi, a);
    normalize<T>(R, x);
}

template <class T>
void bpoly_mul_scalar(
    T& R,
    bpoly<T_elem_t>& x,
    T_elem_t& a)
{
    for (auto& xi : x.coeffs)
        bpoly_mul_scalar<T>(R, xi, a);
    normalize<T>(R, x);
}

template <class T>
void bpoly_mul_scalar(
    T& R,
    bpoly<T_elem_t>& x,
    slong a)
{
    for (auto& xi : x.coeffs)
        upoly_mul_scalar<T>(R, xi, a);
    normalize<T>(R, x);
}

template <class T>
void upoly_divexact_scalar(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    T_elem_t& b)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        R.divexact(x.coeffs[i], a.coeffs[i], b);
    normalize<T>(R, x);
}

template <class T>
void upoly_divexact_scalar(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    slong b)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        R.divexact(x.coeffs[i], a.coeffs[i], b);
    normalize<T>(R, x);
}

template <class T>
void bpoly_mul_scalar(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& a,
    T_elem_t& b)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        upoly_mul_scalar<T>(R, x.coeffs[i], a.coeffs[i], b);
    normalize<T>(R, x);
}


template <class T>
void upoly_shift(
    T& R,
    upoly<T_elem_t>& z,
    slong n)
{
    if (n > 0)
    {
        size_t m = z.coeffs.size();
        z.coeffs.resize(m + n);
        for (size_t i = 0; i < n; i++)
            R.zero(z.coeffs.at(m + i));
        for (size_t i = m; i > 0; i--)
            R.swap(z.coeffs.at(i - 1), z.coeffs.at(i - 1 + n));
    }
    else if (n < 0)
    {
        n = -n;
        size_t m = z.coeffs.size();
        if (n >= m)
        {
            z.coeffs.clear();
            return;
        }

        for (size_t i = 0; i < m - n; i++)
            R.swap(z.coeffs.at(i), z.coeffs.at(i + n));

        z.coeffs.resize(m - n);
    }
}

template<class T>
void bpoly_shift_outside(
    T& R,
    bpoly<T_elem_t>& z,
    slong n)
{
    if (n > 0)
    {
        size_t m = z.coeffs.size();
        z.coeffs.resize(m + n);
        for (size_t i = 0; i < n; i++)
            z.coeffs.at(m + i).coeffs.clear();
        for (size_t i = m; i > 0; i--)
            std::swap(z.coeffs.at(i - 1), z.coeffs.at(i - 1 + n));
    }
    else
    {
        n = -n;
        size_t m = z.coeffs.size();
        if (n >= m)
        {
            z.coeffs.clear();
            return;
        }

        for (size_t i = 0; i < m - n; i++)
            std::swap(z.coeffs[i], z.coeffs[i + n]);

        z.coeffs.resize(m - n);
    }
}

template<class T>
void bpoly_shift_inside(
    T& R,
    bpoly<T_elem_t>& z,
    slong n)
{
    for (auto & zi : z.coeffs)
        upoly_shift<T>(R, zi, n);
}

template<class T>
void shiftx(
    T& R,
    LEquation<T_elem_t>& a,
    slong n)
{
    bpoly_shift_inside(R, a.data.xinside, n);
    bpoly_shift_outside(R, a.data.xoutside, n);
}

template<class T>
void shiftx(
    T& R,
    REquation<T_elem_t>& a,
    slong n)
{
    bpoly_shift_inside(R, a.data.xinside, n);
    bpoly_shift_outside(R, a.data.xoutside, n);
}

template<class T>
void shiftx(
    T& R,
    DEquation<T_elem_t>& a,
    slong n)
{
    bpoly_shift_inside(R, a.data.xinside, n);
    bpoly_shift_outside(R, a.data.xoutside, n);
}

template<class T>
void upoly_taylor_shift(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    T_elem_t& c)
{
    slong n = a.coeffs.size();
    x.coeffs.resize(n);
    if (n < 1)
        return;

    xfmpz_t s;
    std::vector<xfmpz_t> f(n, xfmpz_t(1));
    for (slong i = 0; i < n; i++)
    {
        R.zero(x.coeffs.at(i));
        for (slong j = n - 1; j >= i; j--)
        {
            R.mul(x.coeffs.at(i), c);
            R.madd(x.coeffs.at(i), a.coeffs.at(j), f.at(j).data);
        }
        fmpz_zero(s.data);
        for (slong j = i; j < n; j++)
        {
            fmpz_swap(f.at(j).data, s.data);
            fmpz_add(s.data, s.data, f.at(j).data);
        }
    }
}

template<class T>
void upoly_taylor_shift(
    T& R,
    upoly<T_elem_t>& x,
    upoly<T_elem_t>& a,
    slong c)
{
    slong n = a.coeffs.size();
    x.coeffs.resize(n);
    if (n < 1)
        return;

    xfmpz_t s;
    std::vector<xfmpz_t> f(n, xfmpz_t(1));
    for (slong i = 0; i < n; i++)
    {
        R.zero(x.coeffs.at(i));
        for (slong j = n - 1; j >= i; j--)
        {
            R.mul(x.coeffs.at(i), c);
            R.madd(x.coeffs.at(i), a.coeffs.at(j), f.at(j).data);
        }
        fmpz_zero(s.data);
        for (slong j = i; j < n; j++)
        {
            fmpz_swap(f.at(j).data, s.data);
            fmpz_add(s.data, s.data, f.at(j).data);
        }
    }
}


template<class T>
void bpoly_taylor_shift_inside(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& a,
    T_elem_t& c)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        upoly_taylor_shift<T>(R, x.coeffs[i], a.coeffs[i], c);
}

template<class T>
void bpoly_taylor_shift_inside(
    T& R,
    bpoly<T_elem_t>& x,
    bpoly<T_elem_t>& a,
    slong c)
{
    size_t an = a.coeffs.size();
    x.coeffs.resize(an);
    for (size_t i = 0; i < an; i++)
        upoly_taylor_shift<T>(R, x.coeffs[i], a.coeffs[i], c);
}

/* L, R, D conversion */

template <class T>
void convert(
    T& R,
    REquation<T_elem_t>& z,
    LEquation<T_elem_t>& a)
{
    complete_xoutside<T>(R, a.data);
    auto & A = a.data.xoutside.coeffs;
    auto & Z = a.data.xoutside.coeffs;
    slong n = A.size();
    Z.resize(n);
    for (slong i = 0; i < n; i++)
        upoly_taylor_shift<T>(R, Z.at(i), A.at(i), i);
    z.data.xinside.coeffs.clear();
}

template <class T>
void convert(
    T& R,
    LEquation<T_elem_t>& z,
    REquation<T_elem_t>& a)
{
    complete_xoutside<T>(R, a.data);
    auto & A = a.data.xoutside.coeffs;
    auto & Z = z.data.xoutside.coeffs;
    slong n = A.size();
    Z.resize(n);
    for (slong i = 0; i < n; i++)
        upoly_taylor_shift<T>(R, Z.at(i), A.at(i), -i);
    z.data.xinside.coeffs.clear();
}

void stir2(fmpz_t z, slong n, slong i, slong m)
{
    assert(m == 0);
    arith_stirling_number_2(z, n, i);
}

template <class T>
void convert(
    T& R,
    DEquation<T_elem_t>& z,
    REquation<T_elem_t>& a)
{
    T_elem_t t;
    xfmpz_t f;
    z.data.xinside.coeffs.clear();
    z.data.xoutside.coeffs.clear();
    if (a.data.xinside.coeffs.empty())
    {
        auto & A = a.data.xoutside.coeffs;
        for (slong m = 0; m < A.size(); m++)
        {
            auto & Am = A[m].coeffs;
            for (slong n = 0; n < Am.size(); n++)
            for (slong i = 0; i <= n; i++)
            {
                stir2(f.data, n, i, 0);
                R.mul(t, Am[n], f.data);
                bpoly_add_term<T>(R, z.data.xinside, i, m + i, t);
            }
        }
    }
    else
    {
        auto & A = a.data.xinside.coeffs;
        for (slong n = 0; n < A.size(); n++)
        {
            auto & An = A[n].coeffs;
            for (slong m = 0; m < An.size(); m++)
            for (slong i = 0; i <= n; i++)
            {
                stir2(f.data, n, i, 0);
                R.mul(t, An[m], f.data);
                bpoly_add_term<T>(R, z.data.xinside, i, m + i, t);
            }
        }
    }
}

// input i is current degree
static void _mul_root(std::vector<xfmpz_t>& a, slong i, ulong r)
{
    fmpz_zero(a[i + 1].data);
    for (slong k = 0; k <= i; k++)
    {
        fmpz_submul_ui(a[i + 1].data, a[k].data, r);
        fmpz_swap(a[k].data, a[i + 1].data);
    }
}

template <class T>
void convert(
    T& R,
    REquation<T_elem_t>& z,
    DEquation<T_elem_t>& a)
{
    T_elem_t t;
    z.data.xinside.coeffs.clear();
    z.data.xoutside.coeffs.clear();
    complete_xinside<T>(R, a.data);
    auto & A = a.data.xinside.coeffs;
    std::vector<xfmpz_t> ppow(A.size() + 1);
    fmpz_one(ppow[0].data);
    for (slong i = 0; i < A.size(); i++)
    {
        auto & Ai = A[i].coeffs;
        for (slong j = 0; j < Ai.size(); j++)
        {
            if (j < i)
            {
                if (!R.is_zero(Ai[j]))
                    std::cout << "ignoring non-zero coefficient" << std::endl;
            }
            else
            {
                for (slong k = 0; k <= i; k++)
                {
                    R.mul(t, Ai[j], ppow[k].data);
                    bpoly_add_term<T>(R, z.data.xinside, k, j - i, t);
                }
            }
        }
        _mul_root(ppow, i, i);
    }
}

template <class T>
void convert(
    T& R,
    LEquation<T_elem_t>& z,
    DEquation<T_elem_t>& a)
{
    T_elem_t t;
    z.data.xinside.coeffs.clear();
    z.data.xoutside.coeffs.clear();
    complete_xinside<T>(R, a.data);
    auto & A = a.data.xinside.coeffs;
    std::vector<xfmpz_t> ppow(A.size() + 1);
    for (slong i = 0; i < A.size(); i++)
    {
        auto & Ai = A[i].coeffs;
        for (slong j = 0; j < Ai.size(); j++)
        {
            if (j < i)
            {
                if (!R.is_zero(Ai[j]))
                    std::cout << "ignoring non-zero coefficient" << std::endl;
            }
            else
            {
                fmpz_one(ppow[0].data);
                for (slong k = j - i; k < j; k++)
                    _mul_root(ppow, k - (j - i), k);

                for (slong k = 0; k <= i; k++)
                {
                    R.mul(t, Ai[j], ppow[k].data);
                    bpoly_add_term<T>(R, z.data.xinside, k, j - i, t);
                }
            }
        }
    }
}
    

/* transformations *********/
/*
function transform_rescale(a::Union{LEquation{T}, REquation{T}}, s::T, F) where T
  A = a.data.xoutside.coeffs
  zxoutside = upoly{T}[upoly{T}([A[1+i].coeffs[1+j]*s^i
                       for j in 0:length(A[1+i].coeffs)-1])
                       for i in 0:length(A)-1]
  A = a.data.xinside.coeffs
  zxinside = upoly{T}[upoly{T}([A[1+i].coeffs[1+j]*s^j
                       for j in 0:length(A[1+i].coeffs)-1])
                       for i in 0:length(A)-1]
  return typeof(a)(BiBPoly{T}(bpoly{T}(zxinside), bpoly{T}(zxoutside)))
end
*/
template <class T>
void transform_rescale(
    T& R,
    LEquation<T_elem_t>& z,
    T_elem_t& s)
{
    T_elem_t t;

    auto & A = z.data.xoutside.coeffs;
    for (slong i = 0; i < A.size(); i++)
    for (slong j = 0; j < A[i].coeffs.size(); j++)
    {
        R.pow(t, s, i);
        R.mul(A[i].coeffs[j], t);
    }

    auto & B = z.data.xinside.coeffs;
    for (slong i = 0; i < B.size(); i++)
    for (slong j = 0; j < B[i].coeffs.size(); j++)
    {
        R.pow(t, s, j);
        R.mul(B[i].coeffs[j], t);
    }
}


template <class T>
void transform_shift(
    T& R,
    DEquation<T_elem_t>& z,
    DEquation<T_elem_t>& a,
    T_elem_t& c)
{
    complete_xinside<T>(R, a.data);
    auto & A = a.data.xinside.coeffs;
    auto & Z = z.data.xinside.coeffs;
    slong n = A.size();
    Z.resize(n);
    for (slong i = 0; i < n; i++)
        upoly_taylor_shift<T>(R, Z.at(i), A.at(i), c);
    z.data.xoutside.coeffs.clear();
}

/* equation for pFq(a;b|x) */

template <class T>
void pfq_equ_at_zero(
    T& R,
    LEquation<T_elem_t>& z,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    z.data.xinside.coeffs.clear();
    bpoly<T_elem_t>& Z = z.data.xoutside;
    T_elem_t t;

    Z.coeffs.resize(2);

    Z.coeffs[1].coeffs.resize(1);
    R.set(Z.coeffs[1].coeffs[0], -1);
    for (auto & ai : a)
    {
        R.sub(t, ai, 1);
        upoly_mul_linear<T>(R, Z.coeffs[1], t);
    }

    Z.coeffs[0].coeffs.resize(1);
    R.set(Z.coeffs[0].coeffs[0], 1);
    for (auto & bi : b)
    {
        R.sub(t, bi, 1);
        upoly_mul_linear<T>(R, Z.coeffs[0], t);
    }
    upoly_shift<T>(R, Z.coeffs[0], 1);
}

template <class T>
void pfq_equ_at_zero(
    T& R,
    REquation<T_elem_t>& z,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    z.data.xinside.coeffs.clear();
    bpoly<T_elem_t>& Z = z.data.xoutside;
    T_elem_t t;

    Z.coeffs.resize(2);

    Z.coeffs[1].coeffs.resize(1);
    R.set(Z.coeffs[1].coeffs[0], 1);
    for (auto & ai : a)
        upoly_mul_linear<T>(R, Z.coeffs[1], ai);

    Z.coeffs[0].coeffs.resize(1);
    R.set(Z.coeffs[0].coeffs[0], -1);
    for (auto & bi : b)
    {
        R.sub(t, bi, 1);
        upoly_mul_linear<T>(R, Z.coeffs[0], t);
    }
    upoly_shift<T>(R, Z.coeffs[0], 1);
}


/* equation for pFq(a;b|-1/x) */

template <class T>
void pfq_equ_at_infinity(
    T& R,
    LEquation<T_elem_t>& z,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    z.data.xinside.coeffs.clear();
    bpoly<T_elem_t>& Z = z.data.xoutside;
    T_elem_t t;

    Z.coeffs.resize(2);

    Z.coeffs[1].coeffs.resize(2);
    R.one(Z.coeffs[1].coeffs[0]);
    R.one(Z.coeffs[1].coeffs[1]);
    R.neg(Z.coeffs[1].coeffs[(a.size()^b.size()^1)&1]);
    for (auto & bi : b)
    {
        R.neg(t, bi);
        upoly_mul_linear<T>(R, Z.coeffs[1], t);
    }

    Z.coeffs[0].coeffs.resize(1);
    R.set(Z.coeffs[0].coeffs[0], 1);
    for (auto & ai : a)
    {
        R.neg(t, ai);
        upoly_mul_linear<T>(R, Z.coeffs[0], t);
    }
}

template <class T>
void pfq_equ_at_infinity(
    T& R,
    REquation<T_elem_t>& z,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    z.data.xinside.coeffs.clear();
    bpoly<T_elem_t>& Z = z.data.xoutside;
    T_elem_t t;

    Z.coeffs.resize(2);

    Z.coeffs[0].coeffs.resize(1);
    R.set(Z.coeffs[0].coeffs[0], ((a.size()^b.size())&1) ? 1 : -1);
    for (auto & ai : a)
    {
        R.neg(t, ai);
        upoly_mul_linear<T>(R, Z.coeffs[0], t);
    }

    Z.coeffs[1].coeffs.resize(1);
    R.set(Z.coeffs[1].coeffs[0], -1);
    for (auto & bi : b)
    {
        R.sub(t, 1, bi);
        upoly_mul_linear<T>(R, Z.coeffs[1], t);
    }
    upoly_shift<T>(R, Z.coeffs[1], 1);
}

/* equation for pFq(a;b|1-x) */
template <class T>
void pfq_equ_at_one(
    T& R,
    LEquation<T_elem_t>& z,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    REquation<T_elem_t> eR1;
    DEquation<T_elem_t> eD1, eD2;
    T_elem_t t;

    pfq_equ_at_zero<T>(R, eR1, a, b);
//std::cout << "eR1: " << eR1 << std::endl;
    convert<T>(R, eD1, eR1);
//std::cout << "eD1: " << eD1 << std::endl;
    shiftx<T>(R, eD1, -1);
//std::cout << "eD1: " << eD1 << std::endl;
    R.one(t);
    transform_shift<T>(R, eD2, eD1, t);
//std::cout << "eD2: " << eD2 << std::endl;
    shiftx<T>(R, eD2, b.size());
//std::cout << "eD2: " << eD2 << std::endl;
    convert<T>(R, z, eD2);
//std::cout << "z: " << z << std::endl;
    R.neg(t);
//std::cout << "t: " << t << std::endl;
    transform_rescale<T>(R, z, t);
//std::cout << "z: " << z << std::endl;
}


/* equation for x^((d-1)/2-σ)*exp(-d/x)*pFq(a;b|x^-d) */

template <class T>
void _at_inf_exp_step(
    T& R,
    bpoly<T_elem_t>& z,
    bpoly<T_elem_t>& a,
    slong b,
    T_elem_t c)
{
    slong n = a.coeffs.size() - 1;
    std::vector<T_elem_t> null;
    T_elem_t t;

    z.coeffs.resize(n + 2);
    for (slong i = 0; i <= n + 1; i++)
    {
        std::vector<T_elem_t>& t1 = (i > n)  ? null : a.coeffs.at(i).coeffs;
        std::vector<T_elem_t>& t2 = (i == 0) ? null : a.coeffs.at(i - 1).coeffs;
        slong m = std::max(t1.size(), t2.size());

        z.coeffs[i].resize(m + 1);
        std::vector<T_elem_t>& zi = z.coeffs[i].coeffs;
        for (slong j = 0; j <= m; j++)
        {
            if (j < t1.size())
                R.mul(zi[j], t1[j], b);
            else
                R.zero();

            if (0 <= j - 1 && j - 1 < t1.size())
            {
                R.sub(t, c, j - 1);
                R.madd(zi[j], t1[j - 1], t);
            }

            if (0 <= j - 1 && j - 1 < t1.size())
                R.sub(zi[j], t2[j - 1]);
        }
    }
}

template <class T>
void pfq_equ_at_infinity_exp(
    T& R,
    REquation<T_elem_t>& z,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    T_elem_t t;
    bpoly<T_elem_t> e1, e2, et;
    slong j;

    z.data.xoutside.coeffs.clear();
    slong p = a.size();
    slong q = b.size();
    slong d = q + 1 - p;
    assert(d > 0);

    T_elem_t sigma;
    R.set(sigma, 1 - d);
    R.mul_2exp(t, -1);
    for (auto& bi : b)
        R.add(sigma, bi);
    for (auto& ai : a)
        R.sub(sigma, ai);

    e1.coeffs.resize(1);
    e1.coeffs.at(0).coeffs.resize(1);
    R.set(t, d);
    R.pow(e1.coeffs.at(0).coeffs.at(0), t, d);
    j = 0;
    for (auto& ai : a)
    {
        R.sub(t, j, sigma);
        R.madd(t, ai, d);
        _at_inf_exp_step<T>(R, et, e1, d, t);
        std::swap(e1, et);
        j++;
    }

    e2.coeffs.resize(2);
    e2.coeffs.at(0).coeffs.resize(2);
    R.set(e2.coeffs.at(0).coeffs.at(0), -d);
    R.set(e2.coeffs.at(0).coeffs.at(1), sigma);
    e2.coeffs.at(1).coeffs.resize(2);
    R.zero(e2.coeffs.at(0).coeffs.at(0));
    R.one(e2.coeffs.at(0).coeffs.at(1));
    j = 0;
    for (auto& bi : b)
    {
        R.sub(t, j - d, sigma);
        R.madd(t, bi, d);
        _at_inf_exp_step<T>(R, et, e2, d, t);
        std::swap(e2, et);
        j++;
    }

    bpoly_add<T>(R, z.data.xoutside, e1, e2);
    bpoly_shift_inside<T>(R, z.data.xoutside, -1);
}


/* equation for (1+x)^(-2ai)*pFq(a;b|4x/(1+x)^2) */
/*
# return (c(x)*θ+b(x))*a
function _l_equ_step(a::bpoly{T}, c::upoly{T}, b::upoly{T}, F) where T
  # z_i = b*a_i + c*(x*a_i' + a_{i-1})
  A = a.coeffs
  n = length(A)-1
  z = upoly{T}[i<=n ? mul(b, A[1+i], F) : upoly{T}() for i in 0:n+1]
  for i in 0:n+1
    if i<1
      t = shift!(derivative(A[1+i], F), 1, F)
    elseif i>n
      t = A[1+i-1]
    else
      t = add(shift!(derivative(A[1+i], F), 1, F), A[1+i-1], F)
    end
    z[1+i] = add(z[1+i], mul(c, t, F), F)
  end
  return bpoly{T}(z)
end
*/

template <class T>
void _l_equ_step(
    T& R,
    bpoly<T_elem_t>& z,
    bpoly<T_elem_t>& a,
    upoly<T_elem_t>& c,
    upoly<T_elem_t>& b)
{
    upoly<T_elem_t> t, t2;
    slong m = a.coeffs.size();// n = m-1
    z.coeffs.resize(m + 1);
    for (slong i = 0; i <= m; i++)
    {
        if (i <= m-1)
            upoly_mul<T>(R, z.coeffs.at(i), a.coeffs.at(i), b);
        else
            upoly_zero<T>(R, z.coeffs.at(i));

        if (i < 1)
        {
            upoly_derivative<T>(R, t, a.coeffs.at(i));
            upoly_shift<T>(R, t, 1);
        }
        else if (i > m-1)
        {
            upoly_set<T>(R, t, a.coeffs.at(i - 1));
        }
        else
        {
            upoly_derivative<T>(R, t, a.coeffs.at(i));
            upoly_shift<T>(R, t, 1);
            upoly_add<T>(R, t, a.coeffs.at(i - 1));
        }

        upoly_mul<T>(R, t2, t, c);
        upoly_add<T>(R, z.coeffs.at(i), t2);
    }
}

template <class T>
void pfq_equ_at_zero_quadratic(
    T& R,
    REquation<T_elem_t>& z,
    slong k,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    upoly<T_elem_t> C, B;
    bpoly<T_elem_t> e1, e2, et;
    slong j;

    z.data.xoutside.coeffs.clear();
    C.coeffs.resize(3);
    B.coeffs.resize(3);
    R.set(C.coeffs.at(0), 1);
    R.set(C.coeffs.at(1), 0);
    R.set(C.coeffs.at(2), -1);

    e1.coeffs.resize(2);
    e1.coeffs.at(0).coeffs.resize(1);
    R.set(e1.coeffs.at(0).coeffs.at(0), a[k]);
    e1.coeffs.at(1).coeffs.resize(1);
    R.one(e1.coeffs.at(1).coeffs.at(0));
    j = 1;
    for (slong i = 0; i < a.size(); i++)
    {
        if (i == k)
            continue;
        R.set(B.coeffs[0], a[i]);
        R.sub(B.coeffs[1], a[k], a[i]);
        R.mul(B.coeffs[1], 2);
        R.add(B.coeffs[1], j + 1);
        R.mul(B.coeffs[2], a[k], -2);
        R.add(B.coeffs[2], a[i]);
        R.add(B.coeffs[2], j - 1);
        _l_equ_step<T>(R, et, e1, C, B);
        std::swap(e1, et);
        j += 2;
    }

    e2.coeffs.resize(2);
    e2.coeffs.at(0).coeffs.resize(2);
    R.zero(e2.coeffs.at(0).coeffs.at(1));
    R.mul(e2.coeffs.at(0).coeffs.at(1), a[k], 2);
    e2.coeffs.at(1).coeffs.resize(2);
    R.one(e2.coeffs.at(1).coeffs.at(0));
    R.one(e2.coeffs.at(1).coeffs.at(1));
    j = 1;
    for (auto & bi : b)
    {
        R.sub(B.coeffs.at(0), bi, 1);
        R.sub(B.coeffs.at(1), a[k], bi);
        R.mul(B.coeffs.at(1), 2);
        R.add(B.coeffs.at(1), j + 2);
        R.mul(B.coeffs.at(2), a[k], -2);
        R.add(B.coeffs.at(2), bi);
        R.add(B.coeffs.at(2), j - 1);
        _l_equ_step<T>(R, et, e2, C, B);
        std::swap(e2, et);
        j += 2;
    }

    C.coeffs.resize(2);
    R.one(C.coeffs.at(0));
    R.one(C.coeffs.at(1));
    bpoly_mul_inside<T>(R, et, e2, C);
    bpoly_mul_scalar<T>(R, e1, -4);
    bpoly_shift_inside<T>(R, e1, 1);
    bpoly_add<T>(R, et, e1);

    if (b.empty())
    {
        B.coeffs.resize(2);
        R.set(B.coeffs.at(1), -1);
    }
    else
    {
        B.coeffs.resize(3);
        R.set(B.coeffs.at(2), 1);
        R.set(B.coeffs.at(1), -2);
    }
    R.one(B.coeffs.at(0));
    bpoly_divexact_inside<T>(R, z.data.xinside, et, B);
}


template <class T>
void pfq_equ_at_zero_quadratic(
    T& R,
    LEquation<T_elem_t>& z,
    slong k,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b)
{
    REquation<T_elem_t> t;
    pfq_equ_at_zero_quadratic<T>(R, t, k, a, b);
    convert<T>(R, z, t);
}


/** sum ***********************************************************************/

template <class T>
class BiElem {
public:
    xacb_t approx;
    T_elem_t exact;
    bool is_exact;
    bool have_approx;

    BiElem() : is_exact(false), have_approx(false) {}
    BiElem(const acb_t a) : approx(a), is_exact(false), have_approx(true) {}
    BiElem(const T_elem_t& a) : exact(a), is_exact(true), have_approx(false) {}

    void reset(const acb_t a)
    {
        acb_set(approx.data, a);
        is_exact = false;
        have_approx = true;
    }

    void reset(const T_elem_t& a)
    {
        exact = a;
        is_exact = true;
        have_approx = false;
    }

    xacb_t& get_approx(T& R, slong prec)
    {
        if (!have_approx)
        {
            have_approx = true;
            R.nset(approx.data, exact, prec);
        }
        return approx;
    }
};

template <class T>
class BiFrac {
public:
    xacb_t approx;
    T_elem_t num;
    T_elem_t den;
    bool is_exact;
    bool have_approx;

    BiFrac() : is_exact(false), have_approx(false) {}
    BiFrac(const xacb_t& a) : approx(a), is_exact(false), have_approx(true) {}
    BiFrac(const T_elem_t& n, const T_elem_t& d) : num(n), den(d), is_exact(true), have_approx(false) {}

    void reset(const acb_t a)
    {
        acb_set(approx.data, a);
        is_exact = false;
        have_approx = true;
    }

    void reset(const T_elem_t& a, const T_elem_t& b)
    {
        num = a;
        den = b;
        is_exact = true;
        have_approx = false;
    }

    xacb_t& get_approx(T& R, slong prec)
    {
        if (!have_approx)
        {
            xacb_t n, d;
            R.nset(n.data, num, prec + 1);
            R.nset(d.data, den, prec + 1);
            acb_div(approx.data, n.data, d.data, prec + 1);
            have_approx = true;
        }
        return approx;
    }
};

// matrices of polynomials mod Lambda^tau
template <class T>
class BiMatrix {
public:
    std::vector<xacb_t> approx;    
    std::vector<T_elem_t> num;
    T_elem_t den;
    bool is_exact;
    bool have_approx;

    BiMatrix() : is_exact(false), have_approx(false) {}

    BiMatrix(T& R, slong m, slong n, slong tau) :
        approx(m*n*tau), num(m*n*tau), is_exact(true), have_approx(false)
    {
        R.one(den);
    }

    void resize(slong n)
    {
        approx.resize(n);
        num.resize(n);
    }

    std::vector<xacb_t>& get_approx(T& R, slong prec);
    void force_approx(T& R, slong prec);
    void limit_precision(T& R, slong prec);
};

template <class T>
std::vector<xacb_t>& BiMatrix<T>::get_approx(
    T& R,
    slong prec)
{
    if (!have_approx)
    {
        xacb_t n, d;
        R.nset(d.data, den, prec + 1);

        slong zn = num.size();
        assert(zn == approx.size());
        for (slong i = 0; i < zn; i++)
        {
            R.nset(n.data, num[i], prec + 1);
            acb_div(approx[i].data, n.data, d.data, prec + 1);
        }
        have_approx = true;
    }
    return approx;
}

template <class T>
void BiMatrix<T>::force_approx(
    T& R,
    slong prec)
{
    if (!is_exact)
        return;

    is_exact = false;
    get_approx(R, prec);
}

template <class T>
void BiMatrix<T>::limit_precision(
    T& R,
    slong prec)
{
    if (!is_exact)
        return;

    if (R.nbits(den) > prec)
    {
        force_approx(R, prec);
        return;
    }

    for (const auto& zi : num)
    {
        if (R.nbits(zi) > prec)
        {
            force_approx(R, prec);
            return;
        }
    }
}


template <class T>
void mul_scalar(
    T& R,
    std::vector<T_elem_t>& x,
    T_elem_t& b)
{
    for (auto& xi : x)
        R.mul(xi, b);
}

template <class T>
void mul_scalar(
    T& R,
    std::vector<T_elem_t>& x,
    std::vector<T_elem_t>& a,
    T_elem_t& b)
{
    size_t n = a.size();
    x.resize(n);
    for (size_t i = 0; i < n; i++)
        R.mul(x[i], a[i], b);
}


template <class T>
void madd_scalar(
    T& R,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b,
    T_elem_t& c)
{
    size_t n = a.size();
    assert(b.size() >= n);
    for (size_t i = 0; i < n; i++)
        R.madd(a[i], b[i], c);
}

void mul_scalar(
    std::vector<xacb_t>& a,
    acb_t c,
    slong prec)
{
    for (auto& ai : a)
        acb_mul(ai.data, ai.data, c, prec);
}

void add(
    std::vector<xacb_t>& a,
    std::vector<xacb_t>& b,
    std::vector<xacb_t>& c,
    slong prec)
{
    size_t n = a.size();
    assert(n == b.size());
    assert(n == c.size());
    for (size_t i = 0; i < n; i++)
        acb_add(a[i].data, b[i].data, c[i].data, prec);
}


template <class T>
void mullow(
    T& R,
    BiMatrix<T>& z,  // m x n (x tau)
    BiMatrix<T>& a,  // m x p (x tau)
    BiMatrix<T>& b,  // p x n (x tau)
    slong prec,
    slong m, slong p, slong n, slong tau)
{
    assert(m*n*tau == z.num.size());
    assert(m*n*tau == z.approx.size());
    assert(m*p*tau == a.num.size());
    assert(m*p*tau == a.approx.size());
    assert(p*n*tau == b.num.size());
    assert(p*n*tau == b.approx.size());

#define zIDX(i, j, k) [((i)*n + (j))*tau + (k)]
#define aIDX(i, j, k) [((i)*p + (j))*tau + (k)]
#define bIDX(i, j, k) [((i)*n + (j))*tau + (k)]

//std::cout << "mullow called" << std::endl;
//std::cout << "a.is_exact: " << a.is_exact << std::endl;
//std::cout << "b.is_exact: " << b.is_exact << std::endl;

    if (a.is_exact && b.is_exact)
    {
        z.is_exact = true;
        z.have_approx = false;
        R.mul(z.den, a.den, b.den);

        for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
        for (slong k = 0; k < tau; k++)
            R.zero(z.num zIDX(i,j,k));
        
        for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
        for (slong h = 0; h < p; h++)
        for (slong k = 0; k < tau; k++)
        for (slong l = 0; l <= k; l++)
            R.madd(z.num zIDX(i,j,k), a.num aIDX(i,h,l), b.num bIDX(h,j,k-l));

        return;
    }

    z.is_exact = false;
    z.have_approx = true;
    for (slong i = 0; i < m; i++)
    for (slong j = 0; j < n; j++)
    for (slong k = 0; k < tau; k++)
        acb_zero(z.approx zIDX(i,j,k).data);

    if (!a.is_exact && !b.is_exact)
    {
        for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
        for (slong h = 0; h < p; h++)
        for (slong k = 0; k < tau; k++)
        for (slong l = 0; l <= k; l++)
            acb_addmul(z.approx zIDX(i,j,k).data, a.approx aIDX(i,h,l).data, b.approx bIDX(h,j,k-l).data, prec);
    }
    else
    {
        xacb_t t;

        assert(a.is_exact ^ b.is_exact);

        if (b.is_exact)
        {
            for (slong i = 0; i < m; i++)
            for (slong j = 0; j < n; j++)
            for (slong h = 0; h < p; h++)
            for (slong k = 0; k < tau; k++)
            for (slong l = 0; l <= k; l++)
                R.nmadd(z.approx zIDX(i,j,k).data, a.approx aIDX(i,h,l).data, b.num bIDX(h,j,k-l), prec);

            R.nset(t.data, b.den, prec);
        }
        else
        {
            for (slong i = 0; i < m; i++)
            for (slong j = 0; j < n; j++)
            for (slong h = 0; h < p; h++)
            for (slong k = 0; k < tau; k++)
            for (slong l = 0; l <= k; l++)
                R.nmadd(z.approx zIDX(i,j,k).data, b.approx bIDX(h,j,k-l).data, a.num aIDX(i,h,l), prec);

            R.nset(t.data, a.den, prec);
        }

        for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
        for (slong k = 0; k < tau; k++)
            acb_div(z.approx zIDX(i,j,k).data, z.approx zIDX(i,j,k).data, t.data, prec);
    }

#undef zIDX
#undef aIDX
#undef bIDX
}


template <class T>
void content(
    T& R,
    T_elem_t& g,
    upoly<T_elem_t>& a)
{
    for (auto& ai : a.coeffs)
        R.gcd(g, g, ai);
}

// x = a/b
template <class T>
void numerator(
    T_int_ring_t& ZR,
    upoly<T_int_ring_elem_t>& x,
    T& R,
    upoly<T_elem_t>& a,
    T_elem_t& b)
{
    T_elem_t t;
    T_int_ring_elem_t d;

    size_t n = a.coeffs.size();
    x.coeffs.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        R.divexact(t, a.coeffs[i], b);
        R.numerator_denominator(x.coeffs[i], d, t);
        assert(ZR.is_one(d));
    }
}


template <class T>
class SumBsCtx {
public:
    slong s, tau, delta;
    std::vector<upoly<T_elem_t>> eqnum;
    upoly<T_elem_t> eqden;
    BiFrac<T> z;
    T_elem_t lambda_num;
    T_elem_t lambda_den;
    BiMatrix<T> Sigma;
    BiMatrix<T> M;
    BiMatrix<T> Delta;
    BiMatrix<T> Sigma0;
    BiMatrix<T> M0;
    // temp space
    std::vector<std::pair<BiMatrix<T>, BiMatrix<T>>> tmatrix;
    std::vector<T_elem_t> tnumeval;
    T_elem_t tdeneval;
    // scratch space
    std::vector<T_elem_t> tff;
    T_elem_t t1, t2, t3, tzn, tzd;
    xacb_t ta, ta1, ta2;

    std::vector<xacb_t>& get_sum(T& R, slong prec)
    {
        return Sigma.get_approx(R, prec);
    }

    std::pair<BiMatrix<T>, BiMatrix<T>> pop_tmatrix(T& R)
    {
        if (tmatrix.empty())
        {
            return std::pair<BiMatrix<T>, BiMatrix<T>>(
                               BiMatrix<T>(R, delta, s, tau),
                               BiMatrix<T>(R, s, s, tau));
        }
        else
        {
            std::pair<BiMatrix<T>, BiMatrix<T>> res = std::move(tmatrix.back());
            tmatrix.pop_back();
            return res;
        }
    }

    void push_tmatrix(std::pair<BiMatrix<T>, BiMatrix<T>>&& a)
    {
        tmatrix.push_back(a);
    }
    
};

template <class T>
void sum_bs_der_recursive(
    T& R,
    SumBsCtx<T>& S,
    slong hi, slong lo,
    slong prec);

//  in: sum [lo, mid) is in (S.Σ0, S.M0)
// out: put sum [lo, hi) in (S.Σ, S.M)
template <class T>
void sum_bs_der_continue(
    T& R,
    SumBsCtx<T>& S,
    slong hi, slong mid, slong lo,
    slong prec)
{
    auto pt0 = S.pop_tmatrix(R);
    BiMatrix<T>& Sigma0 = std::get<0>(pt0); 
    BiMatrix<T>& M0 = std::get<1>(pt0);
    std::swap(S.Sigma0, Sigma0);
    std::swap(S.M0, M0);

    assert(0 <= lo);
    assert(lo < mid);
    assert(mid < hi);

    sum_bs_der_recursive<T>(R, S, hi, mid, prec);

    auto pt1 = S.pop_tmatrix(R);
    BiMatrix<T>& Sigma1 = std::get<0>(pt1); 
    BiMatrix<T>& M1 = std::get<1>(pt1);
    std::swap(S.Sigma, Sigma1);
    std::swap(S.M, M1);

//std::cout << "combining [" << lo << ", " << mid << ") with [" << mid << ", " << hi << ")"<< std::endl;
//std::cout << "Sigma1: " << Sigma1.num << "/" << Sigma1.den << std::endl;
//std::cout << "    M1: " << M1.num << "/" << M1.den << std::endl;
//std::cout << "Sigma0: " << Sigma0.num << "/" << Sigma0.den << std::endl;
//std::cout << "    M0: " << M0.num << "/" << M0.den << std::endl;

    mullow<T>(R, S.Delta, Sigma1, M0, prec, S.delta, S.s, S.s, S.tau);
    S.Delta.limit_precision(R, prec);
    if (S.z.is_exact && S.Delta.is_exact && Sigma0.is_exact)
    {
        R.pow(S.tzd, S.z.den, mid - lo);
        R.pow(S.tzn, S.z.num, mid - lo);
        S.Sigma.is_exact = true;
        S.Sigma.have_approx = false;
        R.mul(S.Sigma.den, S.tzd, Sigma0.den);
        R.mul(S.Sigma.den, Sigma1.den);
        R.divexact(S.t2, Sigma0.den, M0.den);
        R.mul(S.t3, S.t2, S.tzn);
        mul_scalar(R, S.Sigma.num, S.Delta.num, S.t3);
        R.mul(S.t3, Sigma1.den, S.tzd);
        madd_scalar(R, S.Sigma.num, Sigma0.num, S.t3);
        S.Sigma.limit_precision(R, prec);
        R.mul(S.Delta.den, S.tzd);
        mul_scalar(R, S.Delta.num, S.tzn);
    }
    else
    {
        S.Delta.force_approx(R, prec);
        acb_pow_si(S.ta.data, S.z.get_approx(R, prec).data, mid - lo, prec);
        mul_scalar(S.Delta.approx, S.ta.data, prec);
        S.Sigma.is_exact = false;
        S.Sigma.have_approx = true;
        add(S.Sigma.approx, Sigma0.get_approx(R, prec), S.Delta.approx, prec);
    }
    mullow<T>(R, S.M, M1, M0, prec, S.s, S.s, S.s, S.tau);
    S.M.limit_precision(R, prec);

//std::cout << "S.Sigma: " << S.Sigma.num << "/" << S.Sigma.den << std::endl;
//std::cout << "    S.M: " << S.M.num << "/" << S.M.den << std::endl;

    S.push_tmatrix(std::move(pt1));

    std::swap(S.Sigma0, Sigma0);
    std::swap(S.M0, M0);
    S.push_tmatrix(std::move(pt0));
}


// out: put sum from [lo, hi) in (S.Σ, S.M)
template <class T>
void sum_bs_der_recursive(
    T & R,
    SumBsCtx<T>& S,
    slong hi, slong lo,
    slong prec)
{
    assert(hi > lo);

//std::cout << "sum_bs_der_recursive called" << std::endl;
//std::cout << "hi: " << hi << std::endl;
//std::cout << "lo: " << lo << std::endl;


    assert(0 <= lo);
    assert(lo < hi);



    if (hi > lo + (S.z.is_exact ? 1 : 1))
    {
        slong mid = (hi + lo)/2;
        sum_bs_der_recursive<T>(R, S, mid, lo, prec);
        std::swap(S.Sigma0, S.Sigma);
        std::swap(S.M0, S.M);
        sum_bs_der_continue<T>(R, S, hi, mid, lo, prec);
        return;
    }

    slong delta = S.delta;
    slong s = S.s;
    slong tau = S.tau;

#define IDX(i0,i1,i2) [((i0)*(s) + (i1))*(tau) + (i2)]

//    S.Sigma.resize(delta*s*tau);
//    S.M.resize(s*s*tau);
//    S.Delta.resize(delta*s*tau);

//std::cout << "S.eqnum: " << S.eqnum << std::endl;
//std::cout << "S.eqden: " << S.eqden << std::endl;

    S.Sigma.is_exact = true;
    S.Sigma.have_approx = false;
    R.one(S.Sigma.den);
    for (slong d = 0; d < delta; d++)
    for (slong j = 0; j < s; j++)
    for (slong k = 0; k < tau; k++)
        R.zero(S.Sigma.num IDX(d,j,k));

    S.M.is_exact = true;
    S.M.have_approx = false;
    R.one(S.M.den);
    for (slong i = 0; i < s; i++)
    for (slong j = 0; j < s; j++)
    for (slong k = 0; k < tau; k++)
        R.zero(S.M.num IDX(i,j,k));

    upoly_evaluate<T>(R, S.M.den, S.eqden, lo);
    for (slong j = 0; j < s; j++)
    for (slong k = 0; k < tau; k++)
        upoly_evaluate<T>(R, S.M.num IDX(0,j,k), S.eqnum[j*tau+k], lo);

    for (slong i = 1; i < s; i++)
        R.set(S.M.num IDX(i,i-1,0), S.M.den);

    // ff = λden^(δ-1)*(Λ+λ+n-0)*(Λ+λ+n-1)*...*(Λ+λ+n-(d-1))
    slong fflen = 1;
    R.pow(S.tff[0], S.lambda_den, delta - 1);

    R.mul(S.Sigma.den, S.tff[0], S.M.den);
    for (slong j = 0; j < s; j++)
    for (slong k = 0; k < tau; k++)
        R.mul(S.Sigma.num IDX(0,j,k), S.tff[0], S.M.num IDX(0,j,k));

    for (slong d = 1; d < delta; d++)
    {
        R.mul(S.t1, S.lambda_den, lo - (d - 1));
        R.add(S.t1, S.lambda_num);

        // ff *= lo-(d-1)+λ + Λ
        for (slong i = 0; i < fflen; i++)
            R.divexact(S.tff[i], S.lambda_den);

        for (slong i = fflen - 1; i > 0; i--)
        {
            R.mul(S.tff[i], S.t1);
            R.madd(S.tff[i], S.tff[i - 1], S.lambda_den);
        }
        R.mul(S.tff[0], S.t1);

        if (fflen < tau)
        {
            fflen++;
            R.pow(S.tff[fflen - 1], S.lambda_den, delta - 1);
        }

        for (slong j = 0; j < s; j++)
        for (slong k = 0; k < tau; k++)
        for (slong l = 0; l < fflen && l <= k; l++)
            R.madd(S.Sigma.num IDX(d,j,k), S.tff[l], S.M.num IDX(0,j,k-l));
    }

    for (slong n = lo + 1; n < hi; n++)
    {
        assert(S.z.is_exact);

        // update M
        upoly_evaluate<T>(R, S.tdeneval, S.eqden, n);
        for (slong i = 0; i < s; i++)
        for (slong k = 0; k < tau; k++)
            upoly_evaluate<T>(R, S.tnumeval[i*tau+k], S.eqnum[i*tau+k], n);

        for (slong j = 0; j < s; j++)
        {
            // update column j
            for (slong k = 0; k < tau; k++)
                R.zero(S.tff[k]);

            for (slong k = 0; k < tau; k++)
            for (slong h = 0; h <= k; h++)
                R.madd(S.tff[k], S.tnumeval[0*tau+h], S.M.num IDX(0,j,k-h));

            for (slong i = s - 1; i >= 1; i--)
            {
                for (slong k = 0; k < tau; k++)
                for (slong h = 0; h <= k; h++)
                    R.madd(S.tff[k], S.tnumeval[i*tau+h], S.M.num IDX(i,j,k-h));

                for (slong k = 0; k < tau; k++)
                    R.mul(S.M.num IDX(i,j,k), S.tdeneval, S.M.num IDX(i-1,j,k));
            }

            for (slong k = 0; k < tau; k++)
                R.swap(S.M.num IDX(0,j,k), S.tff[k]);
        }
        R.mul(S.M.den, S.tdeneval);

        // update Σ
        R.mul(S.t1, S.tdeneval, S.z.den);
        R.mul(S.Sigma.den, S.t1);
        mul_scalar(R, S.Sigma.num, S.t1);
        // Σ.num += zn^(n-lo)*ff*M.num

        // ff = zn^(n-lo)*λd^(δ-1)*(Λ+λ+n-0)*(Λ+λ+n-1)*...*(Λ+λ+n-(d-1))
        fflen = 1;
        R.pow(S.tff[0], S.z.num, n - lo);
        R.pow(S.t1, S.lambda_den, delta - 1);
        R.mul(S.tff[0], S.t1);
        for (slong j = 0; j < s; j++)
        for (slong k = 0; k < tau; k++)
            R.madd(S.Sigma.num IDX(0,j,k), S.tff[0], S.M.num IDX(0,j,k));

        for (slong d = 1; d < delta; d++)
        {
            R.mul(S.t1, S.lambda_den, n - (d - 1));
            R.add(S.t1, S.lambda_num);
            // ff *= n-(d-1)+λ + Λ
            for (slong i = 0; i < fflen; i++)
                R.divexact(S.tff[i], S.lambda_den);

            for (slong i = fflen - 1; i > 0; i--)
            {
                R.mul(S.tff[i], S.t1);
                R.madd(S.tff[i], S.tff[i - 1], S.lambda_den);
            }
            R.mul(S.tff[0], S.t1);

            if (fflen < tau)
            {
                fflen++;
                R.pow(S.tff[fflen - 1], S.z.num, n-lo);
                R.pow(S.t1, S.lambda_den, delta - 1);
                R.mul(S.tff[fflen - 1], S.t1);
            }

            for (slong j = 0; j < s; j++)
            for (slong k = 0; k < tau; k++)
            for (slong l = 0; l < fflen && l <= k; l++)
                R.madd(S.Sigma.num IDX(d,j,k), S.tff[l], S.M.num IDX(0,j,k-l));
        }
    }

#undef IDX

//std::cout << "S.Sigma.num: " << S.Sigma.num << std::endl;
//std::cout << "S.Sigma.den: " << S.Sigma.den << std::endl;

//std::cout << "S.M.num: " << S.Sigma.num << std::endl;
//std::cout << "S.M.den: " << S.Sigma.den << std::endl;

    return;
}


template <class T>
void sum_bs_start(
    T& R, T_int_ring_t& ZR,
    SumBsCtx<T_int_ring_t>& S,
    std::vector<upoly<T_elem_t>>& num, slong s, slong tau,
    upoly<T_elem_t>& den,
    BiElem<T>& z,
    T_elem_t& lambda,
    slong delta,
    slong N1, slong N0,
    slong prec)
{
//std::cout << "sum_bs_start called" << std::endl;
//std::cout << "num: " << num << std::endl;
//std::cout << "den: " << den << std::endl;

    assert(0 <= N0);
    assert(N0 < N1);
    assert(num.size() == s*tau);

    T_elem_t g;

    R.zero(g);
    content<T>(R, g, den);
    for (slong i = 0; i < s*tau; i++)
        content<T>(R, g, num.at(i));

    numerator<T>(ZR, S.eqden, R, den, g);
    S.eqnum.resize(s*tau);
    for (slong i = 0; i < s*tau; i++)
        numerator<T>(ZR, S.eqnum.at(i), R, num.at(i), g);

    R.numerator_denominator(S.lambda_num, S.lambda_den, lambda);

//std::cout << "S.eqnum: " << S.eqnum << std::endl;
//std::cout << "S.eqden: " << S.eqden << std::endl;

    if (z.is_exact)
    {
        S.z.is_exact = true;
        S.z.have_approx = false;
        R.numerator_denominator(S.z.num, S.z.den, z.exact);
    }
    else
    {
        S.z.is_exact = false;
        S.z.have_approx = true;
        S.z.approx = z.approx;
    }

    S.s = s;
    S.tau = tau;
    S.delta = delta;
    S.Sigma.resize(delta*s*tau);
    S.M.resize(s*s*tau);
    S.Delta.resize(delta*s*tau);
    S.Sigma0.resize(delta*s*tau);
    S.M0.resize(s*s*tau);
    S.tmatrix.clear();
    S.tnumeval.resize(s*tau);
    S.tff.resize(tau);

    sum_bs_der_recursive(ZR, S, N1, N0, prec);
}


slong overlap_dist(
    const acb_t a,
    const acb_t Delta,
    const acb_t b,
    slong prec)
{
    if (acb_overlaps(a, b))
        return 0;

    xmag_t m, q;
    acb_get_mag(m.data, a);
    acb_get_mag(q.data, b);
    if (mag_cmp(m.data, q.data) > 0)
        mag_swap(m.data, q.data);
    acb_get_mag(q.data, Delta);
    mag_div(q.data, q.data, m.data);
    if (mag_is_special(q.data))
        return mag_is_zero(q.data) ? 0 : prec;
    fmpz_add_si(&q.data->exp, &q.data->exp, prec);
    if (fmpz_sgn(&q.data->exp) < 0)
        return 0;
    else if (fmpz_cmp_si(&q.data->exp, prec) >= 0)
        return prec;
    else
        return fmpz_get_si(&q.data->exp);
}

slong overlap_dist(
    std::vector<xacb_t>& a,
    std::vector<xacb_t>& Delta,
    std::vector<xacb_t>& b,
    slong length,
    slong prec)
{
    slong dist = 0;
    for (slong i = 0; i < length; i++)
        dist = std::max(dist, overlap_dist(a.at(i).data, Delta.at(i).data, b.at(i).data, prec));
    return dist;
}

template <class T>
slong sum_bs_continue(
    T& R,
    SumBsCtx<T>& S,
    slong hi, slong mid, slong lo,
    slong prec)
{
    std::swap(S.Sigma0, S.Sigma);
    std::swap(S.M0, S.M);
    sum_bs_der_continue<T>(R, S, hi, mid, lo, prec);
    return overlap_dist(S.Sigma.get_approx(R, prec),
                        S.Delta.get_approx(R, prec),
                        S.Sigma0.get_approx(R, prec),
                        S.delta*S.s*S.tau, prec);
}


template <class T>
void sum_bs(
    T& R,
    acb_t res,
    slong N,
    LEquation<T_elem_t>& E,
    std::vector<T_elem_t>& iv,
    T_elem_t& z,
    slong prec)
{
    auto& eq = coeffs_xoutside<T>(R, E);
    slong s = bpoly_degree_outside(eq);

    if (N > s)
    {
        T_int_ring_t ZR;
        SumBsCtx<T_int_ring_t> S;
        std::vector<upoly<T_elem_t>> numeq(s);
        BiElem<T> Z(z);
        T_elem_t lambda;

        for (slong i = 0; i < s; i++)
            upoly_neg<T>(R, numeq.at(i), eq.coeffs.at(1 + i));

        R.zero(lambda);
        sum_bs_start<T>(R, ZR, S, numeq, s, 1, eq.coeffs.at(0),
                        Z, lambda, 1, N, s, prec);

        std::vector<xacb_t>& Sigma = S.get_sum(ZR, prec);

        slong i = 1;
        R.nmul(res, Sigma[s - i].data, iv[i - 1], prec);
        for (i++; i <= s; i++)
        {
            xacb_t t;
            R.nmul(t.data, Sigma[s - i].data, iv[i - 1], prec);
            acb_add(res, res, t.data, prec + 10);
        }
    }
    else
    {
        acb_zero(res);
    }

    for (slong i = std::min(N, s); i > 0; i--)
    {
        R.nmul(res, res, z, prec + 10);
        R.nadd(res, res, iv.at(i - 1), prec + 10);
    }
}   


/* funny stuff for the evaluation of z^a*log(z)^b */

template <class T>
class BranchInfo {
public:
    T_elem_t lambda;
    BiElem<T> z;
    BiElem<T> invz;
    std::vector<xacb_t> logzpow;
    std::vector<xacb_t> zlogzpow;
    bool have_invz;

    BranchInfo(const T_elem_t& _lambda, const T_elem_t& _z) :
        lambda(_lambda), z(_z), have_invz(false)
    {
        logzpow.resize(1);
        acb_one(logzpow[0].data);
    }

    BranchInfo(const T_elem_t& _lambda, const acb_t _z) :
        lambda(_lambda), z(_z), have_invz(false)
    {
        logzpow.resize(1);
        acb_one(logzpow[0].data);
    }

    BranchInfo(const T_elem_t& _lambda, const T_elem_t& _z, const T_elem_t& _invz) :
        lambda(_lambda), z(_z), invz(_invz), have_invz(true)
    {
        logzpow.resize(1);
        acb_one(logzpow[0].data);
    }

    void reset(const T_elem_t& _lambda, const T_elem_t& _z)
    {
        lambda = _lambda;
        z.reset(_z);
        have_invz = false;
        logzpow.resize(1);
        acb_one(logzpow[0].data);
    }

    void reset_lambda(const T_elem_t& _lambda)
    {
        lambda = _lambda;
    }

    xacb_t& get_logz(T& R, slong prec);
    void zpow(T& R, acb_t res, slong i, slong prec);
    void zpowlogpow(T& R, acb_t res, slong i, ulong j, slong prec);
    xacb_t& get_logzpow(T& R, slong j, slong prec);
    xacb_t& get_zlogzpow(T& R, slong j, slong prec);
    void zpowlogpoly_eval(T& R, acb_t res, slong n, std::vector<xacb_t>& e, slong prec, acb_t t1);


};

template<class T>
xacb_t& BranchInfo<T>::get_logz(T& R, slong prec)
{
    if (logzpow.size() > 1)
        return logzpow[1];

    logzpow.resize(2);
    xacb_t& l = logzpow[1];
    if (have_invz)
    {
        acb_log(l.data, invz.get_approx(R, prec).data, prec);
        acb_neg(l.data, l.data);
    }
    else
    {
        acb_log(l.data, z.get_approx(R, prec).data, prec);
    }

    return l;
}

/* z^(λ+i) */
template<class T>
void BranchInfo<T>::zpow(
    T& R,
    acb_t res,
    slong i,
    slong prec)
{
    if (R.is_zero(lambda))
    {
        if (have_invz)
            acb_pow_si(res, invz.get_approx(R, prec).data, -i, prec);
        else
            acb_pow_si(res, z.get_approx(R, prec).data, i, prec);
        return;
    }

    xacb_t t;
    R.nset(t.data, lambda, prec + 60);
    acb_add_si(t.data, t.data, i, prec + 60);
    if (have_invz)
    {
        acb_neg(t.data, t.data);
        acb_pow(res, invz.get_approx(R, prec).data, t.data, prec);
    }
    else
    {
        acb_pow(res, z.get_approx(R, prec).data, t.data, prec);
    }
}

// Given |z| <= zmag and |arg z| <= zarg, return a ball for z^A*log(z)^j/j!
void fancy_pow_at_zero(
    acb_t res,
    arb_t zmag, arb_t zarg,
    acb_t a,
    ulong j,
    slong prec)
{
    assert(arb_is_nonnegative(zmag));
    assert(arb_is_nonnegative(zarg));
    xarb_t l, t, r, aoj, joa;
    if (arb_is_positive(acb_realref(a)))
    {
        if (arb_is_zero(zmag))
        {
            arb_zero(r.data);
        }
        else if (j == 0)
        {
            arb_pow(r.data, zmag, acb_realref(a), prec);
        }
        else
        {
            arb_div_ui(aoj.data, acb_realref(a), j, prec);
            arb_log(l.data, zmag, prec);
            arb_pow(t.data, zmag, aoj.data, prec);
            arb_mul(r.data, t.data, l.data, prec);
            arb_abs(r.data, r.data);
            // if !(-l > j/a)
            arb_ui_div(joa.data, j, acb_realref(a), prec);
            arb_neg(l.data, l.data);
            if (!arb_gt(l.data, joa.data))
            {
                arb_mul_2exp_si(joa.data, joa.data, -1);
                arb_max(r.data, r.data, joa.data, prec);
            }
            arb_addmul(r.data, t.data, zarg, prec);
            arb_pow_ui(r.data, r.data, j, prec);
            arb_fac_ui(t.data, j, prec);
            arb_div(r.data, r.data, t.data, prec);
        }
        arb_mul(t.data, acb_imagref(a), zarg, prec);
        arb_abs(t.data, t.data);
        arb_exp(t.data, t.data, prec);
        arb_mul(r.data, r.data, t.data, prec);
        arb_set(acb_realref(res), r.data);
        arb_swap(acb_imagref(res), r.data);
    }
    else if (j == 0 && acb_is_zero(a))
    {
        acb_one(res);
    }
    else
    {
        acb_indeterminate(res);
    }
}

/* z^(λ+i)*log(z)^j/j! */
template<class T>
void BranchInfo<T>::zpowlogpow(
    T& R,
    acb_t res,
    slong i, ulong j,
    slong prec)
{
    zpow(R, res, i, prec);
    if (j < 1)
        return;

    xacb_t t1, t2;
    xarb_t t3;
    xacb_t& l = get_logz(R, prec);

    if (acb_is_finite(l.data))
    {
        acb_pow_ui(t2.data, l.data, j, prec);
        acb_mul(t1.data, res, t2.data, prec);
        arb_fac_ui(t3.data, j, prec);
        acb_div_arb(res, t1.data, t3.data, prec);
    }
    else
    {
        R.nset(t1.data, lambda, prec + 60);
        arb_add_si(acb_realref(t1.data), acb_realref(t1.data), i, prec + 60);

        if (prec > 100)
            prec = 10 + sqrt(prec + 1000);

        acb_abs(acb_realref(t2.data), z.get_approx(R, prec).data, prec);
        arb_const_pi(acb_imagref(t2.data), prec);
        fancy_pow_at_zero(res, acb_realref(t2.data), acb_imagref(t2.data), t1.data, j, prec);
    }

//std::cout << "res: " << res << std::endl;
}

// log(z)^j/j!
template <class T>
xacb_t& BranchInfo<T>::get_logzpow(T& R, slong j, slong prec)
{
    assert(j >= 0);
    if (j < 1)
        return logzpow.at(0);

    get_logz(R, prec);

    while (j >= logzpow.size())
    {
        ulong n = logzpow.size();
        logzpow.emplace_back();
        acb_pow_ui(logzpow.back().data, logzpow[1].data, n, prec);
        xarb_t f;
        arb_fac_ui(f.data, n, prec);
        acb_div_arb(logzpow.back().data, logzpow.back().data, f.data, prec);
    }

    return logzpow.at(j);
}

// z*log(z)^j/j!
template <class T>
xacb_t& BranchInfo<T>::get_zlogzpow(T& R, slong j, slong prec)
{
    assert(j >= 0);

    // z = get_approx(Z, Φ)
    // j<1 && return z
    xacb_t& z1 = z.get_approx(R, prec);
    if (j < 1)
        return z1;

    xacb_t& l = get_logz(R, prec);
    if (acb_is_finite(l.data))
    {
        while (j > zlogzpow.size())
        {
            zlogzpow.emplace_back();
            ulong n = logzpow.size();
            acb_pow_ui(zlogzpow.back().data, l.data, n, prec);
            acb_mul(zlogzpow.back().data, zlogzpow.back().data, z1.data, prec);
            xarb_t f;
            arb_fac_ui(f.data, n, prec);
            acb_div_arb(zlogzpow.back().data, zlogzpow.back().data, f.data, prec);
        }
    }
    else
    {
        xarb_t r, pi;
        xacb_t one;
        acb_abs(r.data, z1.data, prec);
        arb_const_pi(pi.data, 60);
        acb_one(one.data);
        while (j > zlogzpow.size())
        {
            zlogzpow.emplace_back();
            fancy_pow_at_zero(zlogzpow.back().data, r.data, pi.data, one.data,
                              zlogzpow.size(), prec);
        }
    }

    return zlogzpow.at(j - 1);
}

// z^(λ+n)*Σ_{j<τ}e_j*Λ^(τ-1-j) with Λ^(τ-1-j) replaced by log(z)^j/j!
template <class T>
void BranchInfo<T>::zpowlogpoly_eval(
    T& R,
    acb_t res,
    slong n,
    std::vector<xacb_t>& e,
    slong prec,
    acb_t t1) // temp
{
    slong tau = e.size();
    slong j = tau - 1;
    if (j < 1)
    {
        zpow(R, res, n, prec);
        acb_mul(res, res, e[0].data, prec);
        return;
    }

    bool ok = acb_is_finite(get_logz(R, prec).data);
    xacb_t& p = ok ? get_logzpow(R, j, prec) : get_zlogzpow(R, j, prec);
    acb_mul(res, e[tau-1-j].data, p.data, prec);


    while (--j > 0)
    {
        xacb_t& p2 = ok ? get_logzpow(R, j, prec) : get_zlogzpow(R, j, prec);
        acb_addmul(res, e.at(tau-1-j).data, p2.data, prec);
    }

    if (ok)
        acb_add(res, res, e.at(tau-1).data, prec);
    else
        acb_addmul(res, e.at(tau-1).data, z.get_approx(R, prec).data, prec);

    zpow(R, t1, n - !ok, prec);
    acb_mul(res, res, t1, prec);
}


/* evaluation of series solutions to diff equations */


template <class T>
void next_u_coeffs(
    T& R,
    std::vector<std::vector<std::vector<T_elem_t>>>& u,
    slong & tau,
    LEquation<T_elem_t>& P,
    T_elem_t& lambda,
    std::vector<xfmpz_t>& ns,
    slong & noff,
    slong N)
{
//std::cout << "---- next_u_coeffs called ------------" << std::endl;
//std::cout << "tau: " << tau << std::endl;
//std::cout << "u: " << u << std::endl;
//std::cout << "P: " << P << std::endl;

    slong nu = ns.size();
    T_elem_t lambdaN;
    R.add(lambdaN, lambda, N);
    bpoly<T_elem_t> Pn;
    bpoly_taylor_shift_inside<T>(R, Pn, coeffs_xoutside<T>(R, P), lambdaN);

//std::cout << "Pn: " << Pn << std::endl;

    slong s = bpoly_degree_outside(Pn);

//std::cout << "s: " << s << std::endl;

//std::cout << "iv: " << std::endl;

    std::vector<std::vector<T_elem_t>> uN;
    slong mu = 0;
    while (noff < nu && fmpz_equal_si(ns[noff].data, N))
    {
        uN.push_back(std::vector<T_elem_t>(nu));
        for (slong l = 0; l < nu; l++)
            R.set(uN.at(mu).at(l), slong(l == noff));
        mu++;
        noff++;
    }

//std::cout << "uN: " << uN << std::endl;

//std::cout << "mu: " << mu << std::endl;

    // mu = number of initial conditions at lambda + N = multiplicity of root
    for (slong i = 0; i < mu; i++)
        assert(i >= Pn.coeffs.at(0).coeffs.size() || R.is_zero(Pn.coeffs.at(0).coeffs.at(i)));

//std::cout << "rhs" << std::endl;
    uN.resize(mu + tau);
    for (slong j = 0; j < tau; j++)
    {
        uN.at(mu + j).resize(nu);
        for (slong l = 0; l < nu; l++)
            R.zero(uN.at(mu + j).at(l));
    }

    for (slong i = 0; i < s; i++)
    {
        auto & Pni = Pn.coeffs.at(i + 1).coeffs;
        for (slong k = 0; k < Pni.size(); k++)
        for (slong j = 0; j + k < u.at(i).size(); j++)
        for (slong l = 0; l < nu; l++)
            R.msub(uN.at(mu + j).at(l),
                          u.at(i).at(j + k).at(l),
                          Pni.at(k));
    }

//std::cout << "uN: " << uN << std::endl;

//std::cout << "lhs" << std::endl;
    for (slong i = 1; i <= tau; i++)
    {
        for (slong j = 1; j < i; j++)
        for (slong l = 0; l < nu; l++)
            R.msub(uN.at(mu + tau - i).at(l),
                          uN.at(mu + tau - i + j).at(l),
                          Pn.coeffs.at(0).coeffs.at(mu + j));

        for (slong l = 0; l < nu; l++)
            R.divexact(uN.at(mu + tau - i).at(l),
                       uN.at(mu + tau - i).at(l),
                       Pn.coeffs.at(0).coeffs.at(mu));
    }

//std::cout << "triming zeros" << std::endl;
    slong k = mu;
    for (; k > 0; k--)
    {
        bool allzero = true;
        for (slong l = 0; l < nu; l++)
            allzero = allzero && R.is_zero(uN.at(tau + k - 1).at(l));
        if (!allzero)
            break;
    }
    tau += k;
    uN.resize(tau);

    std::swap(u.at(s - 1), uN);
    for (slong i = s - 1; i > 0; i--)
        std::swap(u.at(i), u.at(i - 1));

//std::cout << "next_u_coeffs returning" << std::endl;
//std::cout << "tau: " << tau << std::endl;
//std::cout << "u: " << u << std::endl;
}

template <class T>
void u_derivative(
    T& R,
    std::vector<std::vector<T_elem_t>>& uN,
    T_elem_t& fac,
    slong tau,
    slong nu)
{
    for (slong j = 0; j < tau; j++)
    for (slong l = 0; l < nu; l++)
    {
        R.mul(uN[j][l], uN[j][l], fac);
        if (j + 1 < tau)
            R.add(uN[j][l], uN[j][l], uN[j + 1][l]);
    }
}

#if 0
template <class T>
void tail_bound(
    T& R,
    std::vector<xarb_poly_t>& res,
    LEquation<T_elem_t>& P,
    SingularLocus<T_elem_t>& T,
    std::vector<xacb_t>& lambda,
    std::vector<xfmpz>& ns,
    slong nsoffs,


    tail_bound<T>(R, Er, P, u, ns, nsoff, alphas, Z, delta, N, tau, nu, prec);
#endif

template <class T>
void eval_basis_forward(
    T& R,
    xacb_mat_t& M,              /* result */
    LEquation<T_elem_t>& P,
    BranchInfo<T>& Z,
    std::vector<xfmpz_t>& ns,
    std::vector<T_elem_t>& alphas,
    slong prec)
{
    slong delta = M.nrows();

//std::cout << "eval_basis_forward called" << std::endl;
//std::cout << "P: " << tostring(P) << std::endl;
//std::cout << "z: " << z.tostring() << std::endl;
////std::cout << "lambda: " << tostring(lambda) << std::endl;
////std::cout << "ns: " << tostring(ns) << std::endl;
//std::cout << "alphas: " << tostring(alphas) << std::endl;
//std::cout << "delta: " << delta << std::endl;
//std::cout << "prec: " << prec << std::endl;

    slong nu = ns.size();   // number of initial conditions
    slong tau = 0;          // strict bound on the power of log(z) thus far
    slong maxN = 20*prec;   // max number of terms to sum
    assert(fmpz_is_zero(ns[0].data));
    assert(M.ncols() == ns.size());

    auto& Ptheta = coeffs_xoutside<T>(R, P);

//std::cout << "Ptheta: " << tostring(Ptheta) << std::endl;

    slong s = bpoly_degree_outside<T_elem_t>(Ptheta);

//std::cout << "s: " << s << std::endl;

    // u is a nested array of the last s solutions
    assert(s > 0 && nu > 0);
    std::vector<std::vector<std::vector<T_elem_t>>> u(s);
    std::vector<std::vector<T_elem_t>> uN;
    std::vector<T_elem_t> uN0;
    T_elem_t lambdaN;
    xacb_t t1, t2, ts;

    slong N = 0;

    // M is the matrix answer
    acb_mat_zero(M.data);

    // the nsoff reads through the ns as we sum past them
    slong nsoff = 0;

    int changed_count = 0;
    while (true)
    {
        bool may_stop = nsoff >= nu;
        assert(N <= maxN);

//std::cout << "N: " << N << std::endl;

        next_u_coeffs<T>(R, u, tau, P, Z.lambda, ns, nsoff, N);
        uN = u[0];

//std::cout << "uN: " << uN << std::endl;

        // add uN*z^N to sum
        bool changed = false;
        for (slong d = 0; d < delta; d++)
        {
            if (d > 0)
            {
                R.add(lambdaN, Z.lambda, N - d + 1);
                u_derivative<T>(R, uN, lambdaN, tau, nu);
            }

            for (slong l = 0; l < nu; l++)
            {
                acb_zero(ts.data);
                for (slong j = 0; j < tau; j++)
                {
                    if (R.is_zero(uN[j][l]))
                        continue;
                    R.nset(t1.data, uN[j][l], prec);
                    Z.zpowlogpow(R, t2.data, N - d, j, prec);
                    acb_addmul(ts.data, t1.data, t2.data, prec);
                }

                acb_add(ts.data, ts.data, &M[d][l], prec);
                changed = changed || !acb_overlaps(ts.data, &M[d][l]);
                acb_swap(&M[d][l], ts.data);
            }
        }

        N++;

        if (!may_stop)
            continue;

        if (!changed)
        {
            if (++changed_count > 10)
                break;
        }
        else
        {
            changed_count = 0;
        }
    }
#if 0
    std::vector<xarb_poly_t> Er;
    tail_bound<T>(R, Er, P, u, ns, nsoff, alphas, Z, delta, N, tau, nu, prec);
#endif
    xfmpz_t fd(1);
    for (slong d = 0; d < delta; d++)
    {
        if (d > 1)
            fmpz_mul_si(fd.data, fd.data, d);

        for (slong l = 0; l < nu; l++)
        {
            if (d > 1)
                acb_div_fmpz(&M[d][l], &M[d][l], fd.data, prec);
#if 0
            if (d < Er[l].data->length)
            {
                arb_add_error(acb_realref(&M[d][l]), Er[l].data->coeffs + d);
                arb_add_error(acb_imagref(&M[d][l]), Er[l].data->coeffs + d);
            }
#endif
        }
    }
}


template <class T>
void eval_basis_bs(
    T& R,
    xacb_mat_t& M,
    LEquation<T_elem_t>& P,
    BranchInfo<T>& Z,
    std::vector<xfmpz_t>& ns,
    std::vector<T_elem_t>& alphas,
    slong prec)
{
    slong delta = M.nrows();

//std::cout << "eval_basis_bs called" << std::endl;
//std::cout << "P: " << P << std::endl;
//std::cout << "z: " << z << std::endl;
//std::cout << "lambda: " << Z.lambda << std::endl;
//std::cout << "ns: " << ns << std::endl;
//std::cout << "alphas: " << alphas << std::endl;
//std::cout << "delta: " << delta << std::endl;
//std::cout << "prec: " << prec << std::endl;

    slong nu = ns.size();   // number of initial conditions
    slong tau = 0;          // strict bound on the power of log(z) thus far
    slong maxN = 20*prec;   // max number of terms to sum
    assert(fmpz_is_zero(ns[0].data));
    assert(M.ncols() == ns.size());

    auto& Ptheta = coeffs_xoutside<T>(R, P);

//std::cout << "Ptheta: " << tostring(Ptheta) << std::endl;

    slong s = bpoly_degree_outside<T_elem_t>(Ptheta);
    assert(s > 0);

    // u is a nested array of the last s solutions
    // each solution is an array of coefficients of log^k(z)/k!  0 <= k < τ
    // each coefficient is an array of ν elements of F
    assert(s > 0 && nu > 0);
    std::vector<std::vector<std::vector<T_elem_t>>> u(s);
    std::vector<std::vector<T_elem_t>> uN;
    std::vector<T_elem_t> uN0(s*nu*1);
    T_elem_t lambdaN;
    xacb_t t1, t2, ts;

    // M is the matrix answer
    acb_mat_zero(M.data);

    // the nsoff reads through the ns as we sum past them
    slong nsoff = 0;

    slong N = 0;
    slong N0 = -1;

    while (true)
    {
        if (N > maxN)
        {
            std::cout<< "internal error: parameters too large" << std::endl;
            abort();
        }

        next_u_coeffs<T>(R, u, tau, P, Z.lambda, ns, nsoff, N);
        uN = u[0];

        // add uN*z^N to sum
        for (slong d = 0; d < delta; d++)
        {
            if (d > 0)
            {
                R.add(lambdaN, Z.lambda, N - d + 1);
                u_derivative<T>(R, uN, lambdaN, tau, nu);
            }

            for (slong l = 0; l < nu; l++)
            {
                acb_zero(ts.data);
                for (slong j = 0; j < tau; j++)
                {
                    if (R.is_zero(uN[j][l]))
                        continue;
                    R.nset(t1.data, uN[j][l], prec);
                    Z.zpowlogpow(R, t2.data, N - d, j, prec);
                    acb_addmul(ts.data, t1.data, t2.data, prec);
                }

                acb_add(&M[d][l], &M[d][l], ts.data, prec);
            }
        }

        N++;

#define IDX(x, y, z) [((x)*nu + (y))*tau + (z)]
        if (nsoff >= ns.size() && N >= delta)
        {
            N0 = N;
            // for u_{N0-1},...,u_{N0-s}, map the 3D array u of F elem's
            // into 2D array uN0 of poly's in Lambda (it's still 3D)
            uN0.resize(s*nu*tau);
            for (slong i = 0; i < s; i++)
            for (slong l = 0; l < nu; l++)
            {
                slong j = 0;
                for (; j < u.at(i).size(); j++)
                    R.set(uN0 IDX(i,l,tau-1-j), u.at(i).at(j).at(l));
                for (; j < tau; j++)
                    R.zero(uN0 IDX(i,l,tau-1-j));
            }

            break;
        }
    }

//std::cout << "N: " << N << std::endl;

//std::cout << "M: " << M << std::endl;
//std::cout << "uN0: " << uN0 << std::endl;


    // set up sum for [N0, N1)
    std::vector<upoly<T_elem_t>> numeq(s*tau);
    std::vector<upoly<T_elem_t>> Pn(s + 1);
    std::vector<upoly<T_elem_t>> denpow(tau + 1);
    upoly<T_elem_t> p1, p2, p3, dPn0;
    for (slong i = 0; i <= s; i++)
        upoly_taylor_shift<T>(R, Pn.at(i), Ptheta.coeffs.at(i), Z.lambda);
    upoly_one<T>(R, denpow.at(0));
    for (slong i = 0; i < tau; i++)
        upoly_mul<T>(R, denpow.at(i + 1), denpow.at(i), Pn.at(0));

    upoly_derivative<T>(R, dPn0, Pn.at(0));
    for (slong j = 0; j < s; j++)
    {
        auto& d = Pn.at(j + 1);
        upoly_neg<T>(R, d);
        upoly_mul<T>(R, numeq.at(j*tau + 0), d, denpow.at(tau - 1));
        for (slong i = 1; i < tau; i++)
        {            
            upoly_derivative<T>(R, p1, d);
            upoly_mul<T>(R, p2, p1, Pn.at(0));
            upoly_divexact_scalar<T>(R, p3, p2, i);            
            upoly_mul<T>(R, p2, d, dPn0);
            upoly_sub<T>(R, d, p3, p2);
            upoly_mul<T>(R, numeq.at(j*tau+i), d, denpow.at(tau - 1 - i));
        }
    }

    // sum [N0, N1)
    slong N1 = N0 + 10;
    T_int_ring_t ZR;
    SumBsCtx<T_int_ring_t> S;

//std::cout << "starting N0 = " << N0 << ", N1 = " << N1 << std::endl;
//std::cout << "numeq: " << numeq << std::endl;
//std::cout << "deneq: " << denpow.at(tau) << std::endl;
    sum_bs_start<T>(R, ZR, S, numeq, s, tau, denpow.at(tau),
                    Z.z, Z.lambda, delta, N1, N0, prec);

//std::cout << " get sum: " << S.get_sum(ZR, prec) << std::endl;

    slong overlap_count = 0;
    slong d1 = prec;
    slong NDelta = N1;
    do {
        slong N2 = N1 + std::min(std::max(NDelta, (N1+100)/64), N1);
        slong d2 = sum_bs_continue<T_int_ring_t>(ZR, S, N2, N1, N0, prec);
        if (d2 < 1)
        {
            overlap_count++;
            NDelta = 10;
        }
        else
        {
            overlap_count = 0;
            if (d1 - d2 <= 0)
                NDelta = N1;
            else
                NDelta = 10 + ((N2 - N1)*d2*1)/(6*(d1 - d2));
        }
        N1 = N2;
        d1 = d2;
    } while (N1 < maxN && overlap_count < 2);

std::cout << "bs used N1 = " << N1 << std::endl;

    // accumulate sums for [N0, N1) and add to sum for [0,N0)
#define IDX2(x, y, z) [((x)*s + (y))*tau + (z)].data
    std::vector<xacb_t>& Sigma1a = S.get_sum(ZR, prec); // d x s x tau

//std::cout << "Sigma1a: " << Sigma1a << std::endl;
//std::cout << "M: " << M << std::endl;

    std::vector<xacb_t> fdz(tau);
    for (slong d = 0; d < delta; d++)
    for (slong l = 0; l < nu; l++)
    {
        // fdz = z^(d-N0-λ)*f^(d)(z) for the l^th initial value
        for (slong i = 0; i < s; i++)
        for (slong k = 0; k < tau; k++)
        for (slong h = 0; h <= k; h++)
        {
            if (i == 0 && h == 0)
            {
                R.nmul(fdz[k].data, Sigma1a IDX2(d,i,h), uN0 IDX(i,l,k-h), prec);
            }
            else
            {
                R.nmul(t1.data, Sigma1a IDX2(d,i,h), uN0 IDX(i,l,k-h), prec);
                acb_add(fdz[k].data, fdz[k].data, t1.data, prec);
            }
        }
//std::cout << "fdz: " << fdz << std::endl;
//std::cout << "N0-d: " << N0-d << std::endl;
        Z.zpowlogpoly_eval(R, t1.data, N0-d, fdz, prec, t2.data);
//std::cout << "log_poly_eval = " << t1 << std::endl;

        acb_add(&M[d][l], &M[d][l], t1.data, prec);
    }

//std::cout << "M: " << M << std::endl;


#undef IDX2

    // compute u_{N1-1},...,u_{N1-s}
    //uN1 = mullow(get_approx(S.M, Φ), uN0, Φ)
    //uN1t = nacb[uN1[i,l,1+τ-j] for i in 1:s, j in 1:τ, l in 1:ν]

    //// uN1 is s x ? x ?
    //mullow(uN1, S.M.get_approx(R, prec), uN0, prec);
    ////uN1t is s x tau x nu with
    ////uN1t[i,j,l] = uN1[i,l,τ-j-1]
    

  // bound sum for [N1, ∞) and add to sum for [0,N1)
/*
  Er = tail_bound(P, L, uN1t, ns, αs, Z, δ, N1, τ, ν, Φ, F)
  for d in 0:δ-1, l in 1:ν
    getindex!(t, M, 1+d,l)
    div!(t, t, factorial(ZZ(d)), Φ)
    add_error!(t, coeff(Er[l], d))
    setindex!(M, t, 1+d,l)
  end
*/

#if 0
    std::vector<xarb_poly_t> Er;
    tail_bound<T>(R, Er, P, u, ns, nsoff, alphas, Z, delta, N, tau, nu, prec);
#endif
    xfmpz_t fd(1);
    for (slong d = 0; d < delta; d++)
    {
        if (d > 1)
            fmpz_mul_si(fd.data, fd.data, d);

        for (slong l = 0; l < nu; l++)
        {
            if (d > 1)
                acb_div_fmpz(&M[d][l], &M[d][l], fd.data, prec);
#if 0
            if (d < Er[l].data->length)
            {
                arb_add_error(acb_realref(&M[d][l]), Er[l].data->coeffs + d);
                arb_add_error(acb_imagref(&M[d][l]), Er[l].data->coeffs + d);
            }
#endif
        }
    }

//std::cout << "M: " << M << std::endl;

#undef IDX1
}


template <class T>
void eval_basis(
    T& R,
    xacb_mat_t& M,
    LEquation<T_elem_t>& P,
    BranchInfo<T>& Z,
    std::vector<xfmpz_t>& ns,
    std::vector<T_elem_t>& alphas,
    slong prec)
{
//std::cout << std::endl << "*******" << std::endl;
//std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1" << std::endl;
//std::cout << "eval basis called" << std::endl;

#if 0

#if 0
    xacb_mat_t M1(M.nrows(), M.ncols());
    eval_basis_forward<T>(R, M1, P, Z, ns, alphas, prec);
#endif
//std::cout << "forward M: " << M1 << std::endl;

    eval_basis_bs<T>(R, M, P, Z, ns, alphas, prec);
//std::cout << "     bs M: " << M << std::endl;

#if 0
    for (slong i = 0; i < M1.nrows(); i++)
    for (slong j = 0; j < M1.ncols(); j++)
    {
        assert(acb_is_finite(&M[i][j]));
    }
#endif

#else

    eval_basis_forward<T>(R, M, P, Z, ns, alphas, prec);

#endif
}


/** hypergeoemtric ************************************************************/

template <class T>
void compute_pfq_near_zero_taylor(
    T& R,
    acb_t res,
    std::vector<T_elem_t> a,
    std::vector<T_elem_t> b,
    T_elem_t z,
    slong prec)
{
    size_t p = a.size();
    size_t q = b.size();
    assert(q + 1 >= p);

    LEquation<T_elem_t> eq;
    pfq_equ_at_zero<T>(R, eq, a, b);

    std::vector<T_elem_t> roots(q + 1);
    R.zero(roots[0]);
    for (size_t i = 0; i < q; i++)
        R.sub(roots[i + 1], slong(1), b[i]);

    std::vector<xfmpz_t> offs(1);
    fmpz_zero(offs[0].data);
    BranchInfo<T> Z(roots[0], z);
    xacb_mat_t e(1, 1);
    eval_basis<T>(R, e, eq, Z, offs, roots, prec + 30);
    acb_swap(res, &e[0][0]);
}

template <class T>
void compute_pfq_anywhere_taylor(
    T& R,
    acb_t res,
    std::vector<T_elem_t> a,
    std::vector<T_elem_t> b,
    T_elem_t z,
    slong prec)
{
    size_t p = a.size();
    size_t q = b.size();
    assert(q + 1 == p);
    slong k = 0;

    LEquation<T_elem_t> eq;
    pfq_equ_at_zero_quadratic<T>(R, eq, k, a, b);

    xacb_t w, s, t, t1;
    R.nset(s.data, z, prec + 60);
    acb_sub_ui(t1.data, s.data, 1, prec + 60);
    acb_neg(t1.data, t1.data);
    acb_sqrt(t.data, t1.data, prec + 60);
    acb_add_ui(s.data, t.data, 1, prec + 60);
    acb_sub_ui(t1.data, t.data, 1, prec + 60);
    acb_neg(t1.data, t1.data);
    acb_div(w.data, t1.data, s.data, prec + 30);

    std::vector<T_elem_t> roots(q + 1);
    R.zero(roots[0]);
    for (size_t i = 0; i < q; i++)
        R.sub(roots[i + 1], slong(1), b[i]);
    std::vector<xfmpz_t> offs(1);
    fmpz_zero(offs[0].data);
    BranchInfo<T> Z(roots[0], w.data);
    xacb_mat_t e(1, 1);
    eval_basis<T>(R, e, eq, Z, offs, roots, prec + 30);

    R.nset(t.data, a[k], prec + 20);
    acb_mul_2exp_si(t.data, t.data, 1);
    acb_neg(t.data, t.data);
    acb_mul_2exp_si(s.data, s.data, -1);
    acb_pow(t1.data, s.data, t.data, prec + 20);
    acb_mul(res, t1.data, &e[0][0], prec + 20);
}

/*
for each s = λ+ns[j] compute the coeff of z^s in
residue Γ(B)/Γ(A)*Γ(A+s)*Γ(-s)/Γ(B+s)*z^-s,  A=(λ+ns,arest), B=(b)
the coeff is returned as the coeffs of 1, log(z), ..., log^(μ-1)(z)/(μ-1)
TODO: currently doesn't work in terminating case
*/
template <class T>
void pfq_initial_values_at_infinity(
    T& R,
    std::vector<xacb_t>& iv,
    T_elem_t& lambda,
    std::vector<xfmpz_t>& ns,
    std::vector<T_elem_t*> arest,
    std::vector<T_elem_t> a,
    std::vector<T_elem_t> b,
    slong prec)
{
    racb_poly_t Rx(1, prec);
    xacb_poly_t s0, s1, s2, st;
    xacb_t f;
    bs_product<racb_poly_t> P(Rx);

    acb_poly_fit_length(st.data, 2);
    acb_one(st.data->coeffs + 0);
    acb_one(st.data->coeffs + 1);
    _acb_poly_set_length(st.data, 2);

    acb_one(f.data);
    for (slong k = 0; k < b.size(); k++)
    {
        R.nset(f.data, b[k], prec);
        acb_gamma(f.data, f.data, prec);
        acb_mul(st.data->coeffs + 0, st.data->coeffs + 0, f.data, prec);
    }
    for (slong k = 0; k < a.size(); k++)
    {
        R.nset(f.data, a[k], prec);
        acb_rgamma(f.data, f.data, prec);
        acb_mul(st.data->coeffs + 1, st.data->coeffs + 1, f.data, prec);
    }
    acb_mul(f.data, st.data->coeffs + 0, st.data->coeffs + 1, prec);

    T_elem_t lambdan, t;
    xfmpz_t d;
    slong m = ns.size();
    slong mu;
    iv.resize(m);
    for (slong j = 0; j < m; j += mu)
    {
        const xfmpz_t& n = ns[j];
        // μ = multiplicity of exponent λ+n
        mu = 1;
        while (j + mu < m && fmpz_equal(ns[j + mu].data, n.data))
            mu++;

        // compute coeffs iv[j], ..., iv[j+μ-1] of z^(λ+n), ..., z^(λ+n)*log^(μ-1)(z)/(μ-1)!
        // there is a pole of order j + mu at s = -λ-n + ε
        Rx.ord = j + mu;
        P.one(Rx);

        // handle poles from s = -λ-n + ε
        ulong sign = 0;
        for (slong k = 0; k < j + mu; k++)
        {
            // Γ(λ+ns[k]+s) = Γ(ns[k]-ns[j] + ε)
            //              = Γ(-d + ε)
            //              = (-1)^d Γ(1+ε)/(ε*(1-ε)_d)
            fmpz_sub(d.data, n.data, ns[k].data);
            assert(fmpz_sgn(d.data) >= 0);
            if (!fmpz_abs_fits_ui(d.data))
                assert(false && "not implemented");
            acb_one(st.data->coeffs + 0);
            acb_one(st.data->coeffs + 1);
            acb_poly_gamma_series(s1.data, st.data, Rx.ord, prec);
            acb_neg(st.data->coeffs + 1, st.data->coeffs + 1);
            ulong dd = fmpz_get_ui(d.data);
            sign += dd;
            acb_poly_rising_ui_series(s2.data, st.data, dd, Rx.ord, prec);
            acb_poly_div_series(s0.data, s1.data, s2.data, Rx.ord, prec);
            P.push(Rx, s0);
        }

        // handle non-poles from s = -λ-n + ε
        for (slong k = j + mu; k < m; k++)
        {
            fmpz_sub(d.data, ns[k].data, n.data);
            assert(fmpz_sgn(d.data) > 0);
            acb_set_fmpz(st.data->coeffs + 0, d.data);
            acb_one(st.data->coeffs + 1);
            acb_poly_gamma_series(s0.data, st.data, Rx.ord, prec);
            P.push(Rx, s0);
        }

        // handle Γ(-s) = Γ(λ+n - ε)
        R.add(lambdan, lambda, n.data);
        R.nset(st.data->coeffs + 0, lambdan, prec);
        acb_set_si(st.data->coeffs + 1, -1);
        acb_poly_gamma_series(s0.data, st.data, Rx.ord, prec);
        P.push(Rx, s0);

        // handle non-poles of Γ(arest[k]+s) = Γ(arest[k]-λ-n + ε)
        for (slong k = 0; k < arest.size(); k++)
        {
            R.sub(t, *arest[k], lambdan);
            R.nset(st.data->coeffs + 0, t, prec);
            acb_one(st.data->coeffs + 1);
            acb_poly_gamma_series(s0.data, st.data, Rx.ord, prec);
            P.push(Rx, s0);
        }

        // handle Γ(b[k]+s) = Γ(b[k]-λ-n + ε)
        for (slong k = 0; k < b.size(); k++)
        {
            R.sub(t, b[k], lambdan);
            R.nset(st.data->coeffs + 0, t, prec);
            acb_one(st.data->coeffs + 1);
            acb_poly_rgamma_series(s0.data, st.data, Rx.ord, prec);
            P.push(Rx, s0);
        }

        // now we have s = ε^-ord*(s[0] + ... s[ord-1]*ε^(ord-1) + O(s^ord))
        //     determine r[j+k] for 0 <= k < μ.
        // Sum[s[0]*ε^j,{j,0,inf}] * Sum[(-1)^k*Log[z]^k/k!*ε^k,{k,0,inf}]
        acb_ptr s = P.product().data->coeffs;
        for (slong k = 0; k < mu; k++)
        {
            // coeff Log[z]^k/k!*ε^(ord-1) has j = ord-1-k
            assert(k < Rx.ord);
            acb_mul(iv[j + k].data, s + Rx.ord - 1 - k, f.data, prec);
            if ((sign + k)&1)
                acb_neg(iv[j + k].data, iv[j + k].data);
        }
    }
}

template <class E>
using tri_tuple = std::tuple<E, std::vector<xfmpz_t>, std::vector<slong>>;


template <class T>
std::vector<tri_tuple<T_elem_t>> partition_mod1(
    T& R,
    std::vector<T_elem_t>& a)
{
    std::vector<tri_tuple<T_elem_t>> r;
    T_elem_t t;
    xfmpz_t k;

    for (slong i = 0; i < a.size(); i++)
    {
        for (slong j = 0; j < r.size(); j++)
        {
            R.sub(t, a[i], std::get<0>(r[j]));
            if (R.is_integer(k, t))
            {
                std::get<1>(r[j]).push_back(std::move(k));
                std::get<2>(r[j]).push_back(i);
                goto continue_outside;
            }
        }
        r.emplace_back(a[i], std::vector<xfmpz_t>{xfmpz_t(0)}, std::vector<slong>{i});
    continue_outside:;
    }

    for (auto & rj : r)
    {
        auto & rj1 = std::get<1>(rj);
        std::sort(rj1.begin(), rj1.end());
        R.add(std::get<0>(rj), std::get<1>(rj)[0].data);
        for (slong l = rj1.size(); l > 0 ; l--)
            fmpz_sub(rj1[l - 1].data, rj1[l - 1].data, rj1[0].data);
    }
    return r;
}

template <class T>
void compute_pfq_near_infinity_taylor(
    T& R,
    acb_t res,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b,
    T_elem_t& z,
    slong prec)
{
//std::cout << "compute_pfq_near_infinity_taylor called" << std::endl;
    slong p = a.size();
    slong q = b.size();
    assert(q + 1 <= p);

    std::vector<T_elem_t *> arest(p);
    std::vector<xacb_t> iv(p);
    xacb_mat_t e(1, 1);

    T_elem_t nz, rz;
    R.neg(nz, z);
    R.inv(rz, nz);
    BranchInfo<T> Z(a[0], rz, nz); // dummy lambda

    LEquation<T_elem_t> eq;
    pfq_equ_at_infinity<T>(R, eq, a, b);

    acb_zero(res);
    std::vector<tri_tuple<T_elem_t>> I(partition_mod1<T>(R, a));
    for (auto & i : I)
    {
        T_elem_t& lambda = std::get<0>(i);
        std::vector<xfmpz_t>& ns = std::get<1>(i);
        std::vector<slong>& indices = std::get<2>(i);
        arest.clear();
        for (slong j = 0; j < p; j++)
            if (std::find(indices.begin(), indices.end(), j) == indices.end())
                arest.push_back(&a[j]);
        pfq_initial_values_at_infinity<T>(R, iv, lambda, ns, arest, a, b, prec + 40);
        e.resize(1, ns.size());
        Z.reset_lambda(lambda);
        eval_basis<T>(R, e, eq, Z, ns, a, prec + 40);
        for (slong j = 0; j < ns.size(); j++)
            acb_addmul(res, &e[0][j], iv[j].data, prec + 40);
    }
}


template <class T>
void compute_pfq_near_one_taylor(
    T& R,
    acb_t res,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b,
    T_elem_t& z,
    slong prec)
{
    slong n = a.size();
    assert(n == b.size() + 1);

    prec += 30;

    T_elem_t w, t;
    R.one(t);
    R.divexact(w, t, std::min(n, slong(8)));

    // 0 -> 1-w
    LEquation<T_elem_t> eq;
    pfq_equ_at_zero<T>(R, eq, a, b);

    std::vector<T_elem_t> roots(n);
    R.zero(roots[0]);
    for (slong i = 1; i < n; i++)
        R.sub(roots[i], slong(1), b[i - 1]);

    std::vector<xfmpz_t> offs(1);
    fmpz_zero(offs[0].data);
    R.sub(t, slong(1), w);
    BranchInfo<T> Z(roots[0], t);

    xacb_mat_t e0(n, 1);
    xacb_mat_t e1(1, n);
    xacb_mat_t e2(n, n);
    xacb_mat_t et(n, 1);
    xacb_mat_t ei(n, 1);

    eval_basis<T>(R, e0, eq, Z, offs, roots, prec);
    for (slong i = 1; i < n; i += 2)
        acb_neg(&e0[i][0], &e0[i][0]);

    T_elem_t onemz;
    R.sub(onemz, 1, z);

    pfq_equ_at_one<T>(R, eq, a, b);

    // roots of indicial equation are sigma and 0:n-2
    R.neg(roots[0], a[n - 1]);
    for (slong i = 0; i < n - 1; i++)
    {
        R.set(roots.at(i + 1), i);
        R.add(roots[0], b[i]);
        R.sub(roots[0], a[i]);
    }

    offs.resize(n);
    for (slong i = 0; i < n - 1; i++)
        fmpz_set_ui(offs[i].data, i);

    bool intcase = R.is_integer(offs[n - 1], roots[0]);
    if (n < 2 || intcase)
    {
        std::sort(offs.begin(), offs.end());

        if (n < 2)
            R.set(t, roots[0]);
        else
            R.set(t, offs[0].data);

        for (slong i = n - 1; i >= 0; i--)
            fmpz_sub(offs[i].data, offs[i].data, offs[0].data);

        Z.reset(t, w);
        eval_basis<T>(R, e2, eq, Z, offs, roots, prec);

        Z.reset(t, onemz);
        eval_basis<T>(R, e1, eq, Z, offs, roots, prec);

        et.resize(1, 1);
    }
    else
    {
        e2.hshrink(n - 1);
        e1.hshrink(n - 1);
        offs.resize(n - 1);
        Z.reset(roots.at(1), w);
        eval_basis<T>(R, e2, eq, Z, offs, roots, prec);

        Z.reset(roots.at(1), onemz);
        eval_basis<T>(R, e1, eq, Z, offs, roots, prec);

        offs.resize(1);
        Z.reset(roots.at(0), w);
        eval_basis<T>(R, et, eq, Z, offs, roots, prec);
        e2.hcat(et.data);

        et.resize(1, 1);
        Z.reset(roots.at(0), onemz);
        eval_basis<T>(R, et, eq, Z, offs, roots, prec);
        e1.hcat(et.data);
    }

    // e1 * e2^-1 * e0
    if (0 == acb_mat_solve(ei.data, e2.data, e0.data, prec))
    {
        acb_indeterminate(res);
        return;
    }
    acb_mat_mul(et.data, e1.data, ei.data, prec);
    acb_swap(res, &et[0][0]);
}


template <class T>
void compute_pfq_plain_terms_horner(
    T& R,
    acb_t res,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b,
    T_elem_t& z,
    slong N,
    slong prec)
{
    T_elem_t t1, t2, t3;
    upoly<T_elem_t> A, B;

    upoly_set<T>(R, A, z);
    for (auto & ai : a)
        upoly_mul_linear<T>(R, A, ai);

    upoly_one<T>(R, B);
    for (auto & bi : b)
        upoly_mul_linear<T>(R, B, bi);

    acb_one(res);
    for (slong i = N - 1; i > 0; i--)
    {
        upoly_evaluate<T>(R, t1, A, i - 1);
        upoly_evaluate<T>(R, t2, B, i - 1);
        R.mul(t2, i);
        R.divexact(t3, t1, t2);
        R.nmul(res, res, t3, prec + 10);
        arb_add_ui(acb_realref(res), acb_realref(res), 1, prec + 10);
    }
}

template <class T>
void compute_pfq_plain_terms_bs(
    T& R,
    acb_t res,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b,
    T_elem_t& z,
    slong N,
    slong prec)
{
    LEquation<T_elem_t> eq;
    pfq_equ_at_zero<T>(R, eq, a, b);
    std::vector<T_elem_t> iv(1);
    R.one(iv[0]);
    sum_bs<T>(R, res, N, eq, iv, z, prec + 10 + 2*FLINT_BIT_COUNT(N));
}

template <class T>
void compute_pfq_plain_terms(
    T& R,
    acb_t res,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b,
    T_elem_t& z,
    slong N,
    slong prec)
{
    if (N < 10)
        compute_pfq_plain_terms_horner(R, res, a, b, z, N, prec);
    else
        compute_pfq_plain_terms_bs(R, res, a, b, z, N, prec);
}


template <class T>
slong pfq_terminates_test(
    T& R,
    std::vector<T_elem_t>& a,
    slong limit)
{
    ulong N = -UWORD(1);
    xfmpz_t i;
    for (auto & ai : a)
    {
        if (!R.is_integer(i, ai))
            continue;
        if (fmpz_sgn(i.data) > 0 || fmpz_cmp_si(i.data, -limit) < 0)
            continue;
        N = std::min(N, ulong(1) - fmpz_get_si(i.data));
    }
    return slong(N);
}

template <class T>
void compute_pfq(
    T& R,
    acb_t res,
    std::vector<T_elem_t>& a,
    std::vector<T_elem_t>& b,
    T_elem_t& z,
    slong prec)
{
    xacb_t zz;
    R.nset(zz.data, z, 60);
    double x = arf_get_d(arb_midref(acb_realref(zz.data)), ARF_RND_DOWN);
    double y = arf_get_d(arb_midref(acb_imagref(zz.data)), ARF_RND_DOWN);
    x = FLINT_MAX(FLINT_MIN(x, 1e5), -1e5);
    y = FLINT_MAX(FLINT_MIN(y, 1e5), -1e5);

    slong N = pfq_terminates_test<T>(R, a, prec);
    if (N > 0)
    {
        compute_pfq_plain_terms<T>(R, res, a, b, z, N, prec);
        return;
    }

    slong p = a.size();
    slong q = b.size();
    slong d = q + 1 - p;
    if (d > 0)
    {
        compute_pfq_near_zero_taylor<T>(R, res, a, b, z, prec);
        return;
    }
    else if (d < 0)
    {
        compute_pfq_near_infinity_taylor<T>(R, res, a, b, z, prec);
        return;
    }
    else
    {
        if (x*x + y*y < 0.7)
            compute_pfq_near_zero_taylor<T>(R, res, a, b, z, prec);
        else if (x*x + y*y > 1.3)
            compute_pfq_near_infinity_taylor<T>(R, res, a, b, z, prec);
        else if (x*x + y*y < 2*x - 0.5) // 2*x + r^2-1
            compute_pfq_near_one_taylor<T>(R, res, a, b, z, prec);
        else
            compute_pfq_anywhere_taylor<T>(R, res, a, b, z, prec);
        return;
    }
}


bool _try_fmpq(acb_t res, er e, slong prec)
{
    std::vector<xfmpq_t> a, b;
    xfmpq_t z;

    if (!from_ex(z, echild(e,3)) || !from_ex(a, echild(e,1)) || !from_ex(b, echild(e,2)))
        return false;

    rfmpq_t r;
    compute_pfq<rfmpq_t>(r, res, a, b, z, prec);
    return true;
}

bool _try_fmpqi(acb_t res, er e, slong prec)
{
    std::vector<xfmpqi_t> a, b;
    xfmpqi_t z;

    if (!from_ex(z, echild(e,3)) || !from_ex(a, echild(e,1)) || !from_ex(b, echild(e,2)))
        return false;

    rfmpqi_t r;
    compute_pfq<rfmpqi_t>(r, res, a, b, z, prec);
    return true;
}

ex ncode_sHypergeometricPFQ(er e, slong prec)
{
//std::cout << "***********************************************" << std::endl;
//std::cout << "ncode_sHypergeometricPFQ(" << prec << "): " << e << std::endl;
    if (!ehas_head_sym_length(e, gs.sym_sHypergeometricPFQ.get(), 3))
        return ecopy(e);

    xacb_t res;
    if (!_try_fmpq(res.data, e, prec) &&
        !_try_fmpqi(res.data, e, prec))
    {
        return ecopy(e);
    }

//std::cout << "xres: " << res << std::endl;
    return emake_cmplx_move(res.data);
}


ex ncode_sMeijerG(er e, slong prec)
{
//std::cout << "ncode_sMeijerG(" << prec << "): " << e << std::endl;
    if (!ehas_head_sym_length(e, gs.sym_sMeijerG.get(), 3))
        return ecopy(e);

    return ecopy(e);
}


// PFQEquation[a, b, x0, x_Symbol, theta_Symbol]
ex dcode_hypPFQEquation(er e)
{
//std::cout << "dcode_hypPFQEquation: " << e << std::endl;
    assert(ehas_head_sym(e, gs.sym_hypPFQEquation.get()));

    if (elength(e) != 5)
        return ecopy(e);

    rfmpq_t r;
    std::vector<xfmpq_t> a, b;

    if (!from_ex(a, echild(e,1)) || !from_ex(b, echild(e,2)))
        return ecopy(e);

    if (eis_zero(echild(e,3)))
    {
        LEquation<xfmpq_t> eq;
        pfq_equ_at_zero<rfmpq_t>(r, eq, a, b);
        complete_xoutside<rfmpq_t>(r, eq.data);
        return emake_list(to_ex(eq.data.xoutside, echild(e,5), echild(e,4)));
    }
    else
    {
        LEquation<xfmpq_t> eq;
        pfq_equ_at_zero<rfmpq_t>(r, eq, a, b);
        complete_xoutside<rfmpq_t>(r, eq.data);
        uex res1(to_ex(eq.data.xoutside, echild(e,5), echild(e,4)));

        REquation<xfmpq_t> eq2;
        pfq_equ_at_infinity<rfmpq_t>(r, eq, a, b);
        convert<rfmpq_t>(r, eq, eq2);
        complete_xoutside<rfmpq_t>(r, eq.data);
        uex res2(to_ex(eq.data.xoutside, echild(e,5), echild(e,4)));

        return emake_list(res1.release(), res2.release());
    }
}
