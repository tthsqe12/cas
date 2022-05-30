#include <cmath>
#include <cfloat>

#include "uex.h"
#include "timing.h"
#include "ex_print.h"
#include "eval.h"
#include "code.h"
#include "hash.h"
#include "arithmetic.h"
#include "polynomial.h"
#include "rings.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_lll.h"
#include "arb_fmpz_poly.h"
#include "xacb_t.h"


enum Algo {Algo_interpolate, Algo_pgcd/*, Algo_newton*/};
static const Algo algo = Algo_pgcd;

/* wp = working precision */

slong binary_wp(const acb_t a, const acb_t b)
{
    slong pa = acb_rel_accuracy_bits(a);
    slong pb = acb_rel_accuracy_bits(b);
    slong p = std::max(pa, pb);
    assert(p < ARF_PREC_EXACT - 100);
    return std::max(p, slong(0)) + 100;
}

slong unary_wp(const acb_t a)
{
    slong p = acb_rel_accuracy_bits(a);
    assert(p < ARF_PREC_EXACT - 100);
    return std::max(p, slong(0)) + 100;
}


/* copy the midpoint and set the radius */
void acb_set_radius(acb_t y, const acb_t x, const mag_t e)
{
    arf_set(arb_midref(acb_realref(y)), arb_midref(acb_realref(x)));
    mag_set(arb_radref(acb_realref(y)), e);
    if (arb_is_zero(acb_imagref(x)))
    {
        arb_zero(acb_imagref(y));
    }
    else
    {
        arf_set(arb_midref(acb_imagref(y)), arb_midref(acb_imagref(x)));
        mag_set(arb_radref(acb_imagref(y)), e);
    }
}

void acb_chop_imag(acb_t x)
{
    if (mag_is_zero(arb_radref(acb_realref(x))))
        mag_set(arb_radref(acb_realref(x)), arb_radref(acb_imagref(x)));
    arb_zero(acb_imagref(x));
}

void arb_limit_precision(arb_t x, slong p)
{
    mag_t t;
    mag_init(t);
    arf_get_mag(t, arb_midref(x));
    mag_mul_2exp_si(t, t, -p);
    if (mag_is_zero(t))
        mag_set_ui_2exp_si(t, 1, -p);
    else if (mag_cmp(t, arb_radref(x)) > 0)
        mag_swap(arb_radref(x), t);
    mag_clear(t);
}


void arb_canonicalize(arb_t x, slong p)
{
    if (!mag_is_zero(arb_radref(x)))
        return;

    arf_get_mag(arb_radref(x), arb_midref(x));
    if (mag_is_zero(arb_radref(x)))
        mag_one(arb_radref(x));
    mag_mul_2exp_si(arb_radref(x), arb_radref(x), -p);
}

void acb_canonicalize(acb_t x, slong p)
{
    if (arb_is_zero(acb_imagref(x)))
    {
        arb_canonicalize(acb_realref(x), p);
    }
    else if (mag_is_zero(arb_radref(acb_realref(x))) && mag_is_zero(arb_radref(acb_imagref(x))))
    {
        arb_canonicalize(acb_realref(x), p);
        arb_canonicalize(acb_imagref(x), p);
    }
}

bool acb_is_canonical(const acb_t x)
{
    if (!acb_is_finite(x))
        return false;

    if (mag_is_zero(arb_radref(acb_realref(x))) && mag_is_zero(arb_radref(acb_imagref(x))))
        return false;

    return true;
}

/* if true, the members of the qbarelem class are valid */
bool newton_test(const fmpz_poly_t f, const acb_t x, mag_t e, slong extra = 0)
{
    assert(acb_is_canonical(x));
    acb_t y; /* TODO: could be shallow copies of other stuff */
    xfmpz_poly_t fp, fpp;
    xacb_t f_eval, fp_eval, fpp_eval;
    xmag_t a0, b0, c, t;
    slong p = unary_wp(x);

    fmpz_poly_derivative(fp.data, f);
    fmpz_poly_derivative(fpp.data, fp.data);

    acb_init(y);

    if (arb_is_zero(acb_imagref(x)))
    {
        arb_fmpz_poly_evaluate_arb(acb_realref(fp_eval.data), fp.data, acb_realref(x), p);
        arb_get_mag_lower(a0.data, acb_realref(fp_eval.data));

        arb_fmpz_poly_evaluate_arb(acb_realref(f_eval.data), f, acb_realref(x), p);
        arb_get_mag(b0.data, acb_realref(f_eval.data));
        mag_div(b0.data, b0.data, a0.data);

        mag_mul_2exp_si(e, b0.data, 1);
        mag_max(e, e, arb_radref(acb_realref(x)));

        arf_set(arb_midref(acb_realref(y)), arb_midref(acb_realref(x)));
        mag_set(arb_radref(acb_realref(y)), e);
        arb_zero(acb_imagref(y));

        arb_fmpz_poly_evaluate_arb(acb_realref(fpp_eval.data), fpp.data, acb_realref(y), p);
        arb_get_mag(c.data, acb_realref(fpp_eval.data));
    }
    else
    {
        arb_fmpz_poly_evaluate_acb(fp_eval.data, fp.data, x, p);
        acb_get_mag_lower(a0.data, fp_eval.data);

        arb_fmpz_poly_evaluate_acb(f_eval.data, f, x, p);
        acb_get_mag(b0.data, f_eval.data);
        mag_div(b0.data, b0.data, a0.data);

        mag_mul_2exp_si(e, b0.data, 1);
        mag_max(e, e, arb_radref(acb_realref(x)));
        mag_max(e, e, arb_radref(acb_imagref(x)));

        arf_set(arb_midref(acb_realref(y)), arb_midref(acb_realref(x)));
        arf_set(arb_midref(acb_imagref(y)), arb_midref(acb_imagref(x)));
        mag_set(arb_radref(acb_realref(y)), e);
        mag_set(arb_radref(acb_imagref(y)), e);
        arb_fmpz_poly_evaluate_acb(fpp_eval.data, fpp.data, y, p);
        acb_get_mag(c.data, fpp_eval.data);
    }

    acb_clear(y);

    mag_mul(c.data, c.data, b0.data);
    mag_div(c.data, c.data, a0.data);

//std::cout << "newton test c: " << c.tostring() << std::endl;

    return mag_is_finite(b0.data) && mag_cmp_2exp_si(c.data, -5 - extra) < 0;
}

void newton_trim(const fmpz_poly_t f, acb_t x, mag_t xe)
{
    xacb_t y;
    xmag_t ye;
    slong p;
    slong lo = 0;
    slong hi = std::max(slong(20), acb_rel_accuracy_bits(x));

again:

    p = (hi + lo)/2;

    if (hi <= 20 || hi - lo < 5)
        return;

    arb_set(acb_realref(y.data), acb_realref(x));
    arb_limit_precision(acb_realref(y.data), p);

    if (arb_is_zero(acb_imagref(x)))
    {
        arb_zero(acb_imagref(y.data));
    }
    else
    {
        arb_set(acb_imagref(y.data), acb_imagref(x));
        arb_limit_precision(acb_imagref(y.data), p);
    }

    if (!newton_test(f, y.data, ye.data, 30))
    {
        lo = p;
        goto again;
    }

    acb_swap(x, y.data);
    mag_swap(xe, ye.data);
    hi = p;
    goto again;
}


/* advance the members of the qbarelem class by one iteration */
void newton_step(const fmpz_poly_t f, acb_t x, mag_t xe)
{
    assert(acb_is_canonical(x));
    xacb_t z;
    xmag_t ze;
    xfmpz_poly_t fp;
    xacb_t t, f_eval, fp_eval;
    slong p = unary_wp(x);

    fmpz_poly_derivative(fp.data, f);

try_again:

    arf_set(arb_midref(acb_realref(z.data)), arb_midref(acb_realref(x)));
    mag_mul_2exp_si(arb_radref(acb_realref(z.data)), arb_radref(acb_realref(x)), -p);
    p = std::min(2*p, slong(ARF_PREC_EXACT));

    if (arb_is_zero(acb_imagref(x)))
    {
        arb_fmpz_poly_evaluate_arb(acb_realref(f_eval.data), f, acb_realref(z.data), p);
        arb_fmpz_poly_evaluate_arb(acb_realref(fp_eval.data), fp.data, acb_realref(z.data), p);
        arb_div(acb_realref(t.data), acb_realref(f_eval.data), acb_realref(fp_eval.data), p);
        arb_sub(acb_realref(z.data), acb_realref(z.data), acb_realref(t.data), p);
    }
    else
    {
        arf_set(arb_midref(acb_imagref(z.data)), arb_midref(acb_imagref(x)));
        mag_mul_2exp_si(arb_radref(acb_imagref(z.data)), arb_radref(acb_imagref(x)), -p);

        arb_fmpz_poly_evaluate_acb(f_eval.data, f, z.data, p);
        arb_fmpz_poly_evaluate_acb(fp_eval.data, fp.data, z.data, p);
        acb_div(t.data, f_eval.data, fp_eval.data, p);
        acb_sub(z.data, z.data, t.data, p);
    }

    if (newton_test(f, z.data, ze.data) && mag_cmp(ze.data, xe) < 0)
    {
        acb_swap(x, z.data);
        mag_swap(xe, ze.data);
        return;
    }

    assert(0 && "shit hit the fan");

    /* something bad happened. try again with greater precision */
    // TODO: and possibly more than one iteration
    p += 10 + p/16;
    goto try_again;
}


class qbarelem {
public:
/*
    The object is valid if mu <= 2^-6, where
        a0 >= |1/f'(x)|
        b0 >= |f(x)/f'(x)|
        e >= max(2*b0, x.rad)
        y := [x.mid +- e]
        c >= |f''(y)|
        mu = 2*a0*b0*c
    If so, newton's method converges to the unique root in y from any starting
    value in x. The 2^-6 can actually be replaced by 1 but is used for sanity.
*/
    xfmpz_poly_t minpoly; // f irreducible, primitive, and positive leadcoeff
    xacb_t location;      // x
    xmag_t epsilon;       // e

    /* expr conversion */
    bool set_ex(er e);
    bool set_exp_ex(er f, er e);
    ex get_ex() const;
//    slong get_index() const;

    /* basic arithmetic */
    bool set(const fmpz_poly_t f, const acb_t x);
    bool set_contains(const fmpz_poly_t f, const acb_t x);
    void set(const qbarelem & a);
    void set(slong a);
    void set(const fmpz_t a);
    void set(const fmpq_t a);
    bool isreal() const;
    bool realpart(const qbarelem & a);
    void neg(const qbarelem & a);
    bool add(const qbarelem & a, const qbarelem & b);
    bool sub(const qbarelem & a, const qbarelem & b);
    bool mul(const qbarelem & a, const qbarelem & b);
    bool div(const fmpz_t a, const fmpz_t b);
    bool div(const qbarelem & a, const qbarelem & b);
    bool pow(const qbarelem & a, slong b, ulong c);
    bool pow(const qbarelem & a, qbarelem & b);
    bool pow(const qbarelem & a, const fmpz_t b);
    bool pow(const qbarelem & a, const fmpq_t b);
    bool set_expipi(const fmpq_t b);
    bool mul_expipi(const qbarelem & a, const fmpq_t b);
    bool apply_moebius(const qbarelem & a, const fmpz_t m11, const fmpz_t m12, const fmpz_t m21, const fmpz_t m22);
    void apply_rat_fxn(const qbarelem & a, const fmpz * num, slong num_length, const fmpz * den, slong den_length);
    void quadraticroot(const fmpz_t a, const fmpz_t b, const fmpz_t c);
    bool circle_root(const qbarelem & a, int sign);

    void swap(qbarelem & a) {
        fmpz_poly_swap(minpoly.data, a.minpoly.data);
        acb_swap(location.data, a.location.data);
        mag_swap(epsilon.data, a.epsilon.data);
    }

    template <typename RefineFunc> bool select_root(fmpz_poly_t p, bool must_factor, RefineFunc T);
    void complete_linear(slong e);
};


std::ostream& operator<<(std::ostream& o, const qbarelem& a)
{
    o << "Root[";
    if (fmpz_poly_degree(a.minpoly.data) > 10)
        o << "degree " << fmpz_poly_degree(a.minpoly.data);
    else
        o << a.minpoly;

    o << ", " << a.location << "]";
    return o;
}



bool qbarelem::set(const fmpz_poly_t f, const acb_t x)
{
    slong the_correct_one = -1;
    xfmpz_poly_factor_t fac;
    fmpz_poly_factor(fac.data, f);

    acb_set(location.data, x);

    if (!acb_is_finite(location.data))
        return false;

    acb_canonicalize(location.data, FLINT_BITS);

    for (slong i = 0; i < fac.data->num; i++)
    {
        if (/*fmpz_poly_length(fac.data->p + i) > mindeg &&*/
            newton_test(fac.data->p + i, location.data, epsilon.data))
        {
            if (the_correct_one >= 0)
                return false;
            the_correct_one = i;
        }
    }

    if (the_correct_one < 0)
        return false;

    fmpz_poly_swap(minpoly.data, fac.data->p + the_correct_one);
    return true;
}

bool qbarelem::set_contains(const fmpz_poly_t f, const acb_t x)
{
    xacb_t t;
    slong prec = unary_wp(x);
    slong the_correct_one = -1;
    xfmpz_poly_factor_t fac;

//printf("set from poly: "); fmpz_poly_print_pretty(f, "#"); flint_printf("\n");

    fmpz_poly_factor(fac.data, f);

    acb_set(location.data, x);

    if (!acb_is_finite(location.data))
        return false;

    acb_canonicalize(location.data, FLINT_BITS);

    for (slong i = 0; i < fac.data->num; i++)
    {
//printf("trying factor: "); fmpz_poly_print_pretty(fac.data->p + i, "#"); flint_printf(" ^ %wd\n", fac.data->exp[i]);
        arb_fmpz_poly_evaluate_acb(t.data, fac.data->p + i, location.data, prec);
        if (!acb_contains_zero(t.data))
            continue;

        if (newton_test(fac.data->p + i, location.data, epsilon.data))
        {
//printf("its good\n");
            if (the_correct_one >= 0)
                return false;
            the_correct_one = i;
        }
    }

    if (the_correct_one < 0)
        return false;

    fmpz_poly_swap(minpoly.data, fac.data->p + the_correct_one);
    return true;
}

void qbarelem::set(const qbarelem & a)
{
    fmpz_poly_set(minpoly.data, a.minpoly.data);
    acb_set(location.data, a.location.data);
}

void qbarelem::complete_linear(slong e)
{
    /* make sure x.rad >= 2^e */
    assert(minpoly.data->length == 2);
    if (mag_cmp_2exp_si(arb_radref(acb_realref(location.data)), e) < 0)
        mag_set_ui_2exp_si(arb_radref(acb_realref(location.data)), 1, e);
    mag_mul_2exp_si(epsilon.data, arb_radref(acb_realref(location.data)), 1);
    arb_zero(acb_imagref(location.data));
}

void qbarelem::set(slong a)
{
    fmpz_poly_fit_length(minpoly.data, 2);
    fmpz_set_si(minpoly.data->coeffs + 0, a);
    fmpz_neg(minpoly.data->coeffs + 0, minpoly.data->coeffs + 0);
    fmpz_one(minpoly.data->coeffs + 1);
    _fmpz_poly_set_length(minpoly.data, 2);

    arf_set_si(arb_midref(acb_realref(location.data)), a);
    mag_set_ui_2exp_si(arb_radref(acb_realref(location.data)), 1, -FLINT_BITS);
    complete_linear(-FLINT_BITS);
}

void qbarelem::set(const fmpz_t a)
{
    fmpz_poly_fit_length(minpoly.data, 2);
    fmpz_neg(minpoly.data->coeffs + 0, a);
    fmpz_one(minpoly.data->coeffs + 1);
    _fmpz_poly_set_length(minpoly.data, 2);

    arf_set_fmpz(arb_midref(acb_realref(location.data)), a);
    mag_set_ui_2exp_si(arb_radref(acb_realref(location.data)), 1, -FLINT_BITS);
    complete_linear(-FLINT_BITS);
}

void qbarelem::set(const fmpq_t a)
{
    div(fmpq_numref(a), fmpq_denref(a));
}

bool qbarelem::div(const fmpz_t a, const fmpz_t b)
{
    if (fmpz_is_zero(b))
        return false;

    fmpz_poly_fit_length(minpoly.data, 2);
    fmpz_neg(minpoly.data->coeffs + 0, a);
    fmpz_set(minpoly.data->coeffs + 1, b);
    _fmpz_poly_set_length(minpoly.data, 2);
    _fmpq_canonicalise(minpoly.data->coeffs + 0, minpoly.data->coeffs + 1);

    slong p = (slong)fmpz_bits(a) - (slong)fmpz_bits(b);
    arb_fmpz_div_fmpz(acb_realref(location.data), a, b, 2*FLINT_BITS);

    complete_linear(p - 2*FLINT_BITS);

    return true;
}

bool qbarelem::pow(const qbarelem & a, const fmpz_t b)
{
//std::cout << "pow called a: " << a.tostring() << std::endl;
//printf("b: "); fmpz_print(b); printf("\n");

    if (fmpz_fits_si(b))
        return pow(a, fmpz_get_si(b), 1);

    if (fmpz_poly_length(a.minpoly.data) > 2)
        return false;

    if (!fmpz_is_one(a.minpoly.data->coeffs + 1))
        return false;

    if (fmpz_is_zero(a.minpoly.data->coeffs + 0))
    {
        set(slong(0));
        return fmpz_sgn(b) > 0;
    }

    if (!fmpz_is_pm1(a.minpoly.data->coeffs + 0))
        return false;

    slong r = 1;
    if (!fmpz_is_one(a.minpoly.data->coeffs + 0) && fmpz_is_odd(b))
        r = -1;
    set(r);
    return true;
}

bool qbarelem::pow(const qbarelem & a, const fmpq_t b)
{
    if (fmpz_fits_si(fmpq_numref(b)) && fmpz_fits_si(fmpq_denref(b)))
        return pow(a, fmpz_get_si(fmpq_numref(b)), fmpz_get_si(fmpq_denref(b)));

    assert(0 && "don't want to do this shit anymore");
    return false;
}



void irred_fmpz_poly_roots(std::vector<qbarelem> & v, fmpz_poly_t f)
{
    slong n = f->length - 1;
    v.resize(n);

    acb_struct * roots = (acb_struct *) malloc(n*sizeof(acb_struct));
    mag_struct * eps = (mag_struct *) malloc(n*sizeof(mag_struct));
    for (slong i = 0; i < n; i++)
    {
        acb_init(roots + i);
        mag_init(eps + i);
        fmpz_poly_set(v[i].minpoly.data, f);
        _fmpz_poly_primitive_part(v[i].minpoly.data->coeffs,
                                  v[i].minpoly.data->coeffs, n + 1);
    }

    slong p = FLINT_BITS;

try_again:

    p += (FLINT_BITS + p)/16;

    arb_fmpz_poly_complex_roots(roots, v[0].minpoly.data, 0, p);

    for (slong i = 0; i < n; i++)
    {
        acb_canonicalize(roots + i, p);
        if (!newton_test(v[0].minpoly.data, roots + i, eps + i))
            goto try_again;
    }

    for (slong i = 0; i < n; i++)
    {
        newton_trim(v[i].minpoly.data, roots + i, eps + i);
        acb_swap(v[i].location.data, roots + i);
        mag_swap(v[i].epsilon.data, eps + i);
        acb_clear(roots + i);
        mag_clear(eps + i);
    }
    free(roots);
    free(eps);
}


/* given a function for refining our location, set us to the right root of p */
template <typename RefineFunc>
bool qbarelem::select_root(fmpz_poly_t p, bool must_factor, RefineFunc refine)
{
//std::cout << "select_root called" << std::endl;
    slong prec;
    xacb_t t;
    xfmpz_poly_factor_t f;
//timeit_t timer;

    if (must_factor)
    {
//timeit_start(timer);
        fmpz_poly_factor(f.data, p);
//timeit_stop(timer);
//flint_printf("factor time: %wd\n", timer->wall);
//std::cout << "done factoring" << std::endl;
    }
    else
    {
        _fmpz_poly_primitive_part(p->coeffs, p->coeffs, p->length);
        fmpz_poly_factor_fit_length(f.data, 1);
        f.data->num = 1;
        fmpz_poly_swap(f.data->p + 0, p);
    }

//timeit_start(timer);

try_again:

//std::cout << "before location: " << location.tostring() << std::endl;

    refine(location.data, f.data->num > 1);

    prec = unary_wp(location.data);

//std::cout << "trying location: " << location.tostring() << std::endl;
//SleepMS(1000);

    if (f.data->num > 1)
    {
        for (slong i = 0; i < f.data->num; i++)
        {
            arb_fmpz_poly_evaluate_acb(t.data, f.data->p + i, location.data, prec);
            if (!acb_contains_zero(t.data))
            {
                fmpz_poly_swap(f.data->p + i, f.data->p + f.data->num - 1);
                f.data->num--;
                i--;
            }
        }

        if (f.data->num > 1)
            goto try_again;
    }

    assert(f.data->num == 1);

    if (!newton_test(f.data->p + 0, location.data, epsilon.data))
        goto try_again;

    fmpz_poly_swap(minpoly.data, f.data->p + 0);

//timeit_stop(timer);
//flint_printf("select time: %wd\n", timer->wall);

//std::cout << "select_root returning" << std::endl;

    return true;
}


/************** (m11*a + m12)/(m21*a + m22) *********/

class moebius_refiner {
    xacb_t a, ta;
    xmag_t ae;
    const fmpz * m11, * m12, * m21, * m22;
    const fmpz_poly_struct * aminpoly;
    bool first;

public:

    moebius_refiner(const qbarelem & A, const fmpz_t M11, const fmpz_t M12,
                                        const fmpz_t M21, const fmpz_t M22)
    {
        m11 = M11;
        m12 = M12;
        m22 = M22;
        m21 = M21;
        first = true;
        aminpoly = A.minpoly.data;
        acb_set(a.data, A.location.data);
        mag_set(ae.data, A.epsilon.data);
    }
    void operator() (acb_t x, bool need_containment)
    {
        xacb_t n, d;

    try_again:

        if (!first)
        {
            newton_step(aminpoly, a.data, ae.data);
        }
        first = false;

        slong p = unary_wp(a.data);

        if (need_containment)
            acb_set_radius(x, a.data, ae.data);
        else
            acb_set(x, a.data);

        acb_mul_fmpz(n.data, x, m11, p);
        acb_add_fmpz(n.data, n.data, m12, p);
        acb_mul_fmpz(d.data, x, m21, p);
        acb_add_fmpz(d.data, d.data, m22, p);
        acb_div(x, n.data, d.data, p);

        if (!acb_is_finite(x))
            goto try_again;
    }
};


bool qbarelem::apply_moebius(const qbarelem & a, const fmpz_t m11, const fmpz_t m12,
                                                 const fmpz_t m21, const fmpz_t m22)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "moebius " << a.tostring() << std::endl;
//printf("m11: "); fmpz_print(m11); printf("\n");
//printf("m12: "); fmpz_print(m12); printf("\n");
//printf("m21: "); fmpz_print(m21); printf("\n");
//printf("m22: "); fmpz_print(m22); printf("\n");

    xfmpz_t n, d;

    if (fmpz_poly_length(a.minpoly.data) <= 2)
    {
        fmpz_mul(n.data, a.minpoly.data->coeffs + 0, m11);
        fmpz_mul(d.data, a.minpoly.data->coeffs + 0, m21);
        fmpz_submul(n.data, a.minpoly.data->coeffs + 1, m12);
        fmpz_submul(d.data, a.minpoly.data->coeffs + 1, m22);
        return div(n.data, d.data);
    }

    fmpz_mul(d.data, m11, m22);
    fmpz_submul(d.data, m12, m21);

    if (fmpz_is_zero(d.data))
        return fmpz_is_zero(m22) ? div(m11, m21) : div(m12, m22);

    xfmpz_poly_t t, u, v, U, V, W;
    slong da = fmpz_poly_degree(a.minpoly.data);

    fmpz_poly_set_coeff_fmpz(u.data, 1, m22);
    fmpz_poly_set_coeff_fmpz(u.data, 0, m12);
    fmpz_neg(u.data->coeffs + 0, u.data->coeffs + 0);

    fmpz_poly_set_coeff_fmpz(v.data, 1, m21);
    fmpz_poly_set_coeff_fmpz(v.data, 0, m11);
    fmpz_neg(v.data->coeffs + 0, v.data->coeffs + 0);
    fmpz_poly_neg(v.data, v.data);

    for (slong i = 0; i <= da; i++)
    {
        fmpz_poly_pow(U.data, u.data, i);
        fmpz_poly_pow(V.data, v.data, da - i);
        fmpz_poly_mul(W.data, U.data, V.data);
        fmpz_poly_scalar_addmul_fmpz(t.data, W.data, a.minpoly.data->coeffs + i);
    }

    select_root(t.data, false, moebius_refiner(a, m11, m12, m21, m22));

    return true;
}

bool qbarelem::set_expipi(const fmpq_t a)
{
    if (!fmpz_fits_si(fmpq_denref(a)))
        return false;

    ulong n = fmpz_get_si(fmpq_denref(a));
    fmpz_poly_cyclotomic(minpoly.data, fmpz_is_even(fmpq_numref(a)) ? n : 2*n);

    slong p = 10;

try_again:

    p += FLINT_BITS;

    arb_sin_cos_pi_fmpq(acb_imagref(location.data), acb_realref(location.data), a, p);
    acb_canonicalize(location.data, p);

    if (!newton_test(minpoly.data, location.data, epsilon.data))
        goto try_again;

    return true;
}

bool qbarelem::mul_expipi(const qbarelem & a, const fmpq_t b)
{
    qbarelem t;
    return t.set_expipi(b) && mul(a, t);
}

/************* a + b ***************/

class add_refiner {
    xacb_t a, b, ta, tb;
    xmag_t ae, be;
    const fmpz_poly_struct * aminpoly, * bminpoly;
    bool first;

public:

    add_refiner(const qbarelem & A, const qbarelem & B)
    {
        first = true;
        aminpoly = A.minpoly.data;
        bminpoly = B.minpoly.data;
        acb_set(a.data, A.location.data);
        acb_set(b.data, B.location.data);
        mag_set(ae.data, A.epsilon.data);
        mag_set(be.data, B.epsilon.data);
    }
    void operator() (acb_t x, bool need_containment)
    {
        if (!first)
        {
            newton_step(aminpoly, a.data, ae.data);
            newton_step(bminpoly, b.data, be.data);
        }

        first = false;

        slong p = binary_wp(a.data, b.data);

        if (need_containment)
        {
            acb_set_radius(ta.data, a.data, ae.data);
            acb_set_radius(tb.data, b.data, be.data);
            acb_add(x, ta.data, tb.data, p);
        }
        else
        {
            acb_add(x, a.data, b.data, p);
        }
    }
};

bool qbarelem::add(const qbarelem & a, const qbarelem & b)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "add " << a.tostring() << std::endl;
//std::cout << "  + " << b.tostring() << std::endl;

    if (fmpz_poly_length(a.minpoly.data) <= 2)
    {
        xfmpz_t m12(a.minpoly.data->coeffs + 0);
        fmpz_neg(m12.data, m12.data);
        return apply_moebius(b, a.minpoly.data->coeffs + 1, m12.data,
                                eget_cint_data(0), a.minpoly.data->coeffs + 1);
    }

    if (fmpz_poly_length(b.minpoly.data) <= 2)
    {
        xfmpz_t m12(b.minpoly.data->coeffs + 0);
        fmpz_neg(m12.data, m12.data);
        return apply_moebius(a, b.minpoly.data->coeffs + 1, m12.data,
                                eget_cint_data(0), b.minpoly.data->coeffs + 1);
    }

//timeit_t timer;
//timeit_start(timer);

    xfmpz_poly_t t;

    if (algo == Algo_pgcd)
    {
		rfmpz_poly_t r;
        sparse_poly<xfmpz_poly_t> A, B;
        set_univar<rfmpz_poly_t>(B, a.minpoly, r);
        set_univar_shift<rfmpz_poly_t>(A, b.minpoly, -1, r);
        resultant<rfmpz_poly_t>(t, A, B, r);
    }
    else
    {
#if 0
        x_fmpz_vector_t xs, ys;

        slong ad = fmpz_poly_degree(a.minpoly.data);
        slong bd = fmpz_poly_degree(b.minpoly.data);
        slong rd = ad*bd;

        _fmpz_vector_fit_length(xs.data, rd + 1);
        _fmpz_vector_fit_length(ys.data, rd + 1);

        fmpz_poly_set(t.data, b.minpoly.data);
        for (slong i = 1; i <= bd; i += 2)
            fmpz_neg(t.data->coeffs + i, t.data->coeffs + i);

        for (slong i = 0; i <= rd; i++)
        {
            fmpz_set_si(xs.data->array + i, i);
            fmpz_poly_resultant(ys.data->array + i, a.minpoly.data, t.data);
            if (i < rd)
                _fmpz_poly_taylor_shift(t.data->coeffs, eget_cint_data(-1), bd + 1);
        }
        fmpz_poly_interpolate_fmpz_vec(t.data, xs.data->array, ys.data->array, rd + 1);
#endif
    }

//timeit_stop(timer);
//flint_printf("poly + time: %wd\n", timer->wall);

    select_root(t.data, true, add_refiner(a, b));

//std::cout << "add return: " << tostring() << std::endl;
    return true;
}



/************* a - b ***************/

class sub_refiner {
    xacb_t a, b, ta, tb;
    xmag_t ae, be;
    const fmpz_poly_struct * aminpoly, * bminpoly;
    bool first;

public:

    sub_refiner(const qbarelem & A, const qbarelem & B)
    {
        first = true;
        aminpoly = A.minpoly.data;
        bminpoly = B.minpoly.data;
        acb_set(a.data, A.location.data);
        acb_set(b.data, B.location.data);
        mag_set(ae.data, A.epsilon.data);
        mag_set(be.data, B.epsilon.data);
    }

    void operator() (acb_t x, bool need_containment)
    {
        if (!first)
        {
            newton_step(aminpoly, a.data, ae.data);
            newton_step(bminpoly, b.data, be.data);
        }

        first = false;

        slong p = binary_wp(a.data, b.data);

        if (need_containment)
        {
            acb_set_radius(ta.data, a.data, ae.data);
            acb_set_radius(tb.data, b.data, be.data);
            acb_sub(x, ta.data, tb.data, p);
        }
        else
        {
            acb_sub(x, a.data, b.data, p);
        }
    }
};


void qbarelem::neg(const qbarelem & a)
{
    if (this != &a)
        set(a);

    for (slong i = minpoly.data->length - 2; i >= 0; i -= 2)
        fmpz_neg(minpoly.data->coeffs + i, minpoly.data->coeffs + i);

    acb_neg(location.data, location.data);
}

bool qbarelem::sub(const qbarelem & a, const qbarelem & b)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "sub " << a.tostring() << " - " << b.tostring() << std::endl;

    if (fmpz_poly_length(a.minpoly.data) <= 2)
    {
        if (!apply_moebius(b, a.minpoly.data->coeffs + 1, a.minpoly.data->coeffs + 0,
                              eget_cint_data(0),          a.minpoly.data->coeffs + 1))
        {
            return false;
        }
        neg(*this);
        return true;
    }

    if (fmpz_poly_length(b.minpoly.data) <= 2)
    {
        return apply_moebius(a, b.minpoly.data->coeffs + 1, b.minpoly.data->coeffs + 0,
                                eget_cint_data(0),          b.minpoly.data->coeffs + 1);
    }

    xfmpz_poly_t t;
	rfmpz_poly_t r;
    sparse_poly<xfmpz_poly_t> A, B;
    set_univar<rfmpz_poly_t>(A, a.minpoly, r);
    set_univar_shift<rfmpz_poly_t>(B, b.minpoly, +1, r);
    resultant<rfmpz_poly_t>(t, A, B, r);

    select_root(t.data, true, sub_refiner(a, b));

//std::cout << "sub return: " << tostring() << std::endl;
    return true;
}



/************* Re[a] ***************/

bool qbarelem::isreal() const
{
    if (minpoly.data->length <= 2 || arb_is_zero(acb_imagref(location.data)))
        return true;

    /* t will always contain a unique root and TODO: can be a shallow object */
    xacb_t a, t;
    xmag_t ae, te;

    acb_set_radius(t.data, location.data, epsilon.data);
    if (!arb_contains_zero(acb_imagref(t.data)))
        return false;

    acb_set(a.data, location.data);
    mag_set(ae.data, epsilon.data);

try_again:

    /* t intersects the real axis */
    arb_zero(acb_imagref(t.data));
    if (newton_test(minpoly.data, t.data, te.data))
        return true;

    newton_step(minpoly.data, a.data, ae.data);

    acb_set_radius(t.data, a.data, ae.data);
    if (!arb_contains_zero(acb_imagref(t.data)))
        return false;

    goto try_again;
}

class realpart_refiner {
    xacb_t a, ta;
    xmag_t ae;
    const fmpz_poly_struct * aminpoly;
    bool first;

public:

    realpart_refiner(const qbarelem & A)
    {
        first = true;
        aminpoly = A.minpoly.data;
        acb_set(a.data, A.location.data);
        mag_set(ae.data, A.epsilon.data);
    }
    void operator() (acb_t x, bool need_containment)
    {
        if (!first)
            newton_step(aminpoly, a.data, ae.data);

        first = false;

        if (need_containment)
            acb_set_radius(x, a.data, ae.data);
        else
            acb_set(x, a.data);

        acb_chop_imag(x);
        acb_mul_2exp_si(x, x, 1);
    }
};

bool qbarelem::realpart(const qbarelem & a)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "realpart " << a.tostring() << std::endl;

    if (a.isreal())
    {
        set(a);
        return true;
    }

//timeit_t timer;
//timeit_start(timer);

    slong ad = fmpz_poly_degree(a.minpoly.data);
    slong td;
    if (mul_si_checked(td, ad, ad - 1))
        throw exception_arithmetic(17);

    td = td/2;

    xfmpz_poly_t t1;

    if (algo == Algo_pgcd)
    {
        xfmpz_poly_series_t t2;

        try {
            slong prec = td + ad;
		    rfmpz_poly_series_t r(prec);
            sparse_poly<xfmpz_poly_series_t> A, B;

            while (true)
            {
                set_univar<rfmpz_poly_series_t>(A, a.minpoly, r);
                set_univar_shift<rfmpz_poly_series_t>(B, a.minpoly, -1, r);
                resultant<rfmpz_poly_series_t>(t2, A, B, r);
                if (t2.absolute_prec() > td)
                    break;
                prec += td - t2.absolute_prec() + 1;
                r.set_prec(prec);
            }
        }
        catch (exception_arithmetic & e)
        {
	        rfmpz_poly_t r;
            sparse_poly<xfmpz_poly_t> A, B;
            set_univar<rfmpz_poly_t>(A, a.minpoly, r);
            set_univar_shift<rfmpz_poly_t>(B, a.minpoly, -1, r);
            resultant<rfmpz_poly_t>(t2.terms, A, B, r);
            t2.start = 0;
            t2.prec = 2*td + ad + 1;
            t2.normalize();
        }

        // divide the terms of t2 by a(x/2)*2^deg(a)
        fmpz_poly_scalar_mul_fmpz(t1.data, a.minpoly.data, a.minpoly.data->coeffs + ad);
        for (slong i = 1; i <= ad; i += 1)
            fmpz_mul_2exp(t1.data->coeffs + ad - i, t1.data->coeffs + ad - i, i);
        fmpz_poly_div_series(t1.data, t2.terms.data, t1.data, td + 1);

        assert((t2.start % 2) == 0);
		fmpz_poly_sqrt_series(t1.data, t1.data, td - t2.start/2 + 1);
        fmpz_poly_shift_left(t1.data, t1.data, t2.start/2);
    }
    else
    {
#if 0
        x_fmpz_vector_t xs, ys;
        xfmpz_poly_t t2, a2;
        xfmpz_t u, res;

        _fmpz_vector_fit_length(xs.data, 2*td + 1);
        _fmpz_vector_fit_length(ys.data, 2*td + 1);

        fmpz_poly_set(t1.data, a.minpoly.data);
        fmpz_poly_scalar_mul_fmpz(a2.data, a.minpoly.data, a.minpoly.data->coeffs + ad);
        for (slong i = 1; i <= ad; i += 1)
        {
            fmpz_mul_2exp(a2.data->coeffs + ad - i, a2.data->coeffs + ad - i, i);
            if (i & 1)
                fmpz_neg(t1.data->coeffs + i, t1.data->coeffs + i);
        }
        fmpz_poly_set(t2.data, t1.data);

        fmpz_zero(xs.data->array + 0);
        fmpz_poly_resultant(res.data, a.minpoly.data, t1.data);
        fmpz_divexact(ys.data->array + 0, res.data, a2.data->coeffs + 0);

        for (slong i = 1; i <= td; i++)
        {
            _fmpz_poly_taylor_shift(t1.data->coeffs, eget_cint_data(-1), ad + 1);
            fmpz_set_si(xs.data->array + 2*i - 1, i);
            fmpz_poly_resultant(res.data, a.minpoly.data, t1.data);
            fmpz_poly_evaluate_fmpz(u.data, a2.data, xs.data->array + 2*i - 1);
            fmpz_divexact(ys.data->array + 2*i - 1, res.data, u.data);

            _fmpz_poly_taylor_shift(t2.data->coeffs, eget_cint_data(+1), ad + 1);
            fmpz_set_si(xs.data->array + 2*i - 0, -i);
            fmpz_poly_resultant(res.data, a.minpoly.data, t2.data);
            fmpz_poly_evaluate_fmpz(u.data, a2.data, xs.data->array + 2*i - 0);
            fmpz_divexact(ys.data->array + 2*i - 0, res.data, u.data);
        }

        // t1 is a square TODO how to get its sqrt better
        fmpz_poly_interpolate_fmpz_vec(t1.data, xs.data->array, ys.data->array, 2*td + 1);
#endif
    }

//timeit_stop(timer);
//flint_printf("poly re time: %wd\n", timer->wall);

    select_root(t1.data, true, realpart_refiner(a));

    apply_moebius(*this, eget_cint_data(1), eget_cint_data(0),
                         eget_cint_data(0), eget_cint_data(2));

//std::cout << "realpart return: " << tostring() << std::endl;

    return true;
}




/************* a*b ***************/

class mul_refiner {
    xacb_t a, b, ta, tb;
    xmag_t ae, be;
    const fmpz_poly_struct * aminpoly, * bminpoly;
    bool first;

public:

    mul_refiner(const qbarelem & A, const qbarelem & B)
    {
        first = true;
        aminpoly = A.minpoly.data;
        bminpoly = B.minpoly.data;
        acb_set(a.data, A.location.data);
        acb_set(b.data, B.location.data);
        mag_set(ae.data, A.epsilon.data);
        mag_set(be.data, B.epsilon.data);
    }

    void operator() (acb_t x, bool need_containment)
    {
        if (!first)
        {
            newton_step(aminpoly, a.data, ae.data);
            newton_step(bminpoly, b.data, be.data);
        }

        first = false;

        slong p = binary_wp(a.data, b.data);

        if (need_containment)
        {
            acb_set_radius(ta.data, a.data, ae.data);
            acb_set_radius(tb.data, b.data, be.data);
            acb_mul(x, ta.data, tb.data, p);
        }
        else
        {
            acb_mul(x, a.data, b.data, p);
        }
    }
};


bool qbarelem::mul(const qbarelem & a, const qbarelem & b)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "mul " << a.tostring() << std::endl;
//std::cout << "  * " << b.tostring() << std::endl;

    if (fmpz_poly_length(a.minpoly.data) <= 2)
    {
        xfmpz_t m11(a.minpoly.data->coeffs + 0);
        fmpz_neg(m11.data, m11.data);
        return apply_moebius(b, m11.data, eget_cint_data(0),
                                eget_cint_data(0), a.minpoly.data->coeffs + 1);
    }

    if (fmpz_poly_length(b.minpoly.data) <= 2)
    {
        xfmpz_t m11(b.minpoly.data->coeffs + 0);
        fmpz_neg(m11.data, m11.data);
        return apply_moebius(a, m11.data, eget_cint_data(0),
                                eget_cint_data(0), b.minpoly.data->coeffs + 1);
    }

    xfmpz_poly_t t;

//timeit_t timer;
//timeit_start(timer);

    if (algo == Algo_pgcd)
    {
		rfmpz_poly_t r;
        sparse_poly<xfmpz_poly_t> A, B;
        set_univar<rfmpz_poly_t>(A, a.minpoly, r);
        set_univar_scale(B, b.minpoly, -1);
        resultant<rfmpz_poly_t>(t, A, B, r);
    }
    else
    {
#if 0
        x_fmpz_vector_t xs, ys;
        xfmpz_t acc;
        slong ad = fmpz_poly_degree(a.minpoly.data);
        slong bd = fmpz_poly_degree(b.minpoly.data);
        slong rd = ad*bd;

        _fmpz_vector_fit_length(xs.data, rd + 1);
        _fmpz_vector_fit_length(ys.data, rd + 1);
        fmpz_poly_fit_length(t.data, bd + 1);
        _fmpz_poly_set_length(t.data, bd + 1);

        for (slong i = 0; i <= rd; i++)
        {
            slong scale = (i & 1) ? -(i/2 + 1) : i/2 + 1;
            fmpz_set_si(xs.data->array + i, scale);
            fmpz_one(acc.data);
            for (slong j = 0; j <= bd; j++)
            {
                fmpz_mul(t.data->coeffs + bd - j, b.minpoly.data->coeffs + j, acc.data);
                fmpz_mul_si(acc.data, acc.data, scale);
            }
            fmpz_poly_resultant(ys.data->array + i, a.minpoly.data, t.data);
        }

        fmpz_poly_interpolate_fmpz_vec(t.data, xs.data->array, ys.data->array, rd + 1);
#endif
    }

//timeit_stop(timer);
//flint_printf("poly * time: %wd\n", timer->wall);

    select_root(t.data, true, mul_refiner(a, b));

//std::cout << "mul return: " << tostring() << std::endl;
    return true;
}


/************* a/b ***************/

class div_refiner {
    xacb_t a, b, ta, tb;
    xmag_t ae, be;

    const fmpz_poly_struct * aminpoly, * bminpoly;
    bool first;

public:

    div_refiner(const qbarelem & A, const qbarelem & B)
    {
        first = true;
        aminpoly = A.minpoly.data;
        bminpoly = B.minpoly.data;
        acb_set(a.data, A.location.data);
        acb_set(b.data, B.location.data);
        mag_set(ae.data, A.epsilon.data);
        mag_set(be.data, B.epsilon.data);
    }
    void operator() (acb_t x, bool need_containment)
    {
    try_again:

        if (!first)
        {
            newton_step(aminpoly, a.data, ae.data);
            newton_step(bminpoly, b.data, be.data);
        }

        first = false;

        slong p = binary_wp(a.data, b.data);

        if (need_containment)
        {
            acb_set_radius(ta.data, a.data, ae.data);
            acb_set_radius(tb.data, b.data, be.data);
            acb_div(x, ta.data, tb.data, p);
        }
        else
        {
            acb_div(x, a.data, b.data, p);
        }

        if (!acb_is_finite(x))
            goto try_again;
    }
};


bool qbarelem::div(const qbarelem & a, const qbarelem & b)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "mul " << a.tostring() << " * " << b.tostring() << std::endl;

    if (fmpz_poly_length(a.minpoly.data) <= 2)
    {
        xfmpz_t m11(a.minpoly.data->coeffs + 0);
        fmpz_neg(m11.data, m11.data);
        return apply_moebius(b, eget_cint_data(0),          m11.data,
                                a.minpoly.data->coeffs + 1, eget_cint_data(0));
    }

    if (fmpz_poly_length(b.minpoly.data) <= 2)
    {
        xfmpz_t m11(b.minpoly.data->coeffs + 0);
        fmpz_neg(m11.data, m11.data);
        return apply_moebius(a, b.minpoly.data->coeffs + 1, eget_cint_data(0),
                                eget_cint_data(0),          m11.data);
    }

    xfmpz_poly_t t;
	rfmpz_poly_t r;
    sparse_poly<xfmpz_poly_t> A, B;
    set_univar<rfmpz_poly_t>(A, a.minpoly, r);
    set_univar_scale(B, b.minpoly, +1);
    resultant<rfmpz_poly_t>(t, A, B, r);

    select_root(t.data, true, div_refiner(a, b));

//std::cout << "div return: " << tostring() << std::endl;
    return true;
}


/************ a^(b/c) **************/
class pow_refiner {
    xacb_t a;
    xmag_t ae;
    const fmpz_poly_struct * aminpoly;
    slong b;
    ulong c;
    xfmpq_t pow;
    bool first;

public:

    pow_refiner(const qbarelem & A, slong B, ulong C)
    {
        b = B;
        c = C;
        first = true;
        aminpoly = A.minpoly.data;
        fmpz_set_si(fmpq_numref(pow.data), B);
        fmpz_set_ui(fmpq_denref(pow.data), C);
        acb_set(a.data, A.location.data);
        mag_set(ae.data, A.epsilon.data);
    }

    void operator() (acb_t x, bool need_containment)
    {
        xacb_t ta, s;
        bool ok;

    try_again:

        if (!first)
            newton_step(aminpoly, a.data, ae.data);

        first = false;

        slong p = unary_wp(a.data);

        /* ta is our ball containing a unique root of aminpoly */
        acb_set_radius(ta.data, a.data, ae.data);

        if (c < 2)
        {
            /* no shenanigans because a -> a^b is continuous */
            if (need_containment)
                acb_pow_si(x, ta.data, b, p);
            else
                acb_pow_si(x, a.data, b, p);
            return;
        }

        if (acb_contains_zero(ta.data))
            goto try_again;

        if (arb_is_zero(acb_imagref(ta.data)))
        {
    try_real:
            assert(arb_is_zero(acb_imagref(ta.data)));

            if (arb_is_positive(acb_realref(ta.data)))
            {
                arb_pow_fmpq(acb_realref(x), acb_realref(ta.data), pow.data, p);
                arb_zero(acb_imagref(x));
            }
            else if (arb_is_negative(acb_realref(ta.data)))
            {
                arb_sin_cos_pi_fmpq(acb_imagref(s.data), acb_realref(s.data), pow.data, p);
                arb_abs(acb_realref(ta.data), acb_realref(ta.data));
                arb_pow_fmpq(acb_realref(ta.data), acb_realref(ta.data), pow.data, p);
                arb_mul(acb_realref(x), acb_realref(s.data), acb_realref(ta.data), p);
                arb_mul(acb_imagref(x), acb_imagref(s.data), acb_realref(ta.data), p);
            }
            else
            {
                /* hmm, ta was supposed to not contain zero */
                goto try_again;
            }
        }
        else if (!arb_contains_zero(acb_imagref(ta.data)) || arb_is_positive(acb_realref(ta.data)))
        {
            /* if this is not good enough this time, eventually it will be */
            acb_root_ui(s.data, ta.data, c, p);
            acb_pow_si(x, s.data, b, p);
        }
        else
        {
            /* intersect ta with real axis and see if this real ball is ok */
            acb_chop_imag(ta.data);
            if (newton_test(aminpoly, ta.data, arb_radref(acb_realref(s.data))))
                goto try_real;
            else
                goto try_again;
        }

        /* may have had trouble if relative precision is too low */
        if (!acb_is_finite(x))
            goto try_again;
    }
};

bool qbarelem::pow(const qbarelem & a, slong b, ulong c)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "pow: " << a.tostring() << "^(" << b << "/" << c << ")" << std::endl;

    assert(c > 0);
    if (fmpz_is_zero(a.minpoly.data->coeffs + 0))
    {
        set(slong(0));
        return b > 0;
    }
    if (b == 0)
    {
        set(slong(1));
        return true;
    }

//timeit_t timer;
//timeit_start(timer);

    xfmpz_poly_t t;
    slong ad = fmpz_poly_degree(a.minpoly.data);
    slong rd = ad*c;
    slong absb = std::abs(b);

    if (b == 1)
    {
        fmpz_poly_inflate(t.data, a.minpoly.data, c);
    }
    else if (b == -1)
    {
        fmpz_poly_reverse(t.data, a.minpoly.data, ad + 1);
        fmpz_poly_inflate(t.data, t.data, c);
    }
    else if (algo == Algo_pgcd)
    {
		rfmpz_poly_t r;
        sparse_poly<xfmpz_poly_t> P, Q;
        Q.fit_length(2);
        Q.length = 2;
        Q.exps[0] = absb;
        Q.exps[1] = 0;
        fmpz_poly_one(Q.coeffs[0].data);
        fmpz_poly_neg(Q.coeffs[1].data, Q.coeffs[0].data);
        fmpz_poly_shift_left(Q.coeffs[b > 0].data, Q.coeffs[b > 0].data, c);
        set_univar<rfmpz_poly_t>(P, a.minpoly, r);
        resultant<rfmpz_poly_t>(t, P, Q, r);
    }
    else
    {
#if 0
        x_fmpz_vector_t xs, ys;
        _fmpz_vector_fit_length(xs.data, rd + 1);
        _fmpz_vector_fit_length(ys.data, rd + 1);
        fmpz_poly_fit_length(t.data, absb + 1);
        _fmpz_poly_set_length(t.data, absb + 1);
        for (slong i = 0; i <= absb; i++)
            fmpz_zero(t.data->coeffs + i);

        for (slong i = 0; i <= rd; i++)
        {
            fmpz_set_si(xs.data->array + i, i + 1);
            if (b > 0)
            {
                fmpz_set_si(t.data->coeffs + absb, -1);
                fmpz_pow_ui(t.data->coeffs + 0, xs.data->array + i, c);
            }
            else
            {
                fmpz_set_si(t.data->coeffs + 0, -1);
                fmpz_pow_ui(t.data->coeffs + absb, xs.data->array + i, c);
            }
            fmpz_poly_resultant(ys.data->array + i, a.minpoly.data, t.data);
        }

        fmpz_poly_interpolate_fmpz_vec(t.data, xs.data->array, ys.data->array, rd + 1);
#endif
    }

//timeit_stop(timer);
//flint_printf("poly ^ time: %wd\n", timer->wall);


    select_root(t.data, c > 1 || absb > 1, pow_refiner(a, b, c));

//std::cout << "pow return: " << tostring() << std::endl;
    return true;
}

void qbarelem::quadraticroot(const fmpz_t a, const fmpz_t b, const fmpz_t c)
{
    xfmpz_t d, r, s;

//std::cout << "quadratic root" << std::endl;
//printf("a: "); fmpz_print(a); printf("\n");
//printf("b: "); fmpz_print(b); printf("\n");
//printf("c: "); fmpz_print(c); printf("\n");

    fmpz_mul(d.data, b, b);
    fmpz_mul(s.data, a, c);
    fmpz_submul_ui(d.data, s.data, 4);
    if (fmpz_sgn(d.data) >= 0)
    {
        fmpz_sqrtrem(s.data, r.data, d.data);
        if (fmpz_is_zero(r.data))
        {
            fmpz_sub(s.data, s.data, b);
            fmpz_mul_2exp(r.data, a, 1);
            div(s.data, r.data);
            return;
        }
    }

    fmpz_poly_fit_length(minpoly.data, 3);
    fmpz_set(minpoly.data->coeffs + 0, c);
    fmpz_set(minpoly.data->coeffs + 1, b);
    fmpz_set(minpoly.data->coeffs + 2, a);
    _fmpz_poly_set_length(minpoly.data, 3);
    _fmpz_poly_primitive_part(minpoly.data->coeffs, minpoly.data->coeffs, 3);

    flint_bitcnt_t d_bits = fmpz_bits(d.data)/2;
    flint_bitcnt_t a_bits = fmpz_bits(a);
    flint_bitcnt_t b_bits = fmpz_bits(b);
    b_bits = std::max(b_bits, d_bits);
    if (b_bits > a_bits + d_bits)
        b_bits -= a_bits + d_bits;

    slong p = b_bits > a_bits + d_bits ? b_bits - (a_bits + d_bits) : b_bits;

try_again:

    p += FLINT_BITS;

    if (fmpz_sgn(d.data) >= 0)
    {
        arb_sqrt_fmpz(acb_realref(location.data), d.data, p);
        arb_sub_fmpz(acb_realref(location.data), acb_realref(location.data), b, p);
        arb_div_fmpz(acb_realref(location.data), acb_realref(location.data), a, p);
        arb_mul_2exp_si(acb_realref(location.data), acb_realref(location.data), -1);
        arb_zero(acb_imagref(location.data));
    }
    else
    {
        fmpz_neg(d.data, d.data);
        arb_sqrt_fmpz(acb_imagref(location.data), d.data, p);
        arb_div_fmpz(acb_imagref(location.data), acb_imagref(location.data), a, p);
        arb_fmpz_div_fmpz(acb_realref(location.data), b, a, p);
        arb_neg(acb_realref(location.data), acb_realref(location.data));
        acb_mul_2exp_si(location.data, location.data, -1);        
    }

    if (!newton_test(minpoly.data, location.data, epsilon.data))
        goto try_again;
}



/********* a +- I*Sqrt[1 - a^2] **************/

class circle_root_refiner {
    xacb_t a;
    xmag_t ae;
    const fmpz_poly_struct * aminpoly;
    int sign;
    bool first;

public:

    circle_root_refiner(const qbarelem & A, int sign_)
    {
        sign = sign_;
        first = true;
        aminpoly = A.minpoly.data;
        acb_set(a.data, A.location.data);
        mag_set(ae.data, A.epsilon.data);
    }

    void operator() (acb_t x, bool need_containment)
    {
        xacb_t ta, s;
        bool ok;

    try_again:

        if (!first)
            newton_step(aminpoly, a.data, ae.data);

        first = false;

        slong p = unary_wp(a.data);

        /* ta is our ball containing a unique root of aminpoly */
        acb_set_radius(ta.data, a.data, ae.data);

        /* nasty shit because the damn function has a discontinuity along Abs[Re[a]] > 1 */

        if (arb_is_zero(acb_imagref(ta.data)))
        {
    try_real:
            assert(arb_is_zero(acb_imagref(ta.data)));

            arb_sqr(acb_realref(s.data), acb_realref(ta.data), p);
            arb_sub_ui(acb_realref(s.data), acb_realref(s.data), 1, p);
            arb_neg(acb_realref(s.data), acb_realref(s.data));

            if (arb_is_positive(acb_realref(s.data)))
            {
                arb_sqrt(acb_imagref(s.data), acb_realref(s.data), p);
                arb_zero(acb_realref(s.data));
            }
            else if (arb_is_negative(acb_realref(s.data)))
            {
                arb_neg(acb_realref(s.data), acb_realref(s.data));
                arb_sqrt(acb_realref(s.data), acb_realref(s.data), p);
                arb_neg(acb_realref(s.data), acb_realref(s.data));
                arb_zero(acb_imagref(s.data));
            }
            else
            {
                goto try_again;
            }
        }
        else if (!arb_contains_zero(acb_imagref(ta.data)))
        {
            /* if this is not good enough this time, eventually it will be */
            acb_sqr(s.data, ta.data, p);
            arb_sub_ui(acb_realref(s.data), acb_realref(s.data), 1, p);
            acb_neg(s.data, s.data);
            acb_sqrt(s.data, s.data, p);
            acb_mul_onei(s.data, s.data);
        }
        else
        {
            /* intersect ta with real axis and see if this real ball is ok */
            acb_chop_imag(ta.data);
            if (newton_test(aminpoly, ta.data, arb_radref(acb_realref(s.data))))
                goto try_real;
            else
                goto try_again;
        }

        if (sign > 0)
            acb_add(x, ta.data, s.data, p);
        else
            acb_sub(x, ta.data, s.data, p);

        /* may have had trouble if relative precision is too low */
        if (!acb_is_finite(x))
            goto try_again;
    }
};


bool qbarelem::circle_root(const qbarelem & a, int sign)
{
//std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
//std::cout << "circleroot: " << a.tostring() << std::endl;

    if (fmpz_poly_length(a.minpoly.data) <= 2)
    {
        xfmpz_t a0(a.minpoly.data->coeffs + 0);
        xfmpz_t a1(a.minpoly.data->coeffs + 1);
        fmpz_mul_2exp(a0.data, a0.data, 1);

        if (fmpz_cmpabs(a.minpoly.data->coeffs + 0, a.minpoly.data->coeffs + 1) >= 0)
        {
            fmpz_neg(a0.data, a0.data);
            fmpz_neg(a1.data, a1.data);
        }

        if (sign < 0)
        {
            fmpz_neg(a0.data, a0.data);
            fmpz_neg(a1.data, a1.data);
        }

        quadraticroot(a1.data, a0.data, a1.data);

//std::cout << "circle root done" << tostring() << std::endl;
        return true;
    }

    xfmpz_poly_t t;

	rfmpz_poly_t r;
    sparse_poly<xfmpz_poly_t> P, Q;
    Q.fit_length(2);
    Q.length = 2;
    Q.exps[0] = 1;
    Q.exps[1] = 0;
    fmpz_poly_zero(Q.coeffs[0].data);
    fmpz_poly_set_coeff_si(Q.coeffs[0].data, 1, -2);
    fmpz_poly_zero(Q.coeffs[1].data);
    fmpz_poly_set_coeff_si(Q.coeffs[1].data, 0, 1);
    fmpz_poly_set_coeff_si(Q.coeffs[1].data, 2, 1);
    set_univar<rfmpz_poly_t>(P, a.minpoly, r);
    resultant<rfmpz_poly_t>(t, P, Q, r);

    select_root(t.data, true, circle_root_refiner(a, sign));

//std::cout << "circleroot return: " << tostring() << std::endl;
    return true;
}




/********* apply rat fxn **************/

class apply_rat_fxn_refiner {
    xacb_t a;
    xmag_t ae;
    const fmpz_poly_struct * aminpoly;
    const fmpz * num; slong num_length;
    const fmpz * den; slong den_length;
    bool first;

public:

    apply_rat_fxn_refiner(const qbarelem & A, const fmpz * num_, slong num_length_,
                                              const fmpz * den_, slong den_length_)
    {
        first = true;
        aminpoly = A.minpoly.data;
        acb_set(a.data, A.location.data);
        mag_set(ae.data, A.epsilon.data);
        num = num_;
        den = den_;
        num_length = num_length_;
        den_length = den_length_;
    }

    void operator() (acb_t x, bool need_containment)
    {
        xacb_t ta, u, v;

    try_again:

        if (!first)
            newton_step(aminpoly, a.data, ae.data);

        first = false;

        slong p = unary_wp(a.data);

        if (need_containment)
        {
            acb_set_radius(ta.data, a.data, ae.data);
            _arb_fmpz_poly_evaluate_acb(u.data, num, num_length, ta.data, p);
            _arb_fmpz_poly_evaluate_acb(v.data, den, den_length, ta.data, p);
        }
        else
        {
            _arb_fmpz_poly_evaluate_acb(u.data, num, num_length, a.data, p);
            _arb_fmpz_poly_evaluate_acb(v.data, den, den_length, a.data, p);
        }

        acb_div(x, u.data, v.data, p);
        if (!acb_is_finite(x))
            goto try_again;
    }
};

void qbarelem::apply_rat_fxn(const qbarelem & a, const fmpz * num, slong num_length,
                                                 const fmpz * den, slong den_length)
{
    xfmpz_poly_t t;

//std::cout << "apply_rat_fxn called" << std::endl;

	rfmpz_poly_t r;
    sparse_poly<xfmpz_poly_t> A, B;

    set_univar<rfmpz_poly_t>(A, a.minpoly, r);

    slong max_length = std::max(num_length, den_length);

    B.fit_length(max_length);
    B.length = 0;
    for (slong i = max_length - 1; i >= 0; i--)
    {
        B.exps[B.length] = i;
        fmpz_poly_zero(B.coeffs[B.length].data);
        if (i < num_length)
            fmpz_poly_set_coeff_fmpz(B.coeffs[B.length].data, 0, num + i);
        fmpz_poly_neg(B.coeffs[B.length].data, B.coeffs[B.length].data);
        if (i < den_length)
            fmpz_poly_set_coeff_fmpz(B.coeffs[B.length].data, 1, den + i);
        B.length += !fmpz_poly_is_zero(B.coeffs[B.length].data);
    }

    resultant<rfmpz_poly_t>(t, A, B, r);

    select_root(t.data, true, apply_rat_fxn_refiner(a, num, num_length, den, den_length));
}



/***************************************************************/

bool eget_fmpq_poly(fmpq_poly_t q, er f, er v)
{
    poly p(1);
    wex vars(emake_node(gs.sym_sList.copy(), ecopy(v)));

    if (!ex_to_polynomial(p, f, vars.get()))
    {
        return false;
    }

    xfmpq_t ai;
    slong ei;

    fmpq_poly_zero(q);

    for (ulong i = 0; i < p.coeffs.size(); i++)
    {
        if (!eget_fmpq(ai.data, p.coeffs[i].get()) ||
            COEFF_IS_MPZ(*p.exps[i].data) ||
            (ei = *p.exps[i].data, ei < 0))
        {
            return false;
        }

        fmpq_poly_set_coeff_fmpq(q, ei, ai.data);
    }

    return true;
}

bool eget_arb(arb_t x, er e)
{
//std::cout << "eget_arb e: " << ex_tostring_full(e) << std::endl;

    if (eis_double(e))
    {
        arb_set_d(x, edouble_number(e));
        arb_limit_precision(x, 53);
        return true;
    }
    else if (eis_real(e))
    {
        arb_set(x, ereal_data(e));
        return true;
    }
    else
    {
        return false;
    }
}

bool eget_acb(acb_t x, er e)
{
//std::cout << "eget_acb e: " << ex_tostring_full(e) << std::endl;

    if (eis_cmplx(e))
    {
        if (!eget_arb(acb_imagref(x), ecmplx_imag(e)))
            return false;
        e = ecmplx_real(e);
    }
    else
    {
        arb_zero(acb_imagref(x));
    }
    return eget_arb(acb_realref(x), e);
}

ex acb_to_ex(const acb_t z)
{
    ex e;

    if (!arb_is_zero(acb_realref(z)))
        e = emake_real_copy(acb_realref(z));
    else
        e = emake_cint(0);

    if (!arb_is_zero(acb_imagref(z)))
        e = emake_cmplx(e, emake_real_copy(acb_imagref(z)));

    return e;
}


ex qbarelem::get_ex() const
{
    std::vector<wex> p;
    wex slot(emake_node(gs.sym_sSlot.copy(), emake_cint(1)));
    uex number;

    if (minpoly.data->length <= 2)
    {
        ex r = emake_rat();
        fmpz_neg(fmpq_numref(erat_data(r)), minpoly.data->coeffs + 0);
        fmpz_set(fmpq_denref(erat_data(r)), minpoly.data->coeffs + 1);
        return efix_rat(r);
    }

    number.setz(acb_to_ex(location.data));

    for (slong i = 0; i < minpoly.data->length; i++)
    {
        if (fmpz_is_zero(minpoly.data->coeffs + i))
            continue;
        ex t;
        if (i > 0)
        {
            t = emake_node(gs.sym_sPower.copy(), slot.copy(), emake_int_si(i));
            if (!fmpz_is_one(minpoly.data->coeffs + i))
            {
                t = emake_node(gs.sym_sTimes.copy(),
                        emake_int_copy(minpoly.data->coeffs + i), t);
            }
        }
        else
        {
            t = emake_int_copy(minpoly.data->coeffs + i);
        }
        p.push_back(wex(t));
    }

    return emake_node(gs.sym_sRoot.copy(), emake_node(gs.sym_sList.copy(),
                emake_node(gs.sym_sFunction.copy(),
                    emake_node(gs.sym_sPlus.copy(), p)
                ),
                number.release()
           ));
}


// for Exp[f*e]  where f is a number
bool qbarelem::set_exp_ex(er f, er e)
{
//std::cout << "set_exp_ex called e: " << ex_tostring_full(e) << std::endl;
//std::cout << "                  f: " << ex_tostring_full(f) << std::endl;

    xfmpq_t a;
    qbarelem t, s;

    if (ehas_head_sym_length(e, gs.sym_sTimes.get(), 2) && eis_finite_number(echild(e,1)))
    {
        wex p(num_Times(f, echild(e,1)));
        er e1 = p.get();
        er e2 = echild(e,2);

        if (eis_zero(ecmplx_real(e1)) && eget_fmpq(a.data, ecmplx_imag(e1)))
        {
            if (eis_sym(e2, gs.sym_sPi.get()))
            {
                return set_expipi(a.data);
            }
            else if (ehas_head_sym_length(e2, gs.sym_sArcCos.get(), 1))
            {
                return t.set_ex(echild(e2,1)) &&
                       t.circle_root(t, +1) &&
                       pow(t, a.data);
            }
            else if (ehas_head_sym_length(e2, gs.sym_sArcSin.get(), 1))
            {
                return t.set_ex(echild(e2,1)) &&
                       t.circle_root(t, -1) &&
                       t.mul_expipi(t, eget_crat_data(1, 2)) &&
                       pow(t, a.data);
            }
            else if (ehas_head_sym_length(e2, gs.sym_sArcTan.get(), 1))
            {
                fmpq_div_2exp(a.data, a.data, 1);
                return t.set_ex(echild(e2,1)) &&
                       t.mul_expipi(t, eget_crat_data(1, 2)) &&
                       s.apply_moebius(t, eget_cint_data(-1), eget_cint_data(1),
                                          eget_cint_data(0), eget_cint_data(1)) &&
                       t.apply_moebius(t, eget_cint_data(+1), eget_cint_data(1),
                                          eget_cint_data(0), eget_cint_data(1)) &&
                       t.pow(t, a.data) &&
                       (fmpq_neg(a.data, a.data), s.pow(s, a.data)) &&
                       mul(s, t);
            }
        }
    }
    

    return false;
}

bool qbarelem::set_ex(er e)
{
    if (eis_int(e))
    {
        set(eint_data(e));
        return true;
    }
    else if (eis_rat(e))
    {
        set(erat_data(e));
        return true;
    }

    qbarelem s, t;

    if (eis_cmplx(e))
    {
        if (!t.set_ex(ecmplx_imag(e)))
            return false;
        s.mul_expipi(t, eget_crat_data(1,2));
        if (!t.set_ex(ecmplx_real(e)))
            return false;
        add(s, t);
        return true;
    }

    if (!eis_node(e))
    {
        return false;
    }

    if (ehas_head_sym(e, gs.sym_sPlus.get()))
    {
        set(slong(0));
        for (ulong i = 0; i < elength(e); i++)
        {
            if (!t.set_ex(echild(e, i + 1)))
                return false;
            swap(s);
            add(s, t);
        }
        return true;
    }
    else if (ehas_head_sym(e, gs.sym_sTimes.get()))
    {
        set(slong(1));
        for (ulong i = 0; i < elength(e); i++)
        {
            if (!t.set_ex(echild(e, i + 1)))
                return false;
            swap(s);
            mul(s, t);
        }
        return true;
    }
    else if (ehas_head_sym_length(e, gs.sym_sPower.get(), 2))
    {
        if (eis_int(echild(e,2)))
        {
            return t.set_ex(echild(e,1)) && pow(t, eint_data(echild(e,2)));
        }
        else if (eis_rat(echild(e,2)))
        {
            return t.set_ex(echild(e,1)) && pow(t, erat_data(echild(e,2)));
        }
        else if (eis_sym(echild(e,1), gs.sym_sE.get()))
        {
            return set_exp_ex(eget_cint(1), echild(e,2));
        }
        else 
        {
            return false;
        }
    }
    else if (ehas_head_sym_length(e, gs.sym_sRe.get(), 1))
    {
        if (!set_ex(echild(e,1)))
            return false;
        return realpart(*this);
    }
    else if (ehas_head_sym_length(e, gs.sym_sExp.get(), 1))
    {
        return set_exp_ex(eget_cint(1), echild(e,1));
    }
    else if (ehas_head_sym_length(e, gs.sym_sCos.get(), 1))
    {
        if (!set_exp_ex(gs.const_i.get(), echild(e,1)))
            return false;
        fmpz num[3] = {1, 0, 1};
        fmpz den[2] = {0, 2};
        apply_rat_fxn(*this, num, 3, den, 2);
        return true;
    }
    else if (ehas_head_sym_length(e, gs.sym_sRoot.get(), 1) &&
             ehas_head_sym_length(echild(e,1), gs.sym_sList.get(), 2))
    {
        er body;
        er func = echild(e,1,1);
        wex var(emake_node(gs.sym_sSlot.copy(), emake_cint(1)));
        xfmpq_poly_t q;
        xacb_t x;

        if (elength(func) == 1)
        {
            body = echild(func, 1);
        }
        else if (elength(func) == 2)
        {
            var.reset(ecopychild(func, 1));
            body = echild(func, 2);            
        }
        else
        {
            return false;
        }

        if (!eget_fmpq_poly(q.data, body, var.get()))
            return false;

        if (!eget_acb(x.data, echild(e,1,2)))
            return false;

        return set(q.zpoly(), x.data);
    }

    return false;
}


ex rootreduce(er e)
{
    qbarelem r;

    if (r.set_ex(e))
    {
        newton_trim(r.minpoly.data, r.location.data, r.epsilon.data);
        return r.get_ex();
    }

    if (!eis_node(e))
        return ecopy(e);

    if (ehas_head_sym(e, gs.sym_sList.get()) || ehas_head_sym(e, gs.sym_sRule.get()))
    {
        uex a; a.init_push_backr(echild(e,0), elength(e));
        bool changed = false;
        for (ulong i = 0; i < elength(e); i++)
        {
            er ei = echild(e, i + 1);
            ex f = rootreduce(ei);
            changed = changed || (etor(f) != ei);
            a.push_back(f);
        }
        return changed ? a.release() : ecopy(e);
    }

    return ecopy(e);
}


ex dcode_sRootReduce(er e)
{
//std::cout << "dcode_sRootReduce: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sRootReduce.get()));

    if (elength(e) != 1)
        return _handle_message_argx1(e);

    return rootreduce(echild(e,1));
}


static ex _galois_group_order_irr3(xfmpq_poly_t &p)
{
    assert(fmpq_poly_degree(p.data) == 3);
    xfmpz_t disc;
    fmpz_poly_discriminant(disc.data, p.zpoly());
    return fmpz_is_square(disc.data) ? emake_cint(3) : emake_cint(6);
}

static ex _galois_group_order_irr4(xfmpq_poly_t &p)
{
    assert(fmpq_poly_degree(p.data) == 4);

    xfmpq_t A, B, C, b0, b1, b2, b3, b3b3, u, v, w;

    /* b3 = c3/(4c4), b2 = c2/(4c4), b1 = c1/(8c4), b0 = c0/(4c4)*/
    fmpq_set_fmpz_frac(b0.data, p.data->coeffs + 0, p.data->coeffs + 4);
    fmpq_set_fmpz_frac(b1.data, p.data->coeffs + 1, p.data->coeffs + 4);
    fmpq_set_fmpz_frac(b2.data, p.data->coeffs + 2, p.data->coeffs + 4);
    fmpq_set_fmpz_frac(b3.data, p.data->coeffs + 3, p.data->coeffs + 4);
    fmpq_div_2exp(b0.data, b0.data, 2);
    fmpq_div_2exp(b1.data, b1.data, 3);
    fmpq_div_2exp(b2.data, b2.data, 2);
    fmpq_div_2exp(b3.data, b3.data, 2);

    fmpq_pow_si(b3b3.data, b3.data, 2);

    /* A = b1 - b2 b3 + b3^3 */
    fmpq_sub(v.data, b3b3.data, b2.data);
    fmpq_mul(u.data, b3.data, v.data);
    fmpq_add(A.data, b1.data, u.data);

    /* B = -2 b2 + 3 b3^2 */
    fmpq_mul_fmpz(u.data, b3b3.data, eint_data(eget_cint(3)));
    fmpq_mul_2exp(v.data, b2.data, 1);
    fmpq_sub(B.data, u.data, v.data);

    /* C = -b0 + b2^2 + 2 b1 b3 - 4 b2 b3^2 + 3 b3^4 */
    fmpq_mul_fmpz(u.data, b3b3.data, eint_data(eget_cint(3)));
    fmpq_mul_2exp(v.data, b2.data, 2);
    fmpq_sub(w.data, u.data, v.data);
    fmpq_mul(u.data, b3.data, w.data);
    fmpq_mul_2exp(v.data, b1.data, 1);
    fmpq_add(w.data, v.data, u.data);
    fmpq_mul(u.data, b3.data, w.data);
    fmpq_pow_si(v.data, b2.data, 2);
    fmpq_sub(w.data, v.data, b0.data);
    fmpq_add(C.data, w.data, u.data);

std::cout << "C: " << std::endl;

    uint8_t Cprog[] = {3,4,0,0,0,0,0,0,0,0,3,5,255,255,255,255,255,255,255,255,1,5,5,0,2,4,4,5,3,5,2,0,0,0,0,0,0,0,1,5,5,1,1,5,5,3,2,4,4,5,3,5,1,0,0,0,0,0,0,0,1,5,5,2,1,5,5,2,2,4,4,5,3,5,252,255,255,255,255,255,255,255,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,3,0,0,0,0,0,0,0,1,5,5,3,1,5,5,3,1,5,5,3,1,5,5,3,2,4,4,5};
    std::vector<xfmpq_t> stack;
    stack.resize(10);
    fmpq_set(stack[0].data, b0.data);
    fmpq_set(stack[1].data, b1.data);
    fmpq_set(stack[2].data, b2.data);
    fmpq_set(stack[3].data, b3.data);
    eval_poly_fmpq(stack, Cprog, sizeof(Cprog));
std::cout << "C: " << stack[4] << std::endl;


    if (fmpq_is_zero(A.data))
    {
        fmpq_pow_si(u.data, B.data, 2);
        fmpq_mul_2exp(v.data, C.data, 2);
        fmpq_sub(w.data, u.data, v.data);
        if (w.is_square())
        {
            return emake_cint(4); // C2 x C2
        }
        else
        {
            fmpq_mul(u.data, w.data, C.data);
            return u.is_square() ? emake_cint(4) : emake_cint(8);
        }
    }
    else
    {
        fmpq_div(u.data, C.data, A.data);
        xfmpq_poly_t q;
        fmpq_poly_set_coeff_fmpz(q.data, 3, eint_data(eget_cint(1)));
        fmpq_poly_set_coeff_fmpq(q.data, 2, u.data);
        fmpq_poly_set_coeff_fmpq(q.data, 1, B.data);
        fmpq_poly_set_coeff_fmpq(q.data, 0, A.data);

        fmpz_poly_factor_t fac;
        fmpz_poly_factor_init(fac);
        fmpz_poly_factor(fac, q.zpoly());
        std::cout << "qfac: " << std::endl;
        fmpz_poly_factor_print(fac);

        if (fac->num == 3)
        {
            fmpz_poly_factor_clear(fac);
            return emake_cint(4); // C2 x C2
        }
        else if (fac->num == 2)
        {
            assert(fmpz_poly_degree(fac->p + 0) == 1);
            assert(fmpz_poly_degree(fac->p + 1) == 2);
            fmpz * fcoeffs = (fac->p + 1)->coeffs;
            xfmpz_t c0c2, c1c1, t1, t2;
            fmpz_mul(c0c2.data, fcoeffs + 0, fcoeffs + 2);
            fmpz_mul(c1c1.data, fcoeffs + 1, fcoeffs + 1);
            fmpz_mul_2exp(t1.data, c0c2.data, 2);
            fmpz_sub(t2.data, c1c1.data, t1.data);
            fmpz_mul(t1.data, c0c2.data, t2.data);
            fmpz_poly_factor_clear(fac);
            return fmpz_is_square(t1.data) ? emake_cint(4) : emake_cint(8);            
        }
        else
        {
            assert(fac->num == 1);
            xfmpz_t disc;
            fmpz_poly_discriminant(disc.data, q.zpoly());
            fmpz_poly_factor_clear(fac);
            return fmpz_is_square(disc.data) ? emake_cint(12) : emake_cint(24);
        }
    }
}

static ex _galois_group_order_irr5(xfmpq_poly_t &p)
{
    assert(fmpq_poly_degree(p.data) == 5);
    fmpq_poly_make_monic(p.data, p.data);

    xfmpq_t disc;
    _fmpz_poly_discriminant(fmpq_numref(disc.data), p.data->coeffs, 6);
    fmpz_pow_ui(fmpq_denref(disc.data), p.data->den, 8);
    fmpq_canonicalise(disc.data);

    bool disc_is_square = disc.is_square();
//std::cout << "disc: " << disc.tostring() << std::endl;

    xfmpq_poly_t q, r;
    xfmpq_t t;
    fmpq_poly_get_coeff_fmpq(t.data, p.data, 4);
    fmpq_div_fmpz(t.data, t.data, eint_data(eget_cint(-5)));
    fmpq_poly_set_coeff_fmpz(r.data, 1, eint_data(eget_cint(1)));
    fmpq_poly_set_coeff_fmpq(r.data, 0, t.data);
    fmpq_poly_compose(q.data, p.data, r.data);

//std::cout << "q: " << q.tostring() << std::endl;

    std::vector<xfmpq_t> stack(20);
    fmpq_poly_get_coeff_fmpq(stack[0].data, q.data, 3); // a2
    fmpq_poly_get_coeff_fmpq(stack[1].data, q.data, 2); // a3
    fmpq_poly_get_coeff_fmpq(stack[2].data, q.data, 1); // a4
    fmpq_poly_get_coeff_fmpq(stack[3].data, q.data, 0); // a5

    /* compute 6th degree resolvent with coefficient rules

        {{6} -> 1, {5} -> 8*a4, {4} -> 2*a2*a3^2 - 6*a2^2*a4 + 
           40*a4^2 - 50*a3*a5, {3} -> -2*a3^4 + 21*a2*a3^2*a4 - 
           40*a2^2*a4^2 + 160*a4^3 - 15*a2^2*a3*a5 - 400*a3*a4*a5 + 
           125*a2*a5^2, {2} -> a2^2*a3^4 - 6*a2^3*a3^2*a4 - 
           8*a3^4*a4 + 9*a2^4*a4^2 + 76*a2*a3^2*a4^2 - 
           136*a2^2*a4^3 + 400*a4^4 - 50*a2*a3^3*a5 + 
           90*a2^2*a3*a4*a5 - 1400*a3*a4^2*a5 + 625*a3^2*a5^2 + 
           500*a2*a4*a5^2, {1} -> -2*a2*a3^6 + 19*a2^2*a3^4*a4 - 
           51*a2^3*a3^2*a4^2 + 3*a3^4*a4^2 + 32*a2^4*a4^3 + 
           76*a2*a3^2*a4^3 - 256*a2^2*a4^4 + 512*a4^5 - 
           31*a2^3*a3^3*a5 - 58*a3^5*a5 + 117*a2^4*a3*a4*a5 + 
           105*a2*a3^3*a4*a5 + 260*a2^2*a3*a4^2*a5 - 
           2400*a3*a4^3*a5 - 108*a2^5*a5^2 - 325*a2^2*a3^2*a5^2 + 
           525*a2^3*a4*a5^2 + 2750*a3^2*a4*a5^2 - 500*a2*a4^2*a5^2 + 
           625*a2*a3*a5^3 - 3125*a5^4, 
         {0} -> a3^8 - 13*a2*a3^6*a4 + a2^5*a3^2*a4^2 + 
           65*a2^2*a3^4*a4^2 - 4*a2^6*a4^3 - 128*a2^3*a3^2*a4^3 + 
           17*a3^4*a4^3 + 48*a2^4*a4^4 - 16*a2*a3^2*a4^4 - 
           192*a2^2*a4^5 + 256*a4^6 - 4*a2^5*a3^3*a5 - 
           12*a2^2*a3^5*a5 + 18*a2^6*a3*a4*a5 + 12*a2^3*a3^3*a4*a5 - 
           124*a3^5*a4*a5 + 196*a2^4*a3*a4^2*a5 + 
           590*a2*a3^3*a4^2*a5 - 160*a2^2*a3*a4^3*a5 - 
           1600*a3*a4^4*a5 - 27*a2^7*a5^2 - 150*a2^4*a3^2*a5^2 - 
           125*a2*a3^4*a5^2 - 99*a2^5*a4*a5^2 - 
           725*a2^2*a3^2*a4*a5^2 + 1200*a2^3*a4^2*a5^2 + 
           3250*a3^2*a4^2*a5^2 - 2000*a2*a4^3*a5^2 - 
           1250*a2*a3*a4*a5^3 + 3125*a2^2*a5^4 - 9375*a4*a5^4}
    */
    fmpq_poly_zero(r.data);
    fmpq_poly_set_coeff_fmpz(r.data, 6, eint_data(eget_cint(1)));
    {
    uint8_t prog5[] = {3,4,0,0,0,0,0,0,0,0,3,5,8,0,0,0,0,0,0,0,1,5,5,2,2,4,4,5};
    uint8_t prog4[] = {3,4,0,0,0,0,0,0,0,0,3,5,250,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,2,2,4,4,5,3,5,2,0,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,1,2,4,4,5,3,5,206,255,255,255,255,255,255,255,1,5,5,1,1,5,5,3,2,4,4,5,3,5,40,0,0,0,0,0,0,0,1,5,5,2,1,5,5,2,2,4,4,5};
    uint8_t prog3[] = {3,4,0,0,0,0,0,0,0,0,3,5,241,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,3,2,4,4,5,3,5,216,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,2,4,4,5,3,5,21,0,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,2,4,4,5,3,5,125,0,0,0,0,0,0,0,1,5,5,0,1,5,5,3,1,5,5,3,2,4,4,5,3,5,254,255,255,255,255,255,255,255,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,2,4,4,5,3,5,112,254,255,255,255,255,255,255,1,5,5,1,1,5,5,2,1,5,5,3,2,4,4,5,3,5,160,0,0,0,0,0,0,0,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5};
    uint8_t prog2[] = {3,4,0,0,0,0,0,0,0,0,3,5,9,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,2,4,4,5,3,5,250,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,2,4,4,5,3,5,1,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,2,4,4,5,3,5,90,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,2,1,5,5,3,2,4,4,5,3,5,120,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,206,255,255,255,255,255,255,255,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,3,2,4,4,5,3,5,76,0,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,2,4,4,5,3,5,244,1,0,0,0,0,0,0,1,5,5,0,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,248,255,255,255,255,255,255,255,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,2,4,4,5,3,5,113,2,0,0,0,0,0,0,1,5,5,1,1,5,5,1,1,5,5,3,1,5,5,3,2,4,4,5,3,5,136,250,255,255,255,255,255,255,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,3,2,4,4,5,3,5,144,1,0,0,0,0,0,0,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5};
    uint8_t prog1[] = {3,4,0,0,0,0,0,0,0,0,3,5,148,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,3,1,5,5,3,2,4,4,5,3,5,117,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,2,1,5,5,3,2,4,4,5,3,5,32,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,225,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,3,2,4,4,5,3,5,205,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,2,4,4,5,3,5,13,2,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,19,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,2,4,4,5,3,5,187,254,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,3,1,5,5,3,2,4,4,5,3,5,4,1,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,3,2,4,4,5,3,5,0,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,254,255,255,255,255,255,255,255,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,2,4,4,5,3,5,105,0,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,3,2,4,4,5,3,5,76,0,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,113,2,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,3,1,5,5,3,1,5,5,3,2,4,4,5,3,5,12,254,255,255,255,255,255,255,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,198,255,255,255,255,255,255,255,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,3,2,4,4,5,3,5,3,0,0,0,0,0,0,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,2,4,4,5,3,5,190,10,0,0,0,0,0,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,160,246,255,255,255,255,255,255,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,3,2,4,4,5,3,5,0,2,0,0,0,0,0,0,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,203,243,255,255,255,255,255,255,1,5,5,3,1,5,5,3,1,5,5,3,1,5,5,3,2,4,4,5};
    uint8_t prog0[] = {3,4,0,0,0,0,0,0,0,0,3,5,229,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,3,1,5,5,3,2,4,4,5,3,5,18,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,2,1,5,5,3,2,4,4,5,3,5,252,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,252,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,3,2,4,4,5,3,5,1,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,2,4,4,5,3,5,157,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,106,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,3,1,5,5,3,2,4,4,5,3,5,196,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,3,2,4,4,5,3,5,48,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,12,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,3,2,4,4,5,3,5,128,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,176,4,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,244,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,3,2,4,4,5,3,5,65,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,2,4,4,5,3,5,43,253,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,96,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,3,2,4,4,5,3,5,64,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,53,12,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,3,1,5,5,3,1,5,5,3,1,5,5,3,2,4,4,5,3,5,243,255,255,255,255,255,255,255,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,2,4,4,5,3,5,131,255,255,255,255,255,255,255,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,3,1,5,5,3,2,4,4,5,3,5,78,2,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,3,2,4,4,5,3,5,240,255,255,255,255,255,255,255,1,5,5,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,30,251,255,255,255,255,255,255,1,5,5,0,1,5,5,1,1,5,5,2,1,5,5,3,1,5,5,3,1,5,5,3,2,4,4,5,3,5,48,248,255,255,255,255,255,255,1,5,5,0,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,1,0,0,0,0,0,0,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,2,4,4,5,3,5,132,255,255,255,255,255,255,255,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,3,2,4,4,5,3,5,17,0,0,0,0,0,0,0,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,178,12,0,0,0,0,0,0,1,5,5,1,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,3,1,5,5,3,2,4,4,5,3,5,192,249,255,255,255,255,255,255,1,5,5,1,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,3,2,4,4,5,3,5,0,1,0,0,0,0,0,0,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,1,5,5,2,2,4,4,5,3,5,97,219,255,255,255,255,255,255,1,5,5,2,1,5,5,3,1,5,5,3,1,5,5,3,1,5,5,3,2,4,4,5};
    eval_poly_fmpq(stack, prog5, sizeof(prog5));
    fmpq_poly_set_coeff_fmpq(r.data, 5, stack[4].data);
    eval_poly_fmpq(stack, prog4, sizeof(prog4));
    fmpq_poly_set_coeff_fmpq(r.data, 4, stack[4].data);
    eval_poly_fmpq(stack, prog3, sizeof(prog3));
    fmpq_poly_set_coeff_fmpq(r.data, 3, stack[4].data);
    eval_poly_fmpq(stack, prog2, sizeof(prog2));
    fmpq_poly_set_coeff_fmpq(r.data, 2, stack[4].data);
    eval_poly_fmpq(stack, prog1, sizeof(prog1));
    fmpq_poly_set_coeff_fmpq(r.data, 1, stack[4].data);
    eval_poly_fmpq(stack, prog0, sizeof(prog0));
    fmpq_poly_set_coeff_fmpq(r.data, 0, stack[4].data);
    }

//std::cout << "r: " << r.tostring() << std::endl;

    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, r.zpoly());
//std::cout << "r fac: " << std::endl;
//fmpz_poly_factor_print(fac);

    xfmpq_t theta;
    bool solvable = false;
    for (slong i = 0; i < fac->num; i++)
    {
        /* look for a linear^1 or linear^6 */
        if ((fac->p + i)->length == 2 && (fac->exp[i] == 1 || fac->exp[i] == 6))
        {
            solvable = true;
            fmpq_set_fmpz_frac(theta.data, (fac->p + i)->coeffs + 0, (fac->p + i)->coeffs + 1);
            fmpq_neg(theta.data, theta.data);
            break;
        }
    }
    fmpz_poly_factor_clear(fac);

    if (!solvable)
    {
        return disc_is_square ? emake_cint(60) : emake_cint(120);
    }

    if (!disc_is_square)
    {
        return emake_cint(20);
    }

    xfmpq_t disc_sqrt;
    fmpz_sqrt(fmpq_numref(disc_sqrt.data), fmpq_numref(disc.data));
    fmpz_sqrt(fmpq_denref(disc_sqrt.data), fmpq_denref(disc.data));

    /* compute 10th degree resolvent with coefficient rules

        {{10} -> 1, {8} -> 3*a2, {7} -> a3, {6} -> 3*a2^2 - 3*a4, 
         {5} -> 2*a2*a3 - 11*a5, {4} -> a2^3 - a3^2 - 2*a2*a4, 
         {3} -> a2^2*a3 - 4*a3*a4 - 4*a2*a5, 
         {2} -> -(a2*a3^2) + a2^2*a4 - 4*a4^2 + 7*a3*a5, 
         {1} -> -a3^3 - a2^2*a5 + 4*a4*a5, 
         {0} -> -(a3^2*a4) + a2*a3*a5 - a5^2}
    */
    fmpq_poly_zero(r.data);
    fmpq_poly_set_coeff_fmpz(r.data, 10, eint_data(eget_cint(1)));
    {
    uint8_t prog8[] = {3,4,0,0,0,0,0,0,0,0,3,5,3,0,0,0,0,0,0,0,1,5,5,0,2,4,4,5};
    uint8_t prog7[] = {3,4,0,0,0,0,0,0,0,0,3,5,1,0,0,0,0,0,0,0,1,5,5,1,2,4,4,5};
    uint8_t prog6[] = {3,4,0,0,0,0,0,0,0,0,3,5,3,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,2,4,4,5,3,5,253,255,255,255,255,255,255,255,1,5,5,2,2,4,4,5};
    uint8_t prog5[] = {3,4,0,0,0,0,0,0,0,0,3,5,2,0,0,0,0,0,0,0,1,5,5,0,1,5,5,1,2,4,4,5,3,5,245,255,255,255,255,255,255,255,1,5,5,3,2,4,4,5};
    uint8_t prog4[] = {3,4,0,0,0,0,0,0,0,0,3,5,1,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,0,2,4,4,5,3,5,254,255,255,255,255,255,255,255,1,5,5,0,1,5,5,2,2,4,4,5,3,5,255,255,255,255,255,255,255,255,1,5,5,1,1,5,5,1,2,4,4,5};
    uint8_t prog3[] = {3,4,0,0,0,0,0,0,0,0,3,5,1,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,1,2,4,4,5,3,5,252,255,255,255,255,255,255,255,1,5,5,0,1,5,5,3,2,4,4,5,3,5,252,255,255,255,255,255,255,255,1,5,5,1,1,5,5,2,2,4,4,5};
    uint8_t prog2[] = {3,4,0,0,0,0,0,0,0,0,3,5,1,0,0,0,0,0,0,0,1,5,5,0,1,5,5,0,1,5,5,2,2,4,4,5,3,5,255,255,255,255,255,255,255,255,1,5,5,0,1,5,5,1,1,5,5,1,2,4,4,5,3,5,7,0,0,0,0,0,0,0,1,5,5,1,1,5,5,3,2,4,4,5,3,5,252,255,255,255,255,255,255,255,1,5,5,2,1,5,5,2,2,4,4,5};
    uint8_t prog1[] = {3,4,0,0,0,0,0,0,0,0,3,5,255,255,255,255,255,255,255,255,1,5,5,0,1,5,5,0,1,5,5,3,2,4,4,5,3,5,255,255,255,255,255,255,255,255,1,5,5,1,1,5,5,1,1,5,5,1,2,4,4,5,3,5,4,0,0,0,0,0,0,0,1,5,5,2,1,5,5,3,2,4,4,5};
    uint8_t prog0[] = {3,4,0,0,0,0,0,0,0,0,3,5,1,0,0,0,0,0,0,0,1,5,5,0,1,5,5,1,1,5,5,3,2,4,4,5,3,5,255,255,255,255,255,255,255,255,1,5,5,1,1,5,5,1,1,5,5,2,2,4,4,5,3,5,255,255,255,255,255,255,255,255,1,5,5,3,1,5,5,3,2,4,4,5};
    eval_poly_fmpq(stack, prog8, sizeof(prog8));
    fmpq_poly_set_coeff_fmpq(r.data, 8, stack[4].data);
    eval_poly_fmpq(stack, prog7, sizeof(prog7));
    fmpq_poly_set_coeff_fmpq(r.data, 7, stack[4].data);
    eval_poly_fmpq(stack, prog6, sizeof(prog6));
    fmpq_poly_set_coeff_fmpq(r.data, 6, stack[4].data);
    eval_poly_fmpq(stack, prog5, sizeof(prog5));
    fmpq_poly_set_coeff_fmpq(r.data, 5, stack[4].data);
    eval_poly_fmpq(stack, prog4, sizeof(prog4));
    fmpq_poly_set_coeff_fmpq(r.data, 4, stack[4].data);
    eval_poly_fmpq(stack, prog3, sizeof(prog3));
    fmpq_poly_set_coeff_fmpq(r.data, 3, stack[4].data);
    eval_poly_fmpq(stack, prog2, sizeof(prog2));
    fmpq_poly_set_coeff_fmpq(r.data, 2, stack[4].data);
    eval_poly_fmpq(stack, prog1, sizeof(prog1));
    fmpq_poly_set_coeff_fmpq(r.data, 1, stack[4].data);
    eval_poly_fmpq(stack, prog0, sizeof(prog0));
    fmpq_poly_set_coeff_fmpq(r.data, 0, stack[4].data);
    }

    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, r.zpoly());
//std::cout << "r fac: " << std::endl;
//fmpz_poly_factor_print(fac);
    assert(fac->num == 2);
    assert((fac->p + 0)->length == 6);
    assert((fac->p + 1)->length == 6);
    assert(fmpz_is_zero((fac->p + 0)->coeffs + 4));
    assert(fmpz_is_zero((fac->p + 1)->coeffs + 4));

    xfmpq_t u, v, w, b[6], c[6];
    for (slong i = 2; i <= 5; i++)
    {
        fmpq_set_fmpz_frac(u.data, (fac->p + 0)->coeffs + 5 - i, (fac->p + 0)->coeffs + 5);
        fmpq_set_fmpz_frac(v.data, (fac->p + 1)->coeffs + 5 - i, (fac->p + 1)->coeffs + 5);
        fmpq_add(b[i].data, u.data, v.data);
        fmpq_sub(c[i].data, u.data, v.data);
        fmpq_div_2exp(b[i].data, b[i].data, 1);
        fmpq_div_2exp(c[i].data, c[i].data, 1);
        fmpq_div(c[i].data, c[i].data, disc_sqrt.data);
    }

    fmpq_set(stack[0].data, b[2].data);
    fmpq_set(stack[1].data, b[3].data);
    fmpq_set(stack[2].data, b[4].data);
    fmpq_set(stack[3].data, b[5].data);
    fmpq_set(stack[4].data, c[2].data);
    fmpq_set(stack[5].data, c[3].data);
    fmpq_set(stack[6].data, c[4].data);
    fmpq_set(stack[7].data, c[5].data);
    fmpq_set(stack[8].data, disc.data);

    /* evaluate A = 
        -81*c2^2(4*b2*c2^2 + 15*c3^2)*disc^2
         - 25*b2^2*(4*b2^3 + 135*b3^2)
         + 45(8*b2^3*c2^2 - 15*b2^2*c3^2 + 180*b3*b2*c2*c3 - 135*b3^2*c2^2)*disc
    */
    xfmpq_t A;
    {
    uint8_t prog[] = {3,9,0,0,0,0,0,0,0,0,3,10,156,255,255,255,255,255,255,255,1,10,10,0,1,10,10,0,1,10,10,0,1,10,10,0,1,10,10,0,2,9,9,10,3,10,104,1,0,0,0,0,0,0,1,10,10,0,1,10,10,0,1,10,10,0,1,10,10,4,1,10,10,4,1,10,10,8,2,9,9,10,3,10,209,242,255,255,255,255,255,255,1,10,10,0,1,10,10,0,1,10,10,1,1,10,10,1,2,9,9,10,3,10,93,253,255,255,255,255,255,255,1,10,10,0,1,10,10,0,1,10,10,5,1,10,10,5,1,10,10,8,2,9,9,10,3,10,164,31,0,0,0,0,0,0,1,10,10,0,1,10,10,1,1,10,10,4,1,10,10,5,1,10,10,8,2,9,9,10,3,10,188,254,255,255,255,255,255,255,1,10,10,0,1,10,10,4,1,10,10,4,1,10,10,4,1,10,10,4,1,10,10,8,1,10,10,8,2,9,9,10,3,10,69,232,255,255,255,255,255,255,1,10,10,1,1,10,10,1,1,10,10,4,1,10,10,4,1,10,10,8,2,9,9,10,3,10,65,251,255,255,255,255,255,255,1,10,10,4,1,10,10,4,1,10,10,5,1,10,10,5,1,10,10,8,1,10,10,8,2,9,9,10};
    eval_poly_fmpq(stack, prog, sizeof(prog));
    fmpq_set(A.data, stack[9].data);
    }
    
    /* evaluate B = 
        c2*(28*b2^2*c2^2 + 90*b4*c2^2 - 135*b3*c3*c2 - 90*b2*c4*c2 + 45*b2*c3^2)*disc
        - 5*b2*(2*b2^3*c2 - 10*b2^2*c4 + 10*b4*b2*c2 + 15*b3*b2*c3 - 45*b3^2*c2)
        - 18*c2^5*disc^2
    */
    xfmpq_t B;
    {
    uint8_t prog[] = {3,9,0,0,0,0,0,0,0,0,3,10,246,255,255,255,255,255,255,255,1,10,10,0,1,10,10,0,1,10,10,0,1,10,10,0,1,10,10,4,2,9,9,10,3,10,50,0,0,0,0,0,0,0,1,10,10,0,1,10,10,0,1,10,10,0,1,10,10,6,2,9,9,10,3,10,181,255,255,255,255,255,255,255,1,10,10,0,1,10,10,0,1,10,10,1,1,10,10,5,2,9,9,10,3,10,206,255,255,255,255,255,255,255,1,10,10,0,1,10,10,0,1,10,10,2,1,10,10,4,2,9,9,10,3,10,28,0,0,0,0,0,0,0,1,10,10,0,1,10,10,0,1,10,10,4,1,10,10,4,1,10,10,4,1,10,10,8,2,9,9,10,3,10,225,0,0,0,0,0,0,0,1,10,10,0,1,10,10,1,1,10,10,1,1,10,10,4,2,9,9,10,3,10,166,255,255,255,255,255,255,255,1,10,10,0,1,10,10,4,1,10,10,4,1,10,10,6,1,10,10,8,2,9,9,10,3,10,45,0,0,0,0,0,0,0,1,10,10,0,1,10,10,4,1,10,10,5,1,10,10,5,1,10,10,8,2,9,9,10,3,10,121,255,255,255,255,255,255,255,1,10,10,1,1,10,10,4,1,10,10,4,1,10,10,5,1,10,10,8,2,9,9,10,3,10,90,0,0,0,0,0,0,0,1,10,10,2,1,10,10,4,1,10,10,4,1,10,10,4,1,10,10,8,2,9,9,10,3,10,238,255,255,255,255,255,255,255,1,10,10,4,1,10,10,4,1,10,10,4,1,10,10,4,1,10,10,4,1,10,10,8,1,10,10,8,2,9,9,10};
    eval_poly_fmpq(stack, prog, sizeof(prog));
    fmpq_set(B.data, stack[9].data);
    }

    /* check if 3*A +- 54*B*disc_sqrt are both squares */

    fmpq_mul_fmpz(u.data, A.data, eint_data(eget_cint(3)));
    fmpq_mul_fmpz(w.data, B.data, eint_data(eget_cint(54)));
    fmpq_mul(v.data, w.data, disc_sqrt.data);

    fmpq_add(w.data, u.data, v.data);
    if (!w.is_square())
        return emake_cint(10);

    fmpq_sub(w.data, u.data, v.data);
    if (!w.is_square())
        return emake_cint(10);

    return emake_cint(5);
}


static ex _galois_group_order_irr(xfmpq_poly_t &p)
{
    slong d = fmpq_poly_degree(p.data);
    switch (d)
    {
        case 5:
            return _galois_group_order_irr5(p);
        case 4:
            return _galois_group_order_irr4(p);
        case 3:
            return _galois_group_order_irr3(p);
        case 2:
            return emake_cint(2);
        case 1:
            return emake_cint(1);
        default:
            return d < 1 ? emake_cint(1) : gs.sym_s$Failed.copy();
    }
}

static ex _galois_group_order(xfmpq_poly_t &p)
{
//std::cout << "_galois_group_order: " << p.tostring() << std::endl;
//std::cout << "_galois_group_order: " << std::endl;
//fmpz_poly_print_pretty(p.zpoly(), "X");
//std::cout << std::endl;

    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, p.zpoly());
//std::cout << "fac: " << std::endl;
//fmpz_poly_factor_print(fac);
    bool irr = fac->num == 1;
    fmpz_poly_factor_clear(fac);
    return (!irr) ? emake_int_si(0) : _galois_group_order_irr(p);
}

#if 0
ex dcode_sGaloisGroup(er e)
{
//std::cout << "dcode_sGaloisGroup: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sGaloisGroup.get()));

    if (elength(e) != 1)
    {
        return ecopy(e);
    }

    er X = echild(e,1);
    if (ehas_head_sym(X, gs.sym_sList.get()))
    {
        size_t n = elength(X);
        xfmpq_poly_t p;
        for (size_t i = 1; i <= n; i++)
        {
            er Xi = echild(X,i);
            if (eis_int(Xi))
            {
                fmpq_poly_set_coeff_fmpz(p.data, i-1, eint_data(Xi));
            }
            else if (eis_rat(Xi))
            {
                fmpq_poly_set_coeff_fmpq(p.data, i-1, erat_data(Xi));
            }
            else
            {
                return ecopy(e);
            }
        }
        return _galois_group_order(p);
    }
    else
    {
        return ecopy(e);
    }
}
#endif



static ex _roots(xfmpq_poly_t &p);

static ex _roots_irr1(const fmpz_poly_t p)
{
    assert(fmpz_poly_degree(p) == 1);
    uex r1(emake_rat());
    fmpq_set_fmpz_frac(erat_data(r1.get()), p->coeffs + 0, p->coeffs + 1);
    fmpq_neg(erat_data(r1.get()), erat_data(r1.get()));
    r1.setz(efix_rat(r1.release()));
    return emake_node(gs.sym_sList.copy(), r1.release());
}

static ex _roots_irr2(const fmpz_poly_t p)
{
    assert(fmpz_poly_degree(p) == 2);
printf("not implemented\n");
abort();
return nullptr;
}



static ex _roots_irr3(const fmpz_poly_t p)
{
    assert(fmpz_poly_degree(p) == 3);

    /* the u1 u2 are roots of
        -(-p2^2 + 3*p1*p3)/(9*p3^2)
        + (-2*p2^3 + 9*p1*p2*p3 - 27*p0*p3^2)/(3*p3*(-p2^2 + 3*p1*p3))*U
        + U^2
    */
printf("not implemented\n");
abort();
return nullptr;
}



static ex _roots_irr4(const fmpz_poly_t p)
{
    assert(fmpz_poly_degree(p) == 4);

printf("not implemented\n");
abort();
return nullptr;
}

static ex _roots_irr5(const fmpz_poly_t p)
{
    assert(fmpz_poly_degree(p) == 5);
    uex r1(emake_cint(0));
    uex r2(emake_cint(0));
    uex r3(emake_cint(0));
    uex r4(emake_cint(0));
    uex r5(emake_cint(0));
    return emake_node(gs.sym_sList.copy(), r1.release(), r2.release(), r3.release(), r4.release(), r5.release());
}


static ex _roots_irr(const fmpz_poly_t p)
{
    slong d = fmpz_poly_degree(p);
    switch (d)
    {
        case 5:
            return _roots_irr5(p);
        case 4:
            return _roots_irr4(p);
        case 3:
            return _roots_irr3(p);
        case 2:
            return _roots_irr2(p);
        case 1:
            return _roots_irr1(p);
        default:
            return d < 1 ? emake_node(gs.sym_sList.copy()) : gs.sym_s$Failed.copy();
    }
}

static ex _roots(xfmpq_poly_t &p)
{
    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, p.zpoly());
    uex res; res.init_push_backr(gs.sym_sList.get(), fmpq_poly_degree(p.data));
    for (slong i = 0; i < fac->num; i++)
    {
        wex irroots(_roots_irr(fac->p + i));
        for (size_t j = 1; j <= elength(irroots.get()); j++)
        {
            for (slong k = 0; k < fac->exp[i]; k++)
            {
                res.push_back(irroots.copychild(j));
            }
        }
    }
    return res.release();
}


#define FMPZ_PTR_SWAP(a, b) \
    do {                    \
        fmpz * __tt__ = b;  \
        b = a;              \
        a = __tt__;         \
    } while (0)


/* f(x) = g(h(x)) */
int _fmpz_poly_decompose(
    fmpz * g,           /* length r + 1 */
    fmpz * h,           /* length s + 1 */
    const fmpz * ff,    /* length r*s + 1 */
    slong r,
    slong s)
{
    int success = 1;
    slong hshift, k, j, i;
    fmpz_t cnum, cden, gg, hden, tden, hdenbar;
    fmpz * t, * tq, * tr;
/*
printf("decompose ff: "); _fmpz_poly_fprint_pretty(stdout, ff, 1 + r*s, "x"); printf("\n");
*/
    assert(r > 1);
    assert(s > 1);

    t = _fmpz_vec_init(r*s + 1);
    tq = _fmpz_vec_init(r*s + 1);
    tr = _fmpz_vec_init(r*s + 1);

    fmpz_init(cnum);
    fmpz_init(cden);
    fmpz_init(gg);
    fmpz_init(hden);
    fmpz_init(tden);
    fmpz_init(hdenbar);

    /* h = x^s */
    k = s; 
    hshift = s;
    fmpz_one(h + hshift);
    fmpz_one(hden);
    for (i = 0; i < hshift; i++)
        fmpz_zero(h + i);
/*
printf("h: ("); _fmpz_poly_fprint_pretty(stdout, h + hshift, 1 + s - hshift, "x"); flint_printf(")*x^%wd", hshift);
printf("/"); fmpz_print(hden); printf("\n");
*/
    for (k = 1; k < s; k++)
    {
        _fmpz_poly_pow(t, h + hshift, 1 + s - hshift, r);
        fmpz_pow_ui(tden, hden, r);

/*
printf("h^r: ("); _fmpz_poly_fprint_pretty(stdout, t, 1 + r*(s - hshift), "x"); flint_printf(")*x^%wd", r*hshift);
printf("/"); fmpz_print(tden); printf("\n");
*/

        if (k <= r*(s - hshift))
        {
            _fmpq_sub(cnum, cden, ff + r*s - k, ff + r*s,
                                  t + r*(s - hshift) - k, tden);
            fmpz_mul_ui(cden, cden, r);
        }
        else
        {
            fmpz_set(cnum, ff + r*s - k);
            fmpz_mul_ui(cden, ff + r*s, r);
        }
/*
printf("c: "); fmpz_print(cnum); printf("/"); fmpz_print(cden); printf("\n");
*/
        if (fmpz_is_zero(cnum))
            continue;

        fmpz_gcd(gg, cnum, cden);
        fmpz_divexact(cnum, cnum, gg);
        fmpz_divexact(cden, cden, gg);
/*
printf("c: "); fmpz_print(cnum); printf("/"); fmpz_print(cden); printf("\n");
*/
        fmpz_gcd(gg, hden, cden);
        fmpz_divexact(hdenbar, hden, gg);
        fmpz_divexact(cden, cden, gg);

        for (i = hshift; i <= s; i++)
            fmpz_mul(h + i, h + i, cden);

        hshift = FLINT_MIN(hshift, s - k);
        fmpz_mul(h + s - k, cnum, hdenbar);
        fmpz_mul(hden, hden, cden);
/*
printf("h: ("); _fmpz_poly_fprint_pretty(stdout, h + hshift, 1 + s - hshift, "x"); flint_printf(")*x^%wd", hshift);
printf("/"); fmpz_print(hden); printf("\n");
*/
    }

/*
printf("f: "); _fmpz_poly_fprint_pretty(stdout, ff, r*s + 1, "x"); printf("\n");
*/
    j = 0;

    fmpz_set(g + 0, ff + 0);

    for (k = 1; k < hshift; k++)
        if (!fmpz_is_zero(ff + k))
            goto fail;

    if (!_fmpz_poly_divrem(t, tr, ff + hshift, r*s + 1 - s*j - hshift,
                                             h + hshift, s + 1 - hshift, 1))
        goto fail;

    for (k = 0; k < r*s + 1 - s*j - hshift; k++)
        if (!fmpz_is_zero(tr + k))
            goto fail;

    for (j = 1; j < r; j++)
    {
        /* t has degree r*s - s*j */
/*
printf("t: "); _fmpz_poly_fprint_pretty(stdout, t, r*s + 1 - s*j, "x"); printf("\n");
*/

        fmpz_swap(g + j, t + 0);

        for (k = 1; k < hshift; k++)
            if (!fmpz_is_zero(t + k))
                goto fail;

/*
printf("num: "); _fmpz_poly_fprint_pretty(stdout, t + hshift, r*s + 1 - s*j - hshift, "x"); printf("\n");
printf("den: "); _fmpz_poly_fprint_pretty(stdout, h + hshift, s + 1 - hshift, "x"); printf("\n");
*/

        if (!_fmpz_poly_divrem(tq, tr, t + hshift, r*s + 1 - hshift - s*j,
                                       h + hshift, s + 1 - hshift, 1))
            goto fail;
/*
printf("quo: "); _fmpz_poly_fprint_pretty(stdout, t2, r*s + 1 - s*(j+1), "x"); printf("\n");
printf("rem: "); _fmpz_poly_fprint_pretty(stdout, tr, r*s + 1 - s*j - hshift, "x"); printf("\n");
*/

        for (k = 0; k < r*s + 1 - hshift - s*j; k++)
            if (!fmpz_is_zero(tr + k))
                goto fail;

        FMPZ_PTR_SWAP(t, tq);
    }

    /* t has degree 0 */
/*
printf("t: "); _fmpz_poly_fprint_pretty(stdout, t, 1, "x"); printf("\n");
*/

    fmpz_swap(g + r, t + 0);

cleanup:

    _fmpz_vec_clear(t, r*s + 1);
    _fmpz_vec_clear(tq, r*s + 1);
    _fmpz_vec_clear(tr, r*s + 1);

    fmpz_clear(cnum);
    fmpz_clear(cden);
    fmpz_clear(gg);
    fmpz_clear(hden);
    fmpz_clear(tden);
    fmpz_clear(hdenbar);

    return success;

fail:

    success = 0;
    goto cleanup;

}

int fmpq_poly_decompose(
    fmpq_poly_t g,
    fmpq_poly_t h,
    const fmpq_poly_t f,
    slong r,
    slong s)
{
    int success;

    if (r <= 1 || s <= 1 || r*s != fmpq_poly_degree(f))
        return 0;

    fmpz_one(g->den);
    fmpq_poly_fit_length(g, r + 1);
    _fmpq_poly_set_length(g, r + 1);

    fmpz_one(h->den);
    fmpq_poly_fit_length(h, s + 1);
    _fmpq_poly_set_length(h, s + 1);

    success = _fmpz_poly_decompose(g->coeffs, h->coeffs, f->coeffs, r, s);

    _fmpq_poly_normalise(g);
    _fmpq_poly_normalise(h);

    if (success)
    {
        fmpq_poly_scalar_div_fmpz(g, g, f->den);
printf("succeess !!!\n");
printf("h: "); fmpq_poly_print_pretty(h, "x"); printf("\n");
printf("g: "); fmpq_poly_print_pretty(g, "x"); printf("\n");

    }

    return success;
}


ex dcode_sDecompose(er e)
{
std::cout << "dcode_sDecompose: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sDecompose.get()));

    if (elength(e) != 2)
        return _handle_message_argx(e, 2);

    xfmpq_poly_t g, h, f;

    if (!eget_fmpq_poly(f.data, echild(e,1), echild(e,2)))
        return ecopy(e);

    fmpq_poly_decompose(g.data, h.data, f.data, 6, 3);

    return ecopy(e);
}






class xfmpz_mat_t {
public:
    fmpz_mat_t data;

    xfmpz_mat_t(slong r, slong c)
    {
        fmpz_mat_init(data, r, c);
    }

    ~xfmpz_mat_t()
    {
        fmpz_mat_clear(data);
    }
};

slong arb_wprec(const arb_t x)
{
    return std::max(slong(FLINT_BITS), arb_rel_accuracy_bits(x));
}


static bool root_approx(fmpz_poly_t f, const arb_t x, slong n)
{
    assert(n >= 1);

    slong try_count = 0;
    slong wp = arb_wprec(x);
    xfmpz_t acc;
    std::vector<xarb_t> v(n + 1);
    xfmpz_mat_t m(n + 1, n + 2);
    fmpz_lll_t fl;
    fmpz_lll_context_init_default(fl);
    fmpz_poly_fit_length(f, n + 1);

    arb_one(v[0].data);
    arb_set(v[1].data, x);
    fmpz_set(acc.data, MAG_EXPREF(arb_radref(v[1].data)));
    for (slong i = 2; i <= n; i++)
    {
        arb_mul(v[i].data, v[i - 1].data, x, wp + 2);
        if (fmpz_cmp(acc.data, MAG_EXPREF(arb_radref(v[i].data))) < 0)
            fmpz_set(acc.data, MAG_EXPREF(arb_radref(v[i].data)));
    }

    fmpz_neg(acc.data, acc.data);
    for (slong i = 0; i <= n; i++)
        arb_mul_2exp_fmpz(v[i].data, v[i].data, acc.data);

    for (try_count = 1; try_count <= 1; try_count++)
    {
        for (slong i = 0; i <= n; i++)
        {
            for (slong j = 0; j <= n; j++)
                fmpz_zero(fmpz_mat_entry(m.data, i, i));
            fmpz_one(fmpz_mat_entry(m.data, i, i));
            // TODO: fail if this is going to throw
            arf_get_fmpz(fmpz_mat_entry(m.data, i, n + 1), arb_midref(v[i].data), ARF_RND_NEAR);
        }

//printf("m: "); fmpz_mat_print_pretty(m.data); printf("\n");

        fmpz_lll(m.data, NULL, fl);

//printf("m: "); fmpz_mat_print_pretty(m.data); printf("\n");

        for (slong i = n; i >= 0; i--)
            fmpz_swap(f->coeffs + i, fmpz_mat_entry(m.data, 0, i));
        f->length = n + 1;
        _fmpz_poly_normalise(f);

        if (fmpz_poly_degree(f) > 0)
            return true;

        for (slong i = 0; i <= n; i++)
            arb_mul_2exp_si(v[i].data, v[i].data, 12*try_count);
    }

    return false;
}

static bool root_approx(fmpz_poly_t f, const acb_t z, slong n)
{
    assert(n >= 1);

    if (arb_is_zero(acb_imagref(z)))
        return root_approx(f, acb_realref(z), n);

    slong try_count = 0;
    slong wp = unary_wp(z);
    xfmpz_t acc;
    std::vector<xacb_t> v(n + 1);
    xfmpz_mat_t m(n + 1, n + 3);
    fmpz_lll_t fl;
    fmpz_lll_context_init_default(fl);
    fmpz_poly_fit_length(f, n + 1);

    acb_one(v[0].data);
    acb_set(v[1].data, z);
    fmpz_set(acc.data, MAG_EXPREF(arb_radref(acb_realref(v[1].data))));
    if (fmpz_cmp(acc.data, MAG_EXPREF(arb_radref(acb_imagref(v[1].data)))) < 0)
        fmpz_set(acc.data, MAG_EXPREF(arb_radref(acb_imagref(v[1].data))));

    for (slong i = 2; i <= n; i++)
    {
        acb_mul(v[i].data, v[i - 1].data, z, wp + 2);
        if (fmpz_cmp(acc.data, MAG_EXPREF(arb_radref(acb_realref(v[i].data)))) < 0)
            fmpz_set(acc.data, MAG_EXPREF(arb_radref(acb_realref(v[i].data))));
        if (fmpz_cmp(acc.data, MAG_EXPREF(arb_radref(acb_imagref(v[i].data)))) < 0)
            fmpz_set(acc.data, MAG_EXPREF(arb_radref(acb_imagref(v[i].data))));
    }

    fmpz_neg(acc.data, acc.data);
    for (slong i = 0; i <= n; i++)
        acb_mul_2exp_fmpz(v[i].data, v[i].data, acc.data);

    for (try_count = 1; try_count <= 10; try_count++)
    {
        for (slong i = 0; i <= n; i++)
            acb_mul_2exp_si(v[i].data, v[i].data, 12*try_count);

        for (slong i = 0; i <= n; i++)
        {
            for (slong j = 0; j <= n; j++)
                fmpz_zero(fmpz_mat_entry(m.data, i, i));
            fmpz_one(fmpz_mat_entry(m.data, i, i));
            // TODO: fail if this is going to throw
            arf_get_fmpz(fmpz_mat_entry(m.data, i, n + 1), arb_midref(acb_realref(v[i].data)), ARF_RND_NEAR);
            arf_get_fmpz(fmpz_mat_entry(m.data, i, n + 2), arb_midref(acb_imagref(v[i].data)), ARF_RND_NEAR);
        }

        fmpz_lll(m.data, NULL, fl);

        for (slong i = n; i >= 0; i--)
            fmpz_swap(f->coeffs + i, fmpz_mat_entry(m.data, 0, i));
        f->length = n + 1;
        _fmpz_poly_normalise(f);

        if (fmpz_poly_degree(f) > 0)
            return true;
    }

    return false;
}


static double weight(fmpz_poly_t f)
{
    double w = f->length;
    for (slong i = 0; i < f->length; i++)
        w += fmpz_bits(f->coeffs + i);
    return w;
}

static bool find_best_poly(qbarelem & s, const acb_t z)
{
    qbarelem r;
    bool success = false;
    xfmpz_poly_t f;

    double best_weight = unary_wp(z);
    for (slong n = 1; n <= 20; n += std::max(n, WORD(8))/8)
    {
//std::cout << "find_best_poly n = " << n << std::endl;
        if (root_approx(f.data, z, n) && r.set_contains(f.data, z))
        {
            double new_weight = weight(r.minpoly.data);
            if (new_weight < best_weight)
            {
                success = true;
                s.swap(r);
                if (new_weight < 0.75 * best_weight)
                    return true;
                best_weight = new_weight;
            }
        }
    }

    return success;
}



ex dcode_sRootApproximant(er e)
{
//std::cout << "dcode_sRootApproximant: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sRootApproximant.get()));

    if (elength(e) < 1 || elength(e) > 3)
        return _handle_message_argb(e, (1 << 0) + (3 << 8));

    xacb_t x;
    if (eget_acb(x.data, echild(e,1)))
    {
        qbarelem r;
        acb_canonicalize(x.data, FLINT_BITS);
        if (find_best_poly(r, x.data))
        {
            newton_trim(r.minpoly.data, r.location.data, r.epsilon.data);
            return r.get_ex();
        }
    }

    return ecopy(e);
}

ex ncode_sRoot(er e, slong prec)
{
    if (!ehas_head_sym(e, gs.sym_sRoot.get()))
        return ecopy(e);

    qbarelem r;
    if (!r.set_ex(e))
        return ecopy(e);

    while (acb_rel_accuracy_bits(r.location.data) <= prec + 5)
        newton_step(r.minpoly.data, r.location.data, r.epsilon.data);

    return acb_to_ex(r.location.data);
}



ex dcode_sRoots(er e)
{
//std::cout << "dcode_sRoots: " << ex_tostring_full(e) << std::endl;
    assert(ehas_head_sym(e, gs.sym_sRoots.get()));

    if (elength(e) != 2)
        return _handle_message_argx(e, 2);

    uex var(emake_node(gs.sym_sList.copy(), ecopychild(e,2)));
    uex eq(ecopychild(e,1));
    if (ehas_head_sym_length(echild(e,1), gs.sym_sEqual.get(), 2))
    {
        eq.setnz(emake_node(gs.sym_sMinus.copy(), ecopychild(e,1,1), ecopychild(e,1,2)));
    }

    poly p(1);
    if (!ex_to_polynomial(p, eq.get(), var.get()))
    {
        _gen_message(gs.sym_sRoots.get(), "neq", "`1` is expected to be a polynomial in the variable `2`.", eq.copy(), var.copychild(1));
        return ecopy(e);
    }

    xfmpq_poly_t q;
    xfmpq_t ai;
    slong ei;

    for (ulong i = 0; i < p.coeffs.size(); i++)
    {
        if (!eget_fmpq(ai.data, p.coeffs[i].get()) ||
            COEFF_IS_MPZ(*p.exps[i].data) ||
            (ei = *p.exps[i].data, ei < 0))
        {
            return ecopy(e);
        }

        fmpq_poly_set_coeff_fmpq(q.data, ei, ai.data);
    }

    xfmpz_poly_factor_t fac;
    fmpz_poly_factor(fac.data, q.zpoly());

    std::vector<wex> l;
    std::vector<qbarelem> v;

    for (slong i = 0; i < fac.data->num; i++)
    {
        irred_fmpz_poly_roots(v, fac.data->p + i);
        for (ulong j = 0; j < v.size(); j++)
        {
            wex rt(v[j].get_ex());
            for (slong k = 0; k < fac.data->exp[i]; k++)
                l.push_back(wex(rt.copy()));
        }
    }

    return emake_node(gs.sym_sList.copy(), l);
}