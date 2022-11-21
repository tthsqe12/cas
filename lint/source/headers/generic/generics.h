#pragma once

#include "fmpz.h"

template <typename Ring_T, typename Elem_T>
struct with {
    const Ring_T* _parent;
    Elem_T _elem;

    with(const Ring_T& p, const Elem_T& e) : _parent(&p), _elem(e) {}
    with(const Ring_T& p, Elem_T&& e) : _parent(&p), _elem(e) {}

    const Ring_T& parent() const {return *_parent;}

    Elem_T& elem() {return _elem;}
    const Elem_T& elem() const {return _elem;}
};

template <typename Ring_T, typename Elem_T>
void generic_pow_binexp(Ring_T& R, Elem_T& x, typename Elem_T::source_t a, ulimb e)
{
    if (e == 0)
    {
        one(R, x);
    }
    else if (e == 1)
    {
        set(R, x, a);
    }
    else if (e == 2)
    {
        sqr(R, x, a);
    }
    else
    {
        Elem_T t, s;
        set(R, t, a);        
        while ((e%2) == 0)
        {
            sqr(R, s, t);
            swap(R, t, s);
            e = e/2;
        }

        set(R, x, t);
        while ((e = e/2) > 0)
        {
            sqr(R, s, t);
            swap(R, t, s);
            if (e&1)
            {
                mul(R, s, x, t);
                swap(R, x, s);
            }
        }
    }
}

// x = a^e mod m
template <typename Ring_T, typename Elem_T>
void generic_powmod_binexp(Ring_T& R,
    Elem_T& x,
    typename Elem_T::source_t a,
    ulimb e,
    typename Elem_T::modulus_t m)
{
    if (e == 0)
    {
        one(R, x);
    }
    else if (e == 1)
    {
        mod(R, x, a, m);
    }
    else if (e == 2)
    {
        sqrmod(R, x, a, m);
    }
    else
    {
        Elem_T t, s;
        mod(R, t, a, m);

        while ((e%2) == 0)
        {
            sqrmod(R, s, t, m);
            swap(R, t, s);
            e = e/2;
        }

        mod(R, x, t, m);
        while ((e = e/2) > 0)
        {
            sqrmod(R, s, t, m);
            swap(R, t, s);
            if (e&1)
            {
                mulmod(R, s, x, t, m);
                swap(R, x, s);
            }
        }
    }
}

// x = a^e mod m
template <typename Ring_T, typename Elem_T>
void generic_powmod_binexp(Ring_T& R,
    Elem_T& x,
    typename Elem_T::source_t a,
    fmpzc e,
    typename Elem_T::modulus_t m)
{
    ulimb small_limb;
    const ulimb* limbs;
    ulimb bits;

    if (e.is_small())
    {
        small_limb = e.small();
        limbs = &small_limb;
        bits = nbits(limbs[0]);
    }
    else
    {
        limbs = e.limbs();
        bits = FLINT_BITS*(e.length()-1) + nbits(limbs[e.length()-1]);
    }

    if (e.is_zero())
    {
        one(R, x);
    }
    else if (e.is_one())
    {
        mod(R, x, a, m);
    }
    else if (e.data == 2)
    {
        sqrmod(R, x, a, m);
    }
    else
    {
        Elem_T t, s;
        mod(R, t, a, m);
        ulimb i = 0;

#define GET_BIT(i) ((limbs[(i)/FLINT_BITS] >> ((i)%LIMB_BITS))&1)
        while (GET_BIT(i) == 0) //while ((n%2) == 0)
        {
            sqrmod(R, s, t, m);
            swap(R, t, s);
            i++; //n = n/2;
        }

        mod(R, x, t, m);
        while (++i < bits) //while ((n = n/2) > 0)
        {
            sqrmod(R, s, t, m);
            swap(R, t, s);
            if (GET_BIT(i)) //if (n&1)
            {
                mulmod(R, s, x, t, m);
                swap(R, x, s);
            }
        }
#undef GET_BIT
    }
}


// f *= a^e  modulo units
template <typename Ring_T, typename PolyElem_T>
static void _generic_append_factor(
    Ring_T& R,
    typename PolyElem_T::product_t& f,
    typename PolyElem_T::source_t a,
    ulimb e)
{
    if (a.length() < 2)
        return;

    ulimb n = f.length();
    f.fit_alloc(n+1);
    f.set_length(n+1);
    f.exp(n) = e;
    unit_normalize(R, f.base(n), a);
}

// f *= a^e  modulo units, maintaining pairwise primeness from the i^th index
template <typename Ring_T, typename PolyElem_T>
static void _generic_append_factor_normalize(
    Ring_T& R,
    typename PolyElem_T::product_t& f,
    PolyElem_T& a,
    ulimb e,
    ulimb i)
{
    FLINT_ASSERT(e > 0);
    FLINT_ASSERT(a.length() > 0);

    PolyElem_T g, lbar, abar;

    while (i < f.length() && a.length() > 0)
    {
        // (g*lbar)^l[i].exp * (g*abar)^e =
        // (lbar)^l[i].exp * (g)^(l[i].exp+e) * (abar)^e
        gcdc(R, g, lbar, abar, f.base(i), a);
        if (g.length() < 2)
        {
            i++;
        }
        else if (lbar.length() < 2)
        {
            swap(R, a, abar);
            swap(R, f.base(i), g);
            f.exp(i) += e;
        }
        else if (abar.length() < 2)
        {
            swap(R, a, g);
            e += f.exp(i);
            unit_normalize(R, f.base(i), lbar);
        }
        else
        {
            swap(R, a, abar);
            _generic_append_factor_normalize<Ring_T, PolyElem_T>(R, f, g, e, i);
        }   
    }

    _generic_append_factor<Ring_T, PolyElem_T>(R, f, a, e);
}

// f *= a^e  modulo units, maintaining pairwise primeness from the i^th index  (or not)
// returned factors are squarefree and pairwise prime
template <typename Ring_T, typename PolyElem_T>
void generic_finite_field_poly_factor_squarefree(
    Ring_T& R,
    typename PolyElem_T::product_t& f,
    PolyElem_T a,
    ulimb e,
    ulimb i)  // 0 for normalization or UWORD_MAX for not
{
    if (a.length() < 2)
        return;

    PolyElem_T u, v, w, g;

again:

    if (a.length() == 2)
    {
        _generic_append_factor_normalize<Ring_T, PolyElem_T>(R, f, a, e, i);
        return;
    }

    derivative(R, g, a);
    gcdc(R, a, w, v, a, g);

    // if the characteristic doesn't fit, let the k <= p-2 condition be always true
    ulimb k, p = fmpz_abs_fits_ui(R.characteristic()) ? fmpz_abs_get_ui(R.characteristic()) : 1;
    for (k = 1; k <= p-2 && !(derivative(R, g, w), sub(R, u, v, g), is_zero(R, u)); k++)
    {
        gcdc(R, g, w, v, w, u);
        _generic_append_factor_normalize<Ring_T, PolyElem_T>(R, f, g, e*k, i);
        divexact(R, u, a, w);
        swap(R, a, u);
    }
    _generic_append_factor_normalize<Ring_T, PolyElem_T>(R, f, w, e*k, i);

    FLINT_ASSERT(!is_zero(R, a));
    FLINT_ASSERT((derivative(R, u, a), is_zero(R, u)));
    if (a.length() < 2)
        return;

    // a^(1/p) must be factored and appended with normalization
    FLINT_ASSERT(p > 1);
    deflate(R, a, p);
    e *= p;
    i = 0;
    goto again;
}



