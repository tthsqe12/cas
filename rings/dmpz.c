#if 0

#include "types.h"
#include "flintarb_wrappers.h"

/*
Integers:

        11xxxxxx    small negative
        11000000    large negative (this one is nullptr and not used)
-2^63:  10xxxxxx    large negative
2^63-1: 01xxxxxx    large positive
        00111111    small positive
        00xxxxxx    small positive

-(2^62-1) <= small <= 2^62-1

2^62 <= large positive

large negative <= -2^62
*/

typedef slong dmpz;
typedef dmpz dmpz_t[1];

struct dmpn_struct {
    slong length;
    slong alloc;
    mp_limb_t d[];
};


//s = +1 -> 01b
//s = -1 -> 10b
#define DPTR_TO_COEFF(x, s) (((ulong) (x) >> 2) | (WORD((3-s)/2) << (FLINT_BITS - 2)))
#define COEFF_TO_DPTR(x) (reinterpret_cast<dmpn_struct*> ((x) << 2))

void dmpz_init(dmpz_t x)
{
    *x = 0;
}

void dmpz_clear(dmpz_t x)
{
    if (COEFF_MIN <= *x && *x <= COEFF_MAX)
        return;

    free(COEFF_TO_DPTR(*x));
}


// assuming x is already big, ensure at least "fit" limbs allocated
dmpn_struct* _dmpn_fit_alloc(dmpn_struct* X, slong fit)
{
    if (X->alloc < fit)
    {
        X->alloc = std::max(2*X->alloc, fit);
        X = (dmpn_struct*) realloc(X, sizeof(dmpn_struct) + X->alloc*sizeof(mp_limb_t));
    }
    return X;
}

// ensure x is big and has at least "fit" limbs allocated
// destroyes the sign of x if x is big
dmpn_struct* _dmpz_promote(dmpz_t x, slong fit, int s)
{
    dmpn_struct* X;
    if (COEFF_MIN <= *x && *x <= COEFF_MAX)
    {
        X = (dmpn_struct*) malloc(sizeof(dmpn_struct) + fit*sizeof(mp_limb_t));
        X->length = 0;
        X->alloc = fit;
    }
    else
    {
        X = _dmpn_fit_alloc(COEFF_TO_DPTR(*x), fit);
    }
    *x = DPTR_TO_COEFF(X, s);
    return X;
}

// ensure x is small
void _dmpz_demote(dmpz_t x)
{
    if (COEFF_MIN <= *x && *x <= COEFF_MAX)
        return;

    free(COEFF_TO_DPTR(*x));
    *x = 0;
}

dmpz _dmpz_normalize_value(dmpn_struct* X, slong len, int s)
{
    assert(len > 0);

    while (len > 1 && X->d[len - 1] == 0)
        len--;

    if (len > 1 || len == 1 && X->d[0] > COEFF_MAX)
        return DPTR_TO_COEFF(X, s);

    dmpz res = s*slong(X->d[0]);
    free(X);
    return res;
}

void dmpz_zero(dmpz_t x)
{
    _dmpz_demote(x);
    *x = 0;
}

void dmpz_one(dmpz_t x)
{
    _dmpz_demote(x);
    *x = 1;
}

void dmpz_set_si(dmpz_t x, slong a)
{
    if (COEFF_MIN <= a && a <= COEFF_MAX)
    {
        _dmpz_demote(x);
        *x = a;
    }
    else
    {
        dmpn_struct* X = _dmpz_promote(x, 2, (a > 0) ? 1 : -1);
        X->d[0] = (a > 0) ? a : -a;
        X->length = 1;
    }
}

void dmpz_set(dmpz_t x, const dmpz_t a)
{
    dmpz sa = *a;
    if (COEFF_MIN <= sa && sa <= COEFF_MAX)
    {
        _dmpz_demote(x);
        *x = sa;
    }
    else
    {
        dmpn_struct* A = COEFF_TO_DPTR(sa);
        slong Alen = A->length;
        dmpn_struct* X = _dmpz_promote(x, Alen, sa > 0 ? 1 : -1);
        flint_mpn_copyi(X->d, A->d, Alen);
        X->length = Alen;
    }
}

void dmpz_neg(dmpz_t x, const dmpz_t a)
{
    dmpz sa = *a;
    if (COEFF_MIN <= sa && sa <= COEFF_MAX)
    {
        _dmpz_demote(x);
        *x = -sa;
    }
    else if (x == a)
    {
        *x ^= UWORD(3) << (FLINT_BITS - 2);
    }
    else
    {
        dmpn_struct* A = COEFF_TO_DPTR(sa);
        slong Alen = A->length;
        dmpn_struct* X = _dmpz_promote(x, Alen, sa > 0 ? -1 : 1);
        flint_mpn_copyi(X->d, A->d, Alen);
        X->length = Alen;
    }
}


// x = (A + b)*s
void _add_small(dmpz_t x, const dmpz_t a, mp_limb_t b, int s)
{
    dmpn_struct* A = COEFF_TO_DPTR(*a);
    slong Alen = A->length;
    dmpn_struct* X = _dmpz_promote(x, Alen + 1, s);
    A = COEFF_TO_DPTR(*a);
    X->d[Alen] = mpn_add_1(X->d, A->d, Alen, b);
    X->length = Alen + X->d[Alen];
}

// x = (A + B)*s
void _add_large(dmpz_t x, const dmpz_t a, const dmpz_t b, int s)
{
    dmpn_struct* X;
    dmpn_struct* A = COEFF_TO_DPTR(*a);
    dmpn_struct* B = COEFF_TO_DPTR(*b);
    slong Alen = A->length;
    slong Blen = B->length;
    if (Alen < Blen)
    {
        X = _dmpz_promote(x, Blen + 1, s);
        A = COEFF_TO_DPTR(*a);
        B = COEFF_TO_DPTR(*b);
        X->d[Blen] = mpn_add(X->d, B->d, Blen, A->d, Alen);
        X->length = Blen + X->d[Blen];
    }
    else
    {
        X = _dmpz_promote(x, Alen + 1, s);
        A = COEFF_TO_DPTR(*a);
        B = COEFF_TO_DPTR(*b);
        X->d[Alen] = mpn_add(X->d, A->d, Alen, B->d, Blen);
        X->length = Alen + X->d[Alen];
    }
}

// x = (A - b)*s
void _sub_small(dmpz_t x, const dmpz_t a, mp_limb_t b, int s)
{
    dmpn_struct* A = COEFF_TO_DPTR(*a);
    slong Alen = A->length;
    dmpn_struct* X = _dmpz_promote(x, Alen, s);
    mpn_sub_1(X->d, A->d, Alen, b);
    *x = _dmpz_normalize_value(X, Alen, s);
}


// x = A - B
void _sub_large(dmpz_t x, const dmpz_t a, const dmpz_t b)
{
    dmpn_struct* X;
    dmpn_struct* A = COEFF_TO_DPTR(*a);
    dmpn_struct* B = COEFF_TO_DPTR(*b);
    slong Alen = A->length;
    slong Blen = B->length;
    if (Alen > Blen)
    {
        X = _dmpz_promote(x, Alen, +1);
        B = COEFF_TO_DPTR(*b);
        mpn_sub(X->d, A->d, Alen, B->d, Blen);
        *x = _dmpz_normalize_value(X, Alen, +1);
    }
    else if (Alen < Blen)
    {
        X = _dmpz_promote(x, Blen, -1);
        A = COEFF_TO_DPTR(*a);
        mpn_sub(X->d, B->d, Blen, A->d, Alen);
        *x = _dmpz_normalize_value(X, Blen, -1);
    }
    else
    {
        int cmp = mpn_cmp(A->d, B->d, Alen);
        if (cmp > 0)
        {
            X = _dmpz_promote(x, Alen, +1);
            mpn_sub_n(X->d, A->d, B->d, Alen);
            *x = _dmpz_normalize_value(X, Alen, +1);
        }
        else if (cmp < 0)
        {
            X = _dmpz_promote(x, Blen, -1);
            mpn_sub_n(X->d, B->d, A->d, Blen);
            *x = _dmpz_normalize_value(X, Blen, -1);
        }
        else
        {
            dmpz_zero(x);
        }
    }
}

void dmpz_add(dmpz_t x, const dmpz_t a, const dmpz_t b)
{
    if (*a > COEFF_MAX)
    {
        if (*b > COEFF_MAX)
            _add_large(x, a, b, +1);
        else if (*b < COEFF_MIN)
            _sub_large(x, a, b);
        else if (*b > 0)
            _add_small(x, a, *b, +1);
        else if (*b < 0)
            _sub_small(x, a, -*b, +1);
        else
            dmpz_set(x, a);
    }
    else if (*a < COEFF_MIN)
    {
        if (*b < COEFF_MIN)
            _add_large(x, a, b, -1);
        else if (*b > COEFF_MAX)
            _sub_large(x, b, a);
        else if (*b < 0)
            _add_small(x, a, *b, -1);
        else if (*b > 0)
            _sub_small(x, b, -*b, -1);
        else
            dmpz_set(x, a);
    }
    else if (*b > COEFF_MAX)
    {
        if (*a > 0)
            _add_small(x, b, *a, +1);
        else if (*a < 0)
            _sub_small(x, b, -*a, +1);
        else
            dmpz_set(x, b);
    }
    else if (*b < COEFF_MAX)
    {
        if (*a < 0)
            _add_small(x, b, *a, -1);
        else if (*a > 0)
            _sub_small(x, b, -*a, -1);
        else
            dmpz_set(x, b);
    }
    else
    {
        dmpz_set_si(x, *a + *b);
    }
}

void dmpz_sub(dmpz_t x, const dmpz_t a, const dmpz_t b)
{
    if (*a > COEFF_MAX)
    {
        if (*b < COEFF_MIN)
            _add_large(x, b, b, +1);
        else if (*b > COEFF_MAX)
            _sub_large(x, a, b);
        else if (*b < 0)
            _add_small(x, a, -*b, +1);
        else if (*b > 0)
            _sub_small(x, a, +*b, +1);
        else
            dmpz_set(x, a);
    }
    else if (*a < COEFF_MIN)
    {
        if (*b > COEFF_MAX)
            _add_large(x, a, b, -1);
        else if (*b < COEFF_MIN)
            _sub_large(x, b, a);
        else if (*b > 0)
            _add_small(x, a, -*b, -1);
        else if (*b < 0)
            _sub_small(x, a, *b, -1);
        else
            dmpz_set(x, a);
    }
    else if (*b > COEFF_MAX)
    {
        if (*a < 0)
            _add_small(x, b, *a, -1);
        else if (*a > 0)
            _sub_small(x, b, -*a, -1);
        else
            dmpz_neg(x, b);
    }
    else if (*b < COEFF_MIN)
    {
        if (*a > 0)
            _add_small(x, b, *a, +1);
        else if (*a < 0)
            _sub_small(x, b, -*a, +1);
        else
            dmpz_neg(x, b);
    }
    else
    {
        dmpz_set_si(x, *a - *b);
    }
}

mywrite(std::ostream& o, const dmpz_t a)
{
    int base = 10;

    if (COEFF_MIN <= *a && *a <= COEFF_MAX)
    {
        o << *a;
        return o;
    }
    else
    {
        TMP_INIT;
        TMP_START;

        dmpn_struct* A = COEFF_TO_DPTR(*a);
        if (*a < COEFF_MIN)
            o << "-";

        slong Alen = A->length;
        size_t alloc_size;
        (alloc_size, A->d, Alen, base);
        unsigned char* res_str = TMP_ARRAY_ALLOC(alloc_size, unsigned char);

        ad = TMP_ARRAY_ALLOC(Alen, mp_limb_t);
        MPN_COPY(ad, A->d, Alen);

        size_t str_size = mpn_get_str(res_str, base, xp, x_size);
        assert(str_size <= alloc_size);

        for (i = 0; i < str_size; i++)
            o << "0123456789abcdefghijklmnopqrstuvwxyz"[res_str[i]];

        TMP_END;
        return o;
    }
}

#endif