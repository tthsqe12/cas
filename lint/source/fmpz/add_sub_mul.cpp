#include "mpn-impl.h"
#include "fmpz.h"

void fmpz_abs(fmpz& x, fmpzc a)
{
    FLINT_ASSERT(fmpz_is_canonical(a));
    if (a.is_negative())
        fmpz_neg(x, a);
    else
        fmpz_set(x, a);
}

void fmpz_neg(fmpz& x)
{
    FLINT_ASSERT(fmpz_is_canonical(x));

    if (x.is_small())
        x.data = -x.data;
    else
        x.data ^= UWORD(3) << (FLINT_BITS - 2);

    FLINT_ASSERT(fmpz_is_canonical(x));
}

void fmpz_neg(fmpz& x, fmpzc a)
{
    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));

    if (a.is_small())
    {
        _fmpz_demote(x, -a.small());
    }
    else if (x.data != a.data)
    {
        fmpn_struct* A = a.ptr();
        slimb Alen = A->length;
        fmpn_struct* X = _fmpz_promote(x, Alen, !a.is_negative());
        MPN_COPY(X->limbs, A->limbs, Alen);
        X->length = Alen;
    }
    else
    {
        fmpz_neg(x);
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));
}


// x = (|A| + b)*s
static void _add_ui(fmpz& x, fmpzc a, ulimb b, int s)
{
    fmpn_struct* A = a.ptr();
    slimb Alen = A->length;
    if (x.data == a.data)
    {
        fmpn_struct* X = x.ptr();
        ulimb out = mpn_add_1(X->limbs, X->limbs, Alen, b);
        if (out != 0)
        {
            X = _fmpn_fit_alloc(X, Alen + 1);
            X->limbs[Alen] = out;
            X->length = Alen + 1;
        }
        x.data = _fmpz_from_ptr(X, s);
    }
    else
    {
        fmpn_struct* X = _fmpz_promote(x, Alen + 1, s);
        X->limbs[Alen] = mpn_add_1(X->limbs, A->limbs, Alen, b);
        X->length = Alen + X->limbs[Alen];
    }
}

// x = (|A| + |B|)*s
static void _add_large(fmpz& x, fmpzc a, fmpzc b, int s)
{
    fmpn_struct* X;
    const ulimb* Alimbs = a.limbs();
    const ulimb* Blimbs = b.limbs();
    slimb Alen = a.length();
    slimb Blen = b.length();
    bool Aaliased = x.data == a.data;
    bool Baliased = x.data == b.data;
    if (Alen < Blen)
    {
        X = _fmpz_promote(x, Blen + 1, s);
        if (Aaliased) Alimbs = X->limbs;
        if (Baliased) Blimbs = X->limbs;
        X->limbs[Blen] = mpn_add(X->limbs, Blimbs, Blen, Alimbs, Alen);
        X->length = Blen + X->limbs[Blen];
    }
    else
    {
        X = _fmpz_promote(x, Alen + 1, s);
        if (Aaliased) Alimbs = X->limbs;
        if (Baliased) Blimbs = X->limbs;
        X->limbs[Alen] = mpn_add(X->limbs, Alimbs, Alen, Blimbs, Blen);
        X->length = Alen + X->limbs[Alen];
    }
}

// x = (|A| - b)*s
static void _sub_ui(fmpz& x, fmpzc a, ulimb b, int s)
{
    fmpn_struct* A = a.ptr();
    slimb Alen = A->length;
    if (Alen > 2)
    {
        fmpn_struct* X;
        ulimb* src;
        if (x.data == a.data)
        {
            X = x.ptr();
            src = X->limbs;
        }
        else
        {
            X = _fmpz_promote(x, Alen, s);
            src = A->limbs;
        }
        mpn_sub_1(X->limbs, A->limbs, Alen, b);
        x.data = _fmpz_normalize_value(X, Alen, s);        
    }
    else if (Alen == 2)
    {
        ulimb a0 = A->limbs[0];
        ulimb a1 = A->limbs[1];
        SUB_DDMMSS(a1,a0, a1,a0, 0,b);
        if (s)
            fmpz_neg_uiui(x, a1,a0);
        else
            fmpz_set_uiui(x, a1,a0);
    }
    else
    {
        ulimb a0 = A->limbs[0];
        if (a0 < b)
            if (s)
                fmpz_set_ui(x, b - a0);
            else
                fmpz_neg_ui(x, b - a0);
        else
            if (s)
                fmpz_neg_ui(x, a0 - b);
            else
                fmpz_set_ui(x, a0 - b);
    }
}

// x = A - B
static void _sub_large(fmpz& x, fmpzc a, fmpzc b)
{
    fmpn_struct* X;
    slimb Alen = a.length();
    slimb Blen = b.length();
    const ulimb* Alimbs = a.limbs();
    const ulimb* Blimbs = b.limbs();
    bool Aaliased = x.data == a.data;
    bool Baliased = x.data == b.data;
    int cmp = mpn_cmp4(Alimbs, Alen, Blimbs, Blen);
    if (cmp > 0)
    {
        X = _fmpz_promote(x, Alen, 0);
        if (Aaliased) Alimbs = X->limbs;
        if (Baliased) Blimbs = X->limbs;
        mpn_sub(X->limbs, Alimbs, Alen, Blimbs, Blen);
        x.data = _fmpz_normalize_value(X, Alen, 0);
    }
    else if (cmp < 0)
    {
        X = _fmpz_promote(x, Blen, 1);
        if (Aaliased) Alimbs = X->limbs;
        if (Baliased) Blimbs = X->limbs;
        mpn_sub(X->limbs, Blimbs, Blen, Alimbs, Alen);
        x.data = _fmpz_normalize_value(X, Blen, 1);
    }
    else
    {
        fmpz_zero(x);
    }
}

void fmpz_add(fmpz& x, fmpzc a, fmpzc b)
{
    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));
    FLINT_ASSERT(fmpz_is_canonical(b));

    if (LIKELY(a.is_small()))
    {
        if (LIKELY(b.is_small()))
        {
            fmpz_set_si(x, a.small() + b.small());
        }
        else
        {
            int sb = b.is_negative();
            int sa = a.is_negative();
            ulimb ua = my_abs(a.small());
            if (UNLIKELY(ua == 0))
                fmpz_set(x, b);
            else if (sa == sb)
                _add_ui(x, b, ua, sb);
            else
                _sub_ui(x, b, ua, sb);
        }
    }
    else
    {
        int sb = b.is_negative();
        int sa = a.is_negative();

        if (LIKELY(b.is_small()))
        {
            ulimb ub = my_abs(b.small());

            if (UNLIKELY(ub == 0))
            {
                fmpz_set(x, a);
            }
            else if (sb == sa)
            {
                _add_ui(x, a, ub, sa);
            }
            else
            {
                _sub_ui(x, a, ub, sa);
            }
        }
        else if (sa == sb)
        {
            _add_large(x, a, b, sa);
        }
        else if (sa)
        {
            _sub_large(x, b, a);
        }
        else
        {
            _sub_large(x, a, b);
        }
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
}

void fmpz_sub(fmpz& x, fmpzc a, fmpzc b)
{
    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));
    FLINT_ASSERT(fmpz_is_canonical(b));

    if (LIKELY(a.is_small()))
    {
        if (LIKELY(b.is_small()))
        {
            fmpz_set_si(x, a.small() - b.small());
        }
        else
        {
            int sb = b.is_negative();
            int sa = a.is_negative();
            ulimb ua = my_abs(a.small());

            if (UNLIKELY(ua == 0))
                fmpz_neg(x, b);
            else if (sa == sb)
                _sub_ui(x, b, ua, !sb);
            else
                _add_ui(x, b, ua, !sb);
        }
    }
    else
    {
        int sb = b.is_negative();
        int sa = a.is_negative();

        if (LIKELY(b.is_small()))
        {
            ulimb ub = my_abs(b.small());

            if (UNLIKELY(ub == 0))
                fmpz_set(x, a);
            else if (sb == sa)
                _sub_ui(x, a, ub, sa);
            else
                _add_ui(x, a, ub, sa);
        }
        else if (sa != sb)
        {
            _add_large(x, a, b, sa);
        }
        else if (sa)
        {
            _sub_large(x, b, a);
        }
        else
        {
            _sub_large(x, a, b);
        }
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
}

void fmpz_add_ui(fmpz& x, fmpzc a, ulimb b)
{
    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));

    if (LIKELY(a.is_small()))
    {
        ulimb Ahi, Alo, hi, lo;
        Alo = a.small();
        Ahi = sign_extend(Alo);
        ADD_SSAAAA(hi, lo, Ahi, Alo, UWORD(0), b);
        fmpz_set_signed_uiui(x, hi, lo);
    }
    else if (a.is_negative())
    {
        _sub_ui(x, a, b, 1);
    }
    else
    {
        _add_ui(x, a, b, 0);
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
}

void fmpz_sub_ui(fmpz& x, fmpzc a, ulimb b)
{
    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));

    if (LIKELY(a.is_small()))
    {
        ulimb Ahi, Alo, hi, lo;
        Alo = a.small();
        Ahi = sign_extend(Alo);
        SUB_DDMMSS(hi, lo, Ahi, Alo, UWORD(0), b);
        fmpz_set_signed_uiui(x, hi, lo);
    }
    else if (a.is_negative())
    {
        _add_ui(x, a, b, 1);
    }
    else
    {
        _sub_ui(x, a, b, 0);
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
}

// x = |A|*|B|*s
static void _mul_ui(fmpz& x, fmpzc a, ulimb b, int s)
{
    fmpn_struct* A = a.ptr();
    fmpn_struct* X;
    slimb Alen = A->length;
    if (x.data != a.data)
    {
        X = _fmpz_promote(x, Alen + 1, s);
        X->limbs[Alen] = mpn_mul_1(X->limbs, A->limbs, Alen, b);
    }
    else
    {
        X = _fmpz_promote(x, Alen + 1, s);
        X->limbs[Alen] = mpn_mul_1(X->limbs, X->limbs, Alen, b);
    }

    //x.data = _fmpz_normalize_value(X, Alen + 1, s);
    X->length = Alen + (X->limbs[Alen] != 0);
}

// x = |A|*|B|*s
static void _mul_large(fmpz& x, fmpzc a, fmpzc b, int s)
{
    fmpn_struct* X;
    fmpn_struct* A = a.ptr();
    fmpn_struct* B = b.ptr();
    slimb Alen = A->length;
    slimb Blen = B->length;
    ulimb Xlen = Alen + Blen;

    if (A == B)
    {
        if (Alen >= SQR_TOOM4_THRESHOLD)
        {
            fmpz_ulimb_realizer xx(x, Alen + Alen, s);
            glb_mpn_ctx.my_mpn_sqr(xx, A->limbs, Alen);
            X = x.ptr();
        }
        else if (UNLIKELY(x.data == a.data))
        {
            tmp_allocator push;
            ulimb* aa = push.recursive_alloc<ulimb>(Alen);
            MPN_COPY(aa, A->limbs, Alen);
            X = _fmpz_fit_destroy(x, Alen + Alen, s);
            my_mpn_sqr_toom_ordered(X->limbs, aa, Alen);
        }
        else
        {
            X = _fmpz_fit_destroy(x, Alen + Alen, s);
            my_mpn_sqr_toom_ordered(X->limbs, A->limbs, Alen);
        }
    }
    else
    {
        if (Alen < Blen)
        {
            std::swap(Alen, Blen);
            std::swap(A, B);
        }

        ulimb* aa = A->limbs, * bb = B->limbs;

        if (Blen >= MUL_TOOM44_THRESHOLD)
        {
            fmpz_ulimb_realizer xx(x, Xlen, s);
            glb_mpn_ctx.my_mpn_mul(xx, aa, Alen, bb, Blen);
            X = x.ptr();
        }
        else if (UNLIKELY(x._ptr() == A))
        {
            tmp_allocator push;
            aa = push.recursive_alloc<ulimb>(Alen);
            MPN_COPY(aa, A->limbs, Alen);
            X = _fmpz_fit_destroy(x, Xlen, s);
            my_mpn_mul_toom_ordered(X->limbs, aa, Alen, bb, Blen);
        }
        else if (UNLIKELY(x._ptr() == B))
        {
            tmp_allocator push;
            bb = push.recursive_alloc<ulimb>(Blen);
            MPN_COPY(bb, B->limbs, Blen);
            X = _fmpz_fit_destroy(x, Xlen, s);
            my_mpn_mul_toom_ordered(X->limbs, aa, Alen, bb, Blen);
        }
        else
        {
            X = _fmpz_fit_destroy(x, Xlen, s);
            my_mpn_mul_toom_ordered(X->limbs, aa, Alen, bb, Blen);
        }
    }

//    x.data = _fmpz_normalize_value(X, Xlen, s);
    X->length = Xlen - (X->limbs[Xlen - 1] == 0);
}

void fmpz_mul(fmpz& x, fmpzc a, fmpzc b)
{
    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));
    FLINT_ASSERT(fmpz_is_canonical(b));

    if (LIKELY(a.is_small()))
    {
        if (LIKELY(b.is_small()))
        {
            ulimb p1, p0;
            SMUL_PPMM(p1, p0, a.small(), b.small());
            fmpz_set_signed_uiui(x, p1, p0);
        }
        else
        {
            int s = a.is_negative()^b.is_negative();

            if (UNLIKELY(a.is_zero()))
                fmpz_zero(x);
            else
                _mul_ui(x, b, my_abs(a.small()), s);
        }
    }
    else
    {
        int s = a.is_negative()^b.is_negative();

        if (b.is_small())
        {
            if (UNLIKELY(b.is_zero()))
                fmpz_zero(x);
            else
                _mul_ui(x, a, my_abs(b.small()), s);
        }
        else
        {
            _mul_large(x, a, b, s);
        }
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
}

void fmpz_mul_ui(fmpz& x, fmpzc a, ulimb b)
{
    FLINT_ASSERT(fmpz_is_canonical(x));
    FLINT_ASSERT(fmpz_is_canonical(a));

    if (UNLIKELY(b == 0))
    {
        fmpz_zero(x);
    }
    else if (a.is_small())
    {
        ulimb p1, p0;
        ulimb ua = my_abs(a.small());
        if (a.is_negative())
        {
            UMUL_PPMM(p1, p0, ua, b);
            fmpz_neg_uiui(x, p1, p0);
        }
        else
        {
            UMUL_PPMM(p1, p0, ua, b);
            fmpz_set_uiui(x, p1, p0);
        }
    }
    else
    {
        _mul_ui(x, a, b, a.is_negative());
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
}

void fmpz_mul_si(fmpz& x, fmpzc a, slimb b)
{
    FLINT_ASSERT(fmpz_is_canonical(a));

    if (a.is_small())
    {
        ulimb p1, p0;
        SMUL_PPMM(p1, p0, a.small(), b);
        fmpz_set_signed_uiui(x, p1, p0);
    }
    else if (UNLIKELY(b == 0))
    {
        fmpz_zero(x);
    }
    else
    {
        _mul_ui(x, a, my_abs(b), (a.is_negative())^(b < 0));
    }

    FLINT_ASSERT(fmpz_is_canonical(x));
}

