#include "flintarb_wrappers.h"

void fmpq_div_si(fmpq_t res, const fmpq_t op, slong x)
{
    fmpz_t y;
    fmpz_init(y);
    fmpz_set_si(y, x);
    fmpq_div_fmpz(res, op, y);
    fmpz_clear(y);
}

bool fmpz_coprime(const fmpz_t a, const fmpz_t b)
{
    fmpz_t g;
    fmpz_init(g);
    fmpz_gcd(g, a, b);
    int res = fmpz_is_one(g);
    fmpz_clear(g);
    return !!res;
}

void fmpz_addmult(fmpz_t x, const fmpz_t a, const fmpz_t b, fmpz_t t)
{
    fmpz_mul(t, a, b);
    fmpz_add(x, x, t);
}

void fmpz_submult(fmpz_t x, const fmpz_t a, const fmpz_t b, fmpz_t t)
{
    fmpz_mul(t, a, b);
    fmpz_sub(x, x, t);
}

void fmpz_neg_inplace(fmpz_t x)
{
    fmpz X = *x;
    if (COEFF_IS_MPZ(X))
        COEFF_TO_PTR(X)->_mp_size *= -1;
    else
        *x = -X;
}

void fmpz_si_sub(fmpz_t x, slong a, const fmpz_t b)
{
    fmpz_sub_si(x, b, a);
    fmpz_neg_inplace(x);
}

#if 0
void fmpz_mpoly_init_set(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    A->coeffs = FLINT_ARRAY_ALLOC(B->length, fmpz);
    for (i = 0; i < B->length; i++)
        fmpz_init_set(A->coeffs + i, B->coeffs + i);

    A->exps = FLINT_ARRAY_ALLOC(B->exps, B->length)
    for (i = 0; i < B->length; i++)
        A->exps[i] = B->exps[i];

    A->length = B->length;
    A->alloc = B->length;
}

void fmpq_mpoly_factor_set(fmpq_mpoly_factor_t res, const fmpq_mpoly_factor_t fac, const fmpq_mpoly_ctx_t ctx)
{
    slong i;

    if (res == fac)
        return;

    fmpq_mpoly_factor_fit_length(res, fac->length, ctx);
    fmpq_set(res->content, fac->content);
    for (i = 0; i < fac->length; i++)
    {
        fmpq_mpoly_set(res->poly + i, fac->poly + i, ctx);
        fmpz_set(res->exp + i, fac->exp + i);
    }
    res->length = fac->length;
}


int nmod_mpoly_is_ui(const nmod_mpoly_t a, const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_is_nmod(a, ctx);
}

void fmpz_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e)
{
    printf("fmpz_pow_fmpz called\n");
    flint_abort();
}

int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e)
{

    if (fmpq_is_zero(b))
    {
        int sgn = fmpz_sgn(e);
        if (sgn > 0)
        {
            fmpq_zero(a);
        }
        else if (sgn == 0)
        {
            fmpq_one(a);
        }
        else
        {
            flint_throw(FLINT_ERROR, "Division by zero in fmpq_pow_fmpz");
        }
        return 1;
    }
    else if (fmpz_is_one(fmpq_denref(b)) && fmpz_is_pm1(fmpq_numref(b)))
    {
        if (!fmpz_is_one(fmpq_numref(b)) && !fmpz_is_even(e))
        {
            fmpz_set_si(fmpq_numref(a), -1);
        }
        else
        {
            fmpz_one(fmpq_numref(a));
        }
        fmpz_one(fmpq_denref(a));
        return 1;
    }
    else if (fmpz_fits_si(e))
    {
        fmpq_pow_si(a, b, fmpz_get_si(e));
        return 1;
    }
    else
    {
        /* overflow */
        return 0;
    }
}

double fmpq_get_d(const fmpq_t a)
{
    double d;
    mpq_t z;
    flint_mpq_init_set_readonly(z, a);
    d = mpq_get_d(z);
    flint_mpq_clear_readonly(z);
    return d;
}

void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
{
    mpz_t F, S, G, H;
    mpz_init(F);
    mpz_init(S);
    mpz_init(G);
    mpz_init(H);
    fmpz_get_mpz(G, g);
    fmpz_get_mpz(H, h);
    mpz_cdiv_qr(F, S, G, H);
    fmpz_set_mpz(f, F);
    fmpz_set_mpz(s, S);
    mpz_clear(F);
    mpz_clear(S);
    mpz_clear(G);
    mpz_clear(H);
}

void fmpz_poly_inflate(fmpz_poly_t result,
                                      const fmpz_poly_t input, ulong inflation)
{
    if (input->length <= 1 || inflation == 1)
    {
        fmpz_poly_set(result, input);
    }
    else if (inflation == 0)
    {
        fmpz_t v;
        fmpz_init_set_ui(v, 1);
        fmpz_poly_evaluate_fmpz(v, input, v);
        fmpz_poly_zero(result);
        fmpz_poly_set_coeff_fmpz(result, 0, v);
        fmpz_clear(v);
    }
    else
    {
        slong i, j, res_length = (input->length - 1) * inflation + 1;

        fmpz_poly_fit_length(result, res_length);

        for (i = input->length - 1; i > 0; i--)
        {
            fmpz_set(result->coeffs + (i * inflation), input->coeffs + i);
            for (j = i * inflation - 1; j > (i - 1) * inflation; j--)
                fmpz_zero(result->coeffs + j);
        }
        fmpz_set(result->coeffs + 0, input->coeffs + 0);
        result->length = res_length;
    }
}


int _fmpq_cmp_fmpz(const fmpz_t p, const fmpz_t q, const fmpz_t r)
{
    int s1, s2, res;
    flint_bitcnt_t bp, bq, br, bs;
    fmpz_t u;

    if (fmpz_is_one(q))
        return fmpz_cmp(p, r);

    s1 = fmpz_sgn(p);
    s2 = fmpz_sgn(r);

    if (s1 != s2)
        return s1 < s2 ? -1 : 1;

    if (s1 == 0)
        return -s2;

    if (s2 == 0)
        return -s1;

    bp = fmpz_bits(p);
    bq = fmpz_bits(q);
    br = fmpz_bits(r);
    bs = 1;

    if (bp + bs + 1 < br + bq)
        return -s1;

    if (bp + bs > br + bq + 1)
        return s1;

    fmpz_init(u);

    fmpz_mul(u, q, r);

    res = fmpz_cmp(p, u);

    fmpz_clear(u);

    return res;
}


int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y)
{
    return _fmpq_cmp_fmpz(fmpq_numref(x), fmpq_denref(x), y);
}



/*
    after initialization one may:

        (1) set length to < 0: terms should not be stored, or
        (2) set want_alt_sum to +-1: terms should be add/sub to alt_sum, or
        (3) set limit to anything >= 0: stop generating terms after limit
        (4) set length to < 0 and set want_alt_sum to +-1: combination of (1) and (2)

    different combinations of these settings may or may not work
*/
void _fmpz_vector_init(_fmpz_vector_t v)
{
    v->array = NULL;
    v->length = 0;
    v->alloc = 0;
    v->limit = WORD_MAX;
    v->want_alt_sum = 0;
    fmpz_init(v->alt_sum);
}


void _fmpz_vector_clear(_fmpz_vector_t v)
{
    slong i;

    for (i = 0; i < v->alloc; i++)
        fmpz_clear(v->array + i);

    if (v->array)
        flint_free(v->array);

    fmpz_clear(v->alt_sum);
}


void _fmpz_vector_fit_length(_fmpz_vector_t v, slong len)
{
    if (len <= v->alloc)
        return;

    if (v->alloc > 0)
    {
        len = FLINT_MAX(len, v->alloc + v->alloc/2);

        v->array = (fmpz *) flint_realloc(v->array, len * sizeof(fmpz));
        FLINT_ASSERT(len > v->alloc);
        flint_mpn_zero((mp_ptr) (v->array + v->alloc), len - v->alloc);
    }
    else
    {
        v->array = (fmpz *) flint_calloc(len, sizeof(fmpz));
    }

    v->alloc = len;
}


void _fmpz_vector_push_back(_fmpz_vector_t v, const fmpz_t a)
{
    if (v->want_alt_sum)
    {
        v->want_alt_sum *= -1;
        if (v->want_alt_sum > 0)
            fmpz_sub(v->alt_sum, v->alt_sum, a);
        else
            fmpz_add(v->alt_sum, v->alt_sum, a);
    }

    if (v->length < 0)
        return;

    _fmpz_vector_fit_length(v, v->length + 1);
    fmpz_set(v->array + v->length, a);
    v->length++;
    FLINT_ASSERT(v->length <= v->limit);
}


void _fmpz_vector_push_back_zero(_fmpz_vector_t v)
{
    v->want_alt_sum *= -1;

    if (v->length < 0)
        return;

    _fmpz_vector_fit_length(v, v->length + 1);
    fmpz_zero(v->array + v->length);
    v->length++;
    FLINT_ASSERT(v->length <= v->limit);
}


void _fmpz_vector_append_ui(_fmpz_vector_t v, const ulong * a, slong n)
{
    slong i;

    if (v->want_alt_sum)
    {
        ulong hi = 0, lo = 0;
        for (i = 0; i + 2 <= n; i += 2)
        {
            add_ssaaaa(hi,lo, hi,lo, 0,a[i]);
            sub_ddmmss(hi,lo, hi,lo, 0,a[i + 1]);
        }

        if (i < n)
        {
            add_ssaaaa(hi,lo, hi,lo, 0,a[i]);
        }

        if (v->want_alt_sum < 0)
        {
            hi = -hi - (lo != 0);
            lo = -lo;
        }

        if (i < n)
        {
            v->want_alt_sum *= -1;
        }

        if (hi == 0)
        {
            fmpz_add_ui(v->alt_sum, v->alt_sum, lo);
        }
        else if (lo != 0 && hi == -UWORD(1))
        {
            fmpz_sub_ui(v->alt_sum, v->alt_sum, -lo);
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_set_signed_uiui(t, hi, lo);
            fmpz_add(v->alt_sum, v->alt_sum, t);
            fmpz_clear(t);
        }
    }

    if (v->length < 0)
        return;

    _fmpz_vector_fit_length(v, v->length + n);
    for (i = 0; i < n; i++)
        fmpz_set_ui(v->array + v->length + i, a[i]);
    v->length += n;

    FLINT_ASSERT(v->length <= v->limit);
}

#endif
