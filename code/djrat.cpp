#include "djrat.h"
#include <unistd.h>

#define BOOL_CHECK(x) (x)

template <class U, class X>
class factored_elem {
public:
    U unit;
    std::vector<X> bases;
    std::vector<fmpz> exps;
    slong length;

    void fit_length(slong alloc)
    {
        if (bases.size() < alloc)
            bases.resize(alloc);
        if (exps.size() < alloc)
            exps.resize(alloc);
    }
};

template <class RU, class RX>
class rfactored_elem {
public:
    RU unit_ring;
    RX terms_ring;

    typedef factored_elem<typename RU::elem_t, typename RX::elem_t> elem_t;

    void mul(elem_t & a, elem_t & b, elem_t & c);
};

template <class RU, class RX>
void rfactored_elem<RU, RX>::mul(
    typename rfactored_elem<RU, RX>::elem_t & a,
    typename rfactored_elem<RU, RX>::elem_t & b,
    typename rfactored_elem<RU, RX>::elem_t & c)
{
    printf("not implemented\n");
}


void fmpq_polyfactor_print_pretty(const fmpq_polyfactor & a, const fmpz_mpoly_ctx_t ctx)
{
fmpq_print(a.sign);
for (slong i = 0; i < a.length; i++)
{
flint_printf("*(");
fmpz_mpoly_print_pretty(a.base + i, NULL, ctx);
flint_printf(")^");
fmpz_print(a.exp + i);
}
}

void nmod_polyfactor_print_pretty(const nmod_polyfactor & a, const nmod_mpoly_ctx_t ctx)
{
flint_printf("%wu", a.sign);
for (slong i = 0; i < a.length; i++)
{
flint_printf("*(");
nmod_mpoly_print_pretty(a.base + i, NULL, ctx);
flint_printf(")^");
fmpz_print(a.exp + i);
}
}

void fmpq_polyfactor::reset(const fmpz_mpoly_ctx_t ctx)
{
	fmpq_one(sign);

	for (size_t i = 0; i < alloc; i++)
	{
		fmpz_mpoly_clear(base + i, ctx);
		fmpz_clear(exp + i);
	}

	if (alloc == 0)
	{
		assert(base == nullptr && exp == nullptr);
	}
	else
	{
		assert(base != nullptr && exp != nullptr);
		free(base);
		free(exp);
	}

	base = nullptr;
	exp = nullptr;
	length = 0;
	alloc = 0;
}
void nmod_polyfactor::reset(const nmod_mpoly_ctx_t ctx)
{
	sign = 1;

	for (size_t i = 0; i < alloc; i++)
	{
		nmod_mpoly_clear(base + i, ctx);
		fmpz_clear(exp + i);
	}

	if (alloc == 0)
	{
		assert(base == nullptr && exp == nullptr);
	}
	else
	{
		assert(base != nullptr && exp != nullptr);
		free(base);
		free(exp);
	}

	base = nullptr;
	exp = nullptr;
	length = 0;
	alloc = 0;
}
void fmpz_mod_polyfactor::reset(const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_one(sign);

	for (size_t i = 0; i < alloc; i++)
	{
        fmpz_mod_mpoly_clear(base + i, ctx);
		fmpz_clear(exp + i);
	}

	if (alloc == 0)
	{
		assert(base == nullptr && exp == nullptr);
	}
	else
	{
		assert(base != nullptr && exp != nullptr);
		free(base);
		free(exp);
	}

	base = nullptr;
	exp = nullptr;
	length = 0;
	alloc = 0;
}


fmpq_polyfactor::~fmpq_polyfactor()
{
	fmpq_clear(sign);

	for (slong i = 0; i < alloc; i++)
	{
		fmpz_mpoly_clear(base + i, NULL);
		fmpz_clear(exp + i);
	}

	if (alloc == 0)
	{
		assert(base == nullptr && exp == nullptr);
	}
	else
	{
		assert(base != nullptr && exp != nullptr);
		free(base);
		free(exp);
	}
}
nmod_polyfactor::~nmod_polyfactor()
{
	for (slong i = 0; i < alloc; i++)
	{
		nmod_mpoly_clear(base + i, NULL);
		fmpz_clear(exp + i);
	}

	if (alloc == 0)
	{
		assert(base == nullptr && exp == nullptr);
	}
	else
	{
		assert(base != nullptr && exp != nullptr);
		free(base);
		free(exp);
	}		
}
fmpz_mod_polyfactor::~fmpz_mod_polyfactor()
{
	fmpz_clear(sign);
	for (slong i = 0; i < alloc; i++)
	{
		fmpz_mod_mpoly_clear(base + i, NULL);
		fmpz_clear(exp + i);
	}

	if (alloc == 0)
	{
		assert(base == nullptr && exp == nullptr);
	}
	else
	{
		assert(base != nullptr && exp != nullptr);
		free(base);
		free(exp);
	}		
}



void fmpq_polyfactor::mul_mpoly(const fmpz_mpoly_t b, const fmpz_mpoly_ctx_t ctx)
{
	if (fmpz_mpoly_is_zero(b, ctx))
	{
        fmpq_zero(sign);
        length = 0;
		return;
	}
	if (!fmpz_mpoly_is_one(b, ctx))
	{
		set_length(length + 1, ctx);
		fmpz_mpoly_set(base + length - 1, b, ctx);
		fmpz_one(exp + length - 1);
	}
}
void nmod_polyfactor::mul_mpoly(const nmod_mpoly_t b, const nmod_mpoly_ctx_t ctx)
{
	if (b->length == 0)
	{
		sign = 0;
		length = 0;
		return;
	}
	sign = nmod_mul(sign, b->coeffs[0], ctx->mod);
	if (!nmod_mpoly_is_ui(b, ctx))
	{
		set_length(length + 1, ctx);
		nmod_mpoly_set(base + length - 1, b, ctx);
		nmod_mpoly_scalar_mul_ui(base + length - 1, base + length - 1, nmod_inv(b->coeffs[0], ctx->mod), ctx);
		fmpz_set_ui(exp + length - 1, 1);
	}
}
void fmpz_mod_polyfactor::mul_mpoly(const fmpz_mod_mpoly_t b, const fmpz_mod_mpoly_ctx_t ctx)
{
	if (b->length == 0)
	{
		fmpz_zero(sign);
		length = 0;
		return;
	}
	fmpz_mod_mul(sign, sign, b->coeffs + 0, ctx->ffinfo);
	if (!fmpz_mod_mpoly_is_fmpz(b, ctx))
	{
		set_length(length + 1, ctx);
		fmpz_mod_mpoly_make_monic(base + length - 1, b, ctx);
		fmpz_one(exp + length - 1);
	}
}

void fmpq_polyfactor::mul_mpoly_pow_fmpz(const fmpz_mpoly_t b, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
{
	xfmpq_t t;
	assert(!fmpz_mpoly_is_zero(b, ctx));
	if (!fmpz_mpoly_is_one(b, ctx))
	{
		set_length(length + 1, ctx);
		fmpz_mpoly_set(base + length - 1, b, ctx);
		fmpz_set(exp + length - 1, e);
	}
}

int nmod_pow_fmpz_checked(mp_limb_t * a, mp_limb_t b, const fmpz_t e, nmod_t ctx)
{
    *a = nmod_pow_fmpz(b, e, ctx);
    if (fmpz_sgn(e) >= 0)
        return 0;
    return n_gcdinv(a, *a, ctx.n) == 1;
}

void nmod_polyfactor::mul_mpoly_pow_fmpz(const nmod_mpoly_t b, const fmpz_t e, const nmod_mpoly_ctx_t ctx)
{
	mp_limb_t t;
	assert(!nmod_mpoly_is_zero(b, ctx));
	BOOL_CHECK(nmod_pow_fmpz_checked(&t, b->coeffs[0], e, ctx->mod));
	sign = nmod_mul(sign, t, ctx->mod);
	if (!nmod_mpoly_is_ui(b, ctx))
	{
		set_length(length + 1, ctx);
		nmod_mpoly_make_monic(base + length - 1, b, ctx);
		fmpz_set(exp + length - 1, e);
	}
}

void fmpz_mod_polyfactor::mul_mpoly_pow_fmpz(const fmpz_mod_mpoly_t b, const fmpz_t e, const fmpz_mod_mpoly_ctx_t ctx)
{
	xfmpz_t t;
	assert(!fmpz_mod_mpoly_is_zero(b, ctx));
	BOOL_CHECK(fmpz_mod_pow_fmpz(t.data, b->coeffs + 0, e, ctx->ffinfo));
	fmpz_mod_mul(sign, sign, t.data, ctx->ffinfo);
	if (!fmpz_mod_mpoly_is_fmpz(b, ctx))
	{
		set_length(length + 1, ctx);
		fmpz_mod_mpoly_make_monic(base + length - 1, b, ctx);
		fmpz_set(exp + length - 1, e);
	}
}


fmpq_polyfactor::fmpq_polyfactor(
    const fmpq_polyfactor & other,
    const fmpz_mpoly_ctx_t ctx) :
	base(nullptr),
	exp(nullptr),
	length(other.length),
	alloc(other.length)
{
	fmpq_init(sign);
	fmpq_set(sign, other.sign);
	if (length > 0)
	{
		base = (fmpz_mpoly_struct *) malloc(length*sizeof(fmpz_mpoly_struct));
		exp = (fmpz *) malloc(length*sizeof(fmpz));
		for (size_t i = 0; i < length; i++)
		{
			fmpz_init_set(exp + i, other.exp + i);
			fmpz_mpoly_init(base + i, ctx);
			fmpz_mpoly_set(base + i, other.base + i, ctx);
		}
	}
}
nmod_polyfactor::nmod_polyfactor(const nmod_polyfactor & other, const nmod_mpoly_ctx_t ctx) :
	base(nullptr),
	exp(nullptr),
	length(other.length),
	alloc(other.length)
{
	sign = other.sign;
	if (length > 0)
	{
		base = (nmod_mpoly_struct *) malloc(length*sizeof(nmod_mpoly_struct));
		exp = (fmpz *) malloc(length*sizeof(fmpz));
		for (size_t i = 0; i < length; i++)
		{
			fmpz_init_set(exp + i, other.exp + i);
			nmod_mpoly_init(base + i, ctx);
			nmod_mpoly_set(base + i, other.base + i, ctx);
		}
	}
}
fmpz_mod_polyfactor::fmpz_mod_polyfactor(const fmpz_mod_polyfactor & other, const fmpz_mod_mpoly_ctx_t ctx) :
	base(nullptr),
	exp(nullptr),
	length(other.length),
	alloc(other.length)
{
	fmpz_set(sign, other.sign);
	if (length > 0)
	{
		base = (fmpz_mod_mpoly_struct *) malloc(length*sizeof(fmpz_mod_mpoly_struct));
		exp = (fmpz *) malloc(length*sizeof(fmpz));
		for (size_t i = 0; i < length; i++)
		{
			fmpz_init_set(exp + i, other.exp + i);
			fmpz_mod_mpoly_init(base + i, ctx);
			fmpz_mod_mpoly_set(base + i, other.base + i, ctx);
		}
	}
}

bool fmpz_mpoly_is_unit_normal(const fmpz_mpoly_t a)
{
    if (a->length < 1)
        return false;

    if (fmpz_sgn(a->coeffs + 0) < 0)
        return false;

    xfmpz_t g(slong(0));
    for (slong i = 0; i < a->length; i++)
    {
        fmpz_gcd(g.data, g.data, a->coeffs + i);
        if (fmpz_is_one(g.data))
            return true;
    }

    return false;
}


bool fmpq_polyfactor::is_nice(const fmpz_mpoly_ctx_t ctx) const
{
	xfmpz_mpoly_t g(ctx);

	/* unit normal bases */
	for (size_t i = 0; i < length; i++)
	    if (!fmpz_mpoly_is_unit_normal(base + i))
		    return false;

	/* rel prime bases */
	for (size_t i = 0; i + 1 < length; i++)
	for (size_t j = i + 1; j < length; j++)
	{
		if (!fmpz_mpoly_gcd(g.data, base + i, base + j, ctx))
            return false;
		if (!fmpz_mpoly_is_one(g.data, ctx))
			return false;
	}

	return true;
}
bool nmod_polyfactor::is_nice(const nmod_mpoly_ctx_t ctx) const
{
	xnmod_mpoly_t g(ctx);

	/* unit normal bases */
	for (size_t i = 0; i < length; i++)
	    if (!nmod_mpoly_is_unit_normal(base + i))
		    return false;

	/* rel prime bases */
	for (size_t i = 0; i + 1 < length; i++)
	for (size_t j = i + 1; j < length; j++)
	{
		nmod_mpoly_gcd(g.data, base + i, base + j, ctx);
		if (!nmod_mpoly_is_one(g.data, ctx))
			return false;
	}

	return true;
}

bool fmpz_mod_polyfactor::is_nice(const fmpz_mod_mpoly_ctx_t ctx) const
{
	xfmpz_mod_mpoly_t g(ctx);

	/* unit normal bases */
	for (size_t i = 0; i < length; i++)
	    if (!fmpz_mod_mpoly_is_unit_normal(base + i))
		    return false;

	/* rel prime bases */
	for (size_t i = 0; i + 1 < length; i++)
	for (size_t j = i + 1; j < length; j++)
	{
		fmpz_mod_mpoly_gcd(g.data, base + i, base + j, ctx);
		if (!fmpz_mod_mpoly_is_one(g.data, ctx))
			return false;
	}

	return true;
}

void fmpq_polyfactor::set_length(slong newlength, const fmpz_mpoly_ctx_t ctx)
{
	if (newlength > alloc)
	{
		slong newalloc = std::max(newlength, 1 + alloc + alloc/2);
        base = FLINT_ARRAY_REALLOC(base, newalloc, fmpz_mpoly_struct);
		exp = FLINT_ARRAY_REALLOC(exp, newalloc, fmpz);

		for (size_t i = alloc; i < newalloc; i++)
		{
			fmpz_mpoly_init(base + i, ctx);
			fmpz_init(exp + i);
		}
		alloc = newalloc;
	}
	length = newlength;
}
void nmod_polyfactor::set_length(slong newlength, const nmod_mpoly_ctx_t ctx)
{
	if (newlength > alloc)
	{
		slong newalloc = std::max(newlength, 1 + alloc + alloc/2);
		if (alloc == 0)
		{
			assert(base == nullptr && exp == nullptr);
			base = (nmod_mpoly_struct *) malloc(newalloc*sizeof(nmod_mpoly_struct));
			exp = (fmpz *) malloc(newalloc*sizeof(fmpz));
		}
		else
		{
			assert(base != nullptr && exp != nullptr);
			base = (nmod_mpoly_struct *) realloc(base, newalloc*sizeof(nmod_mpoly_struct));
			exp = (fmpz *) realloc(exp, newalloc*sizeof(fmpz));
		}

		for (size_t i = alloc; i < newalloc; i++)
		{
			nmod_mpoly_init(base + i, ctx);
			fmpz_init(exp + i);
		}
		alloc = newalloc;
	}
	length = newlength;
}
void fmpz_mod_polyfactor::set_length(slong newlength, const fmpz_mod_mpoly_ctx_t ctx)
{
	if (newlength > alloc)
	{
		slong newalloc = std::max(newlength, 1 + alloc + alloc/2);
		if (alloc == 0)
		{
			assert(base == nullptr && exp == nullptr);
			base = (fmpz_mod_mpoly_struct *) malloc(newalloc*sizeof(fmpz_mod_mpoly_struct));
			exp = (fmpz *) malloc(newalloc*sizeof(fmpz));
		}
		else
		{
			assert(base != nullptr && exp != nullptr);
			base = (fmpz_mod_mpoly_struct *) realloc(base, newalloc*sizeof(fmpz_mod_mpoly_struct));
			exp = (fmpz *) realloc(exp, newalloc*sizeof(fmpz));
		}

		for (size_t i = alloc; i < newalloc; i++)
		{
			fmpz_mod_mpoly_init(base + i, ctx);
			fmpz_init(exp + i);
		}
		alloc = newalloc;
	}
	length = newlength;
}

void fmpq_polyfactor::get_mpoly(fmpz_mpoly_t b, const fmpz_mpoly_ctx_t ctx)
{
	int success;
	xfmpz_mpoly_t t(ctx);
	xfmpz_t e;
	fmpz_mpoly_set_fmpz(b, fmpq_numref(sign), ctx);
	for (size_t i = 0; i < length; i++)
	{
		if (fmpz_sgn(exp + i) > 0)
		{
			fmpz_mpoly_pow_fmpz(t.data, base + i, exp + i, ctx);
			fmpz_mpoly_mul(b, b, t.data, ctx);
		}
		else if (fmpz_sgn(exp + i) < 0)
		{
			fmpz_neg(e.data, exp + i);
			fmpz_mpoly_pow_fmpz(t.data, base + i, e.data, ctx);
			success = fmpz_mpoly_divides(b, b, t.data, ctx);
			assert(success);
		}
	}
    fmpz_mpoly_scalar_divexact_fmpz(b, b, fmpq_denref(sign), ctx);
}
void nmod_polyfactor::get_mpoly(nmod_mpoly_t b, const nmod_mpoly_ctx_t ctx)
{
	int success;
	xnmod_mpoly_t t(ctx);
	xfmpz_t e;
	nmod_mpoly_set_ui(b, sign, ctx);
	for (size_t i = 0; i < length; i++)
	{
		if (fmpz_sgn(exp + i) > 0)
		{
			nmod_mpoly_pow_fmpz(t.data, base + i, exp + i, ctx);
			nmod_mpoly_mul(b, b, t.data, ctx);
		}
		else if (fmpz_sgn(exp + i) < 0)
		{
			fmpz_neg(e.data, exp + i);
			nmod_mpoly_pow_fmpz(t.data, base + i, e.data, ctx);
			success = nmod_mpoly_divides(b, b, t.data, ctx);
			assert(success);
		}
	}
}
void fmpz_mod_polyfactor::get_mpoly(fmpz_mod_mpoly_t b, const fmpz_mod_mpoly_ctx_t ctx)
{
	int success;
	xfmpz_mod_mpoly_t t(ctx);
	xfmpz_t e;
	fmpz_mod_mpoly_set_fmpz(b, sign, ctx);
	for (size_t i = 0; i < length; i++)
	{
		if (fmpz_sgn(exp + i) > 0)
		{
			fmpz_mod_mpoly_pow_fmpz(t.data, base + i, exp + i, ctx);
			fmpz_mod_mpoly_mul(b, b, t.data, ctx);
		}
		else if (fmpz_sgn(exp + i) < 0)
		{
			fmpz_neg(e.data, exp + i);
			fmpz_mod_mpoly_pow_fmpz(t.data, base + i, e.data, ctx);
			success = fmpz_mod_mpoly_divides(b, b, t.data, ctx);
			assert(success); // TODO
		}
	}
}
void fmpq_polyfactor::append_fix_units(
    const fmpq_polyfactor & other,
    const fmpz_mpoly_ctx_t ctx)
{
	xfmpq_t t;
    fmpq_mul(sign, sign, other.sign);
	for (size_t i = 0; i < other.length; i++)
	{
		assert(fmpz_mpoly_is_unit_normal(other.base + i));

		if (!fmpz_mpoly_is_one(other.base + i, ctx) &&
            !fmpz_is_zero(other.exp + i))
		{
			set_length(length + 1, ctx);
			fmpz_mpoly_set(base + length - 1, other.base + i, ctx);
			fmpz_set(exp + length - 1, other.exp + i);
		}
	}
}
void nmod_polyfactor::append_fix_units(const nmod_polyfactor & other, const nmod_mpoly_ctx_t ctx)
{
	mp_limb_t t;
    sign = nmod_mul(sign, other.sign, ctx->mod);
	for (size_t i = 0; i < other.length; i++)
	{
		assert(!nmod_mpoly_is_zero(other.base + i, ctx));
		if (other.base->coeffs[0] != 1)
		{
			BOOL_CHECK(nmod_pow_fmpz_checked(&t, other.base[i].coeffs[0], other.exp + i, ctx->mod));
			sign = nmod_mul(sign, t, ctx->mod);
		}

		if (!nmod_mpoly_is_one(other.base + i, ctx) &&
            !fmpz_is_zero(other.exp + i))
		{
			set_length(length + 1, ctx);
			nmod_mpoly_make_monic(base + length - 1, other.base + i, ctx);
			fmpz_set(exp + length - 1, other.exp + i);
		}
	}
}
void fmpz_mod_polyfactor::append_fix_units(const fmpz_mod_polyfactor & other, const fmpz_mod_mpoly_ctx_t ctx)
{
	xfmpz_t t;
    fmpz_mod_mul(sign, sign, other.sign, ctx->ffinfo);
	for (size_t i = 0; i < other.length; i++)
	{
		assert(!fmpz_mod_mpoly_is_zero(other.base + i, ctx));
		if (other.base->coeffs[0] != 1)
		{
			BOOL_CHECK(fmpz_mod_pow_fmpz(t.data, other.base[i].coeffs + 0, other.exp + i, ctx->ffinfo));
			fmpz_mod_mul(sign, sign, t.data, ctx->ffinfo);
		}

		if (!fmpz_mod_mpoly_is_one(other.base + i, ctx) &&
            !fmpz_is_zero(other.exp + i))
		{
			set_length(length + 1, ctx);
			fmpz_mod_mpoly_make_monic(base + length - 1, other.base + i, ctx);
			fmpz_set(exp + length - 1, other.exp + i);
		}
	}
}



void djrat::normalize()
{
	xfmpz_t g, u, v, e(slong(0));
try_again:
	for (size_t i = 0; i + 1 < length; i++)
	for (size_t j = i + 1; j < length; j++)
	{
		fmpz_gcd(g.data, base + i, base + j);
		if (!fmpz_is_one(g.data))
		{
			fmpz_add(e.data, exp + i, exp + j);
			if (!fmpz_is_zero(e.data))
			{
				set_length(length + 1);
				fmpz_set(base + length - 1, g.data);
				fmpz_set(exp + length - 1, e.data);
			}
			fmpz_divexact(base + i, base + i, g.data);
			fmpz_divexact(base + j, base + j, g.data);
			if (fmpz_is_one(base + j))
				delete_at(j);
			if (fmpz_is_one(base + i))
				delete_at(i);
			goto try_again;
		} 
	}
}

void fmpz_mpoly_unit_normalize(fmpz_t g, fmpz_mpoly_t a)
{
    _fmpz_vec_content(g, a->coeffs, a->length);
    if (fmpz_is_zero(g))
        return;

    if (fmpz_sgn(a->coeffs + 0) < 0)
        fmpz_neg(g, g);

    _fmpz_vec_scalar_divexact_fmpz(a->coeffs, a->coeffs, a->length, g);
}

void fmpq_polyfactor::normalize(const fmpz_mpoly_ctx_t ctx)
{
	xfmpz_mpoly_t g(ctx), u(ctx);
	xfmpz_t e;
    xfmpq_t t;

	for (size_t i = 0; i < length; i++)
    {
        fmpz_one(fmpq_denref(t.data));
        fmpz_mpoly_unit_normalize(fmpq_numref(t.data), base + i);
        BOOL_CHECK(fmpq_pow_fmpz(t.data, t.data, exp + i));
        fmpq_mul(sign, sign, t.data);
		if (fmpz_mpoly_is_one(base + i, ctx))
			delete_at(i, ctx);
    }

try_again:
	for (size_t i = 0; i + 1 < length; i++)
	for (size_t j = i + 1; j < length; j++)
	{
		assert(fmpz_mpoly_is_unit_normal(base + i));
		assert(fmpz_mpoly_is_unit_normal(base + j));
		BOOL_CHECK(fmpz_mpoly_gcd(g.data, base + i, base + j, ctx));
		if (!fmpz_mpoly_is_one(g.data, ctx))
		{
			assert(fmpz_mpoly_is_unit_normal(g.data));
			fmpz_add(e.data, exp + i, exp + j);
			if (!fmpz_is_zero(e.data))
			{
				set_length(length + 1, ctx);
				fmpz_mpoly_set(base + length - 1, g.data, ctx);
				fmpz_set(exp + length - 1, e.data);
			}
			fmpz_mpoly_divides(u.data, base + i, g.data, ctx);
			fmpz_mpoly_swap(base + i, u.data, ctx);
			fmpz_mpoly_divides(u.data, base + j, g.data, ctx);
			fmpz_mpoly_swap(base + j, u.data, ctx);
			assert(fmpz_mpoly_is_unit_normal(base + i));
			assert(fmpz_mpoly_is_unit_normal(base + j));
			if (fmpz_mpoly_is_one(base + j, ctx))
				delete_at(j, ctx);
			if (fmpz_mpoly_is_one(base + i, ctx))
				delete_at(i, ctx);
			goto try_again;
		} 
	}
}
void nmod_polyfactor::normalize(const nmod_mpoly_ctx_t ctx)
{
	xnmod_mpoly_t g(ctx), u(ctx);
	xfmpz_t e;
    mp_limb_t t;

	for (size_t i = 0; i < length; i++)
    {
		assert(!nmod_mpoly_is_zero(base + i, ctx));
        BOOL_CHECK(nmod_pow_fmpz_checked(&t, base[i].coeffs[0], exp + i, ctx->mod));
        sign = nmod_mul(sign, t, ctx->mod);
		nmod_mpoly_make_monic(base + i, base + i, ctx);
		assert(nmod_mpoly_is_unit_normal(base + i));
		if (nmod_mpoly_is_one(base + i, ctx))
			delete_at(i, ctx);
    }

try_again:
	for (size_t i = 0; i + 1 < length; i++)
	for (size_t j = i + 1; j < length; j++)
	{
		assert(nmod_mpoly_is_unit_normal(base + i));
		assert(nmod_mpoly_is_unit_normal(base + j));
		BOOL_CHECK(nmod_mpoly_gcd(g.data, base + i, base + j, ctx));
		if (!nmod_mpoly_is_one(g.data, ctx))
		{
			assert(nmod_mpoly_is_unit_normal(g.data));
			fmpz_add(e.data, exp + i, exp + j);
			if (!fmpz_is_zero(e.data))
			{
				set_length(length + 1, ctx);
				nmod_mpoly_set(base + length - 1, g.data, ctx);
				fmpz_set(exp + length - 1, e.data);
			}
			nmod_mpoly_divides(u.data, base + i, g.data, ctx);
			nmod_mpoly_swap(base + i, u.data, ctx);
			nmod_mpoly_divides(u.data, base + j, g.data, ctx);
			nmod_mpoly_swap(base + j, u.data, ctx);
			assert(nmod_mpoly_is_unit_normal(base + i));
			assert(nmod_mpoly_is_unit_normal(base + j));
			if (nmod_mpoly_is_one(base + j, ctx))
				delete_at(j, ctx);
			if (nmod_mpoly_is_one(base + i, ctx))
				delete_at(i, ctx);
			goto try_again;
		} 
	}
}
void fmpz_mod_polyfactor::normalize(const fmpz_mod_mpoly_ctx_t ctx)
{
	xfmpz_mod_mpoly_t g(ctx), u(ctx);
	xfmpz_t e;
    xfmpz_t t;

	for (size_t i = 0; i < length; i++)
    {
		assert(!fmpz_mod_mpoly_is_zero(base + i, ctx));
        BOOL_CHECK(fmpz_mod_pow_fmpz(t.data, base[i].coeffs + 0, exp + i, ctx->ffinfo));
        fmpz_mod_mul(sign, sign, t.data, ctx->ffinfo);
		fmpz_mod_mpoly_make_monic(base + i, base + i, ctx);
		assert(fmpz_mod_mpoly_is_unit_normal(base + i));
		if (fmpz_mod_mpoly_is_one(base + i, ctx))
			delete_at(i, ctx);
    }

try_again:
	for (size_t i = 0; i + 1 < length; i++)
	for (size_t j = i + 1; j < length; j++)
	{
		assert(fmpz_mod_mpoly_is_unit_normal(base + i));
		assert(fmpz_mod_mpoly_is_unit_normal(base + j));
		BOOL_CHECK(fmpz_mod_mpoly_gcd(g.data, base + i, base + j, ctx));
		if (!fmpz_mod_mpoly_is_one(g.data, ctx))
		{
			assert(fmpz_mod_mpoly_is_unit_normal(g.data));
			fmpz_add(e.data, exp + i, exp + j);
			if (!fmpz_is_zero(e.data))
			{
				set_length(length + 1, ctx);
				fmpz_mod_mpoly_set(base + length - 1, g.data, ctx);
				fmpz_set(exp + length - 1, e.data);
			}
			fmpz_mod_mpoly_divides(u.data, base + i, g.data, ctx);
			fmpz_mod_mpoly_swap(base + i, u.data, ctx);
			fmpz_mod_mpoly_divides(u.data, base + j, g.data, ctx);
			fmpz_mod_mpoly_swap(base + j, u.data, ctx);
			assert(fmpz_mod_mpoly_is_unit_normal(base + i));
			assert(fmpz_mod_mpoly_is_unit_normal(base + j));
			if (fmpz_mod_mpoly_is_one(base + j, ctx))
				delete_at(j, ctx);
			if (fmpz_mod_mpoly_is_one(base + i, ctx))
				delete_at(i, ctx);
			goto try_again;
		} 
	}
}


/*
	input a, g with g|a

	output l as a normalized factorization of s, where a = s * t, where s|g, and (t,g)=1
*/

bool coprime_fac(djrat &l, const fmpz_t a, const fmpz_t g)
{
//std::cout<<" in coprime_fac"<<std::endl;
//std::cout<<"a: "<< a.tostring()<<std::endl;
//std::cout<<"g: "<< g.tostring()<<std::endl;
	l.sign = 1;
	l.length = 0;
	xfmpz_t t(slong(1)), r(slong(1)), rprev(slong(1)), e(slong(1));
	assert(fmpz_divisible(a, g));
	fmpz_divexact(t.data, a, g);
	fmpz_gcd(r.data, t.data, g);
	fmpz_set(rprev.data, g);
	if (!fmpz_is_one(g) && !fmpz_is_zero(e.data))
	{
		l.set_length(l.length + 1);
		fmpz_set(l.base + l.length - 1, g);
		fmpz_set(l.exp + l.length - 1, e.data);
	}
	while (!fmpz_is_one(r.data))
	{
		if (fmpz_cmp(rprev.data, r.data) == 0)
		{
			fmpz_add_ui(l.exp + l.length - 1, l.exp + l.length - 1, 1);
			fmpz_divexact(t.data, t.data, r.data);
			fmpz_set(rprev.data, r.data);
			fmpz_gcd(r.data, t.data, r.data);
		}
		else
		{	
			l.set_length(l.length + 1);
			fmpz_set(l.base + l.length - 1, r.data);
			fmpz_set(l.exp + l.length - 1, e.data);
			fmpz_divexact(t.data,t.data,r.data);
			fmpz_set(rprev.data, r.data);
			fmpz_gcd(r.data, t.data, r.data);
		}
	}
	l.normalize();
	return true;
}


static bool coprime_fac(fmpq_polyfactor &l, fmpz_mpoly_t t, const fmpz_mpoly_t a, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx)
{
//std::cout<<" in coprime_fac"<<std::endl;
//std::cout<<"a: "<< a.tostring()<<std::endl;
//std::cout<<"g: "<< g.tostring()<<std::endl;
	xfmpz_mpoly_t r(ctx), rprev(ctx);
	int success;

    assert(fmpz_mpoly_is_unit_normal(a));
    assert(fmpz_mpoly_is_unit_normal(g));

	fmpz_mpoly_one(t, ctx);
	fmpz_mpoly_one(r.data, ctx);

	fmpz_mpoly_one(rprev.data, ctx);

	l.one();
	success = fmpz_mpoly_divides(t, a, g, ctx);
	assert(success);
	fmpz_mpoly_gcd(r.data, t, g, ctx);
	fmpz_mpoly_set(rprev.data, g, ctx);
	if (!fmpz_mpoly_is_one(g, ctx))
		l.push_back(g, 1, ctx);

	while (!fmpz_mpoly_is_one(r.data, ctx))
	{
		if (fmpz_mpoly_equal(rprev.data, r.data, ctx))
		{
			assert(l.length > 0);
			assert(fmpz_mpoly_equal(l.base + l.length - 1, r.data, ctx));
			fmpz_add_ui(l.exp + l.length - 1, l.exp + l.length - 1, 1);
		}
		else
		{
			l.push_back(r.data, 1, ctx);
		}
		success = fmpz_mpoly_divides(t, t, r.data, ctx);
		assert(success);
		fmpz_mpoly_swap(rprev.data, r.data, ctx);
		fmpz_mpoly_gcd(r.data, t, rprev.data, ctx);
	}
	l.normalize(ctx);
	return true;
}
static bool coprime_fac(nmod_polyfactor &l, nmod_mpoly_t t, const nmod_mpoly_t a, const nmod_mpoly_t g, const nmod_mpoly_ctx_t ctx)
{
//std::cout<<" in coprime_fac"<<std::endl;
//std::cout<<"a: "<< a.tostring()<<std::endl;
//std::cout<<"g: "<< g.tostring()<<std::endl;
	xnmod_mpoly_t r(ctx), rprev(ctx);
	int success;

	assert(nmod_mpoly_is_unit_normal(a));
	assert(nmod_mpoly_is_unit_normal(g));

	nmod_mpoly_one(t, ctx);
	nmod_mpoly_one(r.data, ctx);
	nmod_mpoly_one(rprev.data, ctx);

	l.one();
	success = nmod_mpoly_divides(t, a, g, ctx);
	assert(success);
	nmod_mpoly_gcd(r.data, t, g, ctx);
	nmod_mpoly_set(rprev.data, g, ctx);
	if (!nmod_mpoly_is_one(g, ctx))
		l.push_back(g, 1, ctx);

	while (!nmod_mpoly_is_one(r.data, ctx))
	{
		if (nmod_mpoly_equal(rprev.data, r.data, ctx))
		{
			assert(l.length > 0);
			assert(nmod_mpoly_equal(l.base + l.length - 1, r.data, ctx));
			fmpz_add_ui(l.exp + l.length - 1, l.exp + l.length - 1, 1);
		}
		else
		{
			l.push_back(r.data, 1, ctx);
		}
		success = nmod_mpoly_divides(t, t, r.data, ctx);
		assert(success);
		nmod_mpoly_swap(rprev.data, r.data, ctx);
		nmod_mpoly_gcd(r.data, t, rprev.data, ctx);
	}
	l.normalize(ctx);
	return true;
}
static bool coprime_fac(
    fmpz_mod_polyfactor &l,
    fmpz_mod_mpoly_t t,
    const fmpz_mod_mpoly_t a,
    const fmpz_mod_mpoly_t g,
    const fmpz_mod_mpoly_ctx_t ctx)
{
//std::cout<<" in coprime_fac"<<std::endl;
//std::cout<<"a: "<< a.tostring()<<std::endl;
//std::cout<<"g: "<< g.tostring()<<std::endl;
	xfmpz_mod_mpoly_t r(ctx), rprev(ctx);
	int success;

	assert(fmpz_mod_mpoly_is_unit_normal(a));
	assert(fmpz_mod_mpoly_is_unit_normal(g));

	fmpz_mod_mpoly_one(t, ctx);
	fmpz_mod_mpoly_one(r.data, ctx);
	fmpz_mod_mpoly_one(rprev.data, ctx);

	l.one();
	success = fmpz_mod_mpoly_divides(t, a, g, ctx);
	assert(success);
	fmpz_mod_mpoly_gcd(r.data, t, g, ctx);
	fmpz_mod_mpoly_set(rprev.data, g, ctx);
	if (!fmpz_mod_mpoly_is_one(g, ctx))
		l.push_back(g, 1, ctx);

	while (!fmpz_mod_mpoly_is_one(r.data, ctx))
	{
		if (fmpz_mod_mpoly_equal(rprev.data, r.data, ctx))
		{
			assert(l.length > 0);
			assert(fmpz_mod_mpoly_equal(l.base + l.length - 1, r.data, ctx));
			fmpz_add_ui(l.exp + l.length - 1, l.exp + l.length - 1, 1);
		}
		else
		{
			l.push_back(r.data, 1, ctx);
		}
		success = fmpz_mod_mpoly_divides(t, t, r.data, ctx);
		assert(success);
		fmpz_mod_mpoly_swap(rprev.data, r.data, ctx);
		fmpz_mod_mpoly_gcd(r.data, t, rprev.data, ctx);
	}
	l.normalize(ctx);
	return true;
}


bool mul(djrat & a, djrat & bo, djrat & co)
{
//std::cout<<"start mul: "<<std::endl;
//std::cout<<"****** bo : "<<bo.tostring()<<std::endl;
//std::cout<<"****** co : "<<co.tostring()<<std::endl;

	assert(bo.is_nice());
	assert(co.is_nice());

	a.sign = bo.sign * co.sign;
	a.length = 0;
	if (a.sign == 0)
		return true;

	djrat b(bo), c(co);
	xfmpz_t g, e, t, u, v;
	djrat l;

	a.length = 0;

	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		fmpz_gcd(g.data, b.base + i, c.base + j);
		if (!fmpz_is_one(g.data))
		{
			coprime_fac(l, b.base + i, g.data);
			for (size_t i1 = 0; i1 < l.length; i1++)
			{
				fmpz_mul(e.data, l.exp + i1, b.exp + i);
				a.push_back(l.base + i1, e.data);
			}
			fmpz_one(u.data);
		    for (size_t i1 = 0; i1 < l.length; i1++)
			{
				fmpz_pow_fmpz(t.data, l.base + i1, l.exp + i1);
				fmpz_mul(u.data, u.data, t.data);
			}
			fmpz_divexact(b.base + i, b.base + i, u.data);

			coprime_fac(l, c.base + j, g.data);
			for (size_t i1 = 0; i1 < l.length; i1++)
			{
				fmpz_mul(e.data, l.exp + i1, c.exp + j);
				a.push_back(l.base + i1, e.data);
			}
			fmpz_one(v.data);
		    for (size_t i1 = 0; i1 < l.length; i1++)
			{
				fmpz_pow_fmpz(t.data, l.base + i1, l.exp + i1);
				fmpz_mul(v.data, v.data, t.data);
			}
			fmpz_divexact(c.base + j, c.base + j, v.data);
			a.normalize();
		}
	}
//std::cout<<"**before a : "<<a.tostring()<<std::endl;

	a.append_fix_units(b);
	a.append_fix_units(c);

//std::cout<<"****** a : "<<a.tostring()<<std::endl;
//std::cout<<"****** a sign : "<<a.sign<<std::endl;

	assert(a.is_nice());

	return true;
}

bool fmpq_polyfactor_mul_fmpq(fmpq_polyfactor & a, const fmpq_polyfactor & bo, const fmpq_t co, const fmpz_mpoly_ctx_t ctx)
{
	fmpq_mul(a.sign, bo.sign, co);
	a.length = 0;
	if (fmpq_is_zero(a.sign))
		return true;

	a.append(bo, ctx);
	return true;
}
bool nmod_polyfactor_mul_nmod(nmod_polyfactor & a, const nmod_polyfactor & bo, mp_limb_t co, const nmod_mpoly_ctx_t ctx)
{
	a.sign = nmod_mul(bo.sign, co, ctx->mod);
	a.length = 0;
	if (a.sign == 0)
		return true;

	a.append(bo, ctx);
	return true;
}
bool fmpz_mod_polyfactor_mul_fmpz_mod(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & bo, const fmpz_t co, const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_mod_mul(a.sign, bo.sign, co, ctx->ffinfo);
	a.length = 0;
	if (a.sign == 0)
		return true;

	a.append(bo, ctx);
	return true;
}

bool fmpq_polyfactor_mul(fmpq_polyfactor & a, const fmpq_polyfactor & bo, const fmpq_polyfactor & co, const fmpz_mpoly_ctx_t ctx)
{
//std::cout<<"start mul: "<<std::endl;
//std::cout<<"****** bo : "<<bo.tostring()<<std::endl;
//std::cout<<"****** co : "<<co.tostring()<<std::endl;

	fmpq_mul(a.sign, bo.sign, co.sign);
	a.length = 0;
	if (fmpq_is_zero(a.sign))
		return true;

    if (bo.length == 1 && bo.base[0].length == 1 && fmpz_sgn(bo.exp + 0) >= 0 &&
        co.length == 1 && co.base[0].length == 1 && fmpz_sgn(co.exp + 0) >= 0)
    {
        a.set_length(1, ctx);
        xfmpz_mpoly_t t1(ctx), t2(ctx);
        fmpz_mpoly_pow_fmpz(t1.data, bo.base + 0, bo.exp + 0, ctx);
        fmpz_mpoly_pow_fmpz(t2.data, co.base + 0, co.exp + 0, ctx);
        fmpz_mpoly_mul(a.base + 0, t1.data, t2.data, ctx);
        fmpz_one(a.exp + 0);
        if (fmpz_mpoly_is_zero(a.base + 0, ctx))
            a.zero();
        return true;
    }

    if (bo.length == 0 &&
        co.length == 1 && co.base[0].length == 1 && fmpz_sgn(co.exp + 0) >= 0)
    {
        a.set_length(1, ctx);
        fmpz_mpoly_pow_fmpz(a.base + 0, co.base + 0, co.exp + 0, ctx);
        fmpz_one(a.exp + 0);
        if (fmpz_mpoly_is_zero(a.base + 0, ctx))
            a.zero();
        return true;
    }

    if (co.length == 0 &&
        bo.length == 1 && bo.base[0].length == 1 && fmpz_sgn(bo.exp + 0) >= 0)
    {
        a.set_length(1, ctx);
        fmpz_mpoly_pow_fmpz(a.base + 0, bo.base + 0, bo.exp + 0, ctx);
        fmpz_one(a.exp + 0);
        if (fmpz_mpoly_is_zero(a.base + 0, ctx))
            a.zero();
        return true;
    }


	fmpq_polyfactor l, b(bo, ctx), c(co, ctx);
	xfmpz_mpoly_t g(ctx), t(ctx);
	xfmpz_t e;

    fmpq_one(b.sign);
    fmpq_one(c.sign);

	b.normalize(ctx);
	c.normalize(ctx);
	assert(b.is_nice(ctx));
	assert(c.is_nice(ctx));

	a.length = 0;
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
        fmpz_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
        if (!fmpz_mpoly_is_one(g.data, ctx))
        {
            coprime_fac(l, t.data, b.base + i, g.data, ctx);
            for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, b.exp + i);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mpoly_swap(b.base + i, t.data, ctx);

			coprime_fac(l, t.data, c.base + j, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, c.exp + j);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mpoly_swap(c.base + j, t.data, ctx);

			a.normalize(ctx);
		}
	}
	a.append_fix_units(b, ctx);
	a.append_fix_units(c, ctx);
	assert(a.is_nice(ctx));
	return true;
}
bool nmod_polyfactor_mul(nmod_polyfactor & a, const nmod_polyfactor & bo, const nmod_polyfactor & co, const nmod_mpoly_ctx_t ctx)
{
//flint_printf("nmod_polyfactor mul called\n");
//flint_printf("bo: "); nmod_polyfactor_print_pretty(bo, ctx); flint_printf("\n");
//flint_printf("co: "); nmod_polyfactor_print_pretty(co, ctx); flint_printf("\n");


	a.sign = 1;
	a.length = 0;
	if (bo.sign == 0 || co.sign == 0)
    {
        a.sign = 0;
		return true;
    }

	nmod_polyfactor l, b(bo, ctx), c(co, ctx);
	xnmod_mpoly_t g(ctx), t(ctx);
	xfmpz_t e;

	b.normalize(ctx);
	c.normalize(ctx);
	assert(bo.is_nice(ctx));
	assert(co.is_nice(ctx));

	a.length = 0;
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		nmod_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
		if (!nmod_mpoly_is_one(g.data, ctx))
		{
			coprime_fac(l, t.data, b.base + i, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, b.exp + i);
				a.push_back(l.base + k, e.data, ctx);
			}
			nmod_mpoly_swap(b.base + i, t.data, ctx);

			coprime_fac(l, t.data, c.base + j, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, c.exp + j);
				a.push_back(l.base + k, e.data, ctx);
			}
			nmod_mpoly_swap(c.base + j, t.data, ctx);

			a.normalize(ctx);
		}
	}

	a.append_fix_units(b, ctx);
	a.append_fix_units(c, ctx);
	assert(a.is_nice(ctx));

	return true;
}
bool fmpz_mod_polyfactor_mul(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & bo, const fmpz_mod_polyfactor & co, const fmpz_mod_mpoly_ctx_t ctx)
{
//std::cout<<"start mul: "<<std::endl;
//std::cout<<"****** bo : "<<bo.tostring()<<std::endl;
//std::cout<<"****** co : "<<co.tostring()<<std::endl;

	fmpz_one(a.sign);
	a.length = 0;
	if (fmpz_is_zero(bo.sign) || fmpz_is_zero(co.sign))
    {
        fmpz_zero(a.sign);
		return true;
    }

	fmpz_mod_polyfactor l, b(bo, ctx), c(co, ctx);
	xfmpz_mod_mpoly_t g(ctx), t(ctx);
	xfmpz_t e;

	b.normalize(ctx);
	c.normalize(ctx);
	assert(bo.is_nice(ctx));
	assert(co.is_nice(ctx));

	a.length = 0;
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		fmpz_mod_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
		if (!fmpz_mod_mpoly_is_one(g.data, ctx))
		{
			coprime_fac(l, t.data, b.base + i, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, b.exp + i);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mod_mpoly_swap(b.base + i, t.data, ctx);

			coprime_fac(l, t.data, c.base + j, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, c.exp + j);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mod_mpoly_swap(c.base + j, t.data, ctx);

			a.normalize(ctx);
		}
	}
	a.append_fix_units(b, ctx);
	a.append_fix_units(c, ctx);
	assert(a.is_nice(ctx));
	return true;
}


bool div(djrat & a, djrat & b, djrat & c)
{
	return false;
}


bool fmpq_polyfactor_div(fmpq_polyfactor & a, const fmpq_polyfactor & bo, const fmpq_polyfactor & co, const fmpz_mpoly_ctx_t ctx)
{
	fmpq_div(a.sign, bo.sign, co.sign);
	a.length = 0;
	if (fmpq_is_zero(a.sign))
		return true;

    fmpq_polyfactor l, b(bo, ctx), c(co, ctx);
	xfmpz_mpoly_t g(ctx), t(ctx);
	xfmpz_t e;

	b.normalize(ctx);
	c.normalize(ctx);
	assert(bo.is_nice(ctx));
	assert(co.is_nice(ctx));

	for (size_t j = 0; j < c.length; j++)
		fmpz_neg(c.exp + j, c.exp + j);

	a.length = 0;
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		fmpz_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
		if (!fmpz_mpoly_is_one(g.data, ctx))
		{
			coprime_fac(l, t.data, b.base + i, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, b.exp + i);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mpoly_swap(b.base + i, t.data, ctx);

			coprime_fac(l, t.data, c.base + j, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, c.exp + j);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mpoly_swap(c.base + j, t.data, ctx);

			a.normalize(ctx);
		}
	}
	a.append_fix_units(b, ctx);
	a.append_fix_units(c, ctx);
	assert(a.is_nice(ctx));
	return true;
}
bool nmod_polyfactor_div(nmod_polyfactor & a, const nmod_polyfactor & bo, const nmod_polyfactor & co, const nmod_mpoly_ctx_t ctx)
{
	a.sign = nmod_mul(bo.sign, nmod_inv(co.sign, ctx->mod), ctx->mod);
	a.length = 0;
	if (a.sign == 0)
		return true;

	nmod_polyfactor l, b(bo, ctx), c(co, ctx);
	xnmod_mpoly_t g(ctx), t(ctx);
	xfmpz_t e;

	b.normalize(ctx);
	c.normalize(ctx);
	assert(bo.is_nice(ctx));
	assert(co.is_nice(ctx));

	for (size_t j = 0; j < c.length; j++)
		fmpz_neg(c.exp + j, c.exp + j);

	a.length = 0;
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		nmod_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
		if (!nmod_mpoly_is_one(g.data, ctx))
		{
			coprime_fac(l, t.data, b.base + i, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, b.exp + i);
				a.push_back(l.base + k, e.data, ctx);
			}
			nmod_mpoly_swap(b.base + i, t.data, ctx);

			coprime_fac(l, t.data, c.base + j, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, c.exp + j);
				a.push_back(l.base + k, e.data, ctx);
			}
			nmod_mpoly_swap(c.base + j, t.data, ctx);

			a.normalize(ctx);
		}
	}
	a.append_fix_units(b, ctx);
	a.append_fix_units(c, ctx);
	assert(a.is_nice(ctx));
	return true;
}

bool fmpz_mod_polyfactor_div(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & bo, const fmpz_mod_polyfactor & co, const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_mod_divides(a.sign, bo.sign, co.sign, ctx->ffinfo);
	a.length = 0;
	if (fmpz_is_zero(a.sign))
		return true;

	fmpz_mod_polyfactor l, b(bo, ctx), c(co, ctx);
	xfmpz_mod_mpoly_t g(ctx), t(ctx);
	xfmpz_t e;

	b.normalize(ctx);
	c.normalize(ctx);
	assert(bo.is_nice(ctx));
	assert(co.is_nice(ctx));

	for (size_t j = 0; j < c.length; j++)
		fmpz_neg(c.exp + j, c.exp + j);

	a.length = 0;
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		fmpz_mod_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
		if (!fmpz_mod_mpoly_is_one(g.data, ctx))
		{
			coprime_fac(l, t.data, b.base + i, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, b.exp + i);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mod_mpoly_swap(b.base + i, t.data, ctx);

			coprime_fac(l, t.data, c.base + j, g.data, ctx);
			for (size_t k = 0; k < l.length; k++)
			{
				fmpz_mul(e.data, l.exp + k, c.exp + j);
				a.push_back(l.base + k, e.data, ctx);
			}
			fmpz_mod_mpoly_swap(c.base + j, t.data, ctx);

			a.normalize(ctx);
		}
	}
	a.append_fix_units(b, ctx);
	a.append_fix_units(c, ctx);
	assert(a.is_nice(ctx));
	return true;
}



bool fmpq_polyfactor_equ(const fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpz_mpoly_ctx_t ctx)
{
	fmpq_polyfactor d;
	fmpq_polyfactor_div(d, a, b, ctx);
	return fmpq_is_one(d.sign) && d.length == 0;
}
bool nmod_polyfactor_equ(const nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_mpoly_ctx_t ctx)
{
	nmod_polyfactor d;
	nmod_polyfactor_div(d, a, b, ctx);
	return d.sign == 1 && d.length == 0;
}
bool fmpz_mod_polyfactor_equ(const fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_mod_polyfactor d;
	fmpz_mod_polyfactor_div(d, a, b, ctx);
	return fmpz_is_one(d.sign) && d.length == 0;
}


/*
	input a, b, c

	output a*(b,c), b/(b,c), c/(b,c)
*/
bool gcd_helper(djrat & a, djrat & b, djrat & c)
{
	djrat ct;
	xfmpz_t g, u, v, e(slong(0));

	assert(b.sign != 0);
	assert(c.sign != 0);

try_again: 
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		fmpz_gcd(g.data, b.base + i, c.base + j);
		if (!fmpz_is_one(g.data))
		{
			int cmp = fmpz_cmp(b.exp + i, c.exp + j);
			ct.one();
			fmpz_divexact(b.base + i, b.base + i, g.data);
			fmpz_divexact(c.base + j, c.base + j, g.data);
			if (cmp < 0)
			{
				a.push_back(g.data, b.exp + i);
				fmpz_sub(e.data, c.exp + j, b.exp + i);
				if (fmpz_is_one(b.base + i))
					b.delete_at(i);
				if (fmpz_is_one(c.base + j))
				{
					c.push_back(g.data, e.data);
					c.delete_at(j);
				}
				else
				{
					ct.push_back(g.data, e.data);
					ct.push_back(c.base + j, c.exp + j);
					c.delete_at(j);
					ct.normalize();
					c.append(ct);
				}
			}
			else if (cmp >= 0)
			{
				a.push_back(g.data, c.exp + j);
				fmpz_sub(e.data, b.exp + i, c.exp + j);
				if (fmpz_is_one(c.base + j))
					c.delete_at(j);
				if (fmpz_is_one(b.base + i))
				{
					b.push_back(g.data, e.data);
					b.delete_at(i);
				}
				else if (cmp > 0)
				{
					ct.push_back(g.data, e.data);
					ct.push_back(b.base + i, b.exp + i);
					b.delete_at(i);
					ct.normalize();
					b.append(ct);					
				}
			}
			goto try_again;	
		}
	}

	for (size_t i = 0; i < b.length; i++)
	{
		if (fmpz_sgn(b.exp + i) < 0)
		{
			a.push_back(b.base + i, b.exp + i);
			fmpz_neg(b.exp + i, b.exp + i);
			c.push_back(b.base + i, b.exp + i);
			b.delete_at(i);
			i--;
		}
	}

	for (size_t i = 0; i < c.length; i++)
	{
		if (fmpz_sgn(c.exp + i) < 0)
		{
			a.push_back(c.base + i, c.exp + i);
			fmpz_neg(c.exp + i, c.exp + i);
			b.push_back(c.base + i, c.exp + i);
			c.delete_at(i);
			i--;
		}
	}

	return true;
}

bool gcd_helper(fmpq_polyfactor & a, fmpq_polyfactor & b, fmpq_polyfactor & c, const fmpz_mpoly_ctx_t ctx)
{
	int success;
	xfmpz_mpoly_t g(ctx), u(ctx), v(ctx);
	xfmpz_t e;

	assert(!fmpq_is_zero(b.sign));
	assert(!fmpq_is_zero(c.sign));

    fmpq_gcd(a.sign, b.sign, c.sign);
    fmpq_div(b.sign, b.sign, a.sign);
    fmpq_div(c.sign, c.sign, a.sign);

    a.length = 0;

	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		fmpz_mpoly_gcd_cofactors(g.data, b.base + i, c.base + j, b.base + i, c.base + j, ctx);
		if (fmpz_mpoly_is_one(g.data, ctx))
            continue;

        // {b*g^be, c*g^ce}

        fmpz_sub(e.data, b.exp + i, c.exp + j);
        int sgn = fmpz_sgn(e.data);
        if (sgn >= 0)
        {
            // g^ce*{b*g^(be-ce), c}
            a.push_back(g.data, c.exp + j, ctx);
            if (sgn > 0)
                b.push_back(g.data, e.data, ctx);
        }
        else
        {
            // g^be*{b, c*g^(ce-be)}
            fmpz_neg(e.data, e.data);
            a.push_back(g.data, b.exp + j, ctx);
            c.push_back(g.data, e.data, ctx);            
        }

        if (fmpz_mpoly_is_one(c.base + j, ctx))
            c.delete_at(j, ctx);

        if (fmpz_mpoly_is_one(b.base + i, ctx))
            b.delete_at(i, ctx);
	}

	for (size_t i = 0; i < b.length; i++)
	{
		if (fmpz_sgn(b.exp + i) < 0)
		{
			a.push_back(b.base + i, b.exp + i, ctx);
			fmpz_neg(b.exp + i, b.exp + i);
			c.push_back(b.base + i, b.exp + i, ctx);
			b.delete_at(i, ctx);
			i--;
		}
	}

	for (size_t i = 0; i < c.length; i++)
	{
		if (fmpz_sgn(c.exp + i) < 0)
		{
			a.push_back(c.base + i, c.exp + i, ctx);
			fmpz_neg(c.exp + i, c.exp + i);
			b.push_back(c.base + i, c.exp + i, ctx);
			c.delete_at(i, ctx);
			i--;
		}
	}

	return true;
}
bool gcd_helper(nmod_polyfactor & a, nmod_polyfactor & b, nmod_polyfactor & c, const nmod_mpoly_ctx_t ctx)
{
	int success;
	nmod_polyfactor ct;
	xnmod_mpoly_t g(ctx), u(ctx), v(ctx);
	xfmpz_t e;

	assert(0 != b.sign);
	assert(0 != c.sign);

try_again: 
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		nmod_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
		if (!nmod_mpoly_is_one(g.data, ctx))
		{
			int cmp = fmpz_cmp(b.exp + i, c.exp + j);
			ct.one();
			success = nmod_mpoly_divides(b.base + i, b.base + i, g.data, ctx);
			assert(success);
			success = nmod_mpoly_divides(c.base + j, c.base + j, g.data, ctx);
			assert(success);
			if (cmp < 0)
			{
				a.push_back(g.data, b.exp + i, ctx);
				fmpz_sub(e.data, c.exp + j, b.exp + i);
				if (nmod_mpoly_is_one(b.base + i, ctx))
					b.delete_at(i, ctx);
				if (nmod_mpoly_is_one(c.base + j, ctx))
				{
					c.push_back(g.data, e.data, ctx);
					c.delete_at(j, ctx);
				}
				else
				{
					ct.push_back(g.data, e.data, ctx);
					ct.push_back(c.base + j, c.exp + j, ctx);
					c.delete_at(j, ctx);
					ct.normalize(ctx);
					c.append(ct, ctx);
				}
			}
			else if (cmp >= 0)
			{
				a.push_back(g.data, c.exp + j, ctx);
				fmpz_sub(e.data, b.exp + i, c.exp + j);
				if (nmod_mpoly_is_one(c.base + j, ctx))
					c.delete_at(j, ctx);
				if (nmod_mpoly_is_one(b.base + i, ctx))
				{
					b.push_back(g.data, e.data, ctx);
					b.delete_at(i, ctx);
				}
				else if (cmp > 0)
				{
					ct.push_back(g.data, e.data, ctx);
					ct.push_back(b.base + i, b.exp + i, ctx);
					b.delete_at(i, ctx);
					ct.normalize(ctx);
					b.append(ct, ctx);					
				}
			}
			goto try_again;	
		}
	}

	for (size_t i = 0; i < b.length; i++)
	{
		if (fmpz_sgn(b.exp + i) < 0)
		{
			a.push_back(b.base + i, b.exp + i, ctx);
			fmpz_neg(b.exp + i, b.exp + i);
			c.push_back(b.base + i, b.exp + i, ctx);
			b.delete_at(i, ctx);
			i--;
		}
	}

	for (size_t i = 0; i < c.length; i++)
	{
		if (fmpz_sgn(c.exp + i) < 0)
		{
			a.push_back(c.base + i, c.exp + i, ctx);
			fmpz_neg(c.exp + i, c.exp + i);
			b.push_back(c.base + i, c.exp + i, ctx);
			c.delete_at(i, ctx);
			i--;
		}
	}

	return true;
}
bool gcd_helper(fmpz_mod_polyfactor & a, fmpz_mod_polyfactor & b, fmpz_mod_polyfactor & c, const fmpz_mod_mpoly_ctx_t ctx)
{
	int success;
	fmpz_mod_polyfactor ct;
	xfmpz_mod_mpoly_t g(ctx), u(ctx), v(ctx);
	xfmpz_t e;

	assert(!fmpz_is_zero(b.sign));
	assert(!fmpz_is_zero(c.sign));

try_again: 
	for (size_t i = 0; i < b.length; i++)
	for (size_t j = 0; j < c.length; j++)
	{
		fmpz_mod_mpoly_gcd(g.data, b.base + i, c.base + j, ctx);
		if (!fmpz_mod_mpoly_is_one(g.data, ctx))
		{
			int cmp = fmpz_cmp(b.exp + i, c.exp + j);
			ct.one();
			success = fmpz_mod_mpoly_divides(b.base + i, b.base + i, g.data, ctx);
			assert(success);
			success = fmpz_mod_mpoly_divides(c.base + j, c.base + j, g.data, ctx);
			assert(success);
			if (cmp < 0)
			{
				a.push_back(g.data, b.exp + i, ctx);
				fmpz_sub(e.data, c.exp + j, b.exp + i);
				if (fmpz_mod_mpoly_is_one(b.base + i, ctx))
					b.delete_at(i, ctx);
				if (fmpz_mod_mpoly_is_one(c.base + j, ctx))
				{
					c.push_back(g.data, e.data, ctx);
					c.delete_at(j, ctx);
				}
				else
				{
					ct.push_back(g.data, e.data, ctx);
					ct.push_back(c.base + j, c.exp + j, ctx);
					c.delete_at(j, ctx);
					ct.normalize(ctx);
					c.append(ct, ctx);
				}
			}
			else if (cmp >= 0)
			{
				a.push_back(g.data, c.exp + j, ctx);
				fmpz_sub(e.data, b.exp + i, c.exp + j);
				if (fmpz_mod_mpoly_is_one(c.base + j, ctx))
					c.delete_at(j, ctx);
				if (fmpz_mod_mpoly_is_one(b.base + i, ctx))
				{
					b.push_back(g.data, e.data, ctx);
					b.delete_at(i, ctx);
				}
				else if (cmp > 0)
				{
					ct.push_back(g.data, e.data, ctx);
					ct.push_back(b.base + i, b.exp + i, ctx);
					b.delete_at(i, ctx);
					ct.normalize(ctx);
					b.append(ct, ctx);					
				}
			}
			goto try_again;	
		}
	}

	for (size_t i = 0; i < b.length; i++)
	{
		if (fmpz_sgn(b.exp + i) < 0)
		{
			a.push_back(b.base + i, b.exp + i, ctx);
			fmpz_neg(b.exp + i, b.exp + i);
			c.push_back(b.base + i, b.exp + i, ctx);
			b.delete_at(i, ctx);
			i--;
		}
	}

	for (size_t i = 0; i < c.length; i++)
	{
		if (fmpz_sgn(c.exp + i) < 0)
		{
			a.push_back(c.base + i, c.exp + i, ctx);
			fmpz_neg(c.exp + i, c.exp + i);
			b.push_back(c.base + i, c.exp + i, ctx);
			c.delete_at(i, ctx);
			i--;
		}
	}

	return true;
}


bool gcd(djrat & a, djrat & bo, djrat & co)
{
	djrat b(bo), c(co);
	assert(bo.is_nice());
	assert(co.is_nice());
	a.one();
	if (b.sign == 0)
	{
		if (c.sign == 0)
			a.sign =0;
		else
			a.append(c);
		return true;
	}
	if (c.sign == 0)
	{
		a.append(b);
		return true;
	}
	gcd_helper(a, b, c);
    a.normalize();
	assert(a.is_nice());
	return true;
}



bool fmpq_polyfactor_gcd(fmpq_polyfactor & a, const fmpq_polyfactor & bo, const fmpq_polyfactor & co, const fmpz_mpoly_ctx_t ctx)
{
	fmpq_polyfactor b(bo, ctx), c(co, ctx);
    b.normalize(ctx);
    c.normalize(ctx);

	fmpq_gcd(a.sign, b.sign, c.sign);
	if (fmpq_is_zero(b.sign))
	{
		if (fmpq_is_zero(c.sign))
			fmpq_zero(a.sign);
		else
			a.append(c, ctx);
        fmpq_abs(a.sign, a.sign);
		return true;
	}
	if (fmpq_is_zero(c.sign))
	{
		a.append(b, ctx);
        fmpq_abs(a.sign, a.sign);
		return true;
	}
	gcd_helper(a, b, c, ctx);

    a.normalize(ctx);

	return true;
}
bool nmod_polyfactor_gcd(nmod_polyfactor & a, const nmod_polyfactor & bo, const nmod_polyfactor & co, const nmod_mpoly_ctx_t ctx)
{
	nmod_polyfactor b(bo, ctx), c(co, ctx);
    b.normalize(ctx);
    c.normalize(ctx);
	a.sign = 1;
	if (b.sign == 0)
	{
		if (c.sign == 0)
			a.sign = 0;
		else
			a.append(c, ctx);
		return true;
	}
	if (c.sign == 0)
	{
		a.append(b, ctx);
		return true;
	}
	gcd_helper(a, b, c, ctx);
    a.normalize(ctx);
	return true;
}
bool fmpz_mod_polyfactor_gcd(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & bo, const fmpz_mod_polyfactor & co, const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_mod_polyfactor b(bo, ctx), c(co, ctx);
    b.normalize(ctx);
    c.normalize(ctx);
	fmpz_one(a.sign);
	if (fmpz_is_zero(b.sign))
	{
		if (fmpz_is_zero(c.sign))
			fmpz_zero(a.sign);
		else
			a.append(c, ctx);
		return true;
	}
	if (fmpz_is_zero(c.sign))
	{
		a.append(b, ctx);
		return true;
	}
	gcd_helper(a, b, c, ctx);
    a.normalize(ctx);
	return true;
}

bool fmpq_polyfactor_gcd_fmpq(fmpq_polyfactor & a, const fmpq_polyfactor & bo, const fmpq_t co, const fmpz_mpoly_ctx_t ctx)
{
	fmpq_polyfactor b(bo, ctx), c;
    b.normalize(ctx);
    fmpq_set(c.sign, co);
    c.length = 0;
	fmpq_gcd(a.sign, b.sign, c.sign);
	if (fmpq_is_zero(b.sign))
	{
		if (fmpq_is_zero(c.sign))
			fmpq_zero(a.sign);
		else
			a.append(c, ctx);
        fmpq_abs(a.sign, a.sign);
		return true;
	}
	if (fmpq_is_zero(c.sign))
	{
		a.append(b, ctx);
        fmpq_abs(a.sign, a.sign);
		return true;
	}
	gcd_helper(a, b, c, ctx);
    a.normalize(ctx);
	return true;
}
bool nmod_polyfactor_gcd_nmod(nmod_polyfactor & a, const nmod_polyfactor & bo, mp_limb_t co, const nmod_mpoly_ctx_t ctx)
{
	nmod_polyfactor b(bo, ctx), c;
    b.normalize(ctx);
	c.sign = co;
    c.length = 0;
	a.sign = 1;
	if (b.sign == 0)
	{
		if (c.sign == 0)
			a.sign = 0;
		else
			a.append(c, ctx);
		return true;
	}
	if (c.sign == 0)
	{
		a.append(b, ctx);
		return true;
	}
	gcd_helper(a, b, c, ctx);
    a.normalize(ctx);
	return true;
}
bool fmpz_mod_polyfactor_gcd_fmpz_mod(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & bo, const fmpz_t co, const fmpz_mod_mpoly_ctx_t ctx)
{
	fmpz_mod_polyfactor b(bo, ctx), c;
    b.normalize(ctx);
	fmpz_set(c.sign, co);
    c.length = 0;
	fmpz_one(a.sign);
	if (fmpz_is_zero(b.sign))
	{
		if (fmpz_is_zero(c.sign))
			fmpz_zero(a.sign);
		else
			a.append(c, ctx);
		return true;
	}
	if (fmpz_is_zero(c.sign))
	{
		a.append(b, ctx);
		return true;
	}
	gcd_helper(a, b, c, ctx);
    a.normalize(ctx);
	return true;
}

bool add(djrat & a, djrat & bo, djrat & co)
{
	xfmpz_t u, v, w;
	djrat b(bo), c(co);
	assert(bo.is_nice());
	assert(co.is_nice());
	a.one();
	if (b.sign == 0)
	{
		a.sign = c.sign;
		a.append(c);
		return true;
	}
	if (c.sign == 0)
	{
		a.sign = b.sign;
		a.append(b);
		return true;
	}
	gcd_helper(a, b, c);
	b.get_fmpz(v.data);
	c.get_fmpz(w.data);
	fmpz_add(u.data, v.data, w.data);
	a.mul_fmpz(u.data);
    a.normalize();
	assert(a.is_nice());
	return true;
}

bool fmpq_polyfactor_add(fmpq_polyfactor & a, const fmpq_polyfactor & bo, const fmpq_polyfactor & co, const fmpz_mpoly_ctx_t ctx)
{
    xfmpz_t t1, t2;
    fmpq_gcd_cofactors(a.sign, t1.data, t2.data, bo.sign, co.sign);

    if (bo.length == 1 && fmpz_is_one(bo.exp + 0) &&
        co.length == 1 && fmpz_is_one(co.exp + 0))
    {
        a.set_length(1, ctx);
        fmpz_mpoly_scalar_fmma(a.base + 0, bo.base + 0, t1.data, co.base + 0, t2.data, ctx);
        fmpz_one(a.exp + 0);
        if (fmpz_mpoly_is_zero(a.base + 0, ctx))
            a.zero();
        return true;
    }

    if (bo.length == 1 && fmpz_is_one(bo.exp + 0) &&
        co.length == 0)
    {
        a.set_length(1, ctx);
        fmpz_mpoly_scalar_mul_fmpz(a.base, bo.base + 0, t1.data, ctx);
        fmpz_mpoly_add_fmpz(a.base + 0, a.base + 0, t2.data, ctx);
        fmpz_one(a.exp + 0);
        if (fmpz_mpoly_is_zero(a.base + 0, ctx))
            a.zero();
        return true;
    }

    if (bo.length == 0 &&
        co.length == 1 && fmpz_is_one(co.exp + 0))
    {
        a.set_length(1, ctx);
        fmpz_mpoly_scalar_mul_fmpz(a.base + 0, co.base + 0, t2.data, ctx);
        fmpz_mpoly_add_fmpz(a.base + 0, a.base + 0, t1.data, ctx);
        fmpz_one(a.exp + 0);
        if (fmpz_mpoly_is_zero(a.base + 0, ctx))
            a.zero();
        return true;
    }

	xfmpz_mpoly_t u(ctx), v(ctx), w(ctx);
	fmpq_polyfactor b(bo, ctx), c(co, ctx);

	a.one();
	if (fmpq_is_zero(b.sign))
	{
		fmpq_set(a.sign, c.sign);
		a.append(c, ctx);
		return true;
	}
	if (fmpq_is_zero(c.sign))
	{
		fmpq_set(a.sign, b.sign);
		a.append(b, ctx);
		return true;
	}

	b.normalize(ctx);
	c.normalize(ctx);
	gcd_helper(a, b, c, ctx);

	b.get_mpoly(v.data, ctx);
	c.get_mpoly(w.data, ctx);

	fmpz_mpoly_add(u.data, v.data, w.data, ctx);

	a.mul_mpoly(u.data, ctx);
    a.normalize(ctx);

	return true;
}
bool nmod_polyfactor_add(nmod_polyfactor & a, const nmod_polyfactor & bo, const nmod_polyfactor & co, const nmod_mpoly_ctx_t ctx)
{
	xnmod_mpoly_t u(ctx), v(ctx), w(ctx);

    if (   bo.length == 1 && fmpz_is_one(bo.exp + 0)
        && co.length == 1 && fmpz_is_one(co.exp + 0))
    {
        nmod_mpoly_scalar_mul_ui(u.data, bo.base + 0, bo.sign, ctx);
        nmod_mpoly_scalar_mul_ui(v.data, co.base + 0, co.sign, ctx);
        a.set_length(1, ctx);
		a.sign = 1;
        nmod_mpoly_add(a.base + 0, u.data, v.data, ctx);
        fmpz_one(a.exp + 0);
        a.length = !nmod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

    if (   bo.length == 1 && fmpz_is_one(bo.exp + 0)
        && co.length == 0)
    {
        nmod_mpoly_scalar_mul_ui(u.data, bo.base + 0, bo.sign, ctx);
        a.set_length(1, ctx);
		a.sign = 1;
        nmod_mpoly_add_ui(a.base + 0, u.data, co.sign, ctx);
        fmpz_one(a.exp + 0);
        a.length = !nmod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

    if (   co.length == 1 && fmpz_is_one(co.exp + 0)
        && bo.length == 0)
    {
        nmod_mpoly_scalar_mul_ui(u.data, co.base + 0, co.sign, ctx);
        a.set_length(1, ctx);
        a.sign = 1;
        nmod_mpoly_add_ui(a.base + 0, u.data, bo.sign, ctx);
        fmpz_one(a.exp + 0);
        a.length = !nmod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

	nmod_polyfactor b(bo, ctx), c(co, ctx);

	a.one();
	if (b.sign == 0)
	{
		a.sign = c.sign;
		a.append(c, ctx);
		return true;
	}
	if (c.sign == 0)
	{
		a.sign = b.sign;
		a.append(b, ctx);
		return true;
	}

	b.normalize(ctx);
	c.normalize(ctx);
	gcd_helper(a, b, c, ctx);

	b.get_mpoly(v.data, ctx);
	c.get_mpoly(w.data, ctx);

	nmod_mpoly_add(u.data, v.data, w.data, ctx);

	a.mul_mpoly(u.data, ctx);
    a.normalize(ctx);
	return true;
}
bool fmpz_mod_polyfactor_add(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & bo, const fmpz_mod_polyfactor & co, const fmpz_mod_mpoly_ctx_t ctx)
{
	xfmpz_mod_mpoly_t u(ctx), v(ctx), w(ctx);

    if (   bo.length == 1 && fmpz_is_one(bo.exp + 0)
        && co.length == 1 && fmpz_is_one(co.exp + 0))
    {
        fmpz_mod_mpoly_scalar_mul_fmpz(u.data, bo.base + 0, bo.sign, ctx);
        fmpz_mod_mpoly_scalar_mul_fmpz(v.data, co.base + 0, co.sign, ctx);
        a.set_length(1, ctx);
		fmpz_one(a.sign);
        fmpz_mod_mpoly_add(a.base + 0, u.data, v.data, ctx);
        fmpz_one(a.exp + 0);
        a.length = !fmpz_mod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

    if (   bo.length == 1 && fmpz_is_one(bo.exp + 0)
        && co.length == 0)
    {
        fmpz_mod_mpoly_scalar_mul_fmpz(u.data, bo.base + 0, bo.sign, ctx);
        a.set_length(1, ctx);
		fmpz_one(a.sign);
        fmpz_mod_mpoly_add_fmpz(a.base + 0, u.data, co.sign, ctx);
        fmpz_one(a.exp + 0);
        a.length = !fmpz_mod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

    if (   co.length == 1 && fmpz_is_one(co.exp + 0)
        && bo.length == 0)
    {
        fmpz_mod_mpoly_scalar_mul_fmpz(u.data, co.base + 0, co.sign, ctx);
        a.set_length(1, ctx);
        fmpz_one(a.sign);
        fmpz_mod_mpoly_add_fmpz(a.base + 0, u.data, bo.sign, ctx);
        fmpz_one(a.exp + 0);
        a.length = !fmpz_mod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

	fmpz_mod_polyfactor b(bo, ctx), c(co, ctx);

	a.one();
	if (fmpz_is_zero(b.sign))
	{
		fmpz_set(a.sign, c.sign);
		a.append(c, ctx);
		return true;
	}
	if (fmpz_is_zero(c.sign))
	{
		fmpz_set(a.sign, b.sign);
		a.append(b, ctx);
		return true;
	}

	b.normalize(ctx);
	c.normalize(ctx);
	gcd_helper(a, b, c, ctx);
	b.get_mpoly(v.data, ctx);
	c.get_mpoly(w.data, ctx);
	fmpz_mod_mpoly_add(u.data, v.data, w.data, ctx);
	a.mul_mpoly(u.data, ctx);
    a.normalize(ctx);
	return true;
}



bool fmpq_polyfactor_add_fmpq(
    fmpq_polyfactor & a,
    const fmpq_polyfactor & bo,
    const fmpq_t co,
    const fmpz_mpoly_ctx_t ctx)
{
//std::cout << "fmpq_polyfactor_add_fmpq called" << std::endl;

//std::cout << "bo: " << ex_tostring(etor(ratpoly_get_ex<fmpq_ratpoly>(bo))) << std::endl;
/*
std::cout << "co: ";
fmpq_print(co);
std::cout << std::endl;
*/

    if (bo.length == 1 && fmpz_is_one(bo.exp + 0))
    {
//std::cout << "fmpq_polyfactor_add_fmpq case 1" << std::endl;

        xfmpz_t t1, t2;
        a.set_length(1, ctx);
        fmpq_gcd_cofactors(a.sign, t1.data, t2.data, bo.sign, co);
        fmpz_mpoly_scalar_mul_fmpz(a.base + 0, bo.base + 0, t1.data, ctx);
        fmpz_mpoly_add_fmpz(a.base + 0, a.base + 0, t2.data, ctx);
        fmpz_one(a.exp + 0);
        a.length = !fmpz_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

//std::cout << "fmpq_polyfactor_add_fmpq case 2" << std::endl;


	xfmpz_mpoly_t u(ctx), v(ctx), w(ctx);

	fmpq_polyfactor b(bo, ctx), c;
	fmpq_set(c.sign, co);
	c.length = 0;
	b.normalize(ctx);
	assert(b.is_nice(ctx));
	a.one();
	if (fmpq_is_zero(b.sign))
	{
		fmpq_set(a.sign, c.sign);
		a.append(c, ctx);
		return true;
	}
	if (fmpq_is_zero(c.sign))
	{
		fmpq_set(a.sign, b.sign);
		a.append(b, ctx);
		return true;
	}

	gcd_helper(a, b, c, ctx);

	b.get_mpoly(v.data, ctx);
	c.get_mpoly(w.data, ctx);
	fmpz_mpoly_add(u.data, v.data, w.data, ctx);
	a.mul_mpoly(u.data, ctx);
    a.normalize(ctx);
	assert(a.is_nice(ctx));

	return true;
}
bool nmod_polyfactor_add_nmod(nmod_polyfactor & a, const nmod_polyfactor & bo, mp_limb_t co, const nmod_mpoly_ctx_t ctx)
{
	xnmod_mpoly_t u(ctx), v(ctx), w(ctx);

    if (bo.length == 1 && fmpz_is_one(bo.exp + 0))
    {
        nmod_mpoly_scalar_mul_ui(u.data, bo.base + 0, bo.sign, ctx);
        a.set_length(1, ctx);
        a.sign = 1;
        nmod_mpoly_add_ui(a.base + 0, u.data, co, ctx);
        fmpz_one(a.exp + 0);
        a.length = !nmod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

	nmod_polyfactor b(bo, ctx), c;
	c.sign = co;
	c.length = 0;
	b.normalize(ctx);
	assert(b.is_nice(ctx));
	a.one();
	if (b.sign == 0)
	{
		a.sign = c.sign;
		a.append(c, ctx);
		return true;
	}
	if (c.sign == 0)
	{
		a.sign = b.sign;
		a.append(b, ctx);
		return true;
	}
	gcd_helper(a, b, c, ctx);
	b.get_mpoly(v.data, ctx);
	c.get_mpoly(w.data, ctx);
	nmod_mpoly_add(u.data, v.data, w.data, ctx);
	a.mul_mpoly(u.data, ctx);
    a.normalize(ctx);
	assert(a.is_nice(ctx));
	return true;
}
bool fmpz_mod_polyfactor_add_fmpz_mod(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & bo, const fmpz_t co, const fmpz_mod_mpoly_ctx_t ctx)
{
	xfmpz_mod_mpoly_t u(ctx), v(ctx), w(ctx);

    if (bo.length == 1 && fmpz_is_one(bo.exp + 0))
    {
        fmpz_mod_mpoly_scalar_mul_fmpz(u.data, bo.base + 0, bo.sign, ctx);
        a.set_length(1, ctx);
        fmpz_one(a.sign);
        fmpz_mod_mpoly_add_fmpz(a.base + 0, u.data, co, ctx);
        fmpz_one(a.exp + 0);
        a.length = !fmpz_mod_mpoly_is_zero(a.base + 0, ctx);
        return true;
    }

	fmpz_mod_polyfactor b(bo, ctx), c;
	fmpz_set(c.sign, co);
	c.length = 0;
	b.normalize(ctx);
	assert(b.is_nice(ctx));
	a.one();
	if (fmpz_is_zero(b.sign))
	{
		fmpz_set(a.sign, c.sign);
		a.append(c, ctx);
		return true;
	}
	if (fmpz_is_zero(c.sign))
	{
		fmpz_set(a.sign, b.sign);
		a.append(b, ctx);
		return true;
	}
	gcd_helper(a, b, c, ctx);
	b.get_mpoly(v.data, ctx);
	c.get_mpoly(w.data, ctx);
	fmpz_mod_mpoly_add(u.data, v.data, w.data, ctx);
	a.mul_mpoly(u.data, ctx);
    a.normalize(ctx);
	assert(a.is_nice(ctx));
	return true;
}


bool sub(djrat & a, djrat & bo, djrat & co)
{
	xfmpz_t u, v, w;
	djrat b(bo), c(co);
	assert(bo.is_nice());
	assert(co.is_nice());
	a.one();
	if (b.sign == 0)
	{
		a.sign = -c.sign;
		a.append(c);
		return true;
	}
	if (c.sign == 0)
	{
		a.sign = b.sign;
		a.append(b);
		return true;
	}
	gcd_helper(a, b, c);
	b.get_fmpz(v.data);
	c.get_fmpz(w.data);
	fmpz_sub(u.data, v.data, w.data);
	a.mul_fmpz(u.data);
    a.normalize();
	assert(a.is_nice());
	return true;
}

bool fmpq_polyfactor_sub(fmpq_polyfactor & a, const fmpq_polyfactor & bo, const fmpq_polyfactor & co, const fmpz_mpoly_ctx_t ctx)
{
	xfmpz_mpoly_t u(ctx), v(ctx), w(ctx);
	fmpq_polyfactor b(bo, ctx), c(co, ctx);
	assert(bo.is_nice(ctx));
	assert(co.is_nice(ctx));
	a.one();
	if (fmpq_is_zero(b.sign))
	{
		fmpq_set(a.sign, c.sign);
		a.append(c, ctx);
		return true;
	}
	if (fmpq_is_zero(c.sign))
	{
		fmpq_set(a.sign, b.sign);
		a.append(b, ctx);
		return true;
	}
	gcd_helper(a, b, c, ctx);
	b.get_mpoly(v.data, ctx);
	c.get_mpoly(w.data, ctx);
	fmpz_mpoly_sub(u.data, v.data, w.data, ctx);
	a.mul_mpoly(u.data, ctx);
    a.normalize(ctx);
	assert(a.is_nice(ctx));
	return true;
}

bool fmpq_polyfactor_pow(fmpq_polyfactor & a, const fmpz_t power, const fmpz_mpoly_ctx_t ctx)
{
	if (fmpz_is_zero(power))
	{
		a.one();
        return true;
	}

	BOOL_CHECK(fmpq_pow_fmpz(a.sign, a.sign, power));
	for (slong i = 0; i < a.length; i++)
    {
        fmpz_mul(a.exp + i, a.exp + i, power);
        if (fmpz_mpoly_length(a.base + i, ctx) == 1 && fmpz_sgn(a.exp + i) >= 0)
        {
            fmpz_mpoly_pow_fmpz(a.base + i, a.base + i, a.exp + i, ctx);
            fmpz_one(a.exp + i);
        }
    }

    return true;
}
bool nmod_polyfactor_pow(nmod_polyfactor & a, const fmpz_t power, const nmod_mpoly_ctx_t ctx)
{
	if (fmpz_is_zero(power))
	{
		a.one();
        return true;
	}

	BOOL_CHECK(nmod_pow_fmpz_checked(&a.sign, a.sign, power, ctx->mod));
	for (slong i = 0; i < a.length; i++)
		fmpz_mul(a.exp + i, a.exp + i, power);

    return true;
}
bool fmpz_mod_polyfactor_pow(fmpz_mod_polyfactor & a, const fmpz_t power, const fmpz_mod_mpoly_ctx_t ctx)
{
	if (fmpz_is_zero(power))
	{
		a.one();
        return true;
	}

	BOOL_CHECK(fmpz_mod_pow_fmpz(a.sign, a.sign, power, ctx->ffinfo));
	for (slong i = 0; i < a.length; i++)
		fmpz_mul(a.exp + i, a.exp + i, power);

    return true;
}

bool fmpq_polyfactor_factorize(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpz_mpoly_ctx_t ctx)
{
    xfmpz_mpoly_factor_t fac(ctx);
    xfmpq_t t;
    xfmpz_t s;

	fmpq_set(a.sign, b.sign);
	a.length = 0;
    for (slong i = 0; i < b.length; i++)
    {
        int r = fmpz_mpoly_factor(fac.data, b.base + i, ctx);
        if (r == 0)
            return false;
/*
        r = fmpq_pow_fmpz(t.data, fac.data->constant, b.exp + i);
        if (r == 0)
            return false;

        fmpq_mul(a.sign, a.sign, t.data);
*/
        for (slong j = 0; j < fac.data->num; j++)
        {
            fmpz_mul(s.data, b.exp + i, fac.data->exp + j);
            a.mul_mpoly_pow_fmpz(fac.data->poly + j, s.data, ctx);
        }
    }

    a.normalize(ctx);
    return true;
}
bool nmod_polyfactor_factorize(nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_mpoly_ctx_t ctx)
{
    xnmod_mpoly_factor_t fac(ctx);
    mp_limb_t t;
    xfmpz_t s;

	a.sign = b.sign;
	a.length = 0;
    for (slong i = 0; i < b.length; i++)
    {
        int r = nmod_mpoly_factor(fac.data, b.base + i, ctx);
        if (r == 0)
            return false;

        BOOL_CHECK(nmod_pow_fmpz_checked(&t, fac.data->constant, b.exp + i, ctx->mod));
        a.sign = nmod_mul(a.sign, t, ctx->mod);
        for (slong j = 0; j < fac.data->num; j++)
        {
            fmpz_mul(s.data, b.exp + i, fac.data->exp + j);
            a.mul_mpoly_pow_fmpz(fac.data->poly + j, s.data, ctx);
        }
    }

    a.normalize(ctx);
    return true;
}
bool fmpz_mod_polyfactor_factorize(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_factor_t fac;
    fmpz_mod_mpoly_factor_init(fac, ctx);
    xfmpz_t t;
    xfmpz_t s;

	fmpz_set(a.sign, b.sign);
	a.length = 0;
    for (slong i = 0; i < b.length; i++)
    {
        int r = fmpz_mod_mpoly_factor(fac, b.base + i, ctx);
        if (r == 0)
            return false;

        fmpz_mod_pow_fmpz(t.data, fac->constant, b.exp + i, ctx->ffinfo);
        fmpz_mod_mul(a.sign, a.sign, t.data, ctx->ffinfo);
        for (slong j = 0; j < fac->num; j++)
        {
            fmpz_mul(s.data, b.exp + i, fac->exp + j);
            a.mul_mpoly_pow_fmpz(fac->poly + j, s.data, ctx);
        }
    }

    fmpz_mod_mpoly_factor_clear(fac, ctx);

    a.normalize(ctx);
    return true;
}

bool fmpq_polyfactor_expand_numerator(fmpq_polyfactor & a, const fmpq_polyfactor & b, const fmpz_mpoly_ctx_t ctx)
{
    fmpq_set(a.sign, b.sign);
    a.length = 0;

    if (fmpq_is_zero(b.sign))
        return true;

    xfmpz_mpoly_t num(ctx), t1(ctx), t2(ctx);
    fmpz_mpoly_set_fmpz(num.data, fmpq_numref(b.sign), ctx);
    fmpz_one(fmpq_numref(a.sign));
    for (slong i = 0; i < b.length; i++)
    {
        if (fmpz_sgn(b.exp + i) >= 0)
        {
            if (fmpz_is_one(b.exp + i))
            {
                fmpz_mpoly_mul(t2.data, num.data, b.base + i, ctx);
            }
            else
            {
                fmpz_mpoly_pow_fmpz(t1.data, b.base + i, b.exp + i, ctx);
                fmpz_mpoly_mul(t2.data, num.data, t1.data, ctx);
            }
            fmpz_mpoly_swap(t2.data, num.data, ctx);
        }
        else
        {
            a.push_back(b.base + i, b.exp + i, ctx);
        }
    }
    a.push_back(num.data, 1, ctx);
	return true;
}
bool nmod_polyfactor_expand_numerator(nmod_polyfactor & a, const nmod_polyfactor & b, const nmod_mpoly_ctx_t ctx)
{
    a.sign = b.sign;
    a.length = 0;

    if (b.sign == 0)
        return true;

    xnmod_mpoly_t num(ctx), t1(ctx), t2(ctx);
    nmod_mpoly_set_ui(num.data, b.sign, ctx);
    a.sign = 1;
    for (slong i = 0; i < b.length; i++)
    {
        if (fmpz_sgn(b.exp + i) >= 0)
        {
            if (fmpz_is_one(b.exp + i))
            {
                nmod_mpoly_mul(t2.data, num.data, b.base + i, ctx);
            }
            else
            {
                nmod_mpoly_pow_fmpz(t1.data, b.base + i, b.exp + i, ctx);
                nmod_mpoly_mul(t2.data, num.data, t1.data, ctx);
            }
            nmod_mpoly_swap(t2.data, num.data, ctx);
        }
        else
        {
            a.push_back(b.base + i, b.exp + i, ctx);
        }
    }
    a.push_back(num.data, 1, ctx);
	return true;
}
bool fmpz_mod_polyfactor_expand_numerator(fmpz_mod_polyfactor & a, const fmpz_mod_polyfactor & b, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_set(a.sign, b.sign);
    a.length = 0;

    if (fmpz_is_zero(b.sign))
        return true;

    xfmpz_mod_mpoly_t num(ctx), t1(ctx), t2(ctx);
    fmpz_mod_mpoly_set_fmpz(num.data, b.sign, ctx);
    fmpz_one(a.sign);
    for (slong i = 0; i < b.length; i++)
    {
        if (fmpz_sgn(b.exp + i) >= 0)
        {
            if (fmpz_is_one(b.exp + i))
            {
                fmpz_mod_mpoly_mul(t2.data, num.data, b.base + i, ctx);
            }
            else
            {
                fmpz_mod_mpoly_pow_fmpz(t1.data, b.base + i, b.exp + i, ctx);
                fmpz_mod_mpoly_mul(t2.data, num.data, t1.data, ctx);
            }
            fmpz_mod_mpoly_swap(t2.data, num.data, ctx);
        }
        else
        {
            a.push_back(b.base + i, b.exp + i, ctx);
        }
    }
    a.push_back(num.data, 1, ctx);
	return true;
}


bool fmpq_polyfactor_map(
	fmpq_polyfactor & a,
	const fmpz_mpoly_ctx_t actx,
	const fmpq_polyfactor & b,
	const fmpz_mpoly_ctx_t bctx,
	std::vector<xfmpz_mpoly_t> & values)
{
	fmpz_mpoly_struct ** v = (fmpz_mpoly_struct **) malloc(values.size()*sizeof(fmpz_mpoly_struct *));
	assert(values.size() == bctx->minfo->nvars);
	for (size_t i = 0; i < values.size(); i++)
		v[i] = values[i].data;
	fmpq_set(a.sign, b.sign);
	a.set_length(b.length, actx);
	for (slong i = 0; i < b.length; i++)
	{
		fmpz_mpoly_compose_fmpz_mpoly(a.base + i, b.base + i, v, bctx, actx);
		fmpz_set(a.exp + i, b.exp + i);
	}
	free(v);
	return true;
}

bool nmod_polyfactor_map(
	nmod_polyfactor & a,
	const nmod_mpoly_ctx_t actx,
	const nmod_polyfactor & b,
	const nmod_mpoly_ctx_t bctx,
	std::vector<xnmod_mpoly_t> & values)
{
	nmod_mpoly_struct ** v = (nmod_mpoly_struct **) malloc(values.size()*sizeof(nmod_mpoly_struct *));
	assert(values.size() == bctx->minfo->nvars);
	for (size_t i = 0; i < values.size(); i++)
		v[i] = values[i].data;
	a.sign = b.sign;
	a.set_length(b.length, actx);
	for (slong i = 0; i < b.length; i++)
	{
		nmod_mpoly_compose_nmod_mpoly(a.base + i, b.base + i, v, bctx, actx);
		fmpz_set(a.exp + i, b.exp + i);
	}
	free(v);
	return true;
}

bool fmpz_mod_polyfactor_map(
	fmpz_mod_polyfactor & a,
	const fmpz_mod_mpoly_ctx_t actx,
	const fmpz_mod_polyfactor & b,
	const fmpz_mod_mpoly_ctx_t bctx,
	std::vector<xfmpz_mod_mpoly_t> & values)
{
	fmpz_mod_mpoly_struct ** v = (fmpz_mod_mpoly_struct **) malloc(values.size()*sizeof(fmpz_mod_mpoly_struct *));
	assert(values.size() == bctx->minfo->nvars);
	for (size_t i = 0; i < values.size(); i++)
		v[i] = values[i].data;
	fmpz_set(a.sign, b.sign);
	a.set_length(b.length, actx);
	for (slong i = 0; i < b.length; i++)
	{
		fmpz_mod_mpoly_compose_fmpz_mod_mpoly(a.base + i, b.base + i, v, bctx, actx);
		fmpz_set(a.exp + i, b.exp + i);
	}
	free(v);
	return true;
}



void fmpq_polyfactor_partial_fractions(
    std::vector<fmpq_polyfactor> & v,
    const fmpq_polyfactor & a,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
//std::cout << "fmpq_polyfactor_partial_fractions called" << std::endl;

//    rdense_poly<rfmpz_mpoly> R(ctx);
//    std::vector<dense_poly<xfmpz_mpoly_t>> r;
//    r.push_back(dense_poly<xfmpz_mpoly_t>());

    v.clear();
    for (slong i = 0; i < a.length; i++)
    {
        v.push_back(fmpq_polyfactor());
        fmpq_one(v.back().sign);
        v.back().push_back(a.base + i, -1, ctx);
    }

//std::cout << "fmpq_polyfactor_partial_fractions returning" << std::endl;

}
