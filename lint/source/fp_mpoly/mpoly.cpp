#include "fp_mpoly.h"

std::string monomial_format(
    const ulimb* a,
    ulimb M)
{
    std::stringstream o;
    for (ulimb i = M; i > 0; i--) {
        if (i < M)
            o << "_";
        o << format_hex(a[i-1]);
    }
    return o.str();
}

bool mpoly::is_canonical(const mpoly_ctx& ctxm) const
{
    ulimb M = ctxm.stride(bits());
    const ulimb* exps = data();
    const ulimb l = length();

    if (bits() > FLINT_BITS)
    {
        ulimb mask = ctxm.overflow_mask<true>(bits());
        for (ulimb i = 0; i < l; i++) {
            if (mpoly_monomial_overflowed<true,0>(exps + M*i, M, mask)) {
                std::cerr << "mpoly mp exponents overflowed" << std::endl;
                return false;
            }
        }
    }
    else
    {
        ulimb mask = ctxm.overflow_mask<false>(bits());
        for (ulimb i = 0; i < l; i++) {
            if (mpoly_monomial_overflowed<false,0>(exps + M*i, M, mask)) {
                std::cerr << "mpoly sp exponents overflowed" << std::endl;
                return false;
            }
        }

    }

    for (ulimb i = 1; i < l; i++) {
        if (mpoly_monomial_cmp<0>(exps + M*(i-1), exps + M*i, M) <= 0) {
            std::cerr << "mpoly exponents out of order " << std::endl;            
            std::cerr << "exp[" << i-1 <<"]: " << monomial_format(exps + M*(i-1), M) << std::endl;
            std::cerr << "exp[" << i-0 <<"]: " << monomial_format(exps + M*(i-0), M) << std::endl;
            return false;
        }
    }

    return true;
}

bool fp_mpoly_is_canonical(
    const fp_mpoly_ring& ctx,
    const fp_mpoly& a)
{
    if (!a.m.is_canonical(ctx.m))
        return false;

    for (ulimb i = 0; i < a.length(); i++) {
        if (!fp_is_canonical(ctx.c, a.coeffs() + ctx.c.stride()*i))
            return false;
    }

    for (ulimb i = 0; i < a.length(); i++) {
        if (fp_is_zero(ctx.c, a.coeffs() + ctx.c.stride()*i)) {
            std::cerr << "fp_mpoly has a zero coefficient" << std::endl;
            return false;
        }
    }

    return true;
}

static void _mpoly_ctx_clear(mpoly_ctx& ctxm)
{
    ctxm._lut_var_offsets_sp += ctxm._nvars;
    ctxm._lut_var_shifts_sp += ctxm._nvars;
    delete[] ctxm._lut_var_offsets_sp;
    delete[] ctxm._lut_var_shifts_sp;
    ctxm._lut_var_offsets_sp = nullptr;
    ctxm._lut_var_shifts_sp = nullptr;
    ctxm._nvars = 0;
}

static void _mpoly_ctx_init(mpoly_ctx& ctxm, ulimb n)
{
    FLINT_ASSERT(n > 0);
    for (ulimb i = 1; i <= FLINT_BITS; i++)
        ctxm._lut_stride_sp[i - 1] = (n - 1)/(FLINT_BITS/i) + 1;
    for (ulimb i = 1; i <= FLINT_BITS; i++)
    {
        ulimb j = std::max(i, mpoly_ctx::min_bits);
        while (j < FLINT_BITS && ctxm._lut_stride_sp[j - 1] == ctxm._lut_stride_sp[j])
            j++;
        ctxm._lut_fix_bits_sp[i - 1] = j;
    }

    ctxm._lut_var_offsets_sp = new ulimb[n*FLINT_BITS];
    ctxm._lut_var_shifts_sp = new uint8_t[n*FLINT_BITS];
    ctxm._lut_var_offsets_sp -= n;
    ctxm._lut_var_shifts_sp -= n;
    ctxm._nvars = n;
    for (ulimb i = 1; i <= FLINT_BITS; i++)
    {
        ulimb fpw = FLINT_BITS/i;
        for (ulimb v = 0; v < n; v++)
        {
            ctxm._lut_var_offsets_sp[v+n*i] = (n - 1 - v)/fpw;
            ctxm._lut_var_shifts_sp[v+n*i]  = (n - 1 - v)%fpw*i;
        }
    }

    for (ulimb i = 1; i <= FLINT_BITS; i++)
    {
        ulimb mask = pow2(i - 1);
        for (ulimb j = i; j < FLINT_BITS; j += i)
            mask = (mask << i) + pow2(i - 1);
        ctxm._lut_overflow_mask_sp[i-1] = mask;
    }
}

mpoly_ctx::mpoly_ctx(ulimb n) {
    _mpoly_ctx_init(*this, n);
}

mpoly_ctx::~mpoly_ctx() {
    _mpoly_ctx_clear(*this);
}

void mpoly_ctx::set_nvars(ulimb n) {
    if (_nvars == n)
        return;
    _mpoly_ctx_clear(*this);
    _mpoly_ctx_init(*this, n);
}

template <ulimb ML>
inline void mpoly_monomial_max_sp(
    ulimb* exp1, const ulimb * exp2, const ulimb * exp3,
    ulimb bits,
    ulimb M,
    ulimb mask)
{
    FLINT_ASSERT(ML == 0 || M == ML);
    ulimb i, s, m;
    if constexpr (ML == 0)
    {
        i = 0; do {
            s = mask + exp2[i] - exp3[i];
            m = mask & s;
            m = m - (m >> (bits - 1));
            exp1[i] = exp3[i] + (s & m);
        } while (++i < M);
    }
    else
    {
        for (i = 0; i < ML; i++)
        {
            s = mask + exp2[i] - exp3[i];
            m = mask & s;
            m = m - (m >> (bits - 1));
            exp1[i] = exp3[i] + (s & m);
        }
    }
}

template <ulimb ML>
void mpoly_monomials_max_sp(
    ulimb* p,
    mpolyc A,
    ulimb M,
    ulimb mask)
{
    FLINT_ASSERT(A.len > 0);
    if constexpr (ML == 0)
    {
        mpoly_monomial_set<ML>(p, A.exps + M*0, M);
        for (ulimb i = 1; i < A.len; i++)
            mpoly_monomial_max_sp<ML>(p, p, A.exps + M*i, A.bits, M, mask);
    }
    else
    {
        ulimb tp[ML];
        mpoly_monomial_set<ML>(tp, A.exps + M*0, M);
        for (ulimb i = 1; i < A.len; i++)
            mpoly_monomial_max_sp<ML>(tp, tp, A.exps + M*i, A.bits, M, mask);
        mpoly_monomial_set<ML>(p, tp, M);
    }
}

void mpoly_monomial_max_mp(ulimb * exp1, const ulimb * exp2, const ulimb * exp3,
                                                  ulimb bits, slimb N)
{
    slimb i, j;
    for (i = 0; i < N; i += bits/FLINT_BITS)
    {
        const ulimb * t = exp2;
        for (j = bits/FLINT_BITS - 1; j >= 0; j--)
        {
            if (exp3[i + j] != exp2[i + j])
            {
                if (exp3[i + j] > exp2[i + j])
                    t = exp3;
                break;
            }
        }
        for (j = 0; j < bits/FLINT_BITS; j++)
            exp1[i + j] = t[i + j];
    }
}


template <ulimb ML>
inline void mpoly_monomial_min_sp(
    ulimb* exp1, const ulimb * exp2, const ulimb * exp3,
    ulimb bits,
    ulimb M,
    ulimb mask)
{
    FLINT_ASSERT(ML == 0 || M == ML);
    ulimb i, s, m;
    if constexpr (ML == 0)
    {
        i = 0; do {
            s = mask + exp2[i] - exp3[i];
            m = mask & s;
            m = m - (m >> (bits - 1));
            exp1[i] = exp2[i] - (s & m);
        } while (++i < M);
    }
    else
    {
        for (i = 0; i < ML; i++)
        {
            s = mask + exp2[i] - exp3[i];
            m = mask & s;
            m = m - (m >> (bits - 1));
            exp1[i] = exp2[i] - (s & m);
        }
    }
}

void mpoly_degrees(
    mpoly_ctx& ctxm,
    fmpz* Adegs,
    mpolyc A)
{
    ulimb M = ctxm.stride(A.bits);
    ulimb* p = ctxm.leaf_alloc_ui(M);

    if (UNLIKELY(A.bits > FLINT_BITS))
    {
        mpoly_monomial_set(p, A.exps + M*0, M);
        for (ulimb i = 1; i < A.len; i++)
            mpoly_monomial_max_mp(p, p, A.exps + M*i, A.bits, M);

        ulimb m = A.bits/FLINT_BITS;
        ulimb n = ctxm.nvars();
        ulimb v = 0; do {
            Adegs[v].set_limbs(p + m*(n-1-v), m, 0);
        } while (++v < n);
    }
    else
    {
        ulimb mask = ctxm.overflow_mask_sp(A.bits);
        if (M == 1)
            mpoly_monomials_max_sp<1>(p, A, M, mask);
        else if (M == 2)
            mpoly_monomials_max_sp<2>(p, A, M, mask);
        else if (M == 3)
            mpoly_monomials_max_sp<3>(p, A, M, mask);
        else
            mpoly_monomials_max_sp<0>(p, A, M, mask);

        mask = ctxm.var_mask_sp(A.bits);
        auto shifts = ctxm.var_shifts_sp(A.bits);
        auto offsets = ctxm.var_offsets_sp(A.bits);
        ulimb v = 0; do {
            fmpz_set_ui(Adegs[v], (p[offsets[v]]>>shifts[v])&mask);
        } while (++v < ctxm.nvars());
    }
}


bool mpoly_repacked_down_monomials(
    mpoly_ctx& ctxm,
    ulong* Aexps, ulimb Abits,
    mpolyc B)
{
    FLINT_ASSERT_ALWAYS(false && "not implemented");
}

bool mpoly_repacked_up_monomials(
    mpoly_ctx& ctxm,
    ulong* Aexps, ulimb Abits,
    mpolyc B)
{
    FLINT_ASSERT(Abits > B.bits);

    ulimb n = ctxm.nvars();
    ulimb AM = ctxm.stride(Abits);
    ulimb BM = ctxm.stride(B.bits);

//std::cout << "repacking " << B.bits << "(" << BM << ")" << " -> " << Abits << "(" << AM << ")"  << " with " << ctxm.nvars() << " vars" << std::endl;

    mpoly_zero_monomials(Aexps, B.len, AM);

    if (LIKELY(B.bits < FLINT_BITS))
    {
        auto Boffsets = ctxm.var_offsets_sp(B.bits);
        auto Bshifts = ctxm.var_shifts_sp(B.bits);
        auto Bmask = ctxm.var_mask_sp(B.bits);

        if (Abits < FLINT_BITS)
        {
            auto Aoffsets = ctxm.var_offsets_sp(Abits);
            auto Ashifts = ctxm.var_shifts_sp(Abits);
            for (ulimb i = 0; i < B.len; i++)
            {
                ulimb v = 0; do {
                    ulimb e = ((B.exps + BM*i)[Boffsets[v]] >> Bshifts[v]) & Bmask;
                    (Aexps + AM*i)[Aoffsets[v]] |= e << Ashifts[v];
                } while (++v < n);
            }
        }
        else
        {
            for (ulimb i = 0; i < B.len; i++)
            {
                ulimb v = 0; do {
                    ulimb e = ((B.exps + BM*i)[Boffsets[v]] >> Bshifts[v]) & Bmask;
                    (Aexps + AM*i)[Abits/FLINT_BITS*(n-1-v)] = e;
                } while (++v < n);
            }
        }
    }
    else
    {
        FLINT_ASSERT(Abits >= FLINT_BITS && B.bits >= FLINT_BITS);
        for (ulimb i = 0; i < B.len; i++)
        {
            ulimb v = 0; do { // iterate backwards
                ui_vec_set(Aexps + AM*i + Abits/FLINT_BITS*v,
                           B.exps + BM*i + B.bits/FLINT_BITS*v, B.bits/FLINT_BITS);
            } while (++v < n);
        }
    }

//std::cout << "repacked: " << std::endl;
//for (ulimb k = 0; k < B.len; k++)
//    std::cout << "[" << k << "]: " << monomial_format(B.exps + BM*k, BM) << " -> " << monomial_format(Aexps + AM*k, AM) << std::endl;

    return true;
}


/* return true if the repacking was successful, false if it failed */
bool mpoly_repacked_monomials(
    mpoly_ctx& ctxm,
    ulong* Aexps, ulimb Abits,
    mpolyc B)
{
    if (Abits < B.bits)
        return mpoly_repacked_down_monomials(ctxm, Aexps, Abits, B);

    if (Abits > B.bits)
        return mpoly_repacked_up_monomials(ctxm, Aexps, Abits, B);

    mpoly_copy_monomials(Aexps, B.exps, B.len, ctxm.stride(B.bits));
    return true;
}


bool mpoly::repacked_bits(mpoly_ctx& ctxm, ulimb ebits)
{
    ulimb l = length();

    if (ebits == bits())
        return true;

    if (l < 1)
    {
        set_bits(ebits);
        return true;
    }

    tmp_allocator push;
    ulimb M = ctxm.stride(bits());
    ulimb eM = ctxm.stride(ebits);

    if (ebits > bits())
    {
        mpolyc B(*this);
        B.exps = push.recursive_alloc<ulimb>(l*M);
        mpoly_copy_monomials(B.exps, data(), l, M);
        mpoly_repacked_up_monomials(ctxm, fit_alloc(l+1, eM), ebits, B);
    }
    else
    {
        ulimb* eexps = push.recursive_alloc<ulimb>(l*eM);
        if (!mpoly_repacked_down_monomials(ctxm, eexps, ebits, *this))
            return false;
        mpoly_copy_monomials(data(), eexps, l, M);
    }

    set_bits(ebits);
    return true;
}

void mpoly::push_monomial(mpoly_ctx& ctxm, const ulimb* e) {
    ulimb n = ctxm.nvars();
    ulimb ebits = 1 + ui_vec_max_bits(e, n);
    ulimb l = length();

    if (UNLIKELY(ebits > bits()))
        repacked_bits(ctxm, ctxm.fix_bits(ebits));

    ulimb M = ctxm.stride(bits());
    ulimb* texps = fit_alloc(l+1, M);
    set_length(l + 1);

    mpoly_monomial_zero(texps + M*l, M);
    if (bits() < FLINT_BITS)
    {
        auto offsets = ctxm.var_offsets_sp(bits());
        auto shifts = ctxm.var_shifts_sp(bits());
        ulimb v = 0; do {
            (texps + M*l)[offsets[v]] |= e[v] << shifts[v];
        } while (++v < n);
    }
    else
    {
        ulimb v = 0; do {
            (texps + M*l)[bits()/FLINT_BITS*v] = e[v];
        } while (++v < n);
    }
}

void mpoly::push_monomial(mpoly_ctx& ctxm, const fmpz* e) {
    ulimb n = ctxm.nvars();
    ulimb ebits = 1 + fmpz_vec_max_bits(e, n);
    FLINT_ASSERT(!has_top_bit(ebits-1));
    ulimb l = length();

    if (UNLIKELY(ebits > bits()))
        repacked_bits(ctxm, ctxm.fix_bits(ebits));

    ulimb M = ctxm.stride(bits());
    ulimb* texps = fit_alloc(l+1, M);
    set_length(l + 1);

    mpoly_monomial_zero(texps + M*l, M);
    if (bits() < FLINT_BITS)
    {
        auto offsets = ctxm.var_offsets_sp(bits());
        auto shifts = ctxm.var_shifts_sp(bits());
        ulimb v = 0; do {
            (texps + M*l)[offsets[v]] |= fmpz_get_ui(e[v]) << shifts[v];
        } while (++v < n);
    }
    else
    {
        ulimb m = bits()/FLINT_BITS;
        ulimb v = 0; do {
            fmpz_get_ui_vec(texps + M*l + m*v, m, e[v]);
        } while (++v < n);
    }
}

