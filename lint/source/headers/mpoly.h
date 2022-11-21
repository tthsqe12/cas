#pragma once

#include "fmpz.h"

struct mpoly_ctx {
    static constexpr ulimb min_bits = 4;

    ulimb _nvars;
    ulimb* _lut_var_offsets_sp;
    uint8_t* _lut_var_shifts_sp;
    chunk<fmpz> _leaf_fmpzs;
    chunk<ulimb> _leaf_uis;
    ulimb _lut_stride_sp[FLINT_BITS];
    ulimb _lut_overflow_mask_sp[FLINT_BITS];
    unsigned char _lut_fix_bits_sp[FLINT_BITS];

    mpoly_ctx(ulimb nvars);
    ~mpoly_ctx();
    void set_nvars(ulimb n);
    ulimb nvars() const {FLINT_ASSERT(_nvars > 0); return _nvars;}

    ulimb stride_sp(ulimb bits) const {
        FLINT_ASSERT(0 < bits && bits <= FLINT_BITS);
        return _lut_stride_sp[bits - 1];
    }

    ulimb stride_mp(ulimb bits) const {
        FLINT_ASSERT(0 < bits && (bits % FLINT_BITS) == 0);
        return bits/FLINT_BITS*_nvars;
    }

    ulimb stride(ulimb bits) const {
        return bits <= FLINT_BITS ? stride_sp(bits) : stride_mp(bits);
    }

    ulimb fix_bits_sp(ulimb bits) const {
        FLINT_ASSERT(0 < bits && bits <= FLINT_BITS);
        return _lut_fix_bits_sp[bits-1];
    }

    ulimb fix_bits(ulimb bits) const {
        FLINT_ASSERT(0 < bits);
        return bits <= FLINT_BITS ? fix_bits_sp(bits) : round_up(bits, FLINT_BITS);
    }

    const ulimb* var_offsets_sp(ulimb bits) const {
        FLINT_ASSERT(0 < bits && bits <= FLINT_BITS);
        return _lut_var_offsets_sp + _nvars*(bits);  // not lut shifted
    }

    const uint8_t* var_shifts_sp(ulimb bits) const {
        FLINT_ASSERT(0 < bits && bits <= FLINT_BITS);
        return _lut_var_shifts_sp + _nvars*(bits);   // note lut shifted
    }

    ulimb var_mask_sp(ulimb bits) const {
        FLINT_ASSERT(0 < bits && bits <= FLINT_BITS);
        return (-ulimb(1)) >> (FLINT_BITS - bits);        
    }

    ulimb overflow_mask_sp(ulimb bits) const {
        FLINT_ASSERT(0 < bits && bits <= FLINT_BITS);
        return _lut_overflow_mask_sp[bits-1];
    }

    template <bool MP>
    ulimb overflow_mask(ulimb bits) const {
        if constexpr (MP) {
            FLINT_ASSERT(bits%FLINT_BITS == 0);
            return bits/FLINT_BITS;
        } else {
            return _lut_overflow_mask_sp[bits-1];
        }
    }

    bool write_monomial(std::ostream& o, bool first, const ulimb* a, ulimb M, const char* x) const;

    fmpz* leaf_alloc_fmpz(ulimb n) {return _leaf_fmpzs.fit_alloc(n);}
    ulimb* leaf_alloc_ui(ulimb n) {return _leaf_uis.fit_alloc(n);}
};

template<ulimb N> inline void my_mpn_add_n(ulimb* x, const ulimb* a, const ulimb* b);

template<>
inline void my_mpn_add_n<1>(ulimb* x, const ulimb* a, const ulimb* b) {
    x[0] = a[0] + b[0];
}

template<>
inline void my_mpn_add_n<2>(ulimb* x, const ulimb* a, const ulimb* b) {
    ADD_SSAAAA(x[1], x[0], a[1], a[0], b[1], b[0]);
}

template<>
inline void my_mpn_add_n<3>(ulimb* x, const ulimb* a, const ulimb* b) {
    ADD_S3A3A3(x[2], x[1], x[0], a[2], a[1], a[0], b[2], b[1], b[0]);
}

template<ulimb N> inline void my_mpn_sub_n(ulimb* x, const ulimb* a, const ulimb* b);

template<>
inline void my_mpn_sub_n<1>(ulimb* x, const ulimb* a, const ulimb* b) {
    x[0] = a[0] - b[0];
}

template<>
inline void my_mpn_sub_n<2>(ulimb* x, const ulimb* a, const ulimb* b) {
    SUB_DDMMSS(x[1], x[0], a[1], a[0], b[1], b[0]);
}

template<>
inline void my_mpn_sub_n<3>(ulimb* x, const ulimb* a, const ulimb* b) {
    SUB_D3M3S3(x[2], x[1], x[0], a[2], a[1], a[0], b[2], b[1], b[0]);
}


/*
monomial operation parameters:
    MP: multiprecision or not. true by default for safety.
    ML: the default 0 uses general code for the arbitrary length M
*/
template <bool MP = true, ulimb ML = 0>
inline void mpoly_monomial_add(ulimb* x, const ulimb* a, const ulimb* b, ulimb M) {
    FLINT_ASSERT(ML == 0 || ML == M);
    if constexpr (MP) {
        if constexpr (ML == 0)
            mpn_add_n(x, a, b, M);
        else
            my_mpn_add_n<ML>(x, a, b);
    } else {
        if constexpr (ML == 0) {
            ulimb i = 0; do {
                x[i] = a[i] + b[i];
            } while (++i < M);
        } else {
           for (ulimb i = 0; i < ML; i++)
              x[i] = a[i] + b[i];
        }
    }
}

template <ulimb ML = 0>
inline void mpoly_monomial_set(ulimb* x, const ulimb* a, ulimb M) {
    FLINT_ASSERT(ML == 0 || ML == M);
    if constexpr (ML == 0) {
        ulimb i = 0; do {
            x[i] = a[i];
        } while (++i < M);
    } else {
        for (ulimb i = 0; i < ML; i++)
            x[i] = a[i];
    }
}

template <ulimb ML = 0>
inline void mpoly_monomial_swap(ulimb* x, ulimb* y, ulimb M) {
    FLINT_ASSERT(ML == 0 || ML == M);
    if constexpr (ML == 0) {
        ulimb i = 0; do {
            std::swap(x[i], y[i]);
        } while (++i < M);
    } else {
        for (ulimb i = 0; i < ML; i++)
            std::swap(x[i], y[i]);
    }
}

template <ulimb ML = 0>
inline void mpoly_monomial_zero(ulimb* x, ulimb M) {
    FLINT_ASSERT(ML == 0 || ML == M);
    if constexpr (ML == 0) {
        ulimb i = 0; do {
            x[i] = 0;
        } while (++i < M);
    } else {
       for (ulimb i = 0; i < ML; i++)
          x[i] = 0;
    }
}

// mask is really bits/FLINT_BITS in the MP = true case
template <bool MP, ulimb ML = 0>
bool mpoly_monomial_overflowed(const ulimb* a, ulimb M, ulimb mask) {
    if constexpr (MP) {
        ulong i = mask - 1; do {
            if (has_top_bit(a[i]))
                return true;
            i += mask;
        } while (i < M);
    } else {
        if constexpr (ML == 0) {
            ulimb i = 0; do {
                if ((mask & a[i]) != 0)
                    return true;
            } while (++i < M);
        } else {
            for (ulimb i = 0; i < ML; i++)
                if ((mask & a[i]) != 0)
                    return true;
        }
    }
    return false;
}


template <bool MP, ulimb ML = 0>
inline bool mpoly_monomial_set_overflowed(ulimb* x, const ulimb* a, ulimb M, ulimb mask) {
    FLINT_ASSERT(ML == 0 || ML == M);
    if constexpr (MP) {
        mpoly_monomial_set<ML>(x, a, M);
        return mpoly_monomial_overflowed<MP, ML>(x, M, mask);
    } else {
        ulimb of = 0;
        if constexpr (ML == 0) {
            ulimb i = 0; do {
                of |= mask & a[i];
                x[i] = a[i];
            } while (++i < M);
        } else {
            for (ulimb i = 0; i < ML; i++) {
                of |= mask & a[i];
                x[i] = a[i];
            }
        }
        return of != 0;
    }
}

template <bool MP, ulimb ML = 0>
inline bool mpoly_monomial_sub_overflowed(ulimb* x, const ulimb* a, const ulimb* b, ulimb M, ulimb mask) {
    FLINT_ASSERT(ML == 0 || ML == M);
    if constexpr (MP) {
        if constexpr (ML == 0)
            mpn_sub_n(x, a, b, M);
        else
            my_mpn_sub_n<ML>(x, a, b);
        return mpoly_monomial_overflowed<MP, ML>(x, M, mask);
    } else {
        ulimb of = 0;
        if constexpr (ML == 0) {
            ulimb i = 0; do {
                of |= mask & (a[i] - b[i]);
                x[i] = (a[i] - b[i]);
            } while (++i < M);
        } else {
            for (ulimb i = 0; i < ML; i++) {
                of |= mask & (a[i] - b[i]);
                x[i] = (a[i] - b[i]);
            }
        }
        return of != 0;
    }
}


inline void mpoly_copy_monomials(ulimb * exp1, const ulimb * exp2, slimb len, ulimb M) {
    memcpy(exp1, exp2, M*len*sizeof(ulimb));
}

inline void mpoly_zero_monomials(ulimb * exp1, slimb len, ulimb M) {
    FLINT_ASSERT(len*M > 0);
    memset(exp1, 0, M*len*sizeof(ulimb));
}

template <ulimb ML = 0>
int mpoly_monomial_cmp(const ulimb * a, const ulimb * b, ulimb N)
{
    slimb i = N - 1;
    do {
        if (a[i] != b[i])
            return a[i] > b[i] ? 1 : -1;
    } while (--i >= 0);
    return 0;
}

template <ulimb ML = 0>
bool mpoly_monomial_gt(const ulimb * a, const ulimb * b, ulimb N)
{
    slimb i = N - 1;
    do {
        if (a[i] != b[i])
            return a[i] > b[i];
    } while (--i >= 0);
    return false;
}

template <ulimb ML = 0>
bool mpoly_monomial_equal(const ulimb * a, const ulimb * b, ulimb N)
{
    slimb i = N - 1;
    do {
        if (a[i] != b[i])
            return false;
    } while (--i >= 0);
    return true;
}

struct mpoly : my_vec<ulimb> {
    ulimb _bits = mpoly_ctx::min_bits;

    ulimb bits() const {return _bits;}
    void set_bits(ulimb n) {_bits = n;}

    void zero(mpoly_ctx& ctxm) {
        set_length(0);
        set_bits(ctxm.fix_bits(mpoly_ctx::min_bits));
    }

    bool is_canonical(const mpoly_ctx&) const;
    bool repacked_bits(mpoly_ctx& ctxm, ulimb ebits);
    void push_monomial(mpoly_ctx& ctxm, const ulimb* e);
    void push_monomial(mpoly_ctx& ctxm, const fmpz* e);
};

struct mpolyc {
    ulimb len;
    ulimb bits;
    ulimb* exps;

    mpolyc(mpoly a) : len(a.length()), bits(a.bits()), exps(a.data()) {}
};


void mpoly_degrees(mpoly_ctx& ctxm, fmpz* Adegs, mpolyc A);
bool mpoly_repacked_down_monomials(mpoly_ctx& ctxm, ulong* Aexps, ulimb Abits, mpolyc B);
bool mpoly_repacked_up_monomials(mpoly_ctx& ctxm, ulong* Aexps, ulimb Abits, mpolyc B);
bool mpoly_repacked_monomials(mpoly_ctx& ctxm, ulong* Aexps, ulimb Abits, mpolyc B);

bool mpoly_proved_nonsquare(mpoly_ctx& ctx, mpolyc A);
bool mpoly_proved_irreducible(mpoly_ctx& ctx, mpolyc A, random_state& state);

void fmpz_vec_add(fmpz* x, const fmpz* a, const fmpz* b, ulimb len);
slimb fmpz_vec_max_bits(const fmpz* a, ulimb len);
ulimb ui_vec_max_bits(const ulimb* a, ulimb an);

