#pragma once

#include "packed.h"

template <int vsize>
struct pd_fft_ctx {
    packed<double,vsize> n;
    packed<double,vsize> ninv;
    nmod_t mod[vsize];
    ulong primitive_root[vsize];
    my_vector<packed<double,vsize>> wrevtab;

    void set_mod(const ulong* p);
    void fit_wrevtab(ulong k);
    packed<double,vsize> w2rev(ulong j) {return wrevtab.at(j);}
    packed<double,vsize> w2revinv(ulong j) {
        return (j == 0) ? wrevtab.at(j) : neg(wrevtab.at(j^(saturate_bits(j)>>1)));
    }

    void fft_full_iter(packed<double,vsize>* x, ulong k, ulong j);
    void fft_trunc(packed<double,vsize>* x, ulong k, ulong j, ulong itrunc, ulong otrunc);
    void fft_trunc(packed<double,vsize>* x, ulong depth, ulong itrunc, ulong otrunc) {
        fit_wrevtab(depth);
        fft_trunc(x, depth, 0, itrunc, otrunc);
    }

    template<bool j_is_0> void ifft_full_iter(packed<double,vsize>* x, ulong k, ulong j);
    void ifft_full(packed<double,vsize>* x, ulong k, ulong j);
    void ifft_trunc(packed<double,vsize>* x, ulong k, ulong j, ulong z, ulong nn, bool f);
    void ifft_trunc(packed<double,vsize>* x, ulong depth, ulong trunc) {
        fit_wrevtab(depth);
        ifft_trunc(x, depth, 0, trunc, trunc, false);
    }
};

#define USE_BITS 50

struct mpn_fft_ctx {
    my_vector<packed<double,4>> _two_powers;
    my_vector<packed<double,4>> buffer;
    pd_fft_ctx<4> ctx;

    mpn_fft_ctx() {
#if USE_BITS == 51
        ulong primes[4] = {0x0007ff9000000001,
                           0x0007fe7800000001,
                           0x0007fd8800000001,
                           0x0007fd8000000001};
#elif USE_BITS == 50
        ulong primes[4] = {0x0003fff340000001,
                           0x0003fff300000001,
                           0x0003ffeec0000001,
                           0x0003ffed00000001};
#else
        ulong primes[4] = {0x0000ffdf00000001,
                           0x0000ffcd00000001,
                           0x00007fc800000001,
                           0x00007fbf00000001};
#endif
        ctx.set_mod(primes);
    }

    const packed<double,4>* two_power_table(ulong len)
    {
        if (_two_powers.size() >= len)
            return _two_powers.data();

        _two_powers.resize(len);

        packed<double,4>* d = _two_powers.data();
        packed<double,4> p = ctx.n;
        packed<double,4> two_over_p = add(ctx.ninv, ctx.ninv);
        packed<double,4> t = packed<double,4>(1.0);
        ulong j = 0;
        do {
            d[j] = t;
            t = fnmadd(round(mul(t, two_over_p)), p, add(t,t));
        } while (++j < len);

        return d;
    }

    void my_mpn_mul(ulong* c, const ulong* a, ulong an, const ulong* b, ulong bn);
};

