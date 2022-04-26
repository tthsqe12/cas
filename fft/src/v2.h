#pragma once

#include <cassert>
#include "packed.h"
#include "misc.h"

#define BLK_SZ 256
#define LG_BLK_SZ 8
#define BLK_SHIFT 10
#define VEC_SZ 4

struct fftv2_ctx {
    double* data;
    double p;
    double pinv;
    nmod_t mod;
    ulong primitive_root;
    ulong depth;
    ulong blk_sz;
    double* w2s;
    ulong wtab_depth;

    void init_prime(ulong p);

    ~fftv2_ctx() {free(w2s);}

    fftv2_ctx(const fftv2_ctx& other) = delete;

    fftv2_ctx(fftv2_ctx&& other) noexcept
    {
        data = other.data;
        p = other.p;
        pinv = other.pinv;
        mod = other.mod;
        primitive_root = other.primitive_root;
        depth = other.depth;
        blk_sz = other.blk_sz;
        w2s = other.w2s;
        wtab_depth = other.wtab_depth;
        other.w2s = nullptr;
        other.wtab_depth = 0;
    }

    fftv2_ctx& operator=(fftv2_ctx&& other) noexcept
    {
        if (this != &other) {
            data = other.data;
            p = other.p;
            pinv = other.pinv;
            mod = other.mod;
            primitive_root = other.primitive_root;
            depth = other.depth;
            blk_sz = other.blk_sz;
            w2s = other.w2s;
            wtab_depth = other.wtab_depth;
            other.w2s = nullptr;
            other.wtab_depth = 0;
        }  
        return *this;  
    }

    fftv2_ctx(ulong pp)
    {
        blk_sz = BLK_SZ;
        data = nullptr;
        w2s = nullptr;
        init_prime(pp);
    }

    void set_prime(ulong pp)
    {
        blk_sz = BLK_SZ;
        std::free(w2s);
        w2s = nullptr;
        init_prime(pp);
    }

    void fit_wtab(ulong k);

    inline ulong offset(ulong I)
    {
        return (I << LG_BLK_SZ) + 4*(I >> (BLK_SHIFT+2));
    }

    ulong data_size()
    {
        return offset(ulong(1) << (depth - LG_BLK_SZ));
    }

    void set_data(double* d)
    {
        data = d;
    }
    double* release_data()
    {
        double* d = data;
        data = nullptr;
        return d;
    }

    void set_depth(ulong l)
    {
        depth = l;
        if (l < LG_BLK_SZ)
        {
            std::cout << "depth " << l << " too small" << std::endl;
            abort();
        }
        fit_wtab(l);
    }

    inline double* from_index(ulong I)
    {
        return data + offset(I);
    }

    inline void set_index(ulong i, double x)
    {
        from_index(i/BLK_SZ)[i%BLK_SZ] = x;
    }

    inline double get_index(ulong i)
    {
        return from_index(i/BLK_SZ)[i%BLK_SZ];
    }

    void fft_base(ulong I, ulong j);
    template<int depth, typename std::enable_if<depth == 0>::type* = nullptr> inline void fft_basecase(double* X, ulong j);
    template<int depth, typename std::enable_if<depth == 1>::type* = nullptr> inline void fft_basecase(double* X, ulong j);
    template<int depth, typename std::enable_if<depth == 2>::type* = nullptr> inline void fft_basecase(double* X, ulong j);
    template<int depth, typename std::enable_if<depth == 4>::type* = nullptr> inline void fft_basecase(double* X, ulong j);
    template<int depth, typename std::enable_if<depth >= 6>::type* = nullptr> inline void fft_basecase(double* X, ulong j);

    template<bool j_can_be_0> void ifft_base(ulong I, ulong j);
    template<int depth, bool j_is_0, typename std::enable_if<depth == 0>::type* = nullptr> inline void ifft_basecase(double* X, ulong j, ulong jm);
    template<int depth, bool j_is_0, typename std::enable_if<depth == 1>::type* = nullptr> inline void ifft_basecase(double* X, ulong j, ulong jm);
    template<int depth, bool j_is_0, typename std::enable_if<depth == 2>::type* = nullptr> inline void ifft_basecase(double* X, ulong j, ulong jm);
    template<int depth, bool j_is_0, typename std::enable_if<depth == 4>::type* = nullptr> inline void ifft_basecase(double* X, ulong j, ulong jm);
    template<int depth, bool j_is_0, typename std::enable_if<depth >= 6>::type* = nullptr> inline void ifft_basecase(double* X, ulong j, ulong jm);

    void fft_main(ulong I, ulong S, ulong k, ulong j);
    void fft_main_block(ulong I, ulong S, ulong k, ulong j);
    void fft_trunc(ulong I, ulong S, ulong k, ulong j, ulong itrunc, ulong otrunc);
    void fft_trunc_block(ulong I, ulong S, ulong k, ulong j, ulong itrunc, ulong otrunc);
    inline void fft() {
        fft_main(0, 1, depth - LG_BLK_SZ, 0);
    }
    inline void fft_trunc(ulong itrunc, ulong otrunc) {
        assert(itrunc % BLK_SZ == 0 && otrunc % BLK_SZ == 0);
        fft_trunc(0, 1, depth - LG_BLK_SZ, 0, itrunc/BLK_SZ, otrunc/BLK_SZ);
    }

    void ifft_main(ulong I, ulong S, ulong k, ulong j);
    void ifft_main_block(ulong I, ulong S, ulong k, ulong j);
    void ifft_trunc(ulong I, ulong S, ulong k, ulong j, ulong z, ulong n, bool f);
    void ifft_trunc_block(ulong I, ulong S, ulong k, ulong j, ulong z, ulong n, bool f);
    inline void ifft() {
        ifft_main(0, 1, depth - LG_BLK_SZ, 0);
    }
    inline void ifft_trunc(ulong trunc) {
        assert(trunc % BLK_SZ == 0);
        ifft_trunc(0, 1, depth - LG_BLK_SZ, 0, trunc/BLK_SZ, trunc/BLK_SZ, false);
    }

    void from_mpn(const ulong * a, ulong an, ulong bits);
    void point_mul(const double* b, ulong mm);
};



struct crt_data {
    ulong prime;
    ulong coeff_len;
    ulong nprimes;
    std::vector<ulong> data;

    crt_data(ulong _prime, ulong _coeff_len, ulong _nprimes)
    {
        prime = _prime;
        coeff_len = _coeff_len;
        nprimes = _nprimes;
        data.resize(nprimes*coeff_len + coeff_len + nprimes);
    }

    // return mpn of length coeff_len
    ulong* co_prime(ulong i) {
        assert(i < nprimes);
        return data.data() + i*coeff_len;
    }

    // return mpn of length coeff_len
    ulong* prod_primes() {
        return data.data() + nprimes*coeff_len;
    }

    // the reduction of co_prime(i)
    ulong& co_prime_red(ulong i) {
        assert(i < nprimes);
        return data[nprimes*coeff_len + coeff_len + i];
    }

    // need  ceil(64*bn/bits) <= prod_primes/2^(2*bits)
    // i.e.  (64*bn+bits-1)/bits <= prod_primes/2^(2*bits)
    //       64*bn <= bits*prod_primes/2^(2*bits) - (bits-1)
    ulong find_bound(ulong bits)
    {
        ulong bound = -ulong(1);    
        fmpz_t x;
        fmpz_init(x);
        fmpz_set_ui_array(x, prod_primes(), coeff_len);
        fmpz_mul_ui(x, x, bits);
        fmpz_fdiv_q_2exp(x, x, 2*bits);
        fmpz_sub_ui(x, x, bits - 1);
        fmpz_fdiv_q_2exp(x, x, 6);
        if (fmpz_sgn(x) <= 0)
            bound = 0;
        else if (fmpz_abs_fits_ui(x))
            bound = fmpz_get_ui(x);
        fmpz_clear(x);
        return bound;
    }
};

typedef void (*to_ffts_func)(
        fftv2_ctx* Qffts,
        const ulong* a_, ulong an_, ulong atrunc,
        const packed<double,VEC_SZ>* two_pow);

typedef void (*from_ffts_func)(
        ulong* z, ulong zn, ulong zlen,
        fftv2_ctx* Qffts,
        crt_data* Qcrts,
        ulong bits);

struct profile_entry {
    ulong np;
    ulong bits;
    ulong bn_bound;
    to_ffts_func to_ffts;
    from_ffts_func from_ffts;

    profile_entry(ulong np_, ulong bits_, ulong bn_bound_,
                  to_ffts_func to_ffts_, from_ffts_func from_ffts_) :
            np(np_), bits(bits_), bn_bound(bn_bound_),
            to_ffts(to_ffts_), from_ffts(from_ffts_)
    {};
};


struct mpn_ctx_v2 {
    std::vector<fftv2_ctx> ffts;
    std::vector<crt_data> crts;
    std::vector<my_vector<packed<double,VEC_SZ>>> two_pow_tabs;
    std::vector<profile_entry> profiles;
    double* double_buffer = nullptr;
    ulong double_buffer_alloc = 0;

    mpn_ctx_v2(ulong p);

    ~mpn_ctx_v2() {
        free(double_buffer);
    }

    template <ulong np, ulong bits, ulong n, ulong m> void push_profile();
    const profile_entry* best_profile(ulong an, ulong bn);

    const packed<double,VEC_SZ>* two_pow_table(ulong len, ulong np);

    double* fit_double_buffer(ulong n)
    {
        if (n > double_buffer_alloc)
        {
            std::free(double_buffer);
            double_buffer = reinterpret_cast<double*>(my_aligned_alloc(4096, n*sizeof(double)));
            double_buffer_alloc = n;
        }
        return double_buffer;
    }

    ulong nprimes() {
        return ffts.size();
    }

    void add_prime(ulong p);

    void add_prime()
    {
        ulong p = ffts.back().mod.n;
        do {
            p = next_fft_number(p);
        } while (!n_is_prime(p));
        add_prime(p);
    }

    void my_mpn_mul(ulong* z, const ulong* a, ulong an, const ulong* b, ulong bn);
};

