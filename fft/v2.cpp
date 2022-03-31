#define BLK_SZ 256
#define LG_BLK_SZ 8
#define BLK_SHIFT 10
#define VEC_SZ 4

template<int n> struct is0 {};
template<> struct is0<0> {typedef bool type;};

template<int n> struct is1 {};
template<> struct is1<1> {typedef bool type;};

template<int n> struct is2 {};
template<> struct is2<2> {typedef bool type;};

template<int n> struct is3 {};
template<> struct is3<3> {typedef bool type;};

template<int n> struct is4 {};
template<> struct is4<4> {typedef bool type;};

template<int n> struct is5 {};
template<> struct is5<5> {typedef bool type;};

template<int n> struct is6 {};
template<> struct is6<6> {typedef bool type;};

template<int n> struct is7 {};
template<> struct is7<7> {typedef bool type;};

template<int n> struct is8 {};
template<> struct is8<8> {typedef bool type;};


struct fftv2_ctx {
    double * data;
    double p;
    double pinv;
    nmod_t mod;
    ulong primitive_root;
    ulong depth;
    ulong blk_sz;
    double * w2s;
    ulong wtab_depth;

    void init_prime(ulong p);

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
        return (I << LG_BLK_SZ) + (I >> BLK_SHIFT);
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

    double eval_poly(double a0, double a1, double b)
    {
        double x = a1;
        x = add(a0, mulmod2(x, b, p, pinv));
        return reduce_pm2n_to_pm1n(x, p);
    }

    double eval_poly(double a0, double a1, double a2, double a3, double b)
    {
        double x = a3;
        x = add(a2, mulmod2(x, b, p, pinv));
        x = add(a1, mulmod2(x, b, p, pinv));
        x = add(a0, mulmod2(x, b, p, pinv));
        return reduce_pm2n_to_pm1n(x, p);
    }

    double eval_poly(double* a, ulong an, double b)
    {
        double x = a[--an];
        while (an > 0)
            x = add(a[--an], mulmod2(x, b, p, pinv));
        return reduce_pm2n_to_pm1n(x, p);
    }


    inline void fft_basecase2_1x(double* X, ulong j);
    inline void fft_basecase2_2x(double* X, ulong j);
    inline void fft_basecase2_4x(double* X, ulong j);
    inline void ifft_basecase2_1x(double* X, ulong j);
    inline void ifft_basecase2_2x(double* X, ulong j);
    inline void ifft_basecase2_4x(double* X, ulong j);

    template <int> inline void fft_basecase(double* X, ulong j);
    template <int depth, bool, typename std::enable_if<depth == 0>::type* = nullptr> inline void ifft_basecase(double* X, ulong j);
    template <int depth, bool, typename std::enable_if<depth == 1>::type* = nullptr> inline void ifft_basecase(double* X, ulong j);
    template <int depth, bool, typename std::enable_if<depth == 2>::type* = nullptr> inline void ifft_basecase(double* X, ulong j);
    template <int depth, bool, typename std::enable_if<depth == 4>::type* = nullptr> inline void ifft_basecase(double* X, ulong j);
    template <int depth, bool, typename std::enable_if<depth == 6>::type* = nullptr> inline void ifft_basecase(double* X, ulong j);
    template <int depth, bool, typename std::enable_if<depth == 8>::type* = nullptr> inline void ifft_basecase(double* X, ulong j);
    void fft_base(ulong I, ulong j);
    void ifft_base(ulong I, ulong j);

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

void fftv2_ctx::init_prime(ulong pp)
{
    p = pp;
    pinv = 1.0/p;
    nmod_init(&mod, pp);
    primitive_root = n_primitive_root_prime(pp);

    // fill wtab to a depth of 6 (32 entries: 1, exp(pi*i/2), exp(pi*i/4), exp(pi*i*3/4), ...) 
    wtab_depth = 6;
    ulong N = pow2(wtab_depth-1);
    w2s = reinterpret_cast<double*>(std::malloc(N*sizeof(double)));

    ulong w = nmod_pow_ui(primitive_root, (mod.n - 1)>>wtab_depth, mod);
    ulong wi = 1;
    for (ulong j = 0; j < N; j++)
    {
        w2s[n_revbin(j, wtab_depth)/2] = reduce_0n_to_pmhn(wi, p);
        wi = nmod_mul(wi, w, mod);
    }
}

void fftv2_ctx::fit_wtab(ulong k)
{
    if (wtab_depth >= k)
        return;

    w2s = reinterpret_cast<double*>(std::realloc(w2s, pow2(k)*sizeof(double)));

    while (wtab_depth < k)
    {
        wtab_depth++;
        ulong N = pow2(wtab_depth-1);

        slong ww = nmod_pow_ui(primitive_root, (mod.n - 1)>>wtab_depth, mod);
        packed<double, VEC_SZ> w = reduce_0n_to_pmhn(double(ww), p, 0.5*p);
        packed<double, VEC_SZ> n = p;
        packed<double, VEC_SZ> ninv = pinv;
        packed<double, VEC_SZ> hn = 0.5*p;
        packed<double, VEC_SZ> x0, x1;
        double* wptr = w2s;
        ulong j = 0;
        do {
            x0.load(wptr + j);
            x1.load(wptr + j + VEC_SZ);
            x0 = mulmod2(x0, w, n, ninv);
            x1 = mulmod2(x1, w, n, ninv);
            x0 = reduce_pm1n_to_pmhn(x0, n, hn);
            x1 = reduce_pm1n_to_pmhn(x1, n, hn);
            x0.store(wptr + N/2 + j);
            x1.store(wptr + N/2 + j + VEC_SZ);
            j += 2*VEC_SZ;
        } while (j < N/2);
        assert(j == N/2);
    }
}


#include "basecases.cpp"

void fftv2_ctx::fft_base(ulong I, ulong j)
{
    fft_basecase<LG_BLK_SZ>(from_index(I), j);
}

void fftv2_ctx::ifft_base(ulong I, ulong j)
{
    ifft_basecase<LG_BLK_SZ, true>(from_index(I), j);
}



void fftv2_ctx::fft_main_block(
   ulong I, // starting index
   ulong S, // stride
   ulong k, // transform length 2^k
   ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        // column ffts
        for (ulong a = 0; a < ulong(1)<<k2; a++)
            fft_main_block(I + a*S, S<<k2, k1, j);

        // row ffts
        for (ulong b = 0; b < ulong(1)<<k1; b++)
            fft_main_block(I + (b<<k2)*S, S, k2, (j<<k1) + b);

        return;
    }

    if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
#if LG_BLK_SZ > 4
        PD_RADIX_4_FORWARD_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_4_FORWARD_2X(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
#else
        RADIX_4_FORWARD_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 1)
            RADIX_4_FORWARD(X0[i], X1[i], X2[i], X3[i]);
#endif
    }
    else if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
#if LG_BLK_SZ > 4
        PD_RADIX_2_FORWARD_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_2_FORWARD_2X(X0+i, X1+i, X0+i+4, X1+i+4);
#else
        RADIX_2_FORWARD_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 1)
            RADIX_2_FORWARD(X0+i, X1+i);
#endif
    }
}

void fftv2_ctx::ifft_main_block(
   ulong I, // starting index
   ulong S, // stride
   ulong k, // transform length 2^k
   ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        // row ffts
        for (ulong b = 0; b < ulong(1)<<k1; b++)
            ifft_main_block(I + (b<<k2)*S, S, k2, (j<<k1) + b);

        // column ffts
        for (ulong a = 0; a < ulong(1)<<k2; a++)
            ifft_main_block(I + a*S, S<<k2, k1, j);

        return;
    }

    if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
#if LG_BLK_SZ > 4
        PD_RADIX_4_REVERSE_PARAM(j, true)
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_4_REVERSE_2X(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
#else
        RADIX_4_REVERSE_PARAM(j, true)
        for (ulong i = 0; i < blk_sz; i += 1)
            RADIX_4_REVERSE(X0+i, X1+i, X2+i, X3+i);
#endif
    }
    else if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
#if LG_BLK_SZ > 4
        PD_RADIX_2_REVERSE_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_2_REVERSE_2X(X0+i, X1+i, X0+i+4, X1+i+4);
#else
        RADIX_2_REVERSE_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 1)
            RADIX_2_REVERSE(X0+i, X1+i);
#endif
    }
}

void fftv2_ctx::fft_main(
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k+LG_BLK_SZ)
    ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        for (ulong a = 0; a < ulong(1)<<k2; a++)
            fft_main_block(I + a*S, S<<k2, k1, j);

        for (ulong b = 0; b < ulong(1)<<k1; b++)
            fft_main(I + b*(S<<k2), S, k2, (j<<k1) + b);

        return;
    }

    if (k == 2)
    {
        // k1 = 2; k2 = 0
        fft_main_block(I, S, 2, j);
        fft_base(I+S*0, 4*j+0);
        fft_base(I+S*1, 4*j+1);
        fft_base(I+S*2, 4*j+2);
        fft_base(I+S*3, 4*j+3);
    }
    else if (k == 1)
    {
        // k1 = 1; k2 = 0
        fft_main_block(I, S, 1, j);
        fft_base(I+S*0, 2*j+0);
        fft_base(I+S*1, 2*j+1);
    }
    else
    {
        fft_base(I, j);
    }
}

void fftv2_ctx::ifft_main(
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k+LG_BLK_SZ)
    ulong j)
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        for (ulong b = 0; b < ulong(1)<<k1; b++)
            ifft_main(I + b*(S<<k2), S, k2, (j<<k1) + b);

        for (ulong a = 0; a < ulong(1)<<k2; a++)
            ifft_main_block(I + a*S, S<<k2, k1, j);

        return;
    }

    if (k == 2)
    {
        // k1 = 2; k2 = 0
        ifft_base(I+S*0, 4*j+0);
        ifft_base(I+S*1, 4*j+1);
        ifft_base(I+S*2, 4*j+2);
        ifft_base(I+S*3, 4*j+3);
        ifft_main_block(I, S, 2, j);
    }
    else if (k == 1)
    {
        // k1 = 1; k2 = 0
        ifft_base(I+S*0, 2*j+0);
        ifft_base(I+S*1, 2*j+1);
        ifft_main_block(I, S, 1, j);
    }
    else
    {
        ifft_base(I, j);
    }
}

//////////////// truncated fft /////////////////////////

void fftv2_ctx::fft_trunc_block(
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k)
    ulong j,
    ulong itrunc,
    ulong otrunc)
{
    assert(itrunc <= ulong(1) << k);
    assert(otrunc <= ulong(1) << k);

    if (otrunc < 1)
        return;

    if (itrunc <= 1)
    {
        if (itrunc < 1)
        {
            for (ulong a = 0; a < otrunc; a++)
            {
                double* X0 = from_index(I + S*a);
                packed<double, 4> z;
                z.zero();
                for (ulong i = 0; i < blk_sz; i += 4)
                    z.store(X0 + i);
            }
        }
        else
        {
            double* X0 = from_index(I + S*0);
            for (ulong a = 1; a < otrunc; a++)
            {
                double* X1 = from_index(I + S*a);
                for (ulong i = 0; i < blk_sz; i += 2*VEC_SZ)
                {
                    packed<double, VEC_SZ> x0, x1;
                    x0.load(X0 + i + 0*VEC_SZ);
                    x1.load(X0 + i + 1*VEC_SZ);
                    x0.store(X1 + i + 0*VEC_SZ);
                    x1.store(X1 + i + 1*VEC_SZ);
                }
            }
        }

        return;
    }

    if (itrunc == otrunc && otrunc == (ulong(1) << k))
    {
        fft_main_block(I, S, k, j);
        return;
    }

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = ulong(1) << k2;
        ulong n1 = otrunc >> k2;
        ulong n2 = otrunc & (l2 - 1);
        ulong z1 = itrunc >> k2;
        ulong z2 = itrunc & (l2 - 1);
        ulong n1p = n1 + (n2 != 0);
        ulong z2p = std::min(l2, itrunc);

        // columns
        for (ulong a = 0; a < z2p; a++)
            fft_trunc_block(I + a*S, S << k2, k1, j, z1 + (a < z2), n1p);

        // full rows
        for (ulong b = 0; b < n1; b++)
            fft_trunc_block(I + b*(S << k2), S, k2, (j << k1) + b, z2p, l2);

        // last partial row
        if (n2 > 0)
            fft_trunc_block(I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2);

        return;
    }

    if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
        PD_RADIX_4_FORWARD_PARAM(j)

//        hits[itrunc][otrunc]++;

        if (false){
        }else if (itrunc == 2 && otrunc == 4){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC4(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 2 && otrunc == 3){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC3(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 2 && otrunc == 2){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC2(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 2 && otrunc == 1){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC2_OTRUNC1(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);

        }else if (itrunc == 4 && otrunc == 3){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC4_OTRUNC3(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 4 && otrunc == 2){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC4_OTRUNC2(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 4 && otrunc == 1){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC4_OTRUNC1(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);

        }else if (itrunc == 3 && otrunc == 4){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC3_OTRUNC4(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 3 && otrunc == 3){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC3_OTRUNC3(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 3 && otrunc == 2){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC3_OTRUNC2(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        }else if (itrunc == 3 && otrunc == 1){
            for (ulong i = 0; i < blk_sz; i += 8)
                PD_RADIX_4_FORWARD_2X_ITRUNC3_OTRUNC1(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);

        }else
        {
            std::cout << "oops" << std::endl;
            std::abort();
        }
    }
    else if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        PD_RADIX_2_FORWARD_PARAM(j)

        assert(itrunc == 2 && otrunc == 1);
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_2_FORWARD_2X_ITRUNC2_OTRUNC1(X0+i, X1+i, X0+i+4, X1+i+4);
    }
}


void fftv2_ctx::fft_trunc(
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k + LG_BLK_SZ)
    ulong j,
    ulong itrunc,   // actual trunc is itrunc*BLK_SZ
    ulong otrunc)   // actual trunc is otrunc*BLK_SZ
{
    if (otrunc < 1)
        return;

    if (itrunc < 1)
    {
        for (ulong a = 0; a < otrunc; a++)
        {
            double* X0 = from_index(I + S*a);
            packed<double, 4> z;
            z.zero();
            for (ulong i = 0; i < blk_sz; i += 4)
                z.store(X0 + i);
        }

        return;
    }

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = ulong(1) << k2;
        ulong n1 = otrunc >> k2;
        ulong n2 = otrunc & (l2 - 1);
        ulong z1 = itrunc >> k2;
        ulong z2 = itrunc & (l2 - 1);
        ulong n1p = n1 + (n2 != 0);
        ulong z2p = std::min(l2, itrunc);

        // columns
        for (ulong a = 0; a < z2p; a++)
            fft_trunc_block(I + a*S, S << k2, k1, j, z1 + (a < z2), n1p);

        // full rows
        for (ulong b = 0; b < n1; b++)
            fft_trunc(I + b*(S << k2), S, k2, (j << k1) + b, z2p, l2);

        // last partial row
        if (n2 > 0)
            fft_trunc(I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2);

        return;
    }

    if (k == 2)
    {
        fft_trunc_block(I, S, 2, j, itrunc, otrunc);
        fft_base(I+S*0, 4*j+0);
        if (otrunc > 1) fft_base(I+S*1, 4*j+1);
        if (otrunc > 2) fft_base(I+S*2, 4*j+2);
        if (otrunc > 3) fft_base(I+S*3, 4*j+3);
    }
    else if (k == 1)
    {
        fft_trunc_block(I, S, 1, j, itrunc, otrunc);
        fft_base(I+S*0, 2*j+0);
        if (otrunc > 1) fft_base(I+S*1, 2*j+1);
    }
    else
    {
        fft_base(I, j);
    }
}



void fftv2_ctx::ifft_trunc_block(
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k)
    ulong j,
    ulong z,   // actual trunc is z
    ulong n,   // actual trunc is n
    bool f)
{
    assert(n <= z);
    assert(1 <= z && z <= pow2(k));
    assert(1 <= n+f && n+f <= pow2(k));

    if (!f && z == n && n == pow2(k))
    {
        ifft_main_block(I, S, k, j);
        return;
    }

    if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
        packed<double, VEC_SZ> N = p;
        packed<double, VEC_SZ> Ninv = pinv;

        /*
            common cases (z,n,f) =
               most common  (2,2,true)
               most common  (3,3,false|true)
                            (4,0,true)
                            (4,1,false|true)
                            (4,2,false|true)
                            (4,3,false|true)

            k = 2, z = 4, n = 2, f = true
            [            2                2        -w^2             0]
            [         2//w            -2//w           0          -w^2]
            [1//2*r + 1//2   -1//2*r + 1//2   -1//2*w^2   -1//2*r*w^3]

            k = 2, z = 4, n = 3, f = false
            [      -r + 1           r + 1         2   r*w^3]
            [        2//w           -2//w         0    -w^2]
            [(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]

            k = 2, z = 4, n = 3, f = true
            [      -r + 1           r + 1         2   r*w^3]
            [        2//w           -2//w         0    -w^2]
            [(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
            [          -r               r         1   r*w^3]
        */

        if (z == 2 && n == 2 && f)
        {
        /*
            k = 2, z = 2, n = 2, f = true
            [            2                2]
            [         2//w            -2//w]
            [1//2*r + 1//2   -1//2*r + 1//2]

            {x0, x1} = {2*(x0 + x1), 2*w^-1*(x0 - x1), (x0+x1)/2 + (x0-x1)*i/2}
        */
            ulong mask = saturate_bits(j);
            double W = ((j) == 0) ? -w2s[0] : w2s[(2*(j)  )^(mask)];
            packed<double, VEC_SZ> c0 = 2.0;
            packed<double, VEC_SZ> c1 = reduce_pm1n_to_pmhn(-2.0*W, p);
            packed<double, VEC_SZ> c2 = 0.5 - 0.5*p;
            packed<double, VEC_SZ> c3 = reduce_pm1n_to_pmhn(mulmod2(w2s[1], 0.5 - 0.5*p, p, pinv), p);

            for (ulong i = 0; i < blk_sz; i += 2*VEC_SZ)
            {
                packed<double, VEC_SZ> u0, v0, s0, t0, u1, v1, s1, t1;
                u0.load(X0 + i);
                    u1.load(X0 + i + VEC_SZ);
                v0.load(X1 + i);
                    v1.load(X1 + i + VEC_SZ);
                s0 = add(u0, v0);
                    s1 = add(u1, v1);
                t0 = sub(u0, v0);
                    t1 = sub(u1, v1);
                u0 = mulmod2(s0, c0, N, Ninv);
                    u1 = mulmod2(s1, c0, N, Ninv);
                v0 = mulmod2(t0, c1, N, Ninv);
                    v1 = mulmod2(t1, c1, N, Ninv);
                s0 = mulmod2(s0, c2, N, Ninv);
                    s1 = mulmod2(s1, c2, N, Ninv);
                t0 = mulmod2(t0, c3, N, Ninv);
                    t1 = mulmod2(t1, c3, N, Ninv);
                s0 = add(s0, t0);
                    s1 = add(s1, t1);
                u0.store(X0 + i);
                    u1.store(X0 + i + VEC_SZ);
                v0.store(X1 + i);
                    v1.store(X1 + i + VEC_SZ);
                s0.store(X2 + i);
                    s1.store(X2 + i + VEC_SZ);
            }

            return;
        }
        else if (z == 3 && n == 3 && !f)
        {
        /*
            k = 2, z = 3, n = 3, f = false
            [      -r + 1           r + 1         2]
            [        2//w           -2//w         0]
            [(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]

            {x0, x1, x3} = {        -r*(x0 - x1) + (x0 + x1) + 2*x2,
                            2*w^-1*(x0 - x1),
                             -w^-2*(-r*(x0 - x1) - (x0 + x1) + 2*x2)}

                         = {        2*x2 - r*v0 + u0,
                            2*w^-1*v0,
                             -w^-2*(2*x2 - r*v0 - u0)}
        */
            ulong mask = saturate_bits(j);
            double W = ((j) == 0) ? -w2s[0] : w2s[(2*(j)  )^(mask)];
            double W2 = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)];

            packed<double, VEC_SZ> f0 = w2s[1];                         // r
            packed<double, VEC_SZ> f1 = reduce_pm1n_to_pmhn(-2.0*W, p); // 2*w^-1
            packed<double, VEC_SZ> f2 = 2.0;
            packed<double, VEC_SZ> f3 = W2;                             // -w^-2

            for (ulong i = 0; i < blk_sz; i += 2*VEC_SZ)
            {
                packed<double, VEC_SZ> a0, b0, c0, u0, v0,  a1, b1, c1, u1, v1;
                a0.load(X0 + i);
                    a1.load(X0 + i + VEC_SZ);
                b0.load(X1 + i);
                    b1.load(X1 + i + VEC_SZ);
                c0.load(X2 + i);
                    c1.load(X2 + i + VEC_SZ);
                u0 = add(a0, b0);
                    u1 = add(a1, b1);
                v0 = sub(a0, b0);
                    v1 = sub(a1, b1);
                c0 = reduce_to_pm1n(c0, N, Ninv);
                    c1 = reduce_to_pm1n(c1, N, Ninv);
                u0 = reduce_to_pm1n(u0, N, Ninv);
                    u1 = reduce_to_pm1n(u1, N, Ninv);
                b0 = mulmod2(v0, f1, N, Ninv);
                    b1 = mulmod2(v1, f1, N, Ninv);
                v0 = mulmod2(v0, f0, N, Ninv);
                    v1 = mulmod2(v1, f0, N, Ninv);
                c0 = fmsub(f2, c0, v0);
                    c1 = fmsub(f2, c1, v1);
                a0 = add(c0, u0);
                    a1 = add(c1, u1);
                c0 = sub(c0, u0);                
                    c1 = sub(c1, u1);                
                c0 = mulmod2(c0, f3, N, Ninv);
                    c1 = mulmod2(c1, f3, N, Ninv);
                a0.store(X0 + i);
                    a1.store(X0 + i + VEC_SZ);
                b0.store(X1 + i);
                    b1.store(X1 + i + VEC_SZ);
                c0.store(X2 + i);
                    c1.store(X2 + i + VEC_SZ);
            }

            return;
        }
        else if (z == 3 && n == 3 && f)
        {
        /*
            k = 2, z = 3, n = 3, f = true
            [      -r + 1           r + 1         2]
            [        2//w           -2//w         0]
            [(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]
            [          -r               r         1]

            {x0, x1, x3, x4} = {        -r*(x0 - x1) + (x0 + x1) + 2*x2,
                                2*w^-1*(x0 - x1),
                                 -w^-2*(-r*(x0 - x1) - (x0 + x1) + 2*x2),
                                        -r*(x0 - x1)             +   x2  }
        */
            ulong mask = saturate_bits(j);
            double W = ((j) == 0) ? -w2s[0] : w2s[(2*(j)  )^(mask)];
            double W2 = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)];

            packed<double, VEC_SZ> f0 = w2s[1];                         // r
            packed<double, VEC_SZ> f1 = reduce_pm1n_to_pmhn(-2.0*W, p); // 2*w^-1
            packed<double, VEC_SZ> f2 = 2.0;
            packed<double, VEC_SZ> f3 = W2;                             // -w^-2

            for (ulong i = 0; i < blk_sz; i += 2*VEC_SZ)
            {
                packed<double, VEC_SZ> a0, b0, c0, d0, u0, v0,  a1, b1, c1, d1, u1, v1;
                a0.load(X0 + i);
                    a1.load(X0 + i + VEC_SZ);
                b0.load(X1 + i);
                    b1.load(X1 + i + VEC_SZ);
                c0.load(X2 + i);
                    c1.load(X2 + i + VEC_SZ);
                u0 = add(a0, b0);
                    u1 = add(a1, b1);
                v0 = sub(a0, b0);
                    v1 = sub(a1, b1);
                c0 = reduce_to_pm1n(c0, N, Ninv);
                    c1 = reduce_to_pm1n(c1, N, Ninv);
                u0 = reduce_to_pm1n(u0, N, Ninv);
                    u1 = reduce_to_pm1n(u1, N, Ninv);
                b0 = mulmod2(v0, f1, N, Ninv);
                    b1 = mulmod2(v1, f1, N, Ninv);
                v0 = mulmod2(v0, f0, N, Ninv);
                    v1 = mulmod2(v1, f0, N, Ninv);
                d0 = sub(c0, v0);
                    d1 = sub(c1, v1);
                c0 = fmsub(f2, c0, v0);
                    c1 = fmsub(f2, c1, v1);
                a0 = add(c0, u0);
                    a1 = add(c1, u1);
                c0 = sub(c0, u0);                
                    c1 = sub(c1, u1);                
                c0 = mulmod2(c0, f3, N, Ninv);
                    c1 = mulmod2(c1, f3, N, Ninv);
                a0.store(X0 + i);
                    a1.store(X0 + i + VEC_SZ);
                b0.store(X1 + i);
                    b1.store(X1 + i + VEC_SZ);
                c0.store(X2 + i);
                    c1.store(X2 + i + VEC_SZ);
                d0.store(X3 + i);
                    d1.store(X3 + i + VEC_SZ);
            }

            return;
        }
        else if (z == 4 && n == 1 && !f)
        {
        /*
            k = 2, z = 4, n = 1, f = false
            [4   -w   -w^2   -w^3]

            {x0} = {4*x0 - w*x1 - w^2*(x2 + w*x3)}
        */
            packed<double, VEC_SZ> f0 = 4.0;
            packed<double, VEC_SZ> w  = w2s[2*j];
            packed<double, VEC_SZ> w2 = w2s[j];

            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> a0, b0, c0, d0, a1, b1, c1, d1;
                d0.load(X3 + i);
                a0.load(X0 + i);
                b0.load(X1 + i);
                c0.load(X2 + i);
                d0 = mulmod2(d0, w, N, Ninv);
                a0 = mulmod2(a0, f0, N, Ninv);
                c0 = add(c0, d0);
                b0 = mulmod2(b0, w, N, Ninv);
                a0 = sub(a0, b0);
                c0 = mulmod2(c0, w2, N, Ninv);
                a0 = sub(a0, c0);
                a0.store(X0 + i);
            }

            return;
        }
        else if (z == 4 && n == 0 && f)
        {
        /*
            k = 2, z = 4, n = 0, f = true
            [1//4   1//4*w   1//4*w^2   1//4*w^3]

            {x0} = {1/4*x0 + w/4*x1 + w^2/4*(x2 + w*x3)}
        */
            double ff0 = 0.25 - 0.25*p;
            double ww = w2s[2*j];
            double ww2 = w2s[j];
            packed<double, VEC_SZ> f0 = ff0;
            packed<double, VEC_SZ>  wo4 = reduce_pm1n_to_pmhn(mulmod2(ww, ff0, p, pinv), p);
            packed<double, VEC_SZ> w2o4 = reduce_pm1n_to_pmhn(mulmod2(ww2, ff0, p, pinv), p);
            packed<double, VEC_SZ> w = ww;

            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> a0, b0, c0, d0, a1, b1, c1, d1;
                d0.load(X3 + i);
                a0.load(X0 + i);
                b0.load(X1 + i);
                c0.load(X2 + i);
                d0 = mulmod2(d0, w, N, Ninv);
                a0 = mulmod2(a0, f0, N, Ninv);
                c0 = add(c0, d0);
                b0 = mulmod2(b0, wo4, N, Ninv);
                a0 = add(a0, b0);
                c0 = mulmod2(c0, w2o4, N, Ninv);
                a0 = add(a0, c0);
                a0.store(X0 + i);
            }

            return;            
        }
        else if (z == 4 && n == 2 && !f)
        {
        /*
            k = 2, z = 4, n = 2, f = false
            [   2       2   -w^2      0]
            [2//w   -2//w      0   -w^2]
        */
            ulong mask = saturate_bits(j);
            double mwi = ((j) == 0) ? -w2s[0] : w2s[(2*(j)  )^(mask)];
            packed<double, VEC_SZ> f0 = 2.0;
            packed<double, VEC_SZ> w2 = w2s[j];
            packed<double, VEC_SZ> twowi = reduce_pm1n_to_pmhn(-2.0*mwi, p);

            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> a0, b0, c0, d0, u0, v0;
                a0.load(X0 + i);
                b0.load(X1 + i);
                c0.load(X2 + i);
                d0.load(X3 + i);
                u0 = add(a0, b0);
                v0 = sub(a0, b0);
                u0 = mulmod2(u0, f0, N, Ninv);
                v0 = mulmod2(v0, twowi, N, Ninv);
                c0 = mulmod2(c0, w2, N, Ninv);
                d0 = mulmod2(d0, w2, N, Ninv);
                u0 = sub(u0, c0);
                v0 = sub(v0, d0);
                u0.store(X0 + i);
                v0.store(X1 + i);
            }

            return;            
        }
        else if (z == 4 && n == 1 && f)
        {
        /*
            k = 2, z = 4, n = 1, f = true
            [4        -w   -w^2        -w^3]
            [1   -1//2*w      0   -1//2*w^3]

            2(2*x0 - w/2*(x1 + w^2*x3)) - w^2*x2
                x0 - w/2*(x1 + w^2*x3)
        */
            packed<double, VEC_SZ> f2 = 2.0;
            packed<double, VEC_SZ> w2 = w2s[j];
            packed<double, VEC_SZ> wo2 = reduce_pm1n_to_pmhn(mulmod2(w2s[2*j], 0.5-0.5*p, p, pinv), p);
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> a0, b0, c0, d0, u0;
                a0.load(X0 + i);
                a0 = reduce_to_pm1n(a0, N, Ninv);
                b0.load(X1 + i);
                c0.load(X2 + i);
                d0.load(X3 + i);
                c0 = mulmod2(c0, w2, N, Ninv);
                d0 = mulmod2(d0, w2, N, Ninv);
                b0 = add(b0, d0);
                b0 = mulmod2(b0, wo2, N, Ninv);
                u0 = fmsub(f2, a0, b0);
                b0 = sub(a0, b0);
                a0 = reduce_to_pm1n(add(u0, u0), N, Ninv);
                a0 = sub(a0, c0);
                a0.store(X0 + i);
                b0.store(X1 + i);
            }

            return;
        }

//hits[z][2*n+f]++;

    }


    if (k > 1)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = ulong(1) << k2;
        ulong n1 = n >> k2;
        ulong n2 = n & (l2 - 1);
        ulong z1 = z >> k2;
        ulong z2 = z & (l2 - 1);
        bool fp = n2 + f > 0;
        ulong z2p = std::min(l2, z);
        ulong m = std::min(n2, z2);
        ulong mp = std::max(n2, z2);

        // complete rows
        for (ulong b = 0; b < n1; b++)
            ifft_main_block(I + b*(S << k2), S, k2, (j << k1) + b);

        // rightmost columns
        for (ulong a = n2; a < z2p; a++)
            ifft_trunc_block(I + a*S, S << k2, k1, j, z1 + (a < mp), n1, fp);

        // last partial row
        if (fp)
            ifft_trunc_block(I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2, f);

        // leftmost columns
        for (ulong a = 0; a < n2; a++)
            ifft_trunc_block(I + a*S, S << k2, k1, j, z1 + (a < m), n1 + 1, false);

        return;
    }

    if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);

        packed<double, VEC_SZ> N = p;
        packed<double, VEC_SZ> Ninv = pinv;

        if (n == 1 && z == 2 && f)
        {
            // {x0, x1} = {2*x0 - w*x1, x0 - w*x1}
            packed<double, VEC_SZ> w = w2s[j];
            packed<double, VEC_SZ> c = 2.0;
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> u0, u1, v0, v1;
                u0.load(X0 + i);
                u1.load(X1 + i);
                u0 = reduce_to_pm1n(u0, N, Ninv);
                u1 = mulmod2(u1, w, N, Ninv);
                v0 = fmsub(c, u0, u1);
                v1 = sub(u0, u1);
                v0.store(X0 + i);
                v1.store(X1 + i);
            }
        }
        else if (n == 1 && z == 2 && !f)
        {
            // {x0} = {2*x0 - w*x1}
            packed<double, VEC_SZ> w = w2s[j];
            packed<double, VEC_SZ> c = 2.0;
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> u0, u1, v0, v1;
                u0.load(X0 + i);
                u1.load(X1 + i);
                u0 = reduce_to_pm1n(u0, N, Ninv);
                u1 = mulmod2(u1, w, N, Ninv);
                v0 = fmsub(c, u0, u1);
                v0.store(X0 + i);
            }
        }
        else if (n == 1 && z != 2 && f)
        {
            // {x0, x1} = {2*x0, x0}
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> u0, v0, v1;
                u0.load(X0 + i);
                u0 = reduce_to_pm1n(u0, N, Ninv);
                v0 = add(u0, u0);
                v1 = u0;
                v0.store(X0 + i);
                v1.store(X1 + i);
            }
        }
        else if (n == 1 && z != 2 && !f)
        {
            // {x0} = {2*x0}
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> u0, v0;
                u0.load(X0 + i);
                u0 = reduce_to_pm1n(u0, N, Ninv);
                v0 = add(u0, u0);
                v0.store(X0 + i);
            }
        }
        else if (n == 0 && z == 2 && f)
        {
            // {x0} = {(x0 + w*x1)/2}
            packed<double, VEC_SZ> w = w2s[j];
            packed<double, VEC_SZ> c = 0.5 - 0.5*p;
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> u0, u1;
                u0.load(X0 + i);
                u1.load(X1 + i);
                u1 = mulmod2(u1, w, N, Ninv);
                u0 = mulmod2(add(u0, u1), c, p, pinv);
                u0.store(X0 + i);
            }
        }
        else if (n == 0 && z != 2 && f)
        {
            // {x0} = {x0/2}
            packed<double, VEC_SZ> w = w2s[j];
            packed<double, VEC_SZ> c = 0.5 - 0.5*p;
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double, VEC_SZ> u0, u1;
                u0.load(X0 + i);
                u0 = mulmod2(u0, c, p, pinv);
                u0.store(X0 + i);
            }
        }
        else
        {
            std::cout << "ooops" << std::endl;
            std::abort();
        }

        return;
    }
}

void fftv2_ctx::ifft_trunc(
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k + LG_BLK_SZ)
    ulong j,
    ulong z,   // actual trunc is z*BLK_SZ
    ulong n,   // actual trunc is n*BLK_SZ
    bool f)
{
    assert(n <= z);
    assert(1 <= BLK_SZ*z && BLK_SZ*z <= pow2(k+LG_BLK_SZ));
    assert(1 <= BLK_SZ*n+f && BLK_SZ*n+f <= pow2(k+LG_BLK_SZ));
/*
std::cout << std::endl;
std::cout << "ifft_trunc called k = " << k << std::endl;
std::cout << "j = " << j << ",  ";
std::cout << "z = " << z << ",  ";
std::cout << "n = " << n << ",  ";
std::cout << "f = " << f << std::endl;
*/
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = ulong(1) << k2;
        ulong n1 = n >> k2;
        ulong n2 = n & (l2 - 1);
        ulong z1 = z >> k2;
        ulong z2 = z & (l2 - 1);
        bool fp = n2 + f > 0;
        ulong z2p = std::min(l2, z);
        ulong m = std::min(n2, z2);
        ulong mp = std::max(n2, z2);

        // complete rows
        for (ulong b = 0; b < n1; b++)
            ifft_main(I + b*(S << k2), S, k2, (j << k1) + b);

        // rightmost columns
        for (ulong a = n2; a < z2p; a++)
            ifft_trunc_block(I + a*S, S << k2, k1, j, z1 + (a < mp), n1, fp);

        // last partial row
        if (fp)
            ifft_trunc(I + n1*(S << k2), S, k2, (j << k1) + n1, z2p, n2, f);

        // leftmost columns
        for (ulong a = 0; a < n2; a++)
            ifft_trunc_block(I + a*S, S << k2, k1, j, z1 + (a < m), n1 + 1, false);

        return;
    }

    if (k == 2)
    {
                   ifft_base(I+S*0, 4*j+0);
        if (n > 1) ifft_base(I+S*1, 4*j+1);
        if (n > 2) ifft_base(I+S*2, 4*j+2);
        if (n > 3) ifft_base(I+S*3, 4*j+3);
        ifft_trunc_block(I, S, 2, j, z, n, f);
        if (f) ifft_trunc(I + S*n, S, 0, 4*j+n, 1, 0, f);
        
    }
    else if (k == 1)
    {
                   ifft_base(I+S*0, 2*j+0);
        if (n > 1) ifft_base(I+S*1, 2*j+1);
        ifft_trunc_block(I, S, 1, j, z, n, f);
        if (f) ifft_trunc(I + S*n, S, 0, 2*j+n, 1, 0, f);
    }
    else
    {
        assert(!f);
        ifft_base(I, j);
    }
}


///////////////////////// mpn_mul //////////////////////////////////////


// pointwise mul of self with b and m
void fftv2_ctx::point_mul(const double* b, ulong mm)
{
    double m = mm;
    if (m > 0.5*p)
        m -= p;

    packed<double, 4> M = m;
    packed<double, 4> n = p;
    packed<double, 4> ninv = pinv;
    for (ulong I = 0; I < ulong(1) << (depth - LG_BLK_SZ); I++)
    {
        double* x = from_index(I);
        const double* bx = b + offset(I);
        ulong j = 0;
        do {
            packed<double, 4> x0, x1, x2, x3, b0, b1, b2, b3;
            x0.load(x+j+0);
            x1.load(x+j+4);
            x2.load(x+j+8);
            x3.load(x+j+12);
            b0.load(bx+j+0);
            b1.load(bx+j+4);
            b2.load(bx+j+8);
            b3.load(bx+j+12);
            x0 = mulmod2(x0, M, n, ninv);
            x1 = mulmod2(x1, M, n, ninv);
            x2 = mulmod2(x2, M, n, ninv);
            x3 = mulmod2(x3, M, n, ninv);
            x0 = mulmod2(x0, b0, n, ninv);
            x1 = mulmod2(x1, b1, n, ninv);
            x2 = mulmod2(x2, b2, n, ninv);
            x3 = mulmod2(x3, b3, n, ninv);
            x0.store(x+j+0);
            x1.store(x+j+4);
            x2.store(x+j+8);
            x3.store(x+j+12);
        } while (j += 16, j < blk_sz);
    }
}


// TODO parameterize by bits and number of primes
void fftv2_ctx::from_mpn(
    const ulong * a_,
    ulong an_,
    ulong bits)
{
    bits = 71;
    FLINT_ASSERT(bits >= FLINT_BITS);

    const uint32_t* a = reinterpret_cast<const uint32_t*>(a_);
    ulong an = 2*an_;
    ulong xn = UWORD(1) << depth;
    ulong ttlen = bits + 64;

    double* two_pow = new double[ttlen];

    two_pow[0] = 1;
    ulong tp = 1;
    for (ulong i = 1; i < ttlen; i++)
    {
        tp = nmod_add(tp, tp, mod);
        two_pow[i] = (mod.n - tp < tp) ? -double(mod.n-tp) : double(tp);
    }

    ulong masks[32];
    masks[0] = 0x00000000;
    for (ulong i = 1; i < 32; i++)
        masks[i] = 2*masks[i-1]+1;

    // if i*bits + 32 < 32*an, then aindex is easy
    ulong end_easy = std::min(xn, (32*an - 33)/bits);

    // if i*bits >= 32*an, then aindex is zero
    ulong end_hard = std::min(xn, (32*an + bits - 1)/bits);

    ulong i = 0;

#define CODE(ir) \
    { \
        ulong k = ((i+ir)*bits)/32; \
        ulong j = ((  ir)*bits)%32; \
        double x = double(a[k] >> j); \
        k++; \
        j = 32 - j; \
        while (j + 32 <= bits) \
        { \
            x += mulmod2(double(a[k]), two_pow[j], p, pinv); \
            k++; \
            j += 32; \
        } \
        assert(j <= bits < j + 32); \
        x += mulmod2(double(a[k] & masks[bits - j]), two_pow[j], p, pinv); \
        set_index(i+ir, reduce_to_pm1n(x, p, pinv)); \
    }

    if ((bits % 4) == 0)
    {
        end_easy &= -uint32_t(8);
        for ( ; i < end_easy; i += 8)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
        }
    }
    else if ((bits % 2) == 0)
    {
        end_easy &= -uint32_t(16);
        for ( ; i < end_easy; i += 16)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
            CODE(8);CODE(9);CODE(10);CODE(11);
            CODE(12);CODE(13);CODE(14);CODE(15);
        }
    }
    else
    {
        end_easy &= -uint32_t(32);
        for ( ; i < end_easy; i += 32)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
            CODE(8);CODE(9);CODE(10);CODE(11);
            CODE(12);CODE(13);CODE(14);CODE(15);
            CODE(16);CODE(17);CODE(18);CODE(19);
            CODE(20);CODE(21);CODE(22);CODE(23);
            CODE(24);CODE(25);CODE(26);CODE(27);
            CODE(28);CODE(29);CODE(30);CODE(31);
        }
    }
#undef CODE

    for (; i < end_hard; i++)
    {
        ulong k = (i*bits)/32;
        ulong j = (i*bits)%32;

#define aindex(i) (((i) < an) ? a[i] : UWORD(0))

        double x = double(aindex(k) >> j);
        k++;
        j = 32 - j;
        while (j + 32 <= bits)
        {
            x += mulmod2(double(aindex(k)), two_pow[j], p, pinv);
            k++;
            j += 32;
        }
        x += mulmod2(double(aindex(k) & masks[bits - j]), two_pow[j], p, pinv);

#undef aindex

        set_index(i, reduce_to_pm1n(x, p, pinv));
//std::cout << "x[" << i << "]: " << get_index(i) << std::endl;
    }

    for (; i < xn; i++)
    {
        set_index(i, 0.0);
    }

    delete two_pow;
}


struct crt_data {
    double _prod_primes_d;
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

    double& prod_primes_d() {return _prod_primes_d;}
};


struct mpn_ctx_v2 {
    std::vector<fftv2_ctx> ffts;
    std::vector<crt_data> crts;
    double* double_buffer = nullptr;
    ulong double_buffer_alloc = 0;

    ~mpn_ctx_v2() {
        free(double_buffer);
    }


    double* fit_double_buffer(ulong n)
    {
        if (n > double_buffer_alloc)
        {
            std::free(double_buffer);
            double_buffer = reinterpret_cast<double*>(std::calloc(n, sizeof(double)));
            double_buffer_alloc = n;
        }
        return double_buffer;
    }


    ulong nprimes() {
        return ffts.size();
    }

    mpn_ctx_v2(ulong p) : double_buffer(nullptr), double_buffer_alloc(0)
    {
        ffts.emplace_back(p);
        crts.emplace_back(p, 1, 1);

        crts[0].co_prime_red(0) = 1;
        crts[0].co_prime(0)[0] = 1;
        crts[0].prod_primes()[0] = p;
        crts[0].prod_primes_d() = p;
    }

    void add_prime(ulong p)
    {
        ffts.emplace_back(p);
        ulong len = crts.back().coeff_len;
        ulong* t = new ulong[len + 2];

        t[len + 1] = 0;
        t[len] = mpn_mul_1(t, crts.back().prod_primes(), len, p);
        len += (t[len] != 0);

        mpz_t dummy;
        dummy->_mp_d = t;
        dummy->_mp_size = len;
        dummy->_mp_alloc = len;
        double dd = mpz_get_d(dummy);

        // leave room for one more bit
        if (FLINT_SIGN_EXT(t[len-1]) != 0)
            len += 1;

        crts.emplace_back(p, len, nprimes());

        // set product of primes
        crts.back().prod_primes_d() = dd;
        mpn_copyi(crts.back().prod_primes(), t, len);

        // set cofactors
        for (ulong i = 0; i < nprimes(); i++)
        {
            mpn_divexact_1(crts.back().co_prime(i), t, len, crts[i].prime);
            crts.back().co_prime_red(i) = mpn_mod_1(crts.back().co_prime(i), len, crts[i].prime);
        }

        delete[] t;
    }

    void add_prime()
    {
        ulong p = ffts.back().mod.n;
        do {
            p = next_fft_number(p);
        } while (!n_is_prime(p));
        add_prime(p);
    }
};



void mpn_mul_v2(
    mpn_ctx_v2& Q,
    ulong* z,
    const ulong* a, ulong an,
    const ulong* b, ulong bn)
{
    ulong zn = an + bn;
    ulong bits = 71;
    ulong alen = cld(FLINT_BITS * an, bits);
    ulong blen = cld(FLINT_BITS * bn, bits);
    alen = FLINT_MAX(1, alen);
    alen = FLINT_MAX(1, alen);
    ulong zlen = alen + blen - 1;
    ulong depth = clog2(zlen);
    depth = FLINT_MAX(depth, LG_BLK_SZ);
    ulong np = 4;
timeit_t timer;

    while (Q.nprimes() < np)
        Q.add_prime();

    Q.ffts[0].set_depth(depth);
    ulong zcoeff_len = Q.crts.back().coeff_len;

    ulong* tt = new ulong[zcoeff_len + 1];
    double* abuf = new double[Q.ffts[0].data_size()];
    double* bbuf = new double[Q.ffts[0].data_size()];
    ulong* zbuf = new ulong[zlen*zcoeff_len];

std::cout << "depth: " << depth << std::endl;
std::cout << " alen: " << alen << std::endl;
std::cout << " blen: " << blen << std::endl;

    for (ulong i = 0; i < np; i++)
    {
        fftv2_ctx& ctx = Q.ffts[i];
std::cout << "doing prime " << i << ": " << format_hex(ctx.mod.n) << std::endl;
        ctx.set_depth(depth);

        ctx.set_data(bbuf);
timeit_start(timer);
        ctx.from_mpn(b, bn, bits);
timeit_stop(timer);
std::cout << "b from mpn: " << timer->wall << std::endl;
timeit_start(timer);
        ctx.fft();
timeit_stop(timer);
std::cout << " b forward: " << timer->wall << std::endl;
        ctx.release_data();

        ctx.set_data(abuf);
timeit_start(timer);
        ctx.from_mpn(a, an, bits);
timeit_stop(timer);
std::cout << "a from mpn: " << timer->wall << std::endl;
timeit_start(timer);
        ctx.fft();
timeit_stop(timer);
std::cout << " a forward: " << timer->wall << std::endl;

timeit_start(timer);
        // bake in 2^-depth * (p2*p3*p4)^-1 mod p1
        {
            ulong t1, thi, tlo;
            ulong cop = Q.crts[np - 1].co_prime_red(i);
            thi = cop >> (FLINT_BITS - depth);
            tlo = cop << (depth);
            NMOD_RED2(t1, thi, tlo, ctx.mod);
            ctx.point_mul(bbuf, nmod_inv(t1, ctx.mod));
        }
timeit_stop(timer);
std::cout << " point mul: " << timer->wall << std::endl;

timeit_start(timer);
        ctx.ifft();
timeit_stop(timer);
std::cout << "   reverse: " << timer->wall << std::endl;


timeit_start(timer);
        ulong* mult = Q.crts[np - 1].co_prime(i);
        for (ulong j = 0; j < zlen; j++)
        {
            double x = ctx.get_index(j);
            slong xx = reduce_to_pm1n(x, ctx.p, ctx.pinv);
            ulong y;
            if (xx < 0)
                y = ctx.mod.n + xx;
            else
                y = xx;

            if (i == 0)
                mpn_mul_1(zbuf + j*zcoeff_len, mult, zcoeff_len, y);
            else
                mpn_addmul_1(zbuf + j*zcoeff_len, mult, zcoeff_len, y);
        }

timeit_stop(timer);
std::cout << "       crt: " << timer->wall << std::endl;


        ctx.release_data();
    }

timeit_start(timer);

    mpn_zero(z, zn);

    ulong * limit = Q.crts[np - 1].prod_primes();
    for (ulong i = 0; i < zlen; i++)
    {
        ulong tbits = i*bits;
        ulong* t = zbuf + i*zcoeff_len;

        while (mpn_cmp(t, limit, zcoeff_len) >= 0)
            mpn_sub_n(t, t, limit, zcoeff_len);

        ulong toff = tbits/FLINT_BITS;
        ulong tshift = tbits%FLINT_BITS;

        if (toff >= zn)
            break;

        if (tshift == 0)
        {
            mpn_add_n(z + toff, z + toff, t, FLINT_MIN(zcoeff_len, zn - toff));
        }
        else
        {
            tt[zcoeff_len] = mpn_lshift(tt, t, zcoeff_len, tshift);
            mpn_add_n(z + toff, z + toff, tt, FLINT_MIN(zcoeff_len + 1, zn - toff));
        }
    }

timeit_stop(timer);
std::cout << "     final: " << timer->wall << std::endl;


    delete[] abuf;
    delete[] bbuf;
    delete[] zbuf;
    delete[] tt;
}



template<ulong np, ulong bits>
void mpn_to_ffts(
    mpn_ctx_v2& Q,
    const ulong* a_, ulong an_, ulong atrunc,
    const packed<double, VEC_SZ>* two_pow)
{
    ulong nvs = (np + VEC_SZ - 1)/VEC_SZ;

//std::cout << "mpn_to_ffts called" << std::endl;

    FLINT_ASSERT(bits >= FLINT_BITS);

    const uint32_t* a = reinterpret_cast<const uint32_t*>(a_);
    ulong an = 2*an_;

    packed<double, VEC_SZ> X[nvs];
    packed<double, VEC_SZ> P[nvs];
    packed<double, VEC_SZ> PINV[nvs];

    for (ulong l = 0; l < nvs; l++)
    {
#if VEC_SZ == 4
        P[l] = packed<double, 4>(Q.ffts[4*l+0].p, Q.ffts[4*l+1].p, Q.ffts[4*l+2].p, Q.ffts[4*l+3].p);
        PINV[l] = packed<double, 4>(Q.ffts[4*l+0].pinv, Q.ffts[4*l+1].pinv, Q.ffts[4*l+2].pinv, Q.ffts[4*l+3].pinv);
#else
    #error "unsupported VEC_SZ"
#endif
    }

    // if i*bits + 32 < 32*an, then aindex is easy
    ulong end_easy = std::min(atrunc, (32*an - 33)/bits);

    // if i*bits >= 32*an, then aindex is zero
    ulong end_hard = std::min(atrunc, (32*an + bits - 1)/bits);


    ulong i = 0;

#define CODE(ir)\
    {\
        ulong k = ((i+ir)*bits)/32;\
        ulong j = ((  ir)*bits)%32;\
\
        packed<double, 4> ak = double(a[k] >> j);\
        for (ulong l = 0; l < nvs; l++)\
            X[l] = ak;\
        k++;\
        j = 32 - j;\
        while (j + 32 <= bits)\
        {\
            ak = double(a[k]);\
            for (ulong l = 0; l < nvs; l++)\
                X[l] = add(X[l], mulmod2(ak, two_pow[j*nvs+l], P[l], PINV[l]));\
            k++;\
            j += 32;\
        }\
\
        if ((bits-j) != 0)\
        {\
            ak = double(a[k] << (32-(bits-j)));\
            for (ulong l = 0; l < nvs; l++)\
                X[l] = add(X[l], mulmod2(ak, two_pow[(bits-32)*nvs+l], P[l], PINV[l]));\
        }\
\
        for (ulong l = 0; l < nvs; l++)\
            X[l] = reduce_to_pm1n(X[l], P[l], PINV[l]);\
\
        for (ulong l = 0; l < np; l++)\
            Q.ffts[l].set_index(i+ir, X[l/4].data[l%4]);\
    }

    if ((bits % 4) == 0)
    {
        end_easy &= -ulong(8);
        for ( ; i < end_easy; i += 8)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
        }
    }
    else if ((bits % 2) == 0)
    {
        end_easy &= -ulong(16);
        for ( ; i < end_easy; i += 16)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
            CODE(8);CODE(9);CODE(10);CODE(11);
            CODE(12);CODE(13);CODE(14);CODE(15);
        }
    }
    else
    {
        end_easy &= -ulong(32);
        for ( ; i < end_easy; i += 32)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
            CODE(8);CODE(9);CODE(10);CODE(11);
            CODE(12);CODE(13);CODE(14);CODE(15);
            CODE(16);CODE(17);CODE(18);CODE(19);
            CODE(20);CODE(21);CODE(22);CODE(23);
            CODE(24);CODE(25);CODE(26);CODE(27);
            CODE(28);CODE(29);CODE(30);CODE(31);
        }
    }
#undef CODE

    for (; i < end_hard; i++)
    {
        ulong k = (i*bits)/32;
        ulong j = (i*bits)%32;

#define aindex(i) (((i) < an) ? a[i] : uint32_t(0))

        packed<double, VEC_SZ> ak = double(aindex(k) >> j);
        for (ulong l = 0; l < nvs; l++)
            X[l] = ak;
        k++;
        j = 32 - j;
        while (j + 32 <= bits)
        {
            ak = double(aindex(k));
            for (ulong l = 0; l < nvs; l++)
                X[l] = add(X[l], mulmod2(ak, two_pow[j*nvs+l], P[l], PINV[l]));
            k++;
            j += 32;
        }

        if ((bits-j) != 0)
        {
            ak = double(aindex(k) << (32-(bits-j)));
            for (ulong l = 0; l < nvs; l++)
                X[l] = add(X[l], mulmod2(ak, two_pow[(bits-32)*nvs+l], P[l], PINV[l]));
        }

#undef aindex

        for (ulong l = 0; l < nvs; l++)
            X[l] = reduce_to_pm1n(X[l], P[l], PINV[l]);

        for (ulong l = 0; l < np; l++)
            Q.ffts[l].set_index(i, X[l/VEC_SZ].data[l%VEC_SZ]);
    }

    for (ulong l = 0; l < np; l++)
        for (ulong j = i; j < atrunc; j++)
            Q.ffts[l].set_index(j, 0.0);
}

void fill_two_pow(
    packed<double, VEC_SZ>* two_pow,
    ulong bits,
    std::vector<fftv2_ctx>& Qffts,
    ulong np)
{
    ulong nvs = (np + VEC_SZ - 1)/VEC_SZ;

    packed<double, VEC_SZ>* ps = new packed<double, VEC_SZ>[nvs];
    packed<double, VEC_SZ>* pinvs = new packed<double, VEC_SZ>[nvs];

    for (ulong l = 0; l < nvs; l++)
    {
#if VEC_SZ == 4
        ps[l] = packed<double, 4>(Qffts[4*l+0].p, Qffts[4*l+1].p, Qffts[4*l+2].p, Qffts[4*l+3].p);
        pinvs[l] = packed<double, 4>(Qffts[4*l+0].pinv, Qffts[4*l+1].pinv, Qffts[4*l+2].pinv, Qffts[4*l+3].pinv);
#else
    #error "unsuported VEC_SZ"
#endif
    }

    for (ulong l = 0; l < nvs; l++)
        two_pow[0*nvs+l] = 1;

    for (ulong i = 1; i < bits; i++)
    for (ulong l = 0; l < nvs; l++)
    {
        packed<double, VEC_SZ> t = two_pow[(i-1)*nvs+l];
        two_pow[i*nvs+l] = reduce_to_pm1n(add(t, t), ps[l], pinvs[l]);
    }

    delete[] ps;
    delete[] pinvs;
}

template <ulong np, ulong zcoeff_len>
void mpn_from_ffts(
    ulong* z, ulong zn, ulong zlen,
    mpn_ctx_v2& Q,
    ulong bits)
{
    ulong r[zcoeff_len + 1];
    ulong t[zcoeff_len + 1];
    ulong u, v;

    if (zcoeff_len != Q.crts[np-1].coeff_len)
    {
        std::cout << "oops" << std::endl;
        std::abort();
    }

    mpn_zero(z, zn);

    ulong i = 0;

    // easy if zn-zcoeff_len > floor(i*bits/64)
    ulong end_easy = (zn >= zcoeff_len+1 ? zn - (zcoeff_len+1) : ulong(0))*FLINT_BITS/bits;

    ulong Xs[BLK_SZ*np];

    end_easy &= -BLK_SZ;

    for (; i < end_easy; i += BLK_SZ)
    {
        ulong I = i/BLK_SZ;

        for (ulong l = 0; l < np; l++)
        {
            packed<double, VEC_SZ> P = Q.ffts[l].p;
            packed<double, VEC_SZ> PINV = Q.ffts[l].pinv;
            double* x = Q.ffts[l].from_index(I);
            for (ulong j = 0; j < BLK_SZ; j += 4*VEC_SZ)
            {
                packed<double, VEC_SZ> x0, x1, x2, x3;
                packed<ulong, VEC_SZ> y0, y1, y2, y3;
                x0.load(x + j + 0*VEC_SZ);
                x1.load(x + j + 1*VEC_SZ);
                x2.load(x + j + 2*VEC_SZ);
                x3.load(x + j + 3*VEC_SZ);
                x0 = reduce_to_pm1n(x0, P, PINV);
                x1 = reduce_to_pm1n(x1, P, PINV);
                x2 = reduce_to_pm1n(x2, P, PINV);
                x3 = reduce_to_pm1n(x3, P, PINV);
                x0 = reduce_pm1no_to_0n(x0, P);
                x1 = reduce_pm1no_to_0n(x1, P);
                x2 = reduce_pm1no_to_0n(x2, P);
                x3 = reduce_pm1no_to_0n(x3, P);
                y0 = convert_limited<packed<ulong, VEC_SZ>>(x0);
                y1 = convert_limited<packed<ulong, VEC_SZ>>(x1);
                y2 = convert_limited<packed<ulong, VEC_SZ>>(x2);
                y3 = convert_limited<packed<ulong, VEC_SZ>>(x3);
                y0.store(Xs + l*BLK_SZ + j + 0*VEC_SZ);
                y1.store(Xs + l*BLK_SZ + j + 1*VEC_SZ);
                y2.store(Xs + l*BLK_SZ + j + 2*VEC_SZ);
                y3.store(Xs + l*BLK_SZ + j + 3*VEC_SZ);
            }
        }

        for (ulong j = 0; j < BLK_SZ; j += 1)
        {
            for (ulong l = 0; l < np; l++)
            {
                ulong* mult = Q.crts[np - 1].co_prime(l);
                ulong y = Xs[l*BLK_SZ + j];
                if (l == 0)
                {
                    umul_ppmm(r[1], r[0], mult[0], y);
                    umul_ppmm(r[3], r[2], mult[2], y);
                    umul_ppmm(u, v, mult[1], y);
                    add_sssaaaaaa(r[3],r[2],r[1], r[3],r[2],r[1], 0,u,v);
                }
                else
                {
                    umul_ppmm(t[1], t[0], mult[0], y);
                    umul_ppmm(t[3], t[2], mult[2], y);
                    umul_ppmm(u, v, mult[1], y);
                    add_sssaaaaaa(t[3],t[2],t[1], t[3],t[2],t[1], 0,u,v);
                    add_ssssaaaaaaaa(r[3],r[2],r[1],r[0], r[3],r[2],r[1],r[0], t[3],t[2],t[1],t[0]);
                }
            }

            ulong* limit = Q.crts[np - 1].prod_primes();

        check1:
            if (r[3] < limit[3])
                goto done1;
            if (r[3] > limit[3])
                goto sub1;

            if (r[2] < limit[2])
                goto done1;
            if (r[2] > limit[2])
                goto sub1;

            if (r[1] < limit[1])
                goto done1;
            if (r[1] > limit[1])
                goto sub1;

            if (r[0] < limit[0])
                goto done1;

        sub1:
            sub_ddddmmmmssss(r[3],r[2],r[1],r[0], r[3],r[2],r[1],r[0], limit[3],limit[2],limit[1],limit[0]);
            goto check1;

        done1:

            ulong tbits = (i+j)*bits;
            ulong toff = tbits/FLINT_BITS;
            ulong tshift = tbits%FLINT_BITS;
            assert(zn > zcoeff_len + toff);

            if (tshift == 0)
            {
                add_ssssaaaaaaaa(z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                 z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                 r[3],r[2],r[1],r[0]);
            }
            else
            {
                r[4] =                       r[3] >> (64-tshift);
                r[3] = (r[3] << (tshift)) | (r[2] >> (64-tshift));
                r[2] = (r[2] << (tshift)) | (r[1] >> (64-tshift));
                r[1] = (r[1] << (tshift)) | (r[0] >> (64-tshift));
                r[0] =  r[0] << (tshift);

                add_sssssaaaaaaaaaa(z[toff+4],z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                    z[toff+4],z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                    r[4], r[3],r[2],r[1],r[0]);
            }
        }
    }

    for (; i < zlen; i++)
    {
        for (ulong l = 0; l < np; l++)
        {
            double x = Q.ffts[l].get_index(i);
            x = reduce_to_pm1n(x, Q.ffts[l].p, Q.ffts[l].pinv);
            if (x < 0)
                x += Q.ffts[l].p;
            ulong y = x;

            ulong* mult = Q.crts[np - 1].co_prime(l);
            if (l == 0)
            {
                umul_ppmm(r[1], r[0], mult[0], y);
                umul_ppmm(r[3], r[2], mult[2], y);
                umul_ppmm(u, v, mult[1], y);
                add_sssaaaaaa(r[3],r[2],r[1], r[3],r[2],r[1], 0,u,v);
            }
            else
            {
                umul_ppmm(t[1], t[0], mult[0], y);
                umul_ppmm(t[3], t[2], mult[2], y);
                umul_ppmm(u, v, mult[1], y);
                add_sssaaaaaa(t[3],t[2],t[1], t[3],t[2],t[1], 0,u,v);
                add_ssssaaaaaaaa(r[3],r[2],r[1],r[0], r[3],r[2],r[1],r[0], t[3],t[2],t[1],t[0]);
            }
        }

        ulong* limit = Q.crts[np - 1].prod_primes();

    check:
        if (r[3] < limit[3])
            goto done;
        if (r[3] > limit[3])
            goto sub;

        if (r[2] < limit[2])
            goto done;
        if (r[2] > limit[2])
            goto sub;

        if (r[1] < limit[1])
            goto done;
        if (r[1] > limit[1])
            goto sub;

        if (r[0] < limit[0])
            goto done;

    sub:
        sub_ddddmmmmssss(r[3],r[2],r[1],r[0], r[3],r[2],r[1],r[0], limit[3],limit[2],limit[1],limit[0]);
        goto check;

    done:

        ulong tbits = i*bits;
        ulong toff = tbits/FLINT_BITS;
        ulong tshift = tbits%FLINT_BITS;

        if (toff >= zn)
            break;

        if (tshift == 0)
        {
            if (zn - toff >= zcoeff_len)
            {
                add_ssssaaaaaaaa(z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                 z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                 r[3],r[2],r[1],r[0]);
            }
            else if (zn - toff == 3)
            {
                add_sssaaaaaa(z[toff+2],z[toff+1],z[toff+0],
                              z[toff+2],z[toff+1],z[toff+0],
                              r[2],r[1],r[0]);
            }
            else if (zn - toff == 2)
            {
                add_ssaaaa(z[toff+1],z[toff+0],
                           z[toff+1],z[toff+0],
                           r[1],r[0]);
            }
            else
            {
                z[toff+0] += r[0];
            }
        }
        else
        {
            r[4] =                       r[3] >> (64-tshift);
            r[3] = (r[3] << (tshift)) | (r[2] >> (64-tshift));
            r[2] = (r[2] << (tshift)) | (r[1] >> (64-tshift));
            r[1] = (r[1] << (tshift)) | (r[0] >> (64-tshift));
            r[0] =  r[0] << (tshift);

            if (zn - toff > zcoeff_len)
            {
                add_sssssaaaaaaaaaa(z[toff+4],z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                    z[toff+4],z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                    r[4], r[3],r[2],r[1],r[0]);
            }
            else if (zn - toff == zcoeff_len)
            {
                add_ssssaaaaaaaa(z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                 z[toff+3],z[toff+2],z[toff+1],z[toff+0],
                                 r[3],r[2],r[1],r[0]);
            }
            else if (zn - toff == 3)
            {
                add_sssaaaaaa(z[toff+2],z[toff+1],z[toff+0],
                              z[toff+2],z[toff+1],z[toff+0],
                              r[2],r[1],r[0]);
            }
            else if (zn - toff == 2)
            {
                add_ssaaaa(z[toff+1],z[toff+0],
                           z[toff+1],z[toff+0],
                           r[1],r[0]);
            }
            else
            {
                z[toff+0] += r[0];
            }
        }
    }
}


void mpn_mul_v2p1(
    mpn_ctx_v2& Q,
    ulong* z,
    const ulong* a, ulong an,
    const ulong* b, ulong bn)
{
    ulong zn = an + bn;
    ulong bits = 88;
    ulong alen = cld(FLINT_BITS * an, bits);
    ulong blen = cld(FLINT_BITS * bn, bits);
    alen = FLINT_MAX(1, alen);
    alen = FLINT_MAX(1, alen);
    ulong zlen = alen + blen - 1;
    ulong depth = clog2(zlen);
    depth = FLINT_MAX(depth, LG_BLK_SZ);
    ulong np = 4;

    ulong atrunc = round_up(alen, BLK_SZ);
    ulong btrunc = round_up(blen, BLK_SZ);
    ulong ztrunc = round_up(zlen, BLK_SZ);
//timeit_t timer;

//std::cout << "   zlen: " << zlen << std::endl;
//std::cout << "2^depth: " << (UWORD(1)<<depth) << std::endl;

    while (Q.nprimes() < np)
        Q.add_prime();

    for (ulong l = 0; l < np; l++)
        Q.ffts[l].set_depth(depth);

//std::cout << "bits <= " << 0.5*log2((Q.crts[np - 1].prod_primes_d()/blen)) << std::endl;

    ulong fft_data_size = Q.ffts[0].data_size();
    double* abuf = Q.fit_double_buffer(2*np*fft_data_size);
    double* bbuf = abuf + np*fft_data_size;

    ulong ttlen = bits;
    ulong nvs = (np + VEC_SZ - 1)/VEC_SZ;
    packed<double, VEC_SZ>* two_pow = new packed<double, VEC_SZ>[ttlen*nvs];
    fill_two_pow(two_pow, bits, Q.ffts, np);

//timeit_start(timer);
	for (ulong l = 0; l < np; l++)
        Q.ffts[l].set_data(abuf + l*fft_data_size);
    mpn_to_ffts<4, 88>(Q, a, an, atrunc, two_pow);

	for (ulong l = 0; l < np; l++)
        Q.ffts[l].set_data(bbuf + l*fft_data_size);
    mpn_to_ffts<4, 88>(Q, b, bn, btrunc, two_pow);
//timeit_stop(timer);
//if (timer->wall > 5)
//std::cout << "mod: " << timer->wall << std::endl;

    delete[] two_pow;

//timeit_start(timer);
	for (ulong l = 0; l < np; l++)
        Q.ffts[l].fft_trunc(btrunc, ztrunc);

	for (ulong l = 0; l < np; l++)
    {
        Q.ffts[l].set_data(abuf + l*fft_data_size);
        Q.ffts[l].fft_trunc(atrunc, ztrunc);
        // bake in 2^-depth * (p2*p3*p4)^-1 mod p1
        ulong t1, thi, tlo;
        ulong cop = Q.crts[np - 1].co_prime_red(l);
        thi = cop >> (FLINT_BITS - depth);
        tlo = cop << (depth);
        NMOD_RED2(t1, thi, tlo, Q.ffts[l].mod);
        t1 = nmod_inv(t1, Q.ffts[l].mod);
        Q.ffts[l].point_mul(bbuf + l*fft_data_size, t1);
        Q.ffts[l].ifft_trunc(ztrunc);
    }
//timeit_stop(timer);
//if (timer->wall > 5)
//std::cout << "fft: " << timer->wall << std::endl;

//timeit_start(timer);
    mpn_from_ffts<4, 4>(z, zn, zlen, Q, bits);
//timeit_stop(timer);
//if (timer->wall > 5)
//std::cout << "crt: " << timer->wall << std::endl;

//std::cout << "mpn_mul_v2p1 returning" << std::endl;
}


void mpn_mul_v2(
    ulong* z,
    const ulong* a, ulong an,
    const ulong* b, ulong bn)
{
    ulong zn = an + bn;
    ulong bits = 64;
    ulong alen = cld(FLINT_BITS * an, bits);
    ulong blen = cld(FLINT_BITS * bn, bits);
    alen = FLINT_MAX(1, alen);
    alen = FLINT_MAX(1, alen);
    ulong zlen = alen + blen - 1;
    ulong depth = clog2(zlen);
    depth = FLINT_MAX(depth, LG_BLK_SZ);

    // crt data
    int nprimes = 4;
    // [nprimes]
    ulong primes[4] = {UWORD(0x0003f00000000001),
                       UWORD(0x0002580000000001),
                       UWORD(0x00027c0000000001),
                       UWORD(0x00033c0000000001)};

    // supposed to hold 2*(product of primes)*(number of primes)
    ulong zcoeff_len = 4;

    // [nprimes][zcoeff_len]
    ulong co_primes[4*4] = {
        UWORD(0x0008100000000001), UWORD(0x8000001570500000), UWORD(0x000000000012d53d), UWORD(0x0000000000000000),
        UWORD(0x0009a80000000001), UWORD(0x0000001e8d900000), UWORD(0x00000000001fa3af), UWORD(0x0000000000000000),
        UWORD(0x0009840000000001), UWORD(0x0000001d8b600000), UWORD(0x00000000001dd936), UWORD(0x0000000000000000),
        UWORD(0x0008c40000000001), UWORD(0x00000018d5600000), UWORD(0x000000000016ed56), UWORD(0x0000000000000000),
    };

    // [nprimes]
    ulong co_primes_red[4] = {
        UWORD(0x000382d397829cbd),
        UWORD(0x000253dcc63f1413),
        UWORD(0x000001915c80aefb),
        UWORD(0x00025cb5f55a7ecb),
    };

    // [zcoeff_len]
    ulong prod_primes[4] = {
        UWORD(0x000c000000000001), UWORD(0x800000352f500000), UWORD(0x27a2280000673f78), UWORD(0x000000000000004a),
    };

    // [zcoeff_len + 1]
    ulong tt[5];
/*
std::cout.precision(17);
std::cout << "depth: " << depth << std::endl;
std::cout << "an: " << an << std::endl;
std::cout << "bn: " << bn << std::endl;
std::cout << "zn: " << zn << std::endl;
std::cout << "alen: " << alen << std::endl;
std::cout << "blen: " << blen << std::endl;
std::cout << "zlen: " << zlen << std::endl;
*/
    fftv2_ctx ctx(primes[0]);
    ctx.set_depth(depth);

    double* abuf = new double[ctx.data_size()];
    double* bbuf = new double[ctx.data_size()];
    ulong* zbuf = new ulong[zlen*zcoeff_len];

    for (int i = 0; i < nprimes; i++)
    {
//std::cout << "doing prime: " << ctx.p << std::endl;

        ctx.set_data(bbuf);
        ctx.from_mpn(b, bn, bits);
        ctx.fft();
        ctx.release_data();

        ctx.set_data(abuf);
        ctx.from_mpn(a, an, bits);
        ctx.fft();

        // bake in 2^-depth * (p2*p3*p4)^-1 mod p1
        {
            ulong t1, thi, tlo;
            thi = co_primes_red[i] >> (FLINT_BITS - depth);
            tlo = co_primes_red[i] << (depth);
            NMOD_RED2(t1, thi, tlo, ctx.mod);
//            t1 = ulong(1) << depth;
            ctx.point_mul(bbuf, nmod_inv(t1, ctx.mod));
        }

        ctx.ifft();

        for (ulong j = 0; j < zlen; j++)
        {
            double x = ctx.get_index(j);
            slong xx = reduce_to_pm1n(x, ctx.p, ctx.pinv);
            ulong y;
            if (xx < 0)
                y = ctx.mod.n + xx;
            else
                y = xx;

            if (i == 0)
                mpn_mul_1(zbuf + j*zcoeff_len, co_primes + i*zcoeff_len, zcoeff_len, y);
            else
                mpn_addmul_1(zbuf + j*zcoeff_len, co_primes + i*zcoeff_len, zcoeff_len, y);
        }

        ctx.release_data();

        if (i + 1 < nprimes)
        {
            ctx.set_prime(primes[i + 1]);
            ctx.set_depth(depth);
        }
    }

    mpn_zero(z, zn);

    for (ulong i = 0; i < zlen; i++)
    {
        ulong tbits = i*bits;
        ulong* t = zbuf + i*zcoeff_len;
        while (mpn_cmp(t, prod_primes, zcoeff_len) >= 0)
            mpn_sub_n(t, t, prod_primes, zcoeff_len);

        ulong toff = tbits/FLINT_BITS;
        ulong tshift = tbits%FLINT_BITS;

        if (toff >= zn)
            break;
/*
std::cout << "i: " << i << ", toff: " << toff << ", tshift: " << tshift << ", adding length: " << FLINT_MIN(zcoeff_len, zn - toff) << "  ";
std::cout << "t: " << format_hex(t[0]) << ", " <<
                      format_hex(t[1]) << ", " <<
                      format_hex(t[2]) << ", " <<
                      format_hex(t[3]) << std::endl;
*/

        if (tshift == 0)
        {
            mpn_add_n(z + toff, z + toff, t, FLINT_MIN(zcoeff_len, zn - toff));
//std::cout << "z[" << toff << "]: " << format_hex(z[toff]) << std::endl;
        }
        else
        {
            tt[zcoeff_len] = mpn_lshift(tt, t, zcoeff_len, tshift);
            mpn_add_n(z + toff, z + toff, tt, FLINT_MIN(zcoeff_len + 1, zn - toff));
        }
    }

    delete[] abuf;
    delete[] bbuf;
    delete[] zbuf;
}
