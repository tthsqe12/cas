#define BLK_SZ 256
#define LG_BLK_SZ 8
#define BLK_SHIFT 10

class fftv2_ctx {
public:
    std::vector<double> data;
    std::vector<std::vector<double>> wtab;
    const double * w2s;
    double p;
    double pinv;
    nmod_t mod;
    ulong primitive_root;
    ulong depth;
    ulong blk_sz;

    fftv2_ctx(ulong pp)
    {
        blk_sz = BLK_SZ;
        w2s = nullptr;
        p = pp;
        pinv = 1.0/p;
        nmod_init(&mod, pp);
        primitive_root = n_primitive_root_prime(pp);
        wtab.clear();
    }

    void fit_wtab(ulong k)
    {
        while (wtab.size() <= k)
        {
            ulong l = wtab.size();
            ulong n = UWORD(1) << l;
            // wtab[l] should be filled with n = 2^l th roots of unity
            ulong tabsize = (l == 0) ? n : n/2;
            std::vector<double> tab(tabsize);

            ulong w = nmod_pow_ui(primitive_root, (mod.n - 1)>>l, mod);
            ulong wi = 1;
            double half_p = mul(p, 0.5);

            for (ulong j = 0; j < tabsize; j++)
            {
                double x = wi;
                tab[n_revbin(j, l)/2] = (x >= half_p) ? x - p : x;
                wi = nmod_mul(wi, w, mod);
            }

            wtab.emplace_back(tab);
        }
    }

    inline ulong offset(ulong I)
    {
        return (I << LG_BLK_SZ);// + (I >> BLK_SHIFT);
    }

    void set_depth(ulong l)
    {
        depth = l;
        if (l < LG_BLK_SZ)
        {
            std::cout << "depth too large" << std::endl;
            abort();
        }
        data.resize(offset(ulong(1) << (l - LG_BLK_SZ)));
        fit_wtab(l);
        w2s = wtab[l].data();
    }

    inline double* from_index(ulong I)
    {
        return data.data() + offset(I);
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

    void transform_forward(ulong k, ulong j, ulong I, ulong S);
    void transform_forward_column(ulong k, ulong j, ulong I, ulong S);
    void transform_forward() {transform_forward(depth - LG_BLK_SZ, 0, 0, 1);}

    template <int> inline void transform_forward_basecase(double *X, ulong j);
};

//std::cout << "--- eval poly ---" << std::endl;
//std::cout << x0 << ", " << x1 << ", " << x2 << ", " << x3 << std::endl;
//std::cout << ws[4*j+0] << ", " << ws[4*j+1] << ", " << ws[4*j+2] << ", " << ws[4*j+3] << std::endl;


#define RADIX_4_BUTTERFLY(X0, X1, X2, X3) \
{ \
    double _x0 = X0; \
    double _x1 = X1; \
    double _x2 = X2; \
    double _x3 = X3; \
    double _y0 = add(_x0, mulmod2(_x2, w2, p, pinv)); \
    double _y1 = add(_x1, mulmod2(_x3, w2, p, pinv)); \
    double _y2 = sub(_x0, mulmod2(_x2, w2, p, pinv)); \
    double _y3 = sub(_x1, mulmod2(_x3, w2, p, pinv)); \
    _y0 = reduce_pm2n_to_pm1n(_y0, p); \
    _y2 = reduce_pm2n_to_pm1n(_y2, p); \
    _y1 = mulmod2(_y1, w, p, pinv); \
    _y3 = mulmod2(_y3, iw, p, pinv); \
    X0 = reduce_pm2n_to_pm1n(add(_y0, _y1), p); \
    X1 = reduce_pm2n_to_pm1n(sub(_y0, _y1), p); \
    X2 = reduce_pm2n_to_pm1n(add(_y2, _y3), p); \
    X3 = reduce_pm2n_to_pm1n(sub(_y2, _y3), p); \
}

#define PD_RADIX_4_BUTTERFLY(X0, X1, X2, X3) \
{ \
    PD<4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    _x0.load(X0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
    _x1.load(X1); \
    _x2.load(X2); \
    _x3.load(X3); \
    _y0 = add(_x0, mulmod2(_x2, w2, n, ninv)); \
    _y1 = add(_x1, mulmod2(_x3, w2, n, ninv)); \
    _y2 = sub(_x0, mulmod2(_x2, w2, n, ninv)); \
    _y3 = sub(_x1, mulmod2(_x3, w2, n, ninv)); \
    _y1 = mulmod2(_y1, w, n, ninv); \
    _y3 = mulmod2(_y3, iw, n, ninv); \
    _x0 = add(_y0, _y1); \
    _x1 = sub(_y0, _y1); \
    _x2 = add(_y2, _y3); \
    _x3 = sub(_y2, _y3); \
    _x0.store(X0); \
    _x1.store(X1); \
    _x2.store(X2); \
    _x3.store(X3); \
}

#define PD_RADIX_4_BUTTERFLY_2X(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
{ \
    PD<4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    PD<4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
        _z0 = reduce_to_pm1n(_z0, n, ninv); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _x2.load(X2); \
        _z2.load(Z2); \
    _x3.load(X3); \
        _z3.load(Z3); \
    _y0 = add(_x0, mulmod2(_x2, w2, n, ninv)); \
        _w0 = add(_z0, mulmod2(_z2, w2, n, ninv)); \
    _y1 = add(_x1, mulmod2(_x3, w2, n, ninv)); \
        _w1 = add(_z1, mulmod2(_z3, w2, n, ninv)); \
    _y2 = sub(_x0, mulmod2(_x2, w2, n, ninv)); \
        _w2 = sub(_z0, mulmod2(_z2, w2, n, ninv)); \
    _y3 = sub(_x1, mulmod2(_x3, w2, n, ninv)); \
        _w3 = sub(_z1, mulmod2(_z3, w2, n, ninv)); \
    _y1 = mulmod2(_y1, w, n, ninv); \
        _w1 = mulmod2(_w1, w, n, ninv); \
    _y3 = mulmod2(_y3, iw, n, ninv); \
        _w3 = mulmod2(_w3, iw, n, ninv); \
    _x0 = add(_y0, _y1); \
        _z0 = add(_w0, _w1); \
    _x1 = sub(_y0, _y1); \
        _z1 = sub(_w0, _w1); \
    _x2 = add(_y2, _y3); \
        _z2 = add(_w2, _w3); \
    _x3 = sub(_y2, _y3); \
        _z3 = sub(_w2, _w3); \
    _x0.store(X0); \
        _z0.store(Z0); \
    _x1.store(X1); \
        _z1.store(Z1); \
    _x2.store(X2); \
        _z2.store(Z2); \
    _x3.store(X3); \
        _z3.store(Z3); \
}


void fftv2_ctx::transform_forward_column(
   ulong k, // transform length 2^k
   ulong j,
   ulong I, // starting index
   ulong S) // stride
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        // column ffts
        for (ulong a = 0; a < ulong(1)<<k2; a++)
            transform_forward_column(k1, j, I + a*S, S<<k2);

        // row ffts
        for (ulong b = 0; b < ulong(1)<<k1; b++)
            transform_forward_column(k2, (j<<k1) + b, I + (b<<k2)*S, S);  

        return; 
    }

    if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
#if LG_BLK_SZ > 4
        PD<4> w = w2s[2*j];
        PD<4> w2 = w2s[j];
        PD<4> iw = w2s[2*j+1];
        PD<4> n(p);
        PD<4> ninv(pinv);
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_4_BUTTERFLY_2X(X0+i,   X1+i,   X2+i,   X3+i,
                                    X0+i+4, X1+i+4, X2+i+4, X3+i+4);
#else
        double w = w2s[2*j];
        double w2 = w2s[1*j];
        double iw = w2s[2*j+1];
        for (ulong i = 0; i < blk_sz; i += 4)
            RADIX_4_BUTTERFLY(X0[i], X1[i], X2[i], X3[i]);
#endif
    }
    else if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);

#if LG_BLK_SZ > 4
        PD<4> w(w2s[j]);
        PD<4> n(p);
        PD<4> ninv(pinv);
        for (ulong i = 0; i < blk_sz; i += 8)
        {
            PD<4> x0, x1, y0, y1, z0, z1, w0, w1;
            x0.load(X0+i);
                z0.load(X0+i+4);
            x1.load(X1+i);
                z1.load(X1+i+4);
            x1 = mulmod2(x1, w, n, ninv);
                z1 = mulmod2(z1, w, n, ninv);
            y0 = reduce_pm2n_to_pm1n(add(x0, x1), n);
                w0 = reduce_pm2n_to_pm1n(add(z0, z1), n);
            y1 = reduce_pm2n_to_pm1n(sub(x0, x1), n);
                w1 = reduce_pm2n_to_pm1n(sub(z0, z1), n);
            y0.store(X0+i);
                w0.store(X0+i+4);
            y1.store(X1+i);
                w1.store(X1+i+4);
        }
#else
        double w = w2s[j];
        for (ulong i = 0; i < blk_sz; i++)
        {
            double x0 = X0[i];
            double x1 = X1[i];
            x1 = mulmod2(x1, w, p, pinv);
            X0[i] = reduce_pm2n_to_pm1n(add(x0, x1), p);
            X1[i] = reduce_pm2n_to_pm1n(sub(x0, x1), p);
        }
#endif
    }
}

template <>
inline void fftv2_ctx::transform_forward_basecase<0>(double* X, ulong j)
{
}

template <>
inline void fftv2_ctx::transform_forward_basecase<1>(double* X, ulong j)
{
    double x0 = X[0];
    double x1 = X[1];
    x1 = mulmod2(x1, w2s[1*j], p, pinv);
    X[0] = reduce_pm2n_to_pm1n(add(x0, x1), p);
    X[1] = reduce_pm2n_to_pm1n(sub(x0, x1), p);
}

template <>
inline void fftv2_ctx::transform_forward_basecase<2>(double* X, ulong j)
{
    double w = w2s[2*j];
    double w2 = w2s[1*j];
    double iw = w2s[2*j+1];
    RADIX_4_BUTTERFLY(X[0], X[1], X[2], X[3]);
}

template <>
inline void fftv2_ctx::transform_forward_basecase<4>(double* X, ulong j)
{
    PD<4> w = w2s[2*j];
    PD<4> w2 = w2s[1*j];
    PD<4> iw = w2s[2*j+1];
    PD<4> n(p);
    PD<4> ninv(pinv);
    PD_RADIX_4_BUTTERFLY(X + 4*0, X + 4*1, X + 4*2, X + 4*3);
    transform_forward_basecase<2>(X + 4*0, 4*j+0);
    transform_forward_basecase<2>(X + 4*1, 4*j+1);
    transform_forward_basecase<2>(X + 4*2, 4*j+2);
    transform_forward_basecase<2>(X + 4*3, 4*j+3);
}

template <>
inline void fftv2_ctx::transform_forward_basecase<6>(double* X, ulong j)
{
#if LG_BLK_SZ > 4
    PD<4> w = w2s[2*j];
    PD<4> w2 = w2s[1*j];
    PD<4> iw = w2s[2*j+1];
    PD<4> n(p);
    PD<4> ninv(pinv);
    for (ulong i = 0; i < 16; i += 8)
        PD_RADIX_4_BUTTERFLY_2X(X+i+16*0, X+i+16*1, X+i+16*2, X+i+16*3,
                                X+i+16*0+4, X+i+16*1+4, X+i+16*2+4, X+i+16*3+4);
#else
    double w = w2s[2*j];
    double w2 = w2s[1*j];
    double iw = w2s[2*j+1];
    for (ulong i = 0; i < 16; i += 1)
        RADIX_4_BUTTERFLY(X+i+16*0, X+i+16*1, X+i+16*2, X+i+16*3);
#endif
    transform_forward_basecase<4>(X + 16*0, 4*j+0);
    transform_forward_basecase<4>(X + 16*1, 4*j+1);
    transform_forward_basecase<4>(X + 16*2, 4*j+2);
    transform_forward_basecase<4>(X + 16*3, 4*j+3);
}

template <>
void fftv2_ctx::transform_forward_basecase<8>(double* X, ulong j)
{
    PD<4> w = w2s[2*j];
    PD<4> w2 = w2s[1*j];
    PD<4> iw = w2s[2*j+1];
    PD<4> n(p);
    PD<4> ninv(pinv);
    for (ulong i = 0; i < 64; i += 8)
    {
        PD_RADIX_4_BUTTERFLY_2X(X+i+64*0, X+i+64*1, X+i+64*2, X+i+64*3,
                                X+i+64*0+4, X+i+64*1+4, X+i+64*2+4, X+i+64*3+4);
    }
    transform_forward_basecase<6>(X + 64*0, 4*j+0);
    transform_forward_basecase<6>(X + 64*1, 4*j+1);
    transform_forward_basecase<6>(X + 64*2, 4*j+2);
    transform_forward_basecase<6>(X + 64*3, 4*j+3);
}


void fftv2_ctx::transform_forward(
    ulong k, // transform length 2^(k+LG_BLK_SZ)
    ulong j,
    ulong I, // starting index
    ulong S) // stride
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        for (ulong a = 0; a < ulong(1)<<k2; a++)
            transform_forward_column(k1, j, I + a*S, S<<k2);

        for (ulong b = 0; b < ulong(1)<<k1; b++)
            transform_forward(k2, (j<<k1) + b, I + (b<<k2)*S, S);

        return;
    }

    if (k == 2)
    {
        // k1 = 2; k2 = 0
        transform_forward_column(2, j, I, S);
        transform_forward(0, 4*j+0, I+S*0, S);
        transform_forward(0, 4*j+1, I+S*1, S);
        transform_forward(0, 4*j+2, I+S*2, S);
        transform_forward(0, 4*j+3, I+S*3, S);
    }
    else if (k == 1)
    {
        // k1 = 1; k2 = 0
        transform_forward_column(1, j, I, S);
        transform_forward(0, 2*j+0, I+S*0, S);
        transform_forward(0, 2*j+1, I+S*1, S);
    }
    else
    {
        transform_forward_basecase<LG_BLK_SZ>(from_index(I), j);
    }
}


void test_v2(ulong L)
{
    ulong Xn = ulong(1)<<L;
    flint_rand_t state;
    timeit_t timer;
    ulong new_time, old_time;

    flint_randinit(state);

    if (L < LG_BLK_SZ)
        return;

    {
        double* X = new double[ulong(1)<<L];
        fftv2_ctx ctx(0x03f00000000001ULL);

    std::cout << "------- L: " << L << " --------" << std::endl;

        ctx.set_depth(L);

    //std::cout.precision(17);
    //std::cout << "ctx.p: " << ctx.p << std::endl;

        for (ulong i = 0; i < Xn; i++)
        {
            X[i] = i + 1;
            ctx.set_index(i, X[i]);
        }

    timeit_start(timer);
        ctx.transform_forward();
    timeit_stop(timer);
        //ctx.transform_forward_column(L, 0, 0, 1);

    new_time = FLINT_MAX(1, timer->wall);
    std::cout << "new time: " << new_time << " ms" << std::endl;

        // check answer
        for (int ii = 0; ii < 2 + 1000/L; ii++)
        {
            ulong i = n_randint(state, Xn);
            double y = ctx.eval_poly(X, Xn, (i&1) ? -ctx.w2s[i/2] : ctx.w2s[i/2]);
            if (reduce_pm1n_to_0n(y, ctx.p) !=
                reduce_pm1n_to_0n(reduce_to_pm1n(ctx.get_index(i), ctx.p, ctx.pinv), ctx.p))
            {
                std::cout << "ctx.get_index(" << i << "): " << ctx.get_index(i) << std::endl;
                std::cout << "y: " << y << std::endl;

                std::cout << "oops " << reduce_pm1n_to_0n(ctx.get_index(i), ctx.p) << " != " << reduce_pm1n_to_0n(y, ctx.p) << std::endl;
                std::abort();
            }
        }

        delete[] X;
    }

    {
        PD_fft_ctx<4> ctx;
        PD<4>* y = new PD<4>[Xn];

        ulong primes[4] = {0x0003f00000000001ULL,
                           0x0003dc0000000001ULL,
                           0x0003cf0000000001ULL,
                           0x0003a50000000001ULL};
        ctx.set_mod(primes);

        for (ulong i = 0; i < Xn; i++)
            y[i] = PD<4>(double(i + 1));

        ctx.fit_wtab(L);

    timeit_start(timer);
        PD_fft_short<4>(y, Xn, Xn, L, ctx);
    timeit_stop(timer);

    old_time = FLINT_MAX(1, timer->wall/4);
    std::cout << "old time: " << old_time << " ms" << std::endl;

        delete[] y;
    }

    std::cout << "ratio: " << double(new_time)/double(old_time) << std::endl;

    flint_randclear(state);
}
