#define BLK_SZ 256
#define LG_BLK_SZ 8
#define BLK_SHIFT 10

struct fftv2_ctx {
    std::vector<std::vector<double>> wtab;
    double * data;
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
        data = nullptr;
        p = pp;
        pinv = 1.0/p;
        nmod_init(&mod, pp);
        primitive_root = n_primitive_root_prime(pp);
        wtab.clear();
    }

    void set_prime(ulong pp)
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

//std::cout << "filling depth " << l << " for prime " << format_hex(mod.n) << std::endl;

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
        return (I << LG_BLK_SZ) + (I >> BLK_SHIFT);
    }

    ulong data_size()
    {
        return offset(ulong(1) << (depth - LG_BLK_SZ));
    }

    void set_data(double* d)
    {
        assert(data == nullptr);
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
        assert(data == nullptr);
        depth = l;
        if (l < LG_BLK_SZ)
        {
            std::cout << "depth " << l << " too small" << std::endl;
            abort();
        }
        fit_wtab(l);
        w2s = wtab[l].data();
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

    void transform_forward(ulong k, ulong j, ulong I, ulong S);
    void transform_forward_column(ulong k, ulong j, ulong I, ulong S);
    void transform_forward() {transform_forward(depth - LG_BLK_SZ, 0, 0, 1);}
    template <int> inline void transform_forward_basecase(double *X, ulong j);

    void transform_reverse(ulong k, ulong j, ulong I, ulong S);
    void transform_reverse_column(ulong k, ulong j, ulong I, ulong S);
    void transform_reverse() {transform_reverse(depth - LG_BLK_SZ, 0, 0, 1);}
    template <int> inline void transform_reverse_basecase(double *X, ulong j);


    void from_mpn(const ulong * a, ulong an, ulong bits);
    void point_mul(const double* b, ulong mm);
};

// pointwise mul of self with b and m
void fftv2_ctx::point_mul(const double* b, ulong mm)
{
    double m = mm;
    if (m > 0.5*p)
        m -= p;

    PD<4> M = m;
    PD<4> n = p;
    PD<4> ninv = pinv;
    for (ulong I = 0; I < ulong(1) << (depth - LG_BLK_SZ); I++)
    {
        double* x = from_index(I);
        const double* bx = b + offset(I);
        ulong j = 0;
        do {
            PD<4> x0, x1, x2, x3, b0, b1, b2, b3;
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

ulong saturate_bits(ulong a)
{
    a |= a >> 1;
    a |= a >> 2;
    a |= a >> 4;
    a |= a >> 8;
    a |= a >> 16;
    return a;
}

//std::cout << "--- eval poly ---" << std::endl;
//std::cout << x0 << ", " << x1 << ", " << x2 << ", " << x3 << std::endl;
//std::cout << ws[4*j+0] << ", " << ws[4*j+1] << ", " << ws[4*j+2] << ", " << ws[4*j+3] << std::endl;


/* radix 2 *******************************************************************/

#define RADIX_2_FORWARD_PARAM(j) \
    double w = w2s[j]; \
    double n = p; \
    double ninv = pinv;

#define RADIX_2_FORWARD(X0, X1) \
{ \
    double _x0, _x1, _y0, _y1; \
    _x0 = *(X0); \
    _x1 = *(X1); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
    _x1 = mulmod2(_x1, w, n, ninv); \
    _y0 = add(_x0, _x1); \
    _y1 = sub(_x0, _x1); \
    *(X0) = _y0; \
    *(X1) = _y1; \
}

#define PD_RADIX_2_FORWARD_PARAM(j) \
    PD<4> w = w2s[j]; \
    PD<4> n = p; \
    PD<4> ninv = pinv;

#define PD_RADIX_2_FORWARD_2X(X0, X1, Z0, Z1) \
{ \
    PD<4> _x0, _x1, _y0, _y1, _z0, _z1, _w0, _w1; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
        _z0 = reduce_to_pm1n(_z0, n, ninv); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _x1 = mulmod2(_x1, w, n, ninv); \
        _z1 = mulmod2(_z1, w, n, ninv); \
    _y0 = add(_x0, _x1); \
        _w0 = add(_z0, _z1); \
    _y1 = sub(_x0, _x1); \
        _w1 = sub(_z0, _z1); \
    _y0.store(X0); \
        _w0.store(Z0); \
    _y1.store(X1); \
        _w1.store(Z1); \
}


#define RADIX_2_REVERSE_PARAM(j) \
    ulong mask = saturate_bits(j); \
    double W  = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
    double n = p; \
    double ninv = pinv;

#define RADIX_2_REVERSE(X0, X1) \
{ \
    double _x0, _x1, _y0, _y1; \
    _x0 = *(X0); \
    _x1 = *(X1); \
    _y0 = add(_x0, _x1); \
    _y1 = sub(_x1, _x0); \
    _y0 = reduce_to_pm1n(_y0, n, ninv); \
    _y1 = mulmod2(_y1, W, n, ninv); \
    *(X0) = _y0; \
    *(X1) = _y1; \
}

#define PD_RADIX_2_REVERSE_PARAM(j) \
    ulong mask = saturate_bits(j); \
    PD<4> W  = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
    PD<4> n = p; \
    PD<4> ninv = pinv;

#define PD_RADIX_2_REVERSE_2X(X0, X1, Z0, Z1) \
{ \
    PD<4> _x0, _x1, _y0, _y1; \
    PD<4> _z0, _z1, _w0, _w1; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _y0 = add(_x0, _x1); \
        _w0 = add(_z0, _z1); \
    _y1 = sub(_x1, _x0); \
        _w1 = sub(_z1, _z0); \
    _y0 = reduce_to_pm1n(_y0, n, ninv); \
        _w0 = reduce_to_pm1n(_w0, n, ninv); \
    _y1 = mulmod2(_y1, W, n, ninv); \
        _w1 = mulmod2(_w1, W, n, ninv); \
    _y0.store(X0); \
        _w0.store(Z0); \
    _y1.store(X1); \
        _w1.store(Z1); \
}




/* radix 4 *******************************************************************/

/*
forward butterfly:
    b0 = a0 + w^2*a2 +   w*(a1 + w^2*a3)
    b1 = a0 + w^2*a2 -   w*(a1 + w^2*a3)
    b2 = a0 - w^2*a2 + i*w*(a1 - w^2*a3)
    b3 = a0 - w^2*a2 - i*w*(a1 - w^2*a3)
*/

#define RADIX_4_FORWARD_PARAM(j) \
    double w = w2s[2*j]; \
    double w2 = w2s[1*j]; \
    double iw = w2s[2*j+1]; \
    double n = p; \
    double ninv = pinv;

#define RADIX_4_FORWARD(X0, X1, X2, X3) \
{ \
    double _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    _x0 = *(X0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
    _x1 = *(X1); \
    _x2 = *(X2); \
    _x3 = *(X3); \
    _x2 = mulmod2(_x2, w2, n, ninv); \
    _x3 = mulmod2(_x3, w2, n, ninv); \
    _y0 = add(_x0, _x2); \
    _y1 = add(_x1, _x3); \
    _y2 = sub(_x0, _x2); \
    _y3 = sub(_x1, _x3); \
    _y1 = mulmod2(_y1, w, n, ninv); \
    _y3 = mulmod2(_y3, iw, n, ninv); \
    _x0 = add(_y0, _y1); \
    _x1 = sub(_y0, _y1); \
    _x2 = add(_y2, _y3); \
    _x3 = sub(_y2, _y3); \
    *(X0) = _x0; \
    *(X1) = _x1; \
    *(X2) = _x2; \
    *(X3) = _x3; \
}

#define PD_RADIX_4_FORWARD_PARAM(j) \
    PD<4> w = w2s[2*j]; \
    PD<4> w2 = w2s[1*j]; \
    PD<4> iw = w2s[2*j+1]; \
    PD<4> n(p); \
    PD<4> ninv(pinv);

#define PD_RADIX_4_FORWARD(X0, X1, X2, X3) \
{ \
    PD<4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    _x0.load(X0); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
    _x1.load(X1); \
    _x2.load(X2); \
    _x3.load(X3); \
    _x2 = mulmod2(_x2, w2, n, ninv); \
    _x3 = mulmod2(_x3, w2, n, ninv); \
    _y0 = add(_x0, _x2); \
    _y1 = add(_x1, _x3); \
    _y2 = sub(_x0, _x2); \
    _y3 = sub(_x1, _x3); \
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

#define PD_RADIX_4_FORWARD_2X(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
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
    _x2 = mulmod2(_x2, w2, n, ninv); \
        _z2 = mulmod2(_z2, w2, n, ninv); \
    _x3 = mulmod2(_x3, w2, n, ninv); \
        _z3 = mulmod2(_z3, w2, n, ninv); \
    _y0 = add(_x0, _x2); \
        _w0 = add(_z0, _z2); \
    _y1 = add(_x1, _x3); \
        _w1 = add(_z1, _z3); \
    _y2 = sub(_x0, _x2); \
        _w2 = sub(_z0, _z2); \
    _y3 = sub(_x1, _x3); \
        _w3 = sub(_z1, _z3); \
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


/*
inverse butterfly:
    4*a0 =            (b0 + b1) +        (b2 + b3)
    4*a1 =       w^-1*(b0 - b1) - i*w^-1*(b2 - b3)
    4*a2 = w^-2*(     (b0 + b1) -        (b2 + b3))
    4*a3 = w^-2*(w^-1*(b0 - b1) + i*w^-1*(b2 - b3))

    W  := -w^-1
    W2 := -w^-2
    IW := i*w^-1
*/

#define RADIX_4_REVERSE_PARAM(j) \
    ulong mask = saturate_bits(j); \
    double W  = ((j) == 0) ? -w2s[0] : w2s[(2*(j)  )^(mask)]; \
    double W2 = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
    double IW = ((j) == 0) ?  w2s[1] : w2s[(2*(j)+1)^(mask)]; \
    double n = p; \
    double ninv = pinv;

#define RADIX_4_REVERSE(X0, X1, X2, X3) \
{ \
    double _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    _x0 = *(X0); \
    _x1 = *(X1); \
    _x2 = *(X2); \
    _x3 = *(X3); \
    _y0 = add(_x0, _x1); \
    _y1 = add(_x2, _x3); \
    _y2 = sub(_x0, _x1); \
    _y3 = sub(_x3, _x2); \
    _y2 = mulmod2(_y2, W, n, ninv); \
    _y3 = mulmod2(_y3, IW, n, ninv); \
    _x0 = add(_y0, _y1); \
    _x1 = sub(_y3, _y2); \
    _x2 = sub(_y1, _y0); \
    _x3 = add(_y3, _y2); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
    _x2 = mulmod2(_x2, W2, n, ninv); \
    _x3 = mulmod2(_x3, W2, n, ninv); \
    *(X0) = _x0; \
    *(X1) = _x1; \
    *(X2) = _x2; \
    *(X3) = _x3; \
}

#define PD_RADIX_4_REVERSE_PARAM(j) \
    ulong mask = saturate_bits(j); \
    PD<4> W  = ((j) == 0) ? -w2s[0] : w2s[(2*(j)  )^(mask)]; \
    PD<4> W2 = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
    PD<4> IW = ((j) == 0) ?  w2s[1] : w2s[(2*(j)+1)^(mask)]; \
    PD<4> n(p); \
    PD<4> ninv(pinv); \

#define PD_RADIX_4_REVERSE(X0, X1, X2, X3) \
{ \
    PD<4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    _x0.load(X0); \
    _x1.load(X1); \
    _x2.load(X2); \
    _x3.load(X3); \
    _y0 = add(_x0, _x1); \
    _y1 = add(_x2, _x3); \
    _y2 = sub(_x0, _x1); \
    _y3 = sub(_x3, _x2); \
    _y2 = mulmod2(_y2, W, n, ninv); \
    _y3 = mulmod2(_y3, IW, n, ninv); \
    _x0 = add(_y0, _y1); \
    _x1 = sub(_y3, _y2); \
    _x2 = sub(_y1, _y0); \
    _x3 = add(_y3, _y2); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
    _x2 = mulmod2(_x2, W2, n, ninv); \
    _x3 = mulmod2(_x3, W2, n, ninv); \
    _x0.store(X0); \
    _x1.store(X1); \
    _x2.store(X2); \
    _x3.store(X3); \
}

#define PD_RADIX_4_REVERSE_2X(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
{ \
    PD<4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    PD<4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
    _x0.load(X0); \
        _z0.load(Z0); \
    _x1.load(X1); \
        _z1.load(Z1); \
    _x2.load(X2); \
        _z2.load(Z2); \
    _x3.load(X3); \
        _z3.load(Z3); \
    _y0 = add(_x0, _x1); \
        _w0 = add(_z0, _z1); \
    _y1 = add(_x2, _x3); \
        _w1 = add(_z2, _z3); \
    _y2 = sub(_x0, _x1); \
        _w2 = sub(_z0, _z1); \
    _y3 = sub(_x3, _x2); \
        _w3 = sub(_z3, _z2); \
    _y2 = mulmod2(_y2, W, n, ninv); \
        _w2 = mulmod2(_w2, W, n, ninv); \
    _y3 = mulmod2(_y3, IW, n, ninv); \
        _w3 = mulmod2(_w3, IW, n, ninv); \
    _x0 = add(_y0, _y1); \
        _z0 = add(_w0, _w1); \
    _x1 = sub(_y3, _y2); \
        _z1 = sub(_w3, _w2); \
    _x2 = sub(_y1, _y0); \
        _z2 = sub(_w1, _w0); \
    _x3 = add(_y3, _y2); \
        _z3 = add(_w3, _w2); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
        _z0 = reduce_to_pm1n(_z0, n, ninv); \
    _x2 = mulmod2(_x2, W2, n, ninv); \
        _z2 = mulmod2(_z2, W2, n, ninv); \
    _x3 = mulmod2(_x3, W2, n, ninv); \
        _z3 = mulmod2(_z3, W2, n, ninv); \
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
        PD_RADIX_4_FORWARD_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_4_FORWARD_2X(X0+i,   X1+i,   X2+i,   X3+i,
                                    X0+i+4, X1+i+4, X2+i+4, X3+i+4);
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

void fftv2_ctx::transform_reverse_column(
   ulong k, // transform length 2^k
   ulong j,
   ulong I, // starting index
   ulong S) // stride
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        // row ffts
        for (ulong b = 0; b < ulong(1)<<k1; b++)
            transform_reverse_column(k2, (j<<k1) + b, I + (b<<k2)*S, S);

        // column ffts
        for (ulong a = 0; a < ulong(1)<<k2; a++)
            transform_reverse_column(k1, j, I + a*S, S<<k2);

        return;
    }

    if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
#if LG_BLK_SZ > 4
        PD_RADIX_4_REVERSE_PARAM(j)
        for (ulong i = 0; i < blk_sz; i += 8)
            PD_RADIX_4_REVERSE_2X(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
#else
        RADIX_4_REVERSE_PARAM(j)
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

template <>
inline void fftv2_ctx::transform_forward_basecase<0>(double* X, ulong j)
{
}

template <>
inline void fftv2_ctx::transform_reverse_basecase<0>(double* X, ulong j)
{
}

template <>
inline void fftv2_ctx::transform_forward_basecase<1>(double* X, ulong j)
{
    RADIX_2_FORWARD_PARAM(j)
    RADIX_2_FORWARD(X+0, X+1);
}

template <>
inline void fftv2_ctx::transform_reverse_basecase<1>(double* X, ulong j)
{
    RADIX_2_REVERSE_PARAM(j)
    RADIX_2_REVERSE(X+0, X+1);
}

template <>
inline void fftv2_ctx::transform_forward_basecase<2>(double* X, ulong j)
{
    RADIX_4_FORWARD_PARAM(j)
    RADIX_4_FORWARD(X+0, X+1, X+2, X+3);
}

template <>
inline void fftv2_ctx::transform_reverse_basecase<2>(double* X, ulong j)
{
    RADIX_4_REVERSE_PARAM(j)
    RADIX_4_REVERSE(X+0, X+1, X+2, X+3);
}

template <>
inline void fftv2_ctx::transform_forward_basecase<4>(double* X, ulong j)
{
    PD_RADIX_4_FORWARD_PARAM(j)
    PD_RADIX_4_FORWARD(X+4*0, X+4*1, X+4*2, X+4*3);
    transform_forward_basecase<2>(X+4*0, 4*j+0);
    transform_forward_basecase<2>(X+4*1, 4*j+1);
    transform_forward_basecase<2>(X+4*2, 4*j+2);
    transform_forward_basecase<2>(X+4*3, 4*j+3);
}

template <>
inline void fftv2_ctx::transform_reverse_basecase<4>(double* X, ulong j)
{
    transform_reverse_basecase<2>(X+4*0, 4*j+0);
    transform_reverse_basecase<2>(X+4*1, 4*j+1);
    transform_reverse_basecase<2>(X+4*2, 4*j+2);
    transform_reverse_basecase<2>(X+4*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    PD_RADIX_4_REVERSE(X+4*0, X+4*1, X+4*2, X+4*3);
}

template <>
inline void fftv2_ctx::transform_forward_basecase<6>(double* X, ulong j)
{
    PD_RADIX_4_FORWARD_PARAM(j)
    for (ulong i = 0; i < 16; i += 8)
        PD_RADIX_4_FORWARD_2X(X+i+16*0, X+i+16*1, X+i+16*2, X+i+16*3, X+i+16*0+4, X+i+16*1+4, X+i+16*2+4, X+i+16*3+4);
    transform_forward_basecase<4>(X+16*0, 4*j+0);
    transform_forward_basecase<4>(X+16*1, 4*j+1);
    transform_forward_basecase<4>(X+16*2, 4*j+2);
    transform_forward_basecase<4>(X+16*3, 4*j+3);
}

template <>
inline void fftv2_ctx::transform_reverse_basecase<6>(double* X, ulong j)
{
    transform_reverse_basecase<4>(X+16*0, 4*j+0);
    transform_reverse_basecase<4>(X+16*1, 4*j+1);
    transform_reverse_basecase<4>(X+16*2, 4*j+2);
    transform_reverse_basecase<4>(X+16*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    for (ulong i = 0; i < 16; i += 8)
        PD_RADIX_4_REVERSE_2X(X+i+16*0, X+i+16*1, X+i+16*2, X+i+16*3, X+i+16*0+4, X+i+16*1+4, X+i+16*2+4, X+i+16*3+4);
}

template <>
void fftv2_ctx::transform_forward_basecase<8>(double* X, ulong j)
{
    PD_RADIX_4_FORWARD_PARAM(j)
    for (ulong i = 0; i < 64; i += 8)
        PD_RADIX_4_FORWARD_2X(X+i+64*0, X+i+64*1, X+i+64*2, X+i+64*3, X+i+64*0+4, X+i+64*1+4, X+i+64*2+4, X+i+64*3+4);
    transform_forward_basecase<6>(X+64*0, 4*j+0);
    transform_forward_basecase<6>(X+64*1, 4*j+1);
    transform_forward_basecase<6>(X+64*2, 4*j+2);
    transform_forward_basecase<6>(X+64*3, 4*j+3);
}

template <>
void fftv2_ctx::transform_reverse_basecase<8>(double* X, ulong j)
{
    transform_reverse_basecase<6>(X+64*0, 4*j+0);
    transform_reverse_basecase<6>(X+64*1, 4*j+1);
    transform_reverse_basecase<6>(X+64*2, 4*j+2);
    transform_reverse_basecase<6>(X+64*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    for (ulong i = 0; i < 64; i += 8)
        PD_RADIX_4_REVERSE_2X(X+i+64*0, X+i+64*1, X+i+64*2, X+i+64*3, X+i+64*0+4, X+i+64*1+4, X+i+64*2+4, X+i+64*3+4);
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

void fftv2_ctx::transform_reverse(
    ulong k, // transform length 2^(k+LG_BLK_SZ)
    ulong j,
    ulong I, // starting index
    ulong S) // stride
{
    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        for (ulong b = 0; b < ulong(1)<<k1; b++)
            transform_reverse(k2, (j<<k1) + b, I + (b<<k2)*S, S);

        for (ulong a = 0; a < ulong(1)<<k2; a++)
            transform_reverse_column(k1, j, I + a*S, S<<k2);

        return;
    }

    if (k == 2)
    {
        // k1 = 2; k2 = 0
        transform_reverse(0, 4*j+0, I+S*0, S);
        transform_reverse(0, 4*j+1, I+S*1, S);
        transform_reverse(0, 4*j+2, I+S*2, S);
        transform_reverse(0, 4*j+3, I+S*3, S);
        transform_reverse_column(2, j, I, S);
    }
    else if (k == 1)
    {
        // k1 = 1; k2 = 0
        transform_reverse(0, 2*j+0, I+S*0, S);
        transform_reverse(0, 2*j+1, I+S*1, S);
        transform_reverse_column(1, j, I, S);
    }
    else
    {
        transform_reverse_basecase<LG_BLK_SZ>(from_index(I), j);
    }
}

ulong cld(ulong a, ulong b)
{
    return (a + b - 1)/b;
}

ulong clog2(ulong x)
{
    if (x <= 2)
        return x == 2;
    return FLINT_BITS - __builtin_clzll(x - 1);
}

ulong next_fft_number(ulong p)
{
    ulong bits = FLINT_BIT_COUNT(p);
    ulong l; count_trailing_zeros(l, p - 1);
    ulong q = p - (UWORD(2) << l);
    if (bits < 20)
        std::abort();
    if (FLINT_BIT_COUNT(q) == bits)
        return q;
    if (l < 5)
        return (UWORD(1) << (bits - 2)) + 1;
    return (UWORD(1) << (bits)) - (UWORD(1) << (l - 1)) + 1;   
}

#if 0
void fftv2_ctx::from_mpn(
    const ulong * a,
    ulong an,
    ulong bits)
{
    ulong xn = UWORD(1) << depth;
    ulong ttlen = bits/FLINT_BITS + 1;
    ulong* tt = new ulong[ttlen];

    FLINT_ASSERT(bits >= FLINT_BITS);

    ulong mask = (UWORD(1) << (bits%FLINT_BITS)) - 1;

    for (ulong i = 0; i < xn; i++)
    {
        ulong abits = i*bits;
        ulong aoff = abits/FLINT_BITS;
        ulong ashift = abits%FLINT_BITS;

        if (ashift == 0)
        {
            for (ulong j = 0; j < ttlen; j++)
            {
                tt[j] = (aoff+j < an) ? a[aoff+j] : 0;
            }            
        }
        else
        {
            for (ulong j = 0; j < ttlen; j++)
            {
                ulong alo = (aoff+j+0 < an) ? a[aoff+j+0] : 0;
                ulong ahi = (aoff+j+1 < an) ? a[aoff+j+1] : 0;
                tt[j] = (alo >> ashift) | (ahi << (FLINT_BITS - ashift));
            }
        }

        tt[ttlen-1] &= mask;



        ulong r = mpn_mod_1(tt, ttlen, mod.n);
//std::cout << "x[" << i << "]: " << r << "    tt[0]: " << format_hex(tt[0]) << "  tt[1]: " << format_hex(tt[1]) << std::endl;
        set_index(i, (mod.n - r < r) ? -double(mod.n-r) : double(r));
    }
}
#endif

#if 0
void fftv2_ctx::from_mpn(
    const ulong * a,
    ulong an,
    ulong bits)
{
    ulong xn = UWORD(1) << depth;
    ulong ttlen = bits/FLINT_BITS + 1;
    ulong* tt = new ulong[ttlen];

    FLINT_ASSERT(bits >= FLINT_BITS);

    ulong mask = (UWORD(1) << (bits%FLINT_BITS)) - 1;

    ulong* two_pow = new ulong[ttlen*FLINT_BITS];

    two_pow[0] = 1;
    for (ulong i = 1; i < ttlen*FLINT_BITS; i++)
        two_pow[i] = nmod_add(two_pow[i-1], two_pow[i-1], mod);

#define aindex(i) (((i) < an) ? a[i] : UWORD(0))

#define MADD(hi, lo, a, b) \
{ \
    ulong _a = a; \
    ulong _b = b; \
    ulong _phi, _plo; \
    umul_ppmm(_phi, _plo, _a, _b); \
    add_ssaaaa(hi, lo, hi, lo, _phi, _plo); \
}

    for (ulong i = 0; i < xn; i++)
    {
        ulong abits = i*bits;
        ulong aoff = abits/FLINT_BITS;
        ulong ashift = abits%FLINT_BITS;
        ulong j, c, r, sum_lo, sum_hi;

        sum_lo = sum_hi = 0;
        if (ashift == 0)
        {
            for (j = 0; j < ttlen - 1; j++)
                MADD(sum_hi, sum_lo, aindex(aoff+j), two_pow[FLINT_BITS*j]);
            MADD(sum_hi, sum_lo, mask & aindex(aoff+j), two_pow[FLINT_BITS*(ttlen - 1)]);
        }
        else
        {
            for (j = 0; j < ttlen - 1; j++)
            {
                c = (aindex(aoff+j) >> ashift) | (aindex(aoff+j+1) << (FLINT_BITS - ashift));
                MADD(sum_hi, sum_lo, c, two_pow[FLINT_BITS*j]);
            }
            c = (aindex(aoff+j) >> ashift) | (aindex(aoff+j+1) << (FLINT_BITS - ashift));
            MADD(sum_hi, sum_lo, mask & c, two_pow[FLINT_BITS*j]);
        }

        r = sum_hi;
        sum_hi = 0;
        MADD(sum_hi, sum_lo, r, two_pow[FLINT_BITS]);
        NMOD_RED2(r, sum_hi, sum_lo, mod);

//std::cout << "x[" << i << "]: " << r /* << "    tt[0]: " << format_hex(tt[0]) << "  tt[1]: " << format_hex(tt[1])*/ << std::endl;
        set_index(i, (mod.n - r < r) ? -double(mod.n-r) : double(r));
    }

    delete two_pow;
}
#endif


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
    }

    for (; i < xn; i++)
    {
        set_index(i, 0.0);
    }

    delete two_pow;
}



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
};


struct mpn_ctx_v2 {
    std::vector<fftv2_ctx> ffts;
    std::vector<crt_data> crts;

    ulong nprimes() {
        return ffts.size();
    }

    mpn_ctx_v2(ulong p)
    {
        ffts.emplace_back(p);
        crts.emplace_back(p, 1, 1);

        crts[0].co_prime_red(0) = 1;
        crts[0].co_prime(0)[0] = 1;
        crts[0].prod_primes()[0] = p;
    }

    void add_prime(ulong p)
    {
        ffts.emplace_back(p);
        ulong len = crts.back().coeff_len;
        ulong* t = new ulong[len + 2];

        t[len + 1] = 0;
        t[len] = mpn_mul_1(t, crts.back().prod_primes(), len, p);
        len += (t[len] != 0);

        // leave room for one more bit
        if (FLINT_SIGN_EXT(t[len-1]) != 0)
            len += 1;

        crts.emplace_back(p, len, nprimes());

        // set product of primes
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
        ctx.transform_forward();
timeit_stop(timer);
std::cout << " b forward: " << timer->wall << std::endl;
        ctx.release_data();

        ctx.set_data(abuf);
timeit_start(timer);
        ctx.from_mpn(a, an, bits);
timeit_stop(timer);
std::cout << "a from mpn: " << timer->wall << std::endl;
timeit_start(timer);
        ctx.transform_forward();
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
        ctx.transform_reverse();
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


/*
std::cout << "i: " << i << ", toff: " << toff << ", tshift: " << tshift << ", adding length: " << FLINT_MIN(zcoeff_len, zn - toff) << "  ";
std::cout << "t: " << format_hex(t[0]) << ", " <<
                      format_hex(t[1]) << ", " <<
                      format_hex(t[2]) << ", " <<
                      format_hex(t[3]) << std::endl;
*/

void test_mpn_mul_v2(mpn_ctx_v2& Q, ulong an, ulong bn)
{
    ulong* a = new ulong[an];
    ulong* b = new ulong[bn];
    ulong* z1 = new ulong[an + bn];
    ulong* z2 = new ulong[an + bn];
    timeit_t timer;
    ulong new_time, old_time;

    for (ulong i = 0; i < an; i++)
        a[i] = -UWORD(3);

    for (ulong i = 0; i < bn; i++)
        b[i] = -UWORD(2);

timeit_start(timer);
    if (an < bn)
        mpn_mul(z1, b, bn, a, an);
    else
        mpn_mul(z1, a, an, b, bn);
timeit_stop(timer);
old_time = timer->wall;

timeit_start(timer);
    mpn_mul_v2(Q, z2, a, an, b, bn);
timeit_stop(timer);
new_time = timer->wall;

if (old_time > 30)
{
    std::cout << "----------- " << an << " * " << bn << " words -------------" << std::endl;
    std::cout << "old_time: " << old_time << std::endl;
    std::cout << "new_time: " << new_time << std::endl;
    std::cout << " old/new: " << double(old_time)/double(new_time) << std::endl;
}


    for (ulong i = 0; i < an + bn; i++)
    {
        if (z1[i] != z2[i])
        {
            std::cout << "mismatch!!!!" << std::endl;
            std::cout << "an: " << an << ", bn: " << bn << std::endl;
            std::cout << "z1[" << i << "]: " << format_hex(z1[i]) << std::endl;
            std::cout << "z2[" << i << "]: " << format_hex(z2[i]) << std::endl;
            abort();
        }
    }

    delete[] a;
    delete[] b;
    delete[] z1;
    delete[] z2;
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
        ctx.transform_forward();
        ctx.release_data();

        ctx.set_data(abuf);
        ctx.from_mpn(a, an, bits);
        ctx.transform_forward();

        // bake in 2^-depth * (p2*p3*p4)^-1 mod p1
        {
            ulong t1, thi, tlo;
            thi = co_primes_red[i] >> (FLINT_BITS - depth);
            tlo = co_primes_red[i] << (depth);
            NMOD_RED2(t1, thi, tlo, ctx.mod);
//            t1 = ulong(1) << depth;
            ctx.point_mul(bbuf, nmod_inv(t1, ctx.mod));
        }

        ctx.transform_reverse();

        for (ulong j = 0; j < zlen; j++)
        {
            double x = ctx.get_index(j);
            slong xx = reduce_to_pm1n(x, ctx.p, ctx.pinv);
            ulong y;
            if (xx < 0)
                y = ctx.mod.n + xx;
            else
                y = xx;

if (j < 0)
{
    std::cout << "y[" << j << "]: " << y << std::endl;
}

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

void test_mpn_mul_v2(ulong an, ulong bn)
{
    ulong* a = new ulong[an];
    ulong* b = new ulong[bn];
    ulong* z1 = new ulong[an + bn];
    ulong* z2 = new ulong[an + bn];
    timeit_t timer;
    ulong new_time, old_time;

    for (ulong i = 0; i < an; i++)
        a[i] = -UWORD(1);

    for (ulong i = 0; i < bn; i++)
        b[i] = -UWORD(1);

timeit_start(timer);
    if (an < bn)
        mpn_mul(z1, b, bn, a, an);
    else
        mpn_mul(z1, a, an, b, bn);
timeit_stop(timer);
old_time = timer->wall;

timeit_start(timer);
    mpn_mul_v2(z2, a, an, b, bn);
timeit_stop(timer);
new_time = timer->wall;

if (old_time > 30)
{
    std::cout << "----------- " << an << " * " << bn << " words -------------" << std::endl;
    std::cout << "old_time: " << old_time << std::endl;
    std::cout << "new_time: " << new_time << std::endl;
    std::cout << " old/new: " << double(old_time)/double(new_time) << std::endl;
}


    for (ulong i = 0; i < an + bn; i++)
    {
        if (z1[i] != z2[i])
        {
            std::cout << "mismatch!!!!" << std::endl;
            std::cout << "an: " << an << ", bn: " << bn << std::endl;
            std::cout << "z1[" << i << "]: " << format_hex(z1[i]) << std::endl;
            std::cout << "z2[" << i << "]: " << format_hex(z2[i]) << std::endl;
            abort();
        }
    }

    delete[] a;
    delete[] b;
    delete[] z1;
    delete[] z2;
}


/*
------- depth: 20 --------
new time: 7 ms, 5 ms
   again: 6 ms, 6 ms
old time: 7 ms, 7 ms
 old/new: 1, 1.4
------- depth: 21 --------
new time: 13 ms, 12 ms
   again: 13 ms, 12 ms
old time: 17 ms, 17 ms
 old/new: 1.30769, 1.41667
------- depth: 22 --------
new time: 27 ms, 25 ms
   again: 28 ms, 25 ms
old time: 37 ms, 37 ms
 old/new: 1.37037, 1.48
------- depth: 23 --------
new time: 56 ms, 51 ms
   again: 56 ms, 52 ms
old time: 82 ms, 83 ms
 old/new: 1.46429, 1.62745
------- depth: 24 --------
new time: 111 ms, 104 ms
   again: 111 ms, 104 ms
old time: 179 ms, 181 ms
 old/new: 1.61261, 1.74038
------- depth: 25 --------
new time: 231 ms, 217 ms
   again: 231 ms, 217 ms
old time: 391 ms, 393 ms
 old/new: 1.69264, 1.81106
------- depth: 26 --------
new time: 492 ms, 463 ms
   again: 491 ms, 463 ms
old time: 858 ms, 858 ms
 old/new: 1.7439, 1.85313
------- depth: 27 --------
new time: 1081 ms, 1034 ms
   again: 1083 ms, 1036 ms
Killed
*/
void test_v2(ulong L)
{
    ulong Xn = ulong(1)<<L;
    flint_rand_t state;
    timeit_t timer;
    ulong new_time, new_inv_time, old_time, old_inv_time;

    flint_randinit(state);

    if (L < LG_BLK_SZ)
        return;

    {
        double* X = new double[ulong(1)<<L];
        fftv2_ctx ctx(0x03f00000000001ULL);

    std::cout << "------- depth: " << L << " --------" << std::endl;

        ctx.set_depth(L);
        ctx.set_data(new double[ctx.data_size()]);

//    std::cout.precision(17);
//    std::cout << "ctx.p: " << ctx.p << std::endl;

        for (ulong i = 0; i < Xn; i++)
        {
            X[i] = i + 1;
            ctx.set_index(i, X[i]);
        }

    timeit_start(timer);
        ctx.transform_forward();
    timeit_stop(timer);
    new_time = FLINT_MAX(1, timer->wall);
    std::cout << "new time: " << new_time << " ms, ";

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
                std::cout << "forward oops " << reduce_pm1n_to_0n(ctx.get_index(i), ctx.p) << " != " << reduce_pm1n_to_0n(y, ctx.p) << std::endl;
                std::abort();
            }
        }

        {
            double ti = nmod_pow_ui(nmod_inv(2, ctx.mod), L, ctx.mod);
            if (ti > 0.5*ctx.p)
                ti -= ctx.p;
            for (ulong i = 0; i < Xn; i++)
                ctx.set_index(i, mulmod2(ctx.get_index(i), ti, ctx.p, ctx.pinv));

        timeit_start(timer);
            ctx.transform_reverse();
        timeit_stop(timer);
        new_inv_time = FLINT_MAX(1, timer->wall);
        std::cout << new_inv_time << " ms" << std::endl;

            // check answer
            for (ulong i = 0; i < Xn; i++)
            {
                double y = X[i];
                if (reduce_pm1n_to_0n(y, ctx.p) !=
                    reduce_pm1n_to_0n(reduce_to_pm1n(ctx.get_index(i), ctx.p, ctx.pinv), ctx.p))
                {
                    std::cout << "ctx.get_index(" << i << "): " << ctx.get_index(i) << std::endl;
                    std::cout << "y: " << y << std::endl;
                    std::cout << i << ": reverse oops " << reduce_pm1n_to_0n(ctx.get_index(i), ctx.p) << " != " << reduce_pm1n_to_0n(y, ctx.p) << std::endl;
                    std::abort();
                }
            }

        timeit_start(timer);
            ctx.transform_forward();
        timeit_stop(timer);
        std::cout << "   again: " << timer->wall << " ms";

        timeit_start(timer);
            ctx.transform_reverse();
        timeit_stop(timer);
        std::cout << ", " << timer->wall << " ms" << std::endl;

        }

        delete[] ctx.release_data();
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
    std::cout << "old time: " << old_time << " ms, ";

    timeit_start(timer);
        PD_ifft_short_1<4>(y, Xn, L, ctx);
    timeit_stop(timer);

    old_inv_time = FLINT_MAX(1, timer->wall/4);
    std::cout << old_inv_time << " ms" << std::endl;

        delete[] y;
    }

    std::cout << " old/new: " << double(old_time)/double(new_time) << ", " << double(old_inv_time)/double(new_inv_time) << std::endl;

    flint_randclear(state);
}
