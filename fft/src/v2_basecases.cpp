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
    packed<double,4> w = w2s[j]; \
    packed<double,4> n = p; \
    packed<double,4> ninv = pinv;

#define PD_RADIX_2_FORWARD_2X(X0, X1, Z0, Z1) \
{ \
    packed<double,4> _x0, _x1, _y0, _y1, _z0, _z1, _w0, _w1; \
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

#define PD_RADIX_2_FORWARD_2X_ITRUNC2_OTRUNC1(X0, X1, Z0, Z1) \
{ \
    packed<double,4> _x0, _x1, _y0, _y1, _z0, _z1, _w0, _w1; \
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
    _y0.store(X0); \
        _w0.store(Z0); \
}


#define RADIX_2_REVERSE_PARAM(j) \
    ulong mask = saturate_bits(j); \
    double W  = (UNLIKELY((j) == 0)) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
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
    packed<double,4> W  = (UNLIKELY((j) == 0)) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
    packed<double,4> n = p; \
    packed<double,4> ninv = pinv;

#define PD_RADIX_2_REVERSE_2X(X0, X1, Z0, Z1) \
{ \
    packed<double,4> _x0, _x1, _y0, _y1; \
    packed<double,4> _z0, _z1, _w0, _w1; \
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
    packed<double,4> w = w2s[2*j]; \
    packed<double,4> w2 = w2s[1*j]; \
    packed<double,4> iw = w2s[2*j+1]; \
    packed<double,4> n(p); \
    packed<double,4> ninv(pinv);

#define PD_RADIX_4_FORWARD(X0, X1, X2, X3) \
{ \
    packed<double,4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
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
    packed<double,4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double,4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
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

template <int itrunc, int otrunc>
void fftv2_moth_trunc_block(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j)
{
    packed<double,4> w = Q->w2s[2*j];
    packed<double,4> w2 = Q->w2s[1*j];
    packed<double,4> iw = Q->w2s[2*j+1];
    packed<double,4> N(Q->p);
    packed<double,4> Ninv(Q->pinv);

    ulong i = 0;
    do {
        packed<double,4> x0(0.0), x1(0.0), x2(0.0), x3(0.0), y0, y1, y2, y3;
        packed<double,4> u0(0.0), u1(0.0), u2(0.0), u3(0.0), v0, v1, v2, v3;
        if (0 < itrunc) x0.load(X0+i+0);
            if (0 < itrunc) u0.load(X0+i+4);
        if (0 < itrunc) x0 = reduce_to_pm1n(x0, N, Ninv);
            if (0 < itrunc) u0 = reduce_to_pm1n(u0, N, Ninv);
        if (1 < itrunc) x1.load(X1+i+0);
            if (1 < itrunc) u1.load(X1+i+4);
        if (2 < itrunc) x2.load(X2+i+0);
            if (2 < itrunc) u2.load(X2+i+4);
        if (2 < itrunc) x2 = mulmod2(x2, w2, N, Ninv);
            if (2 < itrunc) u2 = mulmod2(u2, w2, N, Ninv);
        if (3 < itrunc) x3.load(X3+i+0);
            if (3 < itrunc) u3.load(X3+i+4);
        if (3 < itrunc) x3 = mulmod2(x3, w2, N, Ninv);
            if (3 < itrunc) u3 = mulmod2(u3, w2, N, Ninv);
        y0 = (2 < itrunc) ? add(x0, x2) : x0;
            v0 = (2 < itrunc) ? add(u0, u2) : u0;
        y1 = (3 < itrunc) ? add(x1, x3) : x1;
            v1 = (3 < itrunc) ? add(u1, u3) : u1;
        y2 = (2 < itrunc) ? sub(x0, x2) : x0;
            v2 = (2 < itrunc) ? sub(u0, u2) : u0;
        y3 = (3 < itrunc) ? sub(x1, x3) : x1;
            v3 = (3 < itrunc) ? sub(u1, u3) : u1;
        y1 = mulmod2(y1, w, N, Ninv);
            v1 = mulmod2(v1, w, N, Ninv);
        y3 = mulmod2(y3, iw, N, Ninv);
            v3 = mulmod2(v3, iw, N, Ninv);
        x0 = add(y0, y1);
            u0 = add(v0, v1);
        x1 = sub(y0, y1);
            u1 = sub(v0, v1);
        x2 = add(y2, y3);
            u2 = add(v2, v3);
        x3 = sub(y2, y3);
            u3 = sub(v2, v3);        
        if (0 < otrunc) x0.store(X0+i+0);
            if (0 < otrunc) u0.store(X0+i+4);        
        if (1 < otrunc) x1.store(X1+i+0);
            if (1 < otrunc) u1.store(X1+i+4);
        if (2 < otrunc) x2.store(X2+i+0);
            if (2 < otrunc) u2.store(X2+i+4);
        if (3 < otrunc) x3.store(X3+i+0);
            if (3 < otrunc) u3.store(X3+i+4);
    } while (i += 8, i < Q->blk_sz);
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

#define RADIX_4_REVERSE_PARAM(j, jm, j_can_be_0) \
    assert((j) == 0 || ((jm) == (j) ^ (saturate_bits(j)>>1))); \
    double W  = ((j_can_be_0) && UNLIKELY((j) == 0)) ? -w2s[0] : w2s[2*(jm)+1]; \
    double W2 = ((j_can_be_0) && UNLIKELY((j) == 0)) ? -w2s[0] : w2s[(jm)]; \
    double IW = ((j_can_be_0) && UNLIKELY((j) == 0)) ?  w2s[1] : w2s[2*(jm)+0]; \
    double n = p; \
    double ninv = pinv;

#define RADIX_4_REVERSE_PARAM_SURE(j, jm, j_is_0) \
    assert((j) == 0 || ((jm) == (j) ^ (saturate_bits(j)>>1))); \
    double W  = (j_is_zero) ? -w2s[0] : w2s[2*(jm)+1]; \
    double W2 = (j_is_zero) ? -w2s[0] : w2s[(jm)]; \
    double IW = (j_is_zero) ?  w2s[1] : w2s[2*(jm)+0]; \
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

#define PD_RADIX_4_REVERSE_PARAM(j, jm, j_can_be_0) \
    assert((j) == 0 || ((jm) == (j) ^ (saturate_bits(j)>>1))); \
    packed<double,4> W  = ((j_can_be_0) && UNLIKELY((j) == 0)) ? -w2s[0] : w2s[2*(jm)+1]; \
    packed<double,4> W2 = ((j_can_be_0) && UNLIKELY((j) == 0)) ? -w2s[0] : w2s[(jm)]; \
    packed<double,4> IW = ((j_can_be_0) && UNLIKELY((j) == 0)) ?  w2s[1] : w2s[2*(jm)]; \
    packed<double,4> n(p); \
    packed<double,4> ninv(pinv); \

#define PD_RADIX_4_REVERSE_PARAM_SURE(j, jm, j_is_0) \
    assert((j) == 0 || ((jm) == (j) ^ (saturate_bits(j)>>1))); \
    packed<double,4> W  = (j_is_0) ? -w2s[0] : w2s[2*(jm)+1]; \
    packed<double,4> W2 = (j_is_0) ? -w2s[0] : w2s[(jm)]; \
    packed<double,4> IW = (j_is_0) ?  w2s[1] : w2s[2*(jm)]; \
    packed<double,4> n(p); \
    packed<double,4> ninv(pinv); \

#define PD_RADIX_4_REVERSE(X0, X1, X2, X3) \
{ \
    packed<double,4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
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
    _x1.store(X1); \
    _x2 = sub(_y1, _y0); \
    _x3 = add(_y3, _y2); \
    _x0 = reduce_to_pm1n(_x0, n, ninv); \
    _x2 = mulmod2(_x2, W2, n, ninv); \
    _x3 = mulmod2(_x3, W2, n, ninv); \
    _x0.store(X0); \
    _x2.store(X2); \
    _x3.store(X3); \
}

#define PD_RADIX_4_REVERSE_2X(X0, X1, X2, X3, Z0, Z1, Z2, Z3) \
{ \
    packed<double,4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double,4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
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
    _x1.store(X1); \
        _z1.store(Z1); \
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
    _x2.store(X2); \
        _z2.store(Z2); \
    _x3.store(X3); \
        _z3.store(Z3); \
}

template<int depth, typename std::enable_if<depth == 0>::type*>
inline void fftv2_ctx::fft_basecase(double* X, ulong j)
{
}

template<int depth, bool j_can_be_0, typename std::enable_if<depth == 0>::type*>
inline void fftv2_ctx::ifft_basecase(double* X, ulong j, ulong jm)
{
}

template<int depth, typename std::enable_if<depth == 1>::type*>
inline void fftv2_ctx::fft_basecase(double* X, ulong j)
{
    RADIX_2_FORWARD_PARAM(j)
    RADIX_2_FORWARD(X+0, X+1);
}

template<int depth, bool j_can_be_0, typename std::enable_if<depth == 1>::type*>
inline void fftv2_ctx::ifft_basecase(double* X, ulong j, ulong jm)
{
    RADIX_2_REVERSE_PARAM(j)
    RADIX_2_REVERSE(X+0, X+1);
}

template<int depth, typename std::enable_if<depth == 2>::type*>
inline void fftv2_ctx::fft_basecase(double* X, ulong j)
{
    RADIX_4_FORWARD_PARAM(j)
    RADIX_4_FORWARD(X+0, X+1, X+2, X+3);
}


template<int depth, bool j_can_be_0, typename std::enable_if<depth == 2>::type*>
inline void fftv2_ctx::ifft_basecase(double* X, ulong j, ulong jm)
{
    RADIX_4_REVERSE_PARAM(j, jm, j_can_be_0)
    RADIX_4_REVERSE(X+0, X+1, X+2, X+3);
}

template<int depth, typename std::enable_if<depth == 4>::type*>
FORCE_INLINE inline void fftv2_ctx::fft_basecase(double* X, ulong j)
{
#if 0
    PD_RADIX_4_FORWARD_PARAM(j)
    PD_RADIX_4_FORWARD(X+4*0, X+4*1, X+4*2, X+4*3);
    fft_basecase<2>(X+4*0, 4*j+0);
    fft_basecase<2>(X+4*1, 4*j+1);
    fft_basecase<2>(X+4*2, 4*j+2);
    fft_basecase<2>(X+4*3, 4*j+3);
#else
    packed<double,4> w = w2s[2*j];
    packed<double,4> w2 = w2s[1*j];
    packed<double,4> iw = w2s[2*j+1];
    packed<double,4> N(p);
    packed<double,4> Ninv(pinv);
    packed<double,4> x0, x1, x2, x3, y0, y1, y2, y3, w_, w2_, iw_, s, t, u, v;
    x0.load(X+0);
    x0 = reduce_to_pm1n(x0, N, Ninv);
    x1.load(X+4);
    x2.load(X+8);
    x3.load(X+12);
    x2 = mulmod2(x2, w2, N, Ninv);
    x3 = mulmod2(x3, w2, N, Ninv);
    y0 = add(x0, x2);
    y1 = add(x1, x3);
    y2 = sub(x0, x2);
    y3 = sub(x1, x3);
    y1 = mulmod2(y1, w, N, Ninv);
    y3 = mulmod2(y3, iw, N, Ninv);
    x0 = add(y0, y1);
    x1 = sub(y0, y1);
    x2 = add(y2, y3);
    x3 = sub(y2, y3);

//    w_  = packed<double,4>(w2s[2*(4*j+0)+0], w2s[2*(4*j+1)+0], w2s[2*(4*j+2)+0], w2s[2*(4*j+3)+0]);
//    w2_ = packed<double,4>(w2s[1*(4*j+0)+0], w2s[1*(4*j+1)+0], w2s[1*(4*j+2)+0], w2s[1*(4*j+3)+0]);
//    iw_ = packed<double,4>(w2s[2*(4*j+0)+1], w2s[2*(4*j+1)+1], w2s[2*(4*j+2)+1], w2s[2*(4*j+3)+1]);
    u.load_aligned(w2s + 8*j+0);
    v.load_aligned(w2s + 8*j+4);
    w2_.load_aligned(w2s + 4*j);
    w_ = permute<0,2,1,3>(unpacklo(u, v));
    iw_ = permute<0,2,1,3>(unpackhi(u, v));

    // now horizonal action within x0, x1, x2, x3
    //y0 = packed<double,4>(x0[0], x1[0], x2[0], x3[0]);
    //y1 = packed<double,4>(x0[1], x1[1], x2[1], x3[1]);
    //y2 = packed<double,4>(x0[2], x1[2], x2[2], x3[2]);
    //y3 = packed<double,4>(x0[3], x1[3], x2[3], x3[3]);
    //x0 = y0;
    //x1 = y1;
    //x2 = y2;
    //x3 = y3;
    u = unpacklo(x0, x1);
    v = unpackhi(x0, x1);
    s = unpacklo(x2, x3);
    t = unpackhi(x2, x3);
    x0 = permute2<0,2>(u, s);
    x1 = permute2<0,2>(v, t);
    x2 = permute2<1,3>(u, s);
    x3 = permute2<1,3>(v, t);

    x0 = reduce_to_pm1n(x0, N, Ninv);
    x2 = mulmod2(x2, w2_, N, Ninv);
    x3 = mulmod2(x3, w2_, N, Ninv);
    y0 = add(x0, x2);
    y1 = add(x1, x3);
    y2 = sub(x0, x2);
    y3 = sub(x1, x3);
    y1 = mulmod2(y1, w_, N, Ninv);
    y3 = mulmod2(y3, iw_, N, Ninv);
    x0 = add(y0, y1);
    x1 = sub(y0, y1);
    x2 = add(y2, y3);
    x3 = sub(y2, y3);

    //X[ 0+0] = x0[0]; X[ 0+1] = x1[0]; X[ 0+2] = x2[0]; X[ 0+3] = x3[0];
    //X[ 4+0] = x0[1]; X[ 4+1] = x1[1]; X[ 4+2] = x2[1]; X[ 4+3] = x3[1];
    //X[ 8+0] = x0[2]; X[ 8+1] = x1[2]; X[ 8+2] = x2[2]; X[ 8+3] = x3[2];
    //X[12+0] = x0[3]; X[12+1] = x1[3]; X[12+2] = x2[3]; X[12+3] = x3[3];
    u = unpacklo(x0, x1);
    u.store_subset<0,1>(X + 0+0);
    u.store_subset<2,3>(X + 8+0);
    u = unpackhi(x0, x1);
    u.store_subset<0,1>(X + 4+0);
    u.store_subset<2,3>(X + 12+0);
    u = unpacklo(x2, x3);
    u.store_subset<0,1>(X + 0+2);
    u.store_subset<2,3>(X + 8+2);
    u = unpackhi(x2, x3);
    u.store_subset<0,1>(X + 4+2);
    u.store_subset<2,3>(X + 12+2);
#endif
}

template<int depth, bool j_is_0, typename std::enable_if<depth == 4>::type*>
FORCE_INLINE inline void fftv2_ctx::ifft_basecase(double* X, ulong j, ulong jm)
{
#if 0
    ifft_basecase<2>(X+4*0, 4*j+0);
    ifft_basecase<2>(X+4*1, 4*j+1);
    ifft_basecase<2>(X+4*2, 4*j+2);
    ifft_basecase<2>(X+4*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    PD_RADIX_4_REVERSE(X+4*0, X+4*1, X+4*2, X+4*3);
#else
    assert(j_is_0 == (j == 0));
    assert(j == 0 || jm == j^(saturate_bits(j)>>1));

//    double W_0  = ((j) == 0) ? -w2s[0] : w2s[8*jm+7];
//    double IW_0 = ((j) == 0) ?  w2s[1] : w2s[8*jm+6];
//    double W_1  = ((j) == 0) ?  w2s[3] : w2s[8*jm+5];
//    double IW_1 = ((j) == 0) ?  w2s[2] : w2s[8*jm+4];
    packed<double,4> Q0;
    if (j_is_0)
    {
        Q0 = packed<double,4>(-w2s[0], w2s[1], w2s[3], w2s[2]);
    }
    else
    {
        Q0.load_aligned(w2s + 8*jm+4);
        Q0 = permute<3,2,1,0>(Q0);
    }

//    double W_2  = ((j) == 0) ?  w2s[7] : w2s[8*jm+3];
//    double IW_2 = ((j) == 0) ?  w2s[6] : w2s[8*jm+2];
//    double W_3  = ((j) == 0) ?  w2s[5] : w2s[8*jm+1];
//    double IW_3 = ((j) == 0) ?  w2s[4] : w2s[8*jm+0];
    packed<double,4> Q1;
    Q1.load_aligned(w2s + (j_is_0 ? 4 : 8*jm+0));
    Q1 = permute<3,2,1,0>(Q1);

//    double W2_0 = ((j) == 0) ? -w2s[0] : w2s[4*jm+3];
//    double W2_1 = ((j) == 0) ?  w2s[1] : w2s[4*jm+2];
//    double W2_2 = ((j) == 0) ?  w2s[3] : w2s[4*jm+1];
//    double W2_3 = ((j) == 0) ?  w2s[2] : w2s[4*jm+0];
    packed<double,4> Q2;
    if (j_is_0)
        Q2 = permute<3,2,1,0>(Q0);
    else
        Q2.load_aligned(w2s + 4*jm+0);

    packed<double,4> N = p;
    packed<double,4> Ninv = pinv;
    packed<double,4> adxs, bcqr, ADXS, BCQR, eftu, EFTU, ghvw, GHVW;
    packed<double,4> fhuw, egtv, FHUW, EGTV, cbrq, CBRQ, cdrs, CDRS, aAxX;
    packed<double,4> x0, x1, x2, x3, y0, y1, y2, y3;
    adxs = packed<double,4>(X[0+0], X[0+3], X[0+4], X[0+7]);
    bcqr = packed<double,4>(X[0+1], X[0+2], X[0+5], X[0+6]);
    ADXS = packed<double,4>(X[8+0], X[8+3], X[8+4], X[8+7]);
    BCQR = packed<double,4>(X[8+1], X[8+2], X[8+5], X[8+6]);
    eftu = add(adxs, bcqr);
    EFTU = add(ADXS, BCQR);
    ghvw = sub(adxs, bcqr);
    GHVW = sub(ADXS, BCQR);

//    ghvw = mulmod2(ghvw, packed<double,4>(W_0, IW_0, W_1, IW_1), N, Ninv);
    ghvw = mulmod2(ghvw, Q0, N, Ninv);

//    GHVW = mulmod2(GHVW, packed<double,4>(W_2, IW_2, W_3, IW_3), N, Ninv);
    GHVW = mulmod2(GHVW, Q1, N, Ninv);

    fhuw = unpackhi(eftu, ghvw);
    egtv = unpacklo(eftu, ghvw);
    FHUW = unpackhi(EFTU, GHVW);
    EGTV = unpacklo(EFTU, GHVW);
    adxs = add(fhuw, egtv);
    ADXS = add(FHUW, EGTV);
    cbrq = sub(fhuw, egtv);
    CBRQ = sub(FHUW, EGTV);
    cdrs = blend<0,1,0,1>(cbrq, adxs);
    CDRS = blend<0,1,0,1>(CBRQ, ADXS);
//    cdrs = mulmod2(cdrs, packed<double,4>(W2_0, W2_0, W2_1, W2_1), N, Ninv);
//    CDRS = mulmod2(CDRS, packed<double,4>(W2_2, W2_2, W2_3, W2_3), N, Ninv);
    cdrs = mulmod2(cdrs, permute<3,3,2,2>(Q2), N, Ninv);
    CDRS = mulmod2(CDRS, permute<1,1,0,0>(Q2), N, Ninv);

    aAxX = packed<double,4>(unpacklo(adxs, ADXS));
    aAxX = reduce_to_pm1n(aAxX, N, Ninv);
    x0 = packed<double,4>(aAxX[0], cbrq[1], cdrs[0], cdrs[1]);
    x1 = packed<double,4>(aAxX[2], cbrq[3], cdrs[2], cdrs[3]);
    x2 = packed<double,4>(aAxX[1], CBRQ[1], CDRS[0], CDRS[1]);
    x3 = packed<double,4>(aAxX[3], CBRQ[3], CDRS[2], CDRS[3]);
    packed<double,4> W  = j_is_0 ? -w2s[0] : w2s[2*jm+1];
    packed<double,4> IW = j_is_0 ?  w2s[1] : w2s[2*jm+0];
    packed<double,4> W2 = j_is_0 ? -w2s[0] : w2s[jm];
    y0 = add(x0, x1);
    y1 = add(x2, x3);
    y2 = sub(x0, x1);
    y3 = sub(x3, x2);
    y2 = mulmod2(y2, W, N, Ninv);
    y3 = mulmod2(y3, IW, N, Ninv);
    x0 = add(y0, y1);
    x1 = sub(y3, y2);
    x2 = sub(y1, y0);
    x3 = add(y3, y2);
    x0 = reduce_to_pm1n(x0, N, Ninv);
    x2 = mulmod2(x2, W2, N, Ninv);
    x3 = mulmod2(x3, W2, N, Ninv);
    x0.store(X+0);
    x1.store(X+4);
    x2.store(X+8);
    x3.store(X+12);
#endif
}

template<int depth, typename std::enable_if<depth >= 6>::type*>
void fftv2_ctx::fft_basecase(double* X, ulong j)
{
    constexpr ulong l = ulong(1) << (depth-2);
    PD_RADIX_4_FORWARD_PARAM(j)
    for (ulong i = 0; i < l; i += 8)
        PD_RADIX_4_FORWARD_2X(X+i+0*l+0, X+i+1*l+0, X+i+2*l+0, X+i+3*l+0,
                              X+i+0*l+4, X+i+1*l+4, X+i+2*l+4, X+i+3*l+4);
    fft_basecase<depth-2>(X+0*l, 4*j+0);
    fft_basecase<depth-2>(X+1*l, 4*j+1);
    fft_basecase<depth-2>(X+2*l, 4*j+2);
    fft_basecase<depth-2>(X+3*l, 4*j+3);
}

template<int depth, bool j_is_0, typename std::enable_if<depth >= 6>::type*>
void fftv2_ctx::ifft_basecase(double* X, ulong j, ulong jm)
{
    assert(j_is_0 == (j == 0));
    assert(j == 0 || jm == j^(saturate_bits(j)>>1));

    constexpr ulong l = ulong(1) << (depth-2);
    ifft_basecase<depth-2, j_is_0>(X+0*l, 4*j+0, j_is_0 ? 0 : 4*jm+3);
    ifft_basecase<depth-2, false> (X+1*l, 4*j+1, j_is_0 ? 1 : 4*jm+2);
    ifft_basecase<depth-2, false> (X+2*l, 4*j+2, j_is_0 ? 3 : 4*jm+1);
    ifft_basecase<depth-2, false> (X+3*l, 4*j+3, j_is_0 ? 2 : 4*jm+0);
    PD_RADIX_4_REVERSE_PARAM_SURE(j, jm, j_is_0)
    for (ulong i = 0; i < l; i += 8)
        PD_RADIX_4_REVERSE_2X(X+i+0*l+0, X+i+1*l+0, X+i+2*l+0, X+i+3*l+0,
                              X+i+0*l+4, X+i+1*l+4, X+i+2*l+4, X+i+3*l+4);
}

void fftv2_ctx::fft_base(ulong I, ulong j)
{
    fft_basecase<LG_BLK_SZ>(from_index(I), j);
}

template<bool j_can_be_0>
void fftv2_ctx::ifft_base(ulong I, ulong j)
{
    if (j_can_be_0 && j == 0)
        ifft_basecase<LG_BLK_SZ, true>(from_index(I), j, j);
    else
        ifft_basecase<LG_BLK_SZ, false>(from_index(I), j, j^(saturate_bits(j)>>1));
}

///////////////////////////////////////////////////////////////////////////////

#define PD packed<double,VEC_SZ>

static void _missing_helper(int l, int ntrunc, int ztrunc, bool f)
{
    std::cout << std::endl << "function l = " << l << ", n = " << ntrunc << ", z = " << ztrunc << ", f = " << f;
    if (1 <= ztrunc && ztrunc <= l && ntrunc <= ztrunc && 1 <= ntrunc+f && ntrunc+f <= l)
        std::cout << " is not implmented" << std::endl;
    else
        std::cout << " does not exist and should not be called" << std::endl;
    std::abort();
}

static bool _would_be_nice_if_we_had_it_helper(int l, int ntrunc, int ztrunc, bool f)
{
    std::cout << std::endl << "function l = " << l << ", n = " << ntrunc << ", z = " << ztrunc << ", f = " << f;
    if (1 <= ztrunc && ztrunc <= l && ntrunc <= ztrunc && 1 <= ntrunc+f && ntrunc+f <= l)
    {
        std::cout << " is not implmented" << std::endl;
    }
    else
    {
        std::cout << " does not exist and should not be called" << std::endl;
        std::abort();
    }
    return false;
}

/// radix 2 ///

template<int ntrunc, int ztrunc, bool f>
static void radix_2_moth_inv_trunc_block(fftv2_ctx* Q, double* X0, double* X1, ulong j)
{
    _missing_helper(2, ntrunc, ztrunc, f);
}

// {x0, x1} = {2*x0 - w*x1, x0 - w*x1}
template<> void radix_2_moth_inv_trunc_block<1,2,true>(
    fftv2_ctx* Q,
    double* X0, double* X1,
    ulong j)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD w = Q->w2s[j];
    PD c = 2.0;

    ulong i = 0; do {
        PD u0, u1, v0, v1;
        u0.load(X0 + i);
        u1.load(X1 + i);
        u0 = reduce_to_pm1n(u0, N, Ninv);
        u1 = mulmod2(u1, w, N, Ninv);
        v0 = fmsub(c, u0, u1);
        v1 = sub(u0, u1);
        v0.store(X0 + i);
        v1.store(X1 + i);
    } while (i += VEC_SZ, i < Q->blk_sz);
}

// {x0} = {2*x0 - w*x1}
template<> void radix_2_moth_inv_trunc_block<1,2,false>(
    fftv2_ctx* Q,
    double* X0, double* X1,
    ulong j)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD w = Q->w2s[j];
    PD c = 2.0;

    ulong i = 0; do {
        PD u0, u1, v0, v1;
        u0.load(X0 + i);
        u1.load(X1 + i);
        u0 = reduce_to_pm1n(u0, N, Ninv);
        u1 = mulmod2(u1, w, N, Ninv);
        v0 = fmsub(c, u0, u1);
        v0.store(X0 + i);
    } while (i += VEC_SZ, i < Q->blk_sz);
}

// {x0, x1} = {2*x0, x0}
template<> void radix_2_moth_inv_trunc_block<1,1,true>(
    fftv2_ctx* Q,
    double* X0, double* X1,
    ulong j)
{
    PD N = Q->p, Ninv = Q->pinv;

    ulong i = 0; do {
        PD u0, v0, v1;
        u0.load(X0 + i);
        u0 = reduce_to_pm1n(u0, N, Ninv);
        v0 = add(u0, u0);
        v1 = u0;
        v0.store(X0 + i);
        v1.store(X1 + i);
    } while (i += VEC_SZ, i < Q->blk_sz);
}


// {x0} = {2*x0}
template<> void radix_2_moth_inv_trunc_block<1,1,false>(
    fftv2_ctx* Q,
    double* X0, double* X1,
    ulong j)
{
    PD N = Q->p, Ninv = Q->pinv;

    ulong i = 0; do {
        PD u0, v0;
        u0.load(X0 + i);
        u0 = reduce_to_pm1n(u0, N, Ninv);
        v0 = add(u0, u0);
        v0.store(X0 + i);
    } while (i += VEC_SZ, i < Q->blk_sz);
}

// {x0} = {(x0 + w*x1)/2}
template<> void radix_2_moth_inv_trunc_block<0,2,true>(
    fftv2_ctx* Q,
    double* X0, double* X1,
    ulong j)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD w = Q->w2s[j];
    PD c = fnmadd(0.5, Q->p, 0.5);

    ulong i = 0; do {
        PD u0, u1;
        u0.load(X0 + i);
        u1.load(X1 + i);
        u1 = mulmod2(u1, w, N, Ninv);
        u0 = mulmod2(add(u0, u1), c, N, Ninv);
        u0.store(X0 + i);
    } while (i += VEC_SZ, i < Q->blk_sz);
}

// {x0} = {(x0 + w*x1)/2}
template<> void radix_2_moth_inv_trunc_block<0,1,true>(
    fftv2_ctx* Q,
    double* X0, double* X1,
    ulong j)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD w = Q->w2s[j];
    PD c = fnmadd(0.5, Q->p, 0.5);

    ulong i = 0; do {
        PD u0, u1;
        u0.load(X0 + i);
        u0 = mulmod2(u0, c, N, Ninv);
        u0.store(X0 + i);
    } while (i += VEC_SZ, i < Q->blk_sz);
}

/// radix 4 ///

template<int ntrunc, int ztrunc, bool f>
static bool radix_4_moth_inv_trunc_block(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    return _would_be_nice_if_we_had_it_helper(4, ntrunc, ztrunc, f);
}

/*
k = 2, n = 3, z = 4, f = true
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
[          -r               r         1   r*w^3]
*/
template<> bool radix_4_moth_inv_trunc_block<3,4,true>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    PD f0 = Q->w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2
    PD fr = -Q->w2s[2*j+1];     // -r*w
    PD fq = -Q->w2s[j];         // -w^2
    PD fp = reduce_pm1n_to_pmhn(mulmod2(fr, fq, N, Ninv), N);   // r*w^3
    
    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, v0, p0,q0,r0, a1, b1, c1, d1, u1, v1;
        a0.load(X0+i);
        b0.load(X1+i);
        c0.load(X2+i);
        d0.load(X3+i);
        u0 = add(a0, b0);
        v0 = sub(a0, b0);
        p0 = mulmod2(d0, fp, N, Ninv);
        q0 = mulmod2(d0, fq, N, Ninv);
        r0 = mulmod2(d0, fr, N, Ninv);
        c0 = reduce_to_pm1n(c0, N, Ninv);
        u0 = reduce_to_pm1n(u0, N, Ninv);
        b0 = mulmod2(v0, f1, N, Ninv);
        v0 = mulmod2(v0, f0, N, Ninv);
        d0 = sub(c0, v0);
        c0 = fmsub(f2, c0, v0);
        a0 = add(c0, u0);
        c0 = sub(c0, u0);                
        c0 = mulmod2(c0, f3, N, Ninv);
        (a0+p0).store(X0+i);
        (b0+q0).store(X1+i);
        (c0+r0).store(X2+i);
        (d0+p0).store(X3+i);
    } while (i += VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 3, z = 4, f = false
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
*/
template<> bool radix_4_moth_inv_trunc_block<3,4,false>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    PD f0 = Q->w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2
    PD fr = -Q->w2s[2*j+1];     // -r*w
    PD fq = -Q->w2s[j];         // -w^2
    PD fp = reduce_pm1n_to_pmhn(mulmod2(fr, fq, N, Ninv), N);   // r*w^3
    ulong i = 0;
    do {
        PD a0, b0, c0, d0, u0, v0, p0,q0,r0, a1, b1, c1, d1, u1, v1;
        a0.load(X0+i);
        b0.load(X1+i);
        c0.load(X2+i);
        d0.load(X3+i);
        u0 = add(a0, b0);
        v0 = sub(a0, b0);
        p0 = mulmod2(d0, fp, N, Ninv);
        q0 = mulmod2(d0, fq, N, Ninv);
        r0 = mulmod2(d0, fr, N, Ninv);
        c0 = reduce_to_pm1n(c0, N, Ninv);
        u0 = reduce_to_pm1n(u0, N, Ninv);
        b0 = mulmod2(v0, f1, N, Ninv);
        v0 = mulmod2(v0, f0, N, Ninv);
        d0 = sub(c0, v0);
        c0 = fmsub(f2, c0, v0);
        a0 = add(c0, u0);
        c0 = sub(c0, u0);                
        c0 = mulmod2(c0, f3, N, Ninv);
        (a0+p0).store(X0+i);
        (b0+q0).store(X1+i);
        (c0+r0).store(X2+i);
    } while (i += VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 3, z = 3, f = true
[      -r + 1           r + 1         2]
[        2//w           -2//w         0]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2]
[          -r               r         1]

    {x0, x1, x3, x4} = {        -r*(x0 - x1) + (x0 + x1) + 2*x2,
                        2*w^-1*(x0 - x1),
                         -w^-2*(-r*(x0 - x1) - (x0 + x1) + 2*x2),
                                -r*(x0 - x1)             +   x2  }
*/
template<> bool radix_4_moth_inv_trunc_block<3,3,true>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    PD f0 = Q->w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2
    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, v0,  a1, b1, c1, d1, u1, v1;
        a0.load(X0+i);
            a1.load(X0+i+VEC_SZ);
        b0.load(X1+i);
            b1.load(X1+i+VEC_SZ);
        c0.load(X2+i);
            c1.load(X2+i+VEC_SZ);
        v0 = sub(a0, b0);
            v1 = sub(a1, b1);
        mulmod2(v0, f1, N, Ninv).store(X1+i);
            mulmod2(v1, f1, N, Ninv).store(X1+i+VEC_SZ);
        c0 = reduce_to_pm1n(c0, N, Ninv);
            c1 = reduce_to_pm1n(c1, N, Ninv);
        v0 = mulmod2(v0, f0, N, Ninv);
            v1 = mulmod2(v1, f0, N, Ninv);
        sub(c0, v0).store(X3+i);
            sub(c1, v1).store(X3+i+VEC_SZ);
        u0 = reduce_to_pm1n(add(a0, b0), N, Ninv);
            u1 = reduce_to_pm1n(add(a1, b1), N, Ninv);
        c0 = fmsub(f2, c0, v0);
            c1 = fmsub(f2, c1, v1);
        a0 = add(c0, u0);
            a1 = add(c1, u1);
        c0 = sub(c0, u0);                
            c1 = sub(c1, u1);                
        c0 = mulmod2(c0, f3, N, Ninv);
            c1 = mulmod2(c1, f3, N, Ninv);
        a0.store(X0+i);
            a1.store(X0+i+VEC_SZ);
        c0.store(X2+i);
            c1.store(X2+i+VEC_SZ);
    } while (i += 2*VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 3, z = 3, f = false
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
template<> bool radix_4_moth_inv_trunc_block<3,3,false>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD W  = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[jm];
    PD f0 = Q->w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2
    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, v0,  a1, b1, c1, d1, u1, v1;
        a0.load(X0+i);
            a1.load(X0+i+VEC_SZ);
        b0.load(X1+i);
            b1.load(X1+i+VEC_SZ);
        c0.load(X2+i);
            c1.load(X2+i+VEC_SZ);
        v0 = sub(a0, b0);
            v1 = sub(a1, b1);
        mulmod2(v0, f1, N, Ninv).store(X1+i);
            mulmod2(v1, f1, N, Ninv).store(X1+i+VEC_SZ);
        c0 = reduce_to_pm1n(c0, N, Ninv);
            c1 = reduce_to_pm1n(c1, N, Ninv);
        v0 = mulmod2(v0, f0, N, Ninv);
            v1 = mulmod2(v1, f0, N, Ninv);
        u0 = reduce_to_pm1n(add(a0, b0), N, Ninv);
            u1 = reduce_to_pm1n(add(a1, b1), N, Ninv);
        c0 = fmsub(f2, c0, v0);
            c1 = fmsub(f2, c1, v1);
        a0 = add(c0, u0);
            a1 = add(c1, u1);
        c0 = sub(c0, u0);                
            c1 = sub(c1, u1);                
        c0 = mulmod2(c0, f3, N, Ninv);
            c1 = mulmod2(c1, f3, N, Ninv);
        a0.store(X0+i);
            a1.store(X0+i+VEC_SZ);
        c0.store(X2+i);
            c1.store(X2+i+VEC_SZ);
    } while (i += 2*VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 2, z = 4, f = true
[            2                2        -w^2             0]
[         2//w            -2//w           0          -w^2]
[1//2*r + 1//2   -1//2*r + 1//2   -1//2*w^2   -1//2*r*w^3]
*/
template<> bool radix_4_moth_inv_trunc_block<2,4,true>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD W = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD f0 = 2.0;
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = fnmadd(0.5, N, 0.5);    // 1//2
    PD f3 = Q->w2s[1];   // r
    PD f4 = Q->w2s[j];     // w^2
    PD f5 = reduce_pm1n_to_pmhn(mulmod2(f4, Q->w2s[2*j+1], N, Ninv), N);   // r*w^3

    ulong i = 0; do {
        PD a0,b0,u0,v0,s0,t0,g0,h0,p0,q0,r0, u1, v1, s1, t1;
        u0.load(X0+i);
        v0.load(X1+i);
        a0.load(X2+i);
        b0.load(X3+i);
        p0 = mulmod2(a0, f4, N, Ninv);
        q0 = mulmod2(b0, f4, N, Ninv);
        r0 = mulmod2(b0, f5, N, Ninv);
        s0 = reduce_to_pm1n(add(u0, v0), N, Ninv);
        t0 = sub(u0, v0);
        g0 = mulmod2(s0, f0, N, Ninv);
        h0 = mulmod2(t0, f1, N, Ninv);
        t0 = mulmod2(t0, f3, N, Ninv);
        (g0-p0).store(X0+i);
        (h0-q0).store(X1+i);
        mulmod2((s0+t0)-(p0+r0), f2, N, Ninv).store(X2+i);
    } while (i += VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 2, z = 4, f = false
[   2       2   -w^2      0]
[2//w   -2//w      0   -w^2]
*/
template<> bool radix_4_moth_inv_trunc_block<2,4,false>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD mwi = UNLIKELY(j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD w2 = Q->w2s[j];
    PD twowi = reduce_pm1n_to_pmhn(-2.0*mwi, N);

    ulong i = 0; do {
        // BAD
        PD a0, b0, c0, d0, u0, v0, a1, b1, c1, d1, u1, v1;
        a0.load(X0+i);
            a1.load(X0+i+VEC_SZ);
        b0.load(X1+i);
            b1.load(X1+i+VEC_SZ);
        c0.load(X2+i);
            c1.load(X2+i+VEC_SZ);
        d0.load(X3+i);
            d1.load(X3+i+VEC_SZ);
        c0 = mulmod2(c0, w2, N, Ninv);
            c1 = mulmod2(c1, w2, N, Ninv);
        d0 = mulmod2(d0, w2, N, Ninv);
            d1 = mulmod2(d1, w2, N, Ninv);
        u0 = add(a0, b0);
            u1 = add(a1, b1);
        v0 = sub(a0, b0);
            v1 = sub(a1, b1);
        u0 = reduce_to_pm1n(u0+u0, N, Ninv);
            u1 = reduce_to_pm1n(u1+u1, N, Ninv);
        v0 = mulmod2(v0, twowi, N, Ninv);
            v1 = mulmod2(v1, twowi, N, Ninv);
        u0 = sub(u0, c0);
            u1 = sub(u1, c1);
        v0 = sub(v0, d0);
            v1 = sub(v1, d1);
        u0.store(X0+i);
            u1.store(X0+i+VEC_SZ);
        v0.store(X1+i);
            v1.store(X1+i+VEC_SZ);
    } while (i += 2*VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 2, z = 2, f = true
[            2                2]
[         2//w            -2//w]
[1//2*r + 1//2   -1//2*r + 1//2]

{x0, x1} = {2*(x0 + x1), 2*w^-1*(x0 - x1), (x0+x1)/2 + (x0-x1)*i/2}
*/
template<> bool radix_4_moth_inv_trunc_block<2,2,true>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    double W = (j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD c0 = 2.0;
    PD c1 = reduce_pm1n_to_pmhn(-2.0*W, Q->p);
    double ha = fnmadd(0.5, Q->p, 0.5);
    PD c2 = ha;
    PD c3 = reduce_pm1n_to_pmhn(mulmod2(Q->w2s[1], ha, Q->p, Q->pinv), Q->p);

    ulong i = 0; do {
        PD u0, v0, s0, t0, u1, v1, s1, t1;
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
    } while (i += 2*VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 2, z = 2, f = true
[            2                2]
[         2//w            -2//w]
[1//2*r + 1//2   -1//2*r + 1//2]

{x0, x1} = {2*(x0 + x1), 2*w^-1*(x0 - x1), (x0+x1)/2 + (x0-x1)*i/2}
*/
template<> bool radix_4_moth_inv_trunc_block<2,2,false>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    double W = (j == 0) ? -Q->w2s[0] : Q->w2s[2*jm+1];
    PD c0 = 2.0;
    PD c1 = reduce_pm1n_to_pmhn(-2.0*W, Q->p);
    double ha = fnmadd(0.5, Q->p, 0.5);
    PD c2 = ha;
    PD c3 = reduce_pm1n_to_pmhn(mulmod2(Q->w2s[1], ha, Q->p, Q->pinv), Q->p);

    ulong i = 0; do {
        PD u0, v0, s0, t0, u1, v1, s1, t1;
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
        u0.store(X0 + i);
            u1.store(X0 + i + VEC_SZ);
        v0.store(X1 + i);
            v1.store(X1 + i + VEC_SZ);
    } while (i += 2*VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 1, z = 4, f = true
[4        -w   -w^2        -w^3]
[1   -1//2*w      0   -1//2*w^3]
*/
template<> bool radix_4_moth_inv_trunc_block<1,4,true>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD f2 = 2.0;
    PD w2 = Q->w2s[j];
    double ha = fnmadd(0.5, Q->p, 0.5);
    PD wo2 = reduce_pm1n_to_pmhn(mulmod2(Q->w2s[2*j], ha, N, Ninv), N);

    ulong i = 0; do {
        PD a0, b0, c0, d0, u0;
        a0.load(X0+i);
        a0 = reduce_to_pm1n(a0, N, Ninv);
        b0.load(X1+i);
        c0.load(X2+i);
        d0.load(X3+i);
        c0 = mulmod2(c0, w2, N, Ninv);
        d0 = mulmod2(d0, w2, N, Ninv);
        b0 = add(b0, d0);
        b0 = mulmod2(b0, wo2, N, Ninv);
        u0 = fmsub(f2, a0, b0);
        b0 = sub(a0, b0);
        a0 = reduce_to_pm1n(add(u0, u0), N, Ninv);
        a0 = sub(a0, c0);
        a0.store(X0+i);
        b0.store(X1+i);
    } while (i += VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 1, z = 4, f = false
[4   -w   -w^2   -w^3]
*/
template<> bool radix_4_moth_inv_trunc_block<1,4,false>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD f1 = 4.0;
    PD w2 = Q->w2s[j];
    PD w  = Q->w2s[2*j];

    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, a1, b1, c1, d1, u1;
        a0.load(X0+i);
        b0.load(X1+i);
        c0.load(X2+i);
        d0.load(X3+i);
        d0 = mulmod2(d0, w, N, Ninv);
        a0 = mulmod2(a0, f1, N, Ninv);
        b0 = mulmod2(b0, w, N, Ninv);
        a0 = sub(a0, b0);
        c0 = mulmod2(add(c0, d0), w2, N, Ninv);
        sub(a0, c0).store(X0+i);
    } while (i += VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 1, z = 1, f = true
[4]
[1]
*/
template<> bool radix_4_moth_inv_trunc_block<1,1,true>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD f = 4.0;

    ulong i = 0; do {
        PD a0, a1;
        a0.load(X0+i);
            a1.load(X0+i+VEC_SZ);
        reduce_to_pm1n(f*a0, N, Ninv).store(X0+i);
            reduce_to_pm1n(f*a1, N, Ninv).store(X0+i+VEC_SZ);
        reduce_to_pm1n(a0, N, Ninv).store(X1+i);
            reduce_to_pm1n(a1, N, Ninv).store(X1+i+VEC_SZ);
    } while (i += 2*VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
k = 2, n = 1, z = 1, f = false
[4]
*/
template<> bool radix_4_moth_inv_trunc_block<1,1,false>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD f = 4.0;

    ulong i = 0; do {
        PD a0, a1;
        a0.load(X0+i);
            a1.load(X0+i+VEC_SZ);
        reduce_to_pm1n(f*a0, N, Ninv).store(X0+i);
            reduce_to_pm1n(f*a1, N, Ninv).store(X0+i+VEC_SZ);
    } while (i += 2*VEC_SZ, i < Q->blk_sz);

    return true;
}

/*
    k = 2, n = 0, z = 4, f = true
    [1//4   1//4*w   1//4*w^2   1//4*w^3]

    {x0} = {1/4*x0 + w/4*x1 + w^2/4*(x2 + w*x3)}
*/
template<> bool radix_4_moth_inv_trunc_block<0,4,true>(
    fftv2_ctx* Q,
    double* X0, double* X1, double* X2, double* X3,
    ulong j, ulong jm)
{
    PD N = Q->p, Ninv = Q->pinv;
    PD ff0 = fnmadd(0.25, N, 0.25);
    PD ww = Q->w2s[2*j];
    PD ww2 = Q->w2s[j];
    PD f0 = ff0;
    PD  wo4 = reduce_pm1n_to_pmhn(mulmod2(ww, ff0, N, Ninv), N);
    PD w2o4 = reduce_pm1n_to_pmhn(mulmod2(ww2, ff0, N, Ninv), N);
    PD w = ww;

    ulong i = 0; do {
        PD a0, b0, c0, d0, a1, b1, c1, d1;
        d0.load(X3+i);
        a0.load(X0+i);
        b0.load(X1+i);
        c0.load(X2+i);
        d0 = mulmod2(d0, w, N, Ninv);
        a0 = mulmod2(a0, f0, N, Ninv);
        c0 = add(c0, d0);
        b0 = mulmod2(b0, wo4, N, Ninv);
        a0 = add(a0, b0);
        c0 = mulmod2(c0, w2o4, N, Ninv);
        a0 = add(a0, c0);
        a0.store(X0+i);
    } while (i += VEC_SZ, i < Q->blk_sz);

    return true;
}

#undef PD
