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
    packed<double, 4> w = w2s[j]; \
    packed<double, 4> n = p; \
    packed<double, 4> ninv = pinv;

#define PD_RADIX_2_FORWARD_2X(X0, X1, Z0, Z1) \
{ \
    packed<double, 4> _x0, _x1, _y0, _y1, _z0, _z1, _w0, _w1; \
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
    packed<double, 4> _x0, _x1, _y0, _y1, _z0, _z1, _w0, _w1; \
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
    packed<double, 4> W  = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
    packed<double, 4> n = p; \
    packed<double, 4> ninv = pinv;

#define PD_RADIX_2_REVERSE_2X(X0, X1, Z0, Z1) \
{ \
    packed<double, 4> _x0, _x1, _y0, _y1; \
    packed<double, 4> _z0, _z1, _w0, _w1; \
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
    packed<double, 4> w = w2s[2*j]; \
    packed<double, 4> w2 = w2s[1*j]; \
    packed<double, 4> iw = w2s[2*j+1]; \
    packed<double, 4> n(p); \
    packed<double, 4> ninv(pinv);

#define PD_RADIX_4_FORWARD(X0, X1, X2, X3) \
{ \
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
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
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double, 4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
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

#include "itrunc4.cpp"
#include "itrunc3.cpp"
#include "itrunc2.cpp"


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
    packed<double, 4> W  = ((j) == 0) ? -w2s[0] : w2s[(2*(j)  )^(mask)]; \
    packed<double, 4> W2 = ((j) == 0) ? -w2s[0] : w2s[(  (j)  )^(mask>>1)]; \
    packed<double, 4> IW = ((j) == 0) ?  w2s[1] : w2s[(2*(j)+1)^(mask)]; \
    packed<double, 4> n(p); \
    packed<double, 4> ninv(pinv); \

#define PD_RADIX_4_REVERSE(X0, X1, X2, X3) \
{ \
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
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
    packed<double, 4> _x0, _x1, _x2, _x3, _y0, _y1, _y2, _y3; \
    packed<double, 4> _z0, _z1, _z2, _z3, _w0, _w1, _w2, _w3; \
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


template <>
inline void fftv2_ctx::fft_basecase<0>(double* X, ulong j)
{
}

template <>
inline void fftv2_ctx::ifft_basecase<0>(double* X, ulong j)
{
}

template <>
inline void fftv2_ctx::fft_basecase<1>(double* X, ulong j)
{
    RADIX_2_FORWARD_PARAM(j)
    RADIX_2_FORWARD(X+0, X+1);
}

template <>
inline void fftv2_ctx::ifft_basecase<1>(double* X, ulong j)
{
    RADIX_2_REVERSE_PARAM(j)
    RADIX_2_REVERSE(X+0, X+1);
}

template <>
inline void fftv2_ctx::fft_basecase<2>(double* X, ulong j)
{
    RADIX_4_FORWARD_PARAM(j)
    RADIX_4_FORWARD(X+0, X+1, X+2, X+3);
}

template <>
inline void fftv2_ctx::ifft_basecase<2>(double* X, ulong j)
{
    RADIX_4_REVERSE_PARAM(j)
    RADIX_4_REVERSE(X+0, X+1, X+2, X+3);
}

template <>
inline void fftv2_ctx::fft_basecase<4>(double* X, ulong j)
{
    PD_RADIX_4_FORWARD_PARAM(j)
    PD_RADIX_4_FORWARD(X+4*0, X+4*1, X+4*2, X+4*3);
    fft_basecase<2>(X+4*0, 4*j+0);
    fft_basecase<2>(X+4*1, 4*j+1);
    fft_basecase<2>(X+4*2, 4*j+2);
    fft_basecase<2>(X+4*3, 4*j+3);
}

template <>
inline void fftv2_ctx::ifft_basecase<4>(double* X, ulong j)
{
    ifft_basecase<2>(X+4*0, 4*j+0);
    ifft_basecase<2>(X+4*1, 4*j+1);
    ifft_basecase<2>(X+4*2, 4*j+2);
    ifft_basecase<2>(X+4*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    PD_RADIX_4_REVERSE(X+4*0, X+4*1, X+4*2, X+4*3);
}

template <>
inline void fftv2_ctx::fft_basecase<6>(double* X, ulong j)
{
    PD_RADIX_4_FORWARD_PARAM(j)
    for (ulong i = 0; i < 16; i += 8)
        PD_RADIX_4_FORWARD_2X(X+i+16*0, X+i+16*1, X+i+16*2, X+i+16*3, X+i+16*0+4, X+i+16*1+4, X+i+16*2+4, X+i+16*3+4);
    fft_basecase<4>(X+16*0, 4*j+0);
    fft_basecase<4>(X+16*1, 4*j+1);
    fft_basecase<4>(X+16*2, 4*j+2);
    fft_basecase<4>(X+16*3, 4*j+3);
}

template <>
inline void fftv2_ctx::ifft_basecase<6>(double* X, ulong j)
{
    ifft_basecase<4>(X+16*0, 4*j+0);
    ifft_basecase<4>(X+16*1, 4*j+1);
    ifft_basecase<4>(X+16*2, 4*j+2);
    ifft_basecase<4>(X+16*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    for (ulong i = 0; i < 16; i += 8)
        PD_RADIX_4_REVERSE_2X(X+i+16*0, X+i+16*1, X+i+16*2, X+i+16*3, X+i+16*0+4, X+i+16*1+4, X+i+16*2+4, X+i+16*3+4);
}

template <>
void fftv2_ctx::fft_basecase<8>(double* X, ulong j)
{
    PD_RADIX_4_FORWARD_PARAM(j)
    for (ulong i = 0; i < 64; i += 8)
        PD_RADIX_4_FORWARD_2X(X+i+64*0, X+i+64*1, X+i+64*2, X+i+64*3, X+i+64*0+4, X+i+64*1+4, X+i+64*2+4, X+i+64*3+4);
    fft_basecase<6>(X+64*0, 4*j+0);
    fft_basecase<6>(X+64*1, 4*j+1);
    fft_basecase<6>(X+64*2, 4*j+2);
    fft_basecase<6>(X+64*3, 4*j+3);
}

template <>
void fftv2_ctx::ifft_basecase<8>(double* X, ulong j)
{
    ifft_basecase<6>(X+64*0, 4*j+0);
    ifft_basecase<6>(X+64*1, 4*j+1);
    ifft_basecase<6>(X+64*2, 4*j+2);
    ifft_basecase<6>(X+64*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    for (ulong i = 0; i < 64; i += 8)
        PD_RADIX_4_REVERSE_2X(X+i+64*0, X+i+64*1, X+i+64*2, X+i+64*3, X+i+64*0+4, X+i+64*1+4, X+i+64*2+4, X+i+64*3+4);
}

