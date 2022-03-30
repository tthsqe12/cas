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

#if 0
inline void fftv2_ctx::basecase2_special(double* X, ulong j)
{
    double w = w2s[2*j];
    double w2 = w2s[1*j];
    double iw = w2s[2*j+1];
    double n = p;
    double ninv = pinv;

    double W = w2s[2*(j+1)];
    double W2 = w2s[1*(j+1)];
    double IW = w2s[2*(j+1)+1];

    double a = X[0];
    double b = X[1];
    double c = X[2];
    double d = X[3];
    double A = X[4];
    double B = X[5];
    double C = X[6];
    double D = X[7];

    a = reduce_to_pm1n(a, n, ninv);
    A = reduce_to_pm1n(A, n, ninv);

    c = mulmod2(c, w2, n, ninv);
    d = mulmod2(d, w2, n, ninv);
    C = mulmod2(C, W2, n, ninv);
    D = mulmod2(D, W2, n, ninv);

    double f = a - c;
    double g = a + c;
    double h = b - d;
    double q = b + d;
    double F = A - C;
    double G = A + C;
    double H = B - D;
    double Q = B + D;

    h = mulmod2(h, -iw, n, ninv);
    H = mulmod2(H, -IW, n, ninv);
    q = mulmod2(q, -w, n, ninv);
    Q = mulmod2(Q, -W, n, ninv);

    double b0 = g - q;
    double b1 = g + q;
    double b2 = f - h;
    double b3 = f + h;

    double B0 = G - Q;
    double B1 = G + Q;
    double B2 = F - H;
    double B3 = F + H;

    X[0] = b0;
    X[1] = b1;
    X[2] = b2;
    X[3] = b3;
    X[4] = B0;
    X[5] = B1;
    X[6] = B2;
    X[7] = B3;
}
#endif

inline void fftv2_ctx::fft_basecase2_2x(double* X, ulong j)
{
    double w = w2s[2*j];
    double w2 = w2s[1*j];
    double iw = w2s[2*j+1];
    double n = p;
    double ninv = pinv;
    packed<double,4> N = p;
    packed<double,4> Ninv = pinv;
    packed<double,4> U, V, S, T;

    double W = w2s[2*(j+1)];
    double W2 = w2s[1*(j+1)];
    double IW = w2s[2*(j+1)+1];

    double a = X[0];
    double b = X[1];
    double c = X[2];
    double d = X[3];
    double A = X[4];
    double B = X[5];
    double C = X[6];
    double D = X[7];

    a = reduce_to_pm1n(a, n, ninv);
    A = reduce_to_pm1n(A, n, ninv);
    U = packed<double, 4>(b, B, a, A);
    V = mulmod2(packed<double,4>(d, D, c, C), packed<double,4>(w2, W2, w2, W2), N, Ninv);
    S = sub(U, V);
    T = add(U, V);
    V = permute2f128(S, T, 3 + 16*1);
    U = mulmod2(insertf128(T, S, 1), packed<double,4>(-w, -W, -iw, -IW), N, Ninv);
    S = addsub(unpacklo(V, V), unpacklo(U, U));
    T = addsub(unpackhi(V, V), unpackhi(U, U));
    S.store(X + 0);
    T.store(X + 4);
}

inline void fftv2_ctx::fft_basecase2_4x(double* X, ulong j)
{
    packed<double,4> N = p;
    packed<double,4> Ninv = pinv;
    packed<double,4> U0, V0, S0, T0, Q0, R0, U1, V1, S1, T1, Q1, R1;
    double w, w2, iw, W, W2, IW;

//    w2 = w2s[1*j];
//    W2 = w2s[1*(j+1)];
//    Q0 = packed<double,4>(w2, W2, w2, W2);
//    w2 = w2s[1*(j+2)];
//    W2 = w2s[1*(j+3)];
//    Q1 = packed<double,4>(w2, W2, w2, W2);
    Q0.load(w2s + j);
    Q1 = permute2f128(Q0, Q0, 1 + 16*1);
    Q0 = permute2f128(Q0, Q0, 0 + 16*0);

    S0 = packed<double,4>(X[0], X[4], X[8+0], X[8+4]);
    S0 = reduce_to_pm1n(S0, N, Ninv);

//    U0 = packed<double, 4>(X[1], X[5], S[0], S[1]);
    U0 = packed<double, 4>(X[1], X[5], 0.0, 0.0);
    U0 = insertf128(U0, S0, 1);

//    U1 = packed<double, 4>(X[8+1], X[8+5], S[2], S[3]);
    U1 = packed<double, 4>(X[8+1], X[8+5], 0.0, 0.0);
    U1 = blend<0,0,1,1>(U1, S0);

    V0 = packed<double,4>(X[3], X[7], X[2], X[6]);
    V1 = packed<double,4>(X[8+3], X[8+7], X[8+2], X[8+6]);
    V0 = mulmod2(V0, Q0, N, Ninv);
    V1 = mulmod2(V1, Q1, N, Ninv);

//    w = w2s[2*j];
//    iw = w2s[2*j+1];
//    W = w2s[2*(j+1)];
//    IW = w2s[2*(j+1)+1];
//    R0 = neg(packed<double,4>(w, W, iw, IW));
    R0.load(w2s + 2*j);
    R0 = neg(permute<0,2,1,3>(R0));

//    w = w2s[2*(j+2)];
//    W = w2s[2*(j+3)];
//    iw = w2s[2*(j+2)+1];
//    IW = w2s[2*(j+3)+1];
//    R1 = neg(packed<double,4>(w, W, iw, IW));
    R1.load(w2s + 2*(j+2));
    R1 = neg(permute<0,2,1,3>(R1));

    S0 = sub(U0, V0);
    S1 = sub(U1, V1);
    T0 = add(U0, V0);
    T1 = add(U1, V1);
    V0 = permute2f128(S0, T0, 3 + 16*1);
    V1 = permute2f128(S1, T1, 3 + 16*1);
    U0 = mulmod2(insertf128(T0, S0, 1), R0, N, Ninv);
    U1 = mulmod2(insertf128(T1, S1, 1), R1, N, Ninv);
    S0 = addsub(unpacklo(V0, V0), unpacklo(U0, U0));
    S1 = addsub(unpacklo(V1, V1), unpacklo(U1, U1));
    T0 = addsub(unpackhi(V0, V0), unpackhi(U0, U0));
    T1 = addsub(unpackhi(V1, V1), unpackhi(U1, U1));

    S0.store(X + 0);
    T0.store(X + 4);
    S1.store(X + 8+0);
    T1.store(X + 8+4);
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
#if 0
    fft_basecase<2>(X+4*0, 4*j+0);
    fft_basecase<2>(X+4*1, 4*j+1);
    fft_basecase<2>(X+4*2, 4*j+2);
    fft_basecase<2>(X+4*3, 4*j+3);
#elseif 0
    fft_basecase2_2x(X+4*0, 4*j+0);
    fft_basecase2_2x(X+4*2, 4*j+2);
#else
    fft_basecase2_4x(X+4*0, 4*j+0);
#endif
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

    double a, b, c, d, e, f, g, h;
    double x, q, r, s, t, u, v, w;
    a = X[0];
    b = X[1];
    c = X[2];
    d = X[3];
    x = X[4+0];
    q = X[4+1];
    r = X[4+2];
    s = X[4+3];
    e = a + b;
    f = c + d;
    g = a - b;
    h = d - c;
    t = x + q;
    u = r + s;
    v = x - q;
    w = s - r;
    g = mulmod2(g, W_0, n, ninv);
    h = mulmod2(h, IW_0, n, ninv);
    v = mulmod2(v, W_1, n, ninv);
    w = mulmod2(w, IW_1, n, ninv);
    a = e + f;
    b = h - g;
    c = f - e;
    d = h + g;
    x = t + u;
    q = w - v;
    r = u - t;
    s = w + v;
    c = mulmod2(c, W2_0, n, ninv);
    d = mulmod2(d, W2_0, n, ninv);
    r = mulmod2(r, W2_1, n, ninv);
    s = mulmod2(s, W2_1, n, ninv);
    a = reduce_to_xm1n(a, n, ninv);
    x = reduce_to_xm1n(x, n, ninv);
    X[0] = b0;
    X[1] = b;
    X[2] = c;
    X[3] = d;
    X[4+0] = x;
    X[4+1] = q;
    X[4+2] = r;
    X[4+3] = s;
*/

inline void fftv2_ctx::ifft_basecase2_4x(double* io, ulong j)
{
    ulong mask = saturate_bits(j);
    ulong mask1 = saturate_bits(j+1);
    ulong mask2 = saturate_bits(j+2);
    ulong mask3 = saturate_bits(j+3);

    double W_0  = ((j) == 0) ? -w2s[0] : w2s[(2*(j+0)  )^(mask)];
    double W2_0 = ((j) == 0) ? -w2s[0] : w2s[(  (j+0)  )^(mask>>1)];
    double IW_0 = ((j) == 0) ?  w2s[1] : w2s[(2*(j+0)+1)^(mask)];
    double W_1  =                        w2s[(2*(j+1)  )^(mask1)];
    double W2_1 =                        w2s[(  (j+1)  )^(mask1>>1)];
    double IW_1 =                        w2s[(2*(j+1)+1)^(mask1)];
    double W_2  =                        w2s[(2*(j+2)  )^(mask2)];
    double W2_2 =                        w2s[(  (j+2)  )^(mask2>>1)];
    double IW_2 =                        w2s[(2*(j+2)+1)^(mask2)];
    double W_3  =                        w2s[(2*(j+3)  )^(mask3)];
    double W2_3 =                        w2s[(  (j+3)  )^(mask3>>1)];
    double IW_3 =                        w2s[(2*(j+3)+1)^(mask3)];
    double n = p;
    double ninv = pinv;
    double a, b, c, d, e, f, g, h;
    double x, q, r, s, t, u, v, w;
    double A, B, C, D, E, F, G, H;
    double X, Q, R, S, T, U, V, W;

    packed<double,4> N = p;
    packed<double,4> Ninv = pinv;

    a = io[0];
    b = io[1];
    c = io[2];
    d = io[3];
    x = io[4+0];
    q = io[4+1];
    r = io[4+2];
    s = io[4+3];
        A = io[8+0];
        B = io[8+1];
        C = io[8+2];
        D = io[8+3];
        X = io[8+4+0];
        Q = io[8+4+1];
        R = io[8+4+2];
        S = io[8+4+3];

    packed<double,4> adxs = packed<double,4>(a, d, x, s);
    packed<double,4> bcqr = packed<double,4>(b, c, q, r);
        packed<double,4> ADXS = packed<double,4>(A, D, X, S);
        packed<double,4> BCQR = packed<double,4>(B, C, Q, R);

    packed<double,4> eftu = add(adxs, bcqr);
        packed<double,4> EFTU = add(ADXS, BCQR);

    packed<double,4> ghvw = sub(adxs, bcqr);
        packed<double,4> GHVW = sub(ADXS, BCQR);

    ghvw = mulmod2(ghvw, packed<double,4>(W_0, IW_0, W_1, IW_1), N, Ninv);
        GHVW = mulmod2(GHVW, packed<double,4>(W_2, IW_2, W_3, IW_3), N, Ninv);

    packed<double,4> fhuw = unpackhi(eftu, ghvw);
    packed<double,4> egtv = unpacklo(eftu, ghvw);
        packed<double,4> FHUW = unpackhi(EFTU, GHVW);
        packed<double,4> EGTV = unpacklo(EFTU, GHVW);

    adxs = add(fhuw, egtv);
        ADXS = add(FHUW, EGTV);

    packed<double,4> cbrq = sub(fhuw, egtv);
        packed<double,4> CBRQ = sub(FHUW, EGTV);

    packed<double,4> cdrs = blend<0,1,0,1>(cbrq, adxs);
        packed<double,4> CDRS = blend<0,1,0,1>(CBRQ, ADXS);

    cdrs = mulmod2(cdrs, packed<double,4>(W2_0, W2_0, W2_1, W2_1), N, Ninv);
        CDRS = mulmod2(CDRS, packed<double,4>(W2_2, W2_2, W2_3, W2_3), N, Ninv);

    a = adxs[0];
    x = adxs[2];
    b = cbrq[1];
    q = cbrq[3];

    c = cdrs[0];
    d = cdrs[1];
    r = cdrs[2];
    s = cdrs[3];

    a = reduce_to_pm1n(a, n, ninv);
    x = reduce_to_pm1n(x, n, ninv);

        A = ADXS[0];
        X = ADXS[2];
        B = CBRQ[1];
        Q = CBRQ[3];

        C = CDRS[0];
        D = CDRS[1];
        R = CDRS[2];
        S = CDRS[3];

        A = reduce_to_pm1n(A, n, ninv);
        X = reduce_to_pm1n(X, n, ninv);

    io[0] = a;
    io[1] = b;
    io[2] = c;
    io[3] = d;
    io[4+0] = x;
    io[4+1] = q;
    io[4+2] = r;
    io[4+3] = s;

        io[8+0] = A;
        io[8+1] = B;
        io[8+2] = C;
        io[8+3] = D;
        io[8+4+0] = X;
        io[8+4+1] = Q;
        io[8+4+2] = R;
        io[8+4+3] = S;
}

template <>
inline void fftv2_ctx::ifft_basecase<4>(double* X, ulong j)
{
#if 0
    ifft_basecase<2>(X+4*0, 4*j+0);
    ifft_basecase<2>(X+4*1, 4*j+1);
    ifft_basecase<2>(X+4*2, 4*j+2);
    ifft_basecase<2>(X+4*3, 4*j+3);
    PD_RADIX_4_REVERSE_PARAM(j)
    PD_RADIX_4_REVERSE(X+4*0, X+4*1, X+4*2, X+4*3);
#else
    ifft_basecase2_4x(X+4*0, 4*j+0);
    PD_RADIX_4_REVERSE_PARAM(j)
    PD_RADIX_4_REVERSE(X+4*0, X+4*1, X+4*2, X+4*3);
#endif
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

