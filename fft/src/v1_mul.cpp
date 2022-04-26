// a = a*b*m
template<int vsize>
static void pd_mul(
    pd_fft_ctx<vsize>& ctx,
    packed<double,vsize>* Fa,
    packed<double,vsize>* Fb,
    const ulong m__[],
    ulong len)
{
    packed<ulong,vsize> m_;
    m_.load(m__);
    packed<double,vsize> m = convert_limited<packed<double,vsize>>(m_);
    packed<double,vsize> N = ctx.n;
    packed<double,vsize> Ninv = ctx.ninv;

    assert(len > 0);
    assert(len%4 == 0);
    ulong i = 0;
    do {
        packed<double,vsize> x0 = Fa[i+0], x1 = Fa[i+1], x2 = Fa[i+2], x3 = Fa[i+3];
        x0 = mulmod2(x0, m, N, Ninv);
        x1 = mulmod2(x1, m, N, Ninv);
        x2 = mulmod2(x2, m, N, Ninv);
        x3 = mulmod2(x3, m, N, Ninv);
        x0 = mulmod2(x0, Fb[i+0], N, Ninv);
        x1 = mulmod2(x1, Fb[i+1], N, Ninv);
        x2 = mulmod2(x2, Fb[i+2], N, Ninv);
        x3 = mulmod2(x3, Fb[i+3], N, Ninv);
        Fa[i+0] = x0;
        Fa[i+1] = x1;
        Fa[i+2] = x2;
        Fa[i+3] = x3;
    } while (i+=4, i < len);
}

template <ulong MN, ulong vsize>
inline void
crt(ulong res[MN], const packed<double,vsize>& a);

template <>
FORCE_INLINE inline void
crt<4,4>(ulong res[4], const packed<double,4>& a)
{
    packed<ulong,4> ai = convert_limited<packed<ulong,4>>(a);
    ulong a0 = ai[0];
    ulong a1 = ai[1];
    ulong a2 = ai[2];
    ulong a3 = ai[3];

#if USE_BITS == 51
    const ulong f[4]   = {0x001ff91000000001, 0xec853d7f59909fc0, 0x2427742e0a74cd09, 0x0000000000000ff2};
    const ulong f2[4]  = {0x003ff22000000002, 0xd90a7afeb3213f80, 0x484ee85c14e99a13, 0x0000000000001fe4};
    const ulong m0f[3] = {0x0017f98000000001, 0x348ca0bf980dc7c0, 0x0000000001fe606e};
    const ulong m1f[3] = {0x0017fa9800000001, 0xc14cc0bfa9885880, 0x0000000001fea642};
    const ulong m2f[3] = {0x0017fb8800000001, 0xba5340bfb8859780, 0x0000000001fee22c};
    const ulong m3f[3] = {0x0017fb9000000001, 0x3c589cbfb90587c0, 0x0000000001fee42c};
#else
    const ulong f[4]   = {0xf00fffc200000001, 0xbc9747effd180592, 0xf080592e1f40c894, 0x00000000000000ff};
    const ulong f2[4]  = {0xe01fff8400000002, 0x792e8fdffa300b25, 0xe100b25c3e819129, 0x00000000000001ff};
    const ulong m0f[3] = {0x000bffcec0000001, 0x0c7bef8b3e76031f, 0x00000000003ffcec};
    const ulong m1f[3] = {0xf00bffcf00000001, 0x0c57afdd2e780315, 0x00000000003ffcf0};
    const ulong m2f[3] = {0x000bffd340000001, 0x0a3bf3e2be9a028f, 0x00000000003ffd34};
    const ulong m3f[3] = {0xf00bffd500000001, 0x0987b504cea80261, 0x00000000003ffd50};
#endif
    ulong r0,r1,r2,r3, t0,t1,t2,t3;
    ulong u, v;

#if 0
    umul_ppmm(r[1], r[0], m0f[0], a[0]);
    umul_ppmm(r[3], r[2], m0f[2], a[0]);
    umul_ppmm(t[2], t[1], m0f[1], a[0]);

    umul_ppmm(u, v, m1f[0], a[1]);
    add_ssaaaa(r[1],r[0], r[1],r[0], u,v);
    umul_ppmm(u, v, m1f[2], a[1]);
    add_ssaaaa(r[3],r[2], r[3],r[2], u,v);
    umul_ppmm(u, v, m1f[1], a[1]);
    add_ssaaaa(t[2],t[1], t[2],t[1], u,v);

    umul_ppmm(u, v, m2f[0], a[2]);
    add_ssaaaa(r[1],r[0], r[1],r[0], u,v);
    umul_ppmm(u, v, m2f[2], a[2]);
    add_ssaaaa(r[3],r[2], r[3],r[2], u,v);
    umul_ppmm(u, v, m2f[1], a[2]);
    add_ssaaaa(t[2],t[1], t[2],t[1], u,v);

    umul_ppmm(u, v, m3f[0], a[3]);
    add_ssaaaa(r[1],r[0], r[1],r[0], u,v);
    umul_ppmm(u, v, m3f[2], a[3]);
    add_ssaaaa(r[3],r[2], r[3],r[2], u,v);
    umul_ppmm(u, v, m3f[1], a[3]);
    add_ssaaaa(t[2],t[1], t[2],t[1], u,v);
#else
    r0 = _mulx_ulong(m0f[0], a0, r1);
    r2 = _mulx_ulong(m0f[2], a0, r3);
    t1 = _mulx_ulong(m0f[1], a0, t2);

    v = _mulx_ulong(m1f[0], a1, u);
    add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m1f[2], a1, u);
    add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m1f[1], a1, u);
    add_ssaaaa(t2,t1, t2,t1, u,v);

    v = _mulx_ulong(m2f[0], a2, u);
    add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m2f[2], a2, u);
    add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m2f[1], a2, u);
    add_ssaaaa(t2,t1, t2,t1, u,v);

    v = _mulx_ulong(m3f[0], a3, u);
    add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m3f[2], a3, u);
    add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m3f[1], a3, u);
    add_ssaaaa(t2,t1, t2,t1, u,v);
#endif

    add_sssaaaaaa(r3,r2,r1, r3,r2,r1, 0,t2,t1);

    sub_ddddmmmmssss(t3,t2,t1,t0, r3,r2,r1,r0, f2[3],f2[2],f2[1],f2[0]);
    if ((slong)t3 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
        r3 = t3;
    }
    sub_ddddmmmmssss(t3,t2,t1,t0, r3,r2,r1,r0, f[3],f[2],f[1],f[0]);
    if ((slong)t3 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
        r3 = t3;
    }
    res[0] = r0;
    res[1] = r1;
    res[2] = r2;
    res[3] = r3;
}

template <>
FORCE_INLINE inline void
crt<3,4>(ulong res[3], const packed<double,4>& a)
{
    packed<ulong,4> ai = convert_limited<packed<ulong,4>>(a);
    ulong a0 = ai[0];
    ulong a1 = ai[1];
    ulong a2 = ai[2];
    ulong a3 = ai[3];

    const ulong f[3]   = {0x0002ff3300000001, 0x5e893a283e29bc7f, 0x3faea3b0884178cd};
    const ulong f2[3]  = {0x0005fe6600000002, 0xbd1274507c5378fe, 0x7f5d47611082f19a};
    const ulong m0f[3] = {0x0001ff5400000001, 0xda42aad93f17a653, 0x0000000000003fb6};
    const ulong m1f[3] = {0x0001ff6600000001, 0x5602aac93f299dd1, 0x0000000000003fbb};
    const ulong m2f[3] = {0x00027f6b00000001, 0x189bd4aeff001be7, 0x0000000000007f95};
    const ulong m3f[3] = {0x00027f7400000001, 0x15a80fd9ff1218f3, 0x0000000000007f9e};
    ulong r0, r1, r2, t0, t1, t2;
    ulong u, v;

    umul_ppmm(r1,r0, m0f[0], a0);
    umul_ppmm(t2,t1, m0f[1], a0);
    r2 =             m0f[2]* a0;

    umul_ppmm(u,v, m1f[0], a1); add_ssaaaa(r1,r0, r1,r0, u,v);
    umul_ppmm(u,v, m1f[1], a1); add_ssaaaa(t2,t1, t2,t1, u,v);
    r2 +=          m1f[2]* a1;

    umul_ppmm(u,v, m2f[0], a2); add_ssaaaa(r1,r0, r1,r0, u,v);
    umul_ppmm(u,v, m2f[1], a2); add_ssaaaa(t2,t1, t2,t1, u,v);
    r2 +=          m2f[2]* a2;

    umul_ppmm(u,v, m3f[0], a3); add_ssaaaa(r1,r0, r1,r0, u,v);
    umul_ppmm(u,v, m3f[1], a3); add_ssaaaa(t2,t1, t2,t1, u,v);
    r2 +=          m3f[2]* a3;

    add_ssaaaa(r2,r1, r2,r1, t2,t1);

    sub_dddmmmsss(t2,t1,t0, r2,r1,r0, f2[2],f2[1],f2[0]);
    if ((slong)t2 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
    }
    sub_dddmmmsss(t2,t1,t0, r2,r1,r0, f[2],f[1],f[0]);
    if ((slong)t2 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
    }
    res[0] = r0;
    res[1] = r1;
    res[2] = r2;
}


// x += crt << R
template <ulong N, ulong R, ulong vsize>
inline void
crt_shift_add(ulong x[N], const packed<double,vsize>& a);


template <>
FORCE_INLINE inline void
crt_shift_add<4,0,4>(ulong x[4], const packed<double,4>& a)
{
    packed<ulong,4> ai = convert_limited<packed<ulong,4>>(a);
    ulong a0 = ai[0];
    ulong a1 = ai[1];
    ulong a2 = ai[2];
    ulong a3 = ai[3];
    const ulong f[4]   = {0xf00fffc200000001, 0xbc9747effd180592, 0xf080592e1f40c894, 0x00000000000000ff};
    const ulong f2[4]  = {0xe01fff8400000002, 0x792e8fdffa300b25, 0xe100b25c3e819129, 0x00000000000001ff};
    const ulong m0f[3] = {0x000bffcec0000001, 0x0c7bef8b3e76031f, 0x00000000003ffcec};
    const ulong m1f[3] = {0xf00bffcf00000001, 0x0c57afdd2e780315, 0x00000000003ffcf0};
    const ulong m2f[3] = {0x000bffd340000001, 0x0a3bf3e2be9a028f, 0x00000000003ffd34};
    const ulong m3f[3] = {0xf00bffd500000001, 0x0987b504cea80261, 0x00000000003ffd50};
    ulong r0,r1,r2,r3, t0,t1,t2,t3;
    ulong u, v;

    r0 = _mulx_ulong(m0f[0], a0, r1);
    r2 = _mulx_ulong(m0f[2], a0, r3);
    t1 = _mulx_ulong(m0f[1], a0, t2);
    v = _mulx_ulong(m1f[0], a1, u); add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m1f[2], a1, u); add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m1f[1], a1, u); add_ssaaaa(t2,t1, t2,t1, u,v);
    v = _mulx_ulong(m2f[0], a2, u); add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m2f[2], a2, u); add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m2f[1], a2, u); add_ssaaaa(t2,t1, t2,t1, u,v);
    v = _mulx_ulong(m3f[0], a3, u); add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m3f[2], a3, u); add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m3f[1], a3, u); add_ssaaaa(t2,t1, t2,t1, u,v);
    add_sssaaaaaa(r3,r2,r1, r3,r2,r1, 0,t2,t1);

    sub_ddddmmmmssss(t3,t2,t1,t0, r3,r2,r1,r0, f2[3],f2[2],f2[1],f2[0]);
    if ((slong)t3 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
        r3 = t3;
    }
    sub_ddddmmmmssss(t3,t2,t1,t0, r3,r2,r1,r0, f[3],f[2],f[1],f[0]);
    if ((slong)t3 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
        r3 = t3;
    }
    add_ssssaaaaaaaa(x[3],x[2],x[1],x[0], x[3],x[2],x[1],x[0], r3,r2,r1,r0);
}

template <>
FORCE_INLINE inline void
crt_shift_add<4,32,4>(ulong x[4], const packed<double,4>& a)
{
    packed<ulong,4> ai = convert_limited<packed<ulong,4>>(a);
    ulong a0 = ai[0];
    ulong a1 = ai[1];
    ulong a2 = ai[2];
    ulong a3 = ai[3];
    const ulong f[4]   = {0x0000000100000000, 0xfd180592f00fffc2, 0x1f40c894bc9747ef, 0x000000fff080592e};
    const ulong f2[4]  = {0x0000000200000000, 0xfa300b25e01fff84, 0x3e819129792e8fdf, 0x000001ffe100b25c};
    const ulong m0f[3] = {0xc000000100000000, 0x3e76031f000bffce, 0x003ffcec0c7bef8b};
    const ulong m1f[3] = {0x0000000100000000, 0x2e780315f00bffcf, 0x003ffcf00c57afdd};
    const ulong m2f[3] = {0x4000000100000000, 0xbe9a028f000bffd3, 0x003ffd340a3bf3e2};
    const ulong m3f[3] = {0x0000000100000000, 0xcea80261f00bffd5, 0x003ffd500987b504};
    ulong r0,r1,r2,r3, t0,t1,t2,t3;
    ulong u, v;

    r0 = _mulx_ulong(m0f[0], a0, r1);
    r2 = _mulx_ulong(m0f[2], a0, r3);
    t1 = _mulx_ulong(m0f[1], a0, t2);

    v = _mulx_ulong(m1f[0], a1, u); add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m1f[2], a1, u); add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m1f[1], a1, u); add_ssaaaa(t2,t1, t2,t1, u,v);

    v = _mulx_ulong(m2f[0], a2, u); add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m2f[2], a2, u); add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m2f[1], a2, u); add_ssaaaa(t2,t1, t2,t1, u,v);

    v = _mulx_ulong(m3f[0], a3, u); add_ssaaaa(r1,r0, r1,r0, u,v);
    v = _mulx_ulong(m3f[2], a3, u); add_ssaaaa(r3,r2, r3,r2, u,v);
    v = _mulx_ulong(m3f[1], a3, u); add_ssaaaa(t2,t1, t2,t1, u,v);

    add_sssaaaaaa(r3,r2,r1, r3,r2,r1, 0,t2,t1);

    sub_ddddmmmmssss(t3,t2,t1,t0, r3,r2,r1,r0, f2[3],f2[2],f2[1],f2[0]);
    if ((slong)t3 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
        r3 = t3;
    }
    sub_ddddmmmmssss(t3,t2,t1,t0, r3,r2,r1,r0, f[3],f[2],f[1],f[0]);
    if ((slong)t3 >= 0)
    {
        r0 = t0;
        r1 = t1;
        r2 = t2;
        r3 = t3;
    }
    add_ssssaaaaaaaa(x[3],x[2],x[1],x[0], x[3],x[2],x[1],x[0], r3,r2,r1,r0);
}



/********************** crt ***************************
Compile time constants will be denoted with caps.

Set B = bits. The output limbs are x[0], x[1], ....

Assume that the output of the crt fits in M-1 bits, i.e. is less than 2^(M-1).
Further suppose that M >= 128

The i^th coefficient corresponds to output bits [B*i, B*i+M) or output limbs
    x[fdiv(B*i,64)], ..., x[cdiv(B*i+M,64)-1]
the number of such limbs is cdiv(B*i+M,64)-fdiv(B*i,64) which is <= cdiv(M,64)+1.

Let G be any common divisor of 64 and B, let L = B/G and F = 64/G.

The index i can be split as j*F + K, where we consider j fixed for now
and let K range over 0 <= K < F. The output limbs are x[j] for
    fdiv(B*(j*F+K),64) <= j < cdiv(B*(j*F+K)+M,64)
or
    L*j + fdiv(B*K,64) <= j < L*j + cdiv(B*K+M,64)

Therefore, there are N = cdiv(B*(F-1)+M,64) output limbs involved for any fixed j
and K ranging over [0,F). The lower L of these coefficients can be written
out, and the rest carry into the next chunk with incremented iq.
*/
template <bool inbounds, ulong B, ulong G, ulong M, ulong N, ulong vsize, ulong K>
FORCE_INLINE inline bool _from_inner_loop(
    ulong* X, ulong Xn,
    ulong x[N],
    packed<double,vsize>* Y, ulong Yn,
    ulong j,
    packed<double,vsize> n,
    packed<double,vsize> ninv)
{
    constexpr ulong F = 64/G;
    constexpr ulong L = B/G;
    constexpr ulong MN = cdiv(M,64);
    constexpr ulong MR = 1 + (M-1)%64; // number of bits in M[MN-1]
    constexpr ulong R = (B*K)%64;

    static_assert(L <= N);

    if (K >= F)
        return false;

    // shift output limbs for next chunk
#if 0
    if (K == 0)
    {
        for (ulong i = 0; i < N; i++)
            x[i] = i + L < N ? x[i + L] : 0;
    }
#else
    if (K == 0)
    {
        for (ulong i = 0; i < N-L; i++)
            x[i] = x[i + L];

        for (ulong i = N-L; i < N; i++)
            x[i] = 0;
    }
    else
    {
        constexpr ulong K1 = K-1;
        constexpr ulong R1 = (B*K1)%64;
        for (ulong i = ((B*K1)/64) + MN+(MR+R1>64); i < N; i++)
            x[i] = 0;
    }
#endif

    // read input
    if (inbounds || F*j+K < Yn)
    {
#if 1
        ulong z[MN+1]; // only MN limbs defined by crt op
        crt<MN,vsize>(z, reduce_to_0n(Y[F*j+K], n, ninv));
        // shift and add to x[fdiv(B*K,64)], ..., X[cdiv(B*K+M,64)-1]
        if (R > 0)
        {
            z[MN] = z[MN-1] >> ((64-R)%64);
            for (ulong i = MN-1; i > 0; i--)
                z[i] = (z[i] << R) | (z[i-1] >> ((64-R)%64));
            z[0] = z[0] << R;
        }
        _multi_add<MN+(MR+R>64)>(x+((B*K)/64), z);
#else
        crt_shift_add<MN+(MR+R>64),R,vsize>(x+((B*K)/64), reduce_to_0n(Y[F*j+K], n, ninv));
#endif
    }

    // outputs for this value of K
    for (ulong i = (B*K)/64; i < (B*(K+1))/64; i++)
    {
        if (!inbounds && L*j + i >= Xn)
            return true;
        X[L*j + i] = x[i];
    }

    return _from_inner_loop<inbounds, B, G, M, N, vsize, std::min(K+1,F)>(X, Xn, x, Y, Yn, j, n, ninv);
}

template <ulong B, ulong G, ulong M, ulong vsize>
void mpn64_from_pd_fft(
    ulong* X, const ulong Xn,
    packed<double,vsize>* Y, const ulong Yn,
    pd_fft_ctx<vsize>& ctx)
{
    constexpr ulong F = 64/G;
    constexpr ulong L = B/G;
    constexpr ulong N = cdiv(B*(F-1)+M,64);
    packed<double,vsize> n = ctx.n;
    packed<double,vsize> ninv = ctx.ninv;

    ulong x[N];
    for (ulong i = 0; i < N; i++)
        x[i] = 0;

    ulong j = 0;
    ulong easy_end = std::min(Yn/F, Xn/L);

    for (; j < easy_end; j++)
        _from_inner_loop<true, B, G, M, N, vsize, 0>(X, Xn, x, Y, Yn, j, n, ninv);

    // The previous loop was written with the assumption that all reads/writes
    // were in bounds. Now be more careful at the end.
    for (; L*j < Xn; j++)
       if (_from_inner_loop<false, B, G, M, N, vsize, 0>(X, Xn, x, Y, Yn, j, n, ninv))
            break;
}



/*********************** mod *********************************


*/
template <ulong bits, ulong G, int vsize>
void mpn32_to_pd_fft(
    packed<double,vsize>* y,
    const ulong* x_, ulong xn_,
    ulong m,  // = cdiv(32*xn, bits)
    ulong itrunc,
    pd_fft_ctx<vsize>& ctx,
    const packed<double,vsize>* two_pow)
{
    const uint32_t* x = reinterpret_cast<const uint32_t*>(x_);
    ulong xn = 2*xn_;
    ulong L = bits/G;
    ulong F = 32/G;

    assert(m == cdiv(32*xn, bits));
    assert(L*G == bits);
    assert(F*G == 32);

    uint32_t tx[L];
    packed<double,vsize> p = ctx.n;
    packed<double,vsize> pinv = ctx.ninv;

    ulong J;
    for (J = 0; L*J + L <= xn; J++)
    {
        assert(F*J + F <= m);
        /*
            for each 0 <= K < F, y[F*J + K] should come from bits starting from
                bits*(F*J + K) = L*32*J + bits*K
            and use x[L*J + 0], ..., x[L*J + L-1]

            this requires L*32*J + bits*F <= 32*(L*J + L) which is an equality
            because bits*F == 32*L
        */
#define CODE(mult, K, xx) \
{\
    if (mult == 0 && F*J + K >= m)\
        goto done;\
    ulong k = mult*L*J + (bits*K)/32;\
    ulong j = (bits*K)%32;\
    ulong nchunks = 1;\
    packed<double,4> ak = double(xx[k] >> j);\
    packed<double,4> Y = ak;\
    k++;\
    j = 32 - j;\
    while (j + 32 <= bits)\
    {\
        ak = double(xx[k]);\
        Y = add(Y, mulmod2(ak, two_pow[j], p, pinv));\
        nchunks++;\
        k++;\
        j += 32;\
    }\
    if ((bits-j) != 0)\
    {\
        ak = double(xx[k] << (32-(bits-j)));\
        Y = add(Y, mulmod2(ak, two_pow[bits-32], p, pinv));\
        nchunks++;\
    }\
    if (nchunks > 3)\
        Y = reduce_to_pm1n(Y, p, pinv);\
    y[F*J + K] = Y;\
}
        if (F == 1) {
            CODE(1,0,x);
        } else if (F == 2) {
            CODE(1,0,x);CODE(1,1,x);
        } else if (F == 4) {
            CODE(1,0,x);CODE(1,1,x);CODE(1,2,x);CODE(1,3,x);
        } else if (F == 8) {
            CODE(1,0,x);CODE(1,1,x);CODE(1,2,x);CODE(1,3,x);CODE(1,4,x);CODE(1,5,x);CODE(1,6,x);CODE(1,7,x);
        } else if (F == 16) {
            CODE(1,0,x);CODE(1,1,x);CODE(1,2,x);CODE(1,3,x);CODE(1,4,x);CODE(1,5,x);CODE(1,6,x);CODE(1,7,x);
            CODE(1,8,x);CODE(1,9,x);CODE(1,10,x);CODE(1,11,x);CODE(1,12,x);CODE(1,13,x);CODE(1,14,x);CODE(1,15,x);
        } else {
            assert(false);
        }
    }

    /*
        at this point  0 <= xn - l*j < l, and it is possible to know m from
        the value of xn - l*j. Have to produce
            y[F*J + 0], ..., y[min(F*J + F -1, m - 1)]
    */
    assert(xn - L*J < L);

    for (ulong M = 0; M < L; M++)
        tx[M] = (L*J + M < xn) ? x[L*J + M] : 0;

    if (F == 1) {
        CODE(0,0,tx);
    } else if (F == 2) {
        CODE(0,0,tx);CODE(0,1,tx);
    } else if (F == 4) {
        CODE(0,0,tx);CODE(0,1,tx);CODE(0,2,tx);CODE(0,3,tx);
    } else if (F == 8) {
        CODE(0,0,tx);CODE(0,1,tx);CODE(0,2,tx);CODE(0,3,tx);CODE(0,4,tx);CODE(0,5,tx);CODE(0,6,tx);CODE(0,7,tx);
    } else if (F == 16) {
        CODE(0,0,tx);CODE(0,1,tx);CODE(0,2,tx);CODE(0,3,tx);CODE(0,4,tx);CODE(0,5,tx);CODE(0,6,tx);CODE(0,7,tx);
        CODE(0,8,tx);CODE(0,9,tx);CODE(0,10,tx);CODE(0,11,tx);CODE(0,12,tx);CODE(0,13,tx);CODE(0,14,tx);CODE(0,15,tx);
    } else {
        assert(false);
    }

#undef CODE

done:

    while (m < itrunc)
        y[m++] = 0;
}

void mpn_fft_ctx::my_mpn_mul(
    ulong* c,
    const ulong* a, ulong an,
    const ulong* b, ulong bn)
{
    assert(an > 0);
    assert(bn > 0);

    if (an < bn)
    {
        std::swap(a, b);
        std::swap(an, bn);
    }

    ulong bits;
    void (*to_pd_fft)(packed<double,4>*, const ulong*, ulong, ulong, ulong, pd_fft_ctx<4>&, const packed<double,4>*);
    void (*from_pd_fft)(ulong*, ulong, packed<double,4>*, ulong, pd_fft_ctx<4>&);

#if USE_BITS == 51
    ulong m[4] = {0x00058329c9f6b64e, 0x0002def60b4c2386,
                  0x00032d730f79897f, 0x000345fcdf05f1e0};
    if (bn <= 6123)
    {
        to_pd_fft =     &mpn32_to_pd_fft<96,16,4>;
        from_pd_fft = &mpn64_from_pd_fft<96,32,217,4>;
        bits = 96;
    }
    else
    {
        to_pd_fft =     &mpn32_to_pd_fft<64,32,4>;
        from_pd_fft = &mpn64_from_pd_fft<64,32,217,4>;
        bits = 64;
    }
#elif USE_BITS == 50
    ulong m[4] = {0x000209495a71dc4b,
                  0x00014238d1ef2a4a,
                  0x0003eae515fb9ecc,
                  0x00014160108dcea6};
    if (bn <= 382)
    {
        to_pd_fft =     &mpn32_to_pd_fft<96,16,4>;
        from_pd_fft = &mpn64_from_pd_fft<96,32,209,4>;
        bits = 96;
    }
    else if (bn <= 94185)
    {
        to_pd_fft =     &mpn32_to_pd_fft<92,4,4>;
        from_pd_fft = &mpn64_from_pd_fft<92,4,209,4>;
        bits = 92;
    }
    else if (bn <= 23063216)
    {
        to_pd_fft =     &mpn32_to_pd_fft<88,8,4>;
        from_pd_fft = &mpn64_from_pd_fft<88,8,209,4>;
        bits = 88;
    }
    else
    {
        to_pd_fft =     &mpn32_to_pd_fft<64,32,4>;
        from_pd_fft = &mpn64_from_pd_fft<64,32,209,4>;
        bits = 64;
    }
#else
    ulong m[4] = {0x0000d8b548614dc7,
                  0x00003171aecda7bf,
                  0x00002eee585c016c,
                  0x00003fd1d2a52fa3};
    if (false&&bn <= 1431)
    {
        to_pd_fft =     &mpn32_to_pd_fft<90,2,4>;
        from_pd_fft = &mpn64_from_pd_fft<90,2,192,4>;
        bits = 90;
    }
    else if (bn <= 22415)
    {
        to_pd_fft =     &mpn32_to_pd_fft<88,8,4>;
        from_pd_fft = &mpn64_from_pd_fft<88,8,192,4>;
        bits = 88;
    }
    else
    {
        to_pd_fft =     &mpn32_to_pd_fft<64,32,4>;
        from_pd_fft = &mpn64_from_pd_fft<64,32,192,4>;
        bits = 64;
    }
#endif

    ulong a_len = cdiv(64*an, bits);
    ulong b_len = cdiv(64*bn, bits);
    ulong z_len = a_len + b_len - 1;
    ulong a_trunc = round_up(a_len, 32);
    ulong b_trunc = round_up(b_len, 32);
    ulong z_trunc = round_up(z_len, 32);
    ulong k = std::max(ulong(2), clog2(z_trunc));

    ctx.fit_wrevtab(k);
    const packed<double,4>* two_powers = two_power_table(bits + 5);

    for (int i = 0; i < 4; i++)
    {
        ulong thi = m[i] >> (FLINT_BITS - k);
        ulong tlo = m[i] << (k);
        NMOD_RED2(tlo, thi, tlo, ctx.mod[i]);
        m[i] = nmod_inv(tlo, ctx.mod[i]);
    }

    ulong n = pow2(k);
    buffer.resize(2*n);
    packed<double,4>* Fa = buffer.data();
    packed<double,4>* Fb = buffer.data() + n;

    to_pd_fft(Fa, a, an, a_len, a_trunc, ctx, two_powers);
    ctx.fft_trunc(Fa, k, 0, a_trunc, z_trunc);
    to_pd_fft(Fb, b, bn, b_len, b_trunc, ctx, two_powers);
    ctx.fft_trunc(Fb, k, 0, b_trunc, z_trunc);
    pd_mul<4>(ctx, Fa, Fb, m, z_trunc);
    ctx.ifft_trunc(Fa, k, 0, z_trunc, z_trunc, false);
    from_pd_fft(c, an + bn, Fa, z_len, ctx);
}

