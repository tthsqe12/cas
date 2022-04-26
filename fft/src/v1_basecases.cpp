static void _missing_helper(int l, int ntrunc, int ztrunc, bool f)
{
    std::cout << std::endl << "function l = " << l << ", n = " << ntrunc << ", z = " << ztrunc << ", f = " << f;
    if (1 <= ztrunc && ztrunc <= l && ntrunc <= ztrunc && 1 <= ntrunc+f && ntrunc+f <= l)
        std::cout << " is not implmented" << std::endl;
    else
        std::cout << " does not exist and should not be called" << std::endl;
    std::abort();
}

#define PD packed<double,vsize>

template <int vsize, int itrunc, int otrunc>
static void fft_moth_trunc_loop(
    pd_fft_ctx<vsize>& Q,
    PD* x,
    ulong s,
    ulong j,
    ulong len)
{
    assert(len > 0);
    assert(len%2 == 0);
    assert(2*j+1 < Q.wrevtab.size());
    PD n = Q.n;
    PD ninv = Q.ninv;
    PD w2 = Q.wrevtab.data()[j];
    PD  w = Q.wrevtab.data()[2*j];
    PD iw = Q.wrevtab.data()[2*j+1];
    
    ulong i = 0;
    do {
        PD x0(0.0), x1(0.0), x2(0.0), x3(0.0), y0, y1, y2, y3;
        PD u0(0.0), u1(0.0), u2(0.0), u3(0.0), v0, v1, v2, v3;
        if (0 < itrunc) x0 = x[0*s+i];
            if (0 < itrunc) u0 = x[0*s+i+1];
        if (1 < itrunc) x1 = x[1*s+i];
            if (1 < itrunc) u1 = x[1*s+i+1];
        if (2 < itrunc) x2 = x[2*s+i];
            if (2 < itrunc) u2 = x[2*s+i+1];
        if (3 < itrunc) x3 = x[3*s+i];
            if (3 < itrunc) u3 = x[3*s+i+1];
        if (0 < itrunc) x0 = reduce_to_pm1n(x0, n, ninv);
            if (0 < itrunc)u0 = reduce_to_pm1n(u0, n, ninv);
        if (2 < itrunc) x2 = mulmod2(x2, w2, n, ninv);
            if (2 < itrunc) u2 = mulmod2(u2, w2, n, ninv);
        if (3 < itrunc) x3 = mulmod2(x3, w2, n, ninv);
            if (3 < itrunc) u3 = mulmod2(u3, w2, n, ninv);
        y0 = (2 < itrunc) ? add(x0, x2) : x0;
            v0 = (2 < itrunc) ? add(u0, u2) : u0;
        y1 = (3 < itrunc) ? add(x1, x3) : x1;
            v1 = (3 < itrunc) ? add(u1, u3) : u1;
        y2 = (2 < itrunc) ? sub(x0, x2) : x0;
            v2 = (2 < itrunc) ? sub(u0, u2) : u0;
        y3 = (3 < itrunc) ? sub(x1, x3) : x1;
            v3 = (3 < itrunc) ? sub(u1, u3) : u1;
        y1 = mulmod2(y1, w, n, ninv);
            v1 = mulmod2(v1, w, n, ninv);
        y3 = mulmod2(y3, iw, n, ninv);
            v3 = mulmod2(v3, iw, n, ninv);
        x0 = add(y0, y1);
            u0 = add(v0, v1);
        x1 = sub(y0, y1);
            u1 = sub(v0, v1);
        x2 = add(y2, y3);
            u2 = add(v2, v3);
        x3 = sub(y2, y3);
            u3 = sub(v2, v3);
        if (0 < otrunc) x[0*s+i] = x0;
            if (0 < otrunc) x[0*s+i+1] = u0;
        if (1 < otrunc) x[1*s+i] = x1;
            if (1 < otrunc) x[1*s+i+1] = u1;
        if (2 < otrunc) x[2*s+i] = x2;
            if (2 < otrunc) x[2*s+i+1] = u2;
        if (3 < otrunc) x[3*s+i] = x3;
            if (3 < otrunc) x[3*s+i+1] = u3;
    } while (i += 2, i < len);
}

template<int vsize, int ztrunc, int ntrunc, bool f>
static void radix_2_moth_inv_trunc_loop(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    _missing_helper(2, ntrunc, ztrunc, f);
}

#define vsize 4

template<> void radix_2_moth_inv_trunc_loop<vsize,2,2,false>(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    assert(len > 0);
    PD n = Q.n, ninv = Q.ninv;
    PD w = Q.w2revinv(j);
    for (ulong i = 0; i < len; i++)
    {
        // {x0, x1} = {2*x0 - w*x1, x0 - w*x1}
        PD u0, u1;
        u0 = add(x0[i], x1[i]);
        u1 = sub(x0[i], x1[i]);
        x0[i] = reduce_to_pm1n(u0, n, ninv);
        x1[i] = mulmod2(u1, w, n, ninv);
    }
}

template<> void radix_2_moth_inv_trunc_loop<vsize,2,1,true>(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    assert(len > 0);
    PD n = Q.n, ninv = Q.ninv;
    PD w = Q.w2rev(j);
    PD c = 2.0;
    for (ulong i = 0; i < len; i++)
    {
        // {x0, x1} = {2*x0 - w*x1, x0 - w*x1}
        PD u0, u1;
        u0 = x0[i];
        u1 = x1[i];
        u0 = reduce_to_pm1n(u0, n, ninv);
        u1 = mulmod2(u1, w, n, ninv);
        x0[i] = fmsub(c, u0, u1);
        x1[i] = sub(u0, u1);
    }
}

template<> void radix_2_moth_inv_trunc_loop<vsize,2,1,false>(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    assert(len > 0);
    PD n = Q.n, ninv = Q.ninv;
    PD w = Q.w2rev(j);
    PD c = 2.0;
    for (ulong i = 0; i < len; i++)
    {
        // {x0} = {2*x0 - w*x1}
        PD u0, u1;
        u0 = x0[i];
        u1 = x1[i];
        u0 = reduce_to_pm1n(u0, n, ninv);
        u1 = mulmod2(u1, w, n, ninv);
        x0[i] = fmsub(c, u0, u1);
    }
}

template<> void radix_2_moth_inv_trunc_loop<vsize,2,0,true>(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    assert(len > 0);
    PD n = Q.n, ninv = Q.ninv;
    PD w = Q.w2rev(j);
    PD c = fnmadd(0.5, n, 0.5);
    for (ulong i = 0; i < len; i++)
    {
        // {x0} = {(x0 + w*x1)/2}
        PD u0, u1;
        c = fnmadd(0.5, n, 0.5);
        u0 = x0[i];
        u1 = x1[i];
        u1 = mulmod2(u1, w, n, ninv);
        u0 = mulmod2(add(u0, u1), c, n, ninv);
        x0[i] = u0;
    }
}

template<> void radix_2_moth_inv_trunc_loop<vsize,1,1,true>(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    assert(len > 0);
    PD n = Q.n, ninv = Q.ninv;
    PD w = Q.w2rev(j);
    PD c = 2.0;
    for (ulong i = 0; i < len; i++)
    {
        // {x0, x1} = {2*x0, x0}
        PD u0, u1;
        u0 = x0[i];
        u0 = reduce_to_pm1n(u0, n, ninv);
        x0[i] = add(u0, u0);
        x1[i] = u0;
    }
}

template<> void radix_2_moth_inv_trunc_loop<vsize,1,1,false>(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    assert(len > 0);
    PD n = Q.n, ninv = Q.ninv;
    PD w = Q.w2rev(j);
    PD c = 2.0;
    for (ulong i = 0; i < len; i++)
    {
        // {x0} = {2*x0}
        PD u0;
        u0 = x0[i];
        u0 = reduce_to_pm1n(u0, n, ninv);
        x0[i] = add(u0, u0);
    }
}

template<> void radix_2_moth_inv_trunc_loop<vsize,1,0,true>(
    pd_fft_ctx<vsize>& Q,
    PD* x0, PD* x1, ulong len,
    ulong j)
{
    assert(len > 0);
    PD n = Q.n, ninv = Q.ninv;
    PD w = Q.w2rev(j);
    PD c = fnmadd(0.5, n, 0.5);
    for (ulong i = 0; i < len; i++)
    {
        // {x0} = {x0/2}
        PD u0;
        u0 = x0[i];
        u0 = mulmod2(u0, c, n, ninv);
        x0[i] = u0;
    }
}

#undef vsize

template<int vsize, int ntrunc, int ztrunc, bool f>
static void radix_4_moth_inv_trunc_loop(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    _missing_helper(4, ntrunc, ztrunc, f);
}

#define vsize 4

/*
k = 2, n = 4, z = 4, f = false
[     1         1         1         1]
[  1//w     -1//w     -r//w      r//w]
[1//w^2    1//w^2   -1//w^2   -1//w^2]
[1//w^3   -1//w^3    r//w^3   -r//w^3]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,4,4,false>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W  = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -w2s[0] : w2s[jm];
    PD IW = UNLIKELY(j == 0) ?  w2s[1] : w2s[2*jm];
    ulong i = 0; do {
        packed<double,4> x0, x1, x2, x3, y0, y1, y2, y3;
        packed<double,4> u0, u1, u2, u3, v0, v1, v2, v3;
        x0 = X0[i+0];
            u0 = X0[i+1];
        x1 = X1[i+0];
            u1 = X1[i+1];
        x2 = X2[i+0];
            u2 = X2[i+1];
        x3 = X3[i+0];
            u3 = X3[i+1];
        y0 = add(x0, x1);
            v0 = add(u0, u1);
        y1 = add(x2, x3);
            v1 = add(u2, u3);
        y2 = sub(x0, x1);
            v2 = sub(u0, u1);
        y3 = sub(x3, x2);
            v3 = sub(u3, u2);
        y2 = mulmod2(y2, W, N, Ninv);
            v2 = mulmod2(v2, W, N, Ninv);
        y3 = mulmod2(y3, IW, N, Ninv);
            v3 = mulmod2(v3, IW, N, Ninv);
        x0 = add(y0, y1);
            u0 = add(v0, v1);
        x1 = sub(y3, y2);
            u1 = sub(v3, v2);
        X1[i+0] = x1;
            X1[i+1] = u1;
        x2 = sub(y1, y0);
            u2 = sub(v1, v0);
        x3 = add(y3, y2);
            u3 = add(v3, v2);
        x0 = reduce_to_pm1n(x0, N, Ninv);
            u0 = reduce_to_pm1n(u0, N, Ninv);
        x2 = mulmod2(x2, W2, N, Ninv);
            u2 = mulmod2(u2, W2, N, Ninv);
        x3 = mulmod2(x3, W2, N, Ninv);
            u3 = mulmod2(u3, W2, N, Ninv);
        X0[i+0] = x0;
            X0[i+1] = u0;
        X2[i+0] = x2;
            X2[i+1] = u2;
        X3[i+0] = x3;
            X3[i+1] = u3;
    } while (i += 2, i < len);
}

/*
k = 2, n = 3, z = 4, f = true
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
[          -r               r         1   r*w^3]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,3,4,true>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W  = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -w2s[0] : w2s[jm];
    PD f0 = w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2
    PD fr = -w2s[2*j+1];     // -r*w
    PD fq = -w2s[j];         // -w^2
    PD fp = reduce_pm1n_to_pmhn(mulmod2(fr, fq, N, Ninv), N);   // r*w^3
    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, v0, p0,q0,r0, a1, b1, c1, d1, u1, v1;
        a0 = X0[i+0];
        b0 = X1[i+0];
        c0 = X2[i+0];
        d0 = X3[i+0];
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
        X0[i+0] = a0+p0;
        X1[i+0] = b0+q0;
        X2[i+0] = c0+r0;
        X3[i+0] = d0+p0;
    } while (i += 1, i < len);
}

/*
k = 2, n = 3, z = 4, f = false
[      -r + 1           r + 1         2   r*w^3]
[        2//w           -2//w         0    -w^2]
[(r + 1)//w^2   (-r + 1)//w^2   -2//w^2    -r*w]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,3,4,false>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W  = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -w2s[0] : w2s[jm];
    PD f0 = w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2
    PD fr = -w2s[2*j+1];     // -r*w
    PD fq = -w2s[j];         // -w^2
    PD fp = reduce_pm1n_to_pmhn(mulmod2(fr, fq, N, Ninv), N);   // r*w^3
    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, v0, p0,q0,r0, a1, b1, c1, d1, u1, v1;
        a0 = X0[i+0];
        b0 = X1[i+0];
        c0 = X2[i+0];
        d0 = X3[i+0];
        u0 = add(a0, b0);
        v0 = sub(a0, b0);
        p0 = mulmod2(d0, fp, N, Ninv);
        q0 = mulmod2(d0, fq, N, Ninv);
        r0 = mulmod2(d0, fr, N, Ninv);
        c0 = reduce_to_pm1n(c0, N, Ninv);
        u0 = reduce_to_pm1n(u0, N, Ninv);
        b0 = mulmod2(v0, f1, N, Ninv);
        v0 = mulmod2(v0, f0, N, Ninv);
        c0 = fmsub(f2, c0, v0);
        a0 = add(c0, u0);
        c0 = sub(c0, u0);                
        c0 = mulmod2(c0, f3, N, Ninv);
        X0[i+0] = a0+p0;
        X1[i+0] = b0+q0;
        X2[i+0] = c0+r0;
    } while (i += 1, i < len);
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
template<> void radix_4_moth_inv_trunc_loop<vsize,3,3,true>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W  = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -w2s[0] : w2s[jm];
    PD f0 = w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2

    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, v0,  a1, b1, c1, d1, u1, v1;
        a0 = X0[i+0];
            a1 = X0[i+1];
        b0 = X1[i+0];
            b1 = X1[i+1];
        c0 = X2[i+0];
            c1 = X2[i+1];
        v0 = sub(a0, b0);
            v1 = sub(a1, b1);
        X1[i+0] = mulmod2(v0, f1, N, Ninv);
            X1[i+1] = mulmod2(v1, f1, N, Ninv);
        c0 = reduce_to_pm1n(c0, N, Ninv);
            c1 = reduce_to_pm1n(c1, N, Ninv);
        v0 = mulmod2(v0, f0, N, Ninv);
            v1 = mulmod2(v1, f0, N, Ninv);
        X3[i+0] = sub(c0, v0);
            X3[i+1] = sub(c1, v1);
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
        X0[i+0] = a0;
            X0[i+1] = a1;
        X2[i+0] = c0;
            X2[i+1] = c1;
    } while (i += 2, i < len);
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
template<> void radix_4_moth_inv_trunc_loop<vsize,3,3,false>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W  = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD W2 = UNLIKELY(j == 0) ? -w2s[0] : w2s[jm];
    PD f0 = w2s[1];                         // r
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = 2.0;
    PD f3 = W2;                             // -w^-2

    ulong i = 0; do {
        PD a0, b0, c0, u0, v0,  a1, b1, c1, u1, v1;
        a0 = X0[i+0];
            a1 = X0[i+1];
        b0 = X1[i+0];
            b1 = X1[i+1];
        c0 = X2[i+0];
            c1 = X2[i+1];
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
        X1[i+0] = b0;
            X1[i+1] = b1;
        X0[i+0] = a0;
            X0[i+1] = a1;
        X2[i+0] = c0;
            X2[i+1] = c1;
    } while (i += 2, i < len);
}

/*
k = 2, n = 2, z = 4, f = false
[   2       2   -w^2      0]
[2//w   -2//w      0   -w^2]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,2,4,false>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD mwi = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD w2 = w2s[j];
    PD twowi = reduce_pm1n_to_pmhn(-2.0*mwi, N);

    ulong i = 0; do {
        // BAD
        PD a0, b0, c0, d0, u0, v0, a1, b1, c1, d1, u1, v1;
        a0 = X0[i+0];
            a1 = X0[i+1];
        b0 = X1[i+0];
            b1 = X1[i+1];
        c0 = X2[i+0];
            c1 = X2[i+1];
        d0 = X3[i+0];
            d1 = X3[i+1];
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
        X0[i+0] = u0;
            X0[i+1] = u1;
        X1[i+0] = v0;
            X1[i+1] = v1;
    } while (i += 2, i < len);
}

/*
k = 2, n = 2, z = 4, f = true
[            2                2        -w^2             0]
[         2//w            -2//w           0          -w^2]
[1//2*r + 1//2   -1//2*r + 1//2   -1//2*w^2   -1//2*r*w^3]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,2,4,true>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD f0 = 2.0;
    PD f1 = reduce_pm1n_to_pmhn(-2.0*W, N); // 2*w^-1
    PD f2 = fnmadd(0.5, N, 0.5);    // 1//2
    PD f3 = w2s[1];   // r
    PD f4 = w2s[j];     // w^2
    PD f5 = reduce_pm1n_to_pmhn(mulmod2(f4, w2s[2*j+1], N, Ninv), N);   // r*w^3

    ulong i = 0; do {
        PD a0,b0,u0,v0,s0,t0,g0,h0,p0,q0,r0, u1, v1, s1, t1;
        u0 = X0[i+0];
        v0 = X1[i+0];
        a0 = X2[i+0];
        b0 = X3[i+0];
        p0 = mulmod2(a0, f4, N, Ninv);
        q0 = mulmod2(b0, f4, N, Ninv);
        r0 = mulmod2(b0, f5, N, Ninv);
        s0 = reduce_to_pm1n(add(u0, v0), N, Ninv);
        t0 = sub(u0, v0);
        g0 = mulmod2(s0, f0, N, Ninv);
        h0 = mulmod2(t0, f1, N, Ninv);
        t0 = mulmod2(t0, f3, N, Ninv);
        X0[i+0] = g0-p0;
        X1[i+0] = h0-q0;
        X2[i+0] = mulmod2((s0+t0)-(p0+r0), f2, N, Ninv);
    } while (i += 1, i < len);
}

/*
k = 2, n = 2, z = 2, f = true
[            2                2]
[         2//w            -2//w]
[1//2*r + 1//2   -1//2*r + 1//2]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,2,2,true>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD c1 = reduce_pm1n_to_pmhn(-2.0*W, N);
    PD c2 = fnmadd(0.5, N, 0.5);
    PD c3 = reduce_pm1n_to_pmhn(mulmod2(w2s[1], c2, N, Ninv), N);

    ulong i = 0; do {
        PD u0, v0, s0, t0, u1, v1, s1, t1;
        u0 = X0[i+0];
            u1 = X0[i+1];
        v0 = X1[i+0];
            v1 = X1[i+1];
        s0 = add(u0, v0);
            s1 = add(u1, v1);
        t0 = sub(u0, v0);
            t1 = sub(u1, v1);
        u0 = reduce_to_pm1n(s0+s0, N, Ninv);
            u1 = reduce_to_pm1n(s1+s1, N, Ninv);
        v0 = mulmod2(t0, c1, N, Ninv);
            v1 = mulmod2(t1, c1, N, Ninv);
        s0 = mulmod2(s0, c2, N, Ninv);
            s1 = mulmod2(s1, c2, N, Ninv);
        t0 = mulmod2(t0, c3, N, Ninv);
            t1 = mulmod2(t1, c3, N, Ninv);
        s0 = add(s0, t0);
            s1 = add(s1, t1);
        X0[i+0] = u0;
            X0[i+1] = u1;
        X1[i+0] = v0;
            X1[i+1] = v1;
        X2[i+0] = s0;
            X2[i+1] = s1;
    } while (i += 2, i < len);
}

/*
k = 2, n = 2, z = 2, f = false
[   2       2]
[2//w   -2//w]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,2,2,false>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD W = UNLIKELY(j == 0) ? -w2s[0] : w2s[2*jm+1];
    PD c0 = 2.0;
    PD c1 = reduce_pm1n_to_pmhn(-2.0*W, N);

    ulong i = 0; do {
        PD u0, v0, s0, t0, u1, v1, s1, t1;
        u0 = X0[i+0];
            u1 = X0[i+1];
        v0 = X1[i+0];
            v1 = X1[i+1];
        s0 = add(u0, v0);
            s1 = add(u1, v1);
        t0 = sub(u0, v0);
            t1 = sub(u1, v1);
        u0 = mulmod2(s0, c0, N, Ninv);
            u1 = mulmod2(s1, c0, N, Ninv);
        v0 = mulmod2(t0, c1, N, Ninv);
            v1 = mulmod2(t1, c1, N, Ninv);
        X0[i+0] = u0;
            X0[i+1] = u1;
        X1[i+0] = v0;
            X1[i+1] = v1;
    } while (i += 2, i < len);
}

/*
k = 2, n = 1, z = 4, f = true
[4        -w   -w^2        -w^3]
[1   -1//2*w      0   -1//2*w^3]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,1,4,true>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD f2 = 2.0;
    PD w2 = w2s[j];
    PD wo2 = reduce_pm1n_to_pmhn(mulmod2(w2s[2*j], fnmadd(0.5, N, 0.5), N, Ninv), N);

    ulong i = 0; do {
        PD a0, b0, c0, d0, u0;
        a0 = X0[i];
        a0 = reduce_to_pm1n(a0, N, Ninv);
        b0 = X1[i];
        c0 = X2[i];
        d0 = X3[i];
        c0 = mulmod2(c0, w2, N, Ninv);
        d0 = mulmod2(d0, w2, N, Ninv);
        b0 = add(b0, d0);
        b0 = mulmod2(b0, wo2, N, Ninv);
        u0 = fmsub(f2, a0, b0);
        b0 = sub(a0, b0);
        a0 = reduce_to_pm1n(add(u0, u0), N, Ninv);
        a0 = sub(a0, c0);
        X0[i] = a0;
        X1[i] = b0;
    } while (i += 1, i < len);
}

/*
k = 2, n = 1, z = 4, f = false
[4   -w   -w^2   -w^3]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,1,4,false>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD f1 = 4.0;
    PD w2 = w2s[j];
    PD w  = w2s[2*j];

    ulong i = 0; do {
        PD a0, b0, c0, d0, u0, a1, b1, c1, d1, u1;
        a0 = X0[i+0];
        b0 = X1[i+0];
        c0 = X2[i+0];
        d0 = X3[i+0];
        d0 = mulmod2(d0, w, N, Ninv);
        a0 = mulmod2(a0, f1, N, Ninv);
        b0 = mulmod2(b0, w, N, Ninv);
        a0 = sub(a0, b0);
        c0 = mulmod2(add(c0, d0), w2, N, Ninv);
        X0[i+0] = sub(a0, c0);
    } while (i += 1, i < len);
}

/*
k = 2, n = 1, z = 1, f = true
[4]
[1]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,1,1,true>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    PD f1 = 4.0;

    ulong i = 0; do {
        packed<double,4> a0, a1;
        a0 = X0[i+0];
            a1 = X0[i+1];
        X0[i+0] = reduce_to_pm1n(f1*a0, N, Ninv);
            X0[i+1] = reduce_to_pm1n(f1*a1, N, Ninv);
        X1[i+0] = reduce_to_pm1n(a0, N, Ninv);
            X1[i+1] = reduce_to_pm1n(a1, N, Ninv);
    } while (i += 2, i < len);
}

/*
k = 2, n = 1, z = 1, f = false
[4]
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,1,1,false>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    PD f = 4.0;

    ulong i = 0; do {
        PD a0, a1;
        a0 = X0[i+0];
            a1 = X0[i+1];
        X0[i+0] = reduce_to_pm1n(f*a0, N, Ninv);
            X0[i+1] = reduce_to_pm1n(f*a1, N, Ninv);
    } while (i += 2, i < len);
}

/*
    k = 2, n = 0, z = 4, f = true
    [1//4   1//4*w   1//4*w^2   1//4*w^3]

    {x0} = {1/4*x0 + w/4*x1 + w^2/4*(x2 + w*x3)}
*/
template<> void radix_4_moth_inv_trunc_loop<vsize,0,4,true>(
    pd_fft_ctx<vsize>* Q,
    PD* x, ulong s, ulong len,
    ulong j, ulong jm)
{
    assert(len > 0);
    assert(len%2 == 0);
    PD* X0 = x + 0*s, * X1 = x + 1*s, * X2 = x + 2*s, * X3 = x + 3*s;
    PD N = Q->n, Ninv = Q->ninv;
    const PD* w2s = Q->wrevtab.data();
    PD ff0 = fnmadd(0.25, N, 0.25);
    PD ww = w2s[2*j];
    PD ww2 = w2s[j];
    PD f0 = ff0;
    PD  wo4 = reduce_pm1n_to_pmhn(mulmod2(ww, ff0, N, Ninv), N);
    PD w2o4 = reduce_pm1n_to_pmhn(mulmod2(ww2, ff0, N, Ninv), N);
    PD w = ww;

    ulong i = 0; do {
        PD a0, b0, c0, d0, a1, b1, c1, d1;
        d0 = X3[i];
        a0 = X0[i];
        b0 = X1[i];
        c0 = X2[i];
        d0 = mulmod2(d0, w, N, Ninv);
        a0 = mulmod2(a0, f0, N, Ninv);
        c0 = add(c0, d0);
        b0 = mulmod2(b0, wo4, N, Ninv);
        a0 = add(a0, b0);
        c0 = mulmod2(c0, w2o4, N, Ninv);
        a0 = add(a0, c0);
        X0[i] = a0;
    } while (i += 1, i < len);
}

#undef vsize

#undef PD
