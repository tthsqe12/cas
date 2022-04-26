template <int vsize>
void pd_fft_ctx<vsize>::fft_full_iter(packed<double,vsize>* x, ulong k, ulong j)
{
    assert(k >= 3);
    ulong blksz = pow2(k);
    ulong m;
    const packed<double,vsize>* w2s = wrevtab.data();

    for (m = 0; m + 3 < k; m += 2)
    {
        blksz = pow2(k - 2 - m); // blksz >= 2
        ulong jj = j<<m;
        packed<double,vsize>* xx = x;
        packed<double,vsize> p = n;
        packed<double,vsize> pinv = ninv;
        for (ulong nblks = pow2(m); nblks > 0; nblks--)
        {
            packed<double,vsize> w2 = w2s[jj];
            packed<double,vsize>  w = w2s[2*jj+0];
            packed<double,vsize> iw = w2s[2*jj+1];
            for (ulong i = 0; i < blksz; i += 2)
            {
                packed<double,vsize> x0, x1, x2, x3, y0, y1, y2, y3;
                packed<double,vsize> u0, u1, u2, u3, v0, v1, v2, v3;
                x0 = xx[i+0+0*blksz];
                    u0 = xx[i+1+0*blksz];
                x1 = xx[i+0+1*blksz];
                    u1 = xx[i+1+1*blksz];
                x2 = xx[i+0+2*blksz];
                    u2 = xx[i+1+2*blksz];
                x3 = xx[i+0+3*blksz];
                    u3 = xx[i+1+3*blksz];
                x0 = reduce_to_pm1n(x0, p, pinv);
                    u0 = reduce_to_pm1n(u0, p, pinv);
                x2 = mulmod2(x2, w2, p, pinv);
                    u2 = mulmod2(u2, w2, p, pinv);
                x3 = mulmod2(x3, w2, p, pinv);
                    u3 = mulmod2(u3, w2, p, pinv);
                y0 = add(x0, x2);
                    v0 = add(u0, u2);
                y1 = add(x1, x3);
                    v1 = add(u1, u3);
                y2 = sub(x0, x2);
                    v2 = sub(u0, u2);
                y3 = sub(x1, x3);
                    v3 = sub(u1, u3);
                y1 = mulmod2(y1, w, p, pinv);
                    v1 = mulmod2(v1, w, p, pinv);
                y3 = mulmod2(y3, iw, p, pinv);
                    v3 = mulmod2(v3, iw, p, pinv);
                x0 = add(y0, y1);
                    u0 = add(v0, v1);
                x1 = sub(y0, y1);
                    u1 = sub(v0, v1);
                x2 = add(y2, y3);
                    u2 = add(v2, v3);
                x3 = sub(y2, y3);
                    u3 = sub(v2, v3);
                xx[i+0+0*blksz] = x0;
                    xx[i+1+0*blksz] = u0;
                xx[i+0+1*blksz] = x1;
                    xx[i+1+1*blksz] = u1;
                xx[i+0+2*blksz] = x2;
                    xx[i+1+2*blksz] = u2;
                xx[i+0+3*blksz] = x3;
                    xx[i+1+3*blksz] = u3;
            }
            jj++;
            xx += 4*blksz;
        }
        assert(xx == x + pow2(k));
    }
#if 1
    if (m + 3 == k)
    {
        ulong jj = j<<m;
        packed<double,vsize>* xx = x;
        packed<double,vsize> p = n;
        packed<double,vsize> pinv = ninv;
        for (ulong nblks = pow2(m); nblks > 0; nblks--, jj += 1, xx += 8)
        {
            packed<double,vsize> x0, x1, x2, x3, x4, x5, x6, x7;
            packed<double,vsize> y0, y1, y2, y3, y4, y5, y6, y7;
            packed<double,vsize> z0, z1, z2, z3, z4, z5, z6, z7;
            x0 = reduce_to_pm1n(xx[0], p, pinv);
            x1 = reduce_to_pm1n(xx[1], p, pinv);
            x2 = reduce_to_pm1n(xx[2], p, pinv);
            x3 = reduce_to_pm1n(xx[3], p, pinv);
            x4 = mulmod2(xx[4], w2s[jj], p, pinv);
            x5 = mulmod2(xx[5], w2s[jj], p, pinv);
            x6 = mulmod2(xx[6], w2s[jj], p, pinv);
            x7 = mulmod2(xx[7], w2s[jj], p, pinv);
            y0 = add(x0, x4);
            y1 = add(x1, x5);
            y2 = add(x2, x6);
            y3 = add(x3, x7);
            y4 = sub(x0, x4);
            y5 = sub(x1, x5);
            y6 = sub(x2, x6);
            y7 = sub(x3, x7);
            y2 = mulmod2(y2, w2s[2*jj], p, pinv);
            y3 = mulmod2(y3, w2s[2*jj], p, pinv);
            y6 = mulmod2(y6, w2s[2*jj+1], p, pinv);
            y7 = mulmod2(y7, w2s[2*jj+1], p, pinv);
            z0 = add(y0, y2);
            z1 = add(y1, y3);
            z2 = sub(y0, y2);
            z3 = sub(y1, y3);
            z4 = add(y4, y6);
            z5 = add(y5, y7);
            z6 = sub(y4, y6);
            z7 = sub(y5, y7);
            z1 = mulmod2(z1, w2s[4*jj], p, pinv);
            z3 = mulmod2(z3, w2s[4*jj+1], p, pinv);
            z5 = mulmod2(z5, w2s[4*jj+2], p, pinv);
            z7 = mulmod2(z7, w2s[4*jj+3], p, pinv);
            xx[0] = add(z0, z1);
            xx[1] = sub(z0, z1);
            xx[2] = add(z2, z3);
            xx[3] = sub(z2, z3);
            xx[4] = add(z4, z5);
            xx[5] = sub(z4, z5);
            xx[6] = add(z6, z7);
            xx[7] = sub(z6, z7);
        }
    }
    else
#endif
    if (m + 2 == k)
    {
        ulong jj = j<<m;
        packed<double,vsize>* xx = x;
        packed<double,vsize> p = n;
        packed<double,vsize> pinv = ninv;
        for (ulong nblks = pow2(m)/2; nblks > 0; nblks--, jj += 2, xx += 4*2)
        {
            packed<double,vsize> x0, x1, x2, x3, y0, y1, y2, y3;
            packed<double,vsize> u0, u1, u2, u3, v0, v1, v2, v3;
            x0 = xx[0];
            x1 = xx[1];
            x2 = xx[2];
            x3 = xx[3];
                u0 = xx[4+0];
                u1 = xx[4+1];
                u2 = xx[4+2];
                u3 = xx[4+3];
            x0 = reduce_to_pm1n(x0, p, pinv);
                u0 = reduce_to_pm1n(u0, p, pinv);
            x2 = mulmod2(x2, w2s[(jj+0)], p, pinv);
                u2 = mulmod2(u2, w2s[(jj+1)], p, pinv);
            x3 = mulmod2(x3, w2s[(jj+0)], p, pinv);
                u3 = mulmod2(u3, w2s[(jj+1)], p, pinv);
            y0 = add(x0, x2);
                v0 = add(u0, u2);
            y1 = add(x1, x3);
                v1 = add(u1, u3);
            y2 = sub(x0, x2);
                v2 = sub(u0, u2);
            y3 = sub(x1, x3);
                v3 = sub(u1, u3);
            y1 = mulmod2(y1, w2s[2*(jj+0)+0], p, pinv);
                v1 = mulmod2(v1, w2s[2*(jj+1)+0], p, pinv);
            y3 = mulmod2(y3, w2s[2*(jj+0)+1], p, pinv);
                v3 = mulmod2(v3, w2s[2*(jj+1)+1], p, pinv);
            x0 = add(y0, y1);
                u0 = add(v0, v1);
            x1 = sub(y0, y1);
                u1 = sub(v0, v1);
            x2 = add(y2, y3);
                u2 = add(v2, v3);
            x3 = sub(y2, y3);
                u3 = sub(v2, v3);
            xx[0] = x0;
            xx[1] = x1;
            xx[2] = x2;
            xx[3] = x3;
                xx[4+0] = u0;
                xx[4+1] = u1;
                xx[4+2] = u2;
                xx[4+3] = u3;
        }
        assert(xx == x + pow2(k));
    }
    else if (m + 1 == k)
    {
        assert(m >= 2);
        ulong jj = j<<m;
        packed<double,vsize>* xx = x;
        packed<double,vsize> p = n;
        packed<double,vsize> pinv = ninv;
        for (ulong nblks = pow2(m); nblks > 0; nblks--, jj += 1, xx += 2)
        {
            packed<double,vsize> x0, x1;
            x0 = xx[0];
            x1 = xx[1];
            x0 = reduce_to_pm1n(x0, p, pinv);
            x1 = mulmod2(x1, w2s[jj], p, pinv);
            xx[0] = add(x0, x1);
            xx[1] = sub(x0, x1);
        }
        assert(xx == x + pow2(k));
    }
}

template <int vsize>
void pd_fft_ctx<vsize>::fft_trunc(
    packed<double,vsize>* x,
    ulong k,
    ulong j,
    ulong itrunc, ulong otrunc)
{
    assert(itrunc%2 == 0);
    assert(itrunc <= pow2(k));
    assert(otrunc <= pow2(k));

    if (otrunc < 1)
        return;

    if (itrunc < 1)
    {
        for (ulong i = 0; i < otrunc; i++)
            x[i] = 0.0;
        return;
    }
    else if (itrunc == 1)
    {
        packed<double,vsize> x0 = x[0];
        for (ulong i = 1; i < otrunc; i++)
            x[i] = x0;
        return;
    }

    if (k < 1)
        return;

    if (k == 1)
    {
        packed<double,vsize> x0(0.0), x1(0.0);
        if (0 < itrunc) x0 = x[0];
        if (1 < itrunc) x1 = x[1];
        x0 = reduce_to_pm1n(x0, n, ninv);
        x1 = mulmod2(x1, w2rev(j), n, ninv);
        if (0 < otrunc) x[0] = add(x0, x1);
        if (1 < otrunc) x[1] = sub(x0, x1);
        return;
    }
    else if (k == 2)
    {
        packed<double,vsize> x0(0.0), x1(0.0), x2(0.0), x3(0.0), y0, y1, y2, y3;
        if (0 < itrunc) x0 = x[0];
        if (1 < itrunc) x1 = x[1];
        if (2 < itrunc) x2 = x[2];
        if (3 < itrunc) x3 = x[3];
        x0 = reduce_to_pm1n(x0, n, ninv);
        x2 = mulmod2(x2, w2rev(j), n, ninv);
        x3 = mulmod2(x3, w2rev(j), n, ninv);
        y0 = add(x0, x2);
        y1 = add(x1, x3);
        y2 = sub(x0, x2);
        y3 = sub(x1, x3);
        y1 = mulmod2(y1, w2rev(2*j+0), n, ninv);
        y3 = mulmod2(y3, w2rev(2*j+1), n, ninv);
        x0 = add(y0, y1);
        x1 = sub(y0, y1);
        x2 = add(y2, y3);
        x3 = sub(y2, y3);
        if (0 < otrunc) x[0] = x0;
        if (1 < otrunc) x[1] = x1;
        if (2 < otrunc) x[2] = x2;
        if (3 < otrunc) x[3] = x3;
        return;
    }

    ulong l = pow2(k);
    ulong k1 = 2;
    ulong k2 = k - 2;

    if (k <= 10 && itrunc == otrunc && otrunc == l)
    {
        fft_full_iter(x, k, j);
        return;
    }

    ulong l2 = l/4;
    ulong n1 = otrunc >> (k - 2);   /* top two bits of otrunc */
    ulong n2 = otrunc & (l/4 - 1);  /* all but the top two bits */
    ulong z1 = itrunc >> (k - 2);   /* top two bits of itrunc */
    ulong z2 = itrunc & (l/4 - 1);  /* all but the top two bits */
    ulong n1p = n1 + (n2 != 0);
    ulong z2p = std::min(l2, itrunc);

#define IT(ii, oo) fft_moth_trunc_loop<vsize,ii,oo>
#define LOOKUP_IT(ii, oo) tab[(oind) + 4*(ii)]

    static void (*tab[20])(pd_fft_ctx<vsize>&, packed<double,vsize>*,
                           ulong, ulong, ulong) =
                    {IT(0,1), IT(0,2), IT(0,3), IT(0,4),
                     IT(1,1), IT(1,2), IT(1,3), IT(1,4),
                     IT(2,1), IT(2,2), IT(2,3), IT(2,4),
                     IT(3,1), IT(3,2), IT(3,3), IT(3,4),
                     IT(4,1), IT(4,2), IT(4,3), IT(4,4)};

    ulong oind = ((otrunc-1)>>(k-2));
    ulong iind = ((itrunc-1)>>(k-2));
    ulong it = 1 + ((itrunc-1)&(l2-1));

    /* left half */
    LOOKUP_IT(iind+1,oind)(*this, x, l2, j, it);
    /* right half */
    if (l2 > it)
        LOOKUP_IT(iind+0,oind)(*this, x + it, l2, j, l2 - it);

    // full rows
    for (ulong b = 0; b < n1; b++)
        fft_trunc(x + (b << k2), k2, (j << k1) + b, z2p, l2);

    // last partial row
    if (n2 > 0)
        fft_trunc(x + (n1<< k2), k2, (j << k1) + n1, z2p, n2);

#undef LOOKUP_IT
#undef IT

    return;
}





template <int vsize> template<bool j_is_0>
void pd_fft_ctx<vsize>::ifft_full_iter(packed<double,vsize>* x, ulong k, ulong j)
{
    assert(j_is_0 == (j == 0));
    assert(k >= 3);
    const packed<double,vsize>* w2s = wrevtab.data();
    packed<double,vsize> N = n;
    packed<double,vsize> Ninv = ninv;
    ulong m = k;
    ulong jm = j^(saturate_bits(j)>>1);

    if (m%2 == 1)
    {
        ulong jj = j<<(m-1);
        ulong jjm = (jm<<(m-1)) + (pow2(m-1)-1);
        packed<double,vsize> w = j_is_0 ? neg(w2s[0]) : w2s[jjm];
        packed<double,vsize>* xx = x;
        ulong nblks = pow2(m-1);
        do {
            packed<double,vsize> u0, u1, x0, x1;
            x0 = xx[0];
            x1 = xx[1];
            u0 = add(x0, x1);
            u1 = sub(x1, x0);
            xx[0] = reduce_to_pm1n(u0, N, Ninv);
            xx[1] = mulmod2(u1, w, N, Ninv);
            jj++;
            jjm = j_is_0 ? jj^(saturate_bits(jj)>>1) : jjm-1;
            w = w2s[jjm];
            xx += 2;
        } while (nblks -= 1, nblks > 0);
        assert(xx == x + pow2(k));

        m -= 1;
    }
    else
    {
        ulong jj = j<<(m-2);
        ulong jjm = (jm<<(m-2)) + (pow2(m-2)-1);
        packed<double,4> W  = j_is_0 ? -w2s[0] : w2s[2*(jjm)+1];
        packed<double,4> W2 = j_is_0 ? -w2s[0] : w2s[(jjm)];
        packed<double,4> IW = j_is_0 ?  w2s[1] : w2s[2*(jjm)];
        packed<double,vsize>* xx = x;
        ulong nblks = pow2(m-2);
        do {
            packed<double,4> x0, x1, x2, x3, y0, y1, y2, y3;
            x0 = xx[0];
            x1 = xx[1];
            x2 = xx[2];
            x3 = xx[3];
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
            xx[0] = x0;
            xx[1] = x1;
            xx[2] = x2;
            xx[3] = x3;
            jj++;
            jjm = j_is_0 ? jj^(saturate_bits(jj)>>1) : jjm-1;
            W  =  w2s[2*(jjm)+1];
            W2 =  w2s[(jjm)];
            IW =  w2s[2*(jjm)];
            xx += 4;
        } while (nblks -= 1, nblks > 0);
        assert(xx == x + pow2(k));

        m -= 2;
    }

#if 0
    for (; m > 0; m -= 1)
    {
        ulong blksz = pow2(k-m);
        ulong jj = j<<(m-1);
        packed<double,vsize>* xx = x;
        for (ulong nblks = pow2(m-1); nblks > 0; nblks--)
        {
            ulong jjm = jj^(saturate_bits(jj)>>1);
            packed<double,vsize>  w = (jj == 0) ? w2s[0] : neg(w2s[jjm]);
            for (ulong i = 0; i < blksz; i += 1)
            {
                packed<double,vsize> u0, u1, x0, x1;
                x0 = xx[i+0*blksz];
                x1 = xx[i+1*blksz];
                u0 = add(x0, x1);
                u1 = sub(x0, x1);
                xx[i+0*blksz] = reduce_to_pm1n(u0, N, Ninv);
                xx[i+1*blksz] = mulmod2(u1, w, N, Ninv);
            }
            jj++;
            xx += 2*blksz;
        }
        assert(xx == x + pow2(k));
    }
#else
    assert(m%2 == 0);
    for (; m >= 2; m -= 2)
    {
        ulong blksz = pow2(k-m);
        ulong jj = j<<(m-2);
        ulong jjm = (jm<<(m-2)) + (pow2(m-2)-1);
        packed<double,4> W  = j_is_0 ? neg(w2s[0]) : w2s[2*(jjm)+1];
        packed<double,4> W2 = j_is_0 ? neg(w2s[0]) : w2s[(jjm)];
        packed<double,4> IW = j_is_0 ?     w2s[1]  : w2s[2*(jjm)];
        packed<double,vsize>* xx = x;
        ulong nblks = pow2(m-2);
        do {
            ulong i = 0;
            do {
                packed<double,vsize> x0, x1, x2, x3, y0, y1, y2, y3;
                packed<double,vsize> u0, u1, u2, u3, v0, v1, v2, v3;
                x0 = xx[i+0+0*blksz];
                    u0 = xx[i+1+0*blksz];
                x1 = xx[i+0+1*blksz];
                    u1 = xx[i+1+1*blksz];
                x2 = xx[i+0+2*blksz];
                    u2 = xx[i+1+2*blksz];
                x3 = xx[i+0+3*blksz];
                    u3 = xx[i+1+3*blksz];
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
                xx[i+0+0*blksz] = x0;
                    xx[i+1+0*blksz] = u0;
                xx[i+0+1*blksz] = x1;
                    xx[i+1+1*blksz] = u1;
                xx[i+0+2*blksz] = x2;
                    xx[i+1+2*blksz] = u2;
                xx[i+0+3*blksz] = x3;
                    xx[i+1+3*blksz] = u3;
            } while (i+=2, i < blksz);
            jj++;
            jjm = j_is_0 ? jj^(saturate_bits(jj)>>1) : jjm-1;
            W  = w2s[2*(jjm)+1];
            W2 = w2s[(jjm)];
            IW = w2s[2*(jjm)];
            xx += 4*blksz;
        } while (nblks--, nblks > 0);
        assert(xx == x + pow2(k));
    }
#endif
}

template <int vsize>
void pd_fft_ctx<vsize>::ifft_full(packed<double,vsize>* x, ulong k, ulong j)
{
    if (k == 0)
    {
        return;
    }
    else if (k == 1)
    {
        packed<double,vsize> u = add(x[0], x[1]);
        packed<double,vsize> v = sub(x[0], x[1]);
        x[0] = reduce_to_pm1n(u, n, ninv);
        x[1] = mulmod2(v, w2revinv(j), n, ninv);
        return;
    }
    else if (k == 2)
    {
        const packed<double,4>* w2s = wrevtab.data();
        ulong jm = j^(saturate_bits(j)>>1);
        packed<double,4> W  = (UNLIKELY((j) == 0)) ? neg(w2s[0]) : w2s[2*(jm)+1];
        packed<double,4> W2 = (UNLIKELY((j) == 0)) ? neg(w2s[0]) : w2s[(jm)];
        packed<double,4> IW = (UNLIKELY((j) == 0)) ?     w2s[1]  : w2s[2*(jm)];
        packed<double,4> N = n;
        packed<double,4> Ninv = ninv;
        packed<double,4> x0, x1, x2, x3, y0, y1, y2, y3;
        x0 = x[0];
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
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
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
        x[3] = x3;
        return;
    }

    if (k <= 10)
    {
        if (j == 0)
            ifft_full_iter<true>(x, k, j);
        else
            ifft_full_iter<false>(x, k, j);
        return;
    }

    ulong l2 = pow2(k-2);
    ifft_full(x+0*l2, k-2, 4*j+0);
    ifft_full(x+1*l2, k-2, 4*j+1);
    ifft_full(x+2*l2, k-2, 4*j+2);
    ifft_full(x+3*l2, k-2, 4*j+3);

    const packed<double,4>* w2s = wrevtab.data();
    ulong jm = j^(saturate_bits(j)>>1);
    packed<double,4> W  = (UNLIKELY((j) == 0)) ? neg(w2s[0]) : w2s[2*(jm)+1];
    packed<double,4> W2 = (UNLIKELY((j) == 0)) ? neg(w2s[0]) : w2s[(jm)];
    packed<double,4> IW = (UNLIKELY((j) == 0)) ?     w2s[1]  : w2s[2*(jm)];
    packed<double,4> N = n;
    packed<double,4> Ninv = ninv;

    ulong i = 0;
    do {
        packed<double,4> x0, x1, x2, x3, y0, y1, y2, y3;
        x0 = x[i+l2*0];
        x1 = x[i+l2*1];
        x2 = x[i+l2*2];
        x3 = x[i+l2*3];
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
        x[i+l2*0] = x0;
        x[i+l2*1] = x1;
        x[i+l2*2] = x2;
        x[i+l2*3] = x3;
    } while (i += 1, i < l2);

    return;
}

template <int vsize>
void pd_fft_ctx<vsize>::ifft_trunc(
    packed<double,vsize>* x,
    ulong k, // transform length 2^(k)
    ulong j,
    ulong z,   // actual trunc is z
    ulong nn,   // actual trunc is n
    bool f)
{
    assert(nn%4 == 0);
    assert(z%4 == 0);
    assert(nn <= z);
    assert(1 <= z && z <= pow2(k));
    assert(1 <= nn+f && nn+f <= pow2(k));

    if (!f && z == nn && nn == pow2(k))
    {
        ifft_full(x, k, j);
        return;
    }

    assert(k >= 2);

#if 0
    ulong k1 = 1;
    ulong k2 = k - k1;

    ulong l2 = ulong(1) << k2;
    ulong n1 = nn >> k2;
    ulong n2 = nn & (l2 - 1);
    ulong z1 = z >> k2;
    ulong z2 = z & (l2 - 1);
    bool fp = n2 + f > 0;
    ulong z2p = std::min(l2, z);
    ulong m = std::min(n2, z2);
    ulong mp = std::max(n2, z2);

#define IT(zz, nn, ff) radix_2_moth_inv_trunc_loop<vsize,zz,nn,ff>
#define LOOKUP_IT(zz, nn, ff) tab[ulong(ff) + 2*((nn) + 3*((zz)-1))]

    static void (*tab[2*3*2])(pd_fft_ctx<vsize>&, packed<double,4>*, packed<double,4>*, ulong, ulong) =
                    {IT(1,0,false),IT(1,0,true), IT(1,1,false),IT(1,1,true), IT(1,2,false),IT(1,2,true),
                     IT(2,0,false),IT(2,0,true), IT(2,1,false),IT(2,1,true), IT(2,2,false),IT(2,2,true)
                    };

    // complete rows
    for (ulong b = 0; b < n1; b++)
        ifft_full(x + b*l2, k2, (j << k1) + b);

    // rightmost columns
    assert(n2 <= mp);
    assert(mp <= z2p);
    if (mp > n2)
        LOOKUP_IT(z1+1,n1,fp)(*this, x+n2+0*l2, x+n2+1*l2, mp - n2, j);
    if (z2p > mp)
        LOOKUP_IT(z1+0,n1,fp)(*this, x+mp+0*l2, x+mp+1*l2, z2p - mp, j);

    //last partial row
    if (fp)
        ifft_trunc(x + n1*l2, k2, (j << k1) + n1, z2p, n2, f);

    // leftmost columns
    assert(m <= n2);
    if (m > 0)
        LOOKUP_IT(z1+1,n1+1,false)(*this, x+0+0*l2, x+0+1*l2, m, j);
    if (n2 > m)
        LOOKUP_IT(z1+0,n1+1,false)(*this, x+m+0*l2, x+m+1*l2, n2 - m, j);

#undef LOOKUP_IT
#undef IT

#else

    ulong k1 = 2;
    ulong k2 = k - k1;

    ulong l2 = ulong(1) << k2;
    ulong n1 = nn >> k2;
    ulong n2 = nn & (l2 - 1);
    ulong z1 = z >> k2;
    ulong z2 = z & (l2 - 1);
    bool fp = n2 + f > 0;
    ulong z2p = std::min(l2, z);
    ulong m = std::min(n2, z2);
    ulong mp = std::max(n2, z2);

#define IT(nn, zz, ff) radix_4_moth_inv_trunc_loop<vsize,nn,zz,ff>
#define LOOKUP_IT(nn, zz, ff) tab[ulong(ff) + 2*((zz)-1 + 4*(nn))]

    static void (*tab[5*4*2])(pd_fft_ctx<vsize>*, packed<double,4>*, ulong, ulong, ulong, ulong) =
        {IT(0,1,false),IT(0,1,true), IT(0,2,false),IT(0,2,true), IT(0,3,false),IT(0,3,true), IT(0,4,false),IT(0,4,true),
         IT(1,1,false),IT(1,1,true), IT(1,2,false),IT(1,2,true), IT(1,3,false),IT(1,3,true), IT(1,4,false),IT(1,4,true),
         IT(2,1,false),IT(2,1,true), IT(2,2,false),IT(2,2,true), IT(2,3,false),IT(2,3,true), IT(2,4,false),IT(2,4,true),
         IT(3,1,false),IT(3,1,true), IT(3,2,false),IT(3,2,true), IT(3,3,false),IT(3,3,true), IT(3,4,false),IT(3,4,true),
         IT(4,1,false),IT(4,1,true), IT(4,2,false),IT(4,2,true), IT(4,3,false),IT(4,3,true), IT(4,4,false),IT(4,4,true)};

    // complete rows
    for (ulong b = 0; b < n1; b++)
        ifft_full(x + b*l2, k2, (j << k1) + b);

    // rightmost columns
    assert(n2 <= mp);
    assert(mp <= z2p);
    ulong jm = j ^ (saturate_bits(j)>>1);
    if (mp > n2)  LOOKUP_IT(n1,z1+1,fp)(this, x+n2, l2, mp - n2, j, jm);
    if (z2p > mp) LOOKUP_IT(n1,z1+0,fp)(this, x+mp, l2, z2p - mp, j, jm);

    //last partial row
    if (fp)
        ifft_trunc(x + n1*l2, k2, (j << k1) + n1, z2p, n2, f);

    // leftmost columns
    assert(m <= n2);
    if (m > 0)  LOOKUP_IT(n1+1,z1+1,false)(this, x + 0, l2, m, j, jm);
    if (n2 > m) LOOKUP_IT(n1+1,z1+0,false)(this, x + m, l2, n2 - m, j, jm);

#undef LOOKUP_IT
#undef IT

#endif
    return;
}

