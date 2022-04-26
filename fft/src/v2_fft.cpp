void fftv2_ctx::fft_main_block(
   ulong I, // starting index
   ulong S, // stride
   ulong k, // transform length 2^k
   ulong j)
{
    if (k > 8)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        // column ffts
        ulong l2 = pow2(k2);
        ulong a = 0; do {
            fft_main_block(I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        // row ffts
        ulong l1 = pow2(k1);
        ulong b = 0; do {
            fft_main_block(I + (b<<k2)*S, S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        return;
    }
    else if (k > 2)
    {
        ulong k1 = 2;
        ulong k2 = k - k1;
        ulong l2 = pow2(k2);

        // column ffts
        PD_RADIX_4_FORWARD_PARAM(j)
        for (ulong a = 0; a < l2; a++)
        {
            double* X0 = from_index(I+a*S + (S<<k2)*0);
            double* X1 = from_index(I+a*S + (S<<k2)*1);
            double* X2 = from_index(I+a*S + (S<<k2)*2);
            double* X3 = from_index(I+a*S + (S<<k2)*3);
            ulong i = 0; do {
                PD_RADIX_4_FORWARD_2X(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
            } while (i += 8, i < BLK_SZ);
        }
        
        // row ffts
        for (ulong b = 0; b < pow2(k1); b++)
            fft_main_block(I + (b<<k2)*S, S, k2, (j<<k1) + b);

        return;        
    }
    else if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
        PD_RADIX_4_FORWARD_PARAM(j)
        ulong i = 0; do {
            PD_RADIX_4_FORWARD_2X(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        } while (i += 8, i < blk_sz);
    }
    else if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        PD_RADIX_2_FORWARD_PARAM(j)
        ulong i = 0; do {
            PD_RADIX_2_FORWARD_2X(X0+i, X1+i, X0+i+4, X1+i+4);
        } while (i += 8, i < blk_sz);
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
        ulong l1 = pow2(k1);
        ulong b = 0; do {
            ifft_main_block(I + (b<<k2)*S, S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        // column ffts
        ulong l2 = pow2(k2);
        ulong a = 0; do {
            ifft_main_block(I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        return;
    }

    ulong jm = j^(saturate_bits(j)>>1);

    if (k == 2)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        double* X2 = from_index(I + S*2);
        double* X3 = from_index(I + S*3);
        PD_RADIX_4_REVERSE_PARAM(j, jm, true)
        ulong i = 0; do {
            PD_RADIX_4_REVERSE_2X(X0+i, X1+i, X2+i, X3+i, X0+i+4, X1+i+4, X2+i+4, X3+i+4);
        } while(i += 8, i < blk_sz);
    }
    else if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        PD_RADIX_2_REVERSE_PARAM(j)
        ulong i = 0; do {
            PD_RADIX_2_REVERSE_2X(X0+i, X1+i, X0+i+4, X1+i+4);
        } while (i += 8, i < blk_sz);
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

        ulong l2 = pow2(k2);
        ulong a = 0; do {
            fft_main_block(I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        ulong l1 = pow2(k1);
        ulong b = 0; do {
            fft_main(I + b*(S<<k2), S, k2, (j<<k1) + b);
        } while (b++, b < l1);

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

        ulong l1 = pow2(k1);
        ulong b = 0; do {
            ifft_main(I + b*(S<<k2), S, k2, (j<<k1) + b);
        } while (b++, b < l1);

        ulong l2 = pow2(k2);
        ulong a = 0; do {
            ifft_main_block(I + a*S, S<<k2, k1, j);
        } while (a++, a < l2);

        return;
    }

    if (k == 2)
    {
        // k1 = 2; k2 = 0
        ifft_base<true>(I+S*0, 4*j+0);
        ifft_base<false>(I+S*1, 4*j+1);
        ifft_base<false>(I+S*2, 4*j+2);
        ifft_base<false>(I+S*3, 4*j+3);
        ifft_main_block(I, S, 2, j);
    }
    else if (k == 1)
    {
        // k1 = 1; k2 = 0
        ifft_base<true>(I+S*0, 2*j+0);
        ifft_base<false>(I+S*1, 2*j+1);
        ifft_main_block(I, S, 1, j);
    }
    else
    {
        ifft_base<true>(I, j);
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
    assert(itrunc <= pow2(k));
    assert(otrunc <= pow2(k));

    if (otrunc < 1)
        return;

    if (itrunc <= 1)
    {
        if (itrunc < 1)
        {
            for (ulong a = 0; a < otrunc; a++)
            {
                double* X0 = from_index(I + S*a);
                packed<double,VEC_SZ> z;
                z.zero();
                ulong i = 0; do {
                    z.store(X0 + i);
                } while (i += VEC_SZ, i < blk_sz);
            }
        }
        else
        {
            double* X0 = from_index(I + S*0);
            for (ulong a = 1; a < otrunc; a++)
            {
                double* X1 = from_index(I + S*a);
                ulong i = 0; do {
                    packed<double,VEC_SZ> x0, x1;
                    x0.load(X0 + i + 0*VEC_SZ);
                    x1.load(X0 + i + 1*VEC_SZ);
                    x0.store(X1 + i + 0*VEC_SZ);
                    x1.store(X1 + i + 1*VEC_SZ);
                } while (i += 2*VEC_SZ, i < blk_sz);
            }
        }

        return;
    }

    if (itrunc == otrunc && otrunc == pow2(k))
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
#define IT(ii, oo) fftv2_moth_trunc_block<ii,oo>
#define LOOKUP_IT(ii, oo) tab[(oo)-1 + 4*((ii)-2)]

        static void (*tab[3*4])(fftv2_ctx*, double*, double*, double*, double*, ulong) =
                        {IT(2,1), IT(2,2), IT(2,3), IT(2,4),
                         IT(3,1), IT(3,2), IT(3,3), IT(3,4),
                         IT(4,1), IT(4,2), IT(4,3), IT(4,4)};

        LOOKUP_IT(itrunc, otrunc)(this, from_index(I+S*0), from_index(I+S*1),
                                        from_index(I+S*2), from_index(I+S*3), j);

#undef LOOKUP_IT
#undef IT
    }
    else if (k == 1)
    {
        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);
        PD_RADIX_2_FORWARD_PARAM(j)

        assert(itrunc == 2 && otrunc == 1);
        ulong i = 0; do {
            PD_RADIX_2_FORWARD_2X_ITRUNC2_OTRUNC1(X0+i, X1+i, X0+i+4, X1+i+4);
        } while (i += 8, i < blk_sz);
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
            packed<double,4> z;
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
#define IT(nn, zz, ff) radix_4_moth_inv_trunc_block<nn,zz,ff>
#define LOOKUP_IT(nn, zz, ff) tab[ulong(ff) + 2*((zz)-1 + 4*(nn))]

        static bool (*tab[5*4*2])(fftv2_ctx*, double*, double*, double*, double*, ulong, ulong) =
            {IT(0,1,false),IT(0,1,true), IT(0,2,false),IT(0,2,true), IT(0,3,false),IT(0,3,true), IT(0,4,false),IT(0,4,true),
             IT(1,1,false),IT(1,1,true), IT(1,2,false),IT(1,2,true), IT(1,3,false),IT(1,3,true), IT(1,4,false),IT(1,4,true),
             IT(2,1,false),IT(2,1,true), IT(2,2,false),IT(2,2,true), IT(2,3,false),IT(2,3,true), IT(2,4,false),IT(2,4,true),
             IT(3,1,false),IT(3,1,true), IT(3,2,false),IT(3,2,true), IT(3,3,false),IT(3,3,true), IT(3,4,false),IT(3,4,true),
             IT(4,1,false),IT(4,1,true), IT(4,2,false),IT(4,2,true), IT(4,3,false),IT(4,3,true), IT(4,4,false),IT(4,4,true)};

        ulong jm = j ^ (saturate_bits(j)>>1);
        if (LOOKUP_IT(n,z,f)(this, from_index(I+S*0), from_index(I+S*1),
                                   from_index(I+S*2), from_index(I+S*3), j, jm))
        {
            return;
        }

#undef LOOKUP_IT
#undef IT
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
#define IT(nn, zz, ff) radix_2_moth_inv_trunc_block<nn,zz,ff>
#define LOOKUP_IT(nn, zz, ff) tab[ulong(ff) + 2*((zz)-1 + 2*(nn))]

        static void (*tab[3*2*2])(fftv2_ctx*, double*, double*, ulong) =
            {IT(0,1,false),IT(0,1,true), IT(0,2,false),IT(0,2,true),
             IT(1,1,false),IT(1,1,true), IT(1,2,false),IT(1,2,true),
             IT(2,1,false),IT(2,1,true), IT(2,2,false),IT(2,2,true)};

        LOOKUP_IT(n,z,f)(this, from_index(I+S*0), from_index(I+S*1), j);

#undef LOOKUP_IT
#undef IT

#if 0

        double* X0 = from_index(I + S*0);
        double* X1 = from_index(I + S*1);

        packed<double,VEC_SZ> N = p;
        packed<double,VEC_SZ> Ninv = pinv;

        if (n == 1 && z == 2 && f)
        {
            // {x0, x1} = {2*x0 - w*x1, x0 - w*x1}
            packed<double,VEC_SZ> w = w2s[j];
            packed<double,VEC_SZ> c = 2.0;
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double,VEC_SZ> u0, u1, v0, v1;
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
            packed<double,VEC_SZ> w = w2s[j];
            packed<double,VEC_SZ> c = 2.0;
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double,VEC_SZ> u0, u1, v0, v1;
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
                packed<double,VEC_SZ> u0, v0, v1;
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
                packed<double,VEC_SZ> u0, v0;
                u0.load(X0 + i);
                u0 = reduce_to_pm1n(u0, N, Ninv);
                v0 = add(u0, u0);
                v0.store(X0 + i);
            }
        }
        else if (n == 0 && z == 2 && f)
        {
            // {x0} = {(x0 + w*x1)/2}
            packed<double,VEC_SZ> w = w2s[j];
            packed<double,VEC_SZ> c = fnmadd(0.5, p, 0.5);
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double,VEC_SZ> u0, u1;
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
            packed<double,VEC_SZ> w = w2s[j];
            packed<double,VEC_SZ> c = 0.5 - 0.5*p;
            for (ulong i = 0; i < blk_sz; i += VEC_SZ)
            {
                packed<double,VEC_SZ> u0, u1;
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
#endif

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

    if (k > 2)
    {
        ulong k1 = k/2;
        ulong k2 = k - k1;

        ulong l2 = pow2(k2);
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
                   ifft_base<true>(I+S*0, 4*j+0);
        if (n > 1) ifft_base<false>(I+S*1, 4*j+1);
        if (n > 2) ifft_base<false>(I+S*2, 4*j+2);
        if (n > 3) ifft_base<false>(I+S*3, 4*j+3);
        ifft_trunc_block(I, S, 2, j, z, n, f);
        if (f) ifft_trunc(I + S*n, S, 0, 4*j+n, 1, 0, f);
        
    }
    else if (k == 1)
    {
                   ifft_base<true>(I+S*0, 2*j+0);
        if (n > 1) ifft_base<false>(I+S*1, 2*j+1);
        ifft_trunc_block(I, S, 1, j, z, n, f);
        if (f) ifft_trunc(I + S*n, S, 0, 2*j+n, 1, 0, f);
    }
    else
    {
        assert(!f);
        ifft_base<true>(I, j);
    }
}



// pointwise mul of self with b and m
void fftv2_ctx::point_mul(const double* b, ulong mm)
{
    double m = mm;
    if (m > 0.5*p)
        m -= p;

    packed<double,4> M = m;
    packed<double,4> n = p;
    packed<double,4> ninv = pinv;
    for (ulong I = 0; I < pow2(depth - LG_BLK_SZ); I++)
    {
        double* x = from_index(I);
        const double* bx = b + offset(I);
        ulong j = 0;
        do {
            packed<double,4> x0, x1, x2, x3, b0, b1, b2, b3;
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


