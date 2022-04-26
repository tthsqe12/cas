/*
void cpd_fft_ctx::fft_base(complex_packed<double,4>* x, ulong j)
{
    std::complex<double> z0(x[0].real()[0],x[0].imag()[0]),
                         z1(x[0].real()[1],x[0].imag()[1]),
                         z2(x[0].real()[2],x[0].imag()[2]),
                         z3(x[0].real()[3],x[0].imag()[3]);
    radix_4_moth(z0, z1, z2, z3, wrev(4*j), wrev(2*j));
    x[0].real(packed<double,4>(z0.real(), z1.real(), z2.real(), z3.real()));
    x[0].imag(packed<double,4>(z0.imag(), z1.imag(), z2.imag(), z3.imag()));
}

void cpd_fft_ctx::ifft_base(complex_packed<double,4>* x)
{
    std::complex<double> z0(x[0].real()[0],x[0].imag()[0]),
                         z1(x[0].real()[1],x[0].imag()[1]),
                         z2(x[0].real()[2],x[0].imag()[2]),
                         z3(x[0].real()[3],x[0].imag()[3]);
    radix_4_mothi1(z0, z1, z2, z3);
    x[0].real(packed<double,4>(z0.real(), z1.real(), z2.real(), z3.real()));
    x[0].imag(packed<double,4>(z0.imag(), z1.imag(), z2.imag(), z3.imag()));
}
*/

void cpd_fft_ctx::fft_full_block(
    complex_packed<double,4>* x,
    ulong k, // 4 parallel transforms each of length 2^k
    ulong j)
{
    const std::complex<double>* w2r = w2rev_table();
    complex_packed<double,4> w, w2;
    if (k == 0)
    {
    }
    else if (k == 1)
    {
        radix_2_moth_block(x[0], x[1], w2r[j]);
    }
    else if (k == 2)
    {
        radix_4_moth_block(x[0], x[1], x[2], x[3], w2r[2*j], w2r[j]);
    }
    else if (k == 3)
    {
        w = w2r[j];
        complex_packed<double,4> x0 = x[0], x1 = x[1], x2 = x[2], x3 = x[3], x4 = x[4], x5 = x[5], x6 = x[6], x7 = x[7];
        radix_2_moth_block(x0, x4, w);
        radix_2_moth_block(x1, x5, w);
        radix_2_moth_block(x2, x6, w);
        radix_2_moth_block(x3, x7, w);
        radix_4_moth_block(x0, x1, x2, x3, w2r[2*(2*j+0)], w2r[(2*j+0)]);
        x[0] = x0; x[1] = x1; x[2] = x2; x[3] = x3;
        radix_4_moth_block(x4, x5, x6, x7, w2r[2*(2*j+1)], w2r[(2*j+1)]);
        x[4] = x4; x[5] = x5; x[6] = x6; x[7] = x7;
    }
    else
    {
        ulong l2 = pow2(k - 2);
        w = w2r[2*j];
        w2 = w2r[j];
        for (ulong i = 0; i < l2; i++)
            radix_4_moth_block(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
        fft_full_block(x+0*l2, k-2, 4*j+0);
        fft_full_block(x+1*l2, k-2, 4*j+1);
        fft_full_block(x+2*l2, k-2, 4*j+2);
        fft_full_block(x+3*l2, k-2, 4*j+3);
    }
}


void cpd_fft_ctx::fft_full(
    complex_packed<double,4>* x,
    ulong k, // transform length 2^k
    ulong j)
{
    assert(k >= 4);

    // columns: k-2 because 4-wide vectors
    fft_full_block(x, k-2, j);

    // rows of length 4
    ulong i = 0;
    do {
        complex_packed<double,4> t0, t1, t2, t3, w, w2;
        transpose_4x4(t0, t1, t2, t3, x[4*i+0], x[4*i+1], x[4*i+2], x[4*i+3]);
        radix_4_moth_block(t0, t1, t2, t3, w2revpacked(i,1), w2revpacked(i,2));
        transpose_4x4(x[4*i+0], x[4*i+1], x[4*i+2], x[4*i+3], t0, t1, t2, t3);
    } while (++i < pow2(k-4));
}

// input truncation only
void cpd_fft_ctx::fft_trunc_block(
    complex_packed<double,4>* x,
    ulong k,
    ulong j,
    ulong itrunc)
{
    assert(itrunc%4 == 0);
    assert(0 < itrunc);
    assert(itrunc <= pow2(k));

    const std::complex<double>* w2r = w2rev_table();
    complex_packed<double,4> w, w2;

    if (itrunc == pow2(k))
    {
        fft_full_block(x, k, j);
        return;
    }

    assert(k >= 3);

    ulong l2 = pow2(k - 2);
    w = w2r[2*j];
    w2 = w2r[j];

    if (itrunc <= l2)
    {
        for (ulong i = 0; i < itrunc; i++)
            radix_4_moth_trunc_block<1>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
        for (ulong i = itrunc; i < l2; i++)
            radix_4_moth_trunc_block<0>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
    }
    else if (itrunc <= 2*l2)
    {
        itrunc -= l2;
        for (ulong i = 0; i < itrunc; i++)
            radix_4_moth_trunc_block<2>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
        for (ulong i = itrunc; i < l2; i++)
            radix_4_moth_trunc_block<1>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
    }
    else if (itrunc <= 3*l2)
    {
        itrunc -= 2*l2;
        for (ulong i = 0; i < itrunc; i++)
            radix_4_moth_trunc_block<3>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
        for (ulong i = itrunc; i < l2; i++)
            radix_4_moth_trunc_block<2>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
    }
    else
    {
        itrunc -= 3*l2;
        for (ulong i = 0; i < itrunc; i++)
            radix_4_moth_trunc_block<4>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
        for (ulong i = itrunc; i < l2; i++)
            radix_4_moth_trunc_block<3>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
    }

    fft_full_block(x+0*l2, k-2, 4*j+0);
    fft_full_block(x+1*l2, k-2, 4*j+1);
    fft_full_block(x+2*l2, k-2, 4*j+2);
    fft_full_block(x+3*l2, k-2, 4*j+3);
}

void cpd_fft_ctx::fft_trunc(
    complex_packed<double,4>* x,
    ulong k, // transform length 4*2^(k)
    ulong itrunc)
{
    assert(0 < itrunc);
    assert(itrunc <= pow2(k));
    assert(itrunc%16 == 0);

    fft_trunc_block(x, k-2, 0, itrunc/4);

    ulong i = 0;
    do {
        complex_packed<double,4> t0, t1, t2, t3;
        transpose_4x4(t0, t1, t2, t3, x[4*i+0], x[4*i+1], x[4*i+2], x[4*i+3]);
        radix_4_moth_block(t0, t1, t2, t3, w2revpacked(i,1), w2revpacked(i,2));
        transpose_4x4(x[4*i+0], x[4*i+1], x[4*i+2], x[4*i+3], t0, t1, t2, t3);
    } while (++i < pow2(k-4));
}


void cpd_fft_ctx::ifft_full_block(
    complex_packed<double,4>* x,
    ulong k)
{
    if (k == 0)
    {
    }
    else if (k == 1)
    {
        radix_2_moth_1_block(x[0], x[1]);
    }
    else if (k == 2)
    {
        radix_4_moth_inv_1_block(x[0], x[1], x[2], x[3]);
    }
    else if (k == 3)
    {
        complex_packed<double,4> x0 = x[0], x1 = x[1], x2 = x[2], x3 = x[3], x4 = x[4], x5 = x[5], x6 = x[6], x7 = x[7];
        radix_2_moth_1_block(x0, x1);
        radix_2_moth_1_block(x2, x3);
        radix_2_moth_1_block(x4, x5);
        radix_2_moth_1_block(x6, x7);
        radix_4_moth_inv_1_block(x0, x2, x4, x6);
        x[0] = x0; x[2] = x2; x[4] = x4; x[6] = x6;
        radix_4_moth_inv_block(x1, x3, x5, x7, wpow(k,1), wpow(k-1,1));
        x[1] = x1; x[3] = x3; x[5] = x5; x[7] = x7;
    }
    else
    {
        ulong l2 = pow2(k-2);
        ifft_full_block(x+0*l2, k-2);
        ifft_full_block(x+1*l2, k-2);
        ifft_full_block(x+2*l2, k-2);
        ifft_full_block(x+3*l2, k-2);
        const std::complex<double>* wk = wpow_table(k);
        const std::complex<double>* wk1 = wpow_table(k-1);
        for (ulong i = 0; i < l2; i++)
            radix_4_moth_inv_block(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2],
                                   wk[i], wk1[i]);
    }
}


void cpd_fft_ctx::ifft_full(
    complex_packed<double,4>* x,
    ulong k) // transformation length 2^k
{
    assert(k >= 4);

    ulong i = 0;
    do {
        complex_packed<double,4> t0, t1, t2, t3;
        transpose_4x4(t0, t1, t2, t3, x[4*i+0], x[4*i+1], x[4*i+2], x[4*i+3]);
        radix_4_moth_rev_block(t0, t1, t2, t3, w2revpacked(i,1), w2revpacked(i,2), w2revpacked(i,3));
        transpose_4x4(x[4*i+0], x[4*i+1], x[4*i+2], x[4*i+3], t0, t1, t2, t3);
    } while (++i < pow2(k-4));

    ifft_full_block(x, k-2);
}

void cpd_fft_ctx::fft_three_trunc(
    complex_packed<double,4>* x,
    ulong k, // transformation length 3*2^k
    ulong itrunc)
{
    assert(k >= 4);
    assert(itrunc%4 == 0);

    ulong l = pow2(k-2);
    itrunc /= 4;
    assert(itrunc <= 3*l);

    const complex_packed<double,4>* w = wthreepow_table(k-2);

    if (itrunc <= l)
    {
        for (ulong i = 0; i < itrunc; i++)
            radix_3_moth_trunc_block<1>(x[i+0*l], x[i+1*l], x[i+2*l], w[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_3_moth_trunc_block<0>(x[i+0*l], x[i+1*l], x[i+2*l], w[i]);
    }
    else if (itrunc <= 2*l)
    {
        itrunc -= l;
        for (ulong i = 0; i < itrunc; i++)
            radix_3_moth_trunc_block<2>(x[i+0*l], x[i+1*l], x[i+2*l], w[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_3_moth_trunc_block<1>(x[i+0*l], x[i+1*l], x[i+2*l], w[i]);
    }
    else
    {
        itrunc -= 2*l;
        for (ulong i = 0; i < itrunc; i++)
            radix_3_moth_trunc_block<3>(x[i+0*l], x[i+1*l], x[i+2*l], w[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_3_moth_trunc_block<2>(x[i+0*l], x[i+1*l], x[i+2*l], w[i]);
    }

    fft_full(x + 0*l, k, 0);
    fft_full(x + 1*l, k, 0);
    fft_full(x + 2*l, k, 0);
}



void cpd_fft_ctx::ifft_three_full(
    complex_packed<double,4>* x,
    ulong k) // transformation length 3*2^k
{
    assert(k >= 4);
    ulong l = pow2(k-2);

    ifft_full(x + 0*l, k);
    ifft_full(x + 1*l, k);
    ifft_full(x + 2*l, k);

    for (ulong i = 0; i < l; i++)
        radix_3_moth_inv_block(x[i+0*l], x[i+1*l], x[i+2*l], wthreepow(k-2,i));
}


void cpd_fft_ctx::fft_five_trunc(
    complex_packed<double,4>* x,
    ulong k, // transformation length 5*2^k
    ulong itrunc)
{
    assert(k >= 4);
    assert(itrunc%4 == 0);

    ulong l = pow2(k-2);
    itrunc /= 4;
    assert(itrunc <= 5*l);

    const complex_packed<double,4>* w = wfive1pow_table(k-2);
    const complex_packed<double,4>* w2 = wfive2pow_table(k-2);

    if (itrunc <= l)
    {
        for (ulong i = 0; i < itrunc; i++)
            radix_5_moth_trunc_block<1>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_5_moth_trunc_block<0>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
    }
    else if (itrunc <= 2*l)
    {
        itrunc -= l;
        for (ulong i = 0; i < itrunc; i++)
            radix_5_moth_trunc_block<2>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_5_moth_trunc_block<1>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
    }
    else if (itrunc <= 3*l)
    {
        itrunc -= 2*l;
        for (ulong i = 0; i < itrunc; i++)
            radix_5_moth_trunc_block<3>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_5_moth_trunc_block<2>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
    }
    else if (itrunc <= 4*l)
    {
        itrunc -= 3*l;
        for (ulong i = 0; i < itrunc; i++)
            radix_5_moth_trunc_block<4>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_5_moth_trunc_block<3>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
    }
    else
    {
        itrunc -= 4*l;
        for (ulong i = 0; i < itrunc; i++)
            radix_5_moth_trunc_block<5>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
        for (ulong i = itrunc; i < l; i++)
            radix_5_moth_trunc_block<4>(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], w[i], w2[i]);
    }

    fft_full(x + 0*l, k, 0);
    fft_full(x + 1*l, k, 0);
    fft_full(x + 2*l, k, 0);
    fft_full(x + 3*l, k, 0);
    fft_full(x + 4*l, k, 0);
}

void cpd_fft_ctx::ifft_five_full(
    complex_packed<double,4>* x,
    ulong k) // transformation length 5*2^k
{
    assert(k >= 4);
    ulong l = pow2(k-2);

    ifft_full(x + 0*l, k);
    ifft_full(x + 1*l, k);
    ifft_full(x + 2*l, k);
    ifft_full(x + 3*l, k);
    ifft_full(x + 4*l, k);

    for (ulong i = 0; i < l; i++)
        radix_5_moth_inv_block(x[i+0*l], x[i+1*l], x[i+2*l], x[i+3*l], x[i+4*l], wfive1pow(k-2,i), wfive2pow(k-2,i));
}

