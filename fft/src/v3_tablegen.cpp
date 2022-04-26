void cpd_fft_ctx::fit_wtab(ulong depth)
{
    depth = std::max(depth, ulong(4));
    ulong l = pow2(depth);
    if (wrevtab.size() < l)
    {
        wrevtab.resize(l);
        for (ulong j = 0; j < l; j++)
            wrevtab[j] = rev_cispi(j, 1);
    }

    if (w2revtab.size() < l/2)
    {
        w2revtab.resize(l/2);
        for (ulong j = 0; j < l/2; j++)
            w2revtab[j] = wrevtab[2*j];
    }

    if (w2revpackedtab.size() < l/4*3)
    {
        w2revpackedtab.resize(l/4*3);
        for (ulong j = 0; j < l/4; j++)
        for (ulong h = 1; h <= 3; h++)
        {
            std::complex<double> w0 = rev_cispi(4*(4*j+0), h);
            std::complex<double> w1 = rev_cispi(4*(4*j+1), h);
            std::complex<double> w2 = rev_cispi(4*(4*j+2), h);
            std::complex<double> w3 = rev_cispi(4*(4*j+3), h);
            w2revpackedtab[3*j+h-1].real(packed<double,4>(w0.real(), w1.real(), w2.real(), w3.real()));
            w2revpackedtab[3*j+h-1].imag(packed<double,4>(w0.imag(), w1.imag(), w2.imag(), w3.imag()));
        }
    }

    if (wpowtab.size() < l)
    {
        wpowtab.resize(l);
        wpowtab[0] = 1;
        for (ulong k = 1; k <= depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
            wpowtab[pow2(k - 1) + j] = cispi_frac(j, pow2(k - 1));
    }

    if (wtwisttab.size() < l)
    {
        wtwisttab.resize(l);
        wtwisttab[0] = std::complex<double>(0); // should not be used?
        for (ulong k = 1; k < depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
        {
            std::complex<double> z0 = cispi_frac(4*j + 0, pow2(k + 1));
            std::complex<double> z1 = cispi_frac(4*j + 1, pow2(k + 1));
            std::complex<double> z2 = cispi_frac(4*j + 2, pow2(k + 1));
            std::complex<double> z3 = cispi_frac(4*j + 3, pow2(k + 1));
            wtwisttab[pow2(k - 1) + j] = complex_packed<double,4>(
                packed<double,4>(z0.real(), z1.real(), z2.real(), z3.real()),
                packed<double,4>(z0.imag(), z1.imag(), z2.imag(), z3.imag()));
        }
    }

    /////////////////////// three ////////////////////

    if (wthreerev1tab.size() < l)
    {
        wthreerev1tab.resize(l);
        wthreerev1tab[0] = 1;
        for (ulong k = 1; k <= depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
            wthreerev1tab[pow2(k - 1) + j] = cispi_frac(3*n_revbin(j,k-1)+1, 3*pow2(k-1));
    }

    if (wthreerev2tab.size() < l)
    {
        wthreerev2tab.resize(l);
        wthreerev2tab[0] = 1;
        for (ulong k = 1; k <= depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
            wthreerev2tab[pow2(k - 1) + j] = cispi_frac(3*n_revbin(j,k-1)-1, 3*pow2(k-1));
    }

    if (wthreepowtab.size() < 2*l)
    {
        wthreepowtab.resize(2*l);
        wthreepowtab[0].real(1.0);
        wthreepowtab[0].imag(0.0);
        for (ulong k = 0; k <= depth; k++)
        for (ulong j = 0; j < pow2(k); j++)
        {
            std::complex<double> t0 = cispi_frac(2*(4*j+0), 3*pow2(k+2));
            std::complex<double> t1 = cispi_frac(2*(4*j+1), 3*pow2(k+2));
            std::complex<double> t2 = cispi_frac(2*(4*j+2), 3*pow2(k+2));
            std::complex<double> t3 = cispi_frac(2*(4*j+3), 3*pow2(k+2));
            wthreepowtab[pow2(k) + j].real(packed<double,4>(t0.real(),t1.real(),t2.real(),t3.real()));
            wthreepowtab[pow2(k) + j].imag(packed<double,4>(t0.imag(),t1.imag(),t2.imag(),t3.imag()));
        }
    }

    ///////////////////// five /////////////////////////

    if (wfiverev1tab.size() < l)
    {
        wfiverev1tab.resize(l);
        wfiverev1tab[0] = 1;
        for (ulong k = 1; k <= depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
            wfiverev1tab[pow2(k - 1) + j] = cispi_frac(5*n_revbin(j,k-1)+1, 5*pow2(k-1));
    }

    if (wfiverev4tab.size() < l)
    {
        wfiverev4tab.resize(l);
        wfiverev4tab[0] = 1;
        for (ulong k = 1; k <= depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
            wfiverev4tab[pow2(k - 1) + j] = cispi_frac(5*n_revbin(j,k-1)-1, 5*pow2(k-1));
    }

    if (wfiverev2tab.size() < l)
    {
        wfiverev2tab.resize(l);
        wfiverev2tab[0] = 1;
        for (ulong k = 1; k <= depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
            wfiverev2tab[pow2(k - 1) + j] = cispi_frac(5*n_revbin(j,k-1)+2, 5*pow2(k-1));
    }

    if (wfiverev3tab.size() < l)
    {
        wfiverev3tab.resize(l);
        wfiverev3tab[0] = 1;
        for (ulong k = 1; k <= depth; k++)
        for (ulong j = 0; j < pow2(k - 1); j++)
            wfiverev3tab[pow2(k - 1) + j] = cispi_frac(5*n_revbin(j,k-1)-2, 5*pow2(k-1));
    }

    if (wfive1powtab.size() < 2*l)
    {
        wfive1powtab.resize(2*l);
        wfive1powtab[0].real(1.0);
        wfive1powtab[0].imag(0.0);
        for (ulong k = 0; k <= depth; k++)
        for (ulong j = 0; j < pow2(k); j++)
        {
            std::complex<double> t0 = cispi_frac(2*(4*j+0), 5*pow2(k+2));
            std::complex<double> t1 = cispi_frac(2*(4*j+1), 5*pow2(k+2));
            std::complex<double> t2 = cispi_frac(2*(4*j+2), 5*pow2(k+2));
            std::complex<double> t3 = cispi_frac(2*(4*j+3), 5*pow2(k+2));
            wfive1powtab[pow2(k) + j].real(packed<double,4>(t0.real(),t1.real(),t2.real(),t3.real()));
            wfive1powtab[pow2(k) + j].imag(packed<double,4>(t0.imag(),t1.imag(),t2.imag(),t3.imag()));
        }
    }

    if (wfive2powtab.size() < 2*l)
    {
        wfive2powtab.resize(2*l);
        wfive2powtab[0].real(1.0);
        wfive2powtab[0].imag(0.0);
        for (ulong k = 0; k <= depth; k++)
        for (ulong j = 0; j < pow2(k); j++)
        {
            std::complex<double> t0 = cispi_frac(2*(4*j+0)*2, 5*pow2(k+2));
            std::complex<double> t1 = cispi_frac(2*(4*j+1)*2, 5*pow2(k+2));
            std::complex<double> t2 = cispi_frac(2*(4*j+2)*2, 5*pow2(k+2));
            std::complex<double> t3 = cispi_frac(2*(4*j+3)*2, 5*pow2(k+2));
            wfive2powtab[pow2(k) + j].real(packed<double,4>(t0.real(),t1.real(),t2.real(),t3.real()));
            wfive2powtab[pow2(k) + j].imag(packed<double,4>(t0.imag(),t1.imag(),t2.imag(),t3.imag()));
        }
    }
}
