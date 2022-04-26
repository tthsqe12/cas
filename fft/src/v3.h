#pragma once

#include <cassert>
#include "packed.h"
#include "misc.h"

struct cpd_fft_ctx
{
    my_vector<complex_packed<double,4>> cpd_buffer;
    complex_packed<double,4>* complex_packed_double_buffer(ulong n)
    {
        if (cpd_buffer.size() < n)
            cpd_buffer.resize(n);
        return cpd_buffer.data();
    }

    std::vector<std::complex<double>> cd_buffer;
    std::complex<double>* complex_double_buffer(ulong n)
    {
        if (cd_buffer.size() < n)
            cd_buffer.resize(n);
        return cd_buffer.data();
    }

    std::vector<double> d_buffer;
    double* double_buffer(ulong n)
    {
        if (d_buffer.size() < n)
            d_buffer.resize(n);
        return d_buffer.data();
    }

    void fit_wtab(ulong depth);

    // exp(2*pi*i*[0/1, 1/2, 1/4,3/4, 1/8,5/8,3/8,7/18, 1/16,9/16,5/16,13/16,3/16,11/16,7/16,15/16, ...])
    std::vector<std::complex<double>>wrevtab;
    std::complex<double> wrev(ulong j) {
        return wrevtab[j];
    }
    std::complex<double> wrevinv(ulong j) {
        // could also return conj(wtab.at(j))
        return wrevtab[j ^ (saturate_bits(j) >> 1)];
    }
    const std::complex<double>* wrev_table() {
        return wrevtab.data();
    }

    // exp(2*pi*i*[0/1,      1/4,     1/8,    3/8,     1/16,     5/16,      3/16,      7/16,       ...])
    std::vector<std::complex<double>>w2revtab;
    std::complex<double> w2rev(ulong j) {
        return w2revtab[j];
    }
    const std::complex<double>* w2rev_table() {
        return w2revtab.data();
    }

    // w2rev(2*(4*j+{0,1,2,3}))^h
    my_vector<complex_packed<double,4>> w2revpackedtab;
    complex_packed<double,4> w2revpacked(ulong j, ulong h) {
        assert(1 <= h && h <= 3);
        assert(3*j + h - 1 < w2revpackedtab.size());
        return w2revpackedtab[3*j + h - 1];
    }

    // could also do a w4revtab since exp(2*pi*i*(1/4)) is trivial


    // exp(2*pi*i*j/2^k)
    // 0/1, 0/2, 0/4,1/4, 0/8,1/8,2/8,3/8, 0/16,1/16,2/16,3/16,4/16,5/16,6/16,7/16, ...
    std::vector<std::complex<double>> wpowtab;
    std::complex<double> wpow(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k - 1));
        return wpowtab[pow2(k - 1) + j];
    }
    const std::complex<double>* wpow_table(ulong k) {
        assert(k > 0);
        return wpowtab.data() + pow2(k - 1);
    }

    // exp(2*pi*i*(j+(0,1,2,3)/4)/2^k)  0<=j<2^(k-1)
    my_vector<complex_packed<double,4>> wtwisttab;
    const complex_packed<double,4> wtwist(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k - 1));
        return wtwisttab[pow2(k - 1) + j];
    }
    const complex_packed<double,4>* wtwist_table(ulong k) {
        assert(k > 0);
        return wtwisttab.data() + pow2(k - 1);
    }

    ////////////////////// three ////////////////////

    // exp(2*pi*i*(3*n_revbits(j,k-1)+1)/(3*2^k)) for 0 <= j < 2^(k-1)
    std::vector<std::complex<double>> wthreerev1tab;
    std::complex<double> wthreerev1(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k-1));
        return wthreerev1tab.at(pow2(k - 1) + j);
    }
    std::complex<double>* wthreerev1(ulong k) {
        assert(k > 0);
        return wthreerev1tab.data() + pow2(k-1);
    }

    // exp(2*pi*i*(3*n_revbits(j,k-1)-1)/(3*2^k)) for 0 <= j < 2^(k-1)
    std::vector<std::complex<double>> wthreerev2tab;
    std::complex<double> wthreerev2(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k-1));
        return wthreerev2tab[pow2(k - 1) + j];
    }
    std::complex<double>* wthreerev2(ulong k) {
        assert(k > 0);
        return wthreerev2tab.data() + pow2(k-1);
    }

    // exp(2*pi*i*(j+(0,1,2,3)/4)/(3*2^k))
    my_vector<complex_packed<double,4>> wthreepowtab;
    complex_packed<double,4> wthreepow(ulong k, ulong j) {
        assert(j < pow2(k));
        return wthreepowtab[pow2(k) + j];
    }
    complex_packed<double,4>* wthreepow_table(ulong k) {
        return wthreepowtab.data() + pow2(k);
    }

    ////////////////////// five ////////////////////

    // exp(2*pi*i*(5*n_revbits(j,k-1)+1)/(5*2^k)) for 0 <= j < 2^(k-1)
    std::vector<std::complex<double>> wfiverev1tab;
    std::complex<double> wfiverev1(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k-1));
        return wfiverev1tab.at(pow2(k - 1) + j);
    }
    std::complex<double>* wfiverev1(ulong k) {
        assert(k > 0);
        return wfiverev1tab.data() + pow2(k-1);
    }

    // exp(2*pi*i*(5*n_revbits(j,k-1)-1)/(5*2^k)) for 0 <= j < 2^(k-1)
    std::vector<std::complex<double>> wfiverev4tab;
    std::complex<double> wfiverev4(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k-1));
        return wfiverev4tab.at(pow2(k - 1) + j);
    }
    std::complex<double>* wfiverev4(ulong k) {
        assert(k > 0);
        return wfiverev4tab.data() + pow2(k-1);
    }

    // exp(2*pi*i*(5*n_revbits(j,k-1)+2)/(5*2^k)) for 0 <= j < 2^(k-1)
    std::vector<std::complex<double>> wfiverev2tab;
    std::complex<double> wfiverev2(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k-1));
        return wfiverev2tab.at(pow2(k - 1) + j);
    }
    std::complex<double>* wfiverev2(ulong k) {
        assert(k > 0);
        return wfiverev2tab.data() + pow2(k-1);
    }

    // exp(2*pi*i*(5*n_revbits(j,k-1)-2)/(5*2^k)) for 0 <= j < 2^(k-1)
    std::vector<std::complex<double>> wfiverev3tab;
    std::complex<double> wfiverev3(ulong k, ulong j) {
        assert(k > 0);
        assert(j < pow2(k-1));
        return wfiverev3tab.at(pow2(k - 1) + j);
    }
    std::complex<double>* wfiverev3(ulong k) {
        assert(k > 0);
        return wfiverev3tab.data() + pow2(k-1);
    }

    // exp(2*pi*i*(j+(0,1,2,3)/4)/(5*2^k))
    my_vector<complex_packed<double,4>> wfive1powtab;
    complex_packed<double,4> wfive1pow(ulong k, ulong j) {
        assert(j < pow2(k));
        return wfive1powtab[pow2(k) + j];
    }
    complex_packed<double,4>* wfive1pow_table(ulong k) {
        return wfive1powtab.data() + pow2(k);
    }

    // exp(2*pi*i*2*(j+(0,1,2,3)/4)/(5*2^k))
    my_vector<complex_packed<double,4>> wfive2powtab;
    complex_packed<double,4> wfive2pow(ulong k, ulong j) {
        assert(j < pow2(k));
        return wfive2powtab[pow2(k) + j];
    }
    complex_packed<double,4>* wfive2pow_table(ulong k) {
        return wfive2powtab.data() + pow2(k);
    }

    ////////////////////////////////////

    void fft_full_block(complex_packed<double,4>* x, ulong k, ulong j);
    void fft_full(complex_packed<double,4>* x, ulong k, ulong j);
    void fft_base(complex_packed<double,4>* x, ulong j);
    void fft_trunc_block(complex_packed<double,4>* x, ulong k, ulong j, ulong itrunc);
    void fft_trunc(complex_packed<double,4>* x, ulong k, ulong itrunc);
    void ifft_base(complex_packed<double,4>* x);
    void ifft_full(complex_packed<double,4>* x, ulong k);
    void ifft_full_block(complex_packed<double,4>* x, ulong );

    void fft_three_trunc(complex_packed<double,4>* x, ulong k, ulong itrunc);
    void ifft_three_full(complex_packed<double,4>* x, ulong k);

    void fft_five_trunc(complex_packed<double,4>* x, ulong k, ulong itrunc);
    void ifft_five_full(complex_packed<double,4>* x, ulong k);

    double my_mpn_mul(ulong* z, const ulong* a, ulong an, const ulong* b, ulong bn, ulong bits, ulong m);
    void my_mpn_mul(ulong* z, const ulong* a, ulong an, const ulong* b, ulong bn);
};
