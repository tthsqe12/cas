#include "packed.h"
#include "misc.h"
#include "v1.h"
#include "v2.h"
#include "v3.h"

void record_time(
    bool first,
    double len,
    double t1, double t2,
    bool for_latex,
    std::vector<std::pair<double,double>>& v1,
    std::vector<std::pair<double,double>>& v2)
{
    t1 *= 1e6;
    t2 *= 1e6;
    double depth = log2(len);
    if (for_latex)
    {
        v1.emplace_back(depth, 5*len*depth/t1);
        v2.emplace_back(depth, 5*len*depth/t2);
    }
    //else
    {
        std::cout << (first ? "depth" : "     ")
            << format_fixed(depth,3,1) << ":"
            << format_fixed(5*len*depth/t1,4,2)
            << format_fixed(5*len*depth/t2,4,2)
            << std::endl;
    }
}

std::ostream& operator<<(std::ostream& o, const std::vector<std::pair<double,double>>& v)
{
    bool first = true;
    for (auto& vi : v)
    {
        if (!first)
            std::cout << ",";
        first = false;
        std::cout << "(" << vi.first << "," << vi.second << ")";
    }
    return o;
}

#define EXP_FAC 16

void profile_v1_fft(pd_fft_ctx<4>& Q, ulong min_depth, ulong max_depth, bool for_latex)
{
    std::vector<std::pair<double,double>> fft_time, ifft_time;
    ulong t1, t2, t3;
    min_depth = std::max(min_depth, ulong(5));
    max_depth = std::max(max_depth, min_depth + 1);

    t1 = GetMS();
    Q.fit_wrevtab(max_depth);
    t2 = GetMS();
    std::cout << "setup: " << t2 - t1 << " ms" << std::endl;

    packed<double,4>* X = new packed<double,4>[pow2(max_depth)];

    for (ulong L = min_depth; L < max_depth; L++)
    {
        bool first = true;
        for (ulong len = pow2(L); len < 2*pow2(L); len = round_up(1+(EXP_FAC+1)*len/EXP_FAC,16))
        {
            ulong k = clog2(len);
            for (ulong i = 0; i < pow2(k); i++)
                X[i] = 0.0;
            ulong nreps = 1 + 100000000/(len*clog2(len));
            t1 = GetMS();
            for (ulong rep = 0; rep < nreps; rep++)
                Q.fft_trunc(X, k, roundup(len/2,16), len);
            t2 = GetMS();
            for (ulong rep = 0; rep < nreps; rep++)
                Q.ifft_trunc(X, k, 0, len, len, false);
            t3 = GetMS();

            record_time(first, len, double(t2-t1)/(4*nreps), double(t3-t2)/(4*nreps), for_latex, fft_time, ifft_time);
            first = false;
        }
    }

    delete[] X;

    if (for_latex)
    {
        std::cout << fft_time << std::endl;
        std::cout << ifft_time << std::endl;
    }
}

void profile_v2_fft(fftv2_ctx& Q, ulong minL, ulong maxL, bool for_latex)
{
    std::vector<std::pair<double,double>> fft_time, ifft_time;
    ulong t1, t2, t3;
    minL = std::max(minL, ulong(10));
    maxL = std::max(maxL, minL + 1);

    t1 = GetMS();
    Q.set_depth(maxL);
    t2 = GetMS();
    std::cout << "setup: " << t2-t1 << " ms " << std::endl;

    minL = std::max(minL, ulong(LG_BLK_SZ + 1));
    for (ulong L = minL+1; L <= maxL; L++)
    {
        Q.set_depth(L);
        Q.set_data(reinterpret_cast<double*>(my_aligned_alloc(4096, Q.data_size()*sizeof(double))));

        // do 1/2*2^L < otrunc <= 2^L
        ulong otrunc = pow2(L-1);
        while (otrunc < pow2(L))
        {
            otrunc = round_up(1+1+(EXP_FAC+1)*otrunc/EXP_FAC, Q.blk_sz);
            otrunc = std::min(otrunc, pow2(L));
            ulong itrunc = round_up(otrunc/2, Q.blk_sz);
            for (ulong i = 0; i < pow2(L); i++)
                Q.set_index(i, 0);
            ulong nreps = 1 + 300000000/(otrunc*clog2(otrunc));
            t1 = GetMS();
            for (ulong i = 0; i < nreps; i++)
                Q.fft_trunc(itrunc, otrunc);
            t2 = GetMS();
            for (ulong i = 0; i < nreps; i++)
                Q.ifft_trunc(otrunc);
            t3 = GetMS();

            record_time(otrunc == pow2(L), otrunc, double(t2-t1)/(nreps), double(t3-t2)/(nreps), for_latex, fft_time, ifft_time);
        }

        delete[] Q.release_data();
    }

    if (for_latex)
    {
        std::cout << fft_time << std::endl;
        std::cout << ifft_time << std::endl;
    }
}

void profile_v3_fft(cpd_fft_ctx& Q, ulong min_depth, ulong max_depth, bool for_latex)
{
    std::vector<std::pair<double,double>> fft_time, ifft_time;
    ulong t1, t2, t3;
    min_depth = std::max(min_depth, ulong(5));
    max_depth = std::max(max_depth, min_depth + 1);

    t1 = GetMS();
    Q.fit_wtab(max_depth);
    t2 = GetMS();
    std::cout << "setup: " << t2 - t1 << " ms" << std::endl;

    complex_packed<double,4>* A = new complex_packed<double,4>[pow2(max_depth)];

    for (ulong k = min_depth; k < max_depth; k++)
    {
        ulong nreps = 1 + 300000000/(k<<k);

        ulong len = ulong(1) << k;
        for (ulong i = 0; i < len; i++)
        {
            A[i].real(packed<double,4>(0.0));
            A[i].imag(packed<double,4>(0.0));
        }
        t1 = GetMS();
        for (ulong rep = 0; rep < nreps; rep++)
            Q.fft_trunc(A, k, len);
        t2 = GetMS();
        for (ulong rep = 0; rep < nreps; rep++)
            Q.ifft_full(A, k);
        t3 = GetMS();
        record_time(true, len, double(t2-t1)/(nreps), double(t3-t2)/(nreps), for_latex, fft_time, ifft_time);

        len = ulong(5) << (k-2);
        for (ulong i = 0; i < len; i++)
        {
            A[i].real(packed<double,4>(0.0));
            A[i].imag(packed<double,4>(0.0));
        }
        t1 = GetMS();
        for (ulong rep = 0; rep < nreps; rep++)
            Q.fft_five_trunc(A, k-2, len);
        t2 = GetMS();
        for (ulong rep = 0; rep < nreps; rep++)
            Q.ifft_five_full(A, k-2);
        t3 = GetMS();
        record_time(false, len, double(t2-t1)/(nreps), double(t3-t2)/(nreps), for_latex, fft_time, ifft_time);

        len = ulong(3) << (k-1);
        for (ulong i = 0; i < len; i++)
        {
            A[i].real(packed<double,4>(0.0));
            A[i].imag(packed<double,4>(0.0));
        }
        t1 = GetMS();
        for (ulong rep = 0; rep < nreps; rep++)
            Q.fft_three_trunc(A, k-1, len);
        t2 = GetMS();
        for (ulong rep = 0; rep < nreps; rep++)
            Q.ifft_three_full(A, k-1);
        t3 = GetMS();
        record_time(false, len, double(t2-t1)/(nreps), double(t3-t2)/(nreps), for_latex, fft_time, ifft_time);
    }

    delete[] A;

    if (for_latex)
    {
        std::cout << fft_time << std::endl;
        std::cout << ifft_time << std::endl;
    }
}


template <class T>
void print_mul_data(T& Q, int min_lg_bit_size, int max_lg_bit_size, bool for_latex)
{
    bool first = true;
    min_lg_bit_size = std::max(min_lg_bit_size, 12);
    for (int lg_bit_size = min_lg_bit_size; lg_bit_size < max_lg_bit_size; lg_bit_size++)
    {
        ulong t1, t2;
        ulong max_len = pow2(1+lg_bit_size-6);
        ulong* data = new ulong[max_len*2];
        for (ulong i = 0; i < max_len; i++)
            data[i] = -ulong(1+i);
        t1 = GetMS();
        Q.my_mpn_mul(data+max_len, data+max_len/2, max_len/2, data+max_len/2, max_len/2);
        t2 = GetMS();
        if (!for_latex)
            std::cout << "words " << max_len << " precomp: " << t2-t1 << std::endl;
        for (int ci = 0; ci < 10; ci++)
        {
            ulong cn = max_len/2*pow(2.0, ci*1.0/10);
            ulong cbits = cn*64;
            double lgcbits = log2(cbits);
            ulong nsamples = 0;
            double total_time = 0;
            double max_time = 0;
            double min_time = 1.0e100;
            for (ulong an = (cn+1)/2; an <= 3*cn/4; an += 1+an/12)
            {            
                ulong bn = cn - an;
                if (!(bn <= an && an < cn))
                    continue;
                nsamples++;
                ulong* a = data;
                ulong* b = data + an;
                ulong* c = data + max_len;
                ulong nreps = 1 + 30000000/(cn*clog2(cn));
                t1 = GetMS();
                for (ulong rep = 0; rep < nreps; rep++)
                    Q.my_mpn_mul(c, a, an, b, bn);
                t2 = GetMS();
                double time = double(t2 - t1)*1e8/(lgcbits*cbits*nreps);
                total_time += time;
                min_time = std::min(min_time, time);
                max_time = std::max(max_time, time);
            }

            if (nsamples < 1)
                continue;

            if (for_latex)
            {
                if (!first)
                    std::cout << ",";
                std::cout << "(" << lgcbits << "," << min_time << "),";
                std::cout << "(" << lgcbits << "," << total_time/nsamples << "),";
                std::cout << "(" << lgcbits << "," << max_time << ")";
                std::cout << std::flush;
            }
            else
            {
                std::cout << format_fixed(lgcbits,2, 2) << ":"
                          << format_fixed(total_time/nsamples, 3, 3)
                          << format_fixed(min_time, 4, 3)
                          << format_fixed(max_time, 4, 3)
                          << std::endl;
            }
            first = false;
        }

        delete[] data;
    }
}


struct gmp_ctx {
    void my_mpn_mul(ulong* z, const ulong* a, ulong an, const ulong* b, ulong bn) {
        mpn_mul(z, a, an, b, bn);
    }
};

int main(int argc, char *argv[])
{
#if 0
    if(1){
        mpn_fft_ctx Q;
        std::cout << "--- v1 fft/ifft gflops ---" << std::endl;
        profile_v1_fft(Q.ctx, 10, 20, true);
    }

    if(1){
        fftv2_ctx Q(0x0003f00000000001);
        std::cout << "--- v2 fft/ifft gflops ---" << std::endl;
        profile_v2_fft(Q, 10, 20, true);
    }

    if(1){
        cpd_fft_ctx Q;
        std::cout << "--- v3 fft/ifft gflops ---" << std::endl;
        profile_v3_fft(Q, 10, 20, true);
    }
#endif

#if 0
    if(0){
        cpd_fft_ctx Q;
        print_mul_data<cpd_fft_ctx>(Q, 12, 17, false);
    }

    if(0){
        mpn_fft_ctx Q;
        print_mul_data<mpn_fft_ctx>(Q, 12, 25, false);
    }
#endif

#if 1
    char id[] = "coffee";
    if(1){
        gmp_ctx Q;
        std::cout << "\\def\\" << id << "gmp{";
        print_mul_data<gmp_ctx>(Q, 12, 31, true);
        std::cout << "}" << std::endl << std::endl;
    }
    if(1){
        cpd_fft_ctx Q;
        std::cout << "\\def\\" << id << "cpd{";
        print_mul_data<cpd_fft_ctx>(Q, 12, 21, true);
        std::cout << "}" << std::endl << std::endl;
    }
    if(1){
        mpn_fft_ctx Q;
        std::cout << "\\def\\" << id << "pd{";
        print_mul_data<mpn_fft_ctx>(Q, 12, 28, true);
        std::cout << "}" << std::endl << std::endl;
    }
    if(1){
        mpn_ctx_v2 Q(0x0003f00000000001);
        std::cout << "\\def\\" << id << "sd{";
        print_mul_data<mpn_ctx_v2>(Q, 14, 31, true);
        std::cout << "}" << std::endl << std::endl;
    }
#endif

#if 0
    char id[] = "bull";
    if(0){
        gmp_ctx Q;
        std::cout << "\\def\\" << id << "gmp{";
        print_mul_data<gmp_ctx>(Q, 12, 29, true);
        std::cout << "}" << std::endl << std::endl;
    }
    if(0){
        cpd_fft_ctx Q;
        std::cout << "\\def\\" << id << "cpd{";
        print_mul_data<cpd_fft_ctx>(Q, 12, 21, true);
        std::cout << "}" << std::endl << std::endl;
    }
    if(1){
        mpn_fft_ctx Q;
        std::cout << "\\def\\" << id << "pd{";
        print_mul_data<mpn_fft_ctx>(Q, 12, 17, true);
        std::cout << "}" << std::endl << std::endl;
    }
    if(1){
        mpn_ctx_v2 Q(0x0003f00000000001);
        std::cout << "\\def\\" << id << "sd{";
        print_mul_data<mpn_ctx_v2>(Q, 14, 29, true);
        std::cout << "}" << std::endl << std::endl;
    }
#endif

    flint_cleanup();
    return 0;
}
