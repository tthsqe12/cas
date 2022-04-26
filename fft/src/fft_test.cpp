#include "packed.h"
#include "misc.h"
#include "v1.h"
#include "v2.h"
#include "v3.h"


template <class T>
void test_mul(T& Q, ulong minsize, ulong maxsize, ulong nreps, flint_rand_t state)
{
    ulong* a = new ulong[maxsize];
    ulong* b = new ulong[maxsize];
    ulong* c = new ulong[maxsize];
    ulong* d = new ulong[maxsize];

    ulong dgs = n_sizeinbase(nreps, 10);

    std::cout << "mul " << format_fixed(0, dgs) << "/" << format_fixed(nreps, dgs) << std::flush;

    minsize = std::max(ulong(10), minsize);
    maxsize = std::max(minsize, maxsize);
    for (ulong rep = 0; rep < nreps; rep++)
    {
        ulong an = 2 + n_randint(state, maxsize - 4);
        ulong bn = 1 + n_randint(state, std::min(an, maxsize - an));

        for (ulong i = 0; i < maxsize; i++)
        {
            a[i] = n_randlimb(state);
            b[i] = n_randlimb(state);
            c[i] = n_randlimb(state);
            d[i] = n_randlimb(state);
        }

        for (ulong ii = 0; ii < 2*dgs + 5; ii++)
            std::cout << '\b';

        std::cout << "mul " << format_fixed(1+rep, dgs) << "/" << format_fixed(nreps, dgs) << std::flush;

        Q.my_mpn_mul(d, a, an, b, bn);
        mpn_mul(c, a, an, b, bn);
        for (ulong i = 0; i < an + bn; i++)
        {
            if (c[i] != d[i])
            {
                std::cout << std::endl;
                std::cout << "FAILED " << "an: " << an << ", bn: " << bn << std::endl;
                std::cout << "limb[" << i << "] = 0x" << format_hex(d[i]) << " should be 0x" << format_hex(c[i]) << std::endl;
                std::abort();
            }
        }
    }

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    std::cout << std::endl;
}

bool output_length_is_ok(
    cpd_fft_ctx& Q,
    ulong bits, ulong m,
    ulong zn,
    ulong* a, ulong* b, ulong* c, ulong* d,
    flint_rand_t state)
{
    if (zn < 2)
        return false;
    ulong bn = zn/2;
    ulong an = zn - bn;
    while (bn >= zn/8)
    {
        for (ulong i = 0; i < an; i++)
            a[i] |= n_randlimb(state);
        for (ulong i = 0; i < bn; i++)
            b[i] |= n_randlimb(state);

        double error = Q.my_mpn_mul(c, a, an, b, bn, bits, m);
        if (!(error < 0.25))
            return false;
        mpn_mul(d, a, an, b, bn);
        if (mpn_cmp(d, c, an + bn) != 0)
            return false;
        ulong step = 1 + zn/32;
        if (bn <= step)
            return true;
        bn -= step;
        an += step;

        a[n_randint(state, an)] += n_randlimb(state);
        b[n_randint(state, bn)] += n_randlimb(state);
        a[n_randint(state, an)] &= n_randlimb(state);
        b[n_randint(state, bn)] &= n_randlimb(state);
    }

    return true;
}

ulong find_failure_point(cpd_fft_ctx& Q, ulong bits, ulong m, flint_rand_t state)
{
    ulong minlen = 10;
    ulong maxlen = bits >= 19 ? 2000 : bits >= 16 ? 60000 : 200000;

    ulong* a = new ulong[maxlen];
    ulong* b = new ulong[maxlen];
    ulong* c = new ulong[maxlen];
    ulong* d = new ulong[maxlen];

    ulong tmax = maxlen;

    // bsplit a few times for security
    for (int i = 0; i < 10; i++)
    {
        for (ulong j = 0; j < maxlen; j++)
        {
            a[j] = -ulong(2);
            b[j] = -ulong(1);
        }

        ulong tmin = 0;
        while (tmax > tmin + 16)
        {
            ulong tmid = (2*tmax+(2+i)*tmin)/(4+i);
            if (output_length_is_ok(Q, bits, m, tmid, a, b, c, d, state))
                tmin = tmid;
            else
                tmax = tmid;
        }
        for (; tmin < tmax; tmin++)
        {
            if (!output_length_is_ok(Q, bits, m, tmin, a, b, c, d, state))
                break;
        }
        tmax = tmin;
    }
    return tmax;
}


void test_v1_fft(pd_fft_ctx<4>& ctx, ulong minL, ulong maxL, ulong ireps, flint_rand_t state)
{
    ulong irepmul = 40;
    minL = std::max(minL, ulong(3));

    packed<double,4>* X = new packed<double,4>[pow2(maxL)];
    packed<double,4>* F = new packed<double,4>[pow2(maxL)];

    ulong dgs1 = n_sizeinbase(maxL, 10);
    ulong dgs2 = n_sizeinbase(ireps + irepmul*maxL, 10);

    std::cout << "fft " << format_fixed(0, dgs1) << "." << format_fixed(0, dgs2)
              << "/" << format_fixed(maxL, dgs1) << "." << format_fixed(0, dgs2) << std::flush;

    for (ulong L = minL; L <= maxL; L++)
    {
        ulong i;
        ulong Xn = pow2(L);
        ulong nreps = ireps + irepmul*L;
        for (ulong rep = 0; rep < nreps; rep++)
        {
            for (ulong ii = 0; ii < 2*(dgs1 + dgs2) + 7; ii++)
                std::cout << '\b';
            std::cout << "fft " << format_fixed(L, dgs1) << "." << format_fixed(1+rep, dgs2)
                      << "/" << format_fixed(maxL, dgs1) << "." << format_fixed(nreps, dgs2) << std::flush;

            // randomize input data
            for (i = 0; i < Xn; i++)
                F[i] = X[i] = packed<double,4>(i+0, i+1, i+2, i+3);

            // output of fft_trunc is supposed to be eval_poly
            ulong itrunc = round_up(1 + n_randint(state, Xn), 4);
            ulong otrunc = round_up(1 + n_randint(state, Xn), 4);
            ctx.fft_trunc(F, L, itrunc, otrunc);

            for (int check_reps = 0; check_reps < 3; check_reps++)
            {
                i = n_randint(state, otrunc);
                packed<double,4> y = ctx.w2rev(i/2);
                y = (i&1) ? neg(y) : y;
                y = eval_poly_mod(X, itrunc, y, ctx.n, ctx.ninv);
                if (!equal_mod(y, F[i], ctx.n, ctx.ninv))
                {
                    std::cout << std::endl;
                    std::cout << "fft error at index " << i << std::endl;
                    std::cout << "itrunc: " << itrunc << std::endl;
                    std::cout << "otrunc: " << otrunc << std::endl;
                    std::abort();
                }
            }

            // output of ifft_trunc is supposed to be 2^L*input
            ulong trunc = round_up(1 + n_randint(state, Xn), 4);

            for (i = 0; i < trunc; i++)
                F[i] = X[i];

            ctx.fft_trunc(F, L, trunc, trunc);
            ctx.ifft_trunc(F, L, 0, trunc, trunc, false);

            for (int check_reps = 0; check_reps < 3; check_reps++)
            {
                i = n_randint(state, trunc);
                packed<double,4> pow2l = 1.0;
                for (ulong ii = 0; ii < L; ii++)
                    pow2l = reduce_to_pm1n(add(pow2l, pow2l), ctx.n, ctx.ninv);
                if (!equal_mod(F[i], mulmod2(X[i], pow2l, ctx.n, ctx.ninv), ctx.n, ctx.ninv))
                {
                    std::cout << std::endl;
                    std::cout << "ifft error at index " << i << std::endl;
                    std::cout << "trunc: " << trunc << std::endl;
                    std::abort();
                }
            }
        }
    }

    delete[] F;
    delete[] X;

    std::cout << std::endl;
}

void test_v2_fft(fftv2_ctx& ctx, ulong minL, ulong maxL, ulong ireps, flint_rand_t state)
{
    ulong irepmul = 10;
    minL = std::max(minL, ulong(LG_BLK_SZ));

    ulong dgs1 = n_sizeinbase(maxL, 10);
    ulong dgs2 = n_sizeinbase(ireps + irepmul*maxL, 10);

    std::cout << "fft " << format_fixed(0, dgs1) << "." << format_fixed(0, dgs2)
              << "/" << format_fixed(maxL, dgs1) << "." << format_fixed(0, dgs2) << std::flush;

    for (ulong L = minL; L <= maxL; L++)
    {
        ulong i;
        ulong Xn = pow2(L);
        double* X = new double[Xn];

        ctx.set_depth(L);
        ctx.set_data(reinterpret_cast<double*>(my_aligned_alloc(32, ctx.data_size()*sizeof(double))));

        ulong nreps = ireps + irepmul*L;
        for (ulong rep = 0; rep < nreps; rep++)
        {
            for (ulong ii = 0; ii < 2*(dgs1 + dgs2) + 7; ii++)
                std::cout << '\b';
            std::cout << "fft " << format_fixed(L, dgs1) << "." << format_fixed(1+rep, dgs2)
                      << "/" << format_fixed(maxL, dgs1) << "." << format_fixed(nreps, dgs2) << std::flush;

            // randomize input data
            for (i = 0; i < Xn; i++)
                X[i] = i + 1;

            // output of fft_trunc is supposed to be eval_poly
            ulong itrunc = round_up(1 + n_randint(state, Xn), ctx.blk_sz);
            ulong otrunc = round_up(1 + n_randint(state, Xn) , ctx.blk_sz);

            for (i = 0; i < itrunc; i++)
                ctx.set_index(i, X[i]);

            ctx.fft_trunc(itrunc, otrunc);

            for (int check_reps = 0; check_reps < 3; check_reps++)
            {
                i = n_randint(state, otrunc);
                double y = eval_poly_mod(X, itrunc, (i&1) ? -ctx.w2s[i/2] : ctx.w2s[i/2], ctx.p, ctx.pinv);
                if (!equal_mod(y, ctx.get_index(i), ctx.p, ctx.pinv))
                {
                    std::cout << std::endl;
                    std::cout << "fft error at index " << i << std::endl;
                    std::cout << "itrunc: " << itrunc << std::endl;
                    std::cout << "otrunc: " << otrunc << std::endl;
                    std::abort();
                }
            }

            // output of ifft_trunc is supposed to be 2^L*input
            ulong trunc = round_up(1 + n_randint(state, Xn), ctx.blk_sz);
            for (i = 0; i < trunc; i++)
                ctx.set_index(i, X[i]);

            ctx.fft_trunc(trunc, trunc);
            ctx.ifft_trunc(trunc);

            for (int check_reps = 0; check_reps < 3; check_reps++)
            {
                i = n_randint(state, trunc);
                double m = reduce_0n_to_pmhn(nmod_pow_ui(2, L, ctx.mod), ctx.p);
                if (!equal_mod(ctx.get_index(i), mulmod2(X[i], m, ctx.p, ctx.pinv), ctx.p, ctx.pinv))
                {
                    std::cout << std::endl;
                    std::cout << "ifft error at index " << i << std::endl;
                    std::cout << "trunc: " << trunc << std::endl;
                    std::abort();
                }
            }
        }

        delete[] ctx.release_data();
        delete[] X;
    }

    std::cout << std::endl;
}


int main(int argc, char *argv[])
{
    flint_rand_t state;
    flint_randinit(state);

    if(1){
        std::cout << "------------- testing v1 ------------" << std::endl;
        mpn_fft_ctx Q;
        test_v1_fft(Q.ctx, 5, 16, 100, state);
        test_mul(Q, 10, 9000, 5000, state);
        std::cout << "PASSED" << std::endl;
    }
    if(1){
        std::cout << "------------- testing v2 ------------" << std::endl;
        mpn_ctx_v2 Q(0x0003f00000000001);
        test_v2_fft(Q.ffts[0], 14, 20, 20, state);
        test_mul(Q, 10, 30000, 5000, state);
        std::cout << "PASSED" << std::endl;
    }
    if(1){
        std::cout << "------------- testing v3 ------------" << std::endl;
        cpd_fft_ctx Q;
        test_mul(Q, 10, 9000, 5000, state);
        std::cout << "PASSED" << std::endl;

        std::cout << "         22      21      20      19      18      17      16      15" << std::endl;
        for (ulong mult = 1; mult <= 5; mult += 2)
        {
            std::cout << mult << ": " << std::flush;
            for (ulong bits = 22; bits >= 19; bits--)
            {
                flint_randclear(state);
                flint_randinit(state);
                std::cout << format_fixed(find_failure_point(Q, bits, mult, state), 8) << std::flush;
            }
            std::cout << std::endl;
        }
    }

    flint_randclear(state);
    flint_cleanup();
    return 0;
}
