void test_v2_trunc(ulong minL, ulong maxL)
{
    flint_rand_t state;
    flint_randinit(state);

std::cout << "------------- testing fft --------------- " << std::endl;

    fftv2_ctx ctx(0x03f00000000001ULL);

    minL = std::max(minL, ulong(LG_BLK_SZ));
    for (ulong L = minL; L <= maxL; L++)
    {
        std::cout << "-- depth: " << L << " --" << std::endl;

        ulong i;
        ulong Xn = pow2(L);
        double* X = new double[Xn];

        ctx.set_depth(L);
        ctx.set_data(new double[ctx.data_size()]);

        for (ulong rep = 0; rep < 300 + 100*L; rep++)
        {
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
                double y = ctx.eval_poly(X, itrunc, (i&1) ? -ctx.w2s[i/2] : ctx.w2s[i/2]);
                if (reduce_to_0n(y, ctx.p, ctx.pinv) !=
                    reduce_to_0n(ctx.get_index(i), ctx.p, ctx.pinv))
                {
                    std::cout << "fft error at index " << i << std::endl;
                    std::abort();
                }
            }

            // output of ifft_trunc is supposed to be 2^L*input
            ulong trunc = round_up(1 + n_randint(state, Xn), BLK_SZ);

            for (i = 0; i < trunc; i++)
                ctx.set_index(i, X[i]);

            ctx.fft_trunc(trunc, trunc);
            ctx.ifft_trunc(trunc);

            for (int check_reps = 0; check_reps < 3; check_reps++)
            {
                i = n_randint(state, trunc);
                double m = reduce_0n_to_pmhn(nmod_pow_ui(2, L, ctx.mod), ctx.p);
                if (reduce_to_0n(ctx.get_index(i), ctx.p, ctx.pinv) !=
                    reduce_to_0n(mulmod2(X[i], m, ctx.p, ctx.pinv), ctx.p, ctx.pinv))
                {
                    std::cout << "ifft error at index " << i << std::endl;
                    std::abort();
                }
            }
        }

        delete[] ctx.release_data();
        delete[] X;
    }

    flint_randclear(state);
}

void display_dpoints(std::vector<double> v) {
    std::cout << "{";
    for (ulong i = 0; i < v.size(); i += 2)
    {
        if (i > 0)
            std::cout << ", ";
        std::cout << "(" << format_fixed(v[i+0], 2, 2) << ", "
                         << format_fixed(v[i+1], 2, 2) << ")";
    }
    std::cout << "}" << std::endl;
}

void profile_v2_trunc(ulong minL, ulong maxL)
{
    double tmul = 10000000;
    timeit_t timer;
    fftv2_ctx ctx(0x03f00000000001ULL);

    std::vector<double>  fft_trunc_times;
    std::vector<double> ifft_trunc_times;
    std::vector<double>  fft_times;
    std::vector<double> ifft_times;

    double time;

std::cout << "------------- profiling fft --------------- " << std::endl;
std::cout << "*** reports time and (time/(n*log(n)) where n = truncated length ***" << std::endl;

    minL = std::max(minL, ulong(LG_BLK_SZ + 1));
    for (ulong L = minL; L <= maxL; L++)
    {
std::cout << "-- depth " << format_fixed(L, 2) << " -- | --  fft  -- | --  ifft  -- |" << std::endl;

        timeit_start(timer);
        ctx.set_depth(L);
        timeit_stop(timer);
        std::cout << "precomp: "
                  << format_fixed(timer->wall, 4)
                  << " ("
                  << format_fixed(timer->wall*tmul/pow2(L), 2, 2)
                  << ")"
                  << std::endl;

        ctx.set_data(new double[ctx.data_size()]);

        ulong nreps = 1;
        if (L < 24)
            nreps <<= (25 - L)/2;

        // do 1/2*2^L < otrunc <= 2^L
        ulong otrunc = pow2(L-1);
        while (true)
        {
            otrunc = round_up(otrunc + pow2(std::max(L, ulong(5))-5), ctx.blk_sz);
            if (otrunc > pow2(L))
                break;

            // any itrunc <= otrunc will do
            ulong itrunc = round_up(otrunc/2, ctx.blk_sz);

            for (ulong i = 0; i < itrunc; i++)
                ctx.set_index(i, i+1);

            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                ctx.fft_trunc(itrunc, otrunc);
            timeit_stop(timer);
            time = double(timer->wall)/nreps;

            double l = log2(otrunc);
            std::cout << "trunc 2^"
                      << format_fixed(l, 2, 2)
                      << ": "
                      << format_fixed(time, 5)
                      << " ("
                      << format_fixed(time*tmul/(l*otrunc), 1, 2)
                      << ")  | "
                      << std::flush;

            fft_trunc_times.push_back(l);
            fft_trunc_times.push_back(time*tmul/(l*otrunc));


            timeit_start(timer);
            for (ulong i = 0; i < nreps; i++)
                ctx.ifft_trunc(otrunc);
            timeit_stop(timer);
            time = double(timer->wall)/nreps;

            std::cout << format_fixed(time, 5)
                      << " ("
                      << format_fixed(time*tmul/(l*otrunc), 1, 2)
                      << ") |"
                      << std::endl;

            ifft_trunc_times.push_back(l);
            ifft_trunc_times.push_back(time*tmul/(l*otrunc));

            if (nreps == 1)
            {
                double m = reduce_0n_to_pmhn(nmod_pow_ui(2, L, ctx.mod), ctx.p);
                for (ulong i = 0; i < itrunc; i++)
                {
                    if (reduce_to_0n(ctx.get_index(i), ctx.p, ctx.pinv) !=
                        reduce_to_0n(mulmod2(i+1, m, ctx.p, ctx.pinv), ctx.p, ctx.pinv))
                    {
                        std::cout << "oops! error at index " << i << std::endl;
                        std::abort();
                    }
                }
            }
        }

        for (ulong i = 0; i < pow2(L); i++)
            ctx.set_index(i, i+1);

        timeit_start(timer);
        for (ulong i = 0; i < nreps; i++)
            ctx.fft();
        timeit_stop(timer);
        time = double(timer->wall)/nreps;
        double l = L;
        std::cout << " full 2^"
                  << format_fixed(l, 2)
                  << "   : "
                  << format_fixed(time, 5)
                  << " ("
                  << format_fixed(time*tmul/(l*pow2(L)), 1, 2)
                  << ")  | "
                  << std::flush;

        fft_times.push_back(l);
        fft_times.push_back(time*tmul/(l*pow2(L)));

        timeit_start(timer);
        for (ulong i = 0; i < nreps; i++)
            ctx.ifft();
        timeit_stop(timer);
        time = double(timer->wall)/nreps;
        std::cout << format_fixed(time, 5)
                  << " ("
                  << format_fixed(time*tmul/(l*pow2(L)), 1, 2)
                  << ") |"
                  << std::endl;

        ifft_times.push_back(l);
        ifft_times.push_back(time*tmul/(l*pow2(L)));

        delete[] ctx.release_data();
    }

#if 0
    std::cout << "fft_trunc: " << std::endl;
    display_dpoints(fft_trunc_times);
    std::cout << "ifft_trunc: " << std::endl;
    display_dpoints(ifft_trunc_times);
    std::cout << "fft: " << std::endl;
    display_dpoints(fft_times);
    std::cout << "ifft: " << std::endl;
    display_dpoints(ifft_times);
#endif
}


void test_v2_mul(ulong max_len, ulong reps)
{
    std::cout << "--------- testing mpn_mul ---------" << std::endl;

    flint_rand_t state;
    flint_randinit(state);

    mpn_ctx_v2 Q(0x0003f00000000001);

    ulong dgs = n_sizeinbase(reps, 10);

    std::cout << " " << format_fixed(0, dgs) << "/" << format_fixed(reps, dgs) << std::flush;

    for (ulong rep = 0; rep < reps; rep++)
    {
        ulong an = n_randint(state, max_len) + 1;
        ulong bn = n_randint(state, an) + 1;
        ulong zn = an + bn;
        ulong* a = new ulong[an];
        ulong* b = new ulong[bn];
        ulong* z1 = new ulong[zn];
        ulong* z2 = new ulong[zn];

        for (ulong ii = 0; ii < 2*dgs + 2; ii++)
            std::cout << '\b';

        std::cout << " " << format_fixed(1+rep, dgs) << "/" << format_fixed(reps, dgs) << std::flush;

        for (ulong i = 0; i < an; i++)
            a[i] = n_randlimb(state);
        for (ulong i = 0; i < bn; i++)
            b[i] = n_randlimb(state);

        a[an-1] += (a[an-1] == 0);
        b[bn-1] += (b[bn-1] == 0);

        for (ulong i = 0; i < zn; i++)
            z1[i] = n_randlimb(state);

        mpn_mul_v2p1(Q, z1, a, an, b, bn);
        mpn_mul(z2, a, an, b, bn);

        for (ulong i = 0; i < an + bn; i++)
        {
            if (z1[i] != z2[i])
            {
                std::cout << "mpn_mul error" << std::endl;
                std::cout << "an: " << an << ", bn: " << bn << std::endl;
                std::cout << "z1[" << i << "]: " << format_hex(z1[i]) << std::endl;
                std::cout << "z2[" << i << "]: " << format_hex(z2[i]) << std::endl;
                std::abort();
            }
        }
    }

    std::cout << std::endl;

    flint_randclear(state);
}



void profile_v2_mul(ulong max_len, bool use_flint)
{
    std::cout << "--------- profiling mpn_mul ---------" << std::endl;

    mpn_ctx_v2 Q(0x0003f00000000001);
    timeit_t timer;
    flint_rand_t state;
    flint_randinit(state);

    ulong* data = new ulong[2*max_len];

    for (ulong i = 0; i < 2*max_len; i++)
        data[i] = -ulong(1);

    for (ulong zn = 100000; zn <= max_len; zn += 1 + zn*(5 + n_randint(state, 8))/64)
    {
        std::cout << format_fixed(zn, n_sizeinbase(max_len,10)) << ":";

        double ratio_sum = 0;
        double ratio_min = 1000000;
        double ratio_max = 0;

        ulong nreps = 8;
        for (ulong rep = 0; rep < nreps; rep++)
        {
            ulong an = (zn+1)/2 + (rep*zn)/(4*nreps);
            if (an >= zn)
                break;

            ulong bn = zn - an;
            ulong* z = data + max_len;
            ulong* a = data;
            ulong* b = data + an;

            timeit_start(timer);
            mpn_mul_v2p1(Q, z, a, an, b, bn);
            timeit_stop(timer);
            ulong new_time = std::max(slong(1), timer->wall);

            if (z[0]    !=  ulong(1) ||
                z[bn-1] !=  ulong(bn == 1 ? 1 : 0) ||
                z[bn]   != -ulong(an == bn ? 2 : 1) ||
                z[zn-1] != -ulong(bn == 1 ? 2 : 1))
            {
                std::cout << "mpn_mul error" << std::endl;
                std::abort();
            }

            timeit_start(timer);
            if (use_flint)
                flint_mpn_mul_fft_main(z, a, an, b, bn);
            else
                mpn_mul(z, a, an, b, bn);
            timeit_stop(timer);
            ulong old_time = std::max(slong(1), timer->wall);

            double ratio = double(old_time)/double(new_time);
            ratio_sum += ratio;
            ratio_min = std::min(ratio_min, ratio);
            ratio_max = std::max(ratio_max, ratio);

            std::cout << format_fixed(new_time, 5)
                      << "("
                      << format_fixed(ratio, 1, 2)
                      << ")"
                      << std::flush;
        }

        std::cout << " | "
                  << format_fixed(ratio_min, 1, 2)
                  << " "
                  << format_fixed(ratio_sum/nreps, 1, 2)
                  << " "
                  << format_fixed(ratio_max, 1, 2)
                  << std::endl;
    }

    flint_randclear(state);

    delete[] data;
}
