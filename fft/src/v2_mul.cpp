



void mpn_ctx_v2::add_prime(ulong p)
{
    ffts.emplace_back(p);
    ulong len = crts.back().coeff_len;
    ulong* t = new ulong[2*(len + 2)];
    ulong* tt = t + (len + 2);

    t[len + 1] = 0;
    t[len] = mpn_mul_1(t, crts.back().prod_primes(), len, p);

    // leave enough room for (product of primes)*(number of primes)
    len += 2;
    mpn_mul_1(tt, t, len, nprimes());

    while (tt[len-1] == 0)
        len--;

    crts.emplace_back(p, len, nprimes());

    // set product of primes
    mpn_copyi(crts.back().prod_primes(), t, len);

    // set cofactors
    for (ulong i = 0; i < nprimes(); i++)
    {
        mpn_divexact_1(crts.back().co_prime(i), t, len, crts[i].prime);
        crts.back().co_prime_red(i) = mpn_mod_1(crts.back().co_prime(i), len, crts[i].prime);
    }

    delete[] t;
}


const packed<double,VEC_SZ>* mpn_ctx_v2::two_pow_table(ulong len, ulong np)
{
    ulong nvs = cdiv(np, VEC_SZ);

    while (nprimes() < nvs*VEC_SZ)
        add_prime();

    while (two_pow_tabs.size() < nvs)
        two_pow_tabs.emplace_back();

    if (two_pow_tabs[nvs-1].size() >= len*nvs)
        return two_pow_tabs[nvs-1].data();

    two_pow_tabs[nvs-1].resize(len*nvs);
    packed<double,VEC_SZ>* d = two_pow_tabs[nvs-1].data();

    packed<double,VEC_SZ>* ps = new packed<double,VEC_SZ>[2*nvs];
    for (ulong l = 0; l < nvs; l++)
    {
#if VEC_SZ == 4
        ps[2*l+0] = packed<double,4>(ffts[4*l+0].p, ffts[4*l+1].p, ffts[4*l+2].p, ffts[4*l+3].p);
        ps[2*l+1] = mul(2.0, packed<double,4>(ffts[4*l+0].pinv, ffts[4*l+1].pinv, ffts[4*l+2].pinv, ffts[4*l+3].pinv));
#elif VEC_SZ == 2
        ps[2*l+0] = packed<double,2>(ffts[2*l+0].p, ffts[2*l+1].p);
        ps[2*l+1] = mul(2.0, packed<double,2>(ffts[2*l+0].pinv, ffts[2*l+1].pinv));
#else
    #error "unsuported VEC_SZ"
#endif
    }

    for (ulong l = 0; l < nvs; l++)
        d[0*nvs+l] = 1;

    for (ulong i = 1; i < len; i++)
    for (ulong l = 0; l < nvs; l++)
    {
        packed<double,VEC_SZ> t = d[(i-1)*nvs+l];
        packed<double,VEC_SZ> p = ps[2*l+0];
        packed<double,VEC_SZ> two_over_p = ps[2*l+1];
        d[i*nvs+l] = fnmadd(round(mul(t, two_over_p)), p, add(t,t));
    }

    delete[] ps;

    return d;
}


template<ulong np, ulong bits>
void mpn_to_ffts(
    fftv2_ctx* Qffts,
    const ulong* a_, ulong an_, ulong atrunc,
    const packed<double,VEC_SZ>* two_pow)
{
    ulong nvs = (np + VEC_SZ - 1)/VEC_SZ;

    FLINT_ASSERT(bits >= FLINT_BITS);

    const uint32_t* a = reinterpret_cast<const uint32_t*>(a_);
    ulong an = 2*an_;

    packed<double,VEC_SZ> X[nvs];
    packed<double,VEC_SZ> P[nvs];
    packed<double,VEC_SZ> PINV[nvs];

    for (ulong l = 0; l < nvs; l++)
    {
#if VEC_SZ == 4
        P[l] = packed<double,4>(Qffts[4*l+0].p, Qffts[4*l+1].p, Qffts[4*l+2].p, Qffts[4*l+3].p);
        PINV[l] = packed<double,4>(Qffts[4*l+0].pinv, Qffts[4*l+1].pinv, Qffts[4*l+2].pinv, Qffts[4*l+3].pinv);
#else
    #error "unsupported VEC_SZ"
#endif
    }

    // if i*bits + 32 < 32*an, then aindex is easy
    ulong end_easy = std::min(atrunc, (32*an - 33)/bits);

    // if i*bits >= 32*an, then aindex is zero
    ulong end_hard = std::min(atrunc, (32*an + bits - 1)/bits);


    ulong i = 0;

#define CODE(ir)\
    {\
        ulong k = ((i+ir)*bits)/32;\
        ulong j = ((  ir)*bits)%32;\
\
        packed<double,4> ak = double(a[k] >> j);\
        for (ulong l = 0; l < nvs; l++)\
            X[l] = ak;\
        k++;\
        j = 32 - j;\
        while (j + 32 <= bits)\
        {\
            ak = double(a[k]);\
            for (ulong l = 0; l < nvs; l++)\
                X[l] = add(X[l], mulmod2(ak, two_pow[j*nvs+l], P[l], PINV[l]));\
            k++;\
            j += 32;\
        }\
\
        if ((bits-j) != 0)\
        {\
            ak = double(a[k] << (32-(bits-j)));\
            for (ulong l = 0; l < nvs; l++)\
                X[l] = add(X[l], mulmod2(ak, two_pow[(bits-32)*nvs+l], P[l], PINV[l]));\
        }\
\
        for (ulong l = 0; l < nvs; l++)\
            X[l] = reduce_to_pm1n(X[l], P[l], PINV[l]);\
\
        for (ulong l = 0; l < np; l++)\
            Qffts[l].set_index(i+ir, X[l/4].data[l%4]);\
    }

    if ((bits % 4) == 0)
    {
        end_easy &= -ulong(8);
        for ( ; i < end_easy; i += 8)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
        }
    }
    else if ((bits % 2) == 0)
    {
        end_easy &= -ulong(16);
        for ( ; i < end_easy; i += 16)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
            CODE(8);CODE(9);CODE(10);CODE(11);
            CODE(12);CODE(13);CODE(14);CODE(15);
        }
    }
    else
    {
        end_easy &= -ulong(32);
        for ( ; i < end_easy; i += 32)
        {
            CODE(0);CODE(1);CODE(2);CODE(3);
            CODE(4);CODE(5);CODE(6);CODE(7);
            CODE(8);CODE(9);CODE(10);CODE(11);
            CODE(12);CODE(13);CODE(14);CODE(15);
            CODE(16);CODE(17);CODE(18);CODE(19);
            CODE(20);CODE(21);CODE(22);CODE(23);
            CODE(24);CODE(25);CODE(26);CODE(27);
            CODE(28);CODE(29);CODE(30);CODE(31);
        }
    }
#undef CODE

    for (; i < end_hard; i++)
    {
        ulong k = (i*bits)/32;
        ulong j = (i*bits)%32;

#define aindex(i) (((i) < an) ? a[i] : uint32_t(0))

        packed<double,VEC_SZ> ak = double(aindex(k) >> j);
        for (ulong l = 0; l < nvs; l++)
            X[l] = ak;
        k++;
        j = 32 - j;
        while (j + 32 <= bits)
        {
            ak = double(aindex(k));
            for (ulong l = 0; l < nvs; l++)
                X[l] = add(X[l], mulmod2(ak, two_pow[j*nvs+l], P[l], PINV[l]));
            k++;
            j += 32;
        }

        if ((bits-j) != 0)
        {
            ak = double(aindex(k) << (32-(bits-j)));
            for (ulong l = 0; l < nvs; l++)
                X[l] = add(X[l], mulmod2(ak, two_pow[(bits-32)*nvs+l], P[l], PINV[l]));
        }

#undef aindex

        for (ulong l = 0; l < nvs; l++)
            X[l] = reduce_to_pm1n(X[l], P[l], PINV[l]);

        for (ulong l = 0; l < np; l++)
            Qffts[l].set_index(i, X[l/VEC_SZ].data[l%VEC_SZ]);
    }

    for (ulong l = 0; l < np; l++)
        for (ulong j = i; j < atrunc; j++)
            Qffts[l].set_index(j, 0.0);
}



/*
    m = cofac_len
    n = coeff_len

    let c[i] = p[0]*...*p[np-1]/p[i] < 2^(64*m)

    for 0 <= x[i] < p[i] the bound n has
        sum_i c[i]*x[i] < 2^(64*n)
    and
        p[0]*...*p[np-1] < 2^(64*n)

    suppose n = m or n = m + 1

    r[0], r[1], ..., r[n-1]     for even products
          t[1], ..., t[n-1]     for odd products

    odd m = 3, n = m
        r[2] r[1] r[0]
        t[2] t[1]

    odd m = 3, n = m + 1
        r[3] r[2] r[1] r[0]
             t[2] t[1]

    even m = 4, n = m
        r[3] r[2] r[1] r[0]
        t[3] t[2] t[1]

    even m = 4, n = m + 1
        r[4] r[3] r[2] r[1] r[0]
        t[4] t[3] t[2] t[1]
*/

template <ulong n, ulong m, bool first>
FORCE_INLINE static inline void
_big_addmul(ulong r[], ulong t[], ulong C[], ulong y)
{
    for (ulong k = 0; k < n; k += 2)
    {
        if (k + 1 < n)
        {
            assert(k < m);
            if (first)
                _mul(r[k+1],r[k+0], C[k+0], y);
            else
                _madd(r[k+1],r[k+0], C[k+0], y);
        }
        else
        {
            assert(k + 1 == n);
            if (k < m)
            {
                if (first)
                    r[k+0] = C[k+0]*y;
                else
                    r[k+0] += C[k+0]*y;
            }
            else
            {
                if (first)
                    r[k+0] = 0;
            }
        }

        if (k + 2 < n)
        {
            assert(k + 1 < m);
            if (first)
                _mul(t[k+2],t[k+1], C[k+1], y);
            else
                _madd(t[k+2],t[k+1], C[k+1], y);
        }
        else if (k + 1 < n)
        {
            if (k + 1 < m)
            {
                if (first)
                    t[k+1] = C[k+1]*y;
                else
                    t[k+1] += C[k+1]*y;
            }
            else
            {
                if (first)
                    t[k+1] = 0;
            }
        }
    }
}


template <ulong n>
FORCE_INLINE static inline void
_reduce_big_sum(ulong r[], ulong t[], const ulong* limit)
{
    _multi_add<n-1>(r+1, t+1);

check:

    for (ulong k = n; k > 1; k--)
    {
        if (LIKELY(r[k-1] > limit[k-1]))
            goto sub;
        if (r[k-1] < limit[k-1])
            return;
    }

    if (r[0] < limit[0])
        return;
sub:
    _multi_sub<n>(r, limit);
    goto check;
}

template <ulong n, bool easy>
FORCE_INLINE static inline void
_add_to_answer(ulong z[], ulong r[], ulong zn, ulong toff, ulong tshift)
{
    assert(zn > toff);

    if (tshift == 0)
    {
        if (easy || zn - toff >= n)
        {
            _multi_add<n>(z + toff, r);
            return;
        }
    }
    else
    {
        r[n] = r[n-1] >> (64-tshift);
        for (ulong k = n; k >= 2; k--)
            r[k-1] = (r[k-1] << (tshift)) | (r[k-2] >> (64-tshift));
        r[0] =  r[0] << (tshift);

        if (easy || zn - toff > n)
        {
            _multi_add<n + 1>(z + toff, r);
            return;
        }
    }

    zn -= toff;
    unsigned char cf = 0;
    ulong i = 0;
    do {
        cf = _addcarry_ulong(cf, r[i], z[toff + i], z[toff + i]);
    } while (++i < zn);
}


template<ulong np, ulong n, ulong m>
void mpn_from_ffts(
    ulong* z, ulong zn, ulong zlen,
    fftv2_ctx* Qffts,
    crt_data* Qcrts,
    ulong bits)
{
    ulong r[n + 1];
    ulong t[n + 1];

    assert(n == Qcrts[np-1].coeff_len);

    if (n == m + 1)
    {
        for (ulong l = 0; l < np; l++) {
            assert(Qcrts[np - 1].co_prime(l)[m] == 0);
        }
    }
    else
    {
        assert(n == m);
    }

    mpn_zero(z, zn);

    ulong i = 0;

    // easy if zn-n > floor(i*bits/64)
    ulong end_easy = (zn >= n+1 ? zn - (n+1) : ulong(0))*FLINT_BITS/bits;

    ulong Xs[BLK_SZ*np];

    end_easy &= -BLK_SZ;

    for (; i < end_easy; i += BLK_SZ)
    {
        ulong I = i/BLK_SZ;

        for (ulong l = 0; l < np; l++)
        {
            packed<double,VEC_SZ> P = Qffts[l].p;
            packed<double,VEC_SZ> PINV = Qffts[l].pinv;
            double* x = Qffts[l].from_index(I);
            for (ulong j = 0; j < BLK_SZ; j += 4*VEC_SZ)
            {
                packed<double,VEC_SZ> x0, x1, x2, x3;
                packed<ulong,VEC_SZ> y0, y1, y2, y3;
                x0.load(x + j + 0*VEC_SZ);
                x1.load(x + j + 1*VEC_SZ);
                x2.load(x + j + 2*VEC_SZ);
                x3.load(x + j + 3*VEC_SZ);
                x0 = reduce_to_0n(x0, P, PINV);
                x1 = reduce_to_0n(x1, P, PINV);
                x2 = reduce_to_0n(x2, P, PINV);
                x3 = reduce_to_0n(x3, P, PINV);
                y0 = convert_limited<packed<ulong,VEC_SZ>>(x0);
                y1 = convert_limited<packed<ulong,VEC_SZ>>(x1);
                y2 = convert_limited<packed<ulong,VEC_SZ>>(x2);
                y3 = convert_limited<packed<ulong,VEC_SZ>>(x3);
                y0.store(Xs + l*BLK_SZ + j + 0*VEC_SZ);
                y1.store(Xs + l*BLK_SZ + j + 1*VEC_SZ);
                y2.store(Xs + l*BLK_SZ + j + 2*VEC_SZ);
                y3.store(Xs + l*BLK_SZ + j + 3*VEC_SZ);
            }
        }

        for (ulong j = 0; j < BLK_SZ; j += 1)
        {
            ulong l = 0;
            _big_addmul<n,m,true>(r, t, Qcrts[np - 1].co_prime(l), Xs[l*BLK_SZ + j]);
            for (l++; l < np; l++)
                _big_addmul<n,m,false>(r, t, Qcrts[np - 1].co_prime(l), Xs[l*BLK_SZ + j]);

            _reduce_big_sum<n>(r, t, Qcrts[np - 1].prod_primes());

            ulong toff = ((i+j)*bits)/FLINT_BITS;
            ulong tshift = ((i+j)*bits)%FLINT_BITS;

            assert(zn > n + toff);

            _add_to_answer<n, true>(z, r, zn, toff, tshift);
        }
    }

    for (; i < zlen; i++)
    {
        for (ulong l = 0; l < np; l++)
        {
            ulong x = reduce_to_0n(Qffts[l].get_index(i), Qffts[l].p, Qffts[l].pinv);
            if (l == 0)
                _big_addmul<n,m,true>(r, t, Qcrts[np - 1].co_prime(l), x);
            else
                _big_addmul<n,m,false>(r, t, Qcrts[np - 1].co_prime(l), x);
        }

        _reduce_big_sum<n>(r, t, Qcrts[np - 1].prod_primes());

        ulong toff = (i*bits)/FLINT_BITS;
        ulong tshift = (i*bits)%FLINT_BITS;

        if (toff >= zn)
            break;

        _add_to_answer<n,false>(z, r, zn, toff, tshift);
    }
}


template <ulong np, ulong bits, ulong n, ulong m>
void mpn_ctx_v2::push_profile()
{
    while (nprimes() < np)
        add_prime();

    ulong bound = crts[np-1].find_bound(bits);

    // at least exclude bound = 0
    assert(bound > 10);

    profiles.emplace_back(np, bits, bound, &mpn_to_ffts<np,bits>, &mpn_from_ffts<np,n,m>);
}


mpn_ctx_v2::mpn_ctx_v2(ulong p) : double_buffer(nullptr), double_buffer_alloc(0)
{
    ffts.emplace_back(p);
    crts.emplace_back(p, 1, 1);

    crts[0].co_prime_red(0) = 1;
    crts[0].co_prime(0)[0] = 1;
    crts[0].prod_primes()[0] = p;

    // first must always works
    push_profile<4,64, 4,3>();

    // may not always work
    push_profile<3,64, 3,2>();
    push_profile<3,66, 3,2>();
    push_profile<3,68, 3,2>();
    push_profile<3,70, 3,2>();
    push_profile<3,72, 3,2>();
    push_profile<4,88, 4,3>();
    push_profile<4,90, 4,3>();
    push_profile<4,92, 4,3>();
    push_profile<5,112, 4,4>();
    push_profile<5,114, 4,4>();
    push_profile<5,116, 4,4>();
    push_profile<6,136, 5,4>();
    push_profile<6,138, 5,4>();
    push_profile<6,140, 5,4>();
    push_profile<6,142, 5,4>();
    push_profile<7,162, 6,5>();
    push_profile<7,164, 6,5>();
    push_profile<7,166, 6,5>();
    push_profile<8,188, 7,6>();
}

const profile_entry* mpn_ctx_v2::best_profile(ulong an, ulong bn)
{
    ulong i = 0;
    ulong best_i = 0;
    double best_score = 100000000.0*(an + bn);

find_next:

    do {
        i++;
        if (i >= profiles.size())
            return profiles.data() + best_i;
    } while (bn > profiles[i].bn_bound);

maximize_bits:

    assert(i < profiles.size() &&
           bn <= profiles[i].bn_bound);

    while (i+1 < profiles.size() &&
           bn <= profiles[i+1].bn_bound &&
           profiles[i+1].np == profiles[i].np)
    {
        i++;
    }

    ulong np = profiles[i].np;
    ulong bits = profiles[i].bits;
    ulong alen = cdiv(64*an, bits);
    ulong blen = cdiv(64*bn, bits);
    ulong zlen = alen + blen - 1;
    ulong atrunc = round_up(alen, BLK_SZ);
    ulong btrunc = round_up(blen, BLK_SZ);
    ulong ztrunc = round_up(zlen, BLK_SZ);
    ulong depth = std::max(ulong(LG_BLK_SZ), clog2(ztrunc));

    double ratio = double(ztrunc)/double(pow2(depth));
    double score = (1-0.25*ratio)*(1.0/1000000);
    score *= np*depth;
    score *= ztrunc;
    if (score < best_score)
    {
        best_i = i;
        best_score = score;
    }

    goto find_next;
}


void mpn_ctx_v2::my_mpn_mul(
    ulong* z,
    const ulong* a, ulong an,
    const ulong* b, ulong bn)
{
    const profile_entry* P = best_profile(an, bn);
    ulong np = P->np;
    ulong bits = P->bits;
    ulong zn = an + bn;
    ulong alen = cdiv(64*an, bits);
    ulong blen = cdiv(64*bn, bits);
    ulong zlen = alen + blen - 1;
    ulong atrunc = round_up(alen, BLK_SZ);
    ulong btrunc = round_up(blen, BLK_SZ);
    ulong ztrunc = round_up(zlen, BLK_SZ);
    ulong depth = std::max(ulong(LG_BLK_SZ), clog2(ztrunc));

    assert(an > 0);
    assert(bn > 0);
    assert(mpn_cmp_ui_2exp(crts[np-1].prod_primes(), crts[np-1].coeff_len, blen, 2*bits) >= 0);

#define TIME_THIS 0

#if TIME_THIS
timeit_t timer, timer_overall;
std::cout << "\n----------------------------" << std::endl;
#endif

    for (ulong l = 0; l < np; l++)
        ffts[l].set_depth(depth);

#if TIME_THIS
timeit_start(timer_overall);
#endif

    ulong fft_data_size = ffts[0].data_size();
    double* abuf = fit_double_buffer(2*np*fft_data_size);
    double* bbuf = abuf + np*fft_data_size;
    const packed<double,VEC_SZ>* two_pow = two_pow_table(bits+5, np);

#if TIME_THIS
timeit_start(timer);
#endif

	for (ulong l = 0; l < np; l++)
        ffts[l].set_data(abuf + l*fft_data_size);
    P->to_ffts(ffts.data(), a, an, atrunc, two_pow);

	for (ulong l = 0; l < np; l++)
        ffts[l].set_data(bbuf + l*fft_data_size);
    P->to_ffts(ffts.data(), b, bn, btrunc, two_pow);

#if TIME_THIS
timeit_stop(timer);
if (timer->wall > 5)
std::cout << "mod: " << timer->wall << std::endl;
#endif

#if TIME_THIS
timeit_start(timer);
#endif
	for (ulong l = 0; l < np; l++)
    {
        ffts[l].fft_trunc(btrunc, ztrunc);
        ffts[l].set_data(abuf + l*fft_data_size);
        ffts[l].fft_trunc(atrunc, ztrunc);
        ulong t1, thi, tlo;
        ulong cop = crts[np - 1].co_prime_red(l);
        thi = cop >> (FLINT_BITS - depth);
        tlo = cop << (depth);
        NMOD_RED2(t1, thi, tlo, ffts[l].mod);
        t1 = nmod_inv(t1, ffts[l].mod);
        ffts[l].point_mul(bbuf + l*fft_data_size, t1);
        ffts[l].ifft_trunc(ztrunc);
    }
#if TIME_THIS
timeit_stop(timer);
if (timer->wall > 5)
std::cout << "fft: " << timer->wall << std::endl;
#endif

#if TIME_THIS
timeit_start(timer);
#endif

    P->from_ffts(z, zn, zlen, ffts.data(), crts.data(), bits);

#if TIME_THIS
timeit_stop(timer);
if (timer->wall > 5)
std::cout << "crt: " << timer->wall << std::endl;
timeit_stop(timer_overall);
if (timer_overall->wall > 5)
std::cout << "   : " << timer_overall->wall << std::endl;
#endif

#undef TIME_THIS
}

