
template <int vsize>
void pd_fft_ctx<vsize>::set_mod(const ulong* p)
{
    packed<ulong,vsize> t;
    t.load(p);
    n = convert_limited<packed<double,vsize>>(t);
    ninv = div(1.0, n);
    for (int i = 0; i < vsize; i++)
    {
        nmod_init(&mod[i], p[i]);
        primitive_root[i] = n_primitive_root_prime(p[i]);
    }
    wrevtab.clear();
}


template <int vsize>
void pd_fft_ctx<vsize>::fit_wrevtab(ulong depth)
{
    depth = std::max(depth, ulong(6));
    ulong l = pow2(depth-1);
    ulong k = wrevtab.size();
    if (LIKELY(k >= l))
        return;

    wrevtab.resize(l);

    if (k == 0)
    {
        wrevtab[0] = 1.0;
        k = 1;
    }

    while (k < l)
    {
        packed<double,vsize> w;
        packed<double,vsize> N = n;
        packed<double,vsize> Ninv= ninv;
        packed<double,vsize> hN = half(N);
        packed<ulong,vsize> ww;
        ulong www[vsize];

        for (int i = 0; i < vsize; i++)
            www[i] = nmod_pow_ui(primitive_root[i], (mod[i].n - 1)/(4*k), mod[i]);

        ww.load(www);
        w = reduce_0n_to_pmhn(convert_limited<packed<double,vsize>>(ww), N, hN);

        ulong j = 0;
        if (k == 1 || k == 2)
        {
            do {
                wrevtab[j+k] = reduce_pm1n_to_pmhn(mulmod2(wrevtab[j], w, N, Ninv), N, hN);
            } while (j += 1, j < k);
        }
        else
        {
            do {
                packed<double,vsize> x0 = wrevtab[j+0];
                packed<double,vsize> x1 = wrevtab[j+1];
                packed<double,vsize> x2 = wrevtab[j+2];
                packed<double,vsize> x3 = wrevtab[j+3];
                x0 = mulmod2(x0, w, N, Ninv);
                x1 = mulmod2(x1, w, N, Ninv);
                x2 = mulmod2(x2, w, N, Ninv);
                x3 = mulmod2(x3, w, N, Ninv);
                wrevtab[j+0+k] = reduce_pm1n_to_pmhn(x0, N, hN);
                wrevtab[j+1+k] = reduce_pm1n_to_pmhn(x1, N, hN);
                wrevtab[j+2+k] = reduce_pm1n_to_pmhn(x2, N, hN);
                wrevtab[j+3+k] = reduce_pm1n_to_pmhn(x3, N, hN);
            } while (j += 4, j < k);
        }

        k = 2*k;
    }
}
