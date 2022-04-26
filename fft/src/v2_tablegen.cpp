void fftv2_ctx::init_prime(ulong pp)
{
    p = pp;
    pinv = 1.0/p;
    nmod_init(&mod, pp);
    primitive_root = n_primitive_root_prime(pp);

    // fill wtab to a depth of 10 (512 entries: 1, exp(pi*i/2), exp(pi*i/4), exp(pi*i*3/4), ...) 
    wtab_depth = 10;
    ulong N = pow2(wtab_depth-1);

    w2s = reinterpret_cast<double*>(my_aligned_alloc(4096, round_up(N*sizeof(double), 4096)));

    ulong w = nmod_pow_ui(primitive_root, (mod.n - 1)>>wtab_depth, mod);
    ulong wi = 1;
    for (ulong j = 0; j < N; j++)
    {
        w2s[n_revbin(j, wtab_depth)/2] = reduce_0n_to_pmhn(wi, p);
        wi = nmod_mul(wi, w, mod);
    }
}

void fftv2_ctx::fit_wtab(ulong k)
{
    if (wtab_depth >= k)
        return;

    if (w2s == nullptr)
        init_prime(mod.n);

    double* oldw2s = w2s;
    w2s = reinterpret_cast<double*>(my_aligned_alloc(4096, pow2(k)*sizeof(double)));
    std::memcpy(w2s, oldw2s, pow2(wtab_depth-1)*sizeof(double));
    free(oldw2s);

    while (wtab_depth < k)
    {
        wtab_depth++;
        ulong N = pow2(wtab_depth-1);

        slong ww = nmod_pow_ui(primitive_root, (mod.n - 1)>>wtab_depth, mod);
        packed<double, VEC_SZ> w = reduce_0n_to_pmhn(double(ww), p, 0.5*p);
        packed<double, VEC_SZ> n = p;
        packed<double, VEC_SZ> ninv = pinv;
        packed<double, VEC_SZ> hn = half(p);
        packed<double, VEC_SZ> x0, x1;
        double* wptr = w2s;
        ulong j = 0;
        do {
            x0.load_aligned(wptr + j);
            x1.load_aligned(wptr + j + VEC_SZ);
            x0 = mulmod2(x0, w, n, ninv);
            x1 = mulmod2(x1, w, n, ninv);
            x0 = reduce_pm1n_to_pmhn(x0, n, hn);
            x1 = reduce_pm1n_to_pmhn(x1, n, hn);
            x0.store_aligned(wptr + N/2 + j);
            x1.store_aligned(wptr + N/2 + j + VEC_SZ);
            j += 2*VEC_SZ;
        } while (j < N/2);
        assert(j == N/2);
    }
}
