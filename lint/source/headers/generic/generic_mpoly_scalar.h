#pragma once

template <ulimb ML = 0, typename RingT, typename MPolyT>
void _generic_mpoly_add(
    mpoly_ctx& ctxm, RingT& ctxc,
    MPolyT& A,
    typename MPolyT::input_t B,
    typename MPolyT::input_t C)
{
    ulimb M = ctxm.stride(A.bits());
    ulimb* Aexps = A.m.fit_alloc_destroy(M*(B.len + C.len));
    ulimb N = ctxc.stride();
    typename RingT::coeff_t* Acoeffs = A.c.fit_alloc_destroy(N*(B.len + C.len));

    ulimb Ai = 0, Bi = 0, Ci = 0;
    while (Bi < B.len && Ci < C.len)
    {
        int cmp = mpoly_monomial_cmp<ML>(B.exps + M*Bi, C.exps + M*Ci, M);
        if (cmp > 0)
        {
            mpoly_monomial_set<ML>(Aexps + M*Ai, B.exps + M*Bi, M);
            set(ctxc, Acoeffs + N*Ai, B.coeffs + N*Bi);
            Ai++;
            Bi++;
        }
        else if (cmp < 0)
        {
            mpoly_monomial_set<ML>(Aexps + M*Ai, C.exps + M*Ci, M);
            set(ctxc, Acoeffs + N*Ai, C.coeffs + N*Ci);
            Ai++;
            Ci++;
        }
        else
        {
            mpoly_monomial_set<ML>(Aexps + M*Ai, B.exps + M*Bi, M);
            add(ctxc, Acoeffs + N*Ai, B.coeffs + N*Bi, C.coeffs + N*Ci);
            Ai += !is_zero(ctxc, Acoeffs + N*Ai);
            Bi++;
            Ci++;
        }
    }

    while (Bi < B.len)
    {
        mpoly_monomial_set<ML>(Aexps + M*Ai, B.exps + M*Bi, M);
        set(ctxc, Acoeffs + N*Ai, B.coeffs + N*Bi);
        Ai++;
        Bi++;
    }

    while (Ci < C.len)
    {
        mpoly_monomial_set<ML>(Aexps + M*Ai, C.exps + M*Ci, M);
        set(ctxc, Acoeffs + N*Ai, C.coeffs + N*Ci);
        Ai++;
        Ci++;
    }

    A.set_length(Ai);
}


template <ulimb ML = 0, typename RingT, typename MPolyT>
void _generic_mpoly_sub(
    mpoly_ctx& ctxm, RingT& ctxc,
    MPolyT& A,
    typename MPolyT::input_t B,
    typename MPolyT::input_t C)
{
    ulimb M = ctxm.stride(A.bits());
    ulimb* Aexps = A.m.fit_alloc_destroy(M*(B.len + C.len));
    ulimb N = ctxc.stride();
    typename RingT::coeff_t* Acoeffs = A.c.fit_alloc_destroy(N*(B.len + C.len));

    ulimb Ai = 0, Bi = 0, Ci = 0;
    while (Bi < B.len && Ci < C.len)
    {
        int cmp = mpoly_monomial_cmp<ML>(B.exps + M*Bi, C.exps + M*Ci, M);
        if (cmp > 0)
        {
            mpoly_monomial_set<ML>(Aexps + M*Ai, B.exps + M*Bi, M);
            set(ctxc, Acoeffs + N*Ai, B.coeffs + N*Bi);
            Ai++;
            Bi++;
        }
        else if (cmp < 0)
        {
            mpoly_monomial_set<ML>(Aexps + M*Ai, C.exps + M*Ci, M);
            neg(ctxc, Acoeffs + N*Ai, C.coeffs + N*Ci);
            Ai++;
            Ci++;
        }
        else
        {
            mpoly_monomial_set<ML>(Aexps + M*Ai, B.exps + M*Bi, M);
            sub(ctxc, Acoeffs + N*Ai, B.coeffs + N*Bi, C.coeffs + N*Ci);
            Ai += !is_zero(ctxc, Acoeffs + N*Ai);
            Bi++;
            Ci++;
        }
    }

    while (Bi < B.len)
    {
        mpoly_monomial_set<ML>(Aexps + M*Ai, B.exps + M*Bi, M);
        set(ctxc, Acoeffs + N*Ai, B.coeffs + N*Bi);
        Ai++;
        Bi++;
    }

    while (Ci < C.len)
    {
        mpoly_monomial_set<ML>(Aexps + M*Ai, C.exps + M*Ci, M);
        neg(ctxc, Acoeffs + N*Ai, C.coeffs + N*Ci);
        Ai++;
        Ci++;
    }

    A.set_length(Ai);
}

