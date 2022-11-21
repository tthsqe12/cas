#pragma once

//////////////////////////////// heaps /////////////////////////////////////////

// a chain of nodes with the same exponent
struct mpoly_heap_chain
{
   ulimb i, j;
   mpoly_heap_chain* next;
};

/*
    The heap can either store the exponents (of size ML) directly in the
    mpoly_heap_entry_imm. Since these are shuffled while sorting, at some point
    it is more efficient to store pointers as in mpoly_heap_entry_ptr and just
    shuffle the pointers.
*/
template <ulimb ML>
struct mpoly_heap_entry_imm {
    static_assert(ML > 0);
    ulimb exp[ML];
    mpoly_heap_chain* next;
};

struct mpoly_heap_entry_ptr {
   ulimb* exp;
   mpoly_heap_chain* next;
};


#define HEAP_LEFT(i) (2*(i))
#define HEAP_RIGHT(i) (2*(i) + 1)
#define HEAP_PARENT(i) ((i)/2)

// using mpoly_heap_entry_imm for small ML
template <bool MP, ulimb ML>
struct mpoly_heap {
    ulimb next_loc;
    ulimb heap_len; /* heap zero index unused */
    mpoly_heap_entry_imm<ML>* heap;
    mpoly_heap_chain* chain;
    ulimb* store, * store_base;
    ulimb* ind;

    bool is_empty() {return heap_len <= 1;}

    ulimb* top_exp() {
        FLINT_ASSERT(!is_empty());
        return heap[1].exp;
    }

    mpoly_heap(ulimb Blen, ulimb M, tmp_allocator& push)
    {
        FLINT_ASSERT(Blen > 0);
        next_loc = Blen + 4;   /* something bigger than heap can ever be */
        heap = push.recursive_alloc<mpoly_heap_entry_imm<ML>>(Blen + 1);
        chain = push.recursive_alloc<mpoly_heap_chain>(Blen);
        store = store_base = push.recursive_alloc<ulimb>(2*Blen);

        ind = push.recursive_alloc<ulimb>(Blen + 2);
        ind += 1;
        ind[-1] = ind[Blen] = UWORD_MAX;

        for (ulimb i = 0; i < Blen; i++)
            ind[i] = 1;

        /* start with no heap nodes */
        heap_len = 1;
    }

    void start_mul(const ulimb* Bexps, const ulimb* Cexps, ulimb M)
    {
        mpoly_heap_chain* x = chain + 0;
        x->i = 0;
        x->j = 0;
        x->next = NULL;
        heap[1].next = x;
        mpoly_monomial_add<MP, ML>(heap[1].exp, Bexps + M*0, Cexps + M*0, M);
        ind[0] = 2*1 + 0;
        heap_len = 2;
    }

    void start_div(const ulimb* Aexps, ulimb M)
    {
        mpoly_heap_chain* x = chain + 0;
        x->i = -ulimb(1);
        x->j = 0;
        x->next = NULL;
        heap[1].next = x;
        mpoly_monomial_set<ML>(heap[1].exp, Aexps, M);
        heap_len = 2;

        ind[0] = UWORD_MAX;
    }

    mpoly_heap_chain* pop(ulimb M)
    {
        mpoly_heap_chain* x = heap[1].next;
        ulimb i = 1, j = 2, s = --heap_len;
        while (j < s) {
            j += !mpoly_monomial_gt<ML>(heap[j].exp, heap[j + 1].exp, M);
            heap[i] = heap[j];
            i = j;
            j = HEAP_LEFT(j);
        }
        // insert last element into heap[i]
        ulimb exp[ML];
        mpoly_monomial_set<ML>(exp, heap[s].exp, M);
        j = HEAP_PARENT(i);
        while (i > 1 && mpoly_monomial_gt<ML>(exp, heap[j].exp, M)) {
            heap[i] = heap[j];
            i = j;
            j = HEAP_PARENT(j);
        }
        heap[i] = heap[s];
        return x; 
    }

    void insert(mpoly_heap_chain* x, ulimb* exp, ulimb M)
    {
        ulimb i = heap_len, j, n = heap_len;
        if (i != 1 && mpoly_monomial_equal<ML>(exp, heap[1].exp, M)) {
            x->next = heap[1].next;
            heap[1].next = x;
            return;
        }
        if (next_loc < heap_len && mpoly_monomial_equal<ML>(exp, heap[next_loc].exp, M)) {
            x->next = heap[next_loc].next;
            heap[next_loc].next = x;
            return;
        }
        while ((j = HEAP_PARENT(i)) >= 1) {
            int cmp = mpoly_monomial_cmp<ML>(exp, heap[j].exp, M);
            if (cmp == 0) {
                x->next = heap[j].next;
                heap[j].next = x;
                next_loc = j;
                return;
            } else if (cmp > 0) {
                i = j;
            } else {
                break;
            }
        }
        heap_len++;
        while (n > i) {
            heap[n] = heap[HEAP_PARENT(n)];
            n = HEAP_PARENT(n);
        }
        mpoly_monomial_set<ML>(heap[i].exp, exp, M);
        heap[i].next = x;
    }

    void insert_mul(const ulimb* Bexps, ulimb ii, const ulimb* Cexps, ulimb jj, ulimb M)
    {
        mpoly_heap_chain* x = chain + ii;
        x->next = nullptr;
        x->i = ii;
        x->j = jj;
        ind[x->i] = 2*(x->j + 1) + 0;
        ulimb exp[ML];
        mpoly_monomial_add<MP, ML>(exp, Bexps + M*x->i, Cexps + M*x->j, M);
        insert(x, exp, M);
    }

    void insert_div(const ulimb* Aexps, ulimb ii, ulimb jj, ulimb M)
    {
        mpoly_heap_chain* x = chain + 0;
        x->next = nullptr;
        x->i = ii;
        x->j = jj;
        ind[x->i] = 2*(x->j + 1) + 0;
        ulimb exp[ML];
        mpoly_monomial_set<ML>(exp, Aexps + M*x->j, M);
        insert(x, exp, M);
    }

    void process_store_mul(const ulimb* Bexps, const ulimb* Cexps, ulimb Clen, ulimb M);
    ulimb process_store_div(const ulimb* Aexps, ulimb Alen, const ulimb* Bexps, const ulimb* Qexps, ulimb Qi, ulimb M);
};

// using mpoly_heap_entry_ptr for general M
#define ML 0
template <bool MP>
struct mpoly_heap<MP, ML> {
    ulimb next_loc;
    ulimb heap_len; /* heap zero index unused */
    mpoly_heap_entry_ptr* heap;
    mpoly_heap_chain* chain;
    ulimb* store, * store_base;
    ulimb* exps;
    ulimb** exp_list;
    ulimb exp_next;
    ulimb* ind;

    bool is_empty() {return heap_len <= 1;}

    ulimb* top_exp() {
        FLINT_ASSERT(!is_empty());
        return heap[1].exp;
    }

    mpoly_heap(ulimb Blen, ulimb M, tmp_allocator& push)
    {
        FLINT_ASSERT(Blen > 0);
        next_loc = Blen + 4;   // something bigger than heap can ever be
        heap = push.recursive_alloc<mpoly_heap_entry_ptr>(Blen + 1);
        chain = push.recursive_alloc<mpoly_heap_chain>(Blen);
        store = store_base = push.recursive_alloc<ulimb>(2*Blen);
        exps = push.recursive_alloc<ulimb>(Blen*M);
        exp_list = push.recursive_alloc<ulimb*>(Blen);

        ind = push.recursive_alloc<ulimb>(Blen + 2);
        ind += 1;
        ind[-1] = ind[Blen] = UWORD_MAX;

        for (ulimb i = 0; i < Blen; i++)
        {
            exp_list[i] = exps + M*i;
            ind[i] = 1;
        }

        // start with no heap nodes and no exponent vectors in use
        exp_next = 0;
        heap_len = 1;
    }

    void start_mul(const ulimb* Bexps, const ulimb* Cexps, ulimb M)
    {
        mpoly_heap_chain* x = chain + 0;
        x->i = 0;
        x->j = 0;
        x->next = NULL;
        heap[1].next = x;
        heap[1].exp = exp_list[exp_next++];
        mpoly_monomial_add<MP, ML>(heap[1].exp, Bexps + M*0, Cexps + M*0, M);
        ind[0] = 2*1 + 0;
        heap_len = 2;
    }

    void start_div(const ulimb* Aexps, ulimb M)
    {
        mpoly_heap_chain* x = chain + 0;
        x->i = -ulimb(1);
        x->j = 0;
        x->next = NULL;
        heap[1].next = x;
        heap[1].exp = exp_list[exp_next++];
        mpoly_monomial_set<ML>(heap[1].exp, Aexps, M);
        heap_len = 2;

        ind[0] = UWORD_MAX;
    }

    mpoly_heap_chain* pop(ulimb M)
    {
        exp_list[--exp_next] = heap[1].exp;
        mpoly_heap_chain* x = heap[1].next;
        ulimb i = 1, j = 2, s = --heap_len;
        while (j < s) {
            j += !mpoly_monomial_gt<ML>(heap[j].exp, heap[j + 1].exp, M);
            heap[i] = heap[j];
            i = j;
            j = HEAP_LEFT(j);
        }
        // insert last element into heap[i]
        ulimb* exp = heap[s].exp;
        j = HEAP_PARENT(i);
        while (i > 1 && mpoly_monomial_gt<ML>(exp, heap[j].exp, M)) {
            heap[i] = heap[j];
            i = j;
            j = HEAP_PARENT(j);
        }
        heap[i] = heap[s];
        return x; 
    }

    void insert(mpoly_heap_chain* x, ulimb* exp, ulimb M)
    {
        ulimb i = heap_len, j, n = heap_len;
        if (i != 1 && mpoly_monomial_equal<ML>(exp, heap[1].exp, M)) {
            x->next = heap[1].next;
            heap[1].next = x;
            return;
        }
        if (next_loc < heap_len && mpoly_monomial_equal<ML>(exp, heap[next_loc].exp, M)) {
            x->next = heap[next_loc].next;
            heap[next_loc].next = x;
            return;
        }
        while ((j = HEAP_PARENT(i)) >= 1) {
            if (!mpoly_monomial_gt<ML>(exp, heap[j].exp, M))
                break;
            i = j;
        }
        if (j >= 1 && mpoly_monomial_equal<ML>(exp, heap[j].exp, M)) {
            x->next = heap[j].next;
            heap[j].next = x;
            next_loc = j;
            return;
        }
        heap_len++;
        while (n > i) {
            heap[n] = heap[HEAP_PARENT(n)];
            n = HEAP_PARENT(n);
        }
        heap[i].exp = exp;
        heap[i].next = x;
        exp_next++;
    }

    void insert_mul(const ulimb* Bexps, ulimb ii, const ulimb* Cexps, ulimb jj, ulimb M)
    {
        mpoly_heap_chain* x = chain + ii;
        x->next = nullptr;
        x->i = ii;
        x->j = jj;
        ind[x->i] = 2*(x->j + 1) + 0;
        ulimb* exp = exp_list[exp_next];
        mpoly_monomial_add<MP, ML>(exp, Bexps + M*x->i, Cexps + M*x->j, M);
        insert(x, exp, M);
    }

    void insert_div(const ulimb* Aexps, ulimb ii, ulimb jj, ulimb M)
    {
        mpoly_heap_chain* x = chain + 0;
        x->next = nullptr;
        x->i = ii;
        x->j = jj;
        ind[x->i] = 2*(x->j + 1) + 0;
        ulimb* exp = exp_list[exp_next];
        mpoly_monomial_set<ML>(exp, Aexps + M*x->j, M);
        insert(x, exp, M);
    }

    void process_store_mul(const ulimb* Bexps, const ulimb* Cexps, ulimb Clen, ulimb M);
    ulimb process_store_div(const ulimb* Aexps, ulimb Alen, const ulimb* Bexps, const ulimb* Qexps, ulimb Qi, ulimb M);
};
#undef ML

#undef HEAP_LEFT
#undef HEAP_RIGHT
#undef HEAP_PARENT


template <bool MP, ulimb ML>
void mpoly_heap<MP, ML>::process_store_mul(const ulimb* Bexps, const ulimb* Cexps, ulimb Clen, ulimb M)
{
    while (store > store_base) {
        ulimb j = *--store;
        ulimb i = *--store;
        // we have ind[-1] == UWORD_MAX  and  ind[Blen] == UWORD_MAX;
        // should we go right?
        if (ind[i + 1] == 2*j + 1)
            insert_mul(Bexps, i + 1, Cexps, j, M);
        // should we go up? */
        if (j + 1 < Clen && is_odd(ind[i]) && ind[i - 1] >= 2*j + 5)
            insert_mul(Bexps, i, Cexps, j + 1, M);
    }
}

// TODO why must this be coppied with ML=0?
template <bool MP>
void mpoly_heap<MP, 0>::process_store_mul(const ulimb* Bexps, const ulimb* Cexps, ulimb Clen, ulimb M)
{
    while (store > store_base) {
        ulimb j = *--store;
        ulimb i = *--store;
        if (ind[i + 1] == 2*j + 1)
            insert_mul(Bexps, i + 1, Cexps, j, M);
        if (j + 1 < Clen && is_odd(ind[i]) && ind[i - 1] >= 2*j + 5)
            insert_mul(Bexps, i, Cexps, j + 1, M);
    }
}

template <bool MP, ulimb ML>
ulimb mpoly_heap<MP, ML>::process_store_div(const ulimb* Aexps, ulimb Alen, const ulimb* Bexps, const ulimb* Qexps, ulimb Qi, ulimb M)
{
    ulimb s = 0;
    while (store > store_base)
    {
        ulimb j = *--store;
        ulimb i = *--store;
        // FLINT_ASSERT(ind[0] == UWORD_MAX && ind[Blen] == UWORD_MAX);
        if (i + 1 == 0) {
            // take next dividend term
            if (j + 1 < Alen)
                insert_div(Aexps, i, j + 1, M);
        } else {
            // FLINT_ASSERT(0 < i && i < Blen);
            // should we go right?
            if (ind[i + 1] == 2*j + 1)
                insert_mul(Bexps, i + 1, Qexps, j, M);
            // should we go up?
            if (j + 1 == Qi)
                s++;
            else if (is_odd(ind[i]) && ind[i - 1] >= 2*j + 5)
                insert_mul(Bexps, i, Qexps, j + 1, M);
        }
    }
    return s;
}

// TODO ditto
template <bool MP>
ulimb mpoly_heap<MP, 0>::process_store_div(const ulimb* Aexps, ulimb Alen, const ulimb* Bexps, const ulimb* Qexps, ulimb Qi, ulimb M)
{
    ulimb s = 0;
    while (store > store_base)
    {
        ulimb j = *--store;
        ulimb i = *--store;
        if (i + 1 == 0) {
            if (j + 1 < Alen)
                insert_div(Aexps, i, j + 1, M);
        } else {
            if (ind[i + 1] == 2*j + 1)
                insert_mul(Bexps, i + 1, Qexps, j, M);
            if (j + 1 == Qi)
                s++;
            else if (is_odd(ind[i]) && ind[i - 1] >= 2*j + 5)
                insert_mul(Bexps, i, Qexps, j + 1, M);
        }
    }
    return s;
}


template <bool MP = true, ulimb ML = 0, typename RingT, typename MPolyT>
void _generic_mpoly_mul_heap(
    mpoly_ctx& ctxm, RingT& ctxc,
    MPolyT& A,
    typename MPolyT::input_t B,
    typename MPolyT::input_t C,
    tmp_allocator& push)
{
    FLINT_ASSERT(B.len > 0);
    FLINT_ASSERT(C.len > 0);

    ulimb M = ctxm.stride(A.bits());
    ulimb N = ctxc.stride();
    mpoly_heap<MP, ML> H(B.len, M, push);
    typename RingT::dotter_t dotter(ctxc);

    H.start_mul(B.exps, C.exps, M);

    ulimb Ai = 0;
    while (!H.is_empty())
    {
        ulimb* Aexps = A.m.fit_alloc(M*(Ai + 1));
        mpoly_monomial_set<ML>(Aexps + M*Ai, H.top_exp(), M);

        dotter.zero(ctxc);
        do {
            auto x = H.pop(M);
            do {
                FLINT_ASSERT(x->i < B.len); FLINT_ASSERT(x->j < C.len);
                H.ind[x->i] |= 1;
                *H.store++ = x->i;
                *H.store++ = x->j;
                dotter.madd(ctxc, B.coeffs + N*x->i, C.coeffs + N*x->j);
            } while ((x = x->next) != nullptr);
        } while (!H.is_empty() && mpoly_monomial_equal<ML>(H.top_exp(), Aexps + M*Ai, M));
        typename RingT::coeff_t* Acoeffs = A.c.fit_alloc(N*(Ai + 1));
        dotter.reduce(ctxc, Acoeffs + N*Ai);
        Ai += !is_zero(ctxc, Acoeffs + N*Ai);

        H.process_store_mul(B.exps, C.exps, C.len, M);
    }

    A.set_length(Ai);
}

template <bool MP = true, ulimb ML = 0, typename RingT, typename MPolyT>
bool _generic_mpoly_divides_heap(
    mpoly_ctx& ctxm, RingT& ctxc,
    MPolyT& Q,
    typename MPolyT::input_t A,
    typename MPolyT::input_t B,
    tmp_allocator& push)
{
    FLINT_ASSERT(!MP || (Q.bits()%FLINT_BITS == 0));
    FLINT_ASSERT(A.len > 0);
    FLINT_ASSERT(B.len > 0);

    ulimb M = ctxm.stride(Q.bits());
    ulimb N = ctxc.stride();
    mpoly_heap<MP, ML> H(B.len, M, push);
    typename RingT::dotter_t dotter(ctxc);
    typename RingT::divider_t divider(ctxc);
    ulimb s = B.len; // number of terms * (latest quotient) we should put into heap
    ulimb mask = ctxm.overflow_mask<MP>(Q.bits());

    divider.set_divisor(ctxc, B.coeffs + N*0, true); // can throw
    H.start_div(A.exps, M);

    ulimb Qi = 0;
    while (!H.is_empty())
    {
        ulimb* Qexps = Q.m.fit_alloc(M*(Qi + 2));
        ulimb* texp = Qexps + M*(Qi + 1);  // use next for tmp space
        if (mpoly_monomial_set_overflowed<MP, ML>(texp, H.top_exp(), M, mask))
            goto not_exact_division;

        bool lm_not_divides = mpoly_monomial_sub_overflowed<MP, MP>(Qexps + M*Qi, H.top_exp(), B.exps + M*0, M, mask);

        typename RingT::coeff_t* Qcoeffs = Q.c.fit_alloc(N*(Qi + 1));
        dotter.zero(ctxc);
        do {
            auto x = H.pop(M);
            do {
                *H.store++ = x->i;
                *H.store++ = x->j;
                if (x->i + 1 == 0) {
                    FLINT_ASSERT(x->j < A.len);
                    dotter.sub(ctxc, A.coeffs + N*x->j);
                } else {
                    FLINT_ASSERT(0 < x->i && x->i < B.len && x->j < Qi);
                    H.ind[x->i] |= 1;
                    dotter.madd(ctxc, B.coeffs + N*x->i, Qcoeffs + N*x->j);
                }
            } while ((x = x->next) != NULL);
        } while (!H.is_empty() && mpoly_monomial_equal(H.top_exp(), texp, M));
        dotter.reduce(ctxc, Qcoeffs + N*Qi);

        s += H.process_store_div(A.exps, A.len, B.exps, Qexps, Qi, M);

        if (!divider.divides(ctxc, Qcoeffs + N*Qi, Qcoeffs + N*Qi))
            goto not_exact_division;

        if (is_zero(ctxc, Qcoeffs + N*Qi))
            continue;

        if (lm_not_divides)
            goto not_exact_division;

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
            H.insert_mul(B.exps, 1, Qexps, Qi, M);

        s = 1;
        Qi++;
    }

    Q.set_length(Qi);
    return true;

not_exact_division:
    Q.set_length(0);
    return false;
}

