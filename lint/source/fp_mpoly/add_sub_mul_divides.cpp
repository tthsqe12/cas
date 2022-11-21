#include "fp_mpoly.h"
#include "generic/generic_mpoly_scalar.h"
#include "generic/generic_mpoly_mul_divides.h"

void fp_mpoly_set(fp_mpoly_ring& ctx,
    fp_mpoly& A,
    const fp_mpoly& B)
{
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, B));

    if (&A == &B)
        return;

    ulimb len = B.length();
    if (len < 1)
    {
        fp_mpoly_zero(ctx, A);
        return;
    }

    ulimb bits = B.bits();
    A.set_bits(bits);
    A.set_length(len);

    ulimb M = ctx.m.stride(bits);
    ulimb* Aexps = A.m.fit_alloc_destroy(len, M);
    mpoly_copy_monomials(Aexps, B.exps(), len, M);

    ulimb N = ctx.c.stride();
    ulimb* Acoeffs = A.c.fit_alloc_destroy(len, N);
    fp_vec_set(ctx.c, Acoeffs, B.coeffs(), len);

    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, A));
}

void fp_mpoly_neg(fp_mpoly_ring& ctx,
    fp_mpoly& A,
    const fp_mpoly& B)
{
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, B));

    ulimb len = B.length();

    if (&A == &B)
    {
        fp_vec_neg(ctx.c, A.coeffs(), B.coeffs(), len);
        return;
    }

    if (len == 0)
    {
        fp_mpoly_zero(ctx, A);
        return;
    }

    ulimb bits = B.bits();
    A.set_bits(bits);
    A.set_length(len);

    ulimb M = ctx.m.stride(bits);
    ulimb* Aexps = A.m.fit_alloc_destroy(len, M);
    mpoly_copy_monomials(Aexps, B.exps(), len, M);

    ulimb N = ctx.c.stride();
    ulimb* Acoeffs = A.c.fit_alloc_destroy(len, N);
    fp_vec_neg(ctx.c, Acoeffs, B.coeffs(), len);

    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, A));
}

bool fp_mpoly_equal(fp_mpoly_ring& ctx,
    const fp_mpoly& B,
    const fp_mpoly& C)
{
    ulimb N = ctx.c.stride();
    ulimb l = B.length();

    if (B.length() != C.length())
        return false;

    if (l < 1)
        return true;

    for (ulimb i = 0; i < l; i++) {
        if (!fp_equal(ctx.c, B.coeffs() + N*i, C.coeffs() + N*i))
           return false;
    }

    tmp_allocator push;

    ulimb Abits = std::max(B.bits(), C.bits());
    ulimb M = ctx.m.stride(Abits);

    const ulimb* Bexps = B.exps();
    if (Abits != B.bits())
    {
        ulimb* t = push.recursive_alloc<ulimb>(M*l);
        Bexps = t;
        mpoly_repacked_up_monomials(ctx.m, t, Abits, B.m);
    }

    const ulimb* Cexps = C.exps();
    if (Abits != C.bits())
    {
        ulimb* t = push.recursive_alloc<ulimb>(M*l);
        Cexps = t;
        mpoly_repacked_up_monomials(ctx.m, t, Abits, C.m);
    }

    for (ulimb i = 0; i < M*l; i++) {
        if (Bexps[i] != Cexps[i])
            return false;
    }

    return true;
}

void fp_mpoly_add(fp_mpoly_ring& ctx,
    fp_mpoly& A,
    const fp_mpoly& B,
    const fp_mpoly& C)
{
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, B));
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, C));

    FLINT_ASSERT(&A != &B); FLINT_ASSERT(&A != &C);

    if (UNLIKELY(B.length() == 0))
    {
        fp_mpoly_set(ctx, A, C);
        return;
    }
    else if (UNLIKELY(C.length() == 0))
    {
        fp_mpoly_set(ctx, A, B);
        return;
    }

    tmp_allocator push;

    ulimb Abits = std::max(B.bits(), C.bits());
    A.set_bits(Abits);
    ulimb M = ctx.m.stride(Abits);

    fp_mpolyi Bi(B);
    if (Abits != B.bits())
    {
        Bi.exps = push.recursive_alloc<ulimb>(M*Bi.len);
        mpoly_repacked_up_monomials(ctx.m, Bi.exps, Abits, B.m);
    }

    fp_mpolyi Ci(C);
    if (Abits != C.bits())
    {
        Ci.exps = push.recursive_alloc<ulimb>(M*Ci.len);
        mpoly_repacked_up_monomials(ctx.m, Ci.exps, Abits, C.m);
    }

    if (M == 1)
        _generic_mpoly_add<1, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci);
    else if (M == 2)
        _generic_mpoly_add<2, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci);
    else
        _generic_mpoly_add<0, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci);

    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, A));
}



void fp_mpoly_sub(fp_mpoly_ring& ctx,
    fp_mpoly& A,
    const fp_mpoly& B,
    const fp_mpoly& C)
{
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, B));
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, C));

    FLINT_ASSERT(&A != &B); FLINT_ASSERT(&A != &C);

    if (UNLIKELY(B.length() == 0))
    {
        fp_mpoly_neg(ctx, A, C);
        return;
    }
    else if (UNLIKELY(C.length() == 0))
    {
        fp_mpoly_set(ctx, A, B);
        return;
    }

    tmp_allocator push;

    ulimb Abits = std::max(B.bits(), C.bits());
    A.set_bits(Abits);
    ulimb M = ctx.m.stride(Abits);

    fp_mpolyi Bi(B);
    if (Abits != B.bits())
    {
        Bi.exps = push.recursive_alloc<ulimb>(M*Bi.len);
        mpoly_repacked_up_monomials(ctx.m, Bi.exps, Abits, B.m);
    }

    fp_mpolyi Ci(C);
    if (Abits != C.bits())
    {
        Ci.exps = push.recursive_alloc<ulimb>(M*Ci.len);
        mpoly_repacked_up_monomials(ctx.m, Ci.exps, Abits, C.m);
    }

    if (M == 1)
        _generic_mpoly_sub<1, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci);
    else if (M == 2)
        _generic_mpoly_sub<2, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci);
    else
        _generic_mpoly_sub<0, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci);

    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, A));
}

void _fp_mpoly_mul_heap_have_degrees(fp_mpoly_ring& ctx,
    fp_mpoly& A, const fmpz* Adegs,
    const fp_mpoly& B,
    const fp_mpoly& C)
{
    FLINT_ASSERT(B.length() > 0);
    FLINT_ASSERT(C.length() > 0);
    FLINT_ASSERT(&A != &B); FLINT_ASSERT(&A != &C);

    tmp_allocator push;

    ulimb Abits = fmpz_vec_max_bits(Adegs, ctx.m.nvars());
    FLINT_ASSERT(!has_top_bit(Abits));
    Abits = ctx.m.fix_bits(std::max(1 + Abits, std::max(B.bits(), C.bits())));
    ulimb M = ctx.m.stride(Abits);

    fp_mpolyi Bi(B);
    if (Abits != B.bits())
    {
        Bi.exps = push.recursive_alloc<ulimb>(M*Bi.len);
        mpoly_repacked_up_monomials(ctx.m, Bi.exps, Abits, B.m);
    }

    fp_mpolyi Ci(C);
    if (Abits != C.bits())
    {
        Ci.exps = push.recursive_alloc<ulimb>(M*Ci.len);
        mpoly_repacked_up_monomials(ctx.m, Ci.exps, Abits, C.m);
    }

    A.set_bits(Abits);
    if (Bi.len > Ci.len)
        std::swap(Bi, Ci);

    if (Abits > FLINT_BITS)
        _generic_mpoly_mul_heap<true, 0, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci, push);
    else if (M == 1)
        _generic_mpoly_mul_heap<false, 1, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci, push);
    else if (M == 2)
        _generic_mpoly_mul_heap<false, 2, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci, push);
    else
        _generic_mpoly_mul_heap<false, 0, fp_ring, fp_mpoly>(ctx.m, ctx.c, A, Bi, Ci, push);
}

void fp_mpoly_mul(fp_mpoly_ring& ctx,
    fp_mpoly& A,
    const fp_mpoly& B,
    const fp_mpoly& C)
{
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, B));
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, C));

    if (UNLIKELY(B.length() < 1 || C.length() < 1))
    {
        fp_mpoly_zero(ctx, A);
        return;
    }

    fmpz* Bdegs = ctx.m.leaf_alloc_fmpz(2*ctx.m.nvars());
    fmpz* Cdegs = Bdegs + 1*ctx.m.nvars();
    mpoly_degrees(ctx.m, Bdegs, B.m);
    mpoly_degrees(ctx.m, Cdegs, C.m);
    fmpz_vec_add(Bdegs, Bdegs, Cdegs, ctx.m.nvars());
    _fp_mpoly_mul_heap_have_degrees(ctx, A, Bdegs, B, C);

    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, A));
}

bool fp_mpoly_divides(fp_mpoly_ring& ctx,
    fp_mpoly& Q,
    const fp_mpoly& A,
    const fp_mpoly& B)
{
    FLINT_ASSERT(&Q != &A);
    FLINT_ASSERT(&Q != &B);
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, A));
    FLINT_ASSERT(fp_mpoly_is_canonical(ctx, B));

    if (UNLIKELY(A.length() < 1))
    {
        fp_mpoly_zero(ctx, Q);
        return true;
    }
    else if (UNLIKELY(B.length() < 1))
    {
        fp_mpoly_zero(ctx, Q);
        return false;
    }

    tmp_allocator push;

    ulimb Qbits = A.bits();
    ulimb M = ctx.m.stride(Qbits);
    Q.set_bits(Qbits);
    fp_mpolyi Ai(A);
    fp_mpolyi Bi(B);
    if (Qbits != B.bits())
    {
        Bi.exps = push.recursive_alloc<ulimb>(M*Bi.len);
        if (!mpoly_repacked_monomials(ctx.m, Bi.exps, Qbits, B.m))
            return false;
    }

    if (Qbits > FLINT_BITS)
        return _generic_mpoly_divides_heap<true, 0, fp_ring, fp_mpoly>(ctx.m, ctx.c, Q, Ai, Bi, push);
    else if (M == 1)
        return _generic_mpoly_divides_heap<false, 1, fp_ring, fp_mpoly>(ctx.m, ctx.c, Q, Ai, Bi, push);
    else if (M == 2)
        return _generic_mpoly_divides_heap<false, 2, fp_ring, fp_mpoly>(ctx.m, ctx.c, Q, Ai, Bi, push);
    else
        return _generic_mpoly_divides_heap<false, 0, fp_ring, fp_mpoly>(ctx.m, ctx.c, Q, Ai, Bi, push);
}

