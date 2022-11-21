#include "fp_poly.h"


void fp_berlekamp_massey::points_fit_alloc(ulimb l, ulimb N) {

    if (points_alloc >= l*N)
        return;

    ulimb new_points_alloc = std::max(l*N, points_alloc + points_alloc/2);
    ulimb* new_points_data = my_alloc<ulimb>(new_points_alloc);

    for (ulimb i = 0; i < points_length; i++)
    {
        ui_vec_set(new_points_data + new_points_alloc - N*(i+1),
                     points_data + points_alloc - N*(i+1), N);
    }

    my_free(points_data);
    points_data = new_points_data;
    points_alloc = new_points_alloc;
}

fp_berlekamp_massey::~fp_berlekamp_massey()
{
    my_free(points_data);
}

fp_berlekamp_massey::fp_berlekamp_massey()
{
    npoints = 0;
    points_data = nullptr;
    points_alloc = 0;
    points_length = 0;
}

fp_berlekamp_massey::fp_berlekamp_massey(fp_ring& ctx)
{
    npoints = 0;
    points_data = nullptr;
    points_alloc = 0;
    points_length = 0;
    fp_poly_one(ctx, R0);
    fp_poly_one(ctx, V1);
}

void fp_berlekamp_massey::start_over(fp_ring& ctx)
{
    npoints = 0;
    points_length = 0;
    fp_poly_zero(ctx, V0);
    fp_poly_one(ctx, R0);
    fp_poly_one(ctx, V1);
    fp_poly_zero(ctx, R1);
}

std::ostream& fp_berlekamp_massey::write(
    std::ostream& o,
    const fp_ring_base& ctx,
    const char* var) const
{
    o << "bma[";
    fp_poly_write(o, ctx, V1, var);
    o << ", {";
    for (slimb i = 0; i < points_length; i++)
    {
        if (i > 0)
            o << ", ";
        fp_write(o, ctx, point(i, ctx.stride()));
    }
    o << "}]";
    return o;
}

void fp_berlekamp_massey::add_points(
    const fp_ring_base& ctx,
    const ulimb* a,
    ulimb count)
{
    ulimb N = ctx.stride();
    slimb old_length = points_length;
    points_fit_alloc(old_length + count, N);
    for (ulimb i = 0; i < count; i++)
        fp_set(ctx, point(old_length + i, N), a + N*i);
    points_length = old_length + count;
}

void fp_berlekamp_massey::add_zeros(
    fp_ring& ctx,
    ulimb count)
{
    ulimb N = ctx.stride();
    slimb old_length = points_length;
    points_fit_alloc(old_length + count, N);
    for (ulimb i = 0; i < count; i++)
        fp_zero(ctx, point(old_length + i, N));
    points_length = old_length + count;
}

void fp_berlekamp_massey::add_point(fp_ring& ctx,
    const ulimb* a)
{
    ulimb N = ctx.stride();
    slimb old_length = points_length;
    points_fit_alloc(old_length + 1, N);
    FLINT_ASSERT(fp_is_canonical(ctx, a));
    fp_set(ctx, point(old_length, N), a);
    points_length = old_length + 1;
}

void fp_berlekamp_massey::add_point_ui(fp_ring& ctx,
    ulong a)
{
    ulimb N = ctx.stride();
    slimb old_length = points_length;
    points_fit_alloc(old_length + 1, N);
    fp_set_ui(ctx, point(old_length, N), a);
    points_length = old_length + 1;
}

// x = x mod X^k
void fp_poly_truncate_inplace(fp_ring& ctx, fp_poly& x, ulimb k)
{
    if (k >= x.length())
        return;
    x.normalize_length(k, ctx.stride());
}

/* return true if reduction changed the master poly, false otherwise */
bool fp_berlekamp_massey::reduce(fp_ring& ctx)
{
    slimb l, k, queue_len, queue_lo, queue_hi;

    // point(j) for queue_lo <= j < queue_hi needs to be added
    queue_lo = npoints;
    queue_hi = points_length;
    FLINT_ASSERT_ALWAYS(queue_hi >= queue_lo);
    queue_len = queue_hi - queue_lo;
    npoints = queue_hi;

    // rev = point(queue_hi)*X^0 + point(queue_hi-1)*X^1 + ...
    _fp_poly rev;
    rev.coeffs = point(queue_lo + queue_len - 1, ctx.stride());
    rev.len = fp_vec_normalized_length(ctx, rev.coeffs, queue_len);

    /* Ri = Ri * x^queue_len + Vi*rev */
    mul(ctx, qt, V0, rev);
    shift_left(ctx, R0, queue_len);
    add(ctx, R0, qt);
    mul(ctx, qt, V1, rev);
    shift_left(ctx, R1, queue_len);
    add(ctx, R1, qt);

    /* now start reducing R0, R1 */
    if (2*R1.degree() < npoints)
        return false;

    /* one iteration to get deg(R0) >= npoints/2 */
    divrem(ctx, qt, R0, R1);
    swap(ctx, R0, R1);
    submul(ctx, V0, qt, V1);
    swap(ctx, V0, V1);

    l = R0.degree();
    FLINT_ASSERT(npoints <= 2*l && l < npoints);

    k = npoints - l;
    FLINT_ASSERT(0 <= k && k <= l);

    // (l - k)/2 is the expected number of required iterations.
    // Either branch is OK anytime.
    if (l - k < 10)
    {
        while (npoints <= 2*R1.degree())
        {
            divrem(ctx, qt, R0, R1);
            swap(ctx, R0, R1);
            submul(ctx, V0, qt, V1);
            swap(ctx, V0, V1);
        }
    }
    else
    {
        fp_poly m11, m12, m21, m22, r0, r1, t0, t1, rt;

        shift_right(ctx, r0, R0, k);
        shift_right(ctx, r1, R1, k);
        int sgnM = hgcd(ctx, m11, m12, m21, m22, t0, t1, r0, r1);

        /*  multiply M^-1 . {{V0, R0}, {V1 R1}}
                M.{t0,t1} = {r0,r1}
                R0 = r0*x^k + R0'
                R1 = r1*x^k + R1'
        */

        fp_poly v0, v1;
        set(ctx, v0, V0);
        set(ctx, v1, V1);
        mul(ctx, rt, m22, v0);
        mul(ctx, qt, m12, v1);
        sgnM > 0 ? sub(ctx, V0, rt, qt) : sub(ctx, V0, qt, rt);
        mul(ctx, rt, m11, v1);
        mul(ctx, qt, m21, v0);
        sgnM > 0 ? sub(ctx, V1, rt, qt) : sub(ctx, V1, qt, rt);



        fp_poly_truncate_inplace(ctx, R0, k);
        fp_poly_truncate_inplace(ctx, R1, k);

        mul(ctx, rt, m22, R0);
        mul(ctx, qt, m12, R1);
        sgnM > 0 ? sub(ctx, r0, rt, qt) : sub(ctx, r0, qt, rt);

        mul(ctx, rt, m11, R1);
        mul(ctx, qt, m21, R0);
        sgnM > 0 ? sub(ctx, r1, rt, qt) : sub(ctx, r1, qt, rt);
        swap(ctx, R0, r0);
        swap(ctx, R1, r1);

        shift_left(ctx, t0, k);
        shift_left(ctx, t1, k);

        add(ctx, R0, R0, t0);
        add(ctx, R1, R1, t1);
    }

    FLINT_ASSERT(V1.degree() >= 0);
    FLINT_ASSERT(2*V1.degree() <= npoints);
    FLINT_ASSERT(2*R0.degree() >= npoints);
    FLINT_ASSERT(2*R1.degree() <  npoints);

    return true;
}


