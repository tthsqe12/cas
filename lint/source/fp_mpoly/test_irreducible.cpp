#include "limbs.h"
#include "fp_mpoly.h"

static int _point2d_cmp(slimb x1, slimb y1, slimb x2, slimb y2)
{
    if (x1 < x2)
        return -1;
    if (x1 > x2)
        return 1;
    if (y1 < y2)
        return -1;
    if (y1 > y2)
        return 1;
    return 0;
}

struct point2d {
    slimb x;
    slimb y;
    bool operator <(point2d other) {return _point2d_cmp(x, y, other.x, other.y) < 0;}
};

static int point2d_cmp(point2d p, point2d q)
{
    return _point2d_cmp(p.x, p.y, q.x, q.y);
}

/*
    Standing on O and looking at A, is B strictly to the left?
    i.e. (A.y - O.y) * (B.x - O.x) - (A.x - O.x) * (B.y - O.y) < 0
*/
static int _is_ccw(point2d O, point2d A, point2d B)
{
    return si_mat22_det_is_negative(A.y - O.y, A.x - O.x,
                                    B.y - O.y, B.x - O.x);
}

#if FLINT_WANT_ASSERT
static bool point2d_set_is_canonical(const my_vec<point2d>& A)
{
    for (ulimb i = 1; i < A.length(); i++)
        if (point2d_cmp(A[i], A[i - 1]) <= 0)
            return false;

    return true;
}
#endif

static void point2d_set_sort(my_vec<point2d>& A)
{
    std::sort(A.data(), A.data() + A.length());
//    qsort(A->points, A->length, sizeof(point2d), (int(*)(const void*,const void*))point2d_cmp);
}

/*
    P is a sorted array of nP distinct points
    compute the points on a ccw traversal of the convex hull via
        Andrew's monotone chain convex hull algorithm
*/
static slimb convex_hull_ccw(slimb* idxs, const point2d* P, slimb nP)
{
    slimb i, j, k = 0;

    if (nP < 3)
    {
        for (i = 0; i < nP; i++)
            idxs[i] = i;

        return nP;
    }

    for (i = 0; i < nP; i++)
    {
        while (k >= 2 && !_is_ccw(P[idxs[k - 2]], P[idxs[k - 1]], P[i]))
            k--;
        idxs[k++] = i;
    }

    for (i = nP - 1, j = k + 1; i > 0; i--)
    {
        while (k >= j && !_is_ccw(P[idxs[k - 2]], P[idxs[k - 1]], P[i - 1]))
            k--;
        idxs[k++] = i - 1;
    }

    return k - 1;
}

static int _is_in_polygon(
    point2d* V,    /* ccw polygon is V[0] -> ... -> V[nV-1] -> V[0] */
    slimb nV,
    point2d p)
{
    slimb i, a, b, c;
#if FLINT_WANT_ASSERT
    int check;

    i = nV - 1;
    check = !_is_ccw(V[0], V[i], p);
    for (i = nV - 2; i >= 0; i--)
        check &= !_is_ccw(V[i + 1], V[i], p);
#endif

again:

    FLINT_ASSERT(nV >= 3);

    if (nV < 8)
    {
        i = nV - 1;
        if (_is_ccw(V[0], V[i], p))
        {
            FLINT_ASSERT(check == 0);
            return 0;
        }
        for (i = nV - 2; i >= 0; i--)
        {
            if (_is_ccw(V[i + 1], V[i], p))
            {
                FLINT_ASSERT(check == 0);
                return 0;
            }
        }

        FLINT_ASSERT(check == 1);
        return 1;
    }

    a = nV/4;
    b = nV/2;
    c = nV - nV/4;

    if (_is_ccw(V[a], V[0], p))
    {
        V += 0;
        nV = 1 + a;
        goto again;
    }

    if (_is_ccw(V[b], V[a], p))
    {
        V += a;
        nV = 1 + b - a;
        goto again;
    }

    if (_is_ccw(V[c], V[b], p))
    {
        V += b;
        nV = 1 + c - b;
        goto again;
    }

    if (!_is_ccw(V[0], V[c], p))
    {
        FLINT_ASSERT(check == 1);
        return 1;
    }

    if (!_is_ccw(V[nV - 1], V[c], p))
    {
        FLINT_ASSERT(check == !_is_ccw(V[0], V[nV - 1], p));
        return !_is_ccw(V[0], V[nV - 1], p);
    }

    V += c;
    nV = nV - c;

    if (nV >= 3)
        goto again;

    FLINT_ASSERT(nV == 2);
    FLINT_ASSERT(check == 0);
    return 0;
}


/* T = A union ((B + s) intersect V) */
static void point2d_set_merge_shift(
    my_vec<point2d>& T,
    const my_vec<point2d>& A,
    const my_vec<point2d>& B,
    slimb sx, slimb sy,
    point2d * V, slimb nV)
{
    FLINT_ASSERT(&T != &A);
    FLINT_ASSERT(&T != &B);

    ulimb Alen = A.length();
    ulimb Blen = B.length();
    const point2d* Apoints = A.data();
    const point2d* Bpoints = B.data();
          point2d* Tpoints = T.fit_alloc_destroy(Alen + Blen);

    ulimb i = 0, j = 0, k = 0;
    while (i < Alen && j < Blen)
    {
        slimb Bsx = Bpoints[j].x + sx;
        slimb Bsy = Bpoints[j].y + sy;
        int cmp = _point2d_cmp(Apoints[i].x, Apoints[i].y, Bsx, Bsy);

        if (cmp < 0)
        {
            Tpoints[k] = Apoints[i];
            i++;
            k++;
        }
        else if (cmp > 0)
        {
            Tpoints[k].x = Bsx;
            Tpoints[k].y = Bsy;
            j++;
            k += _is_in_polygon(V, nV, Tpoints[k]);
        }
        else
        {
            Tpoints[k] = Apoints[i];
            i++;
            j++;
            k++;
        }

    }

    while (i < Alen)
    {
        Tpoints[k] = Apoints[i];
        i++;
        k++;
    }

    while (j < Blen)
    {
        Tpoints[k].x = Bpoints[j].x + sx;
        Tpoints[k].y = Bpoints[j].y + sy;
        j++;
        k += _is_in_polygon(V, nV, Tpoints[k]);
    }

    T.set_length(k);

    FLINT_ASSERT(point2d_set_is_canonical(T));
}

/* T = A union B */
static void point2d_set_merge(
    my_vec<point2d>& T,
    const my_vec<point2d>& A,
    const my_vec<point2d>& B)
{
    FLINT_ASSERT(&T != &A);
    FLINT_ASSERT(&T != &B);

    ulimb Alen = A.length();
    ulimb Blen = B.length();
    const point2d* Apoints = A.data();
    const point2d* Bpoints = B.data();
          point2d* Tpoints = T.fit_alloc_destroy(Alen + Blen);

    ulimb i = 0, j = 0, k = 0;
    while (i < Alen && j < Blen)
    {
        int cmp = _point2d_cmp(Apoints[i].x, Apoints[i].y,
                               Bpoints[j].x, Bpoints[j].y);
        if (cmp < 0)
        {
            Tpoints[k] = Apoints[i];
            i++;
            k++;
        }
        else if (cmp > 0)
        {
            Tpoints[k] = Bpoints[j];
            j++;
            k += 1;
        }
        else
        {
            Tpoints[k] = Apoints[i];
            i++;
            j++;
            k++;
        }
    }

    while (i < Alen)
    {
        Tpoints[k] = Apoints[i];
        i++;
        k++;
    }

    while (j < Blen)
    {
        Tpoints[k] = Bpoints[j];
        j++;
        k += 1;
    }

    T.set_length(k);

    FLINT_ASSERT(point2d_set_is_canonical(T));
}



#if FLINT_WANT_ASSERT
static bool point2d_set_contains(const my_vec<point2d>& A, slimb x, slimb y)
{
    slimb lo = 0;
    slimb mid;
    slimb hi = A.length();
    const point2d* Apoints = A.data();
    int cmp;

again:

    if (hi - lo < 8)
    {
        for ( ; lo < hi; lo++)
        {
            if (Apoints[lo].x == x && Apoints[lo].y == y)
                return 1;
        }
        return 0;
    }

    mid = lo + (hi - lo)/2;

    cmp = _point2d_cmp(Apoints[mid].x, Apoints[mid].y, x, y);

    if (cmp == 0)
        return 1;

    if (cmp < 0)
        lo = mid;
    else
        hi = mid;

    goto again;
}
#endif

/* is A intersect B empty? */
int point2d_set_disjoint(
    const my_vec<point2d>& A,
    const my_vec<point2d>& B)
{
    const point2d* Apoints = A.data();
    const point2d* Bpoints = B.data();
    slimb Alen = A.length();
    slimb Blen = B.length();
    slimb lo, mid, hi;
    int cmp;
#if FLINT_WANT_ASSERT
    int check = 1;

    for (lo = 0; lo < Blen; lo++)
        check &= !point2d_set_contains(A, Bpoints[lo].x, Bpoints[lo].y);

#endif

again:

    if (Alen < 1 || Blen < 1)
    {
        FLINT_ASSERT(check == 1);
        return 1;
    }

    if (Alen < Blen)
    {
        std::swap(Alen, Blen);
        std::swap(Apoints, Bpoints);
    }

    cmp = point2d_cmp(Bpoints[0], Apoints[0]);

    if (cmp == 0)
    {
        FLINT_ASSERT(check == 0);
        return 0;
    }

    if (cmp < 0)
    {
        Bpoints += 1;
        Blen -= 1;
        goto again;
    }

    /*
        throw out everything from A that is < B[0]
        if A contains B[0], return 0
    */

    lo = 0;
    hi = Alen - 1;

    cmp = point2d_cmp(Bpoints[0], Apoints[hi]);
    if (cmp >= 0)
    {
        FLINT_ASSERT(cmp == check);
        return cmp;
    }

search:

    FLINT_ASSERT(point2d_cmp(Apoints[lo], Bpoints[0]) < 0);
    FLINT_ASSERT(point2d_cmp(Bpoints[0], Apoints[hi]) < 0);

    if (hi - lo < 8)
    {
        for (lo++ ; lo < hi; lo++)
        {
            cmp = point2d_cmp(Bpoints[0], Apoints[lo]);

            if (cmp == 0)
            {
                FLINT_ASSERT(check == 0);
                return 0;
            }

            if (cmp < 0)
                break;
        }

        Apoints += lo;
        Alen -= lo;
        Bpoints += 1;
        Blen -= 1;
        goto again;
    }

    mid = lo + (hi - lo)/2;

    cmp = point2d_cmp(Apoints[mid], Bpoints[0]);

    if (cmp == 0)
    {
        FLINT_ASSERT(check == 0);
        return 0;
    }

    if (cmp < 0)
        lo = mid;
    else
        hi = mid;

    goto search;
}


/*
    ccw polygon is V[0] -> ... -> V[nV-1] -> V[0]
    |verts coordinates| < 2^(FLINT_BITS - 3)

    return: 0 no
            1 yes
           -1 don't know
*/
static int convex_hull_is_indecomposable(
    point2d* V,
    slimb nV,
    ulimb bound,
    my_vec<point2d>& Ai,   /* tmp storage */
    my_vec<point2d>& Aim1,
    my_vec<point2d>& T,
    my_vec<point2d>& S,
    point2d* E,        /* tmp of length nV */
    slimb* Egcd)       /* tmp of length nV */
{
    slimb i, j, k, g, prevx, prevy;

    FLINT_ASSERT(nV >= 3);

    if (nV == 3)
    {
        ulimb gg = my_abs(V[2].x - V[0].x);
        gg = n_gcd(gg, my_abs(V[2].y - V[0].y));
        gg = n_gcd(gg, my_abs(V[1].x - V[0].x));
        gg = n_gcd(gg, my_abs(V[1].y - V[0].y));
        return gg == 1;
    }

    /*
        since |V[i]| < 2^(FLINT_BITS - 3), all |E[i]| < 2^(FLINT_BITS - 2) and
        no addition V[i] + E[j] will overflow
    */
    prevx = V[0].x;
    prevy = V[0].y;
    g = 0;
    ulimb prod = 1;
    for (i = nV - 1; i >= 0; i--)
    {
        E[i].x = prevx - V[i].x;
        E[i].y = prevy - V[i].y;
        prevx = V[i].x;
        prevy = V[i].y;
        Egcd[i] = n_gcd(my_abs(E[i].x), my_abs(E[i].y));
        E[i].x /= Egcd[i];
        E[i].y /= Egcd[i];
        g = n_gcd(g, Egcd[i]);
        if (ui_mul_overflowed(prod, prod, Egcd[i]))
            return -1;
    }

    if (g > 1)
        return 0;

    if (prod > bound)
        return -1;

    /* S = {V[0] + j*E[nV-1]}_j */
    point2d* Spoints = S.fit_alloc(Egcd[nV - 1]);
    for (j = 0; j < Egcd[nV - 1]; j++)
    {
        Spoints[j].x = V[0].x - j*E[nV - 1].x;
        Spoints[j].y = V[0].y - j*E[nV - 1].y;
    }
    S.set_length(Egcd[nV - 1]);
    point2d_set_sort(S);

    /* A_{i-1} is empty */
    Aim1.set_length(0);

    for (i = 0; i < nV - 1; i++)
    {
        point2d* Aipoints = Ai.fit_alloc(Egcd[i]);
        k = 0;
        for (j = 1; j <= Egcd[i]; j++)
        {
            Aipoints[k].x = V[0].x + j*E[i].x;
            Aipoints[k].y = V[0].y + j*E[i].y;
            if (!_is_in_polygon(V, nV, Aipoints[k]))
                break;
            k++;
        }
        Ai.set_length(k);
        point2d_set_sort(Ai);

        if (Aim1.length() > 0)
        {
            point2d_set_merge(T, Ai, Aim1);
            std::swap(Ai, T);
            for (j = 1; j <= Egcd[i]; j++)
            {
                point2d_set_merge_shift(T, Ai, Aim1, j*E[i].x, j*E[i].y, V, nV);
                std::swap(Ai, T);
                if (!point2d_set_disjoint(Ai, S))
                    return 0;
            }
        }
        else
        {
            if (!point2d_set_disjoint(Ai, S))
                return 0;
        }

        std::swap(Aim1, Ai);
    }

    return 1;
}

static void si_rand_vec_primitive(
    random_state& state,
    slimb * v, slimb len,
    ulimb bound)
{
    slimb i, g;

again:

    g = 0;
    for (i = 0; i < len; i++)
    {
        ulimb u = state.get_mod(bound);
        v[i] = state.get_bit() ? -slimb(u) : slimb(u);
        g = n_gcd(g, u);
    }

    if (g == 0)
        goto again;

    if (g == 1)
        return;

    for (i = 0; i < len; i++)
        v[i] /= g;
}


static int _test_indecomposable2(slimb* a, slimb* b, ulimb n)
{
    ulimb g = 0;
    for (ulimb i = 0; i < n; i++)
        g = n_gcd(g, my_abs(a[i] - b[i]));
    return g == 1;
}


static int _test_colinear(slimb* a, slimb* b, slimb* c, slimb n)
{
    fmpz tn, td, sn, sd, g;

    for (ulimb i = 0; i < n; i++)
    {
        fmpz_set_si(sn, a[i]);
        fmpz_sub_si(sn, sn, c[i]);
        fmpz_set_si(sd, a[i]);
        fmpz_sub_si(sd, sd, b[i]);
        fmpz_gcd(g, sn, sd);
        if (fmpz_is_zero(g))
            continue;

        if (sd.is_negative())
            fmpz_neg(g, g);
        fmpz_divexact(sn, sn, g);
        fmpz_divexact(sd, sd, g);

        if (fmpz_is_zero(td))
        {
            fmpz_swap(td, sd);
            fmpz_swap(tn, sn);
        }
        else if (!fmpz_equal(sd, td) || !fmpz_equal(sn, tn))
        {
            return 4;
        }
    }

    if (fmpz_is_zero(td))
        return 0;
    else if (fmpz_sgn(tn) < 0)
        return 1;
    else if (fmpz_cmp(tn, td) > 0)
        return 2;
    else
        return 3;
}

static int _test_indecomposable3(
    slimb * a,
    slimb * b,
    slimb * c,
    slimb n)
{
    slimb i;
    ulimb g;

    switch (_test_colinear(a, b, c, n))
    {
        case 0:
            return 0;
        case 1:
            return _test_indecomposable2(c, b, n);
        case 2:
            return _test_indecomposable2(a, c, n);
        case 3:
            return _test_indecomposable2(a, b, n);
        default:
            break;
    }

    g = 0;
    for (i = 0; i < n; i++)
    {
        g = n_gcd(g, FLINT_ABS(a[i] - b[i]));
        g = n_gcd(g, FLINT_ABS(a[i] - c[i]));
    }
    return g == 1;
}

/*
    Fast Absolute Irreducibility Testing via Newton Polytopes
    Shuhong Gao and Alan G.B. Lauderz

    Aexps is the Alen x nvars exponent matrix with the entry A[i,j]
    at Aexps[i*stride + j]
*/
int _mpoly_test_irreducible(
    slimb * Aexps, slimb stride, slimb Alen,
    slimb nvars,
    slimb tries_left,   /* what they call the "projection bound" */
    random_state& state)
{
    int success;
    slimb i, j, newlen, hull_len;
    ulimb matrix_bound = 2;
    ulimb memory_bound = 1000;
    ulimb max_memory_bound = pow2(20 + FLINT_BITS/8);

    if (Alen < 2 || nvars < 2)
        return false;
    if (Alen == 2)
        return _test_indecomposable2(Aexps + 0*stride, Aexps + 1*stride, nvars);
    if (Alen == 3)
        return _test_indecomposable3(Aexps + 0*stride, Aexps + 1*stride, Aexps + 2*stride, nvars);
    if (tries_left <= 0)
        return false;

    my_vec<point2d> T1, T2, T3, T4;
    chunk<slimb> hull_idxs(10);
    chunk<point2d> hull_points(4);
    tmp_allocator push;

    slimb* rowx = push.recursive_alloc<slimb>(2*nvars);
    slimb* rowy = rowx + nvars;
    slimb* dups = FLINT_ARRAY_ALLOC(Alen, slimb);
    point2d* points = push.recursive_alloc<point2d>(Alen);

again:

    if (--tries_left < 0)
        return false;

    memory_bound = std::min(max_memory_bound, memory_bound/8*9);
    matrix_bound += 1;

    if (nvars == 2)
    {
        tries_left = 0;
        memory_bound = max_memory_bound;

        for (i = 0; i < Alen; i++)
        {
            slimb lox = Aexps[i*stride + 0];
            slimb loy = Aexps[i*stride + 1];

            if (std::min(lox, loy) <= -(WORD(1) << (FLINT_BITS - 3)) ||
                std::max(lox, loy) >= WORD(1) << (FLINT_BITS - 3))
            {
                return false;
            }

            points[i].x = lox;
            points[i].y = loy;
        }
    }
    else
    {
        si_rand_vec_primitive(state, rowx, nvars, matrix_bound);
        si_rand_vec_primitive(state, rowy, nvars, matrix_bound);

        for (i = 0; i < Alen; i++)
        {
            ulimb x2, x1, x0, y2, y1, y0, p2, p1, p0;
            x2 = x1 = x0 = y2 = y1 = y0 = 0;
            for (j = 0; j < nvars; j++)
            {
                SMUL_PPMM(p1, p0, Aexps[i*stride + j], rowx[j]);
                p2 = sign_extend(p1);
                ADD_SSSAAAAAA(x2, x1, x0, x2, x1, x0, p2, p1, p0);

                SMUL_PPMM(p1, p0, Aexps[i*stride + j], rowy[j]);
                p2 = sign_extend(p1);
                ADD_SSSAAAAAA(y2, y1, y0, y2, y1, y0, p2, p1, p0);
            }

            if (x2 != sign_extend(x0) || x1 != sign_extend(x0) ||
                y2 != sign_extend(y0) || y1 != sign_extend(y0))
            {
                goto again;
            }

            points[i].x = x0;
            points[i].y = y0;

            if (points[i].x <= -(WORD(1) << (FLINT_BITS - 3)) ||
                points[i].y <= -(WORD(1) << (FLINT_BITS - 3)) ||
                points[i].x >= WORD(1) << (FLINT_BITS - 3) ||
                points[i].y >= WORD(1) << (FLINT_BITS - 3))
            {
                goto again;
            }
        }
    }

//    qsort(points, Alen, sizeof(point2d), (int(*)(const void*,const void*))point2d_cmp);
    std::sort(points + 0, points + Alen);

    /* delete duplicates and track which are duplicated */
    dups[0] = 0;
    newlen = 1;
    for (i = 1; i < Alen; i++)
    {
        if (point2d_cmp(points[newlen - 1], points[i]) == 0)
        {
            dups[newlen - 1] = 1;
        }
        else
        {
            dups[newlen] = 0;
            points[newlen] = points[i];
            newlen++;
        }
    }

    /* find indices of convex hull */
    hull_len = convex_hull_ccw(hull_idxs.fit_alloc(newlen + 1), points, newlen);
    if (hull_len < 3)
        goto again;

    /* ensure no duplicates on hull */
    hull_points.fit_alloc(hull_len + 1);
    for (i = 0; i < hull_len; i++)
    {
        hull_points[i] = points[hull_idxs[i]];
        if (dups[hull_idxs[i]] != 0)
            goto again;
    }

    /* check indecomposability with a bound on the memory useage */
    success = convex_hull_is_indecomposable(hull_points.data(), hull_len,
                                   memory_bound, T1, T2, T3, T4, points, dups);
    if (success < 1)
    {
        if (success < 0)
            memory_bound = FLINT_MIN(max_memory_bound, memory_bound/8*9);
        goto again;
    }

    return true;
}


bool mpoly_proved_irreducible(
    mpoly_ctx& ctx,
    mpolyc A,
    random_state& state)
{
    slimb n = ctx.nvars();
    slimb i, j;

    if (A.bits > FLINT_BITS || A.len < 2)
        return false;

    tmp_allocator push;
    slimb* max_exps = push.recursive_alloc<slimb>(n + n*A.len);
    slimb* uexps = max_exps + n;

    for (j = 0; j < n; j++)
        max_exps[j] = 0;

    slimb N = ctx.stride_sp(A.bits);
    auto offsets = ctx.var_offsets_sp(A.bits);
    auto shifts = ctx.var_shifts_sp(A.bits);
    ulimb mask = ctx.var_mask_sp(A.bits);
    for (i = 0; i < A.len; i++)
    for (j = 0; j < n; j++)
    {
        uexps[n*i + j] = ((A.exps + N*i)[offsets[j]] >> shifts[j]) & mask;
        max_exps[j] = std::max(max_exps[j], uexps[n*i + j]);
    }

    slimb tdeg = 1;
    for (j = 0; j < n; j++)
    {
        if (si_add_overflowed(tdeg, tdeg, max_exps[j])) {
            tdeg = WORD_MAX;
            break;
        }
    }

    slimb tries = 12 - A.len/tdeg/2;
    return _mpoly_test_irreducible(uexps, n, A.len, n, tries, state);
}

