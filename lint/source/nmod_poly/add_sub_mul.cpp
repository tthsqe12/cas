#include "nmod_poly.h"

inline void ui_vec_zero(ulimb* x, ulimb n)
{
    for (ulimb i = 0; i < n; i++)
        x[i] = 0;
}

inline void ui_vec_set_nz(ulimb* x, const ulimb* a, ulimb n)
{
    FLINT_ASSERT(n > 0);
    ulimb i = 0;
    do {
        x[i] = a[i];
    } while (++i < n);
}

void ui_vec_set(ulimb* x, const ulimb* a, ulimb n)
{
    if (n > 0)
        ui_vec_set_nz(x, a, n);
}

inline bool ui_vec_is_zero(const ulimb* a, ulimb N)
{
    for (ulimb i = 0; i < N; i++)
        if (a[i] != 0)
            return false;
    return true;
}

bool nmod_is_canonical(const nmod_ring_base& R, const ulimb* a)
{
    return mpn_cmp(a, R.modulus_limbs, R.stride) < 0;
}

bool nmod_is_zero(const nmod_ring_base& R, const ulimb* a)
{
    return ui_vec_is_zero(a, R.stride);
}

void nmod_neg(const nmod_ring_base& R, ulimb* x, const ulimb* a)
{
    if (ui_vec_is_zero(a, R.stride))
        ui_vec_zero(x, R.stride);
    else
        mpn_sub_n(x, R.modulus_limbs, a, R.stride);
}

void nmod_add(const nmod_ring_base& R, ulimb* x, const ulimb* a, const ulimb* b)
{
    auto c = mpn_add_n(x, a, b, R.stride);
    if (c != 0 || mpn_cmp(x, R.modulus_limbs, R.stride) >= 0)
        mpn_sub_n(x, x, R.modulus_limbs, R.stride);
}

void nmod_sub(const nmod_ring_base& R, ulimb* x, const ulimb* a, const ulimb* b)
{
    auto c = mpn_sub_n(x, a, b, R.stride);
    if (c != 0)
        mpn_add_n(x, x, R.modulus_limbs, R.stride);
}

void nmod_mul(nmod_ring& R, ulimb* x, const ulimb* a, const ulimb* b)
{
    ulimb N = R.stride;
    my_mpn_mul_n(R.temp2N, a, b, N);
    my_mpn_tdiv_qr(R.tempNp1, x, R.temp2N, 2*N, R.modulus_limbs, N);
}

void nmod_inv(nmod_ring& R, ulimb* x, const ulimb* a)
{
    ulimb N = R.stride;
    ulimb* A = R.temp2N;
    ulimb* M = A + N;
    ulimb* g = R.tempN;
    ulimb* u = R.tempNp1;
    mp_size_t usize;

    MPN_COPY(A, a, N);
    MPN_COPY(M, R.modulus_limbs, N);

    auto gn = mpn_gcdext(g, u, &usize, A, N, M, N);
    if (gn != 1 || g[0] != 1)
        throw with<nmod_ring, fmpz>(R, fmpz(g, gn, 0));

    if (usize < 0)
    {
        usize = -usize;
        FLINT_ASSERT_ALWAYS(usize <= N);
        mpn_sub(x, R.modulus_limbs, N, u, usize);
    }
    else
    {
        FLINT_ASSERT_ALWAYS(usize <= N);
        MPN_COPY(x, u, usize);
        MPN_ZERO(x + usize, N - usize);
    }
}

void nmod_divexact(nmod_ring& R, ulimb* x, const ulimb* a, const ulimb* b)
{
    ulimb* t = R.temp2Np1; // not used by nmod_inv or nmod_mul!
    nmod_inv(R, t, b);
    nmod_mul(R, x, a, t);
}

void nmod_randtest(nmod_ring& R, ulimb* x, rand_state& state)
{
    ulimb N = R.stride;
    ulimb i = 0;
    ulimb* t = R.tempNp1;
    do {
        t[i] = state.get_limb();
    } while (++i < N+1);
    ulimb q[3];
    my_mpn_tdiv_qr(q, x, t, N+1, R.modulus_limbs, N);
}


void nmod_poly::normalize(ulimb n, ulimb N)
{
    ulimb* d = coeffs.data();

    while (n > 0 && ui_vec_is_zero(d + N*(n-1), N))
        n--;

    set_length(n);
}

bool nmod_poly_is_canonical(const nmod_ring_base& R, const nmod_poly& a)
{
    ulimb N = R.stride;
    const ulimb* ad = a.data();
    ulimb an = a.length();
    for (ulimb i = an; i > 0; i--)
    {
        if (!nmod_is_canonical(R, ad + N*(i-1)))
            return false;
        if (i == an && nmod_is_zero(R, ad + N*(i-1)))
            return false;
    }
    return true;
}

void nmod_poly_write(std::ostream& o, const nmod_ring_base& R, const nmod_poly& a, const char* var)
{
    ulimb N = R.stride;
    bool first = true;
    ulimb an = a.length();
    const ulimb* ad = a.data();
    for (ulimb i = an; i > 0; i--)
    {
        if (i < an && ui_vec_is_zero(ad + N*(i-1), N))
            continue;

        if (!first)
            o << " + ";

        first = false;
        mpn_write(o, ad + N*(i-1), N);
        if (i-1 > 0)
        {
            o << "*";
            o << var;
            if (i-1 > 1)
            o << "^" << i-1;
        }
    }
    if (first)
        o << "0";
}

void nmod_poly_randtest(nmod_ring& R, nmod_poly& x, rand_state& state, ulimb len)
{
    ulimb N = R.stride;

    if (len < 1)
    {
        x.set_length(0);
        return;
    }

    ulimb* xd = x.fit_alloc(len, N);
    ulimb i = 0;
    do {
        nmod_randtest(R, xd + i*N, state);
    } while (++i < len);

    if (ui_vec_is_zero(xd + (len-1)*N, N))
        (xd + (len-1)*N)[0] = 1;

    x.set_length(len);
    return;
}

bool nmod_poly_equal(const nmod_ring_base& R, const nmod_poly& a, const nmod_poly& b)
{
    ulimb N = R.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();
    if (an != bn)
        return false;
    for (ulimb i = 0; i < N*an; i++)
        if (ad[i] != bd[i])
            return false;
    return true;
}

void nmod_poly_set(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a)
{
    ulimb N = R.stride;
    ulimb an = a.length();
    ulimb* xd = x.fit_alloc(an, N);
    ui_vec_set(xd, a.data(), N*an);    
}

void nmod_poly_add(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b)
{
    ulimb N = R.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    ulimb xn = std::max(an, bn);
    ulimb* xd = x.fit_alloc(xn, N);
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    xn = std::min(an, bn);
    for (ulimb i = 0; i < xn; i++)
        nmod_add(R, xd + N*i, ad + N*i, bd + N*i);

    if (an == bn)
    {
        x.normalize(an, N);
    }
    else if (an >= bn)
    {
        ui_vec_set_nz(xd + N*bn, ad + N*bn, N*(an - bn));
        x.set_length(an);
    }
    else
    {
        ui_vec_set(xd + N*an, bd + N*an, N*(bn - an));
        x.set_length(bn);
    }
}


void nmod_poly_sub(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b)
{
    ulimb N = R.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    ulimb xn = std::max(an, bn);
    ulimb* xd = x.fit_alloc(xn, N);
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    xn = std::min(an, bn);
    for (ulimb i = 0; i < xn; i++)
        nmod_sub(R, xd + N*i, ad + N*i, bd + N*i);

    if (an == bn)
    {
        x.normalize(an, N);
    }
    else if (an > bn)
    {
        ui_vec_set(xd + N*bn, ad + N*bn, N*(an - bn));
        x.set_length(an);
    }
    else
    {
        ulimb i = an;
        do {
            nmod_neg(R, xd + N*i, bd + N*i);
        } while (++i < bn);

        x.set_length(bn);
    }
}

void nmod_poly_mul(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b)
{
    FLINT_ASSERT(&x != &a); FLINT_ASSERT(&x != &b);

    ulimb N = R.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* ad = a.data();
    const ulimb* bd = b.data();

    if (an < 1 || bn < 1)
    {
        x.set_length(0);
        return;
    }

    ulimb xn = an + bn - 1;
    ulimb* xd = x.fit_alloc(xn, N);

    ulimb* t = R.temp2N;
    ulimb* z = R.temp2Np1;

    for (ulimb i = 0; i < xn; i++)
    {
        ui_vec_zero(z, 2*N+1);

        ulong jstart = i > bn - 1 ? i - (bn - 1) : 0;
        ulong jstop = std::min(i + 1, an);

        for (ulimb j = jstart; j < jstop; j++)
        {
            mpn_mul_n(t, ad + N*j, bd + (i - j)*N, N);
            z[2*N] += mpn_add_n(z, z, t, 2*N);
        }

        my_mpn_tdiv_qr(R.tempNp1, xd + N*i, z, 2*N, R.modulus_limbs, N);
    }

    x.normalize(xn, N);
}

void nmod_poly_divexact(nmod_ring& R,
    nmod_poly& q,
    const nmod_poly& a,
    const nmod_poly& b)
{
    std::cout << "nmod_poly_divexact not implemented" << std::endl;
    std::abort();
}

void nmod_poly_divrem(nmod_ring& R,
    nmod_poly& q,
    nmod_poly& r,
    const nmod_poly& a,
    const nmod_poly& b)
{
    ulimb N = R.stride;
    ulimb an = a.length();
    ulimb bn = b.length();
    const ulimb* bd = b.data();

    if (UNLIKELY(an < bn))
    {
        nmod_poly_set(R, r, a);
        q.set_length(0);
        return;
    }

    ulimb* blcinv = R.tempNb;
    nmod_inv(R, blcinv, bd + bn - 1);

    ulimb* qd = q.fit_alloc(an - (bn - 1), N);

    if (UNLIKELY(bn <= 1))
    {
        FLINT_ASSERT_ALWAYS(bn == 1);
        for (ulimb i = 0; i < an; i++)
            nmod_mul(R, qd + i*N, a.data() + i*N, blcinv);
        r.set_length(0);
        return;
    }

    ulimb* rd = r.fit_alloc(an, N);
    ui_vec_set(rd, a.data(), an*N);

    ulimb* t = R.temp2Np1;
    for (ulimb i = an - 1; i >= bn - 1; i--)
    {
        nmod_mul(R, qd + (i - (bn - 1))*N, rd + (i)*N, blcinv);

        for (ulimb j = 0; j < bn - 1; j++)
        {
            nmod_mul(R, t, bd + j*N, qd + (i - (bn - 1))*N);
            nmod_sub(R, rd + (j+i-(bn-1))*N,
                        rd + (j+i-(bn-1))*N, t);
        }
    }

    r.normalize(bn - 1, N);
    q.normalize(an - (bn - 1), N);
}
