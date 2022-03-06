rfq_nmod::rfq_nmod(const ulong* m, slong d)
{
    t = nmod_neg(nmod_inv(m[d]));

    for (slong i = 0 ; i < d; i++)
        M[i*(d - 1) + 0] = nmod_mul(m[i], t);

    for (j = 1; j < d - 1; j++)
    {
        slong i = 0;
        M[i*(d-1) + j] = nmod_mul(M[i*(d-1) + 0], M[(d-1)*(d-1) + j-1]);
        for (i++; i < d; i++)
            M[i*(d-1) + j] = nmod_addmul(M[(i-1)*(d-1) + (j-1)],
                                         M[i*(d-1) + 0], M[(d-1)*(d-1) + j-1]);
    }
}


void rfq_nmod::mul3_from(const ulong* b, const ulong* c)
{
    UMUL_PPMM(p1, p0, b[0], c[d-1]); p2 = 0;
    for (slong i = 1; i < d; i++)
        UMADD3(p2, p1, p0, b[i], c[d-1-i]);
    P[3*(d-1)+0] = p0; P[3*(d-1)+1] = p1; P[3*(d-1)+2] = p2;

    for (slong j = d - 2; j >= 0; j--)
    {
        // handle P[j] and P[2*d-2-j]
        UMUL_PPMM(p1, p0, b[j], c[0]); p2 = 0;
        UMUL_PPMM(q1, q0, b[d-1], c[d-1-j]); q2 = 0;
        for (slong i = 1; i <= j; i++)
        {
            UMADD3(p2, p1, p0, b[j-i], c[i]);
            UMADD3(q2, q1, q0, b[d-1-i], c[d-1+j+i]);
        }
        P[3*(j)+0] = p0; P[3*(j)+1] = p1; P[3*(j)+2] = p2;
        P[3*(2*d-2-j)+0] = q0; P[3*(2*d-2-j)+1] = q1; P[3*(2*d-2-j)+2] = q2;
    }
}

void rfq_nmod::madd3_from(const ulong* b, const ulong* c)
{
    p0 = P[3*(d-1)+0]; p1 = P[3*(d-1)+1]; p2 = P[3*(d-1)+2];
    for (slong i = 0; i < d; i++)
        MADD3(p2, p1, p0, b[i], c[d-1-i]);
    P[3*(d-1)+0] = p0; P[3*(d-1)+1] = p1; P[3*(d-1)+2] = p2;

    for (slong j = d - 2; j >= 0; j--)
    {
        p0 = P[3*(j)+0] = p0; P[3*(j)+1] = p1; P[3*(j)+2] = p2;
        q0 = P[3*(2*d-2-j)+0]; q1 = P[3*(2*d-2-j)+1]; q2 = P[3*(2*d-2-j)+2];
        for (slong i = 0; i <= j; i++)
        {
            UMADD3(p2, p1, p0, b[j-i], c[i]);
            UMADD3(q2, q1, q0, b[d-1-i], c[d-1+j+i]);
        }
        P[3*(j)+0] = p0; P[3*(j)+1] = p1; P[3*(j)+2] = p2;
        P[3*(2*d-2-j)+0] = q0; P[3*(2*d-2-j)+1] = q1; P[3*(2*d-2-j)+2] = q2;
    }
}

void rfq_nmod::reduce3_to(ulong* a)
{
    for (j = 0; j < d - 1; j++)
        t[j] = nmod_reduce_3(P[3*(j+d)+2], P[3*(j+d)+1], P[3*(j+d)+0]);

    for (i = 0; i < d; i++)
    {
        ulong p2 = P[3*i+2], p1 = P[3*i+1], p0 = P[3*i+0];
        for (j = 0; j < d - 1; j++)
            MADD3(p2, p1, p0, t[j], M[i*(d - 1) + j]);

        a[i] = nmod_reduce_3(p2, p1, p0);
    }
}

template <class T>
void _upoly_make_monic(
    T& R,
    T_elem_t* A, slong Alen)
{
    slong idx = R.elem_stride();

    slong i = Alen - 1;
    if (i < 0)
        return;

    T_elem_t* u = R.leaf_temp(0);
    R.inv(u, A + idx*i);
    R.one(A + idx*i);
    for (slong i = 0; i < Alen - 1; i++)
        R.mul(A + idx*i, u);
}

// reduce A mod B and return new length(A)
template <class T>
slong _upoly_rem_inpl(
    T& R,
    T_elem_t* A, slong Alen,
    T_elem_t* B, slong Blen)
{
    slong idx = R.elem_stride();
    T_elem_t* u = R.leaf_temp(0);
    T_elem_t* q0 = R.leaf_temp(1);
    T_elem_t* q1 = R.leaf_temp(2);

    assert(Blen >= 2);

    while (Alen >= Blen)
    {
        if (Alen == Blen)
        {
            R.div(u, A + idx*(Alen - 1), B + idx*(Blen - 1));
            R.neg(u);

            for (slong i = 0; i < Blen - 1; i++)
                R.madd(A + idx*i, u, B + idx*i);

            Alen -= 1;
        }
        else
        {   // quo(a[n+1]*x^(n+1) + a[n]*x^n + ..., b[n]*x^n + ...) = Q[1]*x + Q[0]
            // Q[1] = a[n+1]/b[n]
            // Q[0] = a[n]/b[n] - Q[1]*b[n-1]/b[n]

            R.inv(u, B + idx*(Blen - 1));
            R.mul(q1, A + idx*(Alen - 1), u);
            R.mul(q0, q1, B + idx*(Blen - 2));
            R.sub(q0, A + idx*(Alen - 2));
            R.mul(q0, u);
            R.neg(q1);

            R.madd(A + idx*(Alen - Blen - 1), q0, B + idx*(0));
            for (slong i = 0; i < Blen - 2; i++)
            {
                R.dot_start(A + idx*(Alen - Blen + i));
                R.dot_madd(q1, B + idx*(i + 0));
                R.dot_madd(q0, B + idx*(i + 1));
                R.dot_finish(A + idx*(Adeg - Bdeg + i));
            }

            Alen -= 2;
        }

        while (Alen > 0 && R.is_zero(A + idx*(Alen - 1)))
            Alen--;
    }

    return Alen;
}


// output Qlen = max(Alen - Blen + 1, 0)
template <class T>
slong _upoly_divrem(
    T& F,
    T_elem_t* Q,
    T_elem_t* R,
    const T_elem_t* A, slong Alen,
    const T_elem_t* B, slong Blen)
{
    
}

template <class T>
slong _upoly_gcd_euclidean_inpl<T>(
    T& R,
    T_elem_t* A, slong Alen,
    T_elem_t* B, slong Blen)
{
    while (Alen > 1 && Blen > 1)
    {
        if (Alen > Blen)
            Alen = _upoly_rem_inpl<T>(R, A, Alen, B, Blen);
        else
            Blen = _upoly_rem_inpl<T>(R, B, Blen, A, Alen);
    }

    if (Alen < 1)
    {
        if (Blen < 1)
            return 0;

        _upoly_make_monic<T>(R, B, Blen);
        return -Blen - 1;
    }

    if (Blen < 1)
    {
        _upoly_make_monic<T>(R, A, Alen);
        return Alen;
    }

    if (Blen <= 1)
    {
        R.one(B[0], d);
        return -1 - 1;        
    }

    assert(Alen == 1);

    R.one(A[0]);
    return 1; 
}


