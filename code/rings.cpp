#if 0
//*********************** ZZ[x] *************************
class rfmpz_poly_t {
public:

	typedef xfmpz_poly_t elem_t;

	int is_zero(const xfmpz_poly_t & a)
	{
	    return fmpz_poly_is_zero(a.data);
	}

	void zero(xfmpz_poly_t & a)
	{
	    fmpz_poly_zero(a.data);
	}

	void swap(xfmpz_poly_t & a, xfmpz_poly_t & b)
	{
	    fmpz_poly_swap(a.data, b.data);
	}

	void set(xfmpz_poly_t & a, const fmpz_t b)
	{
	    fmpz_poly_set_fmpz(a.data, b);
	}

	void set(xfmpz_poly_t & a, const xfmpz_poly_t & b)
	{
	    fmpz_poly_set(a.data, b.data);
	}

	void neg(xfmpz_poly_t & a, const xfmpz_poly_t & b)
	{
	    fmpz_poly_neg(a.data, b.data);
	}

	void add(xfmpz_poly_t & a, const xfmpz_poly_t & b, const xfmpz_poly_t & c)
	{
	    fmpz_poly_add(a.data, b.data, c.data);
	}

	void sub(xfmpz_poly_t & a, const xfmpz_poly_t & b, const xfmpz_poly_t & c)
	{
	    fmpz_poly_sub(a.data, b.data, c.data);
	}

	void mul(xfmpz_poly_t & a, const xfmpz_poly_t & b, const xfmpz_poly_t & c)
	{
	    fmpz_poly_mul(a.data, b.data, c.data);
	}

	void pow_ui(xfmpz_poly_t & a, const xfmpz_poly_t & b, ulong c)
	{
	    fmpz_poly_pow(a.data, b.data, c);
	}

	void divrem(xfmpz_poly_t & a, xfmpz_poly_t & r, const xfmpz_poly_t & b, const xfmpz_poly_t & c)
	{
	    fmpz_poly_divrem(a.data, r.data, b.data, c.data);
	}

	void divexact(xfmpz_poly_t & a, const xfmpz_poly_t & b, const xfmpz_poly_t & c)
	{
		#ifndef NDEBUG
		    xfmpz_poly_t r;
		    fmpz_poly_divrem(a.data, r.data, b.data, c.data);
		    assert(r.data->length == 0);
		#else
		    fmpz_poly_div(a.data, b.data, c.data);
		#endif
	}
};

//************************** ZZ[[x]] *************************
class xfmpz_poly_series_t {
public:
    xfmpz_poly_t terms;
    slong start;
    slong prec;

    bool is_canonical() const
    {
        if (prec < 0 || terms.data->length > prec)
            return false;

        if (terms.data->length <= 0)
            return true;

        if (fmpz_is_zero(terms.data->coeffs + 0))
            return false;

        return true;
    };

    slong absolute_prec() const
    {
        return start + prec;
    }

    void normalize()
    {
        assert(terms.data->length <= prec);
        while (terms.data->length > 0 && fmpz_is_zero(terms.data->coeffs + 0))
        {
            fmpz_poly_shift_right(terms.data, terms.data, 1);
            prec -= 1;
            start += 1;
        }
    };

    std::string tostring() const
    {
        std::string s = "#^" + stdstring_to_string(start) + "*(";
        s += terms.tostring() + " + O(#^" + stdstring_to_string(prec) + "))";
        return s;
    }
};

class rfmpz_poly_series_t {
public:

	slong prec;

	rfmpz_poly_series_t(slong prec_) : prec(prec_) {}

    void set_prec(slong prec_)
    {
        prec = prec_;
    }

	typedef xfmpz_poly_series_t elem_t;

	int is_zero(const xfmpz_poly_series_t & a)
	{
        assert(a.is_canonical());

        if (!fmpz_poly_is_zero(a.terms.data))
            return false;

        // TODO: exact representation so we throw much less often
        throw exception_arithmetic(17);

        return true;
	}

	void zero(xfmpz_poly_series_t & a)
	{
        a.start = 0;
        a.prec = 2*prec;
	    fmpz_poly_zero(a.terms.data);
	}

	void swap(xfmpz_poly_series_t & a, xfmpz_poly_series_t & b)
	{
	    fmpz_poly_swap(a.terms.data, b.terms.data);
        std::swap(a.start, b.start);
        std::swap(a.prec, b.prec);
	}

	void set(xfmpz_poly_series_t & a, const xfmpz_poly_series_t & b)
	{
	    fmpz_poly_set(a.terms.data, b.terms.data);
        a.start = b.start;
        a.prec = b.prec;
	}

	void set(xfmpz_poly_series_t & a, const fmpz_t b)
	{
	    fmpz_poly_set_fmpz(a.terms.data, b);
        a.start = 0;
        a.prec = prec;
	}

	void set(xfmpz_poly_series_t & a, const xfmpz_poly_t & b)
	{
	    fmpz_poly_set(a.terms.data, b.data);
        a.start = 0;
        while (a.terms.data->length > 0 && fmpz_is_zero(a.terms.data->coeffs + 0))
        {
            fmpz_poly_shift_right(a.terms.data, a.terms.data, 1);
            a.start += 1;
        }
        a.prec = prec;
        assert(a.is_canonical());
	}

	void neg(xfmpz_poly_series_t & a, const xfmpz_poly_series_t & b)
	{
	    fmpz_poly_neg(a.terms.data, b.terms.data);
        a.start = b.start;
        a.prec = b.prec;
	}

	void add(xfmpz_poly_series_t & a, const xfmpz_poly_series_t & b, const xfmpz_poly_series_t & c)
	{
        assert(b.is_canonical());
        assert(c.is_canonical());

        slong bs = b.start;
        slong cs = c.start;
        slong bp = b.prec;
        slong cp = c.prec;
        const fmpz_poly_struct * bt = b.terms.data;
        const fmpz_poly_struct * ct = c.terms.data;

        slong as = 0;
        if (bs > cs)
        {
            std::swap(bs, cs);
            std::swap(bp, cp);
            std::swap(bt, ct);
        }

        as = bs;
        cs -= bs;
        slong ap = std::min(bp, cp + cs);

        if (cs > 0)
        {
            fmpz_poly_t t;
            fmpz_poly_init(t);
            fmpz_poly_shift_left(t, ct, cs);
    	    fmpz_poly_add_series(a.terms.data, bt, t, ap);
            fmpz_poly_clear(t);
        }
        else
        {
    	    fmpz_poly_add_series(a.terms.data, bt, ct, ap);
        }
        a.start = as;
        a.prec = ap;
        a.normalize();
        assert(a.is_canonical());

	}

	void sub(xfmpz_poly_series_t & a, const xfmpz_poly_series_t & b, const xfmpz_poly_series_t & c)
	{
        assert(b.is_canonical());
        assert(c.is_canonical());
        bool is_neg = false;
        slong bs = b.start;
        slong cs = c.start;
        slong bp = b.prec;
        slong cp = c.prec;
        const fmpz_poly_struct * bt = b.terms.data;
        const fmpz_poly_struct * ct = c.terms.data;

        slong as = 0;
        if (bs > cs)
        {
            is_neg = true;
            std::swap(bs, cs);
            std::swap(bp, cp);
            std::swap(bt, ct);
        }

        as = bs;
        cs -= bs;
        slong ap = std::min(bp, cp + cs);

        if (cs > 0)
        {
            fmpz_poly_t t;
            fmpz_poly_init(t);
            fmpz_poly_shift_left(t, ct, cs);
            if (is_neg)
        	    fmpz_poly_sub_series(a.terms.data, t, bt, ap);
            else
        	    fmpz_poly_sub_series(a.terms.data, bt, t, ap);
            fmpz_poly_clear(t);
        }
        else
        {
            if (is_neg)
        	    fmpz_poly_sub_series(a.terms.data, ct, bt, ap);
            else
        	    fmpz_poly_sub_series(a.terms.data, bt, ct, ap);
        }

        a.start = as;
        a.prec = ap;
        a.normalize();
        assert(a.is_canonical());
	}

	void mul(xfmpz_poly_series_t & a, const xfmpz_poly_series_t & b, const xfmpz_poly_series_t & c)
	{
        assert(b.is_canonical());
        assert(c.is_canonical());
        slong ap = std::min(b.prec, c.prec);
        fmpz_poly_mullow(a.terms.data, b.terms.data, c.terms.data, ap);
        a.prec = ap;
        if (add_si_checked(a.start, b.start, c.start))
            throw exception_arithmetic(17);
        a.normalize();
        assert(a.is_canonical());
	}

	void pow_ui(xfmpz_poly_series_t & a, const xfmpz_poly_series_t & b, ulong c)
	{
        assert(b.is_canonical());
        fmpz_poly_pow_trunc(a.terms.data, b.terms.data, c, b.prec);        
        a.prec = b.prec;
        if (mul_si_checked(a.start, b.start, c))
            throw exception_arithmetic(17);
        a.normalize();
        assert(a.is_canonical());
	}

	void divrem(xfmpz_poly_series_t & a, xfmpz_poly_series_t & r, const xfmpz_poly_series_t & b, const xfmpz_poly_series_t & c)
	{
		assert(false);
	}

	void divexact(xfmpz_poly_series_t & a, const xfmpz_poly_series_t & b, const xfmpz_poly_series_t & c)
	{
        assert(b.is_canonical());
        assert(c.is_canonical());
        slong ap = std::min(b.prec, c.prec);
        fmpz_poly_div_series(a.terms.data, b.terms.data, c.terms.data, ap);
        a.prec = ap;
        if (sub_si_checked(a.start, b.start, c.start))
            throw exception_arithmetic(17);
        a.normalize();
        assert(a.is_canonical());
	}
};

// ************** R[y]: used for ZZ[x][y] and ZZ[[x]][y] ****************
template <class X> class sparse_poly {
public:
    std::vector<X> coeffs;
    std::vector<ulong> exps;
    slong length;

    void fit_length(slong alloc)
    {
        if (coeffs.size() < alloc)
            coeffs.resize(alloc);
        if (exps.size() < alloc)
            exps.resize(alloc);
    }

    void swap(sparse_poly<X> & a)
    {
        std::swap(coeffs, a.coeffs);
        std::swap(exps, a.exps);
        std::swap(length, a.length);
    }

    std::string tostring() {
        std::string s;
        bool first = true;
        for (slong i = 0; i < length; i++)
        {
            if (!first)
                s.append(" + ");
            s.push_back('(');
            s.append(coeffs[i].tostring());
            s.append(")*x^");
            s += stdstring_to_string(exps[i]);
            first = false;
        }
        if (first)
            s = "0";
        return s;
    }

};


template <class R>
void set_univar(
    sparse_poly<typename R::elem_t> & z,
    const xfmpz_poly_t & a,
    R & r)
{
    fmpz * a_coeffs = a.data->coeffs;
    slong i = a.data->length;
    z.fit_length(i);
    z.length = 0;
    for (i--; i >= 0; i--)
    {
        if (fmpz_is_zero(a_coeffs + i))
            continue;
        r.set(z.coeffs[z.length], a_coeffs + i);
        z.exps[z.length] = i;
        z.length++; // TODO: proper length
    }
}

template <class R>
void set_univar_shift(
    sparse_poly<typename R::elem_t> & z,
    xfmpz_poly_t a, // clobbered
    slong dir,
    R & r)
{
    z.fit_length(a.data->length);
    z.length = a.data->length;
    if (z.length <= 0)
        return;
    
    slong j = 0;
    slong i = z.length - 1; // TODO: proper length
    r.set(z.coeffs[i], a);
    z.exps[i] = j;
    for (i--; i >= 0; i--)
    {
        z.exps[i] = ++j;
        fmpz_poly_derivative(a.data, a.data);
        fmpz_poly_scalar_divexact_si(a.data, a.data, dir*j);
        r.set(z.coeffs[i], a);
    }
}

void set_univar_scale(sparse_poly<xfmpz_poly_t> & t, const xfmpz_poly_t a, slong dir)
{
    t.fit_length(a.data->length);
    fmpz * a_coeffs = a.data->coeffs;
    slong adeg = a.data->length - 1;
    t.length = 0;
    for (slong j = adeg; j >= 0; j--)
    {
        slong i = dir > 0 ? j : adeg - j;
        if (fmpz_is_zero(a_coeffs + i))
            continue;
        fmpz_poly_zero(t.coeffs[t.length].data);
        fmpz_poly_set_coeff_fmpz(t.coeffs[t.length].data, i, a_coeffs + i);
        t.exps[t.length] = j;
        t.length++;
    }
}

void set(sparse_poly<xfmpz_poly_t> & t, const sparse_poly<xfmpz_poly_t> & a)
{
    t.fit_length(a.length);
    for (slong i = a.length - 1; i >= 0; i--)
    {
        t.exps[i] = a.exps[i];
        fmpz_poly_set(t.coeffs[i].data, a.coeffs[i].data);
    }
    t.length = a.length;
}


/*
    A = prem(A, -B)
    C is used for working space
*/
template <class R>
void prem(
	sparse_poly<typename R::elem_t> & A,
	sparse_poly<typename R::elem_t> & B,
	sparse_poly<typename R::elem_t> & C,
	R & r)
{
    slong a_len, b_len, c_len;
    slong a_deg, b_deg;
    ulong * a_exp, * b_exp, * c_exp;
    typename R::elem_t * a_coeff, * b_coeff, * c_coeff;
    slong i, j, delta, delta_org;
    typename R::elem_t u, v;

//std::cout << "prem called" << std::endl;
//std::cout << "A: " << A.tostring() << std::endl;
//std::cout << "B: " << B.tostring() << std::endl;

    assert(A.length > 0);
    assert(B.length > 0);
    assert(A.exps[0] >= B.exps[0]);

    delta_org = A.exps[0] - B.exps[0] + 1;

    A.fit_length(A.length + B.exps[0]);
    C.fit_length(A.length + B.exps[0]);

    b_len = B.length;
    b_deg = B.exps[0];
    b_exp = B.exps.data();
    b_coeff = B.coeffs.data();

looper:

    a_len = A.length;
    a_deg = A.exps[0];
    a_exp = A.exps.data();
    a_coeff = A.coeffs.data();

    c_exp = C.exps.data();
    c_coeff = C.coeffs.data();

    delta = a_deg - b_deg;

    if (a_len == 0 || delta < 0)
        goto done;

    c_len = 0;
    i = 1;
    j = 1;
    while (i < a_len || j < b_len)
    {
        if (i < a_len && (j >= b_len || a_exp[i] > b_exp[j] + delta))
        {
            r.mul(c_coeff[c_len], a_coeff[i], b_coeff[0]);
            r.neg(c_coeff[c_len], c_coeff[c_len]);
            c_exp[c_len++] = a_exp[i++];
        }
        else if (j < b_len && (i >= a_len || b_exp[j] + delta > a_exp[i]))
        {
            r.mul(c_coeff[c_len], a_coeff[0], b_coeff[j]);
            c_exp[c_len++] = b_exp[j++] + delta;
        }
        else
        {
            assert(i < a_len && j < b_len && a_exp[i] == b_exp[j] + delta);
            r.mul(u, a_coeff[i], b_coeff[0]);
            r.mul(v, a_coeff[0], b_coeff[j]);
            r.sub(c_coeff[c_len], v, u);
            c_exp[c_len] = a_exp[i];
            c_len += !r.is_zero(c_coeff[c_len]);
            i++;
            j++;
        }
    }

    C.length = c_len;
    A.swap(C);

    delta_org--;

    goto looper;

done:

    if (delta_org != 0)
    {
        assert(delta_org > 0);
        r.neg(v, b_coeff[0]);
        r.pow_ui(u, v, delta_org);
        for (i = 0; i < A.length; i++)
        {
            r.mul(v, A.coeffs[i], u);
            r.swap(A.coeffs[i], v);
        }
    }

//std::cout << "prem returning" << std::endl;
//std::cout << "A: " << A.tostring() << std::endl;
//std::cout << "B: " << B.tostring() << std::endl;

}


template <class R>
void pgcd(
	sparse_poly<typename R::elem_t> & F,
	sparse_poly<typename R::elem_t> & B,
	sparse_poly<typename R::elem_t> & A,
	R & r)
{
//std::cout << "pgcd called" << std::endl;

    slong i, d, e;
    typename R::elem_t u, v, w, s;
    sparse_poly<typename R::elem_t> C, D;
    sparse_poly<typename R::elem_t> * last;

    assert(B.length > 0);
    assert(A.length > 0);
    assert(B.exps[0] >= A.exps[0]);
    assert(A.exps[0] >= 1);

    i = std::max(B.exps[0], A.exps[0]);
    A.fit_length(i + 1);
    B.fit_length(i + 1);
    C.fit_length(i + 1);
    D.fit_length(i + 1);

    last = &A;

    r.pow_ui(s, A.coeffs[0], B.exps[0] - A.exps[0]);

    prem<R>(B, A, D, r);

looper:

    d = A.exps[0];
    e = B.exps[0];
    if (B.length < 1)
        goto done;

    if (d - e > 1)
    {
        r.pow_ui(u, B.coeffs[0], d - e - 1); 
        r.pow_ui(v, s, d - e - 1);

        for (i = 0; i < B.length; i++)
        {
            r.mul(w, u, B.coeffs[i]);
            r.divexact(C.coeffs[i], w, v);
            C.exps[i] = B.exps[i];
        }
        C.length = B.length;
        r.mul(w, s, A.coeffs[0]);
        r.mul(u, v, w);
    }
    else
    {
        for (i = 0; i < B.length; i++)
        {
            set(C.coeffs[i], B.coeffs[i]);
            C.exps[i] = B.exps[i];
        }
        C.length = B.length;
        r.mul(u, s, A.coeffs[0]);
    }

    last = &C;
    if (e == 0)
    {
        goto done;
    }

    prem<R>(A, B, D);
    for (i = 0; i < A.length; i++)
    {
        r.divexact(B.coeffs[i], A.coeffs[i], u);
        B.exps[i] = A.exps[i];
    }
    B.length = A.length;

    A.swap(C);
    r.set(s, A.coeffs[0]);

    last = &A;

    goto looper;

done:

    F.swap(*last);

    return;
}


template <class R>
void pgcd_ducos(
	sparse_poly<typename R::elem_t> & F,
    sparse_poly<typename R::elem_t> & B,
	sparse_poly<typename R::elem_t> & A,
	R & r)
{
    ulong exp;
    slong i, j, k, d, e;
    slong alpha, n, J, aJ, ae;
    slong a_len, b_len, c_len, d_len, h_len, t_len;
    ulong * a_exp, * b_exp, * c_exp, * d_exp, * h_exp, * t_exp;
    typename R::elem_t * a_coeff, * b_coeff, * c_coeff, * d_coeff, * h_coeff, * t_coeff;
    int iexists, jexists, kexists;
    typename R::elem_t u, v, w, s;
    sparse_poly<typename R::elem_t> C, D, H, T;
    sparse_poly<typename R::elem_t> * last;

    assert(B.length > 0);
    assert(A.length > 0);
    assert(B.exps[0] >= A.exps[0]);
    assert(A.exps[0] >= 1);

    i = std::max(B.exps[0], A.exps[0]);
    A.fit_length(1 + i);
    B.fit_length(1 + i);
    C.fit_length(1 + i);
    D.fit_length(1 + i);
    H.fit_length(1 + i);
    T.fit_length(1 + i);

    last = &A;

    r.pow_ui(s, A.coeffs[0], B.exps[0] - A.exps[0]);

    prem<R>(B, A, D, r);

looper:

    d = A.exps[0];
    e = B.exps[0];
    if (B.length <= 0)
        goto done;

    last = &B;

    if (d - e == 1)
    {
        a_len = A.length;
        a_exp = A.exps.data();
        a_coeff = A.coeffs.data();

        b_len = B.length;
        b_exp = B.exps.data();
        b_coeff = B.coeffs.data();

        d_len = D.length;
        d_exp = D.exps.data();
        d_coeff = D.coeffs.data();

        if (e == 0)
            goto done;

        /* D = (B[e]*A - A[e]*B)/A[d] */
        /*           i        j       */
        i = 1;
        j = 1;
        if (a_len > 1 && a_exp[1] == e)
            i++;
        else
            j = b_len;
        d_len = 0;
        while (i < a_len || j < b_len)
        {
            if (i < a_len && j < b_len && a_exp[i] == b_exp[j])
            {
                r.mul(u, a_coeff[i], b_coeff[0]);
                r.mul(v, a_coeff[1], b_coeff[j]);
                r.sub(w, u, v);
                r.divexact(d_coeff[d_len], w, a_coeff[0]);
                d_exp[d_len] = a_exp[i];
                d_len += !r.is_zero(d_coeff[d_len]);
                i++;
                j++;                
            }
            else if (i < a_len && (j >= b_len || a_exp[i] > b_exp[j]))
            {
                r.mul(u, a_coeff[i], b_coeff[0]);
                r.divexact(d_coeff[d_len], u, a_coeff[0]);
                d_exp[d_len++] = a_exp[i];
                i++;
            }
            else
            {
                assert(j < b_len && (i >= a_len || b_exp[j] > a_exp[i]));
                r.mul(v, a_coeff[1], b_coeff[j]);
                r.divexact(d_coeff[d_len], v, a_coeff[0]);
                r.neg(d_coeff[d_len], d_coeff[d_len]);
                d_exp[d_len++] = b_exp[j];
                j++;
            }
        }
        D.length = d_len;

        /* A = (B[e]*(D - B*x) + B[e-1]*B)/s */
        /*            i    j            k    */
        i = 0;
        if (b_len > 1 && b_exp[1] == e - 1) {
            j = 2;            
            k = 1;
        } else {
            j = 1;
            k = b_len;
        }
        a_len = 0;
        while (i < d_len || j < b_len || k < b_len)
        {
            exp = 0;
            if (i < d_len)
                exp = std::max(exp, d_exp[i]);
            if (j < b_len)
                exp = std::max(exp, b_exp[j] + 1);
            if (k < b_len)
                exp = std::max(exp, b_exp[k]);

            a_exp[a_len] = exp;

            iexists = (i < d_len) && (exp == d_exp[i]);
            jexists = (j < b_len) && (exp == b_exp[j] + 1);
            kexists = (k < b_len) && (exp == b_exp[k]);

            assert(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    r.sub(w, d_coeff[i], b_coeff[j]);
                    r.mul(u, b_coeff[0], w);
                }
                else
                {
                    r.mul(u, b_coeff[0], d_coeff[i]);
                }
                if (kexists)
                {
                    r.mul(v, b_coeff[1], b_coeff[k]);
                    r.add(w, u, v);
                    r.divexact(a_coeff[a_len], w, s);
                }
                else
                {
                    r.divexact(a_coeff[a_len], u, s);
                }
            }
            else if (kexists)
            {
                r.mul(u, b_coeff[1], b_coeff[k]);
                if (jexists)
                {
                    r.mul(v, b_coeff[0], b_coeff[j]);
                    r.sub(w, u, v);
                    r.divexact(a_coeff[a_len], w, s);
                }
                else
                {
                    r.divexact(a_coeff[a_len], u, s);
                }
            }
            else
            {
                r.mul(u, b_coeff[0], b_coeff[j]);
                r.divexact(a_coeff[a_len], u, s);
                r.neg(a_coeff[a_len], a_coeff[a_len]);
            }

            a_len += !r.is_zero(a_coeff[a_len]);

            i += iexists;
            j += jexists;
            k += kexists;
        }
        A.length = a_len;

        /* A <-> B */
        A.swap(B);

        r.set(s, A.coeffs[0]);
        last = &A;
    }
    else
    {
        a_len = A.length;
        a_exp = A.exps.data();
        a_coeff = A.coeffs.data();
        b_len = B.length;
        b_exp = B.exps.data();
        b_coeff = B.coeffs.data();
        c_len = C.length;
        c_exp = C.exps.data();
        c_coeff = C.coeffs.data();
        d_len = D.length;
        d_exp = D.exps.data();
        d_coeff = D.coeffs.data();
        h_len = H.length;
        h_exp = H.exps.data();
        h_coeff = H.coeffs.data();
        t_len = T.length;
        t_exp = T.exps.data();
        t_coeff = T.coeffs.data();

        n = d - e - 1;
        assert(n > 0);

        alpha = 1;
        while (2*alpha <= n)
            alpha = 2*alpha;

        r.set(u, b_coeff[0]);
        n = n - alpha;
        while (alpha > 1)
        {
            alpha = alpha/2;
            r.mul(v, u, u);
            r.divexact(u, v, s);
            if (n >= alpha)
            {
                r.mul(v, u, b_coeff[0]);
                r.divexact(u, v, s);
                n = n - alpha;
            }
        }
        for (i = 0; i < b_len; i++)
        {
            r.mul(v, u, b_coeff[i]);
            r.divexact(c_coeff[i], v, s);
            c_exp[i] = b_exp[i];
        }
        c_len = b_len;
        C.length = c_len;

        last = &C;

        if (e == 0)
            goto done;

        /* H = C - C[e]*x^e */
        for (i = 1; i < c_len; i++)
        {
            r.set(h_coeff[i - 1], c_coeff[i]);
            h_exp[i - 1] = c_exp[i];
        }
        h_len = c_len - 1;
        H.length = h_len;

        /* D = C[e]*A - A[e]*H  (truncated to powers of x < e) */
        i = 0;
        j = h_len;
        ae = a_len;
        while (i < a_len && a_exp[i] >= e)
        {
            if (a_exp[i] == e)
            {
                j = 0;
                ae = i;
            }
            i++;
        }
        d_len = 0;
        while (i < a_len || j < h_len)
        {
            if (i < a_len && j < h_len && a_exp[i] == h_exp[j])
            {
                r.mul(u, a_coeff[i], c_coeff[0]);
                r.mul(v, a_coeff[ae], h_coeff[j]);
                r.sub(d_coeff[d_len], u, v);
                d_exp[d_len] = a_exp[i];
                d_len += !r.is_zero(d_coeff[d_len]);
                i++;
                j++;
            }
            else if (i < a_len && (j >= h_len || a_exp[i] > h_exp[j]))
            {
                r.mul(d_coeff[d_len], a_coeff[i], c_coeff[0]);
                d_exp[d_len++] = a_exp[i];
                i++;
            }
            else
            {
                assert(j < h_len && (i >= a_len || h_exp[j] > a_exp[i]));
                r.mul(d_coeff[d_len], a_coeff[ae], h_coeff[j]);
                r.neg(d_coeff[d_len], d_coeff[d_len]);
                d_exp[d_len++] = h_exp[j];
                j++;
            }
        }
        D.length = d_len;

        for (J = e + 1; J < d; J++)
        {
            if (h_len == 0)
                break;

            /* H = H*x - H[e-1]*B/B[e] */
            if (h_exp[0] == e - 1)
            {
                i = 1;
                j = 1;
                t_len = 0;
                while (i < h_len || j < b_len)
                {
                    if (i < h_len && j < b_len && h_exp[i] + 1 == b_exp[j])
                    {
                        r.mul(u, h_coeff[0], b_coeff[j]);
                        r.divexact(v, u, b_coeff[0]);
                        r.sub(t_coeff[t_len], h_coeff[i], v);
                        t_exp[t_len] = b_exp[j];
                        t_len += !r.is_zero(t_coeff[t_len]);
                        i++;
                        j++;
                    }
                    else if (i < h_len && (j >= b_len || h_exp[i] + 1 > b_exp[j]))
                    {
                        r.swap(t_coeff[t_len], h_coeff[i]);
                        t_exp[t_len++] = h_exp[i] + 1;
                        i++;
                    }
                    else
                    {
                        assert(j < b_len && (i >= h_len || b_exp[j] > h_exp[i] + 1));
                        r.mul(u, h_coeff[0], b_coeff[j]);
                        r.divexact(t_coeff[t_len], u, b_coeff[0]);
                        r.neg(t_coeff[t_len], t_coeff[t_len]);
                        t_exp[t_len++] = b_exp[j];
                        j++;
                    }
                }
                T.length = t_len;

                H.swap(T);
                h_len = H.length;
                h_exp = H.exps.data();
                h_coeff = H.coeffs.data();
                t_len = T.length;
                t_exp = T.exps.data();
                t_coeff = T.coeffs.data();
            }
            else
            {
                assert(h_exp[0] < e - 1);
                for (i = 0; i < h_len; i++)
                    h_exp[i]++;
            }

            /* find coefficient of x^J in A */
            aJ = 0;
            while (aJ < a_len && a_exp[aJ] != J)
                aJ++;
            if (aJ >= a_len)
                continue;

            /* D = D - A[J]*H */
            i = 0;
            j = 0;
            t_len = 0;
            while (i < d_len || j < h_len)
            {
                if (i < d_len && j < h_len && d_exp[i] == h_exp[j])
                {
                    r.mul(u, h_coeff[j], a_coeff[aJ]);
                    r.sub(t_coeff[t_len], d_coeff[i], u);
                    t_exp[t_len] = d_exp[i];
                    t_len += !r.is_zero(t_coeff[t_len]);
                    i++;
                    j++;                
                }
                else if (i < d_len && (j >= h_len || d_exp[i] > h_exp[j]))
                {
                    r.swap(t_coeff[t_len], d_coeff[i]);
                    t_exp[t_len++] = d_exp[i];
                    i++;
                }
                else
                {
                    assert(j < h_len && (i >= d_len || h_exp[j] > d_exp[i]));
                    r.mul(t_coeff[t_len], h_coeff[j], a_coeff[aJ]);
                    r.neg(t_coeff[t_len], t_coeff[t_len]);
                    t_exp[t_len++] = h_exp[j];
                    j++;
                }
            }
            T.length = t_len;

            D.swap(T);
            d_len = D.length;
            d_exp = D.exps.data();
            d_coeff = D.coeffs.data();
            t_len = T.length;
            t_exp = T.exps.data();
            t_coeff = T.coeffs.data();
        }

        /* B = (-1)^(d-e+1) * (B[e]*(D/A[d] - H*x) +  H[e-1]*B)/s */
        i = 0;
        if (h_len > 0 && h_exp[0] == e - 1) {
            j = 1;
            k = 1;
        } else {
            j = 0;
            k = b_len;
        }
        t_len = 0;
        while (i < d_len || j < h_len || k < b_len)
        {
            exp = 0;
            if (i < d_len)
                exp = std::max(exp, d_exp[i]);
            if (j < h_len)
                exp = std::max(exp, h_exp[j] + 1);
            if (k < b_len)
                exp = std::max(exp, b_exp[k]);

            t_exp[t_len] = exp;

            iexists = (i < d_len && exp == d_exp[i]);
            jexists = (j < h_len && exp == h_exp[j] + 1);
            kexists = (k < b_len && exp == b_exp[k]);

            assert(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    r.divexact(u, d_coeff[i], a_coeff[0]);
                    r.sub(w, u, h_coeff[j]);
                    r.mul(u, b_coeff[0], w);
                }
                else
                {
                    r.divexact(u, d_coeff[i], a_coeff[0]);
                    r.mul(u, b_coeff[0], u);
                }
                if (kexists)
                {
                    r.mul(v, h_coeff[0], b_coeff[k]);
                    r.add(w, u, v);
                    r.divexact(t_coeff[t_len], w, s);
                }
                else
                {
                    r.divexact(t_coeff[t_len], u, s);
                }
            }
            else if (kexists)
            {
                r.mul(u, h_coeff[0], b_coeff[k]);
                if (jexists)
                {
                    r.mul(v, b_coeff[0], h_coeff[j]);
                    r.sub(w, u, v);
                    r.divexact(t_coeff[t_len], w, s);
                }
                else
                {
                    r.divexact(t_coeff[t_len], u, s);
                }
            }
            else
            {
                r.mul(u, b_coeff[0], h_coeff[j]);
                r.divexact(t_coeff[t_len], u, s);
                r.neg(t_coeff[t_len], t_coeff[t_len]);
            }

            if (((d - e) & 1) == 0)
                r.neg(t_coeff[t_len], t_coeff[t_len]);

            t_len += !r.is_zero(t_coeff[t_len]);

            i += iexists;
            j += jexists;
            k += kexists;
        }
        T.length = t_len;

        /* B <-> T */
        B.swap(T);
        b_len = T.length;
        b_exp = T.exps.data();
        b_coeff = T.coeffs.data();
        t_len = T.length;
        t_exp = T.exps.data();
        t_coeff = T.coeffs.data();

        /* A <-> C */
        A.swap(C);
        a_len = A.length;
        a_exp = A.exps.data();
        a_coeff = A.coeffs.data();
        c_len = C.length;
        c_exp = C.exps.data();
        c_coeff = C.coeffs.data();

        r.set(s, A.coeffs[0]);

        last = &A;
    }

    goto looper;

done:

    F.swap(*last);
}




template <class R>
void resultant(
	typename R::elem_t & z,
	sparse_poly<typename R::elem_t> & a,
	sparse_poly<typename R::elem_t> & b,
	R & r)
{
//std::cout << "resultant called" << std::endl;
//std::cout << "a: " << a.tostring() << std::endl;
//std::cout << "b: " << b.tostring() << std::endl;

    sparse_poly<typename R::elem_t> t;
    bool neg = false;

    if (a.exps[0] >= b.exps[0])
	{
        pgcd_ducos<R>(t, a, b, r);
	}
    else
	{
        neg = (a.exps[0] & b.exps[0] & 1);
        pgcd_ducos<R>(t, b, a, r);
	}

    if (t.length != 1 || t.exps[0] != 0)
    {
        r.zero(z);
    }
    else
    {
        r.swap(z, t.coeffs[0]);
        if (neg)
            r.neg(z, z);
    }

//std::cout << "resultant returning" << std::endl;
//std::cout << "z: " << z.tostring() << std::endl;
}
#endif