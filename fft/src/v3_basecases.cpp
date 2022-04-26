/*
    z0 = z0 + w*z1
    z1   z0 - w*z1
*/
FORCE_INLINE inline void radix_2_moth(
    std::complex<double>& z0,
    std::complex<double>& z1,
    std::complex<double> w)
{
    double t0, t1, t2, t3;
    double ax = w.real();
    double ay = w.imag();
    double z0x = z0.real();
    double z0y = z0.imag();
    double z1x = z1.real();
    double z1y = z1.imag();
    t0 = fmadd(ax, z1x, z0x);
    t1 = fmadd(ay, z1x, z0y);
    t0 = fnmadd(ay, z1y, t0);
    t1 = fmadd(ax, z1y, t1);
    t2 = fmsub(2.0, z0x, t0);
    t3 = fmsub(2.0, z0y, t1);
    z0.real(t0);
    z0.imag(t1);
    z1.real(t2);
    z1.imag(t3);
}

/*
z0   z0 + w2*z2 + w*(z1 + w2*z3)
z1 = z0 + w2*z2 - w*(z1 + w2*z3)
z2   z0 - w2*z2 + I*w*(z1 - w2*z3)
z3   z0 - w2*z2 - I*w*(z1 - w2*z3)
*/
FORCE_INLINE inline void radix_4_moth(
    std::complex<double>& z0,
    std::complex<double>& z1,
    std::complex<double>& z2,
    std::complex<double>& z3,
    std::complex<double> w,
    std::complex<double> w2)
{
    double z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    double ax, ay, bx, by, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    ax = w2.real();
    ay = w2.imag();
    bx = w.real();
    by = w.imag();
    z0x = z0.real();
    z0y = z0.imag();
    z2x = z2.real();
    z2y = z2.imag();
    z1x = z1.real();
    z1y = z1.imag();
    z3x = z3.real();
    z3y = z3.imag();
    t0 = fmadd(ax, z2x, z0x);
    t1 = fmadd(ay, z2x, z0y);
    t0 = fnmadd(ay, z2y, t0);
    t1 = fmadd(ax, z2y, t1);
    t2 = fmsub(2.0, z0x, t0);
    t3 = fmsub(2.0, z0y, t1);
    t4 = fmadd(ax, z3x, z1x);
    t5 = fmadd(ay, z3x, z1y);
    t4 = fnmadd(ay, z3y, t4);
    t5 = fmadd(ax, z3y, t5);
    t6 = fmsub(2.0, z1x, t4);
    t7 = fmsub(2.0, z1y, t5);
    t8 = fmadd(bx, t4, t0);
    t4 = fmadd(by, t4, t1);
    t9 = fnmadd(by, t6, t2);
    t6 = fmadd(bx, t6, t3);
    t8 = fnmadd(by, t5, t8);
    t5 = fmadd(bx, t5, t4);
    t9 = fnmadd(bx, t7, t9);
    t7 = fnmadd(by, t7, t6);
    t0 = fmsub(2.0, t0, t8);
    t1 = fmsub(2.0, t1, t5);
    t2 = fmsub(2.0, t2, t9);
    t3 = fmsub(2.0, t3, t7);
    z0.real(t8);
    z0.imag(t5);
    z1.real(t0);
    z1.imag(t1);
    z2.real(t9);
    z2.imag(t7);
    z3.real(t2);
    z3.imag(t3);
}

/*
z0   z0 + ~w2*z1 +   ~w*(z2 + ~w2*z3)
z1   z0 - ~w2*z1 - I*~w*(z2 - ~w2*z3)
z2 = z0 + ~w2*z1 -   ~w*(z2 + ~w2*z3)
z3   z0 - ~w2*z1 + I*~w*(z2 - ~w2*z3)
*/
inline void radix_4_mothi(
    std::complex<double>& z0,
    std::complex<double>& z1,
    std::complex<double>& z2,
    std::complex<double>& z3,
    std::complex<double> w,
    std::complex<double> w2)
{
    double z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    double ax, ay, bx, by, m0, m1, m2, m3, m4, m5, m6, m7, m8, m9;
    ax = w2.real();
    ay = w2.imag();
    bx = w.real();
    by = w.imag();
    z0x = z0.real();
    z0y = z0.imag();
    z1x = z2.real();
    z1y = z2.imag();
    z2x = z1.real();
    z2y = z1.imag();
    z3x = z3.real();
    z3y = z3.imag();

    m0 = fmadd(ax, z2x, z0x);
    m1 = fnmadd(ay, z2x, z0y);
    m0 = fmadd(ay, z2y, m0);
    m1 = fmadd(ax, z2y, m1);
    m2 = fmsub(2, z0x, m0);
    m3 = fmsub(2, z0y, m1);
    m4 = fmadd(ax, z3x, z1x);
    m5 = fnmadd(ay, z3x, z1y);
    m4 = fmadd(ay, z3y, m4);
    m5 = fmadd(ax, z3y, m5);
    m6 = fmsub(2, z1x, m4);
    m7 = fmsub(2, z1y, m5);
    m8 = fmadd(bx, m4, m0);
    m4 = fnmadd(by, m4, m1);
    m9 = fmadd(by, m6, m2);
    m6 = fmadd(bx, m6, m3);
    m8 = fmadd(by, m5, m8);
    m5 = fmadd(bx, m5, m4);
    m9 = fnmadd(bx, m7, m9);
    m7 = fmadd(by, m7, m6);
    m0 = fmsub(2, m0, m8);
    m1 = fmsub(2, m1, m5);
    m2 = fmsub(2, m2, m9);
    m3 = fmsub(2, m3, m7);
    z0.real(m8);
    z0.imag(m5);
    z2.real(m0);
    z2.imag(m1);
    z3.real(m9);
    z3.imag(m7);
    z1.real(m2);
    z1.imag(m3);
}

// radix_4_mothi with w = 1
inline void radix_4_mothi1(
    std::complex<double>& z0,
    std::complex<double>& z1,
    std::complex<double>& z2,
    std::complex<double>& z3)
{
    double z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    double s0x, s0y, s1x, s1y, s2x, s2y, s3x, s3y;
    z0x = z0.real();
    z0y = z0.imag();
    z1x = z1.real();
    z1y = z1.imag();
    z2x = z2.real();
    z2y = z2.imag();
    z3x = z3.real();
    z3y = z3.imag();
    s0x = z0x + z1x;
    s0y = z0y + z1y;
    s1x = z0x - z1x;
    s1y = z0y - z1y;
    s2x = z2x + z3x;
    s2y = z2y + z3y;
    s3x = z2x - z3x;
    s3y = z2y - z3y;
    z0.real(s0x + s2x);
    z0.imag(s0y + s2y);
    z2.real(s0x - s2x);
    z2.imag(s0y - s2y);
    z1.real(s1x + s3y);
    z1.imag(s1y - s3x);
    z3.real(s1x - s3y);
    z3.imag(s1y + s3x);
}


#define OPT_MOTH_2
/*
    z0 = z0 + w*z1
    z1   z0 - w*z1
*/
inline void radix_2_moth_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4> w)
{
    packed<double,4> t0, t1, t2, t3;
    packed<double,4> ax = w.real();
    packed<double,4> ay = w.imag();
    packed<double,4> z0x = z0.real();
    packed<double,4> z0y = z0.imag();
    packed<double,4> z1x = z1.real();
    packed<double,4> z1y = z1.imag();
#ifdef OPT_MOTH_2
    t0 = fnmadd(ay, z1y, fmadd(ax, z1x, z0x));
    t1 = fmadd(ax, z1y, fmadd(ay, z1x, z0y));
    z0.real(t0);
    z0.imag(t1);
    z1.real(fmsub(2, z0x, t0));
    z1.imag(fmsub(2, z0y, t1));
#else
    z0.real(fnmadd(ay, z1y,  fmadd(ax, z1x, z0x)));
    z0.imag( fmadd(ay, z1x,  fmadd(ax, z1y, z0y)));
    z1.real( fmadd(ay, z1y, fnmadd(ax, z1x, z0x)));
    z1.imag(fnmadd(ay, z1x, fnmadd(ax, z1y, z0y)));
#endif
}

/*
    z0 = z0 + I*w*z1
    z1   z0 - I*w*z1
*/
inline void radix_2_moth_i_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4> w)
{
    packed<double,4> ax = w.real();
    packed<double,4> ay = w.imag();
    packed<double,4> z0x = z0.real();
    packed<double,4> z0y = z0.imag();
    packed<double,4> z1x = z1.real();
    packed<double,4> z1y = z1.imag();
#ifdef OPT_MOTH_2
    packed<double,4> t0 = fnmadd(ay, z1x, fnmadd(ax, z1y, z0x));
    packed<double,4> t1 = fnmadd(ay, z1y,  fmadd(ax, z1x, z0y));
    z0.real(t0);
    z0.imag(t1);
    z1.real(fmsub(2, z0x, t0));
    z1.imag(fmsub(2, z0y, t1));
#else
    z0.real(fnmadd(ay, z1x, fnmadd(ax, z1y, z0x)));
    z0.imag(fnmadd(ay, z1y,  fmadd(ax, z1x, z0y)));
    z1.real( fmadd(ay, z1x,  fmadd(ax, z1y, z0x)));
    z1.imag( fmadd(ay, z1y, fnmadd(ax, z1x, z0y)));
#endif
}

/*
    z0 = z0 + ~w*z1
    z1   z0 - ~w*z1
*/
inline void radix_2_moth_inv_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4> w)
{
    packed<double,4> ax = w.real();
    packed<double,4> ay = w.imag();
    packed<double,4> z0x = z0.real();
    packed<double,4> z0y = z0.imag();
    packed<double,4> z1x = z1.real();
    packed<double,4> z1y = z1.imag();

#ifdef OPT_MOTH_2
    packed<double,4> t0 =  fmadd(ay, z1y,  fmadd(ax, z1x, z0x));
    packed<double,4> t1 = fnmadd(ay, z1x,  fmadd(ax, z1y, z0y));
    z0.real(t0);
    z0.imag(t1);
    z1.real(fmsub(2, z0x, t0));
    z1.imag(fmsub(2, z0y, t1));
#else
    z0.real( fmadd(ay, z1y,  fmadd(ax, z1x, z0x)));
    z0.imag(fnmadd(ay, z1x,  fmadd(ax, z1y, z0y)));
    z1.real(fnmadd(ay, z1y, fnmadd(ax, z1x, z0x)));
    z1.imag( fmadd(ay, z1x, fnmadd(ax, z1y, z0y)));
#endif
}

// radix_2_moth_block or radix_2_moth_inv_block with w = 1
inline void radix_2_moth_1_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1)
{
    packed<double,4> z0x = z0.real();
    packed<double,4> z0y = z0.imag();
    packed<double,4> z1x = z1.real();
    packed<double,4> z1y = z1.imag();
    z0.real(add(z0x, z1x));
    z0.imag(add(z0y, z1y));
    z1.real(sub(z0x, z1x));
    z1.imag(sub(z0y, z1y));
}


/*
    z0 = z0 - I*w^-1*z1
    z1   z0 + I*w^-1*z1
*/
inline void radix_2_moth_inv_i_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4> w)
{
    packed<double,4> ax = w.real();
    packed<double,4> ay = w.imag();
    packed<double,4> z0x = z0.real();
    packed<double,4> z0y = z0.imag();
    packed<double,4> z1x = z1.real();
    packed<double,4> z1y = z1.imag();
#ifdef OPT_MOTH_2
    packed<double,4> u0 = fnmadd(ay, z1x,  fmadd(ax, z1y, z0x));
    packed<double,4> v0 = fnmadd(ay, z1y, fnmadd(ax, z1x, z0y));
    z0.real(u0);
    z0.imag(v0);
    z1.real(fmsub(2.0, z0x, u0));
    z1.imag(fmsub(2.0, z0y, v0));
#else
    z0.real(fnmadd(ay, z1x,  fmadd(ax, z1y, z0x)));
    z0.imag(fnmadd(ay, z1y, fnmadd(ax, z1x, z0y)));
    z1.real( fmadd(ay, z1x, fnmadd(ax, z1y, z0x)));
    z1.imag( fmadd(ay, z1y,  fmadd(ax, z1x, z0y)));
#endif
}

/*
z0   z0 + w2*z2 +   w*(z1 + w2*z3)
z1 = z0 + w2*z2 -   w*(z1 + w2*z3)
z2   z0 - w2*z2 + I*w*(z1 - w2*z3)
z3   z0 - w2*z2 - I*w*(z1 - w2*z3)
*/
FORCE_INLINE inline void radix_4_moth_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4>& z2,
    complex_packed<double,4>& z3,
    complex_packed<double,4> w,
    complex_packed<double,4> w2)
{
#if 0
    complex_packed<double,4> s0 = z0, s1 = z1, s2 = z2, s3 = z3;
    radix_2_moth_block(s0, s2, w2);
    radix_2_moth_block(s1, s3, w2);
    radix_2_moth_block(s0, s1, w);
    radix_2_moth_i_block(s2, s3, w);
    z0 = s0;
    z1 = s1;
    z2 = s2;
    z3 = s3;
#else
    packed<double,4> z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    packed<double,4> ax, ay, bx, by, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    ax = w2.real();
    ay = w2.imag();
    bx = w.real();
    by = w.imag();
    z0x = z0.real();
    z0y = z0.imag();
    z2x = z2.real();
    z2y = z2.imag();
    z1x = z1.real();
    z1y = z1.imag();
    z3x = z3.real();
    z3y = z3.imag();
    t0 = fmadd(ax, z2x, z0x);
    t1 = fmadd(ay, z2x, z0y);
    t0 = fnmadd(ay, z2y, t0);
    t1 = fmadd(ax, z2y, t1);
    t2 = fmsub(2.0, z0x, t0);
    t3 = fmsub(2.0, z0y, t1);
    t4 = fmadd(ax, z3x, z1x);
    t5 = fmadd(ay, z3x, z1y);
    t4 = fnmadd(ay, z3y, t4);
    t5 = fmadd(ax, z3y, t5);
    t6 = fmsub(2.0, z1x, t4);
    t7 = fmsub(2.0, z1y, t5);
    t8 = fmadd(bx, t4, t0);
    t4 = fmadd(by, t4, t1);
    t9 = fnmadd(by, t6, t2);
    t6 = fmadd(bx, t6, t3);
    t8 = fnmadd(by, t5, t8);
    t5 = fmadd(bx, t5, t4);
    t9 = fnmadd(bx, t7, t9);
    t7 = fnmadd(by, t7, t6);
    t0 = fmsub(2.0, t0, t8);
    t1 = fmsub(2.0, t1, t5);
    t2 = fmsub(2.0, t2, t9);
    t3 = fmsub(2.0, t3, t7);
    z0.real(t8);
    z0.imag(t5);
    z1.real(t0);
    z1.imag(t1);
    z2.real(t9);
    z2.imag(t7);
    z3.real(t2);
    z3.imag(t3);
#endif
}


template<int itrunc>
FORCE_INLINE inline void radix_4_moth_trunc_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4>& z2,
    complex_packed<double,4>& z3,
    complex_packed<double,4> w,
    complex_packed<double,4> w2)
{
    packed<double,4> z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    packed<double,4> ax, ay, bx, by, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    ax = w2.real();
    ay = w2.imag();
    bx = w.real();
    by = w.imag();

    if (itrunc == 0)
    {
        z0.real(0.0);
        z0.imag(0.0);
        z1.real(0.0);
        z1.imag(0.0);
        z2.real(0.0);
        z2.imag(0.0);
        z3.real(0.0);
        z3.imag(0.0);
        return;
    }
    else if (itrunc == 1)
    {
        z0x = z0.real();
        z0y = z0.imag();
        z0.real(z0x);
        z0.imag(z0y);
        z1.real(z0x);
        z1.imag(z0y);
        z2.real(z0x);
        z2.imag(z0y);
        z3.real(z0x);
        z3.imag(z0y);
        return;
    }

    z0x = 0 < itrunc ? z0.real() : 0.0;
    z0y = 0 < itrunc ? z0.imag() : 0.0;
    z2x = 2 < itrunc ? z2.real() : 0.0;
    z2y = 2 < itrunc ? z2.imag() : 0.0;
    z1x = 1 < itrunc ? z1.real() : 0.0;
    z1y = 1 < itrunc ? z1.imag() : 0.0;
    z3x = 3 < itrunc ? z3.real() : 0.0;
    z3y = 3 < itrunc ? z3.imag() : 0.0;
    if (2 < itrunc)
    {
        t0 = fmadd(ax, z2x, z0x);
        t1 = fmadd(ay, z2x, z0y);
        t0 = fnmadd(ay, z2y, t0);
        t1 = fmadd(ax, z2y, t1);
        t2 = fmsub(2.0, z0x, t0);
        t3 = fmsub(2.0, z0y, t1);
    }
    else
    {
        t0 = t2 = z0x;
        t1 = t3 = z0y;
    }
    if (3 < itrunc)
    {
        t4 = fmadd(ax, z3x, z1x);
        t5 = fmadd(ay, z3x, z1y);
        t4 = fnmadd(ay, z3y, t4);
        t5 = fmadd(ax, z3y, t5);
        t6 = fmsub(2.0, z1x, t4);
        t7 = fmsub(2.0, z1y, t5);
    }
    else
    {
        t4 = t6 = z1x;
        t5 = t7 = z1y;
    }
    t8 = fmadd(bx, t4, t0);
    t4 = fmadd(by, t4, t1);
    t9 = fnmadd(by, t6, t2);
    t6 = fmadd(bx, t6, t3);
    t8 = fnmadd(by, t5, t8);
    t5 = fmadd(bx, t5, t4);
    t9 = fnmadd(bx, t7, t9);
    t7 = fnmadd(by, t7, t6);
    t0 = fmsub(2.0, t0, t8);
    t1 = fmsub(2.0, t1, t5);
    t2 = fmsub(2.0, t2, t9);
    t3 = fmsub(2.0, t3, t7);
    z0.real(t8);
    z0.imag(t5);
    z1.real(t0);
    z1.imag(t1);
    z2.real(t9);
    z2.imag(t7);
    z3.real(t2);
    z3.imag(t3);
}

/*
z0   z0 + ~w2*z1 +   ~w*(z2 + ~w2*z3)
z1   z0 - ~w2*z1 - I*~w*(z2 - ~w2*z3)
z2 = z0 + ~w2*z1 -   ~w*(z2 + ~w2*z3)
z3   z0 - ~w2*z1 + I*~w*(z2 - ~w2*z3)
*/
FORCE_INLINE inline void radix_4_moth_inv_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4>& z2,
    complex_packed<double,4>& z3,
    complex_packed<double,4> w,
    complex_packed<double,4> w2)
{
#if 0
    complex_packed<double,4> s0 = z0, s1 = z1, s2 = z2, s3 = z3;
    radix_2_moth_inv_block(s0, s1, w2);
    radix_2_moth_inv_block(s2, s3, w2);
    radix_2_moth_inv_block(s0, s2, w);
    radix_2_moth_inv_i_block(s1, s3, w);
    z0 = s0;
    z1 = s1;
    z2 = s2;
    z3 = s3;
#else
    packed<double,4> z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    packed<double,4> ax, ay, bx, by, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    ax = w2.real();
    ay = w2.imag();
    bx = w.real();
    by = w.imag();
    z0x = z0.real();
    z0y = z0.imag();
    z1x = z2.real();
    z1y = z2.imag();
    z2x = z1.real();
    z2y = z1.imag();
    z3x = z3.real();
    z3y = z3.imag();

    t0 = fmadd(ax, z2x, z0x);
    t1 = fnmadd(ay, z2x, z0y);
    t0 = fmadd(ay, z2y, t0);
    t1 = fmadd(ax, z2y, t1);
    t2 = fmsub(2.0, z0x, t0);
    t3 = fmsub(2.0, z0y, t1);
    t4 = fmadd(ax, z3x, z1x);
    t5 = fnmadd(ay, z3x, z1y);
    t4 = fmadd(ay, z3y, t4);
    t5 = fmadd(ax, z3y, t5);
    t6 = fmsub(2.0, z1x, t4);
    t7 = fmsub(2.0, z1y, t5);
    t8 = fmadd(bx, t4, t0);
    t4 = fnmadd(by, t4, t1);
    t9 = fmadd(by, t6, t2);
    t6 = fmadd(bx, t6, t3);
    t8 = fmadd(by, t5, t8);
    t5 = fmadd(bx, t5, t4);
    t9 = fnmadd(bx, t7, t9);
    t7 = fmadd(by, t7, t6);
    t0 = fmsub(2.0, t0, t8);
    t1 = fmsub(2.0, t1, t5);
    t2 = fmsub(2.0, t2, t9);
    t3 = fmsub(2.0, t3, t7);
    z0.real(t8);
    z0.imag(t5);
    z2.real(t0);
    z2.imag(t1);
    z3.real(t9);
    z3.imag(t7);
    z1.real(t2);
    z1.imag(t3);
#endif
}

// radix_4_mothi with w = 1
inline void radix_4_moth_inv_1(
    std::complex<double>& z0,
    std::complex<double>& z1,
    std::complex<double>& z2,
    std::complex<double>& z3)
{
    double z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    double s0x, s0y, s1x, s1y, s2x, s2y, s3x, s3y;
    z0x = z0.real();
    z0y = z0.imag();
    z1x = z1.real();
    z1y = z1.imag();
    z2x = z2.real();
    z2y = z2.imag();
    z3x = z3.real();
    z3y = z3.imag();
    s0x = z0x + z1x;
    s0y = z0y + z1y;
    s1x = z0x - z1x;
    s1y = z0y - z1y;
    s2x = z2x + z3x;
    s2y = z2y + z3y;
    s3x = z2x - z3x;
    s3y = z2y - z3y;
    z0.real(s0x + s2x);
    z0.imag(s0y + s2y);
    z2.real(s0x - s2x);
    z2.imag(s0y - s2y);
    z1.real(s1x + s3y);
    z1.imag(s1y - s3x);
    z3.real(s1x - s3y);
    z3.imag(s1y + s3x);
}


// radix_4_moth_inv_block with w = 1
FORCE_INLINE inline void radix_4_moth_inv_1_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4>& z2,
    complex_packed<double,4>& z3)
{
    packed<double,4> z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    packed<double,4> s0x, s0y, s1x, s1y, s2x, s2y, s3x, s3y;
    z0x = z0.real();
    z0y = z0.imag();
    z1x = z1.real();
    z1y = z1.imag();
    z2x = z2.real();
    z2y = z2.imag();
    z3x = z3.real();
    z3y = z3.imag();
    s0x = add(z0x, z1x);
    s0y = add(z0y, z1y);
    s1x = sub(z0x, z1x);
    s1y = sub(z0y, z1y);
    s2x = add(z2x, z3x);
    s2y = add(z2y, z3y);
    s3x = sub(z2x, z3x);
    s3y = sub(z2y, z3y);
    z0.real(add(s0x, s2x));
    z0.imag(add(s0y, s2y));
    z2.real(sub(s0x, s2x));
    z2.imag(sub(s0y, s2y));
    z1.real(add(s1x, s3y));
    z1.imag(sub(s1y, s3x));
    z3.real(sub(s1x, s3y));
    z3.imag(add(s1y, s3x));
}

// radix_4_moth_inv_1_block and then multiply by ti by conj(w)^i
FORCE_INLINE inline void radix_4_moth_rev_block(
    complex_packed<double,4>& z0,
    complex_packed<double,4>& z1,
    complex_packed<double,4>& z2,
    complex_packed<double,4>& z3,
    complex_packed<double,4> w1,
    complex_packed<double,4> w2,
    complex_packed<double,4> w3)
{
    packed<double,4> z0x, z0y, z1x, z1y, z2x, z2y, z3x, z3y;
    packed<double,4> s0x, s0y, s1x, s1y, s2x, s2y, s3x, s3y;
    z0x = z0.real();
    z0y = z0.imag();
    z1x = z1.real();
    z1y = z1.imag();
    z2x = z2.real();
    z2y = z2.imag();
    z3x = z3.real();
    z3y = z3.imag();
    s0x = add(z0x, z1x);
    s0y = add(z0y, z1y);
    s1x = sub(z0x, z1x);
    s1y = sub(z0y, z1y);
    s2x = add(z2x, z3x);
    s2y = add(z2y, z3y);
    s3x = sub(z2x, z3x);
    s3y = sub(z2y, z3y);
    z0x = add(s0x, s2x);
    z0y = add(s0y, s2y);
    z2x = sub(s0x, s2x);
    z2y = sub(s0y, s2y);
    z1x = add(s1x, s3y);
    z1y = sub(s1y, s3x);
    z3x = sub(s1x, s3y);
    z3y = add(s1y, s3x);
    z0.real(z0x);
    z0.imag(z0y);
    z2.real(fmadd(z2y, w2.imag(), mul(z2x, w2.real())));
    z2.imag(fmsub(z2y, w2.real(), mul(z2x, w2.imag())));
    z1.real(fmadd(z1y, w1.imag(), mul(z1x, w1.real())));
    z1.imag(fmsub(z1y, w1.real(), mul(z1x, w1.imag())));
    z3.real(fmadd(z3y, w3.imag(), mul(z3x, w3.real())));
    z3.imag(fmsub(z3y, w3.real(), mul(z3x, w3.imag())));
}


/*
a       a +         b +         c
b =  m*(a + om(1/3)*b + om(2/3)*c)
c   ~m*(a + om(2/3)*b + om(1/3)*c)
*/
template <int itrunc>
FORCE_INLINE inline void radix_3_moth_trunc_block(
    complex_packed<double,4>& a,
    complex_packed<double,4>& b,
    complex_packed<double,4>& c,
    complex_packed<double,4> m)
{
    packed<double,4> ax, ay, bx, by, cx, cy, mx, my;
    packed<double,4> Bx, By, Cx, Cy, bxpcx, bypcy, bxmcx, bxmcy, bymcy, tx, ty;
    if (itrunc == 0)
    {
        a.real(0.0);
        a.imag(0.0);
        b.real(0.0);
        b.imag(0.0);
        c.real(0.0);
        c.imag(0.0);
        return;
    }
    ax = 0 < itrunc ? a.real() : 0.0;
    ay = 0 < itrunc ? a.imag() : 0.0;
    bx = 1 < itrunc ? b.real() : 0.0;
    by = 1 < itrunc ? b.imag() : 0.0;
    cx = 2 < itrunc ? c.real() : 0.0;
    cy = 2 < itrunc ? c.imag() : 0.0;
    mx = m.real();
    my = m.imag();
    bxpcx = 2 < itrunc ? add(bx, cx) : bx;
    bypcy = 2 < itrunc ? add(by, cy) : by;
    a.real(1 < itrunc ? add(ax, bxpcx) : ax);
    a.imag(1 < itrunc ? add(ay, bypcy) : ay);
#if 0
    tx = fnmadd(0.86602540378443864676372317, by, ax);
    ty =  fmadd(0.86602540378443864676372317, bx, ay);
    tx =  fmadd(0.86602540378443864676372317, cy, tx);
    ty = fnmadd(0.86602540378443864676372317, cx, ty);
    By = fnmadd(0.5, bypcy, ty);
    Bx = fnmadd(0.5, bxpcx, tx);
    Cx = fnmadd(0.5, bxpcx, fmsub(2.0, ax, tx));
    Cy = fnmadd(0.5, bypcy, fmsub(2.0, ay, ty));
#else
    bxmcx = 2 < itrunc ? sub(bx, cx) : bx;
    bymcy = 2 < itrunc ? sub(by, cy) : by;
    tx = fnmadd(0.5, bxpcx, ax);
    ty = fnmadd(0.5, bypcy, ay);
    Bx = fnmadd(0.86602540378443864676372317, bymcy, tx);
    By =  fmadd(0.86602540378443864676372317, bxmcx, ty);
    Cx =  fmadd(0.86602540378443864676372317, bymcy, tx);
    Cy = fnmadd(0.86602540378443864676372317, bxmcx, ty);
#endif
    b.real(fmsub(mx, Bx, mul(my, By)));
    b.imag(fmadd(mx, By, mul(my, Bx)));
    c.real(fmadd(mx, Cx, mul(my, Cy)));
    c.imag(fmsub(mx, Cy, mul(my, Cx)));
}


/*
a   a +         ~m*b +         m*c
b = a + om(2/3)*~m*b + om(1/3)*m*c
c   a + om(1/3)*~m*b + om(2/3)*m*c
*/
FORCE_INLINE inline void radix_3_moth_inv_block(
    complex_packed<double,4>& a,
    complex_packed<double,4>& b,
    complex_packed<double,4>& c,
    complex_packed<double,4> m)
{
    packed<double,4> ax, ay, bx, by, cx, cy, mx, my;
    packed<double,4> Bx, By, BxpCx, BxmCx, BypCy, BymCy, tx, ty;
    ax = a.real();
    ay = a.imag();
    bx = b.real();
    by = b.imag();
    cx = c.real();
    cy = c.imag();
    mx = m.real();
    my = m.imag();

    tx = fmadd(mx, bx, mul(my, by));
    ty = fmsub(mx, by, mul(my, bx));
    BxmCx =  fmadd(cy, my, fnmadd(cx, mx, tx));
    BymCy = fnmadd(cx, my, fnmadd(cy, mx, ty));
#ifdef OPT_MOTH_2
    BxpCx = fmsub(2.0, tx, BxmCx);
    BypCy = fmsub(2.0, ty, BymCy);
#else
    BxpCx = fnmadd(cy, my, fmadd(cx, mx, tx));
    BypCy =  fmadd(cx, my, fmadd(cy, mx, ty));
#endif
    a.real(add(ax, BxpCx));
    a.imag(add(ay, BypCy));
    tx = fnmadd(0.5, BxpCx, ax);
    ty = fnmadd(0.5, BypCy, ay);
#if 0
    Bx =  fmadd(0.86602540378443864676372317, BymCy, tx);
    By = fnmadd(0.86602540378443864676372317, BxmCx, ty);
    b.real(Bx);
    b.imag(By);
    c.real(fmsub(2.0, tx, Bx));
    c.imag(fmsub(2.0, ty, By));
#else
    b.real( fmadd(0.86602540378443864676372317, BymCy, tx));
    b.imag(fnmadd(0.86602540378443864676372317, BxmCx, ty));
    c.real(fnmadd(0.86602540378443864676372317, BymCy, tx));
    c.imag( fmadd(0.86602540378443864676372317, BxmCx, ty));
#endif
}


FORCE_INLINE inline void cmplx_mul(
    packed<double,4>& zx, packed<double,4>& zy,
    packed<double,4> ax, packed<double,4> ay,
    packed<double,4> bx, packed<double,4> by)
{
    zx = fmsub(ax, bx, mul(ay, by));
    zy = fmadd(ay, bx, mul(ax, by));
}

FORCE_INLINE inline void cmplx_mulc(
    packed<double,4>& zx, packed<double,4>& zy,
    packed<double,4> ax, packed<double,4> ay,
    packed<double,4> bx, packed<double,4> by)
{
    zx = fmadd(ax, bx, mul(ay, by));
    zy = fmsub(ay, bx, mul(ax, by));
}

FORCE_INLINE inline void cmplx_addmul(
    packed<double,4>& zx, packed<double,4>& zy,
    packed<double,4> ax, packed<double,4> ay,
    packed<double,4> bx, packed<double,4> by)
{
    packed<double,4> zzx = zx, zzy = zy;
    zx = fnmadd(ay, by, fmadd(ax, bx, zzx));
    zy =  fmadd(ay, bx, fmadd(ax, by, zzy));

}

/*
A       (a +  om(0/5)*b +  om(0/5)*c +  om(0/5)*d +  om(0/5)*e)
B   m  *(a +  om(1/5)*b +  om(2/5)*c +  om(3/5)*d +  om(4/5)*e)
C = ~m *(a + om(-1/5)*b + om(-2/5)*c + om(-3/5)*d + om(-4/5)*e)
D   m2 *(a +  om(2/5)*b +  om(4/5)*c +  om(6/5)*d +  om(8/5)*e)
E   ~m2*(a + om(-2/5)*b + om(-4/5)*c + om(-6/5)*d + om(-8/5)*e)
*/
template <int itrunc>
FORCE_INLINE inline void radix_5_moth_trunc_block(
    complex_packed<double,4>& a,
    complex_packed<double,4>& b,
    complex_packed<double,4>& c,
    complex_packed<double,4>& d,
    complex_packed<double,4>& e,
    complex_packed<double,4> m,
    complex_packed<double,4> m2)
{
    packed<double,4> ax, ay, bx, by, cx, cy, dx, dy, ex, ey;
    packed<double,4> Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex, Ey;

    if (itrunc == 0)
    {
        a.real(0.0);
        a.imag(0.0);
        b.real(0.0);
        b.imag(0.0);
        c.real(0.0);
        c.imag(0.0);
        d.real(0.0);
        d.imag(0.0);
        e.real(0.0);
        e.imag(0.0);
        return;
    }
    ax = 0 < itrunc ? a.real() : 0.0;
    ay = 0 < itrunc ? a.imag() : 0.0;
    bx = 1 < itrunc ? b.real() : 0.0;
    by = 1 < itrunc ? b.imag() : 0.0;
    cx = 2 < itrunc ? c.real() : 0.0;
    cy = 2 < itrunc ? c.imag() : 0.0;
    dx = 3 < itrunc ? d.real() : 0.0;
    dy = 3 < itrunc ? d.imag() : 0.0;
    ex = 4 < itrunc ? e.real() : 0.0;
    ey = 4 < itrunc ? e.imag() : 0.0;

    a.real(add(add(ax, bx), add(add(cx, dx), ex)));
    a.imag(add(add(ay, by), add(add(cy, dy), ey)));

    Bx = ax; By = ay;
    if (1 < itrunc) cmplx_addmul(Bx,By, bx,by,  0.30901699437494742410,  0.95105651629515357212);
    if (2 < itrunc) cmplx_addmul(Bx,By, cx,cy, -0.80901699437494742410,  0.58778525229247312917);
    if (3 < itrunc) cmplx_addmul(Bx,By, dx,dy, -0.80901699437494742410, -0.58778525229247312917);
    if (4 < itrunc) cmplx_addmul(Bx,By, ex,ey,  0.30901699437494742410, -0.95105651629515357212);

    Cx = ax; Cy = ay;
    if (1 < itrunc) cmplx_addmul(Cx,Cy, bx,by,  0.30901699437494742410, -0.95105651629515357212);
    if (2 < itrunc) cmplx_addmul(Cx,Cy, cx,cy, -0.80901699437494742410, -0.58778525229247312917);
    if (3 < itrunc) cmplx_addmul(Cx,Cy, dx,dy, -0.80901699437494742410,  0.58778525229247312917);
    if (4 < itrunc) cmplx_addmul(Cx,Cy, ex,ey,  0.30901699437494742410,  0.95105651629515357212);

    Dx = ax; Dy = ay;
    if (1 < itrunc) cmplx_addmul(Dx,Dy, bx,by, -0.80901699437494742410,  0.58778525229247312917);
    if (2 < itrunc) cmplx_addmul(Dx,Dy, cx,cy,  0.30901699437494742410, -0.95105651629515357212);
    if (3 < itrunc) cmplx_addmul(Dx,Dy, dx,dy,  0.30901699437494742410,  0.95105651629515357212);
    if (4 < itrunc) cmplx_addmul(Dx,Dy, ex,ey, -0.80901699437494742410, -0.58778525229247312917);

    Ex = ax; Ey = ay;
    if (1 < itrunc) cmplx_addmul(Ex,Ey, bx,by, -0.80901699437494742410, -0.58778525229247312917);
    if (2 < itrunc) cmplx_addmul(Ex,Ey, cx,cy,  0.30901699437494742410,  0.95105651629515357212);
    if (3 < itrunc) cmplx_addmul(Ex,Ey, dx,dy,  0.30901699437494742410, -0.95105651629515357212);
    if (4 < itrunc) cmplx_addmul(Ex,Ey, ex,ey, -0.80901699437494742410,  0.58778525229247312917);

    cmplx_mul (bx,by, Bx,By, m.real(), m.imag());
    cmplx_mulc(cx,cy, Cx,Cy, m.real(), m.imag());
    cmplx_mul (dx,dy, Dx,Dy, m2.real(), m2.imag());
    cmplx_mulc(ex,ey, Ex,Ey, m2.real(), m2.imag());

    b.real(bx); b.imag(by);
    c.real(cx); c.imag(cy);
    d.real(dx); d.imag(dy);
    e.real(ex); e.imag(ey);
}

/*
A       (a +  om(0/5)*b +  om(0/5)*c +  om(0/5)*d +  om(0/5)*e)
B   m  *(a +  om(1/5)*b +  om(2/5)*c +  om(3/5)*d +  om(4/5)*e)
C = m~ *(a + om(-1/5)*b + om(-2/5)*c + om(-3/5)*d + om(-4/5)*e)
D   m2 *(a +  om(2/5)*b +  om(4/5)*c +  om(6/5)*d +  om(8/5)*e)
E   m2~*(a + om(-2/5)*b + om(-4/5)*c + om(-6/5)*d + om(-8/5)*e)

a   A + m~*B + m*C + m2~*D + m2*E
b       -1/5   1/5   -2/5    2/5
c =     -2/5   2/5   -4/5    4/5
d       -3/5   3/5   -6/5    6/5
e       -4/5   4/5   -8/5    8/5
*/
FORCE_INLINE inline void radix_5_moth_inv_block(
    complex_packed<double,4>& a,
    complex_packed<double,4>& b,
    complex_packed<double,4>& c,
    complex_packed<double,4>& d,
    complex_packed<double,4>& e,
    complex_packed<double,4> m,
    complex_packed<double,4> m2)
{
    packed<double,4> ax, ay, bx, by, cx, cy, dx, dy, ex, ey;
    packed<double,4> Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex, Ey;

    Ax = a.real();
    Ay = a.imag();
    cmplx_mulc(Bx,By, b.real(),b.imag(), m.real(),m.imag());
    cmplx_mul (Cx,Cy, c.real(),c.imag(), m.real(), m.imag());
    cmplx_mulc(Dx,Dy, d.real(),d.imag(), m2.real(),m2.imag());
    cmplx_mul (Ex,Ey, e.real(),e.imag(), m2.real(), m2.imag());

    a.real(add(add(Ax, Bx), add(add(Cx, Dx), Ex)));
    a.imag(add(add(Ay, By), add(add(Cy, Dy), Ey)));

    bx = Ax; by = Ay;
    cmplx_addmul(bx,by, Bx,By, -0.80901699437494742410, -0.58778525229247312917);
    cmplx_addmul(bx,by, Cx,Cy, -0.80901699437494742410,  0.58778525229247312917);
    cmplx_addmul(bx,by, Dx,Dy,  0.30901699437494742410,  0.95105651629515357212);
    cmplx_addmul(bx,by, Ex,Ey,  0.30901699437494742410, -0.95105651629515357212);

    cx = Ax; cy = Ay;
    cmplx_addmul(cx,cy, Bx,By,  0.30901699437494742410, +0.95105651629515357212);
    cmplx_addmul(cx,cy, Cx,Cy,  0.30901699437494742410, -0.95105651629515357212);
    cmplx_addmul(cx,cy, Dx,Dy, -0.80901699437494742410,  0.58778525229247312917);
    cmplx_addmul(cx,cy, Ex,Ey, -0.80901699437494742410, -0.58778525229247312917);

    dx = Ax; dy = Ay;
    cmplx_addmul(dx,dy, Bx,By,  0.30901699437494742410, -0.95105651629515357212);
    cmplx_addmul(dx,dy, Cx,Cy,  0.30901699437494742410,  0.95105651629515357212);
    cmplx_addmul(dx,dy, Dx,Dy, -0.80901699437494742410, -0.58778525229247312917);
    cmplx_addmul(dx,dy, Ex,Ey, -0.80901699437494742410,  0.58778525229247312917);

    ex = Ax; ey = Ay;
    cmplx_addmul(ex,ey, Bx,By, -0.80901699437494742410,  0.58778525229247312917);
    cmplx_addmul(ex,ey, Cx,Cy, -0.80901699437494742410, -0.58778525229247312917);
    cmplx_addmul(ex,ey, Dx,Dy,  0.30901699437494742410, -0.95105651629515357212);
    cmplx_addmul(ex,ey, Ex,Ey,  0.30901699437494742410,  0.95105651629515357212);

    // TODO ??
    b.real(dx); b.imag(dy);
    c.real(bx); c.imag(by);
    d.real(ex); d.imag(ey);
    e.real(cx); e.imag(cy);
}

