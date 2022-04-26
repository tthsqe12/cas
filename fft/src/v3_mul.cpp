// break into coefficients from the Kronecker substitution at 2^B
template<ulong B>
void mpn_to_cpd(
    packed<double,4>* z, ulong zn,
    const ulong* a_, ulong an_)
{
    static_assert(0 < B && B < 32);
    assert(zn%32 == 0);
    assert(an_ > 0);
    constexpr ulong stride = sizeof(complex_packed<double,4>)/sizeof(packed<double,4>);
    const uint32_t* a = reinterpret_cast<const uint32_t*>(a_);
    ulong an = 2*an_;

#define CODE(ir) \
        {\
            ulong a0 = *reinterpret_cast<const ulong*>(a+B*i + (B*(4*ir+0))/32);\
            constexpr uint32_t r0 =                            (B*(4*ir+0))%32;\
            ulong a1 = *reinterpret_cast<const ulong*>(a+B*i + (B*(4*ir+1))/32);\
            constexpr uint32_t r1 =                            (B*(4*ir+1))%32;\
            ulong a2 = *reinterpret_cast<const ulong*>(a+B*i + (B*(4*ir+2))/32);\
            constexpr uint32_t r2 =                            (B*(4*ir+2))%32;\
            ulong a3 = *reinterpret_cast<const ulong*>(a+B*i + (B*(4*ir+3))/32);\
            constexpr uint32_t r3 =                            (B*(4*ir+3))%32;\
            r[ir] = packed<ulong,4>(a0,a1,a2,a3);\
            r[ir] = bit_shift_right(r[ir], packed<ulong,4>(r0,r1,r2,r3));\
            r[ir] = bit_and(mask, r[ir]);\
        }
    ulong i = 0;

    // if B*i + fld(B*31,32) + 1 < an, then aindex is in range
    ulong end_easy = std::min(zn/32, std::max(ulong(1), (an-1)/B) - 1);

    // if B*i >= an, then aindex is trivial
    ulong end_hard = std::min(zn/32, cdiv(an, B));

    packed<ulong,4> mask = pow2(B) - 1;
    for (; i < end_easy; i++)
    {
        packed<ulong,4> r[8];
        CODE(0); CODE(1); CODE(2); CODE(3); CODE(4); CODE(5); CODE(6); CODE(7);
        for (ulong ir = 0; ir < 8; ir++)
            z[stride*(8*i+ir)] = convert_limited<packed<double,4>>(r[ir]);
    }
#undef CODE

#define CODE(ir) \
        {\
            ulong q = B*i + (B*(ir))/32;\
            ulong r =       (B*(ir))%32;\
            if (r == 0)\
            {\
                res[ir] = aindex(q) & (pow2(B) - 1);\
            }\
            else if (r + B <= 32)\
            {\
                if (r + B == 32)\
                    res[ir] = (aindex(q) >> r);\
                else\
                    res[ir] = (aindex(q) >> r) & (pow2(B) - 1);\
            }\
            else\
            {\
                res[ir] = (aindex(q) >> r) |\
                             ((aindex(q+1) & (pow2(r+B-32) - 1)) << (32 - r));\
            }\
        }
#define aindex(q)  ((q < an) ? a[q] : uint32_t(0))
    for (; i < end_hard; i++)
    {
        uint32_t res[32];
        CODE(8*0+0);CODE(8*0+1);CODE(8*0+2);CODE(8*0+3);CODE(8*0+4);CODE(8*0+5);CODE(8*0+6);CODE(8*0+7);
        CODE(8*1+0);CODE(8*1+1);CODE(8*1+2);CODE(8*1+3);CODE(8*1+4);CODE(8*1+5);CODE(8*1+6);CODE(8*1+7);
        CODE(8*2+0);CODE(8*2+1);CODE(8*2+2);CODE(8*2+3);CODE(8*2+4);CODE(8*2+5);CODE(8*2+6);CODE(8*2+7);
        CODE(8*3+0);CODE(8*3+1);CODE(8*3+2);CODE(8*3+3);CODE(8*3+4);CODE(8*3+5);CODE(8*3+6);CODE(8*3+7);

        for (ulong k = 0; k < 8; k++)
        {
            packed<ulong,4> rr = packed<ulong,4>(res[4*k+0], res[4*k+1], res[4*k+2], res[4*k+3]);
            z[stride*(8*i+k)] = convert_limited<packed<double,4>>(rr);
        }
    }
#undef aindex
#undef CODE

    for (; i < zn/32; i++)
    {
        for (ulong k = 0; k < 8; k++)
            z[stride*(8*i+k)] = 0.0;
    }
}

inline void mpn_to_cpd(
    packed<double,4>* z, ulong zn,
    const ulong* a_, ulong an_,
    ulong bits)
{
    static void (*tab[22-8+1])(packed<double,4>*, ulong, const ulong*, ulong) {
        mpn_to_cpd<8>,
        mpn_to_cpd<9>,
        mpn_to_cpd<10>,
        mpn_to_cpd<11>,
        mpn_to_cpd<12>,
        mpn_to_cpd<13>,
        mpn_to_cpd<14>,
        mpn_to_cpd<15>,
        mpn_to_cpd<16>,
        mpn_to_cpd<17>,
        mpn_to_cpd<18>,
        mpn_to_cpd<19>,
        mpn_to_cpd<20>,
        mpn_to_cpd<21>,
        mpn_to_cpd<22>};

    assert(8 <= bits); assert(bits <= 22);
    tab[bits-8](z, zn, a_, an_);
}


// reassemble from coefficients from Kronecker substitution at 2^B
// the entries of a should be real integers in [0, 2^M) (hopefully)
template<ulong B>
void mpn_from_cpd(
    ulong* z_, ulong zn_,
    const complex_packed<double,4>* a, ulong an)
{
    static_assert(0 < B && B < 32);
    assert(an%32 == 0);
    uint32_t* z = reinterpret_cast<uint32_t*>(z_);
    ulong zn = 2*zn_;

    ulong x0 = 0, x1 = 0;
    ulong i = 0;
    ulong j = 0;
    ulong end_easy = std::min(an/32, zn/B);
    for (; j < end_easy; j++)
    {
        packed<ulong,4> ax[4], ay[4];
        for (ulong k = 0; k < 4; k++)
        {
            ax[k] = convert_limited<packed<ulong,4>>(round(a[4*j+k].real()));
            ay[k] = convert_limited<packed<ulong,4>>(round(a[4*j+k].imag()));
        }

#define CODE(r)\
        {\
            ulong o = ((r)*B)%32;\
            ulong A = ((r)&1) ? ay[(r)/8][(((r)-1)/2)%4] : ax[(r)/8][((r)/2)%4];\
            x0 += (A & (pow2(32-o)-1)) << o;\
            x1 += A >> (32-o);\
            if (o+B >= 32)\
            {\
                assert(i < zn);\
                z[i] = x0;\
                x0 = (x0 >> 32) + x1;\
                x1 = 0;\
                i++;\
            }\
        }
        CODE(8*0+0);CODE(8*0+1);CODE(8*0+2);CODE(8*0+3);CODE(8*0+4);CODE(8*0+5);CODE(8*0+6);CODE(8*0+7);
        CODE(8*1+0);CODE(8*1+1);CODE(8*1+2);CODE(8*1+3);CODE(8*1+4);CODE(8*1+5);CODE(8*1+6);CODE(8*1+7);
        CODE(8*2+0);CODE(8*2+1);CODE(8*2+2);CODE(8*2+3);CODE(8*2+4);CODE(8*2+5);CODE(8*2+6);CODE(8*2+7);
        CODE(8*3+0);CODE(8*3+1);CODE(8*3+2);CODE(8*3+3);CODE(8*3+4);CODE(8*3+5);CODE(8*3+6);CODE(8*3+7);
#undef CODE

        assert(B*(j + 1) == i);
    }

    for (; j < an/32; j++)
    {
        packed<ulong,4> ax[4], ay[4];
        for (ulong k = 0; k < 4; k++)
        {
            ax[k] = convert_limited<packed<ulong,4>>(round(a[4*j+k].real()));
            ay[k] = convert_limited<packed<ulong,4>>(round(a[4*j+k].imag()));
        }

#define CODE(r)\
        {\
            ulong o = ((r)*B)%32;\
            ulong A = ((r)&1) ? ay[(r)/8][(((r)-1)/2)%4] : ax[(r)/8][((r)/2)%4];\
            x0 += (A & (pow2(32-o)-1)) << o;\
            x1 += A >> (32-o);\
            if (o+B >= 32)\
            {\
                z[i] = x0;\
                x0 = (x0 >> 32) + x1;\
                x1 = 0;\
                if (++i >= zn)\
                    return;\
            }\
        }
        CODE(8*0+0);CODE(8*0+1);CODE(8*0+2);CODE(8*0+3);CODE(8*0+4);CODE(8*0+5);CODE(8*0+6);CODE(8*0+7);
        CODE(8*1+0);CODE(8*1+1);CODE(8*1+2);CODE(8*1+3);CODE(8*1+4);CODE(8*1+5);CODE(8*1+6);CODE(8*1+7);
        CODE(8*2+0);CODE(8*2+1);CODE(8*2+2);CODE(8*2+3);CODE(8*2+4);CODE(8*2+5);CODE(8*2+6);CODE(8*2+7);
        CODE(8*3+0);CODE(8*3+1);CODE(8*3+2);CODE(8*3+3);CODE(8*3+4);CODE(8*3+5);CODE(8*3+6);CODE(8*3+7);
#undef CODE

        assert(B*(j + 1) == i);
    }

    return;
}

void mpn_from_cpd(
    ulong* z_, ulong zn_,
    const complex_packed<double,4>* a, ulong an,
    ulong bits)
{
    static void (*tab[22-8+1])(ulong*, ulong, const complex_packed<double,4>*, ulong) {
        mpn_from_cpd<8>,
        mpn_from_cpd<9>,
        mpn_from_cpd<10>,
        mpn_from_cpd<11>,
        mpn_from_cpd<12>,
        mpn_from_cpd<13>,
        mpn_from_cpd<14>,
        mpn_from_cpd<15>,
        mpn_from_cpd<16>,
        mpn_from_cpd<17>,
        mpn_from_cpd<18>,
        mpn_from_cpd<19>,
        mpn_from_cpd<20>,
        mpn_from_cpd<21>,
        mpn_from_cpd<22>};

    assert(8 <= bits); assert(bits <= 22);
    tab[bits-8](z_, zn_, a, an);
}


// m * (-0.25*im)*(f + conj(g))*(f - conj(g))
FORCE_INLINE inline std::complex<double> mulconj(
    std::complex<double> f,
    std::complex<double> g,
    double m)
{
// -fx^2*i + 2*fx*fy + fy^2*i + gx^2*i + 2*gx*gy - gy^2*i
    double fx = f.real();
    double fy = f.imag();
    double gx = g.real();
    double gy = g.imag();
    double rx = (0.50*m)*fmadd(gx, gy, mul(fx, fy));
    double ry = (0.25*m)*fmadd(fy, fy, fnmadd(fx, fx, fnmadd(gy, gy, mul(gx, gx))));
    return std::complex<double>(rx, ry);
}

template <typename T>
FORCE_INLINE inline void mul_conj(
    T& rx, T& ry,
    T fx, T fy,
    T gx, T gy,
    double m)
{
    rx = (0.50*m)*fmadd(gx, gy, mul(fx, fy));
    ry = (0.25*m)*fmadd(fy, fy, fnmadd(fx, fx, fnmadd(gy, gy, mul(gx, gx))));
}

/*
    Recover the individual ffts, multiply them, and pack down to ifft with
    complex result of half length.
*/
static void _recover_and_mul_and_pack(
    complex_packed<double,4>* A,    // input and output
    ulong l,
    const std::complex<double>* w2revs,
    double recpN)
{
    assert(l >= 16);

    std::complex<double> Z[16]; // for lazy approach

    for (ulong i = 0; i < 16/4; i++)
    for (ulong j = 0; j < 4; j++)
    {
        Z[4*i+j].real(A[i].real()[j]);
        Z[4*i+j].imag(A[i].imag()[j]);
    }

    // 0, 1
    {
        double c1 = Z[0].real()*Z[0].imag()*recpN;
        double c2 = Z[1].real()*Z[1].imag()*recpN;
        double ei = (c1 + c2);
        double oi = (c1 - c2);
        Z[0].real(ei);
        Z[0].imag(oi);
    }

    // [2,4)
    {
        double aix, aiy, ahx, ahy;
        aix = Z[2].real();
        ahx = Z[3].real();
        aiy = Z[2].imag();
        ahy = Z[3].imag();
        double c1x = fmadd(aix, aiy, mul(ahx, ahy));
        double c1y = fmadd(aiy, aiy, fnmadd(aix, aix, fnmadd(ahy, ahy, mul(ahx, ahx))));
        double ei = c1x*(1.0*recpN);
        double oi = c1y*(0.5*recpN);
        Z[1].real(ei);
        Z[1].imag(oi);
    }

    // [4,8), [8,16)
    for (ulong j = 4; j < 16; j *= 2)
    for (ulong i = 0; i < j/2; i += 2)
    {
        std::complex<double> r1, r3, ei, oi, fi, pi;
        double aix, ahx, bix, bhx, aiy, ahy, biy, bhy;
        aix = Z[j+i].real();
        ahx = Z[2*j-1-i].real();
        bix = Z[j+i+1].real();
        bhx = Z[2*j-2-i].real();
        aiy = Z[j+i].imag();
        ahy = Z[2*j-1-i].imag();
        biy = Z[j+i+1].imag();
        bhy = Z[2*j-2-i].imag();
        double r1x = (0.50*recpN)*fmadd(aix, aiy, mul(ahx, ahy));
        double r1y = (0.25*recpN)*fmadd(aiy, aiy, fnmadd(aix, aix, fnmadd(ahy, ahy, mul(ahx, ahx))));
        double r3x = (0.50*recpN)*fmadd(bix, biy, mul(bhx, bhy));
        double r3y = (0.25*recpN)*fmadd(biy, biy, fnmadd(bix, bix, fnmadd(bhy, bhy, mul(bhx, bhx))));
        double eix = r1x + r3x;
        double eiy = r1y + r3y;
        double pix = r1x - r3x;
        double piy = r1y - r3y;
        double qx = w2revs[(j+i)/2].real();
        double qy = w2revs[(j+i)/2].imag();
        double oix =  fmadd(pix, qx, piy*qy);
        double oiy = fnmadd(pix, qy, piy*qx);
        Z[(j+i)/2].real(eix - oiy);
        Z[(j+i)/2].imag(oix + eiy);
        Z[(2*j-i-2)/2].real(eix + oiy);
        Z[(2*j-i-2)/2].imag(oix - eiy);
    }

    for (ulong i = 0; i < 2; i++)
    {
        A[i].real(packed<double,4>(Z[4*i+0].real(), Z[4*i+1].real(), Z[4*i+2].real(), Z[4*i+3].real()));
        A[i].imag(packed<double,4>(Z[4*i+0].imag(), Z[4*i+1].imag(), Z[4*i+2].imag(), Z[4*i+3].imag()));
    }

    // end being lazy

    for (ulong j = 16; j < l; j *= 2)
    for (ulong i = 0; i < j/2; i += 8)
    {
        packed<double,4> zlx, zly, zhx, zhy;
        packed<double,4> aix, ahx, bix, bhx, aiy, ahy, biy, bhy;
        packed<double,4> r1x, r1y, r3x, r3y, eix, eiy, pix, piy, qx, qy, oix, oiy;
        packed<double,4> A0x, A0y, A1x, A1y, B0x, B0y, B1x, B1y;

        A0x = A[(j+i)/4+0].real();
        A0y = A[(j+i)/4+0].imag();
        A1x = A[(j+i)/4+1].real();
        A1y = A[(j+i)/4+1].imag();

        aix = unpacklo(A0x, A1x);
        bix = unpackhi(A0x, A1x);
        aiy = unpacklo(A0y, A1y);
        biy = unpackhi(A0y, A1y);

        B0x = A[(2*j-i)/4-1].real();
        B0y = A[(2*j-i)/4-1].imag();
        B1x = A[(2*j-i)/4-2].real();
        B1y = A[(2*j-i)/4-2].imag();
        B0x = permute2<1,0>(B0x, B0x);
        B0y = permute2<1,0>(B0y, B0y);
        B1x = permute2<1,0>(B1x, B1x);
        B1y = permute2<1,0>(B1y, B1y);
        ahx = unpackhi(B0x, B1x);
        bhx = unpacklo(B0x, B1x);
        ahy = unpackhi(B0y, B1y);
        bhy = unpacklo(B0y, B1y);

        qx = packed<double,4>(w2revs[(j+i)/2+0].real(), w2revs[(j+i)/2+2].real(), w2revs[(j+i)/2+1].real(), w2revs[(j+i)/2+3].real());
        qy = packed<double,4>(w2revs[(j+i)/2+0].imag(), w2revs[(j+i)/2+2].imag(), w2revs[(j+i)/2+1].imag(), w2revs[(j+i)/2+3].imag());

        r1x = mul(0.50*recpN, fmadd(aix, aiy, mul(ahx, ahy)));
        r1y = mul(0.25*recpN, fmadd(aiy, aiy, fnmadd(aix, aix, fnmadd(ahy, ahy, mul(ahx, ahx)))));
        r3x = mul(0.50*recpN, fmadd(bix, biy, mul(bhx, bhy)));
        r3y = mul(0.25*recpN, fmadd(biy, biy, fnmadd(bix, bix, fnmadd(bhy, bhy, mul(bhx, bhx)))));
        eix = add(r1x, r3x);
        eiy = add(r1y, r3y);
        pix = sub(r1x, r3x);
        piy = sub(r1y, r3y);
        oix = fmadd(pix, qx, mul(piy, qy));
        oiy = fnmadd(pix, qy, mul(piy, qx));
        zlx = sub(eix, oiy);
        zly = add(eiy, oix);
        zhx = add(eix, oiy);
        zhy = sub(oix, eiy);

        A[(j+i)/8].real(permute<0,2,1,3>(zlx));
        A[(j+i)/8].imag(permute<0,2,1,3>(zly));
        A[(2*j-i)/8-1].real(permute<3,1,2,0>(zhx));
        A[(2*j-i)/8-1].imag(permute<3,1,2,0>(zhy));
    }
}

static void _recover_and_mul_and_pack(
    complex_packed<double,4>* Z,    // upper half may overlap A
    const complex_packed<double,4>* A,
    const complex_packed<double,4>* B,
    ulong l,
    double recpN,
    const std::complex<double>* tab1,
    const std::complex<double>* tab2)
{
    assert(l >= 16);

    ulong ind;
    std::complex<double> cc, dd, s, m;
    double qx, qy;

#define f1(i) std::complex<double>(A[(i)/4].real()[(i)%4], A[(i)/4].imag()[(i)%4])
#define f2(i) std::complex<double>(B[(i)/4].real()[(i)%4], B[(i)/4].imag()[(i)%4])

    std::complex<double> R1[4], R2[4];
    cc = mulconj(f1(0), f2(0), recpN);
    dd = mulconj(f1(1), f2(1), recpN);
    m = cc - dd;
    s = cc + dd;
    m = conj(tab1[0])*m;
    R1[0].real(s.real() - m.imag());
    R2[0].real(s.real() + m.imag());
    R1[0].imag(s.imag() + m.real());
    R2[0].imag(-s.imag() + m.real());

    for (ulong j = 2; j < 8; j *= 2)
    for (ulong i = 0; i < j; i += 2)
    {
        cc = mulconj(f1(j+(i+0)), f2(2*j-1-(i+0)), recpN);
        dd = mulconj(f1(j+(i+1)), f2(2*j-1-(i+1)), recpN);
        s = cc + dd;
        m = cc - dd;
        ind = (j+i)/2;
//        R1[ind] = s + std::complex<double>(0,1)*conj(tab1[ind])*m;
        qx = tab1[ind].real();
        qy = tab1[ind].imag();
        R1[ind].real(fmadd(qy, m.real(), fnmadd(qx, m.imag(), s.real())));
        R1[ind].imag(fmadd(qy, m.imag(),  fmadd(qx, m.real(), s.imag())));

        ind = (2*j-i-2)/2;
//        R2[ind] = conj(s + std::complex<double>(0,1)*tab2[ind]*m);
        qx = tab2[ind].real();
        qy = tab2[ind].imag();
        R2[ind].real(fnmadd(qy, m.real(), fnmadd(qx, m.imag(), s.real())));
        R2[ind].imag(fnmadd(qx, m.real(),  fmsub(qy, m.imag(), s.imag())));
        
    }

    for (ulong i = 0; i < 4; i += 4)
    {
        Z[(0*l/2+i)/4].real(packed<double,4>(R1[i+0].real(),R1[i+1].real(),R1[i+2].real(),R1[i+3].real()));
        Z[(0*l/2+i)/4].imag(packed<double,4>(R1[i+0].imag(),R1[i+1].imag(),R1[i+2].imag(),R1[i+3].imag()));
        Z[(1*l/2+i)/4].real(packed<double,4>(R2[i+0].real(),R2[i+1].real(),R2[i+2].real(),R2[i+3].real()));
        Z[(1*l/2+i)/4].imag(packed<double,4>(R2[i+0].imag(),R2[i+1].imag(),R2[i+2].imag(),R2[i+3].imag()));
    }

    for (ulong j = 8; j < l; j *= 2)
    for (ulong i = 0; i < j; i += 8)
    {
/*
        std::complex<double> cc_0 = mulconj(f1(j+(i+0)), f2(2*j-1-(i+0)), recpN);
        std::complex<double> dd_0 = mulconj(f1(j+(i+1)), f2(2*j-1-(i+1)), recpN);
        std::complex<double> cc_1 = mulconj(f1(j+(i+2)), f2(2*j-1-(i+2)), recpN);
        std::complex<double> dd_1 = mulconj(f1(j+(i+3)), f2(2*j-1-(i+3)), recpN);
        std::complex<double> cc_2 = mulconj(f1(j+(i+4)), f2(2*j-1-(i+4)), recpN);
        std::complex<double> dd_2 = mulconj(f1(j+(i+5)), f2(2*j-1-(i+5)), recpN);
        std::complex<double> cc_3 = mulconj(f1(j+(i+6)), f2(2*j-1-(i+6)), recpN);
        std::complex<double> dd_3 = mulconj(f1(j+(i+7)), f2(2*j-1-(i+7)), recpN);

        ulong indA = (j+i)/2;   // 0 mod 4
        ulong indB = (2*j-i)/2; // 0 mod 4

        std::complex<double> r1_0, r2_0, r1_1, r2_1, r1_2, r2_2, r1_3, r2_3;

        r1_0 = (cc_0+dd_0) + std::complex<double>(0,1)*conj(tab1[indA])*(cc_0-dd_0);
        r2_3 = conj((cc_0+dd_0) + std::complex<double>(0,1)*tab2[indB-1]*(cc_0-dd_0));

        r1_1 = (cc_1+dd_1) + std::complex<double>(0,1)*conj(tab1[indA+1])*(cc_1-dd_1);
        r2_2 = conj((cc_1+dd_1) + std::complex<double>(0,1)*tab2[indB-2]*(cc_1-dd_1));

        r1_2 = (cc_2+dd_2) + std::complex<double>(0,1)*conj(tab1[indA+2])*(cc_2-dd_2);
        r2_1 = conj((cc_2+dd_2) + std::complex<double>(0,1)*tab2[indB-3]*(cc_2-dd_2));

        r1_3 = (cc_3+dd_3) + std::complex<double>(0,1)*conj(tab1[indA+3])*(cc_3-dd_3);
        r2_0 = conj((cc_3+dd_3) + std::complex<double>(0,1)*tab2[indB-4]*(cc_3-dd_3));

        A[(1*l/2+indA)/4].real(packed<double,4>(r1_0.real(),r1_1.real(),r1_2.real(),r1_3.real()));
        A[(1*l/2+indA)/4].imag(packed<double,4>(r1_0.imag(),r1_1.imag(),r1_2.imag(),r1_3.imag()));
        A[(2*l/2+indB-4)/4].real(packed<double,4>(r2_0.real(),r2_1.real(),r2_2.real(),r2_3.real()));
        A[(2*l/2+indB-4)/4].imag(packed<double,4>(r2_0.imag(),r2_1.imag(),r2_2.imag(),r2_3.imag()));
*/
        packed<double,4> A0x, A0y, A1x, A1y, a_0246x, a_1357x, a_0246y, a_1357y;
        packed<double,4> B0x, B0y, B1x, B1y, b_0246x, b_1357x, b_0246y, b_1357y;
        packed<double,4> cc_0123x, cc_0123y, dd_0123x, dd_0123y;
        packed<double,4> s_0123x, s_0123y, m_0123x, m_0123y;
        packed<double,4> qa0123x, qa0123y, qb1234x, qb1234y;
        packed<double,4> r1_0123x, r1_0123y, r2_3210x, r2_3210y;

        A0x = A[(j+i)/4+0].real();
        A0y = A[(j+i)/4+0].imag();
        A1x = A[(j+i)/4+1].real();
        A1y = A[(j+i)/4+1].imag();
        a_0246x = unpacklo(A0x, A1x);
        a_1357x = unpackhi(A0x, A1x);
        a_0246y = unpacklo(A0y, A1y);
        a_1357y = unpackhi(A0y, A1y);

        B0x = B[(2*j-i)/4-1].real();
        B0y = B[(2*j-i)/4-1].imag();
        B1x = B[(2*j-i)/4-2].real();
        B1y = B[(2*j-i)/4-2].imag();
        B0x = permute2<1,0>(B0x, B0x);
        B0y = permute2<1,0>(B0y, B0y);
        B1x = permute2<1,0>(B1x, B1x);
        B1y = permute2<1,0>(B1y, B1y);
        b_0246x = unpackhi(B0x, B1x);
        b_1357x = unpacklo(B0x, B1x);
        b_0246y = unpackhi(B0y, B1y);
        b_1357y = unpacklo(B0y, B1y);

        mul_conj(cc_0123x,cc_0123y, a_0246x,a_0246y, b_0246x,b_0246y, recpN);
        mul_conj(dd_0123x,dd_0123y, a_1357x,a_1357y, b_1357x,b_1357y, recpN);        

        ulong indA = (j+i)/2;   // 0 mod 4
        ulong indB = (2*j-i)/2; // 0 mod 4
        qa0123x = packed<double,4>(tab1[indA+0].real(), tab1[indA+2].real(), tab1[indA+1].real(), tab1[indA+3].real());
        qa0123y = packed<double,4>(tab1[indA+0].imag(), tab1[indA+2].imag(), tab1[indA+1].imag(), tab1[indA+3].imag());
        qb1234x = packed<double,4>(tab2[indB-1].real(), tab2[indB-3].real(), tab2[indB-2].real(), tab2[indB-4].real());
        qb1234y = packed<double,4>(tab2[indB-1].imag(), tab2[indB-3].imag(), tab2[indB-2].imag(), tab2[indB-4].imag());

        s_0123x = add(cc_0123x, dd_0123x);
        s_0123y = add(cc_0123y, dd_0123y);
        m_0123x = sub(cc_0123x, dd_0123x);
        m_0123y = sub(cc_0123y, dd_0123y);

        r1_0123x = fnmadd(qa0123x, m_0123y,  fmadd(qa0123y, m_0123x, s_0123x));
        r1_0123y =  fmadd(qa0123x, m_0123x,  fmadd(qa0123y, m_0123y, s_0123y));
        r2_3210x = fnmadd(qb1234x, m_0123y, fnmadd(qb1234y, m_0123x, s_0123x));
        r2_3210y = fnmadd(qb1234x, m_0123x,  fmsub(qb1234y, m_0123y, s_0123y));

        Z[(0*l/2+indA)/4].real(permute<0,2,1,3>(r1_0123x));
        Z[(0*l/2+indA)/4].imag(permute<0,2,1,3>(r1_0123y));
        Z[(1*l/2+indB-4)/4].real(permute<3,1,2,0>(r2_3210x));
        Z[(1*l/2+indB-4)/4].imag(permute<3,1,2,0>(r2_3210y));
    }
#undef f1
#undef f2
}



/*
tab0[i] = cispi_frac_ui(2*n_revbin(i,k-1),2^k) for 1 <= i < 2^(k-1)
        = ?, 1/4, 1/8, 3/8, 1/16, 5/16, ...
tab1[i]: add 1/(3*2^k)
tab2[i]: sub 1/(3*2^k)
*/
static inline void _three_recover_and_mul_and_pack(
    complex_packed<double,4>* A,
    ulong l,
    const std::complex<double>* tab0,
    const std::complex<double>* tab1,
    const std::complex<double>* tab2)
{
    double fac = 1.0/slong(3*l);
    _recover_and_mul_and_pack(A, l, tab0, fac);
    _recover_and_mul_and_pack(A+l/8, A+l/4, A+l/4*2, l, fac, tab1, tab2);
}

static inline void _one_recover_and_mul_and_pack(
    complex_packed<double,4>* A,
    ulong l,
    const std::complex<double>* tab0)
{
    double fac = 1.0/slong(l);
    _recover_and_mul_and_pack(A, l, tab0, fac);
}

static inline void _five_recover_and_mul_and_pack(
    complex_packed<double,4>* A,
    ulong l,
    const std::complex<double>* tab0,
    const std::complex<double>* tab1,
    const std::complex<double>* tab4,
    const std::complex<double>* tab2,
    const std::complex<double>* tab3)
{
    double fac = 1.0/slong(5*l);
    _recover_and_mul_and_pack(A+l/8*0, l, tab0, fac);
    _recover_and_mul_and_pack(A+l/8*1, A+l/4, A+l/4*2, l, fac, tab1, tab4);
    _recover_and_mul_and_pack(A+l/8*3, A+l/4*3, A+l/4*4, l, fac, tab2, tab3);
}


double cpd_fft_ctx::my_mpn_mul(
    ulong* z,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    ulong bits,
    ulong m)
{
    ulong zn = an + bn;
    ulong alen = cdiv(64*an, bits);
    ulong blen = cdiv(64*bn, bits);
    ulong zlen = alen + blen - 1;
    ulong atrunc = round_up(alen, 32);
    ulong btrunc = round_up(blen, 32);
    ulong ztrunc = round_up(zlen, 32);
    ulong abtrunc = std::max(atrunc, btrunc);
    ulong k = std::max(ulong(5), clog2(cdiv(ztrunc, m)));
    ulong l = pow2(k);
    fit_wtab(k);

    // temp space
    complex_packed<double,4>* A = complex_packed_double_buffer(m*l/4);

    mpn_to_cpd(&A[0].re, abtrunc, a, an, bits);
    mpn_to_cpd(&A[0].im, abtrunc, b, bn, bits);

    if (m == 5)
    {
        fft_five_trunc(A, k, abtrunc);
        _five_recover_and_mul_and_pack(A, l, w2rev_table(), wfiverev1(k), wfiverev4(k), wfiverev2(k), wfiverev3(k));
        ifft_five_full(A, k-1);
    }
    else if (m == 3)
    {
        fft_three_trunc(A, k, abtrunc);
        _three_recover_and_mul_and_pack(A, l, w2rev_table(), wthreerev1(k), wthreerev2(k));
        ifft_three_full(A, k-1);
    }
    else if (m == 1)
    {
        fft_trunc(A, k, abtrunc);
        _one_recover_and_mul_and_pack(A, l, w2rev_table());
        ifft_full(A, k-1);
    }
    else
    {
        std::cout << "unsupported multiplier " << m << "in cpd_fft_ctx::my_mpn_mul" << std::endl;
        std::abort();
    }

    double max_error = 0;
    for (ulong i = 0; i < ztrunc/8; i++)
    {
        double ai = A[i/4].real()[i%4];
        max_error = std::max(max_error, std::abs(ai - round(ai)));
        ai = A[i/4].imag()[i%4];
        max_error = std::max(max_error, std::abs(ai - round(ai)));
    }

    mpn_from_cpd(z, zn, A, ztrunc, bits);

//std::cout << "max_error: " << max_error << std::endl;

    return max_error;
}

void cpd_fft_ctx::my_mpn_mul(
    ulong* z,
    const ulong* a, ulong an,
    const ulong* b, ulong bn)
{
    ulong zn = an + bn;
    ulong bits;

         if (zn <      30) bits = 22;
    else if (zn <     100) bits = 21;
    else if (zn <     300) bits = 20;
    else if (zn <    1000) bits = 19;
    else if (zn <    3500) bits = 18;
    else if (zn <   10000) bits = 17;
    else if (zn <   35000) bits = 16;
    else if (zn <   90000) bits = 15;
    else if (zn <  200000) bits = 14;
    else if (zn <  400000) bits = 13;
    else if (zn <  800000) bits = 12;
    else if (zn < 1600000) bits = 11;
    else if (zn < 3200000) bits = 10;
    else if (zn < 6400000) bits = 9;
    else                   bits = 8;

    ulong alen = cdiv(64*an, bits);
    ulong blen = cdiv(64*bn, bits);
    ulong zlen = alen + blen - 1;
    ulong atrunc = round_up(alen, 32);
    ulong btrunc = round_up(blen, 32);
    ulong ztrunc = round_up(zlen, 32);
    ulong abtrunc = std::max(atrunc, btrunc);

    // TODO choice of m really should be in a LUT along with choice of bits
    ulong m = 1;
    ulong best_score = (ztrunc-1) << leading_zeros(ztrunc-1);
    for (ulong try_m = 3; try_m <= 5; try_m += 2)
    {
        ulong t = cdiv(ztrunc, try_m) - 1;
        t = t << leading_zeros(t);
        if (t > best_score)
        {
            m = try_m;
            best_score = t;
        }
    }

    ulong k = std::max(ulong(5), clog2(cdiv(ztrunc, m)));
    ulong l = pow2(k);

    complex_packed<double,4>* A = complex_packed_double_buffer(m*l/4);

    fit_wtab(k);

    mpn_to_cpd(&A[0].re, abtrunc, a, an, bits);
    mpn_to_cpd(&A[0].im, abtrunc, b, bn, bits);

    if (m == 5)
    {
        fft_five_trunc(A, k, abtrunc);
        _five_recover_and_mul_and_pack(A, l, w2rev_table(), wfiverev1(k), wfiverev4(k), wfiverev2(k), wfiverev3(k));
        ifft_five_full(A, k-1);
    }
    else if (m == 3)
    {
        fft_three_trunc(A, k, abtrunc);
        _three_recover_and_mul_and_pack(A, l, w2rev_table(), wthreerev1(k), wthreerev2(k));
        ifft_three_full(A, k-1);
    }
    else
    {
        fft_trunc(A, k, abtrunc);
        _one_recover_and_mul_and_pack(A, l, w2rev_table());
        ifft_full(A, k-1);
    }

    mpn_from_cpd(z, zn, A, ztrunc, bits);
}

