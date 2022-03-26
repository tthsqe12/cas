
template<typename T, int SZ>
struct packed {
    packed() = delete;
};

template<>
struct packed<ulong, 4> {
    __m256i data;
    packed() {}
    packed(uint64_t x) : data(_mm256_set1_epi64x(x)) {}
    packed(const uint64_t* a) : data(_mm256_loadu_si256((const __m256i*) a)) {}
    packed(__m256i _data) : data(_data) {}
    packed(uint64_t d0, uint64_t d1, uint64_t d2, uint64_t d3) : data(_mm256_set_epi64x(d3, d2, d1, d0)) { }
    uint64_t operator[](int i) {return data[i];}
    void load(const uint64_t* a) {data = _mm256_loadu_si256((const __m256i*) a);}
    void load_aligned(const uint64_t* a) {data = _mm256_load_si256((const __m256i*) a);}
    void store(uint64_t* a) const {_mm256_storeu_si256((__m256i*) a, data);}
    void store_aligned(uint64_t* a) const {_mm256_store_si256((__m256i*) a, data);}
};

template<>
struct packed<slong, 4> {
    __m256i data;
    packed() {}
    packed(int64_t x) : data(_mm256_set1_epi64x(x)) {}
    packed(__m256i _data) : data(_data) {}
    packed(int64_t d0, int64_t d1, int64_t d2, int64_t d3) : data(_mm256_set_epi64x(d3, d2, d1, d0)) { }
};

template<>
struct packed<double, 4> {
    __m256d data;
    packed() {}
    packed(double x) : data(_mm256_set1_pd(x)) {}
    packed(__m256d _data) : data(_data) {}
    packed(double d0, double d1, double d2, double d3) : data(_mm256_set_pd(d3, d2, d1, d0)) {}
    double operator[](int i) {return data[i];}
    void zero() {data = _mm256_setzero_pd();}
    void load(const double* a) {data = _mm256_loadu_pd(a);}
    void load_aligned(const double* a) {data = _mm256_load_pd(a);}
    void store(double* a) const {_mm256_storeu_pd(a, data);}
    void store_aligned(double* a) const {_mm256_store_pd(a, data);}
};

#define PI4 packed<slong, 4>
#define PU4 packed<ulong, 4>
#define PD4 packed<double, 4>

std::ostream& operator<<(std::ostream& o, const PD4& x)
{
    double y[4];
    x.store(y);
    int64_t z[4];
    bool ok = false;
    for (int v = 0; v < 4; v++)
    {
        z[v] = y[v];
        y[v] -= z[v];
        if (y[v] != 0)
            ok = true;
    }
    o << "{" << z[0] << ", " << z[1] << ", " << z[2] << ", " << z[3] << "}";
    if (ok)
        o << " + <" << y[0]-z[0] << ", " << y[1]-z[1] << ", " << y[2] - z[2] << ", " << y[3] - z[3] << ">";
    return o;
}

inline double blendv(double a, double b, double c) {
    return c >= 0 ? a : b;
}

inline PD4 blendv(PD4 a, PD4 b, PD4 c) {
    return _mm256_blendv_pd(a.data, b.data, c.data);
}

inline double round(double a) {
    return std::rint(a);
}

inline PD4 round(PD4 a) {
    return _mm256_round_pd(a.data, 4);
}

inline double neg(double a) {
    return -a;
}

inline PD4 neg(PD4 a) {
    return _mm256_sub_pd(_mm256_setzero_pd(), a.data);
}

inline double add(double a, double b) {
    return a + b;
}

inline PD4 add(PD4 a, PD4 b) {
    return _mm256_add_pd(a.data, b.data);
}

inline double sub(double a, double b) {
    return a - b;
}

inline PD4 sub(PD4 a, PD4 b) {
    return _mm256_sub_pd(a.data, b.data);
}

inline double mul(double a, double b) {
    return a*b;
}

inline PD4 mul(PD4 a, PD4 b) {
    return _mm256_mul_pd(a.data, b.data);
}

inline double div(double a, double b) {
    return a/b;
}

inline PD4 div(PD4 a, PD4 b) {
    return _mm256_div_pd(a.data, b.data);
}

inline double fmadd(double a, double b, double c) {
    return std::fma(a, b, c);
}

inline PD4 fmadd(PD4 a, PD4 b, PD4 c) {
    return _mm256_fmadd_pd(a.data, b.data, c.data);
}

inline double fmsub(double a, double b, double c) {
    return std::fma(a, b, -c);
}

inline PD4 fmsub(PD4 a, PD4 b, PD4 c) {
    return _mm256_fmsub_pd(a.data, b.data, c.data);
}

inline double fnmadd(double a, double b, double c) {
    return std::fma(-a, b, c);
}

inline PD4 fnmadd(PD4 a, PD4 b, PD4 c) {
    return _mm256_fnmadd_pd(a.data, b.data, c.data);
}

// inputs in the range [0, 2^52)
template <typename T> inline T convert_limited(PU4 a) {}
template <typename T> inline T convert_limited(PD4 a) {}

template <>
inline PD4 convert_limited<PD4>(PU4 a) {
    __m256d t = _mm256_set1_pd(0x0010000000000000);
    __m256i x = _mm256_or_si256(a.data, _mm256_castpd_si256(t));
    return _mm256_sub_pd(_mm256_castsi256_pd(x), t);
}

template <>
inline PU4 convert_limited<PU4>(PD4 a) {
    __m256d t = _mm256_set1_pd(0x0010000000000000);
    __m256d x = _mm256_add_pd(a.data, t);
    return _mm256_xor_si256(_mm256_castpd_si256(x), _mm256_castpd_si256(t));
}

inline PD4 operator>=(PD4& a, PD4 b) {
   return _mm256_cmp_pd(a.data, b.data, _CMP_GE_OQ);
}

inline PD4 operator>(PD4& a, PD4 b) {
   return _mm256_cmp_pd(a.data, b.data, _CMP_GT_OQ);
}

inline PD4 operator<=(PD4& a, PD4 b) {
   return _mm256_cmp_pd(a.data, b.data, _CMP_LE_OQ);
}

inline PD4 operator<(PD4& a, PD4 b) {
   return _mm256_cmp_pd(a.data, b.data, _CMP_LT_OQ);
}

// return in [-n,n] assuming input in [-2*n,2*n]
inline PD4 reduce_pm2n_to_pm1n(PD4 a, PD4 n, PD4 minus_n) {
    return sub(a, blendv(n, minus_n, a));
}

inline PD4 reduce_pm2n_to_pm1n(PD4 a, PD4 n) {
    return reduce_pm2n_to_pm1n(a, n, neg(n));
}

// return in [0,n) assuming input in [-n,n]
inline PD4 reduce_pm1n_to_0n(PD4 a, PD4 n, PD4 minus_n) {
    PD4 t = sub(a, blendv(n, minus_n, a));
    return blendv(t, a, t);
}

// return in [0,n) assuming input in (-n,n)
inline PD4 reduce_pm1no_to_0n(PD4 a, PD4 n) {
    return blendv(a, add(a, n), a);
}

// return in [0,n) assuming input in [-2*n,2*n]
inline PD4 reduce_pm2n_to_0n(PD4 a, PD4 n, PD4 minus_n) {
    a = reduce_pm2n_to_pm1n(a, n, minus_n);
    PD4 t = sub(a, blendv(n, minus_n, a));
    return blendv(t, a, t);
}

inline PD4 reduce_pm2n_to_0n(PD4 a, PD4 n) {
    return reduce_pm2n_to_0n(a, n, neg(n));
}

// return in [-n,n] assuming input in [-2*n,2*n]
inline double reduce_pm2n_to_pm1n(double a, double n, double minus_n) {
    return sub(a, blendv(n, minus_n, a));
}

inline double reduce_pm2n_to_pm1n(double a, double n) {
    return reduce_pm2n_to_pm1n(a, n, neg(n));
}

// return in [0,n) assuming input in [-n,n]
inline double reduce_pm1n_to_0n(double a, double n, double minus_n) {
    double t = sub(a, blendv(n, minus_n, a));
    return blendv(t, a, t);
}

inline double reduce_pm1n_to_0n(double a, double n) {
    return reduce_pm1n_to_0n(a, n, neg(n));
}

// return in [0,n) assuming input in (-n,n)
inline double reduce_pm1no_to_0n(double a, double n) {
    return blendv(a, add(a, n), a);
}

// return in [0,n) assuming input in [-2*n,2*n]
inline double reduce_pm2n_to_0n(double a, double n, double minus_n) {
    a = reduce_pm2n_to_pm1n(a, n, minus_n);
    double t = sub(a, blendv(n, minus_n, a));
    return blendv(t, a, t);
}

inline double reduce_pm2n_to_0n(double a, double n) {
    return reduce_pm2n_to_0n(a, n, neg(n));
}

// reduce a mod n in [-n,n] assuming
inline double reduce_to_pm1n(double a, double n, double ninv)
{
    double q = round(mul(a, ninv));
    return fnmadd(q, n, a);
}

inline PD4 reduce_to_pm1n(PD4 a, PD4 n, PD4 ninv)
{
    PD4 q = round(mul(a, ninv));
    return fnmadd(q, n, a);
}

// return a*b mod n in [-n,n] assuming
// a in [-2*n, 2*n], b in (-n/2, n/2), bon = b/n
// and **all rounds are to some nearest**
// 0 < n < 2^52
inline double mul_mod(double a, double b, double bon, double n) {
    double h = mul(a, b);
    double q = round(mul(a, bon));
    double l = fmsub(a, b, h);
    return add(fnmadd(q, n, h), l);
}

inline PD4 mul_mod(PD4 a, PD4 b, PD4 bon, PD4 n) {
    PD4 h = mul(a, b);
    PD4 q = round(mul(a, bon));
    PD4 l = fmsub(a, b, h);
    return add(fnmadd(q, n, h), l);
}

// return a*b mod n in [-n,n] assuming
// a in [-4*n, 4*n], b in (-n/2, n/2), ninv = 1/n
// and **all rounds are to some nearest**
// 0 < n < 2^50
inline double mulmod2(double a, double b, double n, double ninv) {
    double h = mul(a, b);
    double q = round(mul(h, ninv));
    double l = fmsub(a, b, h);
    return add(fnmadd(q, n, h), l);
}

inline PD4 mulmod2(PD4 a, PD4 b, PD4 n, PD4 ninv) {
    PD4 h = mul(a, b);
    PD4 q = round(mul(h, ninv));
    PD4 l = fmsub(a, b, h);
    return add(fnmadd(q, n, h), l);
}


#define add_sssssaaaaaaaaaa(s4,s3,s2,s1,s0, a4,a3,a2,a1,a0, b4,b3,b2,b1,b0)  \
  __asm__ ("addq %14,%q4\n\tadcq %12,%q3\n\tadcq %10,%q2\n\tadcq %8,%q1\n\tadcq %6,%q0"    \
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "1"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "2"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "3"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "4"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)  \
  __asm__ ("subq %11,%q3\n\tsbbq %9,%q2\n\tsbbq %7,%q1\n\tsbbq %5,%q0"    \
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "1"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "2"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "3"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))


#define _SHIFT_LEFT4(a, k)                  \
t3 = ((a)[3] << (k)) | ((a)[2] >> (64-k));  \
t2 = ((a)[2] << (k)) | ((a)[1] >> (64-k));  \
t1 = ((a)[1] << (k)) | ((a)[0] >> (64-k));  \
t0 = (a)[0] << (k);

#define _SHIFT_LEFT5(a, k)                      \
t4 =                   ((a)[3] >> (64-(k)));    \
t3 = ((a)[3] << (k)) | ((a)[2] >> (64-(k)));    \
t2 = ((a)[2] << (k)) | ((a)[1] >> (64-(k)));    \
t1 = ((a)[1] << (k)) | ((a)[0] >> (64-(k)));    \
t0 =  (a)[0] << (k);


struct format_hex {
    ulong data;
    format_hex(ulong _data) : data(_data) {}
};

std::ostream& operator<<(std::ostream& o, const format_hex& a)
{
    for (int i = 15; i >= 0; i--)
    {
        ulong b = (a.data >> (4*i)) & 0x0f;
        o << char(b > 9 ? 'a' + b - 10 : '0' + b);
    }
    return o;
}


#undef PI4
#undef PU4
#undef PD4
