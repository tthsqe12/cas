#pragma once

#include <cstdlib>
#include <cassert>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
#include <x86intrin.h>
#include <immintrin.h>

#include "arb.h"

#define FORCE_INLINE __attribute__((always_inline))
#define UNLIKELY(x) __builtin_expect((x),0)
#define LIKELY(x)   __builtin_expect((x),1)


#define PI4 packed<slong,4>
#define PU4 packed<ulong,4>
#define PD4 packed<double,4>


template<typename T, int SZ>
struct packed {
    packed() = delete;

    template <int i0, int i1> inline void store_subset(T* a) const;
};

template<>
struct packed<uint32_t,4> {
    __m128i data;
    packed() {}
    packed(__m128i _data) : data(_data) {}
    packed(uint32_t x) : data(_mm_set1_epi32(x)) {}
    packed(const uint32_t* a) : data(_mm_loadu_si128((const __m128i*) a)) {}
    packed(uint32_t d0, uint32_t d1, uint32_t d2, uint32_t d3) : data(_mm_set_epi32(d3, d2, d1, d0)) { }
    uint32_t operator[](int i) const {
        if (i == 0) return _mm_extract_epi32(data, 0);
        if (i == 1) return _mm_extract_epi32(data, 1);
        if (i == 2) return _mm_extract_epi32(data, 2);
        if (i == 3) return _mm_extract_epi32(data, 3);
        std::cout << "oops" << std::endl;
        std::abort();
    }
    void load(const uint32_t* a) {data = _mm_loadu_si128((const __m128i*) a);}
    void load_aligned(const uint32_t* a) {data = _mm_load_si128((const __m128i*) a);}
    void store(uint32_t* a) const {_mm_storeu_si128((__m128i*) a, data);}
    void store_aligned(uint32_t* a) const {_mm_store_si128((__m128i*) a, data);}
};

template<>
struct packed<uint32_t,8> {
    __m256i data;
    packed() {}
    packed(__m256i _data) : data(_data) {}
    packed(uint32_t x) : data(_mm256_set1_epi32(x)) {}
    packed(const uint32_t* a) : data(_mm256_loadu_si256((const __m256i*) a)) {}
    packed(uint32_t d0, uint32_t d1, uint32_t d2, uint32_t d3, uint32_t d4, uint32_t d5, uint32_t d6, uint32_t d7) : data(_mm256_set_epi32(d7, d6, d5, d4, d3, d2, d1, d0)) { }

    uint32_t operator[](int i) const {
        if (i == 0) return _mm256_extract_epi32(data, 0);
        if (i == 1) return _mm256_extract_epi32(data, 1);
        if (i == 2) return _mm256_extract_epi32(data, 2);
        if (i == 3) return _mm256_extract_epi32(data, 3);
        if (i == 4) return _mm256_extract_epi32(data, 4);
        if (i == 5) return _mm256_extract_epi32(data, 5);
        if (i == 6) return _mm256_extract_epi32(data, 6);
        if (i == 7) return _mm256_extract_epi32(data, 7);
        std::cout << "oops" << std::endl;
        std::abort();
    }

    void load(const uint32_t* a) {data = _mm256_loadu_si256((const __m256i*) a);}
    void load_aligned(const uint64_t* a) {data = _mm256_load_si256((const __m256i*) a);}
    void store(uint32_t* a) const {_mm256_storeu_si256((__m256i*) a, data);}
    void store_aligned(uint32_t* a) const {_mm256_store_si256((__m256i*) a, data);}

#ifdef NEW_ALIGNMENT_IS_BUGGY
    void* operator new[] (size_t n) {
        void* p;
        if (posix_memalign(&p, alignof(__m256i), n))
            throw std::bad_alloc();
        return p;
    }
#endif
};

template<>
struct PU4 {
    __m256i data;
    packed() {}
    packed(__m256i _data) : data(_data) {}
    packed(uint64_t x) : data(_mm256_set1_epi64x(x)) {}
    packed(const uint64_t* a) : data(_mm256_loadu_si256((const __m256i*) a)) {}
    packed(uint64_t d0, uint64_t d1, uint64_t d2, uint64_t d3) : data(_mm256_set_epi64x(d3, d2, d1, d0)) { }
    uint64_t operator[](int i) const {
        if (i == 0) return _mm256_extract_epi64(data, 0);
        if (i == 1) return _mm256_extract_epi64(data, 1);
        if (i == 2) return _mm256_extract_epi64(data, 2);
        if (i == 3) return _mm256_extract_epi64(data, 3);
        std::cout << "oops" << std::endl;
        std::abort();        
    }

    void load(const uint64_t* a) {data = _mm256_loadu_si256((const __m256i*) a);}
    void load_aligned(const uint64_t* a) {data = _mm256_load_si256((const __m256i*) a);}
    void store(uint64_t* a) const {_mm256_storeu_si256((__m256i*) a, data);}
    void store_aligned(uint64_t* a) const {_mm256_store_si256((__m256i*) a, data);}

#ifdef NEW_ALIGNMENT_IS_BUGGY
    void* operator new[] (size_t n) {
        void* p;
        if (posix_memalign(&p, alignof(__m256i), n))
            throw std::bad_alloc();
        return p;
    }
#endif
};

template<>
struct PI4 {
    __m256i data;
    packed() {}
    packed(__m256i _data) : data(_data) {}
    packed(int64_t x) : data(_mm256_set1_epi64x(x)) {}
    packed(int64_t d0, int64_t d1, int64_t d2, int64_t d3) : data(_mm256_set_epi64x(d3, d2, d1, d0)) { }

#ifdef NEW_ALIGNMENT_IS_BUGGY
    void* operator new[] (size_t n) {
        void* p;
        if (posix_memalign(&p, alignof(__m256i), n))
            throw std::bad_alloc();
        return p;
    }
#endif
};

template<>
struct PD4 {
    __m256d data;
    packed() {}
    packed(__m256d _data) : data(_data) {}
    packed(double x) : data(_mm256_set1_pd(x)) {}
    packed(double d0, double d1, double d2, double d3) : data(_mm256_set_pd(d3, d2, d1, d0)) {}
    double operator[](int i) const {return data[i];}
    void zero() {data = _mm256_setzero_pd();}
    // unspecified
    void load(const double* a) {data = _mm256_load_pd(a);}
    void store(double* a) const {_mm256_store_pd(a, data);}
    // specified
    void load_unaligned(const double* a) {data = _mm256_loadu_pd(a);}
    void store_unaligned(double* a) const {_mm256_storeu_pd(a, data);}
    void load_aligned(const double* a) {data = _mm256_load_pd(a);}
    void store_aligned(double* a) const {_mm256_store_pd(a, data);}
    template <int i0, int i1> inline void store_subset(double* a) const;

#ifdef NEW_ALIGNMENT_IS_BUGGY
    void* operator new[] (size_t n) {
        void* p;
        if (posix_memalign(&p, alignof(__m256d), n))
            throw std::bad_alloc();
        return p;
    }
#endif
};


template <typename T, int SZ>
struct complex_packed {
    packed<T,SZ> re;
    packed<T,SZ> im;

    complex_packed() {}
    complex_packed(std::complex<T> a) : re(a.real()), im(a.imag()) {}
    complex_packed(packed<T,SZ> x, packed<T,SZ> y) : re(x), im(y) {}

    packed<T,SZ> real() const {return re;}
    packed<T,SZ> imag() const {return im;}
    void real(packed<T,SZ> x) {re = x;}
    void imag(packed<T,SZ> y) {im = y;}

#ifdef NEW_ALIGNMENT_IS_BUGGY
    void* operator new[] (size_t n) {
        void* p;
        if (posix_memalign(&p, alignof(packed<T,SZ>), n))
            throw std::bad_alloc();
        return p;
    }
#endif
};

#ifdef NEW_ALIGNMENT_IS_BUGGY
template<typename T>
struct my_vector {
    T* _data;
    size_t _length;
    size_t _alloc;
    my_vector() : _data(nullptr), _length(0), _alloc(0) {};
    ~my_vector() {free(_data);}

    my_vector(const my_vector<T>& other) = delete;

    my_vector(my_vector<T>&& other) {
        _data = other._data;
        _length = other._length;
        _alloc = other._alloc;
        other._data = nullptr;
        other._length = 0;
        other._alloc = 0;
    }

    T& operator [](size_t i) {return _data[i];}
    T& at(size_t i) {assert(i < _length); return _data[i];}
    T* data() {return _data;}
    size_t size() {return _length;}
    void clear() {_length = 0;}
    void resize(size_t n) {
        if (n > _alloc) {
            size_t newalloc = std::max(n, _alloc + _alloc/2);
            T* newdata = new T[newalloc];
            for (ulong i = 0; i < std::min(_length, n); i++)
                newdata[i] = _data[i];
            free(_data);
            _data = newdata;
            _alloc = newalloc;
        }
        _length = n;
    }
};
#else
#define my_vector std::vector
#endif

#ifdef NEW_ALIGNMENT_IS_BUGGY
inline void* my_aligned_alloc(size_t alignment, size_t n)
{
    void* p;
    if (posix_memalign(&p, alignment, n))
        return nullptr;
    return p;
}

#else
#define my_aligned_alloc std::aligned_alloc
#endif

template <typename T, int SZ>
std::ostream& operator<<(std::ostream& o, const complex_packed<T, SZ>& x)
{
    o << x.real() << " + I*" << x.imag();
    return o;
}


template <typename T, int SZ>
std::ostream& operator<<(std::ostream& o, const packed<T, SZ>& x)
{
    std::cout << "{";
    for (int i = 0; i < SZ; i++)
    {
        if (i > 0)
            std::cout << ", ";
        std::cout << x[i];
    }
    std::cout << "}";
    return o;
}


/* movement */
inline double blendv(double a, double b, bool c) {
    return c ? b : a;
}

inline double blendv(double a, double b, double c) {
    return c >= 0 ? a : b;
}

inline PD4 blendv(PD4 a, PD4 b, PD4 c) {
    return _mm256_blendv_pd(a.data, b.data, c.data);
}

template <int i0, int i1, int i2, int i3>
inline PD4 blend(PD4 a, PD4 b) {
    return _mm256_blend_pd(a.data, b.data, i0 + 2*(i1 + 2*(i2 + 2*i3)));
}

// return {a[0], b[0], a[2], b[2]}
inline PD4 unpacklo(PD4 a, PD4 b) {
    return _mm256_unpacklo_pd(a.data, b.data);
}

// return {a[1], b[1], a[3], b[3]}
inline PD4 unpackhi(PD4 a, PD4 b) {
    return _mm256_unpackhi_pd(a.data, b.data);
}

inline PD4 insertf128(PD4 a, PD4 b, int i) {
    return _mm256_insertf128_pd(a.data, _mm256_castpd256_pd128(b.data), i);
}

// return {v[i0], v[i1]}
// |   v[0]     |    v[1]    |    v[2]    |   v[3]     |
//   a[0], a[1]   a[2], a[3]   b[0], b[1]   b[2], b[3]
template <int i0, int i1>
inline PD4 permute2(PD4 a, PD4 b) {
    return _mm256_permute2f128_pd(a.data, b.data, i0 + 16*i1);
}

// return {a[i0], a[i1], a[i2], a[i3]}
template <int i0, int i1, int i2, int i3>
inline PD4 permute(PD4 a) {
#ifndef AVOID_AVX2
    return _mm256_permute4x64_pd(a.data, i0 + 4*(i1 + 4*(i2 + 4*i3)));
#else
    return PD4(a.data[i0], a.data[i1], a.data[i2], a.data[i3]);
#endif
}

// extract entries at indices i_0, i_1, .. i_{m_1}
// and store them at m consecutive locations
template<>
inline void packed<double,4>::store_subset<0,1>(double* a) const {
    _mm_storeu_pd(a, _mm256_castpd256_pd128(data));
}

template<>
inline void packed<double,4>::store_subset<2,3>(double* a) const {
    _mm_storeu_pd(a, _mm256_extractf128_pd(data, 1));
}



FORCE_INLINE inline void transpose_4x4(
    packed<double,4>& r0,
    packed<double,4>& r1,
    packed<double,4>& r2,
    packed<double,4>& r3,
    packed<double,4> x0,
    packed<double,4> x1,
    packed<double,4> x2,
    packed<double,4> x3)
{
    packed<double,4> t0, t1, t2, t3;
    t0 = unpacklo(x0, x1);
    t1 = unpackhi(x0, x1);
    t2 = unpacklo(x2, x3);
    t3 = unpackhi(x2, x3);
    r0 = permute2<0,2>(t0, t2);
    r1 = permute2<0,2>(t1, t3);
    r2 = permute2<1,3>(t0, t2);
    r3 = permute2<1,3>(t1, t3);
}

FORCE_INLINE inline void transpose_4x4(
    complex_packed<double,4>& r0,
    complex_packed<double,4>& r1,
    complex_packed<double,4>& r2,
    complex_packed<double,4>& r3,
    complex_packed<double,4> x0,
    complex_packed<double,4> x1,
    complex_packed<double,4> x2,
    complex_packed<double,4> x3)
{
    packed<double,4> t0x, t1x, t2x, t3x, t0y, t1y, t2y, t3y;
    transpose_4x4(t0x, t1x, t2x, t3x, x0.real(), x1.real(), x2.real(), x3.real());
    transpose_4x4(t0y, t1y, t2y, t3y, x0.imag(), x1.imag(), x2.imag(), x3.imag());
    r0.real(t0x);
    r0.imag(t0y);
    r1.real(t1x);
    r1.imag(t1y);
    r2.real(t2x);
    r2.imag(t2y);
    r3.real(t3x);
    r3.imag(t3y);
}

/* bit twiddling */

inline PU4 bit_and(PU4 a, PU4 b) {
#ifndef AVOID_AVX2
    return _mm256_and_si256(a.data, b.data);
#else
    return PU4(a[0]&b[0], a[1]&b[1], a[2]&b[2], a[3]&b[3]);
#endif
}

inline PU4 bit_shift_right(PU4 a, packed<uint32_t,4> b) {
    return _mm256_srl_epi64(a.data, b.data);
}

inline PU4 bit_shift_right(PU4 a, PU4 b) {
#ifndef AVOID_AVX2
    return _mm256_srlv_epi64(a.data, b.data);
#else
    return PU4(a[0]>>b[0], a[1]>>b[1], a[2]>>b[2], a[3]>>b[3]);
#endif
}

/* arithmetic */
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
    __m256d mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000));
    return _mm256_xor_pd(a.data, mask);
}

inline PD4 abs(PD4 a) {
    __m256d mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffff));
    return _mm256_and_pd(a.data, mask);
}

inline PD4 max(PD4 a, PD4 b) {
    return _mm256_max_pd(a.data, b.data);
}

inline PD4 min(PD4 a, PD4 b) {
    return _mm256_min_pd(a.data, b.data);
}

inline double add(double a, double b) {
    return a + b;
}

inline PD4 operator+(PD4 a, PD4 b) {
    return _mm256_add_pd(a.data, b.data);
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

inline PD4 addsub(PD4 a, PD4 b) {
    return _mm256_addsub_pd(a.data, b.data);
}

inline double mul(double a, double b) {
    return a*b;
}

inline PD4 mul(PD4 a, PD4 b) {
    return _mm256_mul_pd(a.data, b.data);
}
inline PD4 mul(double a, PD4 b) {
    return mul(PD4(a), b);
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
#ifndef AVOID_AVX2
    return _mm256_fmadd_pd(a.data, b.data, c.data);
#else
    return _mm256_macc_pd(a.data, b.data, c.data);
#endif
}

inline double fmsub(double a, double b, double c) {
    return std::fma(a, b, -c);
}

inline PD4 fmsub(PD4 a, PD4 b, PD4 c) {
#ifndef AVOID_AVX2
    return _mm256_fmsub_pd(a.data, b.data, c.data);
#else
    return _mm256_msub_pd(a.data, b.data, c.data);
#endif
}

inline double fnmadd(double a, double b, double c) {
    return std::fma(-a, b, c);
}

inline PD4 fnmadd(PD4 a, PD4 b, PD4 c) {
#ifndef AVOID_AVX2
    return _mm256_fnmadd_pd(a.data, b.data, c.data);
#else
    return _mm256_nmacc_pd(a.data, b.data, c.data);
#endif
}

inline PD4 fnmadd(double a, PD4 b, double c) {return fnmadd(PD4(a), b, PD4(c));}

inline PD4 half(PD4 b) {return mul(0.5, b);}

inline PD4 operator-(PD4 a) {return neg(a);}
inline PD4 operator-(PD4 a, PD4 b) {return sub(a, b);}
inline PD4 operator*(double a, PD4 b) {return mul(a, b);}
inline PD4 operator*(PD4 a, PD4 b) {return mul(a, b);}


// inputs in the range [0, 2^52)
template <typename T> inline T convert_limited(PU4 a) {}
template <typename T> inline T convert_limited(PD4 a) {}

template <>
inline PD4 convert_limited<PD4>(PU4 a) {
    __m256d t = _mm256_set1_pd(0x0010000000000000);
#ifndef AVOID_AVX2
    __m256i x = _mm256_or_si256(a.data, _mm256_castpd_si256(t));
    return _mm256_sub_pd(_mm256_castsi256_pd(x), t);
#else
    return _mm256_sub_pd(_mm256_or_pd(_mm256_castsi256_pd(a.data), t), t);

#endif
}

template <>
inline PU4 convert_limited<PU4>(PD4 a) {
    __m256d t = _mm256_set1_pd(0x0010000000000000);
    __m256d x = _mm256_add_pd(a.data, t);
#ifndef AVOID_AVX2
    return _mm256_xor_si256(_mm256_castpd_si256(x), _mm256_castpd_si256(t));
#else
    return _mm256_castpd_si256(_mm256_xor_pd(x, t));
#endif
}

inline bool operator==(PD4 a, PD4 b) {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

inline bool operator!=(PD4 a, PD4 b) {
    return a[0] != b[0] || a[1] != b[1] || a[2] != b[2] || a[3] != b[3];
}

inline PD4 operator>=(PD4 a, PD4 b) {
   return _mm256_cmp_pd(a.data, b.data, _CMP_GE_OQ);
}

inline PD4 operator>(PD4 a, PD4 b) {
   return _mm256_cmp_pd(a.data, b.data, _CMP_GT_OQ);
}

inline PD4 operator<=(PD4 a, PD4 b) {
   return _mm256_cmp_pd(a.data, b.data, _CMP_LE_OQ);
}

inline PD4 operator<(PD4 a, PD4 b) {
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

// return a mod n in [-n,n]
template <typename T>
inline T reduce_to_pm1n(T a, T n, T ninv) {
    return fnmadd(round(mul(a, ninv)), n, a);
}

// return a mod n in (-n,n)
template <typename T>
inline T reduce_to_pm1no(T a, T n, T ninv) {
    return fnmadd(round(mul(a, ninv)), n, a);
}

// return a mod n in [0,n)
template <typename T> inline T reduce_to_0n(T a, T n, T ninv) {
    return reduce_pm1no_to_0n(reduce_to_pm1no(a, n, ninv), n);
}

// [0,n] -> [-n/2, n/2]
inline double reduce_0n_to_pmhn(double a, double n) {
    return a > n*0.5 ? a-n : a;
}

template <typename T>
inline T reduce_0n_to_pmhn(T a, T n, T half_n) {
    return blendv(a, sub(a, n), a > half_n);
}


// [-n,n] -> [-n/2, n/2]
inline double reduce_pm1n_to_pmhn(double a, double n) {
    if (a > 0.5*n)
        return a - n;
    else if (a < -0.5*n)
        return a + n;
    else
        return a;
}

template <typename T>
inline T reduce_pm1n_to_pmhn(T a, T n, T half_n) {
    T t = blendv(n, neg(n), a);
    a = blendv(a, sub(a, t), abs(a) > half_n);
    return a;

//    a = blendv(a, add(a,n), a);
//    a = blendv(sub(a,n), a, sub(a, half_n));
//    return a;
}

template <typename T>
inline T reduce_pm1n_to_pmhn(T a, T n) {
    a = blendv(a, add(a, n), a);
    T t = sub(a, n);
    a =  blendv(t, a, add(t, a));
    return a;
}


template <typename T>
inline bool equal_mod(T a, T b, T n, T ninv) {
    return reduce_to_0n(a, n, ninv) == reduce_to_0n(b, n, ninv);
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


#undef PI4
#undef PU4
#undef PD4


#define add_sssssaaaaaaaaaa(s4,s3,s2,s1,s0, a4,a3,a2,a1,a0, b4,b3,b2,b1,b0)  \
  __asm__ ("addq %14,%q4\n\tadcq %12,%q3\n\tadcq %10,%q2\n\tadcq %8,%q1\n\tadcq %6,%q0"    \
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "1"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "2"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "3"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "4"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define add_ssssssaaaaaaaaaaaa(s5,s4,s3,s2,s1,s0, a5,a4,a3,a2,a1,a0, b5,b4,b3,b2,b1,b0)  \
  __asm__ ("addq %17,%q5\nadcq %15,%q4\n\tadcq %13,%q3\n\tadcq %11,%q2\n\tadcq %9,%q1\n\tadcq %7,%q0"    \
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "1"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "2"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "3"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "4"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "5"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define add_sssssssaaaaaaaaaaaaaa(s6,s5,s4,s3,s2,s1,s0, a6,a5,a4,a3,a2,a1,a0, b6,b5,b4,b3,b2,b1,b0)  \
  __asm__ ("addq %20,%q6\nadcq %18,%q5\nadcq %16,%q4\n\tadcq %14,%q3\n\tadcq %12,%q2\n\tadcq %10,%q1\n\tadcq %8,%q0"    \
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a6)), "rme" ((mp_limb_t)(b6)),                 \
         "1"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "2"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "3"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "4"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "5"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "6"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))


#define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)  \
  __asm__ ("subq %11,%q3\n\tsbbq %9,%q2\n\tsbbq %7,%q1\n\tsbbq %5,%q0"    \
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "1"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "2"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "3"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_dddddmmmmmsssss(s4,s3,s2,s1,s0, a4,a3,a2,a1,a0, b4,b3,b2,b1,b0)  \
  __asm__ ("subq %14,%q4\n\tsbbq %12,%q3\n\tsbbq %10,%q2\n\tsbbq %8,%q1\n\tsbbq %6,%q0"    \
       : "=r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "1"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "2"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "3"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "4"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_ddddddmmmmmmssssss(s5,s4,s3,s2,s1,s0, a5,a4,a3,a2,a1,a0, b5,b4,b3,b2,b1,b0)  \
  __asm__ ("subq %17,%q5\nsbbq %15,%q4\n\tsbbq %13,%q3\n\tsbbq %11,%q2\n\tsbbq %9,%q1\n\tsbbq %7,%q0"    \
       : "=r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "1"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "2"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "3"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "4"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "5"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

#define sub_dddddddmmmmmmmsssssss(s6,s5,s4,s3,s2,s1,s0, a6,a5,a4,a3,a2,a1,a0, b6,b5,b4,b3,b2,b1,b0)  \
  __asm__ ("subq %20,%q6\nsbbq %18,%q5\nsbbq %16,%q4\n\tsbbq %14,%q3\n\tsbbq %12,%q2\n\tsbbq %10,%q1\n\tsbbq %8,%q0"    \
       : "=r" (s6), "=&r" (s5), "=&r" (s4), "=&r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a6)), "rme" ((mp_limb_t)(b6)),                 \
         "1"  ((mp_limb_t)(a5)), "rme" ((mp_limb_t)(b5)),                 \
         "2"  ((mp_limb_t)(a4)), "rme" ((mp_limb_t)(b4)),                 \
         "3"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "4"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "5"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "6"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))

inline void _mul(ulong& hi, ulong& lo, ulong y, ulong x)
{
    __uint128_t p = ((__uint128_t) x) * ((__uint128_t) y);
    lo = (ulong) (p);
    hi = (ulong) (p >> 64);
}

inline void _madd(ulong& hi, ulong& lo, ulong y, ulong x)
{
    __uint128_t p = ((__uint128_t) lo) | (((__uint128_t) hi) << 64);
    p += ((__uint128_t) x) * ((__uint128_t) y);
    lo = (ulong) (p);
    hi = (ulong) (p >> 64);
}


#ifndef AVOID_MULX_ETC
inline ulong _mulx_ulong(ulong x, ulong y, ulong& p_hi)
{
    long long unsigned int _p_lo, _p_hi;
    _p_lo = _mulx_u64(static_cast<long long unsigned int>(x),
                      static_cast<long long unsigned int>(y),
                      &_p_hi);
    p_hi = static_cast<ulong>(_p_hi);
    return static_cast<ulong>(_p_lo);
}

inline ulong _addcarry_ulong(unsigned char cf, ulong x, ulong y, ulong& s)
{
    long long unsigned int _s;
    cf = _addcarry_u64(cf, static_cast<long long unsigned int>(x),
                           static_cast<long long unsigned int>(y),
                           &_s);
    s = static_cast<ulong>(_s);
    return cf;
}

inline ulong _addcarryx_ulong(unsigned char cf, ulong x, ulong y, ulong& s)
{
    long long unsigned int _s;
    cf = _addcarryx_u64(cf, static_cast<long long unsigned int>(x),
                            static_cast<long long unsigned int>(y),
                            &_s);
    s = static_cast<ulong>(_s);
    return cf;
}

inline ulong _subborrow_ulong(unsigned char cf, ulong x, ulong y, ulong& s)
{
    long long unsigned int _s;
    cf = _subborrow_u64(cf, static_cast<long long unsigned int>(x),
                           static_cast<long long unsigned int>(y),
                           &_s);
    s = static_cast<ulong>(_s);
    return cf;
}

#else

inline ulong _mulx_ulong(ulong x, ulong y, ulong& p_hi)
{
    long long unsigned int _p_lo, _p_hi;
    umul_ppmm(_p_hi, _p_lo, x, y);
    p_hi = static_cast<ulong>(_p_hi);
    return static_cast<ulong>(_p_lo);
}

inline ulong _addcarry_ulong(unsigned char cf, ulong x, ulong y, ulong& s)
{
    ulong shi, slo;
    add_ssaaaa(shi,slo, 0,x, 0,y);
    add_ssaaaa(shi,slo, shi,slo, 0,cf);
    s = static_cast<ulong>(slo);
    return shi;
}

inline ulong _addcarryx_ulong(unsigned char cf, ulong x, ulong y, ulong& s)
{
    ulong shi, slo;
    add_ssaaaa(shi,slo, 0,x, 0,y);
    add_ssaaaa(shi,slo, shi,slo, 0,cf);
    s = static_cast<ulong>(slo);
    return shi;
}

inline ulong _subborrow_ulong(unsigned char cf, ulong x, ulong y, ulong& s)
{
    ulong shi, slo;
    sub_ddmmss(shi,slo, 0,x, 0,y);
    sub_ddmmss(shi,slo, shi,slo, 0,cf);
    s = static_cast<ulong>(slo);
    return shi;
}

#endif

template<ulong n>
inline void _multi_add(ulong z[], const ulong a[])
{
    unsigned char cf = 0;
    ulong i = 0;
    do {
        cf = _addcarry_ulong(cf, a[i], z[i], z[i]);
    } while (++i < n);
}

template<ulong n>
inline void _multi_sub(ulong z[], const ulong a[])
{
    unsigned char cf = 0;
    ulong i = 0;
    do {
        cf = _subborrow_ulong(cf, z[i], a[i], z[i]);
    } while (++i < n);
}


// these wouldn't be necessary if gcc didn't generate suvh bad code

template<>
inline void _multi_add<2>(ulong z[], const ulong a[])
{
    add_ssaaaa(z[1],z[0],
               z[1],z[0],
               a[1],a[0]);
}

template<>
inline void _multi_add<3>(ulong z[], const ulong a[])
{
    add_sssaaaaaa(z[2],z[1],z[0],
                  z[2],z[1],z[0],
                  a[2],a[1],a[0]);
}

template<>
inline void _multi_add<4>(ulong z[], const ulong a[])
{
    add_ssssaaaaaaaa(z[3],z[2],z[1],z[0],
                     z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0]);
}

template<>
inline void _multi_add<5>(ulong z[], const ulong a[])
{
    add_sssssaaaaaaaaaa(z[4],z[3],z[2],z[1],z[0],
                        z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0]);
}

template<>
inline void _multi_add<6>(ulong z[], const ulong a[])
{
    add_ssssssaaaaaaaaaaaa(z[5],z[4],z[3],z[2],z[1],z[0],
                           z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0]);
}

template<>
inline void _multi_add<7>(ulong z[], const ulong a[])
{
    add_sssssssaaaaaaaaaaaaaa(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}


template<>
inline void _multi_sub<2>(ulong z[], const ulong a[])
{
    sub_ddmmss(z[1],z[0],
               z[1],z[0],
               a[1],a[0]);
}

template<>
inline void _multi_sub<3>(ulong z[], const ulong a[])
{
    sub_dddmmmsss(z[2],z[1],z[0],
                  z[2],z[1],z[0],
                  a[2],a[1],a[0]);
}

template<>
inline void _multi_sub<4>(ulong z[], const ulong a[])
{
    sub_ddddmmmmssss(z[3],z[2],z[1],z[0],
                     z[3],z[2],z[1],z[0],
                     a[3],a[2],a[1],a[0]);
}

template<>
inline void _multi_sub<5>(ulong z[], const ulong a[])
{
    sub_dddddmmmmmsssss(z[4],z[3],z[2],z[1],z[0],
                        z[4],z[3],z[2],z[1],z[0],
                        a[4],a[3],a[2],a[1],a[0]);
}

template<>
inline void _multi_sub<6>(ulong z[], const ulong a[])
{
    sub_ddddddmmmmmmssssss(z[5],z[4],z[3],z[2],z[1],z[0],
                           z[5],z[4],z[3],z[2],z[1],z[0],
                           a[5],a[4],a[3],a[2],a[1],a[0]);
}

template<>
inline void _multi_sub<7>(ulong z[], const ulong a[])
{
    sub_dddddddmmmmmmmsssssss(z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              z[6],z[5],z[4],z[3],z[2],z[1],z[0],
                              a[6],a[5],a[4],a[3],a[2],a[1],a[0]);
}

