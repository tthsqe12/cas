#pragma once

#ifndef _WIN32 // Linux - Unix
  #include <unistd.h>
  #include <time.h> 
  #include <sys/time.h>
  #include <sys/mman.h>
  #include <termios.h>
#else // Windows and MinGW
  #undef _WIN32_WINNT
  #define _WIN32_WINNT 0x0600
  #include <conio.h>
  #include <windows.h>
#endif

#include "packed.h"
#include "flint/mpn_extras.h"
#include "flint/nmod.h"

ulong GetMS();
void SleepMS(ulong ms);
ulong GetUS();
void SleepUS(ulong us);

// pass arguments by constant ref since by value is buggy on some platforms
template<typename T>
T eval_poly_mod(const T* a, ulong an, const T& b, const T& n, const T& ninv)
{
    T x = a[--an];
    while (an > 0)
        x = add(a[--an], mulmod2(x, b, n, ninv));
    return reduce_pm2n_to_pm1n(x, n);
}

std::complex<double> cispi_frac_ui(ulong a, ulong b);
std::complex<double> cispi_frac(slong a, ulong b);
std::complex<double> rev_cispi(ulong j, ulong h);

void out_mpn(const ulong* a, ulong an);
int mpn_cmp_ui_2exp(const ulong* a, ulong an, ulong b, ulong e);

struct format_hex {
    ulong data;
    int size;
    format_hex(uint64_t _data) : data(_data), size(16) {}
    format_hex(uint32_t _data) : data(_data), size(8) {}
};

std::ostream& operator<<(std::ostream& o, const format_hex& a);
std::string format_fixed(double x, ulong l);
std::string format_fixed(double x, ulong l, ulong r);

inline ulong saturate_bits(ulong a)
{
    a |= a >> 1;
    a |= a >> 2;
    a |= a >> 4;
    a |= a >> 8;
    a |= a >> 16;
    return a;
}

inline constexpr ulong cdiv(ulong a, ulong b)
{
    return (a + b - 1)/b;
}

inline constexpr ulong round_up(ulong a, ulong b)
{
    return cdiv(a, b)*b;
}

inline constexpr ulong pow2(int k)
{
    return ulong(1) << k;
}

inline constexpr ulong pow2(ulong k)
{
    return ulong(1) << k;
}

inline constexpr ulong leading_zeros(ulong x)
{
    return __builtin_clzll(x);
}

inline constexpr ulong clog2(ulong x)
{
    if (x <= 2)
        return x == 2;
    return 64 - __builtin_clzll(x - 1);
}

inline constexpr slong clog2(slong xx)
{
    if (xx <= 2)
        return xx == 2;
    ulong x = xx;
    return 64 - __builtin_clzll(x - 1);
}

inline ulong next_fft_number(ulong p)
{
    ulong bits = FLINT_BIT_COUNT(p);
    ulong l; count_trailing_zeros(l, p - 1);
    ulong q = p - (UWORD(2) << l);
    if (bits < 20)
        std::abort();
    if (FLINT_BIT_COUNT(q) == bits)
        return q;
    if (l < 5)
        return (UWORD(1) << (bits - 2)) + 1;
    return (UWORD(1) << (bits)) - (UWORD(1) << (l - 1)) + 1;
}
