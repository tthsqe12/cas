#include <cassert>
#include "misc.h"


ulong GetMS() {
#ifndef _WIN32 // Linux - Unix
  struct timeval tm;
  gettimeofday( &tm, NULL );
  return tm.tv_sec*1000 + tm.tv_usec/1000;
#else // Windows and MinGW
    return GetTickCount64();
#endif
}

void SleepMS(ulong ms) {
#ifndef _WIN32 // Linux - Unix
    usleep(1000*ms);
#else // Windows and MinGW
    Sleep(ms);
#endif
}

ulong GetUS() {
#ifndef _WIN32 // Linux - Unix
  struct timeval tm;
  gettimeofday( &tm, NULL );
  return tm.tv_sec*1000000 + tm.tv_usec;
#else // Windows and MinGW
    return GetTickCount64()*1000;
#endif
}

void SleepUS(ulong us) {
#ifndef _WIN32 // Linux - Unix
    usleep(us);
#else // Windows and MinGW
    Sleep(us/1000);
#endif
}


slong fft_round_up(slong xxn, slong k)
{
    ulong n = UWORD(1) << k;
    ulong xn = xxn <= 0 ? 1 : xxn;
    xn = ((xn+((UWORD(1) << 4)-1)) >> 4) << 4; 

    if (k >= 10)
    {
        if (xn > n - (n >> 4))
            xn = n;
    }
    else
    {
        if (xn > n - (n >> 3))
            xn = n;
    }

   return slong(xn);
}


// exp(i*pi*a/b)
std::complex<double> cispi_frac_ui(ulong a, ulong b)
{
    arb_t x, y;
    fmpq_t f;
    arb_init(x);
    arb_init(y);
    fmpq_init(f);
    fmpq_set_ui(f, a, b);
    arb_sin_cos_pi_fmpq(y, x, f, 100);
    std::complex<double> z(arf_get_d(arb_midref(x), ARF_RND_NEAR),
                           arf_get_d(arb_midref(y), ARF_RND_NEAR));
    fmpq_clear(f);
    arb_clear(y);
    arb_clear(x);

    return z;
}

std::complex<double> cispi_frac(slong a, ulong b)
{
    if (a >= 0)
        return cispi_frac_ui(a, b);
    else
        return conj(cispi_frac_ui(-a, b));
}

std::complex<double> rev_cispi(ulong j, ulong h)
{
    if (j == 0)
        return 1;
    ulong n = FLINT_BIT_COUNT(j);
    return cispi_frac(h*n_revbin(j, n), pow2(n - 1));
}




void out_mpn(const ulong* a, ulong an)
{
    if (an == 0)
    {
        std::cout << "0";
        return;
    }

    for (ulong i = 0; i < an; i++)
    {
        if (i > 0)
            std::cout << " ";
        std::cout << format_hex(a[i]);
    }
}

// cmp(a, b*2^e)
int mpn_cmp_ui_2exp(const ulong* a, ulong an, ulong b, ulong e)
{
    ulong q = e/FLINT_BITS;
    ulong r = e%FLINT_BITS;
    ulong x, b0, b1;

    // allow one-off normalized input
    if (a[an-1] == 0)
        an--;

    assert(a[an-1] != 0);

    // b*2^e = (b*2^r       )*2^(64*q)
    //       = (b0 + b1*2^64)*2^(64*q)
    if (r == 0)
    {
        b0 = b;
        b1 = 0;
    }
    else
    {
        b0 = b << r;
        b1 = b >> (64 - r);
    }

    //      check words [q+2,infty)
    // then check words [q+1, 64*q+128) against b1
    // then check words [q, q+1) against b0
    // then check words [0, q)

    if (an > q + 2)
        return 1;

    x = (q+1 < an) ? a[q+1] : 0;
    if (x != b1)
        return x > b1 ? 1 : -1;

    x = (q < an) ? a[q] : 0;
    if (x != b0)
        return x > b0 ? 1 : -1;

    q = std::min(q, an);
    while (q > 0)
    {
        q--;
        if (a[q] != 0)
            return 1;
    }

    return 0;
}


std::ostream& operator<<(std::ostream& o, const format_hex& a)
{
    for (int i = a.size-1; i >= 0; i--)
    {
        ulong b = (a.data >> (4*i)) & 0x0f;
        o << char(b > 9 ? 'a' + b - 10 : '0' + b);
    }
    return o;
}

std::string format_fixed(double x, ulong l)
{
    std::string s(l, ' ');
    ulong y = std::abs(round(x));
    while (l > 0)
    {
        s[--l] = '0' + (y%10);
        y = y/10;
        if (y == 0)
            break;
    }

    return s;
}

std::string format_fixed(double x, ulong l, ulong r)
{
    std::string s(l+r+1, ' ');
    s[l] = '.';
    ulong y = std::abs(round(x*n_pow(10, r)));
    while (r > 0)
    {
        --r;
        s[l+1+r] = '0' + (y%10);
        y = y/10;
    }
    while (l > 0)
    {
        --l;
        s[l] = '0' + (y%10);
        y = y/10;
        if (y == 0)
            break;
    }

    return s;
}

