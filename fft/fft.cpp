#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<string>
#include<memory>
#include<vector>
#include<stack>
#include<list>
#include<map>
#include<set>
#include<assert.h>
#include<string.h>
#include<cmath>
#include<new>

#include "flint/mpn_extras.h"
#include "flint/profiler.h"
#include "flint/fft.h"

#include "timing.h"
#include <immintrin.h>

ulong hits[9][9];

#include "packed.h"
#include "v1.cpp"
#include "v2.cpp"
#include "v2_test.cpp"

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



int main(int argc, char *argv[])
{
// graph
#if 1
    test_v2_trunc(10, 14);
    profile_v2_trunc(16, 27);
#endif

// candlestick
#if 0
    test_v2_trunc(10, 14);

    for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
        hits[i][j] = 0;

    profile_v2_trunc(16, 25);

    for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
        if (hits[i][j] != 0)
            std::cout << "[" << i << "][" << j << "] = " << hits[i][j] << std::endl;

    test_v2_mul(100000, 100);

    for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
        hits[i][j] = 0;

    profile_v2_mul(1000000, false);

    for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
        if (hits[i][j] != 0)
            std::cout << "[" << i << "][" << j << "] = " << hits[i][j] << std::endl;
#endif


#if 0
    test_v2_trunc(10, 14);
    profile_v2_trunc(17, 22);
    test_v2_mul(100000, 1000);
    profile_v2_mul(10000000, false);
#endif

// flint graph
#if 0
    profile_flint_trunc(12, 18);
#endif

    ulong maxlen = 5700000;

#if 0
std::cout << std::endl;
std::cout << "***********************************" << std::endl;
std::cout << "*********     old        **********" << std::endl;
std::cout << "***********************************" << std::endl;

    {
        mpn_fft_ctx ZZfft;
        ulong* a = new ulong[maxlen];
        ulong* b = new ulong[maxlen];
        ulong* c = new ulong[2*maxlen];
        ulong* d = new ulong[2*maxlen];

        for (int i = 0; i < maxlen; i++)
        {
            a[i] = -UWORD(1);
            b[i] = -UWORD(2);
        }

        my_PD_mpn_mul(ZZfft, c, a, maxlen, b, maxlen);

        for (slong an = maxlen/2; an <= maxlen; an += maxlen/2)
        for (slong bn = maxlen/2; bn <= an; bn += maxlen/2)
        {
            int reps = 1;
ulong t1 = GetMS();
            for (int i = 0; i < reps; i++) {
//                flint_mpn_mul_fft_main(d, a, an, b, bn);
                mpn_mul(d, a, an, b, bn);
            }

ulong t2 = GetMS();

            for (int i = 0; i < reps; i++)
                my_PD_mpn_mul(ZZfft, c, a, an, b, bn);

ulong t3 = GetMS();
std::cout << "** " << an << "*" << bn << ": ";
std::cout << " ratio: " << double(t2 - t1)/double(t3 - t2);
std::cout << "    old time: " << t2 - t1;
std::cout << " new time: " << t3 - t2;
std::cout << std::endl;

            bool matches = true;
            for (slong i = 0; i < an + bn; i++)
            {
                matches = matches && (c[i] == d[i]);
            }
            if (!matches)
                std::cout << "does not match!!!!" << std::endl;
        }

        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
    }
#endif


#if 0
    {
        ulong p = (UWORD(1)<<(50-1)) + 1;
        for (int j = 0; j < 20; j++)
        {
            std::cout << format_hex(p) << std::endl;
            p = next_fft_number(p);
        }
    }
#endif

#if 0
    for (ulong an = 100; an <= 2000; an += 100)
    {
std::cout << "an: " << an << std::endl;
        for (ulong bn = 100; bn <= 2000; bn += 100)
            test_mpn_mul_v2(an, bn);
    }

    {
        mpn_ctx_v2 ctx;
        for (ulong an = 1000000; an <= 2000000; an += 1000000)
        {
            for (ulong bn = 1000000; bn <= 2000000; bn += 1000000)
            {
                test_mpn_mul_v2(ctx, an, bn);
                test_mpn_mul_v2(ctx, an, bn);
            }
        }
    }
#endif

#if 0
    {
        mpn_fft_ctx ZZfft;
        slong maxlen = 20000;
        ulong* a = new ulong[maxlen];
        ulong* b = new ulong[maxlen];
        ulong* c = new ulong[2*maxlen];
        ulong* d = new ulong[2*maxlen];

        for (int i = 0; i < maxlen; i++)
        {
            a[i] = -UWORD(1);
            b[i] = -UWORD(2);
        }

        for (slong an = 3; an <= 50; an++)
        for (slong bn = 3; bn <= an; bn++)
        {
std::cout << "****************  " << an << "*" << bn << "  **************" << std::endl;
            mpn_mul(d, a, an, b, bn);
            my_PD_mpn_mul(ZZfft, c, a, an, b, bn);
            bool matches = true;
            for (slong i = 0; i < an + bn; i++)
            {
//                std::cout <<    "c[" << i << "]: 0x" << format_hex(c[i]);
//                std::cout <<  "  d[" << i << "]: 0x" << format_hex(d[i]);
//                std::cout << "    matches: " << (c[i] == d[i]) << std::endl;
                matches = matches && (c[i] == d[i]);
            }
if (!matches){
std::cout << "does not match!!!" << std::endl;
abort();
}
        }

        for (slong an = 4000; an <= maxlen; an += 2000)
        for (slong bn = 4000; bn <= an; bn += 2000)
        {
            int reps = 30;
ulong t1 = GetMS();
            for (int i = 0; i < reps; i++) {
//                flint_mpn_mul_fft_main(d, a, an, b, bn);
                mpn_mul(d, a, an, b, bn);
            }

ulong t2 = GetMS();

            for (int i = 0; i < reps; i++)
                my_PD_mpn_mul(ZZfft, c, a, an, b, bn);

ulong t3 = GetMS();
std::cout << "** " << an << "*" << bn << ": ";
std::cout << " ratio: " << double(t2 - t1)/double(t3 - t2);
std::cout << "    old time: " << t2 - t1;
std::cout << " new time: " << t3 - t2;
std::cout << std::endl;

            bool matches = true;
            for (slong i = 0; i < an + bn; i++)
            {
                matches = matches && (c[i] == d[i]);
            }
            if (!matches)
                std::cout << "does not match!!!!" << std::endl;
        }

/*
        int reps = 1000;
        slong smalln = 800;

ulong t1 = GetMS();
        for (int i = 0; i < reps; i++)
            flint_mpn_mul_fft_main(d, a, smalln, b, smalln);
ulong t2 = GetMS();
        for (int i = 0; i < reps; i++)
            my_PD_mpn_mul(ZZfft, c, a, smalln, b, smalln);
ulong t3 = GetMS();

std::cout << "** small: ";
std::cout << " ratio: " << double(t2 - t1)/double(t3 - t2);
std::cout << " gmp time: " << t2 - t1;
std::cout << "  my time: " << t3 - t2;
std::cout << std::endl;
*/

        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
    }

    if (0) {
        ulong n = (UWORD(1) << 63) - 1;
        nmod_poly_t a, b, c;
        nmod_poly_init(a,n);
        nmod_poly_init(b,n);
        nmod_poly_init(c,n);
        for (slong i = 1000000; i >= 0; i--)
        {
            nmod_poly_set_coeff_ui(a, i, 2*i*i-1);
            nmod_poly_set_coeff_ui(b, i, 3*i*i*i-2);
        }
std::cout << "deg a: " << nmod_poly_degree(a) << std::endl;
std::cout << "deg b: " << nmod_poly_degree(b) << std::endl;
ulong t1 = GetMS();
        nmod_poly_mul(c, a, b);
ulong t2 = GetMS();
std::cout << "poly mul: " << t2-t1 << std::endl;

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }
#endif

    flint_cleanup();

    std::cout << "PASS all good!" << std::endl;
    return 0;
}
