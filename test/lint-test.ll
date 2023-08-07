#define NATIVE_PACKED_BITS 256

#include <typeinfo>
#include "fmpz.h"
#include "test_helpers.h"
#include "timing.h"
#include "fp_poly.h"
#include "fp_mpoly.h"
#include "arb.h"
#include "mpn-impl.h"
#include "double_extras.h"
#include "asm_inlines.h"
#include <complex>

#include "generic/generic_spoly_gcd.h"
#include "generic/generic_poly_series.h"
#include "generic/generic_poly_divrem.h"

#include "fft_cpd.h"
#include "fft_cpd_basecases.h"

// number of terms k+1 should be a multiple of 8
// convert from z0 + z1/u + ... zk/u^k to l0. l1 l2 l3 l4 ... format

// first make sure each 0 <= z1 ... zk < u using rounding up
// then write z0 in twos complement

// reverse is similar: make sure |z1|, ..., |zk| <= u/2 and z0 gets the rest;


std::ostream& operator<<(std::ostream& o, const std::vector<std::complex<double>>& v)
{
   o << "{";
   for (auto x : v)
      o << x << ", ";
   o << "}";
   return o;
}

uword test_multiplier = 20;

fmpz_ring ZZ;

void dpoly_square_classical(
   std::vector<std::complex<double>>& x,
   const std::vector<std::complex<double>>& a)
{
   uword n = a.size();
   if (n < 1) {
      x.clear();
      return;
   }
   x.resize(2*n-1);
   for (uword i = 0; i < 2*n-1; i++) {
      x[i] = 0;
      for (uword j = 0; j < n; j++)
         if (i-j < n)
            x[i] += a[j]*a[i-j];
   }
}

double fetch_real(const std::vector<std::complex<double>>& a, uword i) {
   return (i < a.size()) ? a[i].real() : 0.0;
}

double fetch_imag(const std::vector<std::complex<double>>& a, uword i) {
   return (i < a.size()) ? a[i].imag() : 0.0;
}

/////// tables //////////

std::vector<std::complex<double>>cis2pirev_tab;   // e([0/1, 1/4, 1/8,  3/8, 1/16, 5/16, 3/16, 7/16,  ...])
std::vector<std::complex<double>>cis2pi1_tab[20]; //[k][i] = cis2pi_frac(i,1<<k)  k<=d, i<2^(k-1)
std::vector<std::complex<double>>cis2pi3_tab[20]; //[k][i] = cis2pi_frac(i,3<<k)  k<=d, i<1*2^k
std::vector<std::complex<double>>cis2pi5_tab[20]; //[k][i] = cis2pi_frac(i,5<<k)  k<=d, i<2^(k+1)


void fit_tab(uword depth)
{
   uword len;
   for (uword k = 0; k <= depth; k++)
   {
      len = (k<1) ? 0 : ui_mulpow2(1, k-1);
      cis2pi1_tab[k].resize(len);
      for (uword i = 0; i < len; i++)
         cis2pi1_tab[k][i] = my_cis2pi_frac(i,1<<k);

      len = ui_mulpow2(3, k);
      cis2pi3_tab[k].resize(len);
      for (uword i = 0; i < len; i++)
         cis2pi3_tab[k][i] = my_cis2pi_frac(i,3<<k);

      len = ui_mulpow2(5, k+1);
      cis2pi5_tab[k].resize(len);
      for (uword i = 0; i < len; i++)
         cis2pi5_tab[k][i] = my_cis2pi_frac(i,5<<k);
   }

   len = ui_pow2(depth-1);
   cis2pirev_tab.resize(len);
   for (uword i = 0; i < len; i++)
      cis2pirev_tab[i] = my_cis2pi_frac(n_revbin(i, depth-1), 2*len);
}

///////////////////// fft //////////////////

template<int V>
FORCE_NOINLINE void my_fft_one_block(
   complex_packed<double,V>* x,
   uword k, // V parallel transforms each of length 2^k
   uword j)
{
   const std::complex<double>* w2r = cis2pirev_tab.data();
   complex_packed<double,V> w, w2;
   if (k == 0)
   {
   }
   else if (k == 1)
   {
      radix_2_moth_block<V>(x[0], x[1], w2r[j]);
   }
   else if (k == 2)
   {
      radix_4_moth_block<V>(x[0], x[1], x[2], x[3], w2r[2*j], w2r[j]);
   }
   else if (k == 3)
   {
      w = w2r[j];
      complex_packed<double,V> x0 = x[0], x1 = x[1], x2 = x[2], x3 = x[3], x4 = x[4], x5 = x[5], x6 = x[6], x7 = x[7];
      radix_2_moth_block<V>(x0, x4, w);
      radix_2_moth_block<V>(x1, x5, w);
      radix_2_moth_block<V>(x2, x6, w);
      radix_2_moth_block<V>(x3, x7, w);
      radix_4_moth_block<V>(x0, x1, x2, x3, w2r[2*(2*j+0)], w2r[(2*j+0)]);
      x[0] = x0; x[1] = x1; x[2] = x2; x[3] = x3;
      radix_4_moth_block<V>(x4, x5, x6, x7, w2r[2*(2*j+1)], w2r[(2*j+1)]);
      x[4] = x4; x[5] = x5; x[6] = x6; x[7] = x7;
   }
   else
   {
      uword l2 = pow2(k - 2);
      w = w2r[2*j];
      w2 = w2r[j];
      for (uword i = 0; i < l2; i++)
         radix_4_moth_block<V>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
      my_fft_one_block<V>(x+0*l2, k-2, 4*j+0);
      my_fft_one_block<V>(x+1*l2, k-2, 4*j+1);
      my_fft_one_block<V>(x+2*l2, k-2, 4*j+2);
      my_fft_one_block<V>(x+3*l2, k-2, 4*j+3);
   }
}

template<uint32_t k, int V>
FORCE_NOINLINE void tmpl_fft_one_block(
   complex_packed<double,V>* x,
   uword j)
{
   const std::complex<double>* w2r = cis2pirev_tab.data();
   complex_packed<double,V> w, w2;
   if constexpr (k == 0)
   {
   }
   else if constexpr (k == 1)
   {
      radix_2_moth_block<V>(x[0], x[1], w2r[j]);
   }
   else if constexpr (k == 2)
   {
      radix_4_moth_block<V>(x[0], x[1], x[2], x[3], w2r[2*j], w2r[j]);
   }
   else if constexpr (k == 3)
   {
      w = w2r[j];
      complex_packed<double,V> x0 = x[0], x1 = x[1], x2 = x[2], x3 = x[3], x4 = x[4], x5 = x[5], x6 = x[6], x7 = x[7];
      radix_2_moth_block<V>(x0, x4, w);
      radix_2_moth_block<V>(x1, x5, w);
      radix_2_moth_block<V>(x2, x6, w);
      radix_2_moth_block<V>(x3, x7, w);
      radix_4_moth_block<V>(x0, x1, x2, x3, w2r[2*(2*j+0)], w2r[(2*j+0)]);
      x[0] = x0; x[1] = x1; x[2] = x2; x[3] = x3;
      radix_4_moth_block<V>(x4, x5, x6, x7, w2r[2*(2*j+1)], w2r[(2*j+1)]);
      x[4] = x4; x[5] = x5; x[6] = x6; x[7] = x7;
   }
   else
   {
      constexpr uword l2 = ui_pow2(k - 2);
      w = w2r[2*j];
      w2 = w2r[j];
      for (uword i = 0; i < l2; i++)
         radix_4_moth_block<V>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w2);
      tmpl_fft_one_block<k-2,V>(x+0*l2, 4*j+0);
      tmpl_fft_one_block<k-2,V>(x+1*l2, 4*j+1);
      tmpl_fft_one_block<k-2,V>(x+2*l2, 4*j+2);
      tmpl_fft_one_block<k-2,V>(x+3*l2, 4*j+3);
   }
}


template <int V>
FORCE_NOINLINE void my_ifft_one_block_sqr(
   complex_packed<double,V>* x,
   uword k) // 4 transformations of length 2^k
{
   if (k == 0)
   {
      x[0] = sqr(x[0]);
   }
   else if (k == 1)
   {
      complex_packed<double,V> x0 = sqr(x[0]), x1 = sqr(x[1]);
      radix_2_moth_1_block(x0, x1);
      x[0] = x0;
      x[1] = x1;
   }
   else if (k == 2)
   {
      complex_packed<double,V> x0 = sqr(x[0]), x1 = sqr(x[1]), x2 = sqr(x[2]), x3 = sqr(x[3]);
      radix_4_moth_inv_1_block(x0, x1, x2, x3);
      x[0] = x0;
      x[1] = x1;
      x[2] = x2;
      x[3] = x3;
   }
   else if (k == 3)
   {
      complex_packed<double,V> x0 = sqr(x[0]), x1 = sqr(x[1]), x2 = sqr(x[2]), x3 = sqr(x[3]), x4 = sqr(x[4]), x5 = sqr(x[5]), x6 = sqr(x[6]), x7 = sqr(x[7]);
      radix_2_moth_1_block<V>(x0, x1);
      radix_2_moth_1_block<V>(x2, x3);
      radix_2_moth_1_block<V>(x4, x5);
      radix_2_moth_1_block<V>(x6, x7);
      radix_4_moth_inv_1_block<V>(x0, x2, x4, x6);
      x[0] = x0; x[2] = x2; x[4] = x4; x[6] = x6;
      radix_4_moth_inv_block<V>(x1, x3, x5, x7, cis2pi1_tab[3][1], cis2pi1_tab[2][1]);
      x[1] = x1; x[3] = x3; x[5] = x5; x[7] = x7;
   }
   else
   {
      uword l2 = pow2(k-2);
      my_ifft_one_block_sqr<V>(x+0*l2, k-2);
      my_ifft_one_block_sqr<V>(x+1*l2, k-2);
      my_ifft_one_block_sqr<V>(x+2*l2, k-2);
      my_ifft_one_block_sqr<V>(x+3*l2, k-2);
      uword i = 0; do {
         auto w = cis2pi1_tab[k][i];
         auto w1 = cis2pi1_tab[k-1][i];
         radix_4_moth_inv_block<V>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w1);
      } while (i++, i < l2);
   }
}

template <uint32_t k, int V>
FORCE_NOINLINE void tmpl_ifft_one_block_sqr(
   complex_packed<double,V>* x)
{
   if constexpr (k == 0)
   {
      x[0] = sqr(x[0]);
   }
   else if constexpr (k == 1)
   {
      complex_packed<double,V> x0 = sqr(x[0]), x1 = sqr(x[1]);
      radix_2_moth_1_block(x0, x1);
      x[0] = x0;
      x[1] = x1;
   }
   else if constexpr (k == 2)
   {
      complex_packed<double,V> x0 = sqr(x[0]), x1 = sqr(x[1]), x2 = sqr(x[2]), x3 = sqr(x[3]);
      radix_4_moth_inv_1_block(x0, x1, x2, x3);
      x[0] = x0;
      x[1] = x1;
      x[2] = x2;
      x[3] = x3;
   }
   else if constexpr (k == 3)
   {
      complex_packed<double,V> x0 = sqr(x[0]), x1 = sqr(x[1]), x2 = sqr(x[2]), x3 = sqr(x[3]), x4 = sqr(x[4]), x5 = sqr(x[5]), x6 = sqr(x[6]), x7 = sqr(x[7]);
      radix_2_moth_1_block<V>(x0, x1);
      radix_2_moth_1_block<V>(x2, x3);
      radix_2_moth_1_block<V>(x4, x5);
      radix_2_moth_1_block<V>(x6, x7);
      radix_4_moth_inv_1_block<V>(x0, x2, x4, x6);
      x[0] = x0; x[2] = x2; x[4] = x4; x[6] = x6;
      radix_4_moth_inv_block<V>(x1, x3, x5, x7, cis2pi1_tab[3][1], cis2pi1_tab[2][1]);
      x[1] = x1; x[3] = x3; x[5] = x5; x[7] = x7;
   }
   else
   {
      constexpr uword l2 = ui_pow2(k-2);
      tmpl_ifft_one_block_sqr<k-2,V>(x+0*l2);
      tmpl_ifft_one_block_sqr<k-2,V>(x+1*l2);
      tmpl_ifft_one_block_sqr<k-2,V>(x+2*l2);
      tmpl_ifft_one_block_sqr<k-2,V>(x+3*l2);
      uword i = 0; do {
         auto w = cis2pi1_tab[k][i];
         auto w1 = cis2pi1_tab[k-1][i];
         radix_4_moth_inv_block<V>(x[i+0*l2], x[i+1*l2], x[i+2*l2], x[i+3*l2], w, w1);
      } while (i++, i < l2);
   }
}

template <int V>
FORCE_NOINLINE void my_fft_half_block(
   complex_packed<double,V>* x,
   uword m,
   uword k) // four parallel transformations of length m*2^k
{
   ASSERT(k >= 1);
   uword l = pow2(k);
   std::complex<double> one(1,0);

   if (m == 3) {
      for (uword i = 0; i < l/2; i++) {
         uword j = i+l/2;
         auto w = cis2pi3_tab[k][i];
         radix_3_moth_trunc_block<2,V>(x[i+0*l], x[i+1*l], x[i+2*l], w);
         w = cis2pi3_tab[k][j];
         radix_3_moth_trunc_block<1,V>(x[j+0*l], x[j+1*l], x[j+2*l], w);
      }
      for (uword mm = 0; mm < m; mm++)
         my_fft_one_block(x + mm*l, k, 0);
   } else if (m == 5) {
      uword i = 0;
      for (uword i = 0; i < l/2; i++) {
         auto w  = cis2pi5_tab[k][i];
         auto w2 = cis2pi5_tab[k-1][i];
         radix_5_moth_trunc_block<3,V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], w, w2);
      }
      for (uword i = 0; i < l/2; i++) {
         uword j = i+l/2;
         auto w  = cis2pi5_tab[k][j];
         auto w2 = cis2pi5_tab[k-1][j];
         radix_5_moth_trunc_block<2,V>(x[j+0*l],x[j+1*l],x[j+2*l],x[j+3*l],x[j+4*l], w, w2);
      }
      for (uword mm = 0; mm < m; mm++)
         my_fft_one_block(x + mm*l, k, 0);
   } else {
      uword l2 = pow2(k-1);
      for (uword i = 0; i < l2; i++)
         x[i+1*l2] = x[i+0*l2];
      my_fft_one_block<V>(x+0*l2, k-1, 0);
      my_fft_one_block<V>(x+1*l2, k-1, 1);
   }
}

template <uint32_t m, uint32_t k, int V>
FORCE_NOINLINE void tmpl_fft_half_block(
   complex_packed<double,V>* x)
{
   ASSERT(k >= 1);
   constexpr uword l = ui_pow2(k);
   std::complex<double> one(1,0);

   if constexpr (m == 3) {
      for (uword i = 0; i < l/2; i++) {
         uword j = i+l/2;
         auto w = cis2pi3_tab[k][i];
         radix_3_moth_trunc_block<2,V>(x[i+0*l], x[i+1*l], x[i+2*l], w);
         w = cis2pi3_tab[k][j];
         radix_3_moth_trunc_block<1,V>(x[j+0*l], x[j+1*l], x[j+2*l], w);
      }
      for (uword mm = 0; mm < m; mm++)
         tmpl_fft_one_block<k,V>(x + mm*l, 0);
   } else if constexpr (m == 5) {
      uword i = 0;
      for (uword i = 0; i < l/2; i++) {
         auto w  = cis2pi5_tab[k][i];
         auto w2 = cis2pi5_tab[k-1][i];
         radix_5_moth_trunc_block<3,V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], w, w2);
      }
      for (uword i = 0; i < l/2; i++) {
         uword j = i+l/2;
         auto w  = cis2pi5_tab[k][j];
         auto w2 = cis2pi5_tab[k-1][j];
         radix_5_moth_trunc_block<2,V>(x[j+0*l],x[j+1*l],x[j+2*l],x[j+3*l],x[j+4*l], w, w2);
      }
      for (uword mm = 0; mm < m; mm++)
         tmpl_fft_one_block<k,V>(x + mm*l, 0);
   } else {
      constexpr uword l2 = ui_pow2(k-1);
      for (uword i = 0; i < l2; i++)
         x[i+1*l2] = x[i+0*l2];
      tmpl_fft_one_block<k-1,V>(x+0*l2, 0);
      tmpl_fft_one_block<k-1,V>(x+1*l2, 1);
   }
}

template <int V>
FORCE_NOINLINE void my_fft_block(
   complex_packed<double,V>* x,
   uword m,
   uword k) // V parallel transformations of length m*2^k
{
   uword l = pow2(k);
   std::complex<double> one(1,0);

   if (m == 3) {
      uword i = 0;
      radix_3_moth_trunc_block<3,V>(x[i+0*l], x[i+1*l], x[i+2*l], one);
      for (i++; i < l; i++) {
         auto w = cis2pi3_tab[k][i];
         radix_3_moth_trunc_block<3,V>(x[i+0*l], x[i+1*l], x[i+2*l], w);
      }
   } else if (m == 5) {
      uword i = 0;
      radix_5_moth_trunc_block<5,V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], one, one);
      for (i++; i < l; i++) {
         auto w  = cis2pi5_tab[k][i];
         auto w2 = cis2pi5_tab[k-1][i];
         radix_5_moth_trunc_block<5,V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], w, w2);
      }
   }

   for (uword mm = 0; mm < m; mm++)
      my_fft_one_block<V>(x + mm*l, k, 0);
}

template <uint32_t m, uint32_t k, int V>
FORCE_NOINLINE void tmpl_fft_block(
   complex_packed<double,V>* x)
{
   constexpr uword l = ui_pow2(k);
   std::complex<double> one(1,0);

   if constexpr (m == 3) {
      uword i = 0;
      radix_3_moth_trunc_block<3,V>(x[i+0*l], x[i+1*l], x[i+2*l], one);
      for (i++; i < l; i++) {
         auto w = cis2pi3_tab[k][i];
         radix_3_moth_trunc_block<3,V>(x[i+0*l], x[i+1*l], x[i+2*l], w);
      }
   } else if constexpr (m == 5) {
      uword i = 0;
      radix_5_moth_trunc_block<5,V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], one, one);
      for (i++; i < l; i++) {
         auto w  = cis2pi5_tab[k][i];
         auto w2 = cis2pi5_tab[k-1][i];
         radix_5_moth_trunc_block<5,V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], w, w2);
      }
   }

   for (uword mm = 0; mm < m; mm++)
      tmpl_fft_one_block<k,V>(x + mm*l, k, 0);
}

// V parallel transformations of length m*2^k
template <int V>
FORCE_NOINLINE void my_ifft_block_sqr(
   complex_packed<double,V>* x,
   uword m,
   uword k)
{
   uword l = pow2(k);
   std::complex<double> one(1,0);

   for (uword mm = 0; mm < m; mm++)
      my_ifft_one_block_sqr<V>(x + mm*l, k);

   if (m == 3) {
      uword i = 0;
      radix_3_moth_inv_block<V>(x[i+0*l], x[i+1*l], x[i+2*l], one);
      for (i++; i < l; i++) {
         auto w = cis2pi3_tab[k][i];
         radix_3_moth_inv_block<V>(x[i+0*l], x[i+1*l], x[i+2*l], w);
      }
   } else if (m == 5) {
      uword i = 0;
      radix_5_moth_inv_block<V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], one, one);
      for (i++; i < l; i++) {
         auto w = cis2pi5_tab[k][i];
         auto w2 = cis2pi5_tab[k-1][i];
         radix_5_moth_inv_block<V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], w, w2);
      }
   }
}

template <uint32_t m, uint32_t k, int V>
FORCE_NOINLINE void tmpl_ifft_block_sqr(
   complex_packed<double,V>* x)
{
   constexpr uword l = ui_pow2(k);
   std::complex<double> one(1,0);

   for (uword mm = 0; mm < m; mm++)
      tmpl_ifft_one_block_sqr<k,V>(x + mm*l);

   if constexpr (m == 3) {
      uword i = 0;
      radix_3_moth_inv_block<V>(x[i+0*l], x[i+1*l], x[i+2*l], one);
      for (i++; i < l; i++) {
         auto w = cis2pi3_tab[k][i];
         radix_3_moth_inv_block<V>(x[i+0*l], x[i+1*l], x[i+2*l], w);
      }
   } else if constexpr (m == 5) {
      uword i = 0;
      radix_5_moth_inv_block<V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], one, one);
      for (i++; i < l; i++) {
         auto w = cis2pi5_tab[k][i];
         auto w2 = cis2pi5_tab[k-1][i];
         radix_5_moth_inv_block<V>(x[i+0*l],x[i+1*l],x[i+2*l],x[i+3*l],x[i+4*l], w, w2);
      }
   }
}


template <uint32_t mult, uint32_t depth, int V>
FORCE_NOINLINE void tmpl_mandel_iter(
   complex_packed<double,V>* z,
   const complex_packed<double,V>* c)
{
// (z7 + z6/u + ... + z0/u^7)^2 = s14 + s13/u + ... + s1/u^13 + s0/u^14
// (z7*u^7 + z6*u^6 + ... + z1*u + z0)^2 = s14*u^14 + s13*u^13 + ... + s1*u + s0
// write 8 units s14,s13,s12,s11,s10,s9,s8,s7
   constexpr uword len = ui_mulpow2(mult, depth);
   tmpl_fft_half_block<mult,depth,V>(z);
   tmpl_ifft_block_sqr<mult,depth,V>(z);
   packed<double,V> x, y, t; x = zero<double,V>(); y = zero<double,V>();
   packed<double,V> f = broadcast<double,V>(1.0/len);
   packed<double,V> u = broadcast<double,V>(65536.0);
   packed<double,V> uinv = broadcast<double,V>(1.0/65536.0);
   // the first len/2-1 coeffs just propogate rounding
   for (uword i = 0; i < len/2-1; i++) {
      x = add(x, round(mul(f, z[i].real())));
      y = add(y, round(mul(f, z[i].imag())));
      t = round(mul(x, uinv)); x = t;
      t = round(mul(y, uinv)); y = t;
   }
   // the next len/2 coeffs contribute to the answer
   for (uword i = 0; i < len/2; i++) {
      x = add(add(x,c[i].real()), round(mul(f, z[len/2-1+i].real())));
      y = add(add(y,c[i].imag()), round(mul(f, z[len/2-1+i].imag())));
      // x = t*u + r
      t = round(mul(x, uinv)); z[i].real(fnmadd(t, u, x)); x = t;
      t = round(mul(y, uinv)); z[i].imag(fnmadd(t, u, y)); y = t;
   }
   // the last coeff was zero, and we ignore the carry out   
}

template <int V>
FORCE_NOINLINE void mandel_iter(
   complex_packed<double,V>* z,
   const complex_packed<double,V>* c,
   uint32_t mult, uint32_t depth)
{
// (z7 + z6/u + ... + z0/u^7)^2 = s14 + s13/u + ... + s1/u^13 + s0/u^14
// (z7*u^7 + z6*u^6 + ... + z1*u + z0)^2 = s14*u^14 + s13*u^13 + ... + s1*u + s0
// write 8 units s14,s13,s12,s11,s10,s9,s8,s7

   if (depth <= 4) {
      static void (*tab[3*4])(complex_packed<double,V>* z, const complex_packed<double,V>*) =
         {  tmpl_mandel_iter<1,1,V>, tmpl_mandel_iter<3,1,V>, tmpl_mandel_iter<5,1,V>,
            tmpl_mandel_iter<1,2,V>, tmpl_mandel_iter<3,2,V>, tmpl_mandel_iter<5,2,V>,
            tmpl_mandel_iter<1,3,V>, tmpl_mandel_iter<3,3,V>, tmpl_mandel_iter<5,3,V>,
            tmpl_mandel_iter<1,4,V>, tmpl_mandel_iter<3,4,V>, tmpl_mandel_iter<5,4,V>
         };
      return tab[3*(depth-1)+(mult-1)/2](z, c);
   }


   sword len = ui_mulpow2(mult, depth);
   my_fft_half_block(z, mult, depth);
   my_ifft_block_sqr(z, mult, depth);
   packed<double,V> x, y, t; x = zero<double,V>(); y = zero<double,V>();
   packed<double,V> f = broadcast<double,V>(1.0/len);
   packed<double,V> u = broadcast<double,V>(65536.0);
   packed<double,V> uinv = broadcast<double,V>(1.0/65536.0);
   // the first len/2-1 coeffs just propogate rounding
   for (uword i = 0; i < len/2-1; i++) {
      x = add(x, round(mul(f, z[i].real())));
      y = add(y, round(mul(f, z[i].imag())));
      t = round(mul(x, uinv)); x = t;
      t = round(mul(y, uinv)); y = t;
   }
   // the next len/2 coeffs contribute to the answer
   for (uword i = 0; i < len/2; i++) {
      x = add(add(x,c[i].real()), round(mul(f, z[len/2-1+i].real())));
      y = add(add(y,c[i].imag()), round(mul(f, z[len/2-1+i].imag())));
      // x = t*u + r
      t = round(mul(x, uinv)); z[i].real(fnmadd(t, u, x)); x = t;
      t = round(mul(y, uinv)); z[i].imag(fnmadd(t, u, y)); y = t;
   }
   // the last coeff was zero, and we ignore the carry out   
}



void dpoly_square_fft(
   std::vector<std::complex<double>>& x,
   const std::vector<std::complex<double>>& a)
{
   uword n = a.size();

// (z0 + z1/u + z2/u^2 + ... + z7/u^7)^2 = w0 + w1/u + ... + w14/u^14
// (z0*u^7 + z1*u^6 + z2*u^5 + ... + z7*u^0)^2 = w0*u^14 + ... + w14*u^0  mod u^16-1

   for (uword mult = 1; mult <= 5; mult += 2) {
      uword depth = 3;
      uword convlen = ui_mulpow2(mult, depth);
   
      global_fft_cpd_ctx.fit_wtab(10);
      x.resize(ui_mulpow2(mult,depth));
   
      tmp_allocator push;
      complex_packed<double,4>* z = push.recursive_alloc<complex_packed<double,4>>(convlen);
   
      for (uword i = 0; i < convlen; i++) {
         z[i].real(packed<double,4>(fetch_real(a,i),0.0,0.0,0.0));
         z[i].imag(packed<double,4>(fetch_imag(a,i),0.0,0.0,0.0));
      }
   
      my_fft_block(z, mult, depth); // four parallel transformation length mult*2^depth
      my_ifft_block_sqr(z, mult, depth);
   
      for (uword i = 0; i < convlen; i++) {
         z[i].real( mul(1.0/convlen, z[i].real()));
         z[i].imag( mul(1.0/convlen, z[i].imag()));
         x[i] = std::complex<double>(std::round(z[i].real()[0]), std::round(z[i].imag()[0]));
      }
      std::cout << "x: " << x << "\n";
   }
}

void profile_fft(uword mult, uword depth)
{
   uword len = ui_mulpow2(mult, depth);
   tmp_allocator push;
   constexpr int V = 8;
   double* z = push.recursive_aligned_alloc<double,64>(len*2*V);
   double* c = push.recursive_aligned_alloc<double,64>(len*2*V);

   for (int j = 0; j < V; j++) {
      for (sword i = 0; i < len/2; i++) {
         c[2*V*i + 0*V + j] = (i+2)*(+j+3);
         c[2*V*i + 1*V + j] = (i+1)*(-j-2);
      }
   }

   complex_packed<double,V>* Z = reinterpret_cast<complex_packed<double,V>*>(z);
   complex_packed<double,V>* C = reinterpret_cast<complex_packed<double,V>*>(c);
   for (sword i = 0; i < len/2; i++)
      Z[i] = C[i];

   uword t1 = get_ms();
   uword nreps = 100 + 100000000/len;
   for (uword reps = 0; reps < nreps; reps++)
      mandel_iter(Z, C, mult, depth);

/*
   4 cores / 8 threads  with schoolbook mpn arithmetic
   qwords| MIPS
   ------+------
       2 | 380
       3 | 280
       4 | 210
       6 | 120
       8 |  75

fft:
length  8 | mips 51 | qwords 2
length 10 | mips 34 | qwords 3
length 12 | mips 29 | qwords 3
length 16 | mips 20 | qwords 5
length 20 | mips 14 | qwords 6
length 24 | mips 12 | qwords 7
length 32 | mips 9 | qwords 10
length 40 | mips 6 | qwords 12
length 48 | mips 5 | qwords 15

*/

   uword t2 = get_ms();
   std::cout << "length " << len/2 << " | mips " << (nreps*V)/(1000*(t2-t1)) << " | qwords " << (20*(len/2)/64) << std::endl;
}

void test_fft()
{
   {
      std::vector<std::complex<double>> x, a;
      uword n = 2;
      a.resize(n);
      for (uword i = 0; i < n; i++)
         a[i] = std::complex<double>(2*i+1,3*i+3);
   
      dpoly_square_classical(x, a);
   std::cout << "a: " << a << "\n";
   std::cout << "x: " << x << "\n";
   
      dpoly_square_fft(x, a);
   }

#if 1
   uword depth_limit = 10;
   global_fft_cpd_ctx.fit_wtab(10);

   tmp_allocator push;
   complex_packed<double,4>* a = push.recursive_alloc<complex_packed<double,4>>(ui_mulpow2(5, depth_limit));
   complex_packed<double,4>* b = push.recursive_alloc<complex_packed<double,4>>(ui_mulpow2(5, depth_limit));
   complex_packed<double,4>* c = push.recursive_alloc<complex_packed<double,4>>(ui_mulpow2(5, depth_limit));

   for (uword depth = 1; depth <= depth_limit; depth++) { 
   for (uword mult = 1; mult <= 5; mult += 2) {
      // l = convolution length
      uword l = ui_mulpow2(mult, depth);
      // a = rand mod X^l-1
      for (sword i = 0; i < l; i++) {
         a[i].real(packed<double,4>(+3*i+0, +4*i+1, -7*i+2, +6*i-7));
         a[i].imag(packed<double,4>(-5*i+2, -2*i+1, +4*i+2, -1*i-9));
      }
      // b = a*a mod X^l-1
      for (uword i = 0; i < l; i++) {
         packed<double,4> x(0), y(0);
         for (uword j = 0; j < l; j++) {
            // += a[j]*a[i-j]
            packed<double,4> ax = a[j].real();
            packed<double,4> ay = a[j].imag();
            packed<double,4> bx = a[(i-j+l)%l].real();
            packed<double,4> by = a[(i-j+l)%l].imag();
            x = fmadd(ax, bx, x);
            x = fnmadd(ay, by, x);
            y = fmadd(ax, by, y);
            y = fmadd(ay, bx, y);
         }
         b[i].real(x); b[i].imag(y);
      }
      // c = a^2 mod X^l-1
      for (uword i = 0; i < l; i++) {
         c[i].real(mul(1.0, a[i].real()));
         c[i].imag(mul(1.0, a[i].imag()));
      }
      my_fft_block(c, mult, depth);
      my_ifft_block_sqr(c, mult, depth);
      for (uword i = 0; i < l; i++) {
         c[i].real(mul(1.0/l, c[i].real()));
         c[i].imag(mul(1.0/l, c[i].imag()));
         for (int j = 0; j < 4; j++) {
            TEST(abs(c[i].real()[j] - b[i].real()[j]) < 1e-3, i, j, c[i].real()[j], b[i].real()[j]);
            TEST(abs(c[i].imag()[j] - b[i].imag()[j]) < 1e-3, i, j, c[i].imag()[j], b[i].imag()[j]);
         }
      }

      if (depth < 1)
         continue;

      // a = rand top half zero
      for (sword i = 0; i < l/2; i++) {
         a[i].real(packed<double,4>(+3*i+0, +4*i+1, -7*i+2, +6*i-7));
         a[i].imag(packed<double,4>(-5*i+2, -2*i+1, +4*i+2, -1*i-9));
      }
      for (sword i = l/2; i < l; i++) {
         a[i].real(packed<double,4>(0,0,0,0));
         a[i].imag(packed<double,4>(0,0,0,0));
      }
      // b = a*a mod
      for (uword i = 0; i < l; i++) {
         packed<double,4> x(0), y(0);
         for (uword j = 0; j <= i; j++) {
            if (j >= l/2 || i-j >= l/2)
               continue;
            // += a[j]*a[i-j]
            packed<double,4> ax = a[j].real();
            packed<double,4> ay = a[j].imag();
            packed<double,4> bx = a[i-j].real();
            packed<double,4> by = a[i-j].imag();
            x = fmadd(ax, bx, x);
            x = fnmadd(ay, by, x);
            y = fmadd(ax, by, y);
            y = fmadd(ay, bx, y);
         }
         b[i].real(x); b[i].imag(y);
      }
      // c = a^2 mod X^l-1
      for (uword i = 0; i < l; i++) {
         c[i].real(a[i].real());
         c[i].imag(a[i].imag());
      }
      my_fft_half_block(c, mult, depth);
      my_ifft_block_sqr(c, mult, depth);
      for (uword i = 0; i < l; i++) {
         c[i].real(mul(1.0/l, c[i].real()));
         c[i].imag(mul(1.0/l, c[i].imag()));
         for (int j = 0; j < 4; j++) {
            TEST(abs(c[i].real()[j] - b[i].real()[j]) < 1e-3, i, j, c[i].real()[j], b[i].real()[j]);
            TEST(abs(c[i].imag()[j] - b[i].imag()[j]) < 1e-3, i, j, c[i].imag()[j], b[i].imag()[j]);
         }
      }

   }
std::cout << "depth: " << depth << std::endl;

   }
#endif
}


void test_fmpq(random_state& state)
{
   fmpq_t a;

   fmpz_set_ui(a->num, 1);
   fmpz_set_ui(a->den, 30);

   double d = fmpq_round_d(a, ARF_RND_NEAR);

std::cout << "d: " << format_base_form(d, 3) << std::endl;

}

void test_arb()
{
   arb x, y;
   fmpz u;

   for (sword p = 192; p <= 192; p += 64)
{
/*
   for (uword k = 1; k <= 16; k++)
   {
      arb_set_ui(&y, 10);
      arb_log(&y, &y, p);
      arb_set_ui(&x, 10*k);
      arb_mul_2exp_si(&x, &x, -4);
      arb_log(&x, &x, p);
      arb_div(&x, &x, &y, p);
      arb_mul_2exp_si(&x, &x, 64+1);
      arf_get_fmpz(u, &x.mid, ARF_RND_CEIL);
      std::cout << k << ": " << format_hex(u) << std::endl;
   }
*/
   arb_set_si(&x, 2);
   arb_set_si(&y, 10);
   arb_log(&x, &x, p);
   arb_log(&y, &y, p);
   arb_div(&x, &x, &y, p);
//    write(std::cout, &x) << "  or  "; write_base_form(std::cout, &x, 10); std::cout << std::endl;

   arb_set_si(&x, 5);
//    write(std::cout, &x) << "  or  "; write_base_form(std::cout, &x, 10); std::cout << std::endl;

   arb_set_si(&x, 5);
//    write(std::cout, &x) << "  or  "; write_base_form(std::cout, &x, 10); std::cout << std::endl;

   arb_log(&x, &x, p);
//    write(std::cout, &x) << "  or  "; write_base_form(std::cout, &x, 10); std::cout << std::endl;

   arb_exp(&x, &x, p);
//    write(std::cout, &x) << "  or  "; write_base_form(std::cout, &x, 10); std::cout << std::endl;

   arb_rsqrt(&x, &x, p);
//    write_base_form(std::cout, &x, 10); std::cout << std::endl;

   uword f = 1;
   for (int j = 10; j > 0; j--)
   {
      f *= 10;
      arb_mul_ui(&y, &x, f, p);
//        write_base_form(std::cout, &y, 10); std::cout << std::endl;
   }

   f = 1;
   for (int j = 10; j > 0; j--)
   {
      f *= 10;
      arb_div_ui(&y, &x, f, p);
//        write_base_form(std::cout, &y, 10); std::cout << std::endl;
   }

}
}

void test_mpoly(random_state& state)
{
   fmpz p(2);
   fp_mpoly_ring R(p, 8);
   fp_mpoly f, g, h;

TESTSET_BEGIN("mpoly", ": proved irreducible ....")

//    set(R, f, "x0^0+x1^1+x2^2+x3^3+x4^4");
//    TEST(mpoly_proved_irreducible(R.m, f.m, state));

   set(R, f, "x0^809*x1^75*x2^384*x3^324*x4^74*x5^788*x6^83*x7^414+"
            "x0^805*x1^343*x2^595*x3^246*x4^32*x5^90*x6^473*x7^591+"
            "x0^718*x1^108*x2^680*x3^368*x4^358*x7^276+"
            "x0^683*x1^533*x3^649*x4^619*x5^136*x6^223*x7^610+"
            "x1^617*x2^777*x3^799*x4^443*x5^545*x6^166*x7^216+"
            "x0^485*x1^646*x2^424*x3^265*x4^416*x5^400*x6^278+"
            "x0^336*x1^149*x2^361*x3^691*x4^629*x5^282*x6^530*x7^259+"
            "x0^266*x2^258*x4^422*x5^637*x6^244*x7^236+"
            "x0^74*x1^812*x2^162*x3^417*x4^71*x5^188*x6^258*x7^637+"
            "x0^37*x1^604*x2^94*x3^474*x5^853*x6^521*x7^250");
   TEST(mpoly_proved_irreducible(R.m, f.m, state));

   R.m.set_nvars(2);
   set(R, f, "x1^2639+x0^4432*x1^2436+x0^400*x1^1827+x0^1300");
   TEST(mpoly_proved_irreducible(R.m, f.m, state));

   set(R, f, "x1^5481+x0^2*x1^5477+x0^4167*x1^4366+x0^2700");
   TEST(mpoly_proved_irreducible(R.m, f.m, state));

   for (uword i = 100*test_multiplier; i > 0; i--)
   {
      fmpz_next_prime(p, p, false);
      R.c.set_modulus(p);
      R.m.set_nvars(1 + state.get_mod(5));
      random_upto_length_and_degree_bits(state, R, g, 2+state.get_mod(10), 5+state.get_mod(5));
      random_upto_length_and_degree_bits(state, R, h, 2+state.get_mod(10), 5+state.get_mod(5));
      if (g.length() < 2 || h.length() < 2)
         continue;
      mul(R, f, g, h);
      TEST(!mpoly_proved_irreducible(R.m, f.m, state), with(R, f), with(R, g), with(R, h));
   }

TESTSET_END()

}


void check_gcd(fp_mpoly_ring& R,
   fp_mpoly& g,
   fp_mpoly& abar,
   fp_mpoly& bbar,
   fp_mpoly& a,
   fp_mpoly& b,
   const fp_mpoly& gdiv)
{
   fp_mpoly ca, cb, cg, u, v, w;

   gcdc(R, g, abar, bbar, a, b);
   TEST(is_canonical(R, g));
   TEST(is_canonical(R, abar));
   TEST(is_canonical(R, bbar));
   TEST(divides(R, ca, g, gdiv));

   mul(R, ca, g, abar);
   mul(R, cb, g, bbar);
   TEST(equal(R, ca, a) && equal(R, cb, b), with(R,g), with(R,abar), with(R,bbar), with(R,a), with(R,b));

   if (is_zero(R, g))
   {
      TEST(is_zero(R, a) && is_zero(R, b));
      return;
   }

   TEST(is_one(R.c, g.coeffs()), with(R, g));

   set(R, u, b);
   gcdc(R, u, v, w, a, u);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, v, b);
   gcdc(R, u, v, w, a, v);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, w, b);
   gcdc(R, u, v, w, a, w);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, u, a);
   gcdc(R, u, v, w, u, b);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, v, a);
   gcdc(R, u, v, w, v, b);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, w, a);
   gcdc(R, u, v, w, w, b);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, u, a);
   set(R, v, b);
   gcdc(R, u, v, w, u, v);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, v, a);
   set(R, u, b);
   gcdc(R, u, v, w, v, u);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, u, a);
   set(R, w, b);
   gcdc(R, u, v, w, u, w);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, w, a);
   set(R, u, b);
   gcdc(R, u, v, w, w, u);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, v, a);
   set(R, w, b);
   gcdc(R, u, v, w, v, w);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   set(R, w, a);
   set(R, v, b);
   gcdc(R, u, v, w, w, v);
   TEST(equal(R, g, u) && equal(R, abar, v) && equal(R, bbar, w));

   gcdc(R, cg, ca, cb, abar, bbar);
   TEST(is_canonical(R, cg));
   TEST(is_one(R, cg));
   TEST(equal(R, ca, abar) && equal(R, cb, bbar));

   gcdc(R, cg, abar, bbar, abar, bbar);
   TEST(is_canonical(R, cg));
   TEST(equal(R, ca, abar) && equal(R, cb, bbar));
}

void check_factor(fp_mpoly_ring& R,
   fp_mpoly& a,
   uword lower,
   uword upper)
{
   fp_mpoly q, t;
   fp_mpoly_product g, h;
   fmpz omega;

   factor(R, g, a);

   for (uword i = 0; i < g.length(); i++)
   {
      TEST(g.base(i).length() > 0);
      TEST(is_one(R.c, g.base(i).coeffs()));
   }

   fmpz_zero(omega);
   for (uword i = 0; i < g.length(); i++)
      fmpz_add(omega, omega, g.exp(i));
   TEST(fmpz_cmp_ui(omega, lower) >= 0);
   TEST(fmpz_cmp_ui(omega, upper) <= 0);

   expand(R, t, g);
   TEST(divides(R, q, a, t));
   TEST(is_constant(R, q));

   for (uword i = 0; i < g.length(); i++)
   {
      factor(R, h, g.base(i));
      TEST(h.length() == 1 && fmpz_is_one(h.exp(0)));
   }
}



template <typename RingT>
void test_prime_field_mpoly(random_state& state)
{
   fmpz m(2);
   typename RingT::mpoly_ring_type R(m, 1);
   typename RingT::mpoly_type a, b, c, d, e, f, g, q, r, s, t;

TESTSET_BEGIN(typeid(typename RingT::mpoly_type).name(), ": add sub mul divides ....")

   for (uword i = 10*test_multiplier; i > 0; i--)
   {
      random_of_bits(state, ZZ, m, 2 + state.get_mod(200));
      R.c.set_modulus(m);
      R.m.set_nvars(1 + state.get_mod(20));
      for (uword j = 10; j > 0; j--)
      {
         random_upto_length_and_degree_bits(state, R, a, 20, 1+state.get_mod(150));
         random_upto_length_and_degree_bits(state, R, b, 20, 1+state.get_mod(150));
         random_upto_length_and_degree_bits(state, R, c, 20, 1+state.get_mod(150));

         add(R, d, b, c);
         mul(R, e, a, d);
         mul(R, f, a, b);
         mul(R, g, a, c);
         add(R, d, f, g);
         TEST(equal(R, d, e), with(R, a), with(R, b), with(R, c));

         sub(R, d, b, c);
         mul(R, e, a, d);
         mul(R, f, a, b);
         mul(R, g, a, c);
         sub(R, d, f, g);
         TEST(equal(R, d, e), with(R, a), with(R, b), with(R, c));

         add(R, a, b, c);
         sub(R, d, a, c);
         TEST(equal(R, d, b), with(R, a), with(R, b), with(R, c), with(R, d));

         sub(R, e, d, b);
         TEST(is_zero(R, e), with(R, e));

         try {
            random_upto_length_and_degree_bound(state, R, a, 20, state.get_mod(10));
            random_upto_length_and_degree_bound(state, R, b, 20, state.get_mod(10));

            mul(R, c, a, b);
            TEST(divides(R, d, c, b), R.c, with(R, a), with(R, b));
         }
         catch (const parent_and_elem<fp_ring, fmpz>& ex)
         {
            TEST(&ex.parent() == &R.c, &ex.parent(), &R.c);
            TEST(fmpz_cmp_si(ex.elem(), 1) > 0, ex.elem());
            TEST(fmpz_cmp(R.c.characteristic(), ex.elem()) > 0, R.c.characteristic(), ex.elem());
            TEST(fmpz_divisible(R.c.characteristic(), ex.elem()), R.c.characteristic(), ex.elem());
         }
      }
   }

TESTSET_END()

TESTSET_BEGIN(typeid(typename RingT::mpoly_type).name(), ": gcd monomial ....")

   for (uword i = 10*test_multiplier; i > 0; i--)
   {
      random_of_bits(state, ZZ, m, 2 + state.get_mod(200));
      fmpz_next_prime(m, m, false);
      R.c.set_modulus(m);
      R.m.set_nvars(1 + state.get_mod(20));
      for (uword j = 5; j > 0; j--)
      {
         random_upto_length_and_degree_bits(state, R, t, 1, 2 + state.get_mod(10));
         random_upto_length_and_degree_bits(state, R, a, 1, 2 + state.get_mod(10));
         random_upto_length_and_degree_bits(state, R, b, state.get_mod(15), 2 + state.get_mod(10));
         mul(R, a, t);
         mul(R, b, t);
         if (state.get_bit())
            swap(R, a, b);
         check_gcd(R, g, q, r, a, b, t);

         random_upto_length_and_degree_bits(state, R, t, state.get_mod(15), 2 + state.get_mod(10));
         random_upto_length_and_degree_bits(state, R, a, 1, 2 + state.get_mod(10));
         random_upto_length_and_degree_bits(state, R, b, 1, 2 + state.get_mod(10));
         mul(R, a, t);
         mul(R, b, t);
         check_gcd(R, g, q, r, a, b, t);
      }
   }

TESTSET_END()

TESTSET_BEGIN(typeid(typename RingT::mpoly_type).name(), ": gcd .............")

   R.c.set_modulus(fmpz(101));
   R.m.set_nvars(2);
   set(R, a, "x0^3-x1^3");
   set(R, b, "x0^2-x1^2");
   set(R, t, "x0-x1");
   check_gcd(R, g, q, r, a, b, t);

   for (uword i = 1*test_multiplier; i > 0; i--)
   {
      random_of_bits(state, ZZ, m, 2 + state.get_mod(200));
      fmpz_next_prime(m, m, false);
      R.c.set_modulus(m);
      R.m.set_nvars(1 + state.get_mod(3));

      uword degbound = 1 + 8/R.m.nvars();
      for (uword j = 4; j > 0; j--)
      {
         random_upto_length_and_degree_bound(state, R, t, 1 + state.get_mod(12), degbound);
         random_upto_length_and_degree_bound(state, R, a, state.get_mod(12), degbound);
         random_upto_length_and_degree_bound(state, R, b, state.get_mod(12), degbound);
/*
std::cout << "-------------------------------------" << std::endl;
std::cout << "nvars: " << R.m.nvars() << std::endl;
std::cout << "coefficient ring: " << R.c << std::endl;
std::cout << "a: " << with(R, a) << std::endl;
std::cout << "b: " << with(R, b) << std::endl;
std::cout << "t: " << with(R, t) << std::endl;
*/
         mul(R, a, t);
         mul(R, b, t);
         check_gcd(R, g, q, r, a, b, t);
      }
   }

TESTSET_END()

TESTSET_BEGIN(typeid(typename RingT::mpoly_type).name(), ": factor ..........")

   R.c.set_modulus(fmpz(101));
   R.m.set_nvars(2);
   set(R, a, "x0^2-x1^2");
   check_factor(R, a, 2, 2);
   set(R, a, "x0^3-1");
   check_factor(R, a, 2, 2);

TESTSET_END()

}

template <typename RingT>
void test_squarefree_factorization(
   RingT& R,
   typename RingT::poly_type::source_type a,
   uword omega_min)
{
   typename RingT::poly_type::product_type f;
   typename RingT::poly_type p, q;

   factor_squarefree(R, f, a);
   shift_right(R, p, a, a.degree());
   uword omega = 0;
   for (uword i = 0; i < f.length(); i++)
   {
      omega += f.exp(i);
      TEST(f.exp(i) > 0, f.exp(i));
      pow(R, q, f.base(i), f.exp(i));
      mul(R, p, q);
      TEST(is_monic(R, f.base(i)), R, with(R, f.base(i)));
      TEST(is_squarefree(R, f.base(i)), R, with(R, f.base(i)), with(R, a), with(R, f));
   }
   TEST(equal(R, p, a), R, with(R, p), with(R, a));
   TEST(omega >= omega_min, R, with(R, a), with(R, f));
}

template <typename RingT>
void test_factorization(
   RingT& R,
   typename RingT::poly_type::source_type a,
   uword omega_min)
{
   typename RingT::poly_type::product_type f;
   typename RingT::poly_type p, q;

   factor(R, f, a);
   shift_right(R, p, a, a.degree());
   uword omega = 0;
   for (uword i = 0; i < f.length(); i++)
   {
      omega += f.exp(i);
      TEST(f.exp(i) > 0, f.exp(i));
      pow(R, q, f.base(i), f.exp(i));
      mul(R, p, q);
      TEST(is_monic(R, f.base(i)), R, with(R, f.base(i)));
      TEST(is_irreducible(R, f.base(i)), R, with(R, f.base(i)), with(R, a), with(R, f));
   }
   TEST(equal(R, p, a), R, with(R, a), with(R, f));
   TEST(omega >= omega_min, R, with(R, a), with(R, f));
}


template <typename RingT>
void test_prime_field_poly(random_state& state)
{
   fmpz m(2);
   RingT R(m);
   typename RingT::poly_type a, b, c, d, e, f, g, q, r, s, t;

   {
      R.set_modulus(101);
      set(R, a, "1-x");
      generic_poly_inv_series<RingT>(R, b, a, 10);
      generic_poly_inv_series<RingT>(R, c, b, 10);
      TEST(equal(R, a, c), with(R, a), with(R, c));

      set(R, a, "1-x+x^2 + x^15");
      generic_poly_inv_series<RingT>(R, b, a, 20);
      generic_poly_inv_series<RingT>(R, c, b, 20);
      TEST(equal(R, a, c), with(R, a), with(R, c));

      set(R, a, "1-x+x^2+3*x^3+4*x^4+5*x^5+6*x^6+x^7");
      generic_poly_invert<RingT>(R, b, a, 6);
      mul(R, c, a, b);
   }


TESTSET_BEGIN(typeid(typename RingT::poly_type).name(), ": add sub mul ........")

   for (uword i = 0; i < 500*test_multiplier; i++)
   {
      random_of_bits(state, ZZ, m, 2 + state.get_mod(250));
      R.set_modulus(m);

      random_of_length(state, R, a, state.get_mod(20));
      random_of_length(state, R, b, state.get_mod(20));
      sub(R, c, a, b);
      sub(R, d, b, a);
      add(R, c, c, d);
      TEST(is_zero(R, c), with(R, c));

      random_of_length(state, R, a, state.get_mod(20));
      random_of_length(state, R, b, state.get_mod(20));
      set(R, e, a);
      add(R, a, a, b);
      sub(R, a, a, b);
      TEST(equal(R, e, a));

      sub(R, a, c, b);
      sub(R, c, b, c);
      neg(R, c, c);
      TEST(equal(R, a, c));

      zero(R, a);
      random_of_length(state, R, b, state.get_mod(20));
      add(R, c, a, b);
      sub(R, d, c, b);
      TEST(is_zero(R, d));
      TEST(is_zero(R, a));

      set(R, e, a);
      add(R, a, a, b);
      sub(R, a, a, b);
      TEST(equal(R, e, a));

      sub(R, e, c, b);
      sub(R, t, b, c);
      neg(R, t);
      TEST(equal(R, t, e));

      random_of_length(state, R, a, state.get_mod(20));
      random_of_length(state, R, b, state.get_mod(20));
      random_of_length(state, R, c, state.get_mod(20));
      add(R, t, b, c);
      mul(R, d, a, t);
      mul(R, t, a, b);
      mul(R, s, a, c);

      add(R, e, t, s);
      TEST(equal(R, d, e), with(R, d), with(R, e));

      random_of_length(state, R, a, state.get_mod(20));
      random_of_length(state, R, b, state.get_mod(20));
      add(R, t, a, b);
      mul(R, d, t, t);
      mul(R, e, a, a);
      sqr(R, t, b);
      add(R, e, e, t);
      mul(R, t, a, b);
      add(R, e, e, t);
      add(R, e, e, t);
      TEST(equal(R, d, e), with(R, d), with(R, e));

      random_of_length(state, R, a, state.get_mod(20));
      random_of_length(state, R, b, state.get_mod(20));
      sub(R, t, a, b);
      mul(R, d, t, t);
      mul(R, e, a, a);
      mul(R, t, b, b);
      add(R, e, e, t);
      mul(R, t, a, b);
      sub(R, e, e, t);
      sub(R, e, e, t);
      TEST(equal(R, d, e), with(R, d), with(R, e));
   }

TESTSET_END()

TESTSET_BEGIN(typeid(typename RingT::poly_type).name(), ": div mod gcd ........")

   for (uword i = 100*test_multiplier; i > 0; i--)
   {
      random_of_bits(state, ZZ, m, 2 + state.get_mod(300));
      R.set_modulus(m);

      for (uword j = 15; j > 0; j--)
      {
         try {
            random_of_length(state, R, a, state.get_mod(20));
            random_of_length(state, R, b, state.get_mod(20));
            divrem(R, q, r, a, b);
            mul(R, c, q, b);
            add(R, c, c, r);
            TEST(equal(R, a, c), with(R, q), with(R, r), with(R, a), with(R, b));
            TEST(b.degree() < 0 || r.degree() < b.degree());

            random_of_length(state, R, a, state.get_mod(10));
            random_of_length(state, R, b, state.get_mod(10));
            for (uword k = state.get_mod(10); k > 0; k--)
            {
               random_of_length(state, R, q, state.get_mod(10));
               mul(R, f, b, q);
               add(R, a, a, f);
               std::swap(a, b);
            }
            random_of_length(state, R, c, state.get_mod(5));
            set(R, s, a);
            set(R, t, b);
            mul(R, a, s, c);
            mul(R, b, t, c);

            gcd(R, g, a, b);
            TEST(is_canonical(R, g), with(R, g), with(R, a), with(R, b));

//ha we have gcd(s*c, t*c) = g in ZZ/4ZZ with
//    g = x^2 + 3*x + 2
//    s = 2*x^7 + 3*x^6 + 2*x^4 + 3*x + 3
//    t = 2*x^7 + 3*x^6 + 2*x^5 + 3*x^4 + x^3 + x^2 + 3*x + 2
//    c = 2*x^3 + 3*x^2 + 3*x + 2
            TEST(divides(R, q, g, c), R, with(R, g), with(R, c));
            TEST(divides(R, q, a, g));
            TEST(divides(R, q, b, g));

            fp_poly d1, s1;
            gcdx(R, d, s, t, a, b);
            TEST(equal(R, g, d), with(R, g), with(R, d));
            gcdinv(R, d1, s1, a, b);
            TEST(equal(R, d, d1), with(R, d), with(R, d1));
            TEST(equal(R, s, s1), with(R, s), with(R, s1));
            mul(R, e, s, a);
            mul(R, f, t, b);
            add(R, e, e, f);
            TEST(equal(R, e, d), with(R, e), with(R, d));
            TEST(a.degree() < 0 || s.degree() < 0 || t.degree() + g.degree() < a.degree(), with(R, g), with(R, s), with(R, t), with(R, a), with(R, b));
            TEST(b.degree() < 0 || t.degree() < 0 || s.degree() + g.degree() < b.degree(), with(R, g), with(R, s), with(R, t), with(R, a), with(R, b));

            if (a.length() < b.length())
               swap(R, a, b);

            if (a.length() > b.length())
            {
               fp_poly m11, m12, m21, m22;

               int sgnM = hgcd(R, m11, m12, m21, m22, s, t, a, b);

               mul(R, c, m11, s);
               mul(R, d, m12, t);
               add(R, c, c, d);
               mul(R, d, m21, s);
               mul(R, e, m22, t);
               add(R, d, d, e);
               TEST(equal(R, c, a), with(R, c), with(R, a));
               TEST(equal(R, d, b), with(R, d), with(R, b));

               mul(R, c, m11, m22);
               mul(R, d, m12, m21);
               sub(R, c, c, d);
               TEST(equal(R, c, sgnM), with(R, c), sgnM);

               // deg(s) >= ceil(deg(a)/2) > deg(t)
               TEST(s.length() > a.length()/2, s.length(), a.length());
               TEST(a.length()/2 >= t.length(), a.length(), t.length());
            }
         }
         catch (const parent_and_elem<fp_ring, fmpz>& ex)
         {
            TEST(&ex.parent() == &R, &ex.parent(), &R);
            TEST(fmpz_cmp_si(ex.elem(), 1) > 0, ex.elem());
            TEST(fmpz_cmp(R.characteristic(), ex.elem()) > 0, R.characteristic(), ex.elem());
            TEST(fmpz_divisible(R.characteristic(), ex.elem()), R.characteristic(), ex.elem());
         }
      }
   }

TESTSET_END()

TESTSET_BEGIN(typeid(typename RingT::poly_type).name(), ": bma ................")

   fp_berlekamp_massey B1, B2;

   for (uword i = 0; i < 20*test_multiplier; i++)
   {
      random_of_bits(state, ZZ, m, 1 + state.get_mod(300));
      fmpz_next_prime(m, m, false);
      R.set_modulus(m);

      B1.start_over(R);
      B2.start_over(R);

      // check intermediate reductions match
      for (uword k = 0; k < 10; k++)
      {
         B1.add_point_ui(R, state.get_limb());
         B1.add_zeros(R, state.get_mod(5));
         B1.add_point_ui(R, state.get_limb());
         if (state.get_bit())
            B1.reduce(R);
      }

      for (uword k = 0; k < B1.point_count(); k++)
         B2.add_point(R, B1.point(k, R.stride()));
      B2.reduce(R);
      B1.reduce(R);
      TEST(equal(R, B1.v_poly(), B2.v_poly()), R, with(R, B1.v_poly()), with(R, B2.v_poly()));

      //  Check berlekamp-massey does its job - 2k coeffcients of
      //      u     a1    a2          a(2k)
      //     --- = --- + --- + ... + ------ + ...
      //      v     x    x^2         x^(2k)
      //  should be sufficient to reconstruct a divisor of v
      for (uword k = 1; k < 15; k++)
      {
         fp_poly u, v;

         // deg(u) < deg(v), deg(v) = k
         random_upto_length(state, R, u, k);
         random_monic_of_length(state, R, v, k + 1);

         // q has enough coefficients of expansion of u/v at infty
         shift_left(R, s, u, 2*k);
         divrem(R, q, r, s, v);
         B1.start_over(R);
         for (uword l = 2*k; l > 0; l--)
         {
            if (l-1 < q.length())
               B1.add_point(R, q.coeff(l-1, R.stride()));
            else
               B1.add_zeros(R, 1);
         }

         B1.reduce(R);
         TEST(divides(R, q, v, B1.v_poly()), R, with(R, v), with(R, B1.v_poly()));
         TEST(B1.r_poly().degree() < B1.v_poly().degree());
      }
   }

TESTSET_END()

TESTSET_BEGIN(typeid(typename RingT::poly_type).name(), ": factor .............")

   for (uword i = 0; i < 3*test_multiplier; i++)
   {
      random_of_bits(state, ZZ, m, 2 + state.get_mod(200));
      fmpz_next_prime(m, m, false);
      R.set_modulus(m);

      random_of_length(state, R, a, 2 + state.get_mod(8));
      uword omega_min = 1;
      while (a.length() < 25)
      {
         random_of_length(state, R, b, 2 + state.get_mod(8));
         uword k = 1 + state.get_mod(3);
         omega_min += k;
         pow(R, c, b, k);
         mul(R, a, c);
         test_factorization<RingT>(R, a, omega_min);
      }
   }

TESTSET_END()

TESTSET_BEGIN(typeid(typename RingT::poly_type).name(), ": input output .......")

   R.set_modulus(11);

   random_of_length(state, R, a, 5);
   random_of_length(state, R, b, 4);
   pow(R, c, a, 2);
   pow(R, d, b, 3);
   mul(R, q, c, d);
   test_squarefree_factorization(R, q, 2+3);

   set(R, q, "(x+1)^4*(x+2)^5");
   TEST(!is_irreducible(R, q), with(R, q));
   test_squarefree_factorization<RingT>(R, q, 9);

   set(R, q, "(x+1)*(x^2+1)");
   TEST(!is_irreducible(R, q), with(R, q));
   test_factorization<RingT>(R, q, 2);

   set(R, q, "(x+1)*(x+2)*(x^2+8*x+10)*(x^2+10*x+8)*(x^2+7*x+5)");
   TEST(!is_irreducible(R, q), with(R, q));
   test_factorization<RingT>(R, q, 5);

   set(R, q, "(x^2+8*x+10)^2");
   TEST(!is_irreducible(R, q), with(R, q));

   set(R, q, "(x^2+8*x+10)");
   TEST(is_irreducible(R, q), with(R, q));

   R.set_modulus(2);
   set(R, q, "x^3*(x+1)^4*(x^2+x+1)^5*(x^3+x+1)^6");
   test_factorization<RingT>(R, q, 3+4+5+6);

   for (uword i = 0; i < 100*test_multiplier; i++)
   {
      std::stringstream o;
      random_of_bits(state, ZZ, m, 2 + state.get_mod(1000));
      R.set_modulus(m);
      random_of_length(state, R, a, state.get_mod(20));
      o << with(R, a);
      set(R, b, o.str().c_str());
      TEST(equal(R, a, b), with(R, a), with(R, b));
   }
   try {
      set(R, q, "x+");
      TEST(false && "unreachable");
   }
   catch (const char*& oops)
   {
      TEST(strcmp(oops, "syntax_error") == 0);
   }

TESTSET_END()
}


void test_fmpz(random_state& state)
{
   fmpz a, b, c, d, e, s, t, q, r;

TESTSET_BEGIN("fmpz", ": add sub mul ........")

   for (uword max_bits = 200; max_bits <= ui_pow2(20); max_bits += 1 + max_bits/2)
   for (uword i = 0; i < 1+200000*test_multiplier/max_bits; i++)
   {
      random_of_bits(state, ZZ, a, state.get_mod(max_bits));
      random_of_bits(state, ZZ, b, state.get_mod(max_bits));
      add(ZZ, c, a, b);
      sub(ZZ, d, c, b);
      TEST(equal(ZZ, d, a));

      set(ZZ, e, a);
      add(ZZ, a, a, b);
      sub(ZZ, a, a, b);
      TEST(equal(ZZ, e, a));

      sub(ZZ, a, c, b);
      sub(ZZ, c, b, c);
      neg(ZZ, c, c);
      TEST(equal(ZZ, a, c));

      zero(ZZ, a);
      random_of_bits(state, ZZ, b, state.get_mod(max_bits));
      add(ZZ, c, a, b);
      sub(ZZ, d, c, b);
      TEST(is_zero(ZZ, d));
      TEST(is_zero(ZZ, a));

      set(ZZ, e, a);
      add(ZZ, a, a, b);
      sub(ZZ, a, a, b);
      TEST(equal(ZZ, e, a));

      sub(ZZ, e, c, b);
      sub(ZZ, t, b, c);
      neg(ZZ, t);
      TEST(equal(ZZ, t, e));

      random_of_bits(state, ZZ, a, state.get_mod(max_bits));
      random_of_bits(state, ZZ, b, state.get_mod(max_bits));
      random_of_bits(state, ZZ, c, state.get_mod(max_bits));
      add(ZZ, t, b, c);
      mul(ZZ, d, a, t);
      mul(ZZ, t, a, b);
      mul(ZZ, s, a, c);

      add(ZZ, e, t, s);
      TEST(equal(ZZ, d, e), d, e);

      random_of_bits(state, ZZ, a, state.get_mod(max_bits));
      random_of_bits(state, ZZ, b, state.get_mod(max_bits));
      add(ZZ, t, a, b);
      sqr(ZZ, d, t);
      mul(ZZ, e, a, a);
      mul(ZZ, t, b, b);
      add(ZZ, e, t);
      mul(ZZ, t, a, b);
      add(ZZ, e, t);
      add(ZZ, e, t);
      TEST(equal(ZZ, d, e), d, e);

      random_of_bits(state, ZZ, a, state.get_mod(max_bits));
      random_of_bits(state, ZZ, b, state.get_mod(max_bits));
      sub(ZZ, t, a, b);
      sqr(ZZ, d, t);
      mul(ZZ, e, a, a);
      sqr(ZZ, t, b);
      add(ZZ, e, t);
      mul(ZZ, t, a, b);
      sub(ZZ, e, t);
      sub(ZZ, e, t);
      TEST(equal(ZZ, d, e), d, e);
   }

TESTSET_END()

#if 1
TESTSET_BEGIN("fmpz", ": div mod ............")

   for (uword max_bits = 200; max_bits <= ui_pow2(19); max_bits += 1 + max_bits/2)
   for (uword i = 0; i < 1+200000*test_multiplier/max_bits; i++)
   {
      random_of_bits(state, ZZ, a, state.get_mod(max_bits));
      random_of_bits(state, ZZ, b, state.get_mod(max_bits)+1);

      mul(ZZ, c, a, b);
      divexact(ZZ, d, c, b);
      TEST(equal(ZZ, a, d), a, d);

      mul(ZZ, c, a, b);
      divexact(ZZ, c, c, b);
      TEST(equal(ZZ, a, c), a, c);

      fmpz_tdiv_qr(q, r, a, b);
      TEST(fmpz_cmpabs(r, b) < 0);
      mul(ZZ, t, q, b);
      add(ZZ, t, t, r);
      TEST(equal(ZZ, t, a), t, a);
      TEST(r.sign()*a.sign() >= 0, r, a, b);

      fmpz_fdiv_qr(q, r, a, b);
      TEST(fmpz_cmpabs(r, b) < 0, r, b);
      mul(ZZ, t, q, b);
      add(ZZ, t, t, r);
      TEST(equal(ZZ, t, a), q, r, a, b);
      TEST(r.sign()*b.sign() >= 0, r, a, b);

      fmpz_cdiv_qr(q, r, a, b);
      TEST(fmpz_cmpabs(r, b) < 0);
      mul(ZZ, t, q, b);
      add(ZZ, t, t, r);
      TEST(equal(ZZ, t, a), t, a);
      TEST(r.sign()*b.sign() <= 0, q, r, a, b);
   }

TESTSET_END()

TESTSET_BEGIN("fmpz", ": gcd ................")

   for (uword max_bits = 200; max_bits <= ui_pow2(16); max_bits += 1 + max_bits/2)
   for (uword i = 0; i < 1 + 200000*test_multiplier/max_bits; i++)
   {
      if (i % 2)
      {
         random_of_bits(state, ZZ, a, state.get_mod(max_bits));
         random_of_bits(state, ZZ, b, state.get_mod(max_bits));
         random_of_bits(state, ZZ, e, state.get_mod(max_bits));
         mul_2exp(ZZ, a, a, state.get_mod(100));
         mul_2exp(ZZ, b, b, state.get_mod(100));
         mul_2exp(ZZ, e, e, state.get_mod(100));
         mul(ZZ, a, a, e);
         mul(ZZ, b, b, e);
      }
      else
      {
         random_of_bits(state, ZZ, e, 1 + state.get_mod(20 + max_bits/20));
         zero(ZZ, a);
         set(ZZ, b, e);
         while (nbits(a) < max_bits && nbits(b) < max_bits)
         {
            random_of_bits(state, ZZ, c, state.get_mod(20 + max_bits/20));
            mul(ZZ, d, c, b);
            add(ZZ, a, d);
            swap(ZZ, a, b);
         }
      }

      gcd(ZZ, d, a, b);
      TEST(divisible(ZZ, a, d), a, b, d);
      TEST(divisible(ZZ, b, d), a, b, d);
      TEST(divisible(ZZ, d, e), d, e);
      if (!is_zero(ZZ, d))
      {
         divexact(ZZ, a, a, d);
         divexact(ZZ, b, b, d);
         gcd(ZZ, e, a, b);
         TEST(is_one(ZZ, e), e, a, b);
      }

      random_of_bits(state, ZZ, a, state.get_mod(max_bits));
      random_of_bits(state, ZZ, b, state.get_mod(max_bits));
      random_of_bits(state, ZZ, c, state.get_mod(max_bits));
      random_of_bits(state, ZZ, e, state.get_mod(max_bits));
      mul(ZZ, a, a, e);
      mul(ZZ, b, b, e);
      mul(ZZ, c, c, e);
      gcd(ZZ, d, a, b, c);
      TEST(divisible(ZZ, a, d), a, b, c, d);
      TEST(divisible(ZZ, b, d), a, b, c, d);
      TEST(divisible(ZZ, c, d), a, b, c, d);
      TEST(divisible(ZZ, d, e), d, e);
      if (!is_zero(ZZ, d))
      {
         divexact(ZZ, a, a, d);
         divexact(ZZ, b, b, d);
         divexact(ZZ, c, c, d);
         gcd(ZZ, e, a, b, c);
         TEST(is_one(ZZ, e), e, a, b, c);
      }
   }

   for (uword max_bits = 20; max_bits <= ui_pow2(15); max_bits += 1 + max_bits/2)
   for (uword i = 0; i < 1 + 20000*test_multiplier/max_bits; i++)
   {
      if (i % 2)
      {
         random_of_bits(state, ZZ, a, state.get_mod(max_bits));
         random_of_bits(state, ZZ, b, state.get_mod(max_bits));
         random_of_bits(state, ZZ, e, state.get_mod(max_bits));
         mul_2exp(ZZ, a, a, state.get_mod(100));
         mul_2exp(ZZ, b, b, state.get_mod(100));
         mul_2exp(ZZ, e, e, state.get_mod(100));
         mul(ZZ, a, a, e);
         mul(ZZ, b, b, e);
      }
      else
      {
         random_of_bits(state, ZZ, e, 1 + state.get_mod(20 + max_bits/20));
         zero(ZZ, a);
         set(ZZ, b, e);
         while (nbits(a) < max_bits && nbits(b) < max_bits)
         {
            random_of_bits(state, ZZ, c, state.get_mod(20 + max_bits/20));
            mul(ZZ, d, c, b);
            add(ZZ, a, a, d);
            swap(ZZ, a, b);
         }
      }
      gcdx(ZZ, d, s, t, a, b);
      TEST(divisible(ZZ, b, d), b, d);
      TEST(divisible(ZZ, b, d), d, e);
      mul(ZZ, q, a, s);
      mul(ZZ, r, b, t);
      add(ZZ, q, r);
      TEST(equal(ZZ, q, d), d, s, t, a, b);
   }

TESTSET_END()

TESTSET_BEGIN("fmpz", ": primes .............")

   for (uword i = 0; i < 5000*test_multiplier; i++)
   {
      random_of_bits(state, ZZ, a, state.get_mod(1500)+2);
      random_of_bits(state, ZZ, b, state.get_mod(1500)+2);
      mul(ZZ, c, a, b);
      TEST(!fmpz_is_probable_prime(c), a, b);
   }

TESTSET_END()

TESTSET_BEGIN("fmpz", ": input output .......")

   for (uword bits = 2; bits < 100000*test_multiplier; bits += 1 + bits/16)
   {
      std::stringstream o;
      random_of_bits(state, ZZ, a, bits);
      o << a;
      b.set_c_str(o.str().c_str());
      TEST(equal(ZZ, a, b), a, b);
   }

TESTSET_END()
#endif
}



extern "C" {
   std::complex<double> my_cis2pi_frac(sword a, sword b);
}


/*
cnorm[x_]:=1/2(1+Erf[x/Sqrt[2]]);
relerr[a_,b_]:=Abs[a-b]/Max[Abs[a],Abs[b]];
<< FunctionApproximations`

f[x_] = x*(mmlist[[2,1]]/.x->x^2);
Table[relerr[cnorm[x],f[x]],{x,N[-1,30],N[1,30],2^-6}]


mmlist = MiniMaxApproximation[(cnorm[Sqrt[x]]-1/2)/Sqrt[x], {x, {10^-15, 1}, 4, 4},WorkingPrecision->30]
{{1.00000000000000000000000000000*10^-15,0.0292305695345921437251722626141,0.113699874277154120922294311163,0.244001062428549069456872698306,0.405339755094518381530192301552,0.578902572037679744280411286294,0.743810966720942619462477999261,0.879583839221869181529353695684,0.968867447130072264760090449468,1.00000000000000000000000000000},
{(0.398942280401432671967506317593+0.0242770520839500197327775480615 x+0.00357628924737828104575497487818 x^2+0.0000550328822551643311115885318363 x^3+1.43162831811074734139722489263*10^-6 x^4)
/(1+0.227520211844982736553773950853 x+0.0218844630533032192598007160720 x^2+0.00107354266691392553483307336620 x^3+0.0000231924011957566413906256612629 x^4),1.49706863266813104570786526893*10^-17}}

mmlist = MiniMaxApproximation[cnorm[x], {x, {1,Sqrt[5]}, 7, 7},WorkingPrecision->30]
{{1.00000000000000000000000000000,1.01226538075438489495726818470,1.04863194599927966147491951182,1.10782192806509219372656841499,1.18774196606034808890716762627,1.28553543572287505994943198621,1.39764924574690640856351839870,1.51990680663218916112702590167,1.64758023200528589206533879549,1.77546395332493843603142825199,1.89796872172554254609994977898,2.00927386533422053730535214338,2.10358507641815695917899445322,2.17552992775548752848124968119,2.22067409838346254865729471281,2.23606797749978969640917366873},{(0.499998824085616158823924205124+0.520103823175521671596495580219 x+0.101136322156088087077718146739 x^2-0.00358188259512522426332825241504 x^3+0.0174175756030943271926816097366 x^4+0.00609428810670336664837911010447 x^5-0.00119839875863545361524170790446 x^6+0.000394653238766413875205839165013 x^7)/(1+0.242299585286215558640526602925 x+0.00905368589951736524399529640027 x^2+0.118290204612757364569152284230 x^3-0.0267504298291202772494479768904 x^4+0.0140120531877952910805626977733 x^5-0.00191846034224973269770278741823 x^6+0.000421129927632370081710404825236 x^7),2.92560466454817734663603997048*10^-17}}

mmlist = MiniMaxApproximation[(1-cnorm[x])/Exp[-1/2*x^2], {x, {Sqrt[5],8+1/2}, 6, 6},WorkingPrecision->30]
{{2.23606797749978969640917366873,2.28032484250310128017592116307,2.41436813373067240141794630668,2.64193536017794786773939523005,2.96889598908902154952821207851,3.40236524581944588468198474138,3.94868442750833422467088332918,4.60936828231433557213991613127,5.37400305515510463753581802774,6.20985831191355236629632440698,7.05081607207042260469860263704,7.79351233954416294880267997466,8.31256269500075408010931182485,8.50000000000000000000000000000},{(0.499995747012825136822562135422+0.535353678313639053651989655839 x+0.278178003821268015136607175806 x^2+0.0842340709824248563578004218590 x^3+0.0149048904068757366053385321737 x^4+0.00131352116385916283693265841370 x^5+3.52093296442626990311440433034*10^-11 x^6)/(1+1.86854689841039619802245759303 x+1.54735700543080909827822808509 x^2+0.734572265583459906820580853880 x^3+0.214441130378615430811179389473 x^4+0.0373607854718780428236728089347 x^5+0.00329251594513778102306435732839 x^6),5.62057648441586986977096568556*10^-16}}

mmlist = MiniMaxApproximation[cnorm[x], {x, {-2,-1}, 7, 7},WorkingPrecision->30]
{{-2.00000000000000000000000000000,-1.98968239481199984708797064818,-1.95910366735360255530328335902,-1.90938321761777172014351395834,-1.84237449528240821641509823048,-1.76063831060068563448074886989,-1.66739185956488525961031839700,-1.56642576872169791071047675650,-1.46198478156410438039954964064,-1.35861157726505337484285489029,-1.26095642361998356603220871157,-1.17355803140074187311320739970,-1.10060435293087663672804439402,-1.04568721063632819683450522673,-1.01157093073841770107369592862,-1.00000000000000000000000000000},{(0.500000479674713047960298217499+0.364162299972498300817899551200 x+0.0200367443156036883575726226228 x^2-0.0595001625243434310079460657655 x^3-0.0250900635425580417406651282947 x^4-0.00452483787507044253887217131538 x^5-0.000386548072424271473474770623061 x^6-0.0000122882789027113745235406971513 x^7)/(1-0.0695699209107777662730855933368 x+0.0955345531696435822533871927313 x^2-0.0623844858268367144291148549564 x^3-0.00993253616994806137996448403858 x^4-0.00875937890802753513330020889159 x^5-0.00109672645166158707718814778213 x^6-0.000303535021920076703631358404870 x^7),2.16942545691135755327817794616*10^-17}}

mmlist = MiniMaxApproximation[(cnorm[x]), {x, {-3,-2}, 7, 7},WorkingPrecision->40]
{{-3.000000000000000000000000000000000000000,-2.990430979143158369736446998754545554583,-2.961988185479236514179096131310140910649,-2.915474451173349930490316947101804676355,-2.852258371264448783131648103746128779120,-2.774310293273339108845887077463979260862,-2.684234047067025324371610769714323020547,-2.585278294438305816193604664728831022342,-2.481310524025581632264928863155399576514,-2.376738112520520190472981126378572363211,-2.276363895684003555114853018443023120460,-2.185167930420281333913012222307098323101,-2.108013836816493305537870856994419877492,-2.049290611300419874447205470255159133216,-2.012521811932646656040319259036293363635,-2.000000000000000000000000000000000000000},
{(0.5000272052634430144374276382858466005875+0.7026946614641770915093512445773582055591 x+0.4292886404379185678042417596888852436654 x^2+0.1478070866573307398101863446173223260535 x^3+0.03097898176022996582758832827879041741233 x^4+0.003952693846305842985569833648500469771112 x^5+0.0002842927236327892172729502902811384453710 x^6+8.891678849330812012775906582645614578425*10^-6 x^7)
/(1+0.6073091279263267564355511834691489080986 x+0.3738693389633952304114169372025940446839 x^2+0.1307706943994043226995854142018856210128 x^3+0.03981275781123116329028630949077172752658 x^4+0.007814252419509082154700686175222777497832 x^5+0.001133295729212619404972186472699704009562 x^6+0.00007646283706930613856477617723455695355712 x^7),4.149517202893365702278797654488761053387*10^-17}}

mmlist = MiniMaxApproximation[(cnorm[x])/Exp[-1/2*x^2], {x, {-9,-3}, 6,6},WorkingPrecision->60]
{{-9.00000000000000000000000000000000000000000000000000000000000,-8.83272307481869700967292936279944049508322652939218277435093,-8.36576785227670631724902405160664646959859745574917754251699,-7.68786563380989233896740269092679904511377277489481510393791,-6.90601838731266846070062785929275057781487347751072508154618,-6.11349588475585701046506219758629687837147551182646070870261,-5.37479076280326181002536881688424717662154458379369066728325,-4.72592297042395370079799472377344341513139080258952179740089,-4.18205756694651905533615101269899539933408310617006950069899,-3.74595724096216178040465132968384752065955146422165753695445,-3.41443822986720554153739966199321898571258320740510908382040,-3.18244952087179838899489513367046993888600824044856405395819,-3.04534153203908047344870640688118105479653698030748331550524,-3.00000000000000000000000000000000000000000000000000000000000},
{(0.499980970297891012776749712977838328496094337510055847749244-0.555967971027661368474134336778404009336191082594056817229524 x+0.298658083351315058744560418384265283586597620363916744578155 x^2-0.0939854270027801368801692159770319023307085194739260635508334 x^3+0.0173886053360974089999868730925341461382565017381597663180374 x^4-0.00166232626990265084261588858214248132257947887767015473250521 x^5+1.35679311883112384685161517635830988730636952158465779467859*10^-11 x^6)
/(1-1.90964129559407800511844764169197382265857045602560703113695 x+1.62140601258396136906870801866431493225103252430218351740549 x^2-0.792166850997619519046075020567612684399560237199538261984232 x^3+0.239756037982848283322999777283276814105425766461280147083825 x^4-0.0435866589598823234670435491216236933889123041472803736637494 x^5+0.00416683687383734231335071931372141986591214883811483182410396 x^6),-4.57324066688671939673345545620476167414577159660296713222945*10^-17}}

mmlist = MiniMaxApproximation[(x*cnorm[x])/Exp[-1/2*x^2], {x, {-38,-9}, 6,5},WorkingPrecision->60] 
{{-38.0000000000000000000000000000000000000000000000000000000000,-36.3737882920588509310028020883792819286710711172038838466131,-32.2415153190794303597522449233477669141466325845929225374920,-27.1528066702834448214027840726546554766064162367974188121051,-22.3508884850264598229375016066237156767140205727557375518317,-18.3893553379911910285396221286622712766397313009554094210077,-15.3414534888589626847909261624867013230979215896109169699812,-13.0828652992180740380071583316711951018716698608336265480106,-11.4524294447069329140514672805596621156294386393614080072629,-10.3126505939447950992948018536782399218104488884497276775449,-9.56308549050447756026020538430176201754143874245042701339950,-9.13781834504231768775157749620461330905384372144565338162081,-9.00000000000000000000000000000000000000000000000000000000000},
{(-1.42298670320543621309264026578670165724599032878166817867300-3.70185571706124598750538863763052890284479839964569901338729 x+1.99547061425098914745484654827013948497994492428961358101670 x^2-3.33505056618593269877474917288060675240517742893506426639096 x^3+0.487585214685351266418297775658826893033648074042177240236148 x^4-0.302182040872990143277697925092467197085062531068371805669341 x^5-1.11579873338046773584051341302790196357393408839588091301281*10^-13 x^6)/(1+16.1235235012955598214886090344731803127527217056695185864888 x-6.22411460326542449336535809630051850603507966206337269066680 x^2+9.11718968298034514568654219520594031832811684729394675436631 x^3-1.22219489212309550970637204218357629245567424948836797113933 x^4+0.757458047673258365521143965662215195327527607555919905806763 x^5),3.23791534228809291134937316803378873326080170602256986955273*10^-17}} 



$MinMachineNumber 2.22e-308

gaussint[-37] = 5.72e-300
gaussint[-38] = 2.88e-316
*/



template<typename T, int n, int stride>
T horner_eval_rec(const T* cs, T c, T x) {
   if constexpr (n < 1)
      return c;
   else
      return horner_eval_rec<T,n-1,stride>(cs + stride, fmadd(x, c, cs[0]), x);
}

template<typename T, int n, int stride>
T horner_eval(const T* cs, T x) {
   if constexpr(n < 1)
      return T(0);
   else
      return horner_eval_rec<T,n-1,stride>(cs+stride, cs[0], x);
}

template<typename T, int n>
T horner(const T* cs, T x) {
   return horner_eval<T,n,1>(cs, x);
}

template<typename T, int n>
T horner2(const T* cs, T x, T x2) {
   T evenpart, oddpart;
   if constexpr (n%2 == 0) {
      evenpart = horner_eval<T,n/2,2>(cs+1, x2);
       oddpart = horner_eval<T,n/2,2>(cs+0, x2);
   } else {
      evenpart = horner_eval<T,(n+1)/2,2>(cs+0, x2);
       oddpart = horner_eval<T,(n-1)/2,2>(cs+1, x2);
   }
   return fmadd(oddpart, x, evenpart);
}

static packed<double,2> T1[5] = {
   packed<double,2>(1.43162831811074734139722489263E-6   , 0.0000231924011957566413906256612629),
   packed<double,2>(0.0000550328822551643311115885318363 , 0.00107354266691392553483307336620),
   packed<double,2>(0.00357628924737828104575497487818   , 0.0218844630533032192598007160720),
   packed<double,2>(0.0242770520839500197327775480615    , 0.227520211844982736553773950853),
   packed<double,2>(0.398942280401432671967506317593     , 1.0),
};

static packed<double,2> T2[8] = {
   packed<double,2>(0.000394653238766413875205839165013 , 0.000421129927632370081710404825236),
   packed<double,2>(-0.00119839875863545361524170790446 , -0.00191846034224973269770278741823),
   packed<double,2>(0.00609428810670336664837911010447  , 0.0140120531877952910805626977733),
   packed<double,2>(0.0174175756030943271926816097366   , -0.0267504298291202772494479768904),
   packed<double,2>(-0.00358188259512522426332825241504 , 0.118290204612757364569152284230),
   packed<double,2>(0.101136322156088087077718146739    , 0.00905368589951736524399529640027),
   packed<double,2>(0.520103823175521671596495580219    , 0.242299585286215558640526602925),
   packed<double,2>(0.499998824085616158823924205124    , 1.0),
};

static packed<double,2> T3[7] = {
   packed<double,2>(-3.52093296442626990311440433034e-11 , 0.00329251594513778102306435732839),
   packed<double,2>(-0.00131352116385916283693265841370  , 0.0373607854718780428236728089347),
   packed<double,2>(-0.0149048904068757366053385321737   , 0.214441130378615430811179389473),
   packed<double,2>(-0.0842340709824248563578004218590   , 0.734572265583459906820580853880),
   packed<double,2>(-0.278178003821268015136607175806    , 1.54735700543080909827822808509 ),
   packed<double,2>(-0.535353678313639053651989655839    , 1.86854689841039619802245759303 ),
   packed<double,2>(-0.499995747012825136822562135422    , 1.0 ),
};

static packed<double,2> T4[8] = {
   packed<double,2>(-0.0000122882789027113745235406971513 , -0.000303535021920076703631358404870   ),
   packed<double,2>(-0.000386548072424271473474770623061  , -0.00109672645166158707718814778213   ),
   packed<double,2>(-0.00452483787507044253887217131538   , -0.00875937890802753513330020889159   ),
   packed<double,2>(-0.0250900635425580417406651282947    , -0.00993253616994806137996448403858   ),
   packed<double,2>(-0.0595001625243434310079460657655    , -0.0623844858268367144291148549564   ),
   packed<double,2>(0.0200367443156036883575726226228     , 0.0955345531696435822533871927313   ),
   packed<double,2>(0.364162299972498300817899551200      , -0.0695699209107777662730855933368   ),
   packed<double,2>(0.500000479674713047960298217499      , 1.0   ),
};


static packed<double,2> T5[8] = {
   packed<double,2>(8.891678849330812012775906582645614578425E-6 , 0.00007646283706930613856477617723455695355712),
   packed<double,2>(0.0002842927236327892172729502902811384453710   , 0.001133295729212619404972186472699704009562),
   packed<double,2>(0.003952693846305842985569833648500469771112    , 0.007814252419509082154700686175222777497832),
   packed<double,2>(0.03097898176022996582758832827879041741233     , 0.03981275781123116329028630949077172752658),
   packed<double,2>(0.1478070866573307398101863446173223260535      , 0.1307706943994043226995854142018856210128),
   packed<double,2>(0.4292886404379185678042417596888852436654      , 0.3738693389633952304114169372025940446839),
   packed<double,2>(0.7026946614641770915093512445773582055591      , 0.6073091279263267564355511834691489080986 ),
   packed<double,2>(0.5000272052634430144374276382858466005875      , 1.0    ),
};

// 1.49e-18 
static packed<double,2> T1_to_2[8] = {
   packed<double,2>(0.0002758945479707823974, 0.00028593403288013499834),
   packed<double,2>(-0.00063856177824806612406, -0.00097362459574962795361),
   packed<double,2>(0.0039717710448780453451, 0.007984764480314199896),
   packed<double,2>(0.014816795195927067171, -0.0074503836005480758024),
   packed<double,2>(0.0038487708093279040159, 0.054262967863203988667),
   packed<double,2>(0.071985020985342096307, 0.10815812074256124728),
   packed<double,2>(0.42141310785043729006, 0.044933445803134179322),
   packed<double,2>(0.4999996092149238051, 1.),
};

// 2.4967e-19 
static packed<double,2> T2_to_3[8] = {
   packed<double,2>(-0.000073344567813344635627, -0.000082556195268299153039),
   packed<double,2>(0.00087842004234247627455, 0.0011714311430496478821),
   packed<double,2>(-0.0039952625182962965542, -0.0080484670321535666073),
   packed<double,2>(0.0088724353456808434303, 0.040478922376938874875),
   packed<double,2>(0.017020688323489819477, -0.13302378749671349753),
   packed<double,2>(-0.056572712040268286385, 0.37703465981196528762),
   packed<double,2>(0.091763112669811394042, -0.61445356031450586951),
   packed<double,2>(0.49998268514587597539, 1.),
};

// -3.3979e-19
static packed<double,2> Te3_to_10[8] = {
   packed<double,2>(1.7270625120325394616e-13, 0.00092712044999684058744),
   packed<double,2>(-0.00036986756431373905668, 0.012291547481018709102),
   packed<double,2>(-0.0049036171207845836766, 0.085249627719192768147),
   packed<double,2>(-0.033639839384918177279, 0.36766686274985971499),
   packed<double,2>(-0.14177368676002865086, 1.0406957380611832602),
   packed<double,2>(-0.38228614772497846262, 1.8998985937802478687),
   packed<double,2>(-0.62587750953312411524, 2.0496679083115780626),
   packed<double,2>(-0.50000274715153076135, 1.),
}; 

// -8.3129e-18
static packed<double,2> Te2_to_9[8] = {
   packed<double,2>(7.9288661922271320799e-13, 0.00063925795620326458577),
   packed<double,2>(-0.00025502709725860926017, 0.0094517495515727622791),
   packed<double,2>(-0.0037706995709418668266, 0.069949482292557856199),
   packed<double,2>(-0.027650856328813710824, 0.31815591097938899294),
   packed<double,2>(-0.12315371480831391442, 0.94075373373706946421),
   packed<double,2>(-0.34818537217806044599, 1.783690561741582071),
   packed<double,2>(-0.59576667926715572245, 1.9894206779315979699),
   packed<double,2>(-0.5000002250079655289, 1.),
};

// 1.947e-19
static packed<double,2> Tem4_to_m1[8] = {
   packed<double,2>(3.6261569275188610059e-11, -0.00026786900228293985465),
   packed<double,2>(0.00010686629378713327388, 0.004940605637456463353),
   packed<double,2>(-0.0019709617320100679539, -0.042783597967720134734),
   packed<double,2>(0.016962278229217860721, 0.22121551923043560142),
   packed<double,2>(-0.086269114416847391547, -0.7291621581375166917),
   packed<double,2>(0.27426343107170390137, 1.5199234514665179614),
   packed<double,2>(-0.52311879900327086896, -1.8441221611951403308),
   packed<double,2>(0.50000000011932526069, 1.),
};

// 2.4026e-21
static packed<double,2> Tem3_to_m1[8] = {
   packed<double,2>(7.4928822150106415609e-11, -0.0002263720710905964461),
   packed<double,2>(0.000090313168938319448292, 0.0043404356288538093427),
   packed<double,2>(-0.0017314899702133361243, -0.038744991327229186074),
   packed<double,2>(0.015368207097106376408, 0.20547431251924899545),
   packed<double,2>(-0.080223177506494775524, -0.69222933579026740944),
   packed<double,2>(0.2611336237897270574, 1.4709526082241171413),
   packed<double,2>(-0.50888664234973558816, -1.8156578459532313573),
   packed<double,2>(0.50000000002063121704, 1.),
};

// 3.3979e-19
static packed<double,2> Tem10_to_m3[8] = {
   packed<double,2>(1.7270625120325394616e-13, -0.00092712044999684058744),
   packed<double,2>(0.00036986756431373905668, 0.012291547481018709102),
   packed<double,2>(-0.0049036171207845836766, -0.085249627719192768147),
   packed<double,2>(0.033639839384918177279, 0.36766686274985971499),
   packed<double,2>(-0.14177368676002865086, -1.0406957380611832602),
   packed<double,2>(0.38228614772497846262, 1.8998985937802478687),
   packed<double,2>(-0.62587750953312411524, -2.0496679083115780626),
   packed<double,2>(0.50000274715153076135, 1.),
};

// -2.125e-21
static packed<double,2> Tem30_to_m10[8] = {
   packed<double,2>(9.2418640140077299593e-17, 0.023680755369163037647),
   packed<double,2>(-0.0094472545485792304555, -0.071213173049780216519),
   packed<double,2>(0.028409945653764268677, 0.47866643319091255304),
   packed<double,2>(-0.18151302367353650272, -0.78762438587225384979),
   packed<double,2>(0.28580673185131337386, 1.9169692458498049354),
   packed<double,2>(-0.60214126208540244431, -1.3000477437020902405),
   packed<double,2>(0.28966495608200491063, 1.),
   packed<double,2>(-0.065205400187472230042, 0),
};


// -9.0817e-21 
static packed<double,2> Tem37_to_m10[8] = {
   packed<double,2>(3.1731014978267615543e-17, 0.023769724171374703516),
   packed<double,2>(-0.0094827479654324240891, -0.065552029294590859925),
   packed<double,2>(0.026151476052855241668, 0.46930576252868136995),
   packed<double,2>(-0.1777431630559614097, -0.71684080433876216784),
   packed<double,2>(0.2598266336853098743, 1.8624728646124899274),
   packed<double,2>(-0.58424133339750590564, -1.1710483804975333045),
   packed<double,2>(0.25966176666432516633, 1.),
   packed<double,2>(-0.075261703001719193474, 0),
};

// -1.0627e-20
static packed<double,2> Tem38_to_m10[8] = {
   packed<double,2>(2.7774903860616780814e-17, 0.023764894891223400864),
   packed<double,2>(-0.009480821361397161786, -0.064863187811687378365),
   packed<double,2>(0.025876668060716849691, 0.46802291038106014868),
   packed<double,2>(-0.1772333057066222612, -0.70840804969281230074),
   packed<double,2>(0.2567372589487876441, 1.8557489246421754028),
   packed<double,2>(-0.58206488534966707045, -1.1559258857135235533),
   packed<double,2>(0.25616824331987456929, 1.),
   packed<double,2>(-0.076442597452770179415, 0),
};

static packed<double,2> Tzero[8] = {
   packed<double,2>(0,0),
   packed<double,2>(0,0),
   packed<double,2>(0,0),
   packed<double,2>(0,0),
   packed<double,2>(0,0),
   packed<double,2>(0,0),
   packed<double,2>(0,0),
   packed<double,2>(0,1),
};
 

#include <math.h>
#include <stdint.h>

#include <math.h>
#include <stdint.h>
#define attribute_hidden

static inline uint64_t
asuint64 (double f)
{
  union
  {
    double f;
    uint64_t i;
  } u = {f};
  return u.i;
}
static inline double
asdouble (uint64_t i)
{
  union
  {
    uint64_t i;
    double f;
  } u = {i};
  return u.f;
}

#define EXP_TABLE_BITS 7
#define EXP_POLY_ORDER 5
extern const struct exp_data
{
  double invln2N;
  double shift;
  double negln2hiN;
  double negln2loN;
  double poly[4]; /* Last four coefficients.  */
  uint64_t tab[2*(1 << EXP_TABLE_BITS)];
} __exp_data attribute_hidden;


#define N (1 << EXP_TABLE_BITS)
#define InvLn2N __exp_data.invln2N
#define NegLn2hiN __exp_data.negln2hiN
#define NegLn2loN __exp_data.negln2loN
#define T __exp_data.tab
#define C2 __exp_data.poly[5 - EXP_POLY_ORDER]
#define C3 __exp_data.poly[6 - EXP_POLY_ORDER]
#define C4 __exp_data.poly[7 - EXP_POLY_ORDER]
#define C5 __exp_data.poly[8 - EXP_POLY_ORDER]

#define N (1 << EXP_TABLE_BITS)
const struct exp_data __exp_data = {
// N/ln2
.invln2N = 0x1.71547652b82fep0 * N,
// -ln2/N
.negln2hiN = -0x1.62e42fefa0000p-8,
.negln2loN = -0x1.cf79abc9e3b3ap-47,
.poly = {
0x1.ffffffffffdbdp-2,
0x1.555555555543cp-3,
0x1.55555cf172b91p-5,
0x1.1111167a4d017p-7,
},
// 2^(k/N) ~= H[k]*(1 + T[k]) for int k in [0,N)
// tab[2*k] = asuint64(T[k])
// tab[2*k+1] = asuint64(H[k]) - (k << 52)/N
.tab = {
0x0, 0x3ff0000000000000,
0x3c9b3b4f1a88bf6e, 0x3feff63da9fb3335,
0xbc7160139cd8dc5d, 0x3fefec9a3e778061,
0xbc905e7a108766d1, 0x3fefe315e86e7f85,
0x3c8cd2523567f613, 0x3fefd9b0d3158574,
0xbc8bce8023f98efa, 0x3fefd06b29ddf6de,
0x3c60f74e61e6c861, 0x3fefc74518759bc8,
0x3c90a3e45b33d399, 0x3fefbe3ecac6f383,
0x3c979aa65d837b6d, 0x3fefb5586cf9890f,
0x3c8eb51a92fdeffc, 0x3fefac922b7247f7,
0x3c3ebe3d702f9cd1, 0x3fefa3ec32d3d1a2,
0xbc6a033489906e0b, 0x3fef9b66affed31b,
0xbc9556522a2fbd0e, 0x3fef9301d0125b51,
0xbc5080ef8c4eea55, 0x3fef8abdc06c31cc,
0xbc91c923b9d5f416, 0x3fef829aaea92de0,
0x3c80d3e3e95c55af, 0x3fef7a98c8a58e51,
0xbc801b15eaa59348, 0x3fef72b83c7d517b,
0xbc8f1ff055de323d, 0x3fef6af9388c8dea,
0x3c8b898c3f1353bf, 0x3fef635beb6fcb75,
0xbc96d99c7611eb26, 0x3fef5be084045cd4,
0x3c9aecf73e3a2f60, 0x3fef54873168b9aa,
0xbc8fe782cb86389d, 0x3fef4d5022fcd91d,
0x3c8a6f4144a6c38d, 0x3fef463b88628cd6,
0x3c807a05b0e4047d, 0x3fef3f49917ddc96,
0x3c968efde3a8a894, 0x3fef387a6e756238,
0x3c875e18f274487d, 0x3fef31ce4fb2a63f,
0x3c80472b981fe7f2, 0x3fef2b4565e27cdd,
0xbc96b87b3f71085e, 0x3fef24dfe1f56381,
0x3c82f7e16d09ab31, 0x3fef1e9df51fdee1,
0xbc3d219b1a6fbffa, 0x3fef187fd0dad990,
0x3c8b3782720c0ab4, 0x3fef1285a6e4030b,
0x3c6e149289cecb8f, 0x3fef0cafa93e2f56,
0x3c834d754db0abb6, 0x3fef06fe0a31b715,
0x3c864201e2ac744c, 0x3fef0170fc4cd831,
0x3c8fdd395dd3f84a, 0x3feefc08b26416ff,
0xbc86a3803b8e5b04, 0x3feef6c55f929ff1,
0xbc924aedcc4b5068, 0x3feef1a7373aa9cb,
0xbc9907f81b512d8e, 0x3feeecae6d05d866,
0xbc71d1e83e9436d2, 0x3feee7db34e59ff7,
0xbc991919b3ce1b15, 0x3feee32dc313a8e5,
0x3c859f48a72a4c6d, 0x3feedea64c123422,
0xbc9312607a28698a, 0x3feeda4504ac801c,
0xbc58a78f4817895b, 0x3feed60a21f72e2a,
0xbc7c2c9b67499a1b, 0x3feed1f5d950a897,
0x3c4363ed60c2ac11, 0x3feece086061892d,
0x3c9666093b0664ef, 0x3feeca41ed1d0057,
0x3c6ecce1daa10379, 0x3feec6a2b5c13cd0,
0x3c93ff8e3f0f1230, 0x3feec32af0d7d3de,
0x3c7690cebb7aafb0, 0x3feebfdad5362a27,
0x3c931dbdeb54e077, 0x3feebcb299fddd0d,
0xbc8f94340071a38e, 0x3feeb9b2769d2ca7,
0xbc87deccdc93a349, 0x3feeb6daa2cf6642,
0xbc78dec6bd0f385f, 0x3feeb42b569d4f82,
0xbc861246ec7b5cf6, 0x3feeb1a4ca5d920f,
0x3c93350518fdd78e, 0x3feeaf4736b527da,
0x3c7b98b72f8a9b05, 0x3feead12d497c7fd,
0x3c9063e1e21c5409, 0x3feeab07dd485429,
0x3c34c7855019c6ea, 0x3feea9268a5946b7,
0x3c9432e62b64c035, 0x3feea76f15ad2148,
0xbc8ce44a6199769f, 0x3feea5e1b976dc09,
0xbc8c33c53bef4da8, 0x3feea47eb03a5585,
0xbc845378892be9ae, 0x3feea34634ccc320,
0xbc93cedd78565858, 0x3feea23882552225,
0x3c5710aa807e1964, 0x3feea155d44ca973,
0xbc93b3efbf5e2228, 0x3feea09e667f3bcd,
0xbc6a12ad8734b982, 0x3feea012750bdabf,
0xbc6367efb86da9ee, 0x3fee9fb23c651a2f,
0xbc80dc3d54e08851, 0x3fee9f7df9519484,
0xbc781f647e5a3ecf, 0x3fee9f75e8ec5f74,
0xbc86ee4ac08b7db0, 0x3fee9f9a48a58174,
0xbc8619321e55e68a, 0x3fee9feb564267c9,
0x3c909ccb5e09d4d3, 0x3feea0694fde5d3f,
0xbc7b32dcb94da51d, 0x3feea11473eb0187,
0x3c94ecfd5467c06b, 0x3feea1ed0130c132,
0x3c65ebe1abd66c55, 0x3feea2f336cf4e62,
0xbc88a1c52fb3cf42, 0x3feea427543e1a12,
0xbc9369b6f13b3734, 0x3feea589994cce13,
0xbc805e843a19ff1e, 0x3feea71a4623c7ad,
0xbc94d450d872576e, 0x3feea8d99b4492ed,
0x3c90ad675b0e8a00, 0x3feeaac7d98a6699,
0x3c8db72fc1f0eab4, 0x3feeace5422aa0db,
0xbc65b6609cc5e7ff, 0x3feeaf3216b5448c,
0x3c7bf68359f35f44, 0x3feeb1ae99157736,
0xbc93091fa71e3d83, 0x3feeb45b0b91ffc6,
0xbc5da9b88b6c1e29, 0x3feeb737b0cdc5e5,
0xbc6c23f97c90b959, 0x3feeba44cbc8520f,
0xbc92434322f4f9aa, 0x3feebd829fde4e50,
0xbc85ca6cd7668e4b, 0x3feec0f170ca07ba,
0x3c71affc2b91ce27, 0x3feec49182a3f090,
0x3c6dd235e10a73bb, 0x3feec86319e32323,
0xbc87c50422622263, 0x3feecc667b5de565,
0x3c8b1c86e3e231d5, 0x3feed09bec4a2d33,
0xbc91bbd1d3bcbb15, 0x3feed503b23e255d,
0x3c90cc319cee31d2, 0x3feed99e1330b358,
0x3c8469846e735ab3, 0x3feede6b5579fdbf,
0xbc82dfcd978e9db4, 0x3feee36bbfd3f37a,
0x3c8c1a7792cb3387, 0x3feee89f995ad3ad,
0xbc907b8f4ad1d9fa, 0x3feeee07298db666,
0xbc55c3d956dcaeba, 0x3feef3a2b84f15fb,
0xbc90a40e3da6f640, 0x3feef9728de5593a,
0xbc68d6f438ad9334, 0x3feeff76f2fb5e47,
0xbc91eee26b588a35, 0x3fef05b030a1064a,
0x3c74ffd70a5fddcd, 0x3fef0c1e904bc1d2,
0xbc91bdfbfa9298ac, 0x3fef12c25bd71e09,
0x3c736eae30af0cb3, 0x3fef199bdd85529c,
0x3c8ee3325c9ffd94, 0x3fef20ab5fffd07a,
0x3c84e08fd10959ac, 0x3fef27f12e57d14b,
0x3c63cdaf384e1a67, 0x3fef2f6d9406e7b5,
0x3c676b2c6c921968, 0x3fef3720dcef9069,
0xbc808a1883ccb5d2, 0x3fef3f0b555dc3fa,
0xbc8fad5d3ffffa6f, 0x3fef472d4a07897c,
0xbc900dae3875a949, 0x3fef4f87080d89f2,
0x3c74a385a63d07a7, 0x3fef5818dcfba487,
0xbc82919e2040220f, 0x3fef60e316c98398,
0x3c8e5a50d5c192ac, 0x3fef69e603db3285,
0x3c843a59ac016b4b, 0x3fef7321f301b460,
0xbc82d52107b43e1f, 0x3fef7c97337b9b5f,
0xbc892ab93b470dc9, 0x3fef864614f5a129,
0x3c74b604603a88d3, 0x3fef902ee78b3ff6,
0x3c83c5ec519d7271, 0x3fef9a51fbc74c83,
0xbc8ff7128fd391f0, 0x3fefa4afa2a490da,
0xbc8dae98e223747d, 0x3fefaf482d8e67f1,
0x3c8ec3bc41aa2008, 0x3fefba1bee615a27,
0x3c842b94c3a9eb32, 0x3fefc52b376bba97,
0x3c8a64a931d185ee, 0x3fefd0765b6e4540,
0xbc8e37bae43be3ed, 0x3fefdbfdad9cbe14,
0x3c77893b4d91cd9d, 0x3fefe7c1819e90d8,
0x3c5305c14160cc89, 0x3feff3c22b8f71f1,
},
};

inline FORCE_INLINE double my_exp(double x)
{
  uint64_t ki, idx, top, sbits;
  double_t kd, z, r, r2, scale, tail, tmp;
  /* exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].  */
  /* x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].  */
  z = InvLn2N * x;
  kd = std::round(z);
  ki = int(kd);
  r = fmadd(kd, NegLn2hiN, x);
  r = fmadd(kd, NegLn2loN, r);
  /* 2^(k/N) ~= scale * (1 + tail).  */
  tail  = asdouble(T[2*(ki % N) + 0]);
  scale = asdouble(T[2*(ki % N) + 1] + (ki << (52 - EXP_TABLE_BITS)));
  /* exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).  */
  r2 = r*r;
  tmp = fmadd(r2*r2, fmadd(C5,r,C4), fmadd(r2, fmadd(C3,r,C2), r + tail));
  return fmadd(scale, tmp, scale);
}


struct entry {
   double off;
   double expfac;
   const packed<double,2>* tab;
   entry(double off_, double expfac_, const packed<double,2>* tab_) : off(off_), expfac(expfac_), tab(tab_) {}
};
static entry LUT[41] = {
   entry(0.0, -0.5, Tem38_to_m10), // -20
   entry(0.0, -0.5, Tem38_to_m10), // -19
   entry(0.0, -0.5, Tem38_to_m10), // -18
   entry(0.0, -0.5, Tem38_to_m10), // -17
   entry(0.0, -0.5, Tem38_to_m10), // -16
   entry(0.0, -0.5, Tem38_to_m10), // -15
   entry(0.0, -0.5, Tem38_to_m10), // -14
   entry(0.0, -0.5, Tem38_to_m10), // -13
   entry(0.0, -0.5, Tem38_to_m10), // -12
   entry(0.0, -0.5, Tem38_to_m10), // -11
   entry(0.0, -0.5, Tem38_to_m10), // -10
   entry(0.0, -0.5, Tem10_to_m3 ), // - 9
   entry(0.0, -0.5, Tem10_to_m3 ), // - 8
   entry(0.0, -0.5, Tem10_to_m3 ), // - 7
   entry(0.0, -0.5, Tem10_to_m3 ), // - 6
   entry(0.0, -0.5, Tem10_to_m3 ), // - 5
   entry(0.0, -0.5, Tem10_to_m3 ), // - 4
   entry(0.0, -0.5, Tem10_to_m3 ), // - 3
   entry(0.0, -0.5, Tem3_to_m1  ), // - 2
   entry(0.0, -0.5, Tem3_to_m1  ), // - 1
   entry(0.0,  0.0, nullptr     ), // - 0
   entry(0.0,  0.0, T1_to_2     ), //   1
   entry(1.0, -0.5, Te2_to_9    ), //   2
   entry(1.0, -0.5, Te2_to_9    ), //   3
   entry(1.0, -0.5, Te2_to_9    ), //   4
   entry(1.0, -0.5, Te2_to_9    ), //   5
   entry(1.0, -0.5, Te2_to_9    ), //   6
   entry(1.0, -0.5, Te2_to_9    ), //   7
   entry(1.0, -0.5, Te2_to_9    ), //   8
   entry(1.0,  0.0, Tzero       ), //   9
   entry(1.0,  0.0, Tzero       ), //  10
   entry(1.0,  0.0, Tzero       ), //  11
   entry(1.0,  0.0, Tzero       ), //  12
   entry(1.0,  0.0, Tzero       ), //  13
   entry(1.0,  0.0, Tzero       ), //  14
   entry(1.0,  0.0, Tzero       ), //  15
   entry(1.0,  0.0, Tzero       ), //  16
   entry(1.0,  0.0, Tzero       ), //  17
   entry(1.0,  0.0, Tzero       ), //  18
   entry(1.0,  0.0, Tzero       ), //  19
   entry(1.0,  0.0, Tzero       ), //  20
};

/*
   Rn(x) stands for some rational function of x of degree n-1 (2*n total cofficients).
   The approximations for double precision in the various ranges look like
      [-inf,-38]: gaussint(x) = 0
      [-38,-10]:  gaussint(x) = exp(-1/2*x^2)*R8(x)
      [-10,-3]:   gaussint(x) = exp(-1/2*x^2)*R8(x)
      [-3,-1]:    gaussint(x) = exp(-1/2*x^2)*R8(x)
      [-1,1]:     gaussint(x) = 1/2 + x*R5(x^2)
      [1,2]:      gaussint(x) = R8(x)
      [2,9]:      gaussint(x) = 1 - exp(-1/2*x^2)*R8(x)
      [9,inf]:    gaussint(x) = 1

   It seems necessary to include the factors exp(-1/2*x^2) when they are
   present because they smooth out function being approximated and greaty increase
   the effectiveness of the minimax approximations.

   In the range [-20,-1] U [1,20] all approximations are combined into one
   lookup table, and the exp is calculated via a fast inline version
   because it is nowhere near underflow.

   In the range [-38,-37] the exp starts to underflow. This is a non-issue
   because the corresponding R8(x) is less than 1 as 
      R8(x) ~= 1/(sqrt(2*pi)*x)

   Accuracy: Sampling all rational inputs in [-20,20] with a denominator dividing 64
   yields a max relative error of 5.16e-16 at x = -17.375.
         the fmadd-evaluated formula gives 6.380314775457178e-68 = 0x1.b85a67826b46ap-224
         and correctly rounded value is    6.380314775457182e-68 = 0x1.b85a67826b46ep-224
*/
FORCE_NOINLINE double gaussint(double x)
{
   double x2 = x*x;
   sword n = static_cast<sword>(x); // rounding toward zero needed
   if (LIKELY(x2 <= 400)) {
      ASSERT(-20 <= n && n <= 20);
      double earg = static_cast<double>(LUT[n+20].expfac)*x2;
      double off = static_cast<double>(LUT[n+20].off);
      auto lut = LUT[n+20].tab;
      if (LIKELY(n != 0)) {
         double e = my_exp(earg);
         auto f = horner2<packed<double,2>,8>(lut, x, x2);
         return fmadd(e, f[0]/f[1], off);
      } else {
         auto f = horner<packed<double,2>,5>(T1, x2);
         return fmadd(x, f[0]/f[1], 0.5);
      }
   } else {
      if (x > 0) {
         return 1.0;
      } else if (x < 0) {
         if (x > -38) {
            auto f = horner2<packed<double,2>,8>(Tem38_to_m10, x, x2);
            return std::exp(-0.5*x2)*(f[0]/f[1]);  // std::exp because close to underflow
         } else {
            return 0.0;
         }
      } else {
         return x;
      }
   }
}

double gaussint_ref(double x) {
   double xs = x*0.707106781186547524401;
   if (x > 1)
      return 1.0 - 0.5*std::erfc(xs);
   else if (x < -1)
      return 0.5*std::erfc(-xs);
   else
      return 0.5 + 0.5*std::erf(xs);
}

double relerr(double a, double b) {
   return std::abs(a-b)/std::max(std::abs(a), std::abs(b));
}

double testf(random_state& state, double (*f)(double), uword nreps = 10000000)
{
   tmp_allocator push;
   double* r = push.recursive_alloc<double>(nreps);
   for (uword i = 0; i < nreps; i++) {
      double a = -20.0;
      double b = +20.0;
      double x = state.get_uniform_unit<double>();
      r[i] = 0.5*((a+b) + (b-a)*x);
   }

   uword t1 = get_timestamp();
   for (uword i = nreps; i > 0; i--) {
      r[i-1] = f(r[i-1]);
   }
   uword t2 = get_timestamp();

   return sword(t2 - t1)/sword(nreps);
}

double testtab[] = {2.7536241186062336951e-89,3.7662263334836240645e-89,5.1499428190388342654e-89,7.0403237620598904797e-89,9.6222600364827875694e-89,1.3147881495778655039e-88,1.7960925318584341315e-88,2.4529905085539114303e-88,3.3493247875065261e-88,4.5720702514275450532e-88,6.2396863856681483134e-88,8.5134750101441807859e-88,1.1613020851210218879e-87,1.5837177462446084115e-87,2.1592583024383053332e-87,2.9432398234314227876e-87,4.0108917631137030177e-87,5.4645005309425636981e-87,7.4431066486513431594e-87,1.0135665065943082219e-86,1.3798901032686610472e-86,1.878153091078421741e-86,2.5557108680879568397e-86,3.4768559002856967379e-86,4.728853992980909886e-86,6.4301241877427025857e-86,8.741320803665996653e-86,1.1880343570679967338e-85,1.6142663174063522179e-85,2.192883772200104488e-85,2.9781755005028509583e-85,4.0437025727695737022e-85,5.4891154756604099475e-85,7.4493739465076358037e-85,1.0107213281096316672e-84,1.3709996099143926539e-84,1.859248668465094732e-84,2.5207621804034156109e-84,3.4168075450595728869e-84,4.6302390004764620946e-84,6.2730759925285082706e-84,8.4967323431016443588e-84,1.1505820019971231742e-83,1.5576772063224661818e-83,2.1082960927163342773e-83,2.8528569387523383076e-83,3.8594252602525989859e-83,5.2198680878544283819e-83,7.0581465785834792445e-83,9.5414871804700076344e-83,1.2895426847310043856e-82,1.7424070785580980706e-82,2.353736192453177599e-82,3.1787779986463436945e-82,4.2919718660641557761e-82,5.7935904276998913839e-82,7.8186715447061467094e-82,1.0549027021341289573e-81,1.4229383919473699703e-81,1.9189075283028613531e-81,2.5871180974001240644e-81,3.4871667845986716854e-81,4.6991947434872856334e-81,6.3309433362057489093e-81,8.5272239526309765105e-81,1.1482624733502767563e-80,1.5458556403337501493e-80,2.0806113293674640204e-80,2.799672584154638522e-80,3.7663249901290332298e-80,5.0655032883650991722e-80,6.8111692822513274356e-80,9.1561945111588756685e-80,1.2305594780306917152e-79,1.6534249003163862635e-79,2.2210615938249653682e-79,2.9828471292198971043e-79,4.0049366976535319417e-79,5.3759420895849973069e-79,7.2145255246497344953e-79,9.6795514791342034313e-79,1.2983654826817848462e-78,1.7411370479981956773e-78,2.3343352088481755954e-78,3.1288716891133317812e-78,4.1928232637463995883e-78,5.6171966627501621073e-78,7.5236220833031172308e-78,1.0074617976746963338e-77,1.3487283082693851091e-77,1.8051555502370300485e-77,2.4154556153908936409e-77,3.231303606276285521e-77,4.3216614538555932801e-77,5.7785382696354350936e-77,7.7246626074064151977e-77,1.0323698689563289609e-76,1.3793846764306790766e-76,1.8425944116840324625e-76,2.4607550148053844829e-76,3.2854985895182178319e-76,4.3855944534367624695e-76,5.8526146610607594231e-76,7.808465062917017308e-76,1.0415393709124496536e-75,1.3889288234811220695e-75,1.8517338694597799414e-75,2.4681493199726017575e-75,3.2889598502530222864e-75,4.3816733191034873171e-75,5.8360056077400683989e-75,7.7711571017035671155e-75,1.0345463677570100615e-74,1.3769193894190166181e-74,1.8321513852648959691e-74,2.4372971454902203586e-74,3.2415288395428497446e-74,4.3100825258047140427e-74,5.729484754128485245e-74,7.6144727625890760488e-74,1.0117154785399613579e-73,1.3439132613483795478e-73,1.7847540373774250566e-73,2.3696259425623591124e-73,3.145397290215618652e-73,4.1741255816872454868e-73,5.5379598713163184776e-73,7.3456188394970550305e-73,9.7409489189371504826e-73,1.291422756664683255e-72,1.7117088537705789988e-72,2.2682224917615319142e-72,3.004939290146891603e-72,3.9799722724755555649e-72,5.2700979352278974708e-72,6.9767252597517049242e-72,9.2337657417417730654e-72,1.2218007311206189815e-71,1.6162785387980673571e-71,2.1375994544943497959e-71,2.8263813070037469773e-71,3.736194470541266817e-71,4.9376749225174976525e-71,6.5239378542833451829e-71,8.6177013091948126816e-71,1.1380657686048830966e-70,1.5025799373764473732e-70,1.9833625724502389758e-70,2.6173448327816989447e-70,3.4531391455728440528e-70,4.5547182532389365258e-70,6.0062480077797457136e-70,7.9184341572394480923e-70,1.0436855358306085931e-69,1.3752901298489972835e-69,1.8118126799806771054e-69,2.3863084058000589129e-69,3.1422023977094414151e-69,4.1365287697925449797e-69,5.4441771980057627572e-69,7.1634587662350358454e-69,9.4233991962637169401e-69,1.2393293198606707427e-68,1.6295217863584046514e-68,2.1420417471123528154e-68,2.8150752159527030169e-68,3.6986768534089442806e-68,4.8584427547151758952e-68,6.3803147754571814798e-68,8.3768634519114269418e-68,1.0995502540265110104e-67,1.4429226373712320312e-67,1.8930641363039957634e-67,2.4830297723183077167e-67,3.2560633387522878531e-67,4.2687240521137581064e-67,5.5949683949048849452e-67,7.3314776421613048213e-67,9.604609504792336631e-67,1.2579466763193610279e-66,1.647172465243011129e-66,2.1563052416041464671e-66,2.8221216267103162278e-66,3.6926278371158788959e-66,4.8304731958223797703e-66,6.3173965548022032596e-66,8.2600167490633690141e-66,1.0797371479651941057e-65,1.4110730134605852064e-65,1.8436363269527449444e-65,2.4082155493076007759e-65,3.1449213209263639131e-65,4.1059962020989062896e-65,5.3594673852482205521e-65,6.9938941099953107671e-65,9.1245363314209739988e-65,1.1901368143523794411e-64,1.5519485059087454935e-64,2.0232616918120308143e-64,2.6370668829016189569e-64,3.4362485211887880364e-64,4.4765384893209556527e-64,5.8303471778566935697e-64,7.5917316230908929871e-64,9.8828366054712001244e-64,1.2862243346110297862e-63,1.6735788076544631414e-63,2.1770577695685493363e-63,2.8313142815440512382e-63,3.6812943719306315675e-63,4.7852800145634569955e-63,6.218827459393573391e-63,8.079863097391510122e-63,1.049527544317681562e-62,1.3629440375044204925e-62,1.7695243875530393636e-62,2.2968330394155251073e-62,2.9805516194683920541e-62,3.8668579272014643879e-62,5.0154988145274544418e-62,6.5037583576689000735e-62,8.4315808982561427988e-62,1.0928183420103342311e-61,1.4160588116635117465e-61,1.83446300316473111e-61,2.3759154208470926323e-61,3.0764321983366932278e-61,3.9825209902056437627e-61,5.1542223897425251755e-61,6.6690287231849460717e-61,8.626932071316279455e-61,1.1156925301530464217e-60,1.4425372266545481286e-60,1.8646781918796991731e-60,2.4097673306364897691e-60,3.1134413741279943564e-60,4.0216163009681090437e-60,5.1934377875065093828e-60,6.7050743213400754428e-60,8.6545924380084591168e-60,1.1168221242476203505e-59,1.4408400121155699862e-59,1.8584114665160547069e-59,2.3964168122449427704e-59,3.0894217880359330152e-59,3.9818639672761250592e-59,5.130858061454176605e-59,6.6097943932983638152e-59,8.5129536553876437445e-59,1.0961423507485076713e-58,1.4110682384644269184e-58,1.8160317901267828633e-58,2.3366477582341720686e-58,3.0057816159331397684e-58,3.8655916312747566726e-58,4.9701431242702892909e-58,6.3887544005380872813e-58,8.2102780806053434736e-58,1.0548577967130695336e-57,1.3549533559717528556e-57,1.739999609785131164e-57,2.233923700678503194e-57,2.8673578734153854659e-57,3.6795087966566584623e-57,4.7205453256663497711e-57,6.0546477347908803731e-57,7.7639006972704078615e-57,9.9532622510032161504e-57,1.2756904569893195128e-56,1.6346303253977229004e-56,2.0940554745840210603e-56,2.6819532326552218908e-56,3.4340657492721407802e-56,4.3960274947072842898e-56,5.6260888098567395128e-56,7.1985853271117868919e-56,9.2083564293161253678e-56,1.1776370953614259242e-55,1.5056888215062040493e-55,1.9246571094128765749e-55,2.459608043368449817e-55,3.1424822694363566526e-55,4.0139703784687316641e-55,5.1258971416366740563e-55,6.5442519142864861747e-55,8.3530393745278443663e-55,1.065917144759156372e-54,1.3598680395108560938e-54,1.7344607917938700513e-54,2.2117019898896871623e-54,2.8195716361204307027e-54,3.5936357986823390307e-54,4.579092251821903729e-54,5.8333647150319089919e-54,7.4293919973528075248e-54,9.4597971517079707767e-54,1.2042170773396668471e-53,1.5325764518319550496e-53,1.949996915106238104e-53,2.4805050221470059418e-53,3.1545739179365502534e-53,4.0108435133125833081e-53,5.0982971055175700685e-53,6.4790148758154011937e-53,8.2316562905314155529e-53,1.0455863266029367139e-52,1.3277826168978334365e-52,1.6857317997929893727e-52,2.1396581809565285708e-52,2.7151556860494359484e-52,3.4446055371335625187e-52,4.3689660991088017129e-52,5.5400320371518537533e-52,7.0232851333632415461e-52,8.901490820459535564e-52,1.1279234369658349643e-51,1.4288640811087078974e-51,1.8096585692822496995e-51,2.2913782990990226701e-51,2.9006235995972850175e-51,3.6709661993127508858e-51,4.6447658515371300926e-51,5.8754575994667432879e-51,7.4304308957463053357e-51,9.3946528240874358979e-51,1.1875226600878751588e-50,1.5007125359259316513e-50,1.8960402436237748751e-50,2.3949256121695641252e-50,3.0243422996937320192e-50,3.8182494479347884281e-50,4.8193902110879551673e-50,6.0815506120148439663e-50,7.6723958394785289946e-50,9.6770307143657778906e-50,1.2202468109482323264e-49,1.5383235465068449506e-49,1.9388407520500541626e-49,2.4430425876456976867e-49,3.0776156622347635032e-49,3.8760750514910227565e-49,4.8805011691323118196e-49,6.1437157644323225724e-49,7.7320073901737038133e-49,9.7285442408027829602e-49,1.2237646653996830192e-48,1.5390134486226883759e-48,1.9350018113187996457e-48,2.432286858167744064e-48,3.0566285701436608489e-48,3.840298665766472699e-48,4.8237167076883784447e-48,6.0574947644152207796e-48,7.6049920121814658704e-48,9.5455056619479333925e-48,1.1978255721404550505e-47,1.5027359842544487485e-47,1.8848042713091332278e-47,2.3634384423082767748e-47,2.9628986828741323078e-47,3.7135029683861072057e-47,4.6531303605059777824e-47,5.829095022677565384e-47,7.3004817936437999887e-47,9.1410562728837819868e-47,1.1442889835262069943e-46,1.4320874089022449449e-46,1.7918341613061332141e-46,2.2414062326936399914e-46,2.8030950000045726725e-46,3.5046894214412297523e-46,4.3808233294748792368e-46,5.4746507931743543662e-46,6.8399289165585165936e-46,8.5436065133556273194e-46,1.0669040725329624461e-45,1.3319992910664114326e-45,1.6625591354382103733e-45,2.0746493191582424121e-45,2.5882533420310207378e-45,3.2282217524992798201e-45,4.0254499132372236053e-45,5.0183389103225114733e-45,6.2546072241275373887e-45,7.7935368191928002544e-45,9.7087571278782708884e-45,1.209169488242082468e-44,1.5055847981346326096e-44,1.8742078900824964207e-44,2.3325169232167703041e-44,2.9021933778183095518e-44,3.610126277857209877e-44,4.4896547340839899411e-44,5.5821049821470999302e-44,6.9386912398101924346e-44,8.6228659134801592807e-44,1.0713224652335862241e-43,1.3307096345383816796e-43,1.6524978448143288541e-43,2.0516015320176239815e-43,2.5464763159739570043e-43,3.1599541611100079572e-43,3.9202741616230448326e-43,4.8623544558882327395e-43,6.0293612853201779787e-43,7.4746441335424782674e-43,9.2641217612460374833e-43,1.1479223463174752412e-42,1.4220513840677808672e-42,1.7612158816378449969e-42,2.1807426754187585839e-42,2.6995462903636170892e-42,3.340962981784774743e-42,4.1337773170701541062e-42,5.1134854300305531639e-42,6.3238491233269362261e-42,7.8188073056578912157e-42,9.6648263386904615036e-42,1.1943789355818455446e-41,1.4756547260582779002e-41,1.8227281845412811254e-41,2.2508865425348874688e-41,2.7789442941571430439e-41,3.4300513346830572774e-41,4.2326849303180303166e-41,5.2218670320349596464e-41,6.4406577513183898247e-41,7.9419871906293981064e-41,9.7909017232703136579e-41,1.2067317803099211333e-40,1.4869397133048853191e-40,1.8317682359717825927e-40,2.2560163396857890795e-40,2.7778482240190467241e-40,3.4195530213214742096e-40,4.2084747796381036584e-40,5.1781505752190133123e-40,6.3697029857326959068e-40,7.8335433462452329331e-40,9.6314546299817554858e-40,1.1839137928519191058e-39,1.4549324939151997426e-39,1.787558131483605239e-39,2.195695306171616458e-39,2.6963641434018883413e-39,3.3103932254657306305e-39,4.0632654834059098963e-39,4.9861505552387224361e-39,6.1171643995498796823e-39,7.5029058062573090278e-39,9.2003302051581188208e-39,1.1279034244818495995e-38,1.3824040492682328268e-38,1.6939190889533392842e-38,2.0751281001392188271e-38,2.5415095525901571864e-38,3.1119539989087361135e-38,3.809510539590283149e-38,4.6622953327485512296e-38,5.7045970488425370787e-38,6.9782216270608677677e-38,8.5341277263163391311e-38,1.0434415212059208527e-37,1.2754742281487163866e-37,1.5587262888811992081e-37,1.904419557428939938e-37,2.3262158332354799479e-37,2.8407432628383707712e-37,3.4682354119616476761e-37,4.233306929545630712e-37,5.1658947623446517986e-37,6.3023999675810378388e-37,7.6870725294701756137e-37,9.3736904756221862811e-37,1.1427595327535455293e-36,1.3928158886455247421e-36,1.697177201039138365e-36,2.067546493221628988e-36,2.5181291467073023907e-36,3.0661636959502549761e-36,3.7325642988777133772e-36,4.5426981841917659454e-36,5.5273261968939978455e-36,6.7237403726279083495e-36,8.1771394701711018668e-36,9.9422918211867000554e-36,1.2085545007163764688e-35,1.4687254093576554518e-35,1.784471485872796787e-35,2.1675706151255652074e-35,2.6322766797893599168e-35,3.1958358083957139966e-35,3.8791093609198041529e-35,4.7073255319481531181e-35,5.7109858974716940648e-35,6.926958572697236876e-35,8.3997960636334176589e-35,1.0183323598561655713e-34,1.2342552973057985294e-34,1.495598804068333158e-34,1.8118401297841544678e-34,2.1944176984744095677e-34,2.6571335280163300392e-34,3.216637513378980926e-34,3.8930100807404207899e-34,4.710463018214514162e-34,5.6981822406196119463e-34,6.8913409786119039156e-34,8.332317549552754811e-34,1.0072158651315057422e-33,1.2172337238851219722e-33,1.4706863756874546068e-33,1.7764821120776789977e-33,2.1453407732039311052e-33,2.5901589425211081997e-33,3.1264481097094107985e-33,3.7728602492030690313e-33,4.5518180957070171246e-33,5.4902707826994188318e-33,6.6205995511135628874e-33,7.9817030592877487432e-33,9.6202975819752935354e-33,1.1592474254326828726e-32,1.396556370875586839e-32,1.6820368221143052176e-32,2.0253833128324921265e-32,2.4382243158307504348e-32,2.9345045852433198939e-32,3.5309423958859932586e-32,4.2475762137827208016e-32,5.1084181252086901721e-32,6.1422346764321576225e-32,7.3834797328669167012e-32,8.8734086733593168972e-32,1.0661408833734466489e-31,1.280658777060154184e-31,1.5379668829845260555e-31,1.8465252908964712798e-31,2.2164516476613195357e-31,2.6598429185826117352e-31,3.1911590179682452642e-31,3.8276800901475214727e-31,4.5900514432631508379e-31,5.5029327737605438034e-31,6.5957714461136750791e-31,7.9037233004804079054e-31,9.4687488576950285274e-31,1.1340918002789864238e-30,1.3579962404589949943e-30,1.6257122246465673567e-30,1.9457342510666637356e-30,2.3281884322119713001e-30,2.7851429007937065463e-30,3.3309766909317966735e-30,3.9828179997936382729e-30,4.7610647476877055547e-30,5.6900027349909626672e-30,6.798539508633455664e-30,8.1210753774347496899e-30,9.6985369465928574965e-30,1.1579603185686417654e-29,1.3822159529583282581e-29,1.6495021988259143119e-29,1.9679980886988583798e-29,2.3474222881278247536e-29,2.7993200536669754795e-29,3.3374031320768851444e-29,3.977952266291616155e-29,4.7402937194715979418e-29,5.6473632860973805675e-29,6.7263736836225372558e-29,8.0096040743869728517e-29,9.5353338312560303841e-29,1.1348946620968979085e-28,1.3504235541107094668e-28,1.6064945532717536391e-28,1.9106595744986757112e-28,2.2718632119849124989e-28,2.7006969391058129772e-28,3.2096992186122819941e-28,3.8137097252562240497e-28,4.5302873320714363767e-28,5.3802032136921402876e-28,6.3880224190088506768e-28,7.5827896122632516205e-28,8.9988374361959265705e-28,1.0676739183071078709e-27,1.2664431251165092987e-27,1.5018535311328052281e-27,1.7805915322383528051e-27,2.110551064615989048e-27,2.5010493675586882988e-27,2.9630808780943585853e-27,3.5096159208320410137e-27,4.1559520073834124367e-27,4.9201269070909047935e-27,5.8234042277126743358e-27,6.8908440897973778162e-27,8.1519736367661152034e-27,9.6415746467100419923e-27,1.1400608462803911688e-26,1.3477301908329593007e-26,1.5928421882655551674e-26,1.8820771042851147105e-26,2.2232942474536280465e-26,2.6257377676144159231e-26,3.1002779675518013462e-26,3.6596941844086012077e-26,4.3190063178092303465e-26,5.0958632718399506313e-26,6.0109979659060780474e-26,7.0887601874213907064e-26,8.3577404449306504654e-26,9.8515001773413903485e-26,1.1609426234208618951e-25,1.3677730522400718783e-25,1.6110619184357563166e-25,1.8971659711865949423e-25,2.2335379098837562412e-25,2.6289131603675404771e-25,3.0935281050138933496e-25,3.6393749988503634647e-25,4.2804996632328119231e-25,5.0333490472648138072e-25,5.9171769073656178503e-25,6.9545172029758819014e-25,8.1717363711766063381e-25,9.599677459515507705e-25,1.1274411204527677237e-24,1.323811158949185796e-24,1.5540076252277647374e-24,1.8237915404428059025e-24,2.1398936737196289872e-24,2.5101758211489043672e-24,2.9438185751689486693e-24,3.451539879796828085e-24,4.0458493544704729476e-24,4.7413441650303496302e-24,5.5550531414741490274e-24,6.5068369080211525201e-24,7.619853024160526066e-24,8.9210965615904447072e-24,1.0442028191082502249e-23,1.2219303759657977736e-23,1.4295621541576481812e-23,1.6720705892035703914e-23,1.9552448972562696176e-23,2.2858235612140870236e-23,2.6716480287532248167e-23,3.1218409729804991081e-23,3.6470129883497467284e-23,4.2595021965591553367e-23,4.9736519314709646537e-23,5.8061324725911419525e-23,6.7763137193047887718e-23,7.9066967612434647506e-23,9.2234135249394181485e-23,1.0756805087485894653e-22,1.2542090872021723362e-22,1.4620142809320145732e-22,1.7038380701082068747e-22,1.9851807495549177461e-22,2.3124206032708511375e-22,2.6929522089589358921e-22,3.1353462318874763957e-22,3.6495339998311041372e-22,4.2470206476676916542e-22,4.9411311908760526415e-22,5.7472945424883824557e-22,6.6833712403538012027e-22,7.770031514990511413e-22,9.0311913189920099492e-22,1.0494515075362607493e-21,1.2191995205375166775e-21,1.4160619990654651175e-21,1.6443143036698162041e-21,1.9088969567343640576e-21,2.2155177027629505882e-21,2.5707690046932005577e-21,2.9822632761862595798e-21,3.4587884872335800397e-21,4.0104871665106872421e-21,4.6490622656026408198e-21,5.3880138553995535576e-21,6.2429112025809131746e-21,7.2317054343867722681e-21,8.3750887544380714959e-21,9.6969070344491538404e-21,1.122463359132798266e-20,1.2989913083508600655e-20,1.5029185743879838063e-20,1.7384403631503056621e-20,2.0103852255611309169e-20,2.3243092831741545195e-20,2.6866042603744333502e-20,3.1046213143718663348e-20,3.5868129366497283884e-20,4.1428955213674537982e-20,4.7840355628233941143e-20,5.5230628615747385543e-20,6.3747145941225602322e-20,7.3559146420613100296e-20,8.4860931921735715329e-20,9.7875523192172870406e-20,1.1285884059538406477e-19,1.3010448389095242599e-19,1.4994919548660062846e-19,1.7277910328458185617e-19,1.9903685253085042662e-19,2.2922975116438247154e-19,2.6393907029622427316e-19,3.0383066089469541809e-19,3.4966706982139007807e-19,4.0232136339874911549e-19,4.6279289508544661233e-19,5.3222528603264391488e-19,6.1192692379126954061e-19,7.0339432580098193602e-19,8.0833876115066666255e-19,9.2871657717477232598e-19,1.0667637375474858003e-18,1.2250351465685001506e-18,1.4064494113263724927e-18,1.6143397805346065643e-18,1.8525120973640817054e-18,2.1253107150098219453e-18,2.4376934496811696007e-18,2.7953167880499122959e-18,3.2046327270207596575e-18,3.6729988053533166422e-18,4.2088030918072204152e-18,4.8216061260730898064e-18,5.5223020701296050893e-18,6.3233016225815386962e-18,7.2387395811965817635e-18,8.2847103139990239201e-18,9.4795348222033183542e-18,1.0844063554935935042e-17,1.2402019672780564723e-17,1.4180388062181988768e-17,1.6209856084057607457e-17,1.8525312807007280755e-17,2.1166414338830591989e-17,2.4178223841472180968e-17,2.7611935907252616192e-17,3.1525696203117875469e-17,3.5985528671233795164e-17,4.106638412711216508e-17,4.6853325841330381009e-17,5.3442869650996983651e-17,6.0944498348355937506e-17,6.9482372565296592333e-17,7.919726314642477341e-17,9.0248733115919252816e-17,1.0281760083502848576e-16,1.1710871986283027359e-16,1.3335411542308186756e-16,1.5181652230074356778e-16,1.7279337450564592119e-16,1.9662130321756201278e-16,2.2368120644441050467e-16,2.5440396157010670808e-16,2.8927686063722590044e-16,3.2885085790963938824e-16,3.7374873011130015499e-16,4.2467426187257848014e-16,4.8242258248316495474e-16,5.4789179521539443687e-16,6.2209605742717841235e-16,7.0618028858330207856e-16,8.0143670447405791972e-16,9.0932339951258424896e-16,1.0314852253362169947e-15,1.1697772433328206463e-15,1.3262910615047228373e-15,1.5033844026525641441e-15,1.7037142916328732075e-15,1.9302742948864479817e-15,2.1864362960706803378e-15,2.4759973480352417779e-15,2.8032322041968334875e-15,3.1729522023036343782e-15,3.5905712514183292111e-15,4.0621797595587663469e-15,4.5946274357785954602e-15,5.195616007610049193e-15,5.8738030139034052695e-15,6.6389179654828509761e-15,7.5018923131336120483e-15,8.4750048258282593259e-15,9.572044163545130403e-15,1.0808490630465748394e-14,1.2201719317899234759e-14,1.3771227094329743358e-14,1.553888617512127473e-14,1.7529227309514161537e-14,1.9769755960774854e-14,2.229130523020482087e-14,2.5128429691020393762e-14,2.8319844758118643785e-14,3.1908916729108962278e-14,3.5944209195850430887e-14,4.0480092149741948049e-14,4.5577420794382140367e-14,5.1304291842786688617e-14,5.7736885920639461147e-14,6.4960405630323749944e-14,7.3070119861810333441e-14,8.2172526075843372513e-14,9.2386643543215120184e-14,1.0384545191327865576e-13,1.1669749101840814784e-13,1.3110863951335308415e-13,1.472640918152211233e-13,1.6537055486856833873e-13,1.8565868852985160949e-13,2.0838581586720694312e-13,2.3383893242805330403e-13,2.6233804656345224993e-13,2.9423988624191975438e-13,3.2994201146651548761e-13,3.6988737546037443184e-13,4.1456938224331767484e-13,4.6453749312505392883e-13,5.204034400316781493e-13,5.8284810950857648986e-13,6.5262916775566686076e-13,7.3058950420605198153e-13,8.1766657901766044647e-13,9.149027684758831881e-13,1.023456811776201475e-12,1.1446164730485354628e-12,1.2798125438858350044e-12,1.4306343241423417043e-12,1.598846732474427056e-12,1.7864092131205221678e-12,1.9954966218778487863e-12,2.2285222922646015873e-12,2.4881635026006916672e-12,2.7773895863544933101e-12,3.0994929517572154306e-12,3.4581233025652678998e-12,3.8573253801549928226e-12,4.3015805780813278145e-12,4.7958528140589781865e-12,5.3456390812882952603e-12,5.9570251414268948362e-12,6.6367468656042665709e-12,7.3922577780178224195e-12,8.2318034091900680151e-12,9.1645031232925383409e-12,1.020044014646945587e-11,1.1350760591273835392e-11,1.2627782346649196837e-11,1.4045114783879381989e-11,1.5617790317158512797e-11,1.7362408953520568751e-11,1.9297297071481001189e-11,2.1442681781602537137e-11,2.3820882346082881911e-11,2.6456520269214492411e-11,2.9376749817093723474e-11,3.2611510884237823163e-11,3.6193806297859729956e-11,4.0160005838591178083e-11,4.45501794606615341e-11,4.9408462416255272888e-11,5.4783455229408999403e-11,6.0728661725902389577e-11,6.7302968608796068078e-11,7.4571170376347963727e-11,8.2604543711906546919e-11,9.1481475836086101718e-11,1.0128815170228035391e-10,1.121193053397344175e-10,1.2407904110651388212e-10,1.3728173111051324235e-10,1.5185299559306220346e-10,1.6793077364985859541e-10,1.8566649229124820522e-10,2.0522634252189388816e-10,2.2679267185249535204e-10,2.5056550344757771685e-10,2.7676419296789392935e-10,3.0562923508842658684e-10,3.3742423266840233849e-10,3.7243804262357172427e-10,4.1098711370905734515e-10,4.5341803266952844889e-10,5.0011029655893621151e-10,5.5147933048160384486e-10,6.0797977156764575619e-10,6.7010904167652938546e-10,7.3841123313166358495e-10,8.1348133373533399943e-10,8.9596981940684144492e-10,9.865876450376981407e-10,1.0861116665772313767e-9,1.1953905299616736614e-9,1.315351065292267705e-9,1.4470052276663568783e-9,1.5914576292839718842e-9,1.7499137109060444707e-9,1.9236886044449862641e-9,2.1142167424408471014e-9,2.3230622744345209801e-9,2.551930354812487831e-9,2.8026793715854913265e-9,3.0773341907976775275e-9,3.3781004968656327915e-9,3.7073803151423572501e-9,4.0677888094147716721e-9,4.4621724539016118731e-9,4.8936286866497691392e-9,5.365527159061134966e-9,5.8815327046503206156e-9,6.445630159069510023e-9,7.0621511729752581307e-9,7.7358031694902638006e-9,8.471700608870017738e-9,9.2753987345608213957e-9,1.0152929987175230962e-8,1.1110843286058999028e-8,1.2156246392127915766e-8,1.329685158056388978e-8,1.4541024867829999523e-8,1.5897839054349576791e-8,1.7377130862152886383e-8,1.8989562465887719384e-8,2.07466877358812509e-8,2.2661023533496910386e-8,2.4746126421922018275e-8,2.7016675179823067326e-8,2.9488559531092654866e-8,3.2178975531265566719e-8,3.510652808018532516e-8,3.8291341061244284186e-8,4.1755175640091525471e-8,4.5521557290198718601e-8,4.9615912149194103598e-8,5.4065713348522279478e-8,5.8900637999870186672e-8,6.4152735565029494723e-8,6.9856608381558874667e-8,7.6049605164887142511e-8,8.2772028358485350873e-8,9.0067356257562549756e-8,9.7982480888540797399e-8,1.065679626864794951e-7,1.1587830307579147043e-7,1.2597223612617251669e-7,1.369130405258062746e-7,1.4876887318776628501e-7,1.6161312588328054078e-7,1.7552480637731963354e-7,1.9058894562799124252e-7,2.0689703270164974189e-7,2.2454747915064343004e-7,2.4364611470041139831e-7,2.6430671619740013034e-7,2.8665157187919391167e-7,3.1081208314354484873e-7,3.369294061138538775e-7,3.6515513542530385986e-7,3.9565203278849394636e-7,4.2859480302628878223e-7,4.6417092042489828109e-7,5.0258150839216821847e-7,5.4404227557491632113e-7,5.8878451175312290538e-7,6.3705614700211166579e-7,6.891228777947688982e-7,7.4526936390458350911e-7,8.0580050016708166796e-7,8.7104276736231520674e-7,9.4134566669467785549e-7,1.0170832425687031713e-6,1.0986556985908760113e-6,1.1864911119680967289e-6,1.2810472517235019479e-6,1.3828135064100919022e-6,1.49231292727226123e-6,1.6101043930850913963e-6,1.7367849031913455195e-6,1.8729920055567094861e-6,2.0194063669751372615e-6,2.1767544928783636475e-6,2.345811604536913664e-6,2.5274046817844210212e-6,2.722415679752911828e-6,2.9317849284740415792e-6,3.156514724580204785e-6,3.3976731247300604017e-6,3.6563979507854058118e-6,3.9339010171805394158e-6,4.2314725913513058091e-6,4.5504860985289220255e-6,4.8924030826534157691e-6,5.2587784356230157198e-6,5.6512659075690381388e-6,6.071623911330598909e-6,6.5217216347996970351e-6,7.0035454753146790574e-6,7.5192058107985795587e-6,8.0709441228680756994e-6,8.6611404876784849927e-6,9.2923214508200212846e-6,9.967168303140002239e-6,0.000010688525774934420469,0.000011459411166529745151,0.00001228302393386145192,0.000013162755748248969124,0.00001410220105016680142,0.000015105168117417808295,0.000016175690668726172002,0.000017318040024383616155,0.00001853673784620199306,0.000019836569479647433367,0.000021222597921654760479,0.00002270017843824465265,0.000024274973856688847993,0.000025952970557589211166,0.000027740495192853312311,0.000029644232156160821455,0.000031671241833119921254,0.000033828979658909425936,0.000036125316011788613228,0.00003856855697143108651,0.000041167465971599351906,0.000043931286377221191283,0.00004686976501645521491,0.000049993176698838969458,0.000053312349751096345086,0.000056838692602639360983,0.000060584221453230207317,0.000064561589055670083169,0.000068784114646749199484,0.000073265815060024524303,0.000078021437054285552077,0.000083066490891820605942,0.000088417285200803867818,0.000094090963156283324692,0.00010010554001435888442,0.00010647994203419473213,0.00011323404682250717192,0.00012038872513510525576,0.00012796588416993289911,0.00013598851238586431234,0.00014448072588123576744,0.00015346781636575224608,0.00016297630075898360227,0.00017303397244815770307,0.00018366995423736373063,0.00019491475301959355599,0.00020680031620226892919,0.00021936008991602325604,0.00023262907903552503635,0.0002466439090400417116,0.00026144288974024382273,0.00027706608089643614884,0.0002935553597519710628,0.00031095449050404592986,0.0003293091957324072834,0.00034866722980467610843,0.00036907845427506730867,0.0003905949152911988914,0.00041327092302146925839,0.00043716313311312105788,0.00046233063018860429559,0.00048883501338519695102,0.00051674048394003550863,0.00054611393481974809588,0.00057702504239076704292,0.00060954636012312358986,0.00064375341431709336611,0.00067972480183846561486,0.00071754228984445068228,0.00075729091747831808207,0.00079905909950677085141,0.00084293873186981063777,0.00088902529910843205933,0.00093741798363090480051,0.00098821977677365846372,0.0010415375916078796373,0.0010974823774378646204,0.001156169235931946876,0.0012177175388214381393,0.0012822510470974892265,0.0013498980316300945267,0.0014207913951276369571,0.0014950687953494024029,0.0015728727694773894565,0.001654350859547507298,0.0017396557388338980305,0.001828945339073646532,0.0019223829784125582824,0.002020137489946001681,0.0021223833507220346619,0.0022293008110671751632,0.002341076024088241023,0.0024579011751966876073,0.0025799746114948219875,0.002707500970856182468,0.0028406913105252540319,0.0029797632350545567543,0.0031249410233900089929,0.0032764557549083443492,0.0034345454342032657803,0.0035994551144099671843,0.0037714370188506582168,0.0039507506607768085312,0.0041376629609770001283,0.004332448363012558626,0.0045353889458365441349,0.0047467745335452386041,0.0049669028020049880651,0.0051960793820911645917,0.0054346179592701239237,0.0056828403692493708054,0.0059410766894157257619,0.006209665325776135167,0.0064889530951109021828,0.0067792953020445608492,0.0070810558107353907856,0.0073946071108806973216,0.0077203303777314831097,0.0080586155258070341749,0.0084098612559972576383,0.008774475095738361686,0.0091528734319456822301,0.0095454815363861564596,0.0099527335831721421444,0.010375072658058003871,0.010812950759221153707,0.011266828789210063988,0.011737176537743183412,0.012224472655044703153,0.012729204615405755451,0.013251868670662900529,0.013792969793289685268,0.014353021608801654747,0.014932546317180480212,0.015532074603028848233,0.016152145534174448402,0.01679330644844881259,0.017456112828374907947,0.018141128163506277527,0.018848923800170164568,0.019580078778377455304,0.020335179655673435736,0.021114820317715283316,0.021919601775374906859,0.0227501319481792072,0.023607025433914055715,0.024490903264233274266,0.025402392646129643489,0.026342126689141459372,0.027310744118185392028,0.028308888971924364529,0.02933721028659785079,0.030396361765261375051,0.031487001432402064565,0.032609791273917842858,0.033765396862469232365,0.034954486968234739425,0.036177733155123396422,0.037435809362521208083,0.038729391472671962445,0.040059156863817090419,0.041425783949243957934,0.042829951702417115974,0.04427233916839257909,0.04575362496174111297,0.047274486751232745106,0.048835600731561226576,0.050437641082413922541,0.052081279415219547733,0.053767184207933242703,0.055496020228245657071,0.057268447945629916099,0.05908512293266754391,0.060946695256121546281,0.062853808858251863047,0.064807100928895228063,0.066807201268858066004,0.068854731645197351635,0.070950305138990299273,0.07309452548621927831,0.075287986412423404339,0.077531270961792774056,0.079824950821405234578,0.082169585641328843105,0.084565722351335719946,0.087013894474994765924,0.089514621441931643676,0.092068407899064451367,0.094675743021642587549,0.097337099824934358903,0.10005293447742586092,0.10282368561640950884,0.10564977366685525769,0.10853160016447097194,0.11146954708387053618,0.11446397617277908736,0.11751522829321414915,0.12062362277058941369,0.12378945675169440144,0.12701300457250819566,0.13029451713680885461,0.13363422130654191658,0.13703231930491159719,0.14048898813315680503,0.14400437900197094344,0.14757861677851959987,0.15121179945000362926,0.15490399760470679743,0.15865525393145705141,0.16246558273841861648,0.16633496949211847846,0.17026337037759539274,0.17425071188054237046,0.1782968903922946329,0.18240177183849430805,0.18656519133224068144,0.19078695285251062562,0.19506682894860794062,0.19940456047137276796,0.20379985633185302292,0.20825239328810895951,0.21276181576078957615,0.21732773567808563269,0.22194973235062862446,0.22662735237686819933,0.23136010957942226179,0.23614748497285444223,0.24098892676329278134,0.24588385038026145383,0.25083163854105420091,0.25583164134793392879,0.26088317641839773589,0.26598552904870053231,0.27113795241078349204,0.27633966778270591481,0.28158986481263075648,0.2868877018163652025,0.29223230610840829832,0.29762277436640790689,0.30305817302887922923,0.30853753872598689636,0.3140598787431423152,0.31962417151711762604,0.32522936716432740356,0.33087438804087920536,0.33655812933394434262,0.34227945968395091558,0.34803722183705232157,0.35383023332727620563,0.35965728718771128007,0.3655171526900426897,0.37140857611170073789,0.37733028152984291317,0.3832809716413453576,0.38925932860793729097,0.39526401492557053677,0.40129367431707627576,0.40734693264712256468,0.41342239885844908401,0.41951866592832009891,0.42563431184410280688,0.43176790059684617809,0.43791798319170513781,0.44408309867402656087,0.45026177516988710702,0.45645253093984848046,0.46265387544467330042,0.46886431042172447063,0.4750823309707527788,0.48130642664776148169,0.4875350825656228727,0.49376678050011031717,0.5,0.50623321949988968283,0.5124649174343771273,0.51869357335223851831,0.5249176690292472212,0.53113568957827552937,0.53734612455532669958,0.54354746906015151954,0.54973822483011289298,0.55591690132597343913,0.56208201680829486219,0.56823209940315382191,0.57436568815589719312,0.58048133407167990109,0.58657760114155091599,0.59265306735287743532,0.59870632568292372424,0.60473598507442946323,0.61074067139206270903,0.6167190283586546424,0.62266971847015708683,0.62859142388829926211,0.6344828473099573103,0.64034271281228871993,0.64616976667272379437,0.65196277816294767843,0.65772054031604908442,0.66344187066605565738,0.66912561195912079464,0.67477063283567259644,0.68037582848288237396,0.6859401212568576848,0.69146246127401310364,0.69694182697112077077,0.70237722563359209311,0.70776769389159170168,0.7131122981836347975,0.71841013518736924352,0.72366033221729408519,0.72886204758921650796,0.73401447095129946769,0.73911682358160226411,0.74416835865206607121,0.74916836145894579909,0.75411614961973854617,0.75901107323670721866,0.76385251502714555777,0.76863989042057773821,0.77337264762313180067,0.77805026764937137554,0.78267226432191436731,0.78723818423921042385,0.79174760671189104049,0.79620014366814697708,0.80059543952862723204,0.80493317105139205938,0.80921304714748937438,0.81343480866775931856,0.81759822816150569195,0.8217031096077053671,0.82574928811945762954,0.82973662962240460726,0.83366503050788152154,0.83753441726158138352,0.84134474606854294859,0.84509600239529320257,0.84878820054999637074,0.85242138322148040013,0.85599562099802905656,0.85951101186684319497,0.86296768069508840281,0.86636577869345808342,0.86970548286319114539,0.87298699542749180434,0.87621054324830559856,0.87937637722941058631,0.88248477170678585085,0.88553602382722091264,0.88853045291612946382,0.89146839983552902806,0.89435022633314474231,0.89717631438359049116,0.89994706552257413908,0.9026629001750656411,0.90532425697835741245,0.90793159210093554863,0.91048537855806835632,0.91298610552500523408,0.91543427764866428005,0.91783041435867115689,0.92017504917859476542,0.92246872903820722594,0.92471201358757659566,0.92690547451378072169,0.92904969486100970073,0.93114526835480264837,0.933192798731141934,0.93519289907110477194,0.93714619114174813695,0.93905330474387845372,0.94091487706733245609,0.9427315520543700839,0.94450397977175434293,0.9462328157920667573,0.94791872058478045227,0.94956235891758607746,0.95116439926843877342,0.95272551324876725489,0.95424637503825888703,0.95572766083160742091,0.95717004829758288403,0.95857421605075604207,0.95994084313618290958,0.96127060852732803755,0.96256419063747879192,0.96382226684487660358,0.96504551303176526058,0.96623460313753076763,0.96739020872608215714,0.96851299856759793544,0.96960363823473862495,0.97066278971340214921,0.97169111102807563547,0.97268925588181460797,0.97365787331085854063,0.97459760735387035651,0.97550909673576672573,0.97639297456608594429,0.9772498680518207928,0.97808039822462509314,0.97888517968228471668,0.97966482034432656426,0.9804199212216225447,0.98115107619982983543,0.98185887183649372247,0.98254388717162509205,0.98320669355155118741,0.9838478544658255516,0.98446792539697115177,0.98506745368281951979,0.98564697839119834525,0.98620703020671031473,0.98674813132933709947,0.98727079538459424455,0.98777552734495529685,0.98826282346225681659,0.98873317121078993601,0.98918704924077884629,0.98962492734194199613,0.99004726641682785786,0.99045451846361384354,0.99084712656805431777,0.99122552490426163831,0.99159013874400274236,0.99194138447419296583,0.99227966962226851689,0.99260539288911930268,0.99291894418926460921,0.99322070469795543915,0.99351104690488909782,0.99379033467422386483,0.99405892331058427424,0.99431715963075062919,0.99456538204072987608,0.99480392061790883541,0.99503309719799501193,0.9952532254664547614,0.99546461105416345587,0.99566755163698744137,0.99586233703902299987,0.99604924933922319147,0.99622856298114934178,0.99640054488559003282,0.99656545456579673422,0.99672354424509165565,0.99687505897660999101,0.99702023676494544325,0.99715930868947474597,0.99729249902914381753,0.99742002538850517801,0.99754209882480331239,0.99765892397591175898,0.99777069918893282484,0.99787761664927796534,0.99797986251005399832,0.99807761702158744172,0.99817105466092635347,0.99826034426116610197,0.9983456491404524927,0.99842712723052261054,0.9985049312046505976,0.99857920860487236304,0.99865010196836990547,0.99871774895290251077,0.99878228246117856186,0.99884383076406805312,0.99890251762256213538,0.99895846240839212036,0.99901178022322634154,0.9990625820163690952,0.99911097470089156794,0.99915706126813018936,0.99920094090049322915,0.99924270908252168192,0.99928245771015554932,0.99932027519816153439,0.99935624658568290663,0.99939045363987687641,0.99942297495760923296,0.9994538860651802519,0.99948325951605996449,0.99951116498661480305,0.9995376693698113957,0.99956283686688687894,0.99958672907697853074,0.99960940508470880111,0.99963092154572493269,0.99965133277019532389,0.99967069080426759272,0.99968904550949595407,0.99970644464024802894,0.99972293391910356385,0.99973855711025975618,0.99975335609095995829,0.99976737092096447496,0.99978063991008397674,0.99979319968379773107,0.99980508524698040644,0.99981633004576263627,0.9998269660275518423,0.9998370236992410164,0.99984653218363424775,0.99985551927411876423,0.99986401148761413569,0.9998720341158300671,0.99987961127486489474,0.99988676595317749283,0.99989352005796580527,0.99989989445998564112,0.99990590903684371668,0.99991158271479919613,0.99991693350910817939,0.99992197856294571445,0.99992673418493997548,0.9999312158853532508,0.99993543841094432992,0.99993941577854676979,0.99994316130739736064,0.99994668765024890365,0.99995000682330116103,0.99995313023498354479,0.99995606871362277881,0.99995883253402840065,0.99996143144302856891,0.99996387468398821139,0.99996617102034109057,0.99996832875816688008,0.99997035576784383918,0.99997225950480714669,0.99997404702944241079,0.99997572502614331115,0.99997729982156175535,0.99997877740207834524,0.99998016343052035257,0.99998146326215379801,0.99998268195997561638,0.99998382430933127383,0.99998489483188258219,0.9999858977989498332,0.99998683724425175103,0.99998771697606613855,0.99998854058883347025,0.99998931147422506558,0.99999003283169686,0.99999070767854917998,0.99999133885951232152,0.99999192905587713192,0.99999248079418920142,0.99999299645452468532,0.9999934782783652003,0.9999939283760886694,0.99999434873409243096,0.99999474122156437698,0.99999510759691734658,0.99999544951390147108,0.99999576852740864869,0.99999606609898281946,0.99999634360204921459,0.99999660232687526994,0.9999968434852754198,0.99999706821507152596,0.99999727758432024709,0.99999747259531821558,0.99999765418839546309,0.99999782324550712164,0.99999798059363302486,0.99999812700799444329,0.99999826321509680865,0.99999838989560691491,0.99999850768707272774,0.99999861718649358991,0.9999987189527482765,0.9999988135088880319,0.99999890134430140912,0.9999989829167574313,0.99999905865433330532,0.99999912895723263768,0.99999919419949983292,0.99999925473063609542,0.99999931087712220523,0.99999936294385299789,0.99999941121548824688,0.99999945595772442508,0.99999949741849160783,0.9999995358290795751,0.99999957140519697371,0.99999960434796721151,0.9999996348448645747,0.99999966307059388615,0.99999968918791685646,0.99999971334842812081,0.9999997356932838026,0.99999975635388529959,0.99999977545252084936,0.99999979310296729835,0.99999980941105437201,0.99999982447519362268,0.99999983838687411672,0.99999985123112681223,0.99999986308695947419,0.99999987402776387383,0.99999988412169692421,0.99999989343203731352,0.99999990201751911146,0.99999990993264374244,0.99999991722797164151,0.99999992395039483511,0.99999993014339161844,0.99999993584726443497,0.99999994109936200013,0.99999994593428665148,0.99999995038408785081,0.9999999544784427098,0.99999995824482435991,0.99999996170865893876,0.99999996489347191981,0.99999996782102446873,0.99999997051144046891,0.99999997298332482018,0.99999997525387357808,0.9999999773389764665,0.99999997925331226412,0.99999998101043753411,0.99999998262286913785,0.99999998410216094565,0.99999998545897513217,0.99999998670314841944,0.99999998784375360787,0.99999998888915671394,0.99999998984707001282,0.99999999072460126544,0.99999999152829939113,0.99999999226419683051,0.99999999293784882702,0.99999999355436984093,0.99999999411846729535,0.99999999463447284094,0.99999999510637131335,0.9999999955378275461,0.99999999593221119059,0.99999999629261968486,0.99999999662189950313,0.9999999969226658092,0.99999999719732062841,0.99999999744806964519,0.99999999767693772557,0.99999999788578325756,0.99999999807631139556,0.99999999825008628909,0.99999999840854237072,0.99999999855299477233,0.99999999868464893471,0.99999999880460947004,0.99999999891388833342,0.99999999901341235496,0.99999999910403018059,0.99999999918651866626,0.99999999926158876687,0.99999999932989095832,0.99999999939202022843,0.99999999944852066952,0.9999999994998897034,0.9999999995465819673,0.9999999995890128863,0.9999999996275619574,0.9999999996625757673,0.9999999996943707649,0.999999999723235807,0.9999999997494344966,0.9999999997732073281,0.9999999997947736575,0.9999999998143335077,0.9999999998320692264,0.9999999998481470044,0.9999999998627182689,0.9999999998759209589,0.9999999998878806947,0.9999999998987118483,0.9999999999085185242,0.9999999999173954563,0.9999999999254288296,0.9999999999326970314,0.9999999999392713383,0.9999999999452165448,0.9999999999505915376,0.9999999999554498205,0.9999999999598399942,0.9999999999638061937,0.9999999999673884891,0.9999999999706232502,0.9999999999735434797,0.9999999999761791177,0.9999999999785573182,0.9999999999807027029,0.999999999982637591,0.9999999999843822097,0.9999999999859548852,0.9999999999873722177,0.9999999999886492394,0.9999999999897995599,0.9999999999908354969,0.9999999999917681966,0.9999999999926077422,0.9999999999933632531,0.9999999999940429749,0.9999999999946543609,0.9999999999952041472,0.9999999999956984194,0.9999999999961426746,0.9999999999965418767,0.999999999996900507,0.9999999999972226104,0.9999999999975118365,0.9999999999977714777,0.9999999999980045034,0.9999999999982135908,0.9999999999984011533,0.9999999999985693657,0.9999999999987201875,0.9999999999988553835,0.9999999999989765432,0.9999999999990850972,0.9999999999991823334,0.9999999999992694105,0.9999999999993473708,0.9999999999994171519,0.9999999999994795966,0.9999999999995354625,0.9999999999995854306,0.9999999999996301126,0.999999999999670058,0.9999999999997057601,0.999999999999737662,0.9999999999997661611,0.9999999999997916142,0.9999999999998143413,0.9999999999998346294,0.9999999999998527359,0.9999999999998688914,0.9999999999998833025,0.9999999999998961545,0.9999999999999076134,0.9999999999999178275,0.9999999999999269299,0.9999999999999350396,0.9999999999999422631,0.9999999999999486957,0.9999999999999544226,0.9999999999999595199,0.9999999999999640558,0.9999999999999680911,0.9999999999999716802,0.9999999999999748716,0.9999999999999777087,0.9999999999999802302,0.9999999999999824708,0.9999999999999844611,0.9999999999999862288,0.9999999999999877983,0.9999999999999891915,0.999999999999990428,0.999999999999991525,0.9999999999999924981,0.9999999999999933611,0.9999999999999941262,0.9999999999999948044,0.9999999999999954054,0.9999999999999959378,0.9999999999999964094,0.999999999999996827,0.9999999999999971968,0.999999999999997524,0.9999999999999978136,0.9999999999999980697,0.9999999999999982963,0.9999999999999984966,0.9999999999999986737,0.9999999999999988302,0.9999999999999989685,0.9999999999999990907,0.9999999999999991986,0.9999999999999992938,0.9999999999999993779,0.9999999999999994521,0.9999999999999995176,0.9999999999999995753,0.9999999999999996263,0.9999999999999996711,0.9999999999999997107,0.9999999999999997456,0.9999999999999997763,0.9999999999999998034,0.9999999999999998272,0.9999999999999998482,0.9999999999999998666,0.9999999999999998829,0.9999999999999998972,0.9999999999999999098,0.9999999999999999208,0.9999999999999999305,0.9999999999999999391,0.9999999999999999466,0.9999999999999999531,0.9999999999999999589,0.999999999999999964,0.9999999999999999685,0.9999999999999999724,0.9999999999999999758,0.9999999999999999788,0.9999999999999999815,0.9999999999999999838,0.9999999999999999858,0.9999999999999999876,0.9999999999999999892,0.9999999999999999905,0.9999999999999999917,0.9999999999999999928,0.9999999999999999937,0.9999999999999999945,0.9999999999999999952,0.9999999999999999958,0.9999999999999999963,0.9999999999999999968,0.9999999999999999972,0.9999999999999999976,0.9999999999999999979,0.9999999999999999981,0.9999999999999999984,0.9999999999999999986,0.9999999999999999988,0.9999999999999999989,0.9999999999999999991,0.9999999999999999992,0.9999999999999999993,0.9999999999999999994,0.9999999999999999995,0.9999999999999999995,0.9999999999999999996,0.9999999999999999997,0.9999999999999999997,0.9999999999999999997,0.9999999999999999998,0.9999999999999999998,0.9999999999999999998,0.9999999999999999999,0.9999999999999999999,0.9999999999999999999,0.9999999999999999999,0.9999999999999999999,0.9999999999999999999,0.9999999999999999999,0.9999999999999999999,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

int main(int, char**)
{
   random_state_test state;

   {
std::cout << "gaussint err on [-20,10]: " << std::endl;
      double maxerr = 0;
      int i = 0;
      for (double x = -20.0; x <= -10.0; x += 1.0/64.0, i++) {
         double p1 = gaussint(x);
         double p2 = testtab[i];
         double err = relerr(p1, p2);
         if (err > maxerr) {
            maxerr = err;
            std::cout << "at " << x << " rel err " << err << std::endl;
            std::cout << "              my: " << format(p1) << " or " << format_hex(p1) << std::endl;
            std::cout << "             ref: " << format(p2) << " or " << format_hex(p2) << std::endl;
         }
      }

std::cout << "exp err: " << std::endl;
      maxerr = 0;
      for (double x = -20.0; x <= 1.0; x += 1.0/(15*1024.0)) {
         double p1 = std::exp(x);
         double p2 = my_exp(x);
         double err = relerr(p1, p2);
         if (err > maxerr) {
            maxerr = err;
            std::cout << "at " << format(8*1024*x) << "/(8*1024)  rel err " << err << std::endl;
            std::cout << "              my: " << format_hex(p1) << std::endl;
            std::cout << "             ref: " << format_hex(p2) << std::endl;
         }
      }
std::cout << "exp err done " << std::endl;
   }

   {
      double t1 = 0, t2 = 0;
      for (int i = 0; i < 10; i++) {
         t1 += testf(state, gaussint);
         t2 += testf(state, gaussint_ref);
         std::cout << " my: " << format_fixed(t1/(i+1),2,2) <<
                      ", ref: " << format_fixed(t2/(i+1),2,2) <<
                      ", ratio: " << format_fixed(t2/t1,2,2) << std::endl;
      }
   }

   if (1) {
      uword t1 = get_timestamp();
      double x = 0.4;
      double y = 2.000000000000001;
      for (int i = 0; i < 100000; i++) {
         x = y-x;
         x = y-x;
         x = y-x;
      }
      uword t2 = get_timestamp();
      std::cout << "cycles: " << sword(t2-t1)/double(3*100000) << std::endl;
      std::cout << "x: " << format(x) << std::endl;
   }

   return 0;

try {
   fit_tab(10);
   test_fft();
   profile_fft(1,3+1);
   profile_fft(5,1+1);
   profile_fft(3,2+1);
   profile_fft(1,3+2);
   profile_fft(5,1+2);
   profile_fft(3,2+2);
   profile_fft(1,3+3);
   profile_fft(5,1+3);
   profile_fft(3,2+3);

//   test_fmpz(state);
/*
   test_fmpq(state);
   test_arb();
   test_prime_field_mpoly<fp_ring>(state);
   test_prime_field_poly<fp_ring>(state);
   test_mpoly(state);
*/

   if (1) {
      uword nreps = 1000000;
      std::complex<double> res(0);

      uword t1 = get_timestamp();
      for (uword i = 0; i < nreps; i++) {
         res += d_cis2pi_frac(1267, 5, 11);
      }
      uword t2 = get_timestamp();
      for (uword i = 0; i < nreps; i++) {
         res += my_cis2pi_frac(1267, 5<<11);
      }
      uword t3 = get_timestamp();
      std::cout << "form 1: " << (double)(t2-t1)/(nreps) << std::endl;
      std::cout << "form 2: " << (double)(t3-t2)/(nreps) << std::endl;
   }


   if (0) {
      uword nreps = 100000;
      uword n = 200;
      uword* a = new uword[n];
      uword* b = new uword[n];
      uword* c = new uword[n];
      uword m1 = 0x7F0FFFFFFFFFFFFF;
      uword m2 = 0x7F0FFFFFF7FFFFFF;

      for (int i = 0; i < n; i++)
         b[i] = c[i] = -uword(1);

      uword t1 = get_timestamp();
      for (uword i = 0; i < nreps; i++) {
         mpn_addaddmul_1msb0(a, b, c, n-(i%8), m1, m2);
      }
      uword t2 = get_timestamp();
      for (uword i = 0; i < nreps; i++) {
         my_mpn_mul_1(a, b, n-(i%8), m1);
         my_mpn_addmul_1(a, c, n-(i%8), m2);
      }
      uword t3 = get_timestamp();
      std::cout << "form 1: " << (double)(t2-t1)/(n*nreps) << std::endl;
      std::cout << "form 2: " << (double)(t3-t2)/(n*nreps) << std::endl;
   }

   if (0) {
      uword c[100];
      uword d[100];
      uword a = -uword(1);
      uword b = -uword(1);
      for (int n = 1; n < 20; n++)
      {
         for (int i = 0; i < n; i++)
            c[i] = d[i] = -uword(1);

         c[n] = mpn_mul_1c(c, c, n, a, b);

         d[n] = mpn_mul_1(d, d, n, a);
         d[n] += mpn_add_1(d, d, n, b);

         for (int i = 0; i <= n; i++)
            TEST(c[i] == d[i], i);
      }
   }

   if (0) {
      fmpz a, b, q, r;
      for (uword i = 19; i <= 25; i++)
      {
         fmpz_fib_ui(a, ui_pow2(i+1));
         fmpz_fib_ui(b, ui_pow2(i)+100);
         fmpz_tdiv_qr(q, r, a, b);
         uword t1 = get_ms();
         fmpz_tdiv_qr(q, r, a, b);
         uword t2 = get_ms();
         std::cout << "bits " << fmpz_bits(a) << "/" << fmpz_bits(b) << ": " << t2-t1 << std::endl;
      }
   }

   if (0) {
      fmpz a, b, c;
      for (uword i = 5; i < 23; i++)
      {
         fmpz_fib_ui(a, ui_pow2(i));
         fmpz_fib_ui(b, ui_pow2(i)+1);
         gcd(ZZ, c, a, b);
         uword t1 = get_ms();
         gcd(ZZ, c, a, b);
         uword t2 = get_ms();
         std::cout << "bits " << fmpz_bits(a) << " gcd: " << t2-t1 << std::endl;
      }
   }

#if 0
   if (0) {

      uword max_len = 3000;

      uword* aa = new uword[max_len];
      uword* zz = new uword[2*max_len];

      for (uword i = 0; i < max_len; i++)
         aa[i] = ~i;

      for (int jj = 0; jj < 2; jj++)
      for (uword an = 1; an < max_len; an += 1 + an/8)
      {
         uword n = an*128;
         uword nreps = 1 + 3000000000/(n*clog2(n));
         uword t1 = get_ms();
         for (uword j = 0; j < nreps; j++)
            my_mpn_sqr(zz, aa, an);
         uword t2 = get_ms();
         std::cout << format_fixed(log2(n),2,2) << ": " << double(t2-t1)*1e8/(nreps*log2(n)*n) << std::endl;
      }

      delete[] aa;
      delete[] zz;
   }
#endif

   return 0;
}
catch (const char* e)
{
   std::cout << "oops: " << e << std::endl;
   return 1;
}
}


#if 0
#include <map>
   {
      using mat22 = std::tuple<sword,sword,sword,sword>;
      std::map<mat22,uword> T, S[10];

      auto r = T.insert_or_assign(mat22(1,0,0,1), 0);
      TEST(r.second);
      r = S[0].insert_or_assign(mat22(1,0,0,1), 0);
      TEST(r.second);

      for (int len = 1; len <= 3; len++) {
         for (auto& m : S[len-1]) {
// 1 -> -2, -2 -> 2
            for (sword q = 1; q <= 2; q = -q-(q>0)) {
               TEST (!(q == 0 || q == -1));
               sword m11 = std::get<0>(m.first);
               sword m12 = std::get<1>(m.first);
               sword m21 = std::get<2>(m.first);
               sword m22 = std::get<3>(m.first);
               uword mhist = m.second;

//std::cout << "q: " << q << ", m11="<<m11<< ", m12="<<m12<<", m21="<<m21<<", m22="<<m22;


               // m11 m12  q 1
               // m21 m22  1 0
               sword a11 = m11*q+m12;
               sword a12 = m11;
               sword a21 = m21*q+m22;
               sword a22 = m21;
               if (a21 < 0 || (a21 == 0 && a11 < 0)) {
                  a11 = -a11;
                  a12 = -a12;
                  a21 = -a21;
                  a22 = -a22;
               }
               uword ahist = mhist + ((q&255) << (8*(len-1)));

               r = T.insert_or_assign(mat22(a11,a12,a21,a22), ahist);
if (!r.second)
std::cout << "duplicate";
std::cout << "len: " << " ahist: " << format_hex(ahist) << " a: " << a11 << ", " << a12 << ", " << a21 << ", " << a22 << std::endl;

               r = S[len].insert_or_assign(mat22(a11,a12,a21,a22), ahist);
//               TEST(r.second, len, a11,a12,a21,a22);
            }
         }
      }
return 0;
   }
#endif
