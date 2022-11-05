#include <typeinfo>
#include "fmpz.h"
#include "test_helpers.h"
#include "timing.h"
#include "nmod_poly.h"


template <typename RingType>
void test_finite_field_poly(RingType& R, rand_state& state)
{
    typename RingType::poly_t a, b, c, d, e, q, r;
    typename RingType::poly_product_t f;

    std::cout << typeid(typename RingType::poly_t).name() << "_sub_neg_mul..." << std::flush;

    randtest(R, a, state, 8);
    randtest(R, b, state, 4);
    sub(R, c, a, b);
    sub(R, d, b, a);
    add(R, c, d);
    TEST1(is_zero(R, c), with(R, c));

    std::cout << "PASS" << std::endl;

    std::cout << typeid(typename RingType::poly_t).name() << "_poly_divrem..." << std::flush;

    randtest(R, a, state, 8);
    randtest(R, b, state, 4);

    try {
        divrem(R, q, r, a, b);
        mul(R, c, q, b);
        add(R, c, r);
        TEST4(equal(R, a, c), with(R, q), with(R, r), with(R, a), with(R, b));
    }
    catch (const with<RingType, typename RingType::elem_t>& oops)
    {
        std::cout << std::endl;
        std::cout << "zero divisor: " << oops << std::endl;
        std::cout << "     in ring: " << oops.parent() << std::endl;
    }

    std::cout << "PASS" << std::endl;

    std::cout << "factoring: " << std::endl;

    randtest(R, a, state, 5);
    randtest(R, b, state, 4);
    pow(R, c, a, 2);
    pow(R, d, b, 3);
    mul(R, q, c, d);
    factor_squarefree(R, f, q);
    std::cout << with(R, q) << " = " << with(R, f) << std::endl;

    set(R, q, "(x+1)^4*(x+2)^5");
    factor_squarefree(R, f, q);
    std::cout << with(R, q) << " = " << with(R, f) << std::endl;

    std::cout << "parsing: " << std::endl;

    try {
        set(R, q, "x+");
        TEST(false && "unreachable");
    }
    catch (const char*& oops)
    {
        std::cout << oops << std::endl;
    }
    
}

void test_fmpz(rand_state& state)
{
    fmpz_ring ZZ;
    fmpz a, b, c, d, e, s, t, q, r;

    std::cout << "fmpz_add_sub_neg_mul..." << std::flush;

    for (ulimb max_bits = 200; max_bits <= pow2(18); max_bits += 1 + max_bits/2)
    for (ulimb i = 0; i < 2+10000000/max_bits; i++)
    {
        randtest(ZZ, a, state, max_bits);
        randtest(ZZ, b, state, max_bits);
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
        randtest(ZZ, b, state, max_bits);
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

        randtest(ZZ, a, state, max_bits);
        randtest(ZZ, b, state, max_bits);
        randtest(ZZ, c, state, max_bits);
        add(ZZ, t, b, c);
        mul(ZZ, d, a, t);
        mul(ZZ, t, a, b);
        mul(ZZ, s, a, c);

        add(ZZ, e, t, s);
        TEST2(equal(ZZ, d, e), d, e);

        randtest(ZZ, a, state, max_bits);
        randtest(ZZ, b, state, max_bits);
        add(ZZ, t, a, b);
        mul(ZZ, d, t, t);
        mul(ZZ, e, a, a);
        mul(ZZ, t, b, b);
        add(ZZ, e, e, t);
        mul(ZZ, t, a, b);
        add(ZZ, e, e, t);
        add(ZZ, e, e, t);
        TEST2(equal(ZZ, d, e), d, e);

        randtest(ZZ, a, state, max_bits);
        randtest(ZZ, b, state, max_bits);
        sub(ZZ, t, a, b);
        mul(ZZ, d, t, t);
        mul(ZZ, e, a, a);
        mul(ZZ, t, b, b);
        add(ZZ, e, e, t);
        mul(ZZ, t, a, b);
        sub(ZZ, e, e, t);
        sub(ZZ, e, e, t);
        TEST2(equal(ZZ, d, e), d, e);
    }

    std::cout << "PASS" << std::endl;

    std::cout << "fmpz_div..." << std::flush;

    for (ulimb max_bits = 200; max_bits <= pow2(19); max_bits += 1 + max_bits/2)
    for (ulimb i = 0; i < 2+10000000/max_bits; i++)
    {
        randtest(ZZ, a, state, max_bits);
        randtest_not_zero(ZZ, b, state, max_bits);

        mul(ZZ, c, a, b);
        divexact(ZZ, d, c, b);
        TEST2(equal(ZZ, a, d), a, d);

        mul(ZZ, c, a, b);
        divexact(ZZ, c, c, b);
        TEST2(equal(ZZ, a, c), a, c);

        fmpz_tdiv_qr(q, r, a, b);
        TEST(fmpz_cmpabs(r, b) < 0);
        mul(ZZ, t, q, b);
        add(ZZ, t, t, r);
        TEST2(equal(ZZ, t, a), t, a);
        TEST3(r.sign()*a.sign() >= 0, r, a, b);

        fmpz_fdiv_qr(q, r, a, b);
        TEST2(fmpz_cmpabs(r, b) < 0, r, b);
        mul(ZZ, t, q, b);
        add(ZZ, t, t, r);
        TEST4(equal(ZZ, t, a), q, r, a, b);
        TEST3(r.sign()*b.sign() >= 0, r, a, b);

        fmpz_cdiv_qr(q, r, a, b);
        TEST(fmpz_cmpabs(r, b) < 0);
        mul(ZZ, t, q, b);
        add(ZZ, t, t, r);
        TEST2(equal(ZZ, t, a), t, a);
        TEST4(r.sign()*b.sign() <= 0, q, r, a, b);
    }

    std::cout << "PASS" << std::endl;

    std::cout << "fmpz_gcd..." << std::flush;

    for (ulimb max_bits = 200; max_bits <= pow2(16); max_bits += 1 + max_bits/2)
    for (ulimb i = 0; i < 2 + 3000000/max_bits; i++)
    {
        if (i % 2)
        {
            randtest(ZZ, a, state, max_bits);
            randtest(ZZ, b, state, max_bits);
            randtest(ZZ, e, state, max_bits);
            mul_2exp(ZZ, a, a, state.get_mod(100));
            mul_2exp(ZZ, b, b, state.get_mod(100));
            mul_2exp(ZZ, e, e, state.get_mod(100));
            mul(ZZ, a, a, e);
            mul(ZZ, b, b, e);
        }
        else
        {
            randtest_not_zero(ZZ, e, state, 20 + max_bits/20);
            zero(ZZ, a);
            set(ZZ, b, e);
            while (fmpz_bits(a) < max_bits && fmpz_bits(b) < max_bits)
            {
                randtest(ZZ, c, state, 20 + max_bits/20);
                mul(ZZ, d, c, b);
                add(ZZ, a, a, d);
                swap(ZZ, a, b);
            }
        }
        gcd(ZZ, d, a, b);
        TEST3(divisible(ZZ, a, d), a, b, d);
        TEST3(divisible(ZZ, b, d), a, b, d);
        TEST2(divisible(ZZ, d, e), d, e);
        if (!is_zero(ZZ, d))
        {
            divexact(ZZ, a, a, d);
            divexact(ZZ, b, b, d);
            gcd(ZZ, e, a, b);
            TEST3(is_one(ZZ, e), e, a, b);
        }

        randtest(ZZ, a, state, max_bits);
        randtest(ZZ, b, state, max_bits);
        randtest(ZZ, c, state, max_bits);
        randtest(ZZ, e, state, max_bits);
        mul(ZZ, a, a, e);
        mul(ZZ, b, b, e);
        mul(ZZ, c, c, e);
        gcd(ZZ, d, a, b, c);
        TEST4(divisible(ZZ, a, d), a, b, c, d);
        TEST4(divisible(ZZ, b, d), a, b, c, d);
        TEST4(divisible(ZZ, c, d), a, b, c, d);
        TEST2(divisible(ZZ, d, e), d, e);
        if (!is_zero(ZZ, d))
        {
            divexact(ZZ, a, a, d);
            divexact(ZZ, b, b, d);
            divexact(ZZ, c, c, d);
            gcd(ZZ, e, a, b, c);
            TEST4(is_one(ZZ, e), e, a, b, c);
        }
    }

    std::cout << "PASS" << std::endl;

    std::cout << "fmpz_gcdx..." << std::flush;

    for (ulimb max_bits = 20; max_bits <= pow2(15); max_bits += 1 + max_bits/2)
    for (ulimb i = 0; i < 2 + 500000/max_bits; i++)
    {
        if (i % 2)
        {
            randtest(ZZ, a, state, max_bits);
            randtest(ZZ, b, state, max_bits);
            randtest(ZZ, e, state, max_bits);
            mul_2exp(ZZ, a, a, state.get_mod(100));
            mul_2exp(ZZ, b, b, state.get_mod(100));
            mul_2exp(ZZ, e, e, state.get_mod(100));
            mul(ZZ, a, a, e);
            mul(ZZ, b, b, e);
        }
        else
        {
            randtest_not_zero(ZZ, e, state, 20 + max_bits/20);
            zero(ZZ, a);
            set(ZZ, b, e);
            while (fmpz_bits(a) < max_bits && fmpz_bits(b) < max_bits)
            {
                randtest(ZZ, c, state, 20 + max_bits/20);
                mul(ZZ, d, c, b);
                add(ZZ, a, a, d);
                swap(ZZ, a, b);
            }
        }
        gcdx(ZZ, d, s, t, a, b);
        TEST2(divisible(ZZ, b, d), b, d);
        TEST2(divisible(ZZ, b, d), d, e);
        mul(ZZ, q, a, s);
        mul(ZZ, r, b, t);
        add(ZZ, q, q, r);
        TEST5(equal(ZZ, q, d), d, s, t, a, b);
    }

    std::cout << "PASS" << std::endl;
}


int main(int, char**)
{
    rand_state state;

    nmod_ring R(fmpz(3));
    test_finite_field_poly(R, state);

    test_fmpz(state);

#if 0
    if (0) {
        for (ulimb i = 19; i <= 25; i++)
        {
            fmpz_fib_ui(a, pow2(i+1));
            fmpz_fib_ui(b, pow2(i)+100);
            fmpz_tdiv_qr(q, r, a, b);
            ulimb t1 = GetMS();
            fmpz_tdiv_qr(q, r, a, b);
            ulimb t2 = GetMS();
            std::cout << "bits " << fmpz_bits(a) << "/" << fmpz_bits(b) << ": " << t2-t1 << std::endl;
        }
    }

    if (0) {
        for (ulimb i = 5; i < 25; i++)
        {
            fmpz_fib_ui(a, pow2(i));
            fmpz_fib_ui(b, pow2(i)+1);
            gcd(ZZ, c, a, b);
            ulimb t1 = GetMS();
            gcd(ZZ, c, a, b);
            ulimb t2 = GetMS();
            std::cout << "bits " << fmpz_bits(a) << " gcd: " << t2-t1 << std::endl;
        }
    }

    if (0) {

        ulimb max_len = 3000;

        ulimb* aa = new ulimb[max_len];
        ulimb* zz = new ulimb[2*max_len];

        for (ulimb i = 0; i < max_len; i++)
            aa[i] = ~i;

        for (int jj = 0; jj < 2; jj++)
        for (ulimb an = 1; an < max_len; an += 1 + an/8)
        {
            ulimb n = an*128;
            ulimb nreps = 1 + 3000000000/(n*clog2(n));
            ulimb t1 = GetMS();
            for (ulimb j = 0; j < nreps; j++)
                my_mpn_sqr(zz, aa, an);
            ulimb t2 = GetMS();
            std::cout << format_fixed(log2(n),2,2) << ": " << double(t2-t1)*1e8/(nreps*log2(n)*n) << std::endl;
        }

        delete[] aa;
        delete[] zz;
    }
#endif

    return 0;
}

