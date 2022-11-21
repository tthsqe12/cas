#pragma once

#include "types.h"
#include "limbs.h"

struct random_state {
    ulimb cur_limb;
    unsigned int bits_left = 0;

    virtual ulimb next_limb() {FLINT_ASSERT_ALWAYS(false && "don't use this"); return 0;}
    virtual void seed(const uint8_t* data, size_t len) {}

    bool get_bit();
    ulimb get_limb();
    ulimb get_bits(unsigned int count);
    ulimb get_mod(ulimb m);
};

struct random_state_xorshift64 : random_state {
    ulimb s;

    ulimb next_limb() {
        ulimb x = s;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        return s = x;
    }

    void seed(const ulimb* data, ulimb len) {
        s = UWORD(0xd2c98b5bec27f9ae);
        for (ulimb i = 0; i < len; i++)
            s += data[i];
        s += (s == 0);
    }

    random_state_xorshift64() {seed(nullptr, 0);}
};

struct random_state_xoshiro256starstar : random_state {
    ulimb s[4];

    ulimb next_limb() {
        ulimb res = rotate_left(s[1]*5, 7)*9;
        ulimb t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotate_left(s[3], 45);
        return res;
    }

    void seed(const ulimb* data, ulimb len) {
        s[0] = UWORD(0x3852b46af8a5fc1c);
        s[1] = UWORD(0x24cc941b22b1176a);
        s[2] = UWORD(0x5f51122dccce526d);
        s[3] = UWORD(0xccb6c5550042970d);
        for (ulimb i = 0; i < len; i++)
            s[i%4] += data[i];
        s[0] += (s[0] | s[1] | s[2] | s[3]) == 0;
    }

    random_state_xoshiro256starstar() {seed(nullptr, 0);}
};


struct random_state_test : random_state {
    ulimb s[4];

    ulimb next_limb() {
        ulimb res = ((s[0] & s[2]) & (rotate_left(s[1]*5, 7)*9)) ^ sign_extend(s[3]);
        ulimb t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotate_left(s[3], 45);
        return res;
    }

    void seed(const ulimb* data, ulimb len) {
        s[0] = UWORD(0x3852b46af8a5fc1c);
        s[1] = UWORD(0x24cc941b22b1176a);
        s[2] = UWORD(0x5f51122dccce526d);
        s[3] = UWORD(0xccb6c5550042970d);
        for (ulimb i = 0; i < len; i++)
            s[i%4] += data[i];
        s[0] += (s[0] | s[1] | s[2] | s[3]) == 0;
    }

    random_state_test() {seed(nullptr, 0);}
};

struct random_state_cell5120 : random_state {
    unsigned int off;
    ulimb D[80];

    ulimb next_limb();
    void seed(const ulimb* data, ulimb len);

    random_state_cell5120();
    void startover();
    void shuffle();
};

