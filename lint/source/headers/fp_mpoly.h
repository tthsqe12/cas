#pragma once

#include "buffers.h"
#include "fp.h"
#include "generic/generics.h"
#include "mpoly.h"

struct fp_mpoly_ring {
    mpoly_ctx m;
    fp_ring c;

    fp_mpoly_ring(fmpz&& mo, ulimb nvars) : c(std::move(mo)), m(nvars) {}
    fp_mpoly_ring(fmpzc mo, ulimb nvars) : c(mo), m(nvars) {}
};

struct fp_mpolyi;

struct fp_mpoly {
    typedef fp_mpolyi input_t;

    mpoly m;
    chunk<ulimb> c;

    ulimb length() const {return m.length();}
    void set_length(ulimb n) {m.set_length(n);}
    ulimb bits() const {return m.bits();}
    void set_bits(ulimb n) {m.set_bits(n);}
    ulimb* exps() {return m.data();}
    const ulimb* exps() const {return m.data();}
    ulimb* coeffs() {return c.data();}
    const ulimb* coeffs() const {return c.data();}
};

struct fp_mpolyi {
    ulimb len;
    ulimb* exps;
    const ulimb* coeffs;
    fp_mpolyi(const fp_mpoly& a) : len(a.length()), exps((ulimb*)a.exps()), coeffs(a.coeffs()) {}
};

inline bool fp_mpoly_is_zero(fp_mpoly_ring& ctx, const fp_mpoly& a) {return a.length() == 0;}
inline void fp_mpoly_zero(fp_mpoly_ring& ctx, fp_mpoly& x) {x.m.zero(ctx.m);}

void fp_mpoly_sort_terms(
    fp_mpoly_ring& ctx,
    fp_mpoly& X);

void fp_mpoly_random_up_to_length_and_degree_bound(
    random_state& state,
    fp_mpoly_ring& ctx,
    fp_mpoly& X,
    ulimb l,
    ulimb d);

void fp_mpoly_random_up_to_length_and_degree_bits(
    random_state& state,
    fp_mpoly_ring& ctx,
    fp_mpoly& X,
    ulimb l,
    ulimb dbits);


inline void random_up_to_length_and_degree_bound(
    random_state& state,
    fp_mpoly_ring& ctx,
    fp_mpoly& X,
    ulimb l,
    ulimb d) {fp_mpoly_random_up_to_length_and_degree_bound(state, ctx, X, l, d);}

inline void random_up_to_length_and_degree_bits(
    random_state& state,
    fp_mpoly_ring& ctx,
    fp_mpoly& X,
    ulimb l,
    ulimb dbits) {fp_mpoly_random_up_to_length_and_degree_bits(state, ctx, X, l, dbits);}


std::ostream& fp_mpoly_write(std::ostream& o, const fp_mpoly_ring& R, const fp_mpoly& a, const char* var);

inline std::ostream& operator<<(std::ostream& o, const with<fp_mpoly_ring, fp_mpoly>& e) {
    return fp_mpoly_write(o, e.parent(), e.elem(), "x");
}

bool fp_mpoly_is_canonical(const fp_mpoly_ring& ctx, const fp_mpoly& a);


void fp_mpoly_combine_like_terms(fp_mpoly_ring& ctx, fp_mpoly& X);
void fp_mpoly_sort_terms(fp_mpoly_ring& ctx, fp_mpoly& X);



bool fp_mpoly_equal(fp_mpoly_ring& ctx, const fp_mpoly& B, const fp_mpoly& C);
void fp_mpoly_add(fp_mpoly_ring& ctx, fp_mpoly& A, const fp_mpoly& B, const fp_mpoly& C);
void fp_mpoly_sub(fp_mpoly_ring& ctx, fp_mpoly& A, const fp_mpoly& B, const fp_mpoly& C);
void fp_mpoly_mul(fp_mpoly_ring& ctx, fp_mpoly& A, const fp_mpoly& B, const fp_mpoly& C);
bool fp_mpoly_divides(fp_mpoly_ring& ctx, fp_mpoly& Q, const fp_mpoly& A, const fp_mpoly& B);

inline bool equal(fp_mpoly_ring& ctx, const fp_mpoly& a, const fp_mpoly& b) {return fp_mpoly_equal(ctx, a, b);}
inline bool is_zero(fp_mpoly_ring& ctx, const fp_mpoly& a) {return fp_mpoly_is_zero(ctx, a);}
inline void zero(fp_mpoly_ring& ctx, fp_mpoly& x) {fp_mpoly_zero(ctx, x);}
inline void add(fp_mpoly_ring& ctx, fp_mpoly& x, const fp_mpoly& a, const fp_mpoly& b) {fp_mpoly_add(ctx, x, a, b);}
inline void sub(fp_mpoly_ring& ctx, fp_mpoly& x, const fp_mpoly& a, const fp_mpoly& b) {fp_mpoly_sub(ctx, x, a, b);}
inline void mul(fp_mpoly_ring& ctx, fp_mpoly& x, const fp_mpoly& a, const fp_mpoly& b) {fp_mpoly_mul(ctx, x, a, b);}
inline bool divides(fp_mpoly_ring& ctx, fp_mpoly& x, const fp_mpoly& a, const fp_mpoly& b) {return fp_mpoly_divides(ctx, x, a, b);}

