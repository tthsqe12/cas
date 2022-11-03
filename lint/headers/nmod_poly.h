#pragma once

#include "types.h"
#include "ulimb_extras.h"
#include "fmpz.h"

template <typename S, typename T>
struct with {
    const S* _parent;
    T _elem;

    with(const S& p, const T& e) : _parent(&p), _elem(e) {}
    with(const S& p, T&& e) : _parent(&p), _elem(e) {}

    const S& parent() const {return *_parent;}
    T& elem() {return _elem;}
    const T& elem_const() const {return _elem;}
};

// there is an passed parameter stride N which is the size of each element
template <class T>
struct vec_strided {
    T* _data;
    size_t _alloc;

    vec_strided() : _data(nullptr), _alloc(0) {};
    ~vec_strided() {std::free(_data);}

    vec_strided(const vec_strided<T>& other) {
        _alloc = other._alloc;
        _data = reinterpret_cast<T*>(std::malloc(_alloc*sizeof(T)));
        for (ulimb i = 0; i < _alloc; i++)
            _data[i] = other._data[i]; // TODO placement new?
    }

    vec_strided(vec_strided<T>&& other) {
        _data = other._data;
        _alloc = other._alloc;
        other._data = nullptr;
        other._alloc = 0;
    }

    T& operator [](size_t i) const {return _data[i];}
    T* at(size_t i, size_t N) const {assert(i < _alloc); return _data + i*N;}
    T* data() const {return _data;}

    T* fit_alloc(size_t n, size_t N) {
        if (n*N <= _alloc)
            return _data;
        size_t newalloc = std::max(n*N, _alloc + _alloc/2);
        T* newdata = reinterpret_cast<T*>(std::malloc(newalloc*sizeof(T)));
        for (ulimb i = 0; i < _alloc; i++)
            newdata[i] = _data[i]; // TODO placement new?
        std::free(_data);
        // TODO new elements?
        _alloc = newalloc;
        _data = newdata;
        return _data;
    }

    T* fit_alloc_destroy(size_t n, size_t N) {
        if (n*N <= _alloc)
            return _data;
        size_t newalloc = std::max(n*N, _alloc + _alloc/2);
        T* newdata = reinterpret_cast<T*>(std::realloc(_data, newalloc*sizeof(T)));
        // TODO new elements?
        _alloc = newalloc;
        _data = newdata;
        return _data;
    }
};


struct nmod_poly {
    ulimb len;
    vec_strided<ulimb> coeffs;

    slimb degree() const {return slimb(len)-1;}
    ulimb length() const {return len;}
    void set_length(ulimb n) {len = n;}
    ulimb* data() const {return coeffs.data();}
    ulimb* fit_alloc(ulimb n, ulimb N) {return coeffs.fit_alloc(n, N);}
    ulimb* fit_alloc_destroy(ulimb n, ulimb N) {return coeffs.fit_alloc_destroy(n, N);}
    void normalize(ulimb len, ulimb N);
};

// immutable thread safe base
struct nmod_ring_base {
    fmpzc modulus;
    ulimb stride;
    const ulimb* modulus_limbs;

    void complete_modulus_data() {
        FLINT_ASSERT_ALWAYS(modulus.cmp_small(1) > 0);
        stride = modulus.size();
        modulus_limbs = modulus.is_ptr() ? modulus.limbs() :
                                 reinterpret_cast<const ulimb*>(&modulus.data);
    }

    nmod_ring_base(fmpz&& m) : modulus(m.steal_data()) {
        complete_modulus_data();
    }

    nmod_ring_base(fmpzc m) : modulus(m.copy_abs()) {
        complete_modulus_data();
    }

    nmod_ring_base(const nmod_ring_base& R) : modulus(R.modulus.copy_abs()) {
        complete_modulus_data();
    }

    ~nmod_ring_base() {
        modulus.clear();
    }
};

// mutable non thread safe
struct nmod_ring : nmod_ring_base {
    typedef nmod_poly poly_elem_t;

    ulimb* tempN;
    ulimb* tempNb;
    ulimb* tempNp1;
    ulimb* temp2N;
    ulimb* temp2Np1;

    void complete_tmp_data() {
        ulimb N = stride;
        ulimb* p = my_alloc<ulimb>(N+N+N+1+2*N+2*N+1);
        tempN    = p; p += N;
        tempNb   = p; p += N;
        tempNp1  = p; p += N+1;
        temp2N   = p; p += 2*N;
        temp2Np1 = p; p += 2*N+1;
    }

    nmod_ring(fmpz&& m) : nmod_ring_base(std::move(m)) {
        complete_tmp_data();
    }

    nmod_ring(fmpzc m) : nmod_ring_base(m) {
        complete_tmp_data();
    }

    ~nmod_ring() {
        my_free(tempN);
    }
};


void nmod_add(ulimb* x, const ulimb* a, const ulimb* b, const nmod_ring_base& R);
void nmod_sub(ulimb* x, const ulimb* a, const ulimb* b, const nmod_ring_base& R);
void nmod_mul(ulimb* x, const ulimb* a, const ulimb* b, const nmod_ring_base& R);



bool nmod_poly_is_canonical(const nmod_ring_base& R, const nmod_poly& a);
void nmod_poly_write(std::ostream& o, const nmod_ring_base& R, const nmod_poly& a, const char* var);

inline std::ostream& operator<<(std::ostream& o, const with<nmod_ring, nmod_poly>& e) {
    nmod_poly_write(o, e.parent(), e.elem_const(), "x");
    return o;
}

inline std::ostream& operator<<(std::ostream& o, const with<nmod_ring, fmpz>& e) {
    o << e.elem_const();
    return o;
}

inline std::ostream& operator<<(std::ostream& o, const nmod_ring_base& R) {
    o << "Integers modulo " << R.modulus;
    return o;
}


void nmod_poly_randtest(nmod_ring& R, nmod_poly& x, rand_state& state, ulimb len);

inline void randtest(nmod_ring& R, nmod_poly& x, rand_state& state, ulimb len) {nmod_poly_randtest(R, x, state, len);}

inline bool nmod_poly_is_zero(const nmod_ring_base& R, const nmod_poly& a) {return a.length() == 0;}
bool nmod_poly_equal(const nmod_ring_base& R, const nmod_poly& a, const nmod_poly& b);

inline bool is_zero(nmod_ring& R, const nmod_poly& a) {return nmod_poly_is_zero(R, a);}
inline bool equal(nmod_ring& R, const nmod_poly& a, const nmod_poly& b) {return nmod_poly_equal(R, a, b);}

// get_set.cpp
void nmod_poly_set(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a);

inline void set(nmod_ring& R, nmod_poly& x, const nmod_poly& a) {nmod_poly_set(R, x, a);}

// add_sub_mul_div.cpp
void nmod_poly_add(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_sub(const nmod_ring_base& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_mul(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_divexact(nmod_ring& R, nmod_poly& q, const nmod_poly& a, const nmod_poly& b);
void nmod_poly_divrem(nmod_ring& R, nmod_poly& q, nmod_poly& r, const nmod_poly& a, const nmod_poly& b);

inline void add(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_add(R, x, a, b);}
inline void sub(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_sub(R, x, a, b);}
inline void mul(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_mul(R, x, a, b);}
inline void divexact(nmod_ring& R, nmod_poly& x, const nmod_poly& a, const nmod_poly& b) {nmod_poly_divexact(R, x, a, b);}
inline void divrem(nmod_ring& R, nmod_poly& q, nmod_poly& r, const nmod_poly& a, const nmod_poly& b) {nmod_poly_divrem(R, q, r, a, b);}


