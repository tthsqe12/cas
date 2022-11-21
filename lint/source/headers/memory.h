#pragma once
#include "buffers.h"

struct ulimb_realizer {
    virtual ulimb* get_limbs() = 0;
};


struct leaf_allocator {
    void* buffer;
    size_t alloc;

    leaf_allocator() : buffer(nullptr), alloc(0) {}

    leaf_allocator(leaf_allocator&& other) = delete;
    leaf_allocator(const leaf_allocator& other) = delete;
    leaf_allocator& operator=(const leaf_allocator& other) = delete;

    template <typename T>
    T* get_space(size_t n) {
        n *= sizeof(T);
        if (UNLIKELY(n > alloc)) {
            std::free(buffer);
            n = round_up(n, 4096);
            buffer = my_aligned_alloc(4096, n);
            alloc = n;
        }
        return reinterpret_cast<T*>(buffer);
    }
};

struct recursive_allocator {
    size_t length;
    std::vector<std::pair<void*, size_t>> stack;

    recursive_allocator() : length(0) {}

    void* alloc(size_t n);

    size_t in_use()
    {
        return length;
    }

    void restore_in_use(size_t n)
    {
        assert(n <= length);
        length = n;
    }

    ~recursive_allocator()
    {
        for (auto& i : stack)
            std::free(std::get<0>(i));
    }
};

struct tmp_allocator_internal {
    leaf_allocator leaf;
    recursive_allocator rec;

    recursive_allocator& rec_buffer() {return rec;}
};

extern tmp_allocator_internal tls_tmp_allocator;

template<typename T>
inline T* leaf_alloc(size_t n) {
    return tls_tmp_allocator.leaf.get_space<T>(n);
}

class tmp_allocator
{
    size_t _length;
    tmp_allocator_internal* _a;

public:
    tmp_allocator() {
        _a = &tls_tmp_allocator;
        _length = _a->rec.in_use();
    }

    ~tmp_allocator() {
        _a->rec.restore_in_use(_length);
    }

    template <typename T>
    inline T* recursive_alloc(size_t n) {
        return reinterpret_cast<T*>(_a->rec.alloc(n*sizeof(T)));
    }

    template <typename T>
    inline T* recursive_alloc(size_t n, size_t m) {
        return reinterpret_cast<T*>(_a->rec.alloc(n*m*sizeof(T)));
    }
};

