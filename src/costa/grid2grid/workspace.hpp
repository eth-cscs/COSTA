#pragma once
#include <memory>
#include <omp.h>

namespace costa {

enum CommType {send, recv};

namespace memory {

template <typename T>
struct workspace {
    workspace() = default;

    workspace(int block_dim, int max_threads=omp_get_max_threads())
        : block_dim(block_dim)
        , max_threads(max_threads) {
        transpose_buffer = std::vector<T>(block_dim * max_threads);
    }

    void resize_buffer(CommType type, std::size_t size) {
        switch(type) {
	    case costa::CommType::send: send_buffer.resize(size);
			 break;
	    case costa::CommType::recv: recv_buffer.resize(size);
			 break;
        }
    }

    T* buffer_ptr(CommType type) {
        switch(type) {
	    case costa::CommType::send: return send_buffer.data();
	    case costa::CommType::recv: return recv_buffer.data();
        }
        return nullptr;
    }

    int block_dim = 0;
    int max_threads = 0;
    std::vector<T> transpose_buffer;
    std::vector<T> send_buffer;
    std::vector<T> recv_buffer;
};

template <typename Scalar>
using costa_context = std::unique_ptr<workspace<Scalar>>;

template <typename Scalar>
costa_context<Scalar> make_costa_context() {
    return std::make_unique<workspace<Scalar>>(256);
}

// Meyer's singleton, thread-safe in C++11, but not in C++03.
// The thread-safety is guaranteed by the standard in C++11:
//     If control enters the declaration concurrently
//     while the variable is being initialized,
//     the concurrent execution shall wait
//     for completion of the initialization
template <typename Scalar>
workspace<Scalar>* get_costa_context_instance() {
    static costa_context<Scalar> ctxt = make_costa_context<Scalar>();
    return ctxt.get();
}

}}
