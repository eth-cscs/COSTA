#include <complex>

#include <costa/grid2grid/memory_buffer.hpp>
#include <costa/grid2grid/communication_data.hpp>

#include <umpire/ResourceManager.hpp>
#include <umpire/strategy/QuickPool.hpp>
#include <umpire/strategy/ThreadSafeAllocator.hpp>

namespace costa {
namespace memory {

umpire::Allocator& get_memory_allocator() {
    auto& rm = umpire::ResourceManager::getInstance();
    static auto pool = rm.makeAllocator<umpire::strategy::QuickPool>
        ("pool", rm.getAllocator("HOST"));

    /*
    static auto thread_safe_pool = rm.makeAllocator
        <umpire::strategy::ThreadSafeAllocator>("thread_safe_pool", pool);
    return thread_safe_pool;
    */
    return pool;
}

/*
template <typename T>
umpire::Allocator& memory_buffer<T>::get_typed_memory_allocator() {
    auto& allocator = get_memory_allocator();
    static umpire::TypedAllocator<T> typed_allocator{allocator};

    return typed_allocator;
}
*/

template <typename T>
memory_buffer<T>::memory_buffer() : memory_buffer(0) {}

template <typename T>
memory_buffer<T>::memory_buffer(std::size_t size): size_(size), ptr_(nullptr), allocated_(true) {
    if (size == 0) return;

    std::size_t buffer_size = static_cast<std::size_t>(size_) * sizeof(T);
    ptr_ = static_cast<T*>(get_memory_allocator().allocate(buffer_size));
}

// template <typename T>
// memory_buffer<T>::memory_buffer(const memory_buffer&) = delete;

// move constructor
template <typename T>
memory_buffer<T>::memory_buffer(memory_buffer&& rhs) : size_(rhs.size_), ptr_(rhs.ptr_), allocated_(rhs.allocated_) {
    rhs.ptr_ = nullptr;
    rhs.size_ = 0;
    rhs.allocated_ = false;
}

// template <typename T>
// memory_buffer<T>& memory_buffer<T>::operator=(const memory_buffer<T>&) = delete;

// move assignment operator
template <typename T>
memory_buffer<T>& memory_buffer<T>::operator=(memory_buffer<T>&& rhs) {
    deallocate();

    size_ = rhs.size_;
    ptr_ = rhs.ptr_;
    allocated_ = rhs.allocated_;

    rhs.size_ = 0;
    rhs.ptr_ = nullptr;
    rhs.allocated_ = false;

    return *this;
}

template <typename T>
memory_buffer<T>::~memory_buffer() {
    deallocate();
}

/// Returns a pointer to the underlying memory at a given index.
///
/// @param index index of the position,
/// @pre @p index < @p size.
template <typename T>
T* memory_buffer<T>::operator()(std::size_t index) {
    assert(index < size_ && index >= 0);
    return ptr_ + index;
}

template <typename T>
const T* memory_buffer<T>::operator()(std::size_t index) const {
    assert(index < size_ && index >= 0);
    return ptr_ + index;
}

/// Returns a pointer to the underlying memory.
/// If @p size == 0 a @c nullptr is returned.
template <typename T>
T* memory_buffer<T>::operator()() {
    return ptr_;
}
template <typename T>
const T* memory_buffer<T>::operator()() const {
    return ptr_;
}

/// Returns the number of elements of type @c T allocated.
template <typename T>
std::size_t memory_buffer<T>::size() const {
    return size_;
}

template <typename T>
void memory_buffer<T>::deallocate() {
    if (allocated_) {
        get_memory_allocator().deallocate(ptr_);
        allocated_ = false;
    }
}
}
}

template class costa::memory::memory_buffer<double>;
template class costa::memory::memory_buffer<float>;
template class costa::memory::memory_buffer<std::complex<float>>;
template class costa::memory::memory_buffer<std::complex<double>>;

template class costa::memory::memory_buffer<costa::message<float>>;
template class costa::memory::memory_buffer<costa::message<double>>;
template class costa::memory::memory_buffer<costa::message<std::complex<float>>>;
template class costa::memory::memory_buffer<costa::message<std::complex<double>>>;

