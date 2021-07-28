#pragma once

#include <memory>
#include <cstdlib>

namespace costa {
namespace memory {

template <typename T>
class memory_buffer {
public:
    memory_buffer();

    memory_buffer(std::size_t size);

    memory_buffer(const memory_buffer&) = delete;

    // move constructor
    memory_buffer(memory_buffer&& rhs);

    memory_buffer& operator=(const memory_buffer&) = delete;

    // move assignment operator
    memory_buffer& operator=(memory_buffer&& rhs);

    ~memory_buffer();

    /// Returns a pointer to the underlying memory at a given index.
    ///
    /// @param index index of the position,
    /// @pre @p index < @p size.
    T* operator()(std::size_t index);

    const T* operator()(std::size_t index) const;

    /// Returns a pointer to the underlying memory.
    /// If @p size == 0 a @c nullptr is returned.
    T* operator()();
    const T* operator()() const;

    /// Returns the number of elements of type @c T allocated.
    std::size_t size() const;

    // umpire::Allocator& get_memory_allocator();
    // umpire::Allocator& get_typed_memory_allocator();
private:
    void deallocate();

    std::size_t size_;
    T* ptr_;
    bool allocated_;
};
}
}
