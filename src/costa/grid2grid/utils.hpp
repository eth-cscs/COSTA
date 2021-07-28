#pragma once
#include <costa/grid2grid/profiler.hpp>
#include <costa/grid2grid/grid_layout.hpp>
#include <costa/grid2grid/grid_cover.hpp>
#include <costa/grid2grid/communication_data.hpp>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <queue>
#include <tuple>

namespace costa {

std::vector<std::vector<int>> topology_cost(MPI_Comm comm);

namespace utils {

bool if_should_transpose(const char src_ordering,
                         const char dest_ordering,
                         const char trans);

std::unordered_map<int, int> rank_to_comm_vol_for_block(
        const assigned_grid2D& g_init,
        const block_coordinates &b_coord,
        grid_cover &g_cover,
        const assigned_grid2D& g_final);

template <typename T>
communication_data<T> prepare_to_send(
                                      costa::memory::memory_buffer<message<T>>& messages_buffer,
                                      grid_layout<T> &init_layout,
                                      grid_layout<T> &final_layout,
                                      int rank,
                                      const T alpha, const T beta,
                                      bool transpose, bool conjugate);

template <typename T>
communication_data<T> prepare_to_send(
                                      costa::memory::memory_buffer<message<T>>& messages_buffer,
                                      std::vector<layout_ref<T>>& from,
                                      std::vector<layout_ref<T>>& to,
                                      int rank,
                                      const T* alpha, const T* beta,
                                      bool* transpose,
                                      bool* conjugate);
template <typename T> 
communication_data<T> prepare_to_recv(
                                      costa::memory::memory_buffer<message<T>>& messages_buffer,
                                      grid_layout<T> &final_layout,
                                      grid_layout<T> &init_layout,
                                      int rank,
                                      const T alpha, const T beta,
                                      const bool transpose, const bool conjugate);

template <typename T>
communication_data<T> prepare_to_recv(costa::memory::memory_buffer<message<T>>& messages_buffer,
                                      std::vector<layout_ref<T>>& to,
                                      std::vector<layout_ref<T>>& from,
                                      int rank,
                                      const T* alpha, const T* beta,
                                      bool* transpose,
                                      bool* conjugate);
} // namespace utils
} // namespace costa

