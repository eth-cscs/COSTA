#include <mpi.h>
#include <costa/layout.hpp>
#include <costa/grid2grid/transformer.hpp>
#include <complex>
#include <chrono>
#include <numeric>

auto grid_layout_in(int size, int rank, std::vector<int>& K_rank, int N, std::vector<std::complex<double>>& data_in)
{
    std::vector<int> rowsplit(size + 1);
    rowsplit[0] = 0;
    for (int i = 0; i < size; i++) {
        rowsplit[i + 1] = rowsplit[i] + K_rank[i];
    }
    std::vector<int> colsplit({0, size * N});
    std::vector<int> owners(size);
    for (int i = 0; i < size; i++) {
        owners[i] = i;
    }
    costa::block_t localblock;
    localblock.data = data_in.data();
    localblock.ld = K_rank[rank];
    localblock.row = rank;
    localblock.col = 0;

    return costa::custom_layout<std::complex<double>>(size, 1, rowsplit.data(), colsplit.data(),
            owners.data(), 1, &localblock, 'C');
}

auto grid_layout_out(int size, int rank, int K, int N, std::vector<std::complex<double>>& data_out)
{
    std::vector<int> rowsplit({0, K});
    std::vector<int> colsplit(size + 1);
    colsplit[0] = 0;
    for (int i = 0; i < size; i++) {
        colsplit[i + 1] = colsplit[i] + N;
    }
    std::vector<int> owners(size);
    for (int i = 0; i < size; i++) {
        owners[i] = i;
    }
    costa::block_t localblock;
    localblock.data = data_out.data();
    localblock.ld = K;
    localblock.row = 0;
    localblock.col = rank;

    return costa::custom_layout<std::complex<double>>(1, size, rowsplit.data(), colsplit.data(),
            owners.data(), 1, &localblock, 'C');
}

using time_point_t = std::chrono::high_resolution_clock::time_point;

inline std::chrono::high_resolution_clock::time_point time_now()
{
    return std::chrono::high_resolution_clock::now();
}

inline double time_interval(std::chrono::high_resolution_clock::time_point t0)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(time_now() - t0).count();
}

int main()
{
    int rank, size;

    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> K_rank({113965, 113967}); // long dimension per rank
    int N{200};
    int K = std::accumulate(K_rank.begin(), K_rank.end(), 0);

    std::vector<std::complex<double>> data_in(K_rank[rank] * N * size);
    std::vector<std::complex<double>> data_out(K * N);

    auto layout_in = grid_layout_in(size, rank, K_rank, N, data_in);
    auto layout_out = grid_layout_out(size, rank, K, N, data_out);

    for (int i = 0; i < 10; i++) {
        auto t0 = time_now();
        costa::transform(layout_in, layout_out, 'N', std::complex<double>(1, 0), std::complex<double>(0, 0), MPI_COMM_WORLD);
        if (rank == 0) {
            auto t = time_interval(t0);
            std::cout << "[costa::transform] throughput: "
                  << 16.0 * size * K * N / std::pow(2.0, 30) / t << " Gb/sec" << std::endl;
        }
    }

    MPI_Finalize();

    return 0;
}

