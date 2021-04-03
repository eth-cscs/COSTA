#include "../gtest.h"
#include <costa/grid2grid/memory_utils.hpp>
#include <costa/grid2grid/threads_workspace.hpp>
#include <vector>
#include <cmath>

TEST(copy2D, row_major) {
    // this is the input
    std::vector<int> in = {
         9,  1,  1, -1, // [0]
         7,  3,  4, -1, // [1]
         5,  5,  1, -1, // [2]
         9,  2,  3, -1, // [3]
         7,  6,  5, -1, // [4]
         2,  2,  4, -1, // [5]
         3,  7,  4, -1, // [6]
         3,  8,  1, -1  // [7]
    };

    // this is the correct result we are expecting
    std::vector<int> result = {
         9,  1,  1, -1, -1, // [0]
         7,  3,  4, -1, -1, // [1]
         5,  5,  1, -1, -1, // [2]
         9,  2,  3, -1, -1, // [3]
         7,  6,  5, -1, -1, // [4]
         2,  2,  4, -1, -1, // [5]
         3,  7,  4, -1, -1, // [6]
         3,  8,  1, -1, -1  // [7]
    };

    std::vector<int> out(result.size());

    costa::memory::threads_workspace<int> workspace(256);

    int n_rows = 8; int n_cols = 3;
    int in_stride = 4;
    int out_stride = 5;

    costa::memory::copy_and_transform(n_rows, n_cols,
                                      in.data(), 
                                      in_stride, // in-stride
                                      false, // => row-major
                                      out.data(),
                                      out_stride, // out-stride
                                      false, // => row-major
                                      false, // => no transpose
                                      false, // => no conjugate
                                      1, 0, // alpha, beta
                                      workspace
                                      );

    // ignore stride when checking the correctness
    std::cout << "Result should be: " << std::endl;
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            std::cout << result[i * out_stride + j] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "Output is: " << std::endl;
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            std::cout << out[i * out_stride + j] << ", ";
        }
        std::cout << std::endl;
    }
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            EXPECT_TRUE(result[i * out_stride + j] == out[i * out_stride + j]);
        }
    }
}

TEST(copy2D, col_major) {
    // this is the input
    std::vector<int> in = {
         9,  1,  1, -1, // [0]
         7,  3,  4, -1, // [1]
         5,  5,  1, -1, // [2]
         9,  2,  3, -1, // [3]
         7,  6,  5, -1, // [4]
         2,  2,  4, -1, // [5]
         3,  7,  4, -1, // [6]
         3,  8,  1, -1  // [7]
    };

    // this is the correct result we are expecting
    std::vector<int> result = {
         9,  1,  1, -1, -1, // [0]
         7,  3,  4, -1, -1, // [1]
         5,  5,  1, -1, -1, // [2]
         9,  2,  3, -1, -1, // [3]
         7,  6,  5, -1, -1, // [4]
         2,  2,  4, -1, -1, // [5]
         3,  7,  4, -1, -1, // [6]
         3,  8,  1, -1, -1  // [7]
    };

    std::vector<int> out(result.size());

    costa::memory::threads_workspace<int> workspace(256);

    int n_rows = 3; int n_cols = 8;
    int in_stride = 4;
    int out_stride = 5;

    costa::memory::copy_and_transform(n_rows, n_cols,
                                      in.data(), 
                                      in_stride, // in-stride
                                      true, // => col-major
                                      out.data(),
                                      out_stride, // out-stride
                                      true, // => col-major
                                      false, // => no transpose
                                      false, // => no conjugate
                                      1, 0, // alpha, beta
                                      workspace
                                      );

    // ignore stride when checking the correctness
    std::cout << "Result should be: " << std::endl;
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            std::cout << result[j * out_stride + i] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "Output is: " << std::endl;
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            std::cout << out[j * out_stride + i] << ", ";
        }
        std::cout << std::endl;
    }
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            EXPECT_TRUE(result[j * out_stride + i] == out[j * out_stride + i]);
        }
    }
}

TEST(transpose, row_to_col_major) {
    // this is the input
    std::vector<int> in = {
         9,  1,  1, -1, // [0]
         7,  3,  4, -1, // [1]
         5,  5,  1, -1, // [2]
         9,  2,  3, -1, // [3]
         7,  6,  5, -1, // [4]
         2,  2,  4, -1, // [5]
         3,  7,  4, -1, // [6]
         3,  8,  1, -1  // [7]
    };

    // this is the correct result we are expecting
    std::vector<int> result = {
        9, 7, 5, 9, 7, 2, 3, 3, -1, -1,
        1, 3, 5, 2, 6, 2, 7, 8, -1, -1,
        1, 4, 1, 3, 5, 4, 4, 1, -1, -1
    };

    std::vector<int> out(result.size());

    costa::memory::threads_workspace<int> workspace(256);

    int n_rows = 8; int n_cols = 3;
    int in_stride = 4;
    int out_stride = 10;

    // row -> col major ordering = transpose
    costa::memory::copy_and_transform(n_rows, n_cols,
                                      in.data(), 
                                      in_stride, // in-stride
                                      false, // => row-major
                                      out.data(),
                                      out_stride, // out-stride
                                      true, // => col-major
                                      false, // => no transpose
                                      false, // => no conjugate
                                      1, 0, // alpha, beta
                                      workspace
                                      );

    // ignore stride when checking the correctness
    std::cout << "Result should be: " << std::endl;
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            std::cout << result[j * out_stride + i] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "Output is: " << std::endl;
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            std::cout << out[j * out_stride + i] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "Finished the output" << std::endl;
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            EXPECT_TRUE(result[j * out_stride + i] == out[j * out_stride + i]);
        }
    }
}

