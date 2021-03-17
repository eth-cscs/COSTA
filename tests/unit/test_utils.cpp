#include "../gtest.h"
#include <vector>
#include <conflux/lu/conflux_opt.hpp>
#include <conflux/lu/utils.hpp>
#include <conflux/lu/matrix_view.hpp>
#include <cmath>

TEST(push_pivots_up, row_major) {
    // this is the input
    std::vector<int> in = {
        9, 1, 1, 9, 5, 5, 3, 3, // [0]
        7, 3, 4, 5, 2, 4, 5, 2, // [1]
        5, 5, 1, 3, 3, 9, 1, 9, // [2]
        9, 2, 3, 9, 5, 2, 2, 9, // [3]
        7, 6, 5, 7, 8, 1, 4, 4, // [4]
        2, 2, 4, 6, 5, 2, 6, 5, // [5]
        3, 7, 4, 4, 4, 7, 1, 7, // [6]
        3, 8, 1, 6, 4, 8, 7, 8  // [7]
    };

    // this is the correct result we are expecting
    std::vector<int> result = {
        5, 5, 1, 3, 3, 9, 1, 9, // [2]
        7, 3, 4, 5, 2, 4, 5, 2, // [1]
        2, 2, 4, 6, 5, 2, 6, 5, // [5]
        9, 2, 3, 9, 5, 2, 2, 9, // [3]
        7, 6, 5, 7, 8, 1, 4, 4, // [4]
        9, 1, 1, 9, 5, 5, 3, 3, // [0]
        3, 7, 4, 4, 4, 7, 1, 7, // [6]
        3, 8, 1, 6, 4, 8, 7, 8  // [7]
    };


    // dimension
    int n = (int) sqrt(in.size());

    // the starting row for pushing pivots
    int first_non_pivot_row = 0;

    // matrix layout (row/col major)
    auto layout = conflux::order::row_major;


    // current (local) pivots
    // curPivots[0] is the number of local pivots
    // (i.e. the size of the array)
    std::vector<int> curPivots = {3, 2, 1, 5};

    // temporary matrix
    std::vector<int> in_temp(result.size());

    // input matrix
    conflux::matrix_view in_mat(in.data(), n, n, layout);
    conflux::matrix_view result_mat(result.data(), n, n, layout);

    std::cout << "Input matrix:" << std::endl;
    std::cout << in_mat.to_string() << std::endl;

    // in-place pushing pivots up
    conflux::push_pivots_up<int>(in, in_temp, 
                        n, n, layout, 
                        curPivots, first_non_pivot_row);

    std::vector<bool> correct_row(n, true);

    std::cout << "Output matrix:" << std::endl;
    std::cout << in_mat.to_string() << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "Target Matrix:" << std::endl;
    std::cout << result_mat.to_string() << std::endl;
    std::cout << "===================" << std::endl;

    // the output matrix and the result matrix should match
    EXPECT_TRUE(in_mat == result_mat);
}

TEST(permute_rows, row_major) {
    // input matrix
    std::vector<std::vector<double>> in_mat = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            1, 2, 3,
            4, 5, 6,
            7, 8, 9
        },
        // rectangular case
        {
            1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12,
        },
        // rectangular case - strided
        {
            1, 2, 3, -1, -1,
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            10, 11, 12, -1, -1
        },
        // rectangular case - strided
        // (with col-offset)
        {
            1, 2, 3, -1, -1,
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            10, 11, 12, -1, -1
        },
        // rectangular case - strided
        // (with col-offset and different target strides)
        {
            1, 2, 3, -1, -1,
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            10, 11, 12, -1, -1
        },
    };
    // this is the target (correct) result
    std::vector<std::vector<double>> result = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            4, 5, 6,
            7, 8, 9,
            1, 2, 3
        },
        // rectangular case
        {
            4, 5, 6,
            7, 8, 9,
            1, 2, 3,
            10, 11, 12
        },
        // rectangular case - strided
        {
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            1, 2, 3, -1, -1,
            10, 11, 12, -1, -1
        },
        // rectangular case - strided
        // (with col-offset)
        {
            5, 6, -1, -1, -1,
            8, 9, -1, -1, -1,
            2, 3, -1, -1, -1,
            11, 12, -1, -1, -1
        },
        // rectangular case - strided
        // (with col-offset and different target stride)
        {
            5, 6, -1, -1, -1, -1, -1, -1,
            8, 9, -1, -1, -1, -1, -1, -1,
            2, 3, -1, -1, -1, -1, -1, -1,
            11, 12, -1, -1, -1, -1, -1, -1
        }
    };

    // permutation: i-th row -> perm[i]-th row
    // 0th row -> 2nd row
    // 1st row -> 0th row
    // 2nd row -> 1st row
    // 3rd row -> 3rd row
    std::vector<std::vector<int>> perm = {
        {0},
        {0},
        {2, 0, 1},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3}
    };

    std::vector<int> col_offset = {0, 0, 0, 0, 0, 1, 1};
    std::vector<int> n_rows = {1, 1, 3, 4, 4, 4, 4};
    std::vector<int> n_cols = {1, 1, 3, 3, 3, 3, 3};
    std::vector<int> strides = {1, 2, 3, 3, 5, 5, 5};
    std::vector<int> strides_offset = {0, 0, 0, 0, 0, 0, 3};

    // output matrix
    std::vector<std::vector<double>> out_mat(result.size());

    // row major data layout
    auto layout = conflux::order::row_major;

    for (int i = 0; i < in_mat.size(); ++i) {
        // flush the output
        out_mat[i] = std::vector<double>(result[i].size());

        // input matrix view
        conflux::matrix_view<double> in(in_mat[i].data(),
                                        n_rows[i], n_cols[i],
                                        strides[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == in.layout());
        }

        // output matrix view
        conflux::matrix_view<double> out(out_mat[i].data(),
                                        n_rows[i], n_cols[i]-col_offset[i],
                                        strides[i] + strides_offset[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == out.layout());
        }

        // permute rows
        conflux::permute_rows<double>(in, out, perm[i]);

        // target matrix view (used only for output)
        conflux::matrix_view<double> res(result[i].data(),
                                         n_rows[i], 
                                         n_cols[i]-col_offset[i], 
                                         strides[i] + strides_offset[i],
                                         layout);

        std::cout << "Input matrix:" << std::endl;
        std::cout << in.to_string() << std::endl;
        std::cout << "Output matrix:" << std::endl;
        std::cout << out.to_string() << std::endl;
        std::cout << "-------------------" << std::endl;
        std::cout << "Target Matrix:" << std::endl;
        std::cout << res.to_string() << std::endl;
        std::cout << "===================" << std::endl;

        EXPECT_TRUE(out == res);
    }
}

TEST(permute_rows, col_major) {
    // input matrix
    std::vector<std::vector<double>> in_mat = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            1, 4, 7,
            2, 5, 8,
            3, 6, 9
        },
        // rectangular case
        {
            1, 5, 9,
            2, 6, 10,
            3, 7, 11,
            4, 8, 12
        },
        // rectangular case - strided
        {
            1, 5, 9, -1, -1,
            2, 6, 10, -1, -1,
            3, 7, 11, -1, -1,
            4, 8, 12, -1, -1
        },
        // rectangular case - strided (with col-offset)
        {
            1, 5, 9, -1, -1,
            2, 6, 10, -1, -1,
            3, 7, 11, -1, -1,
            4, 8, 12, -1, -1
        },
        // rectangular case - strided 
        // (with col-offset, and different target stride)
        {
            1, 5, 9, -1, -1,
            2, 6, 10, -1, -1,
            3, 7, 11, -1, -1,
            4, 8, 12, -1, -1
        }
    };
    // this is the target (correct) result
    std::vector<std::vector<double>> result = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            4, 7, 1,
            5, 8, 2,
            6, 9, 3
        },
        // rectangular case
        {
            5, 9,  1,
            6, 10, 2,
            7, 11, 3,
            8, 12, 4
        },
        // rectangular case - strided
        {
            5, 9, 1, -1, -1,
            6, 10, 2, -1, -1,
            7, 11, 3, -1, -1,
            8, 12, 4, -1, -1
        },
        // rectangular case - strided (with col-offset)
        {
            6, 10, 2, -1, -1,
            7, 11, 3, -1, -1,
            8, 12, 4, -1, -1
        },
        // rectangular case - strided 
        // (with col-offset, and different target stride)
        {
            6, 10, 2, -1, -1, -1, -1, -1,
            7, 11, 3, -1, -1, -1, -1, -1,
            8, 12, 4, -1, -1, -1, -1, -1
        }
    };

    // permutation: i-th row -> perm[i]-th row
    // 0th row -> 2nd row
    // 1st row -> 0th row
    // 2nd row -> 1st row
    // 3rd row -> 3rd row
    std::vector<std::vector<int>> perm = {
        {0},
        {0},
        {2, 0, 1},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3}
    };

    std::vector<int> col_offset = {0, 0, 0, 0, 0, 1, 1};
    std::vector<int> n_rows = {1, 1, 3, 3, 3, 3, 3};
    std::vector<int> n_cols = {1, 1, 3, 4, 4, 4, 4};
    std::vector<int> strides = {1, 2, 3, 3, 5, 5, 5};
    std::vector<int> strides_offset = {0, 0, 0, 0, 0, 0, 3};

    // output matrix
    std::vector<std::vector<double>> out_mat(result.size());

    // row major data layout
    auto layout = conflux::order::col_major;

    for (int i = 0; i < in_mat.size(); ++i) {
        // flush the output
        out_mat[i] = std::vector<double>(result[i].size());

        // input matrix view
        conflux::matrix_view<double> in(in_mat[i].data(),
                                        n_rows[i], n_cols[i],
                                        strides[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == in.layout());
        }

        // output matrix view
        conflux::matrix_view<double> out(out_mat[i].data(),
                                        n_rows[i], n_cols[i]-col_offset[i],
                                        strides[i] + strides_offset[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == out.layout());
        }

        // permute rows
        conflux::permute_rows<double>(in, out, perm[i]);

        // target matrix view (used only for output)
        conflux::matrix_view<double> res(result[i].data(),
                                         n_rows[i], 
                                         n_cols[i]-col_offset[i], 
                                         strides[i] + strides_offset[i],
                                         layout);

        std::cout << "Input matrix:" << std::endl;
        std::cout << in.to_string() << std::endl;
        std::cout << "Output matrix:" << std::endl;
        std::cout << out.to_string() << std::endl;
        std::cout << "-------------------" << std::endl;
        std::cout << "Target Matrix:" << std::endl;
        std::cout << res.to_string() << std::endl;
        std::cout << "===================" << std::endl;

        EXPECT_TRUE(out == res);
    }
}

TEST(inverse_permute_rows, row_major) {
    // input matrix
    std::vector<std::vector<double>> in_mat = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            1, 2, 3,
            4, 5, 6,
            7, 8, 9
        },
        // rectangular case
        {
            1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12,
        },
        // rectangular case - strided
        {
            1, 2, 3, -1, -1,
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            10, 11, 12, -1, -1
        },
        // rectangular case - strided
        // (with col-offset)
        {
            1, 2, 3, -1, -1,
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            10, 11, 12, -1, -1
        },
        // rectangular case - strided
        // (with col-offset and different target strides)
        {
            1, 2, 3, -1, -1,
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            10, 11, 12, -1, -1
        },
    };
    // this is the target (correct) result
    std::vector<std::vector<double>> result = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            4, 5, 6,
            7, 8, 9,
            1, 2, 3
        },
        // rectangular case
        {
            4, 5, 6,
            7, 8, 9,
            1, 2, 3,
            10, 11, 12
        },
        // rectangular case - strided
        {
            4, 5, 6, -1, -1,
            7, 8, 9, -1, -1,
            1, 2, 3, -1, -1,
            10, 11, 12, -1, -1
        },
        // rectangular case - strided
        // (with col-offset)
        {
            5, 6, -1, -1, -1,
            8, 9, -1, -1, -1,
            2, 3, -1, -1, -1,
            11, 12, -1, -1, -1
        },
        // rectangular case - strided
        // (with col-offset and different target stride)
        {
            5, 6, -1, -1, -1, -1, -1, -1,
            8, 9, -1, -1, -1, -1, -1, -1,
            2, 3, -1, -1, -1, -1, -1, -1,
            11, 12, -1, -1, -1, -1, -1, -1
        }
    };

    // permutation: i-th row -> perm[i]-th row
    // 0th row -> 2nd row
    // 1st row -> 0th row
    // 2nd row -> 1st row
    // 3rd row -> 3rd row
    std::vector<std::vector<int>> perm = {
        {0},
        {0},
        {2, 0, 1},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3}
    };

    std::vector<int> col_offset = {0, 0, 0, 0, 0, 1, 1};
    std::vector<int> n_rows = {1, 1, 3, 4, 4, 4, 4};
    std::vector<int> n_cols = {1, 1, 3, 3, 3, 3, 3};
    std::vector<int> strides = {1, 2, 3, 3, 5, 5, 5};
    std::vector<int> strides_offset = {0, 0, 0, 0, 0, 0, 3};

    // output matrix
    std::vector<std::vector<double>> out_mat(result.size());
    std::vector<std::vector<double>> in2_mat(in_mat.size());

    // row major data layout
    auto layout = conflux::order::row_major;

    for (int i = 0; i < in_mat.size(); ++i) {
        // flush the output
        out_mat[i] = std::vector<double>(result[i].size());
        in2_mat[i] = std::vector<double>(in_mat[i].size());

        // input matrix view
        conflux::matrix_view<double> in(in_mat[i].data(),
                                        n_rows[i], n_cols[i],
                                        strides[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == in.layout());
        }

        // output matrix view
        conflux::matrix_view<double> out(out_mat[i].data(),
                                        n_rows[i], n_cols[i]-col_offset[i],
                                        strides[i] + strides_offset[i],
                                        layout);

        // input matrix view
        conflux::matrix_view<double> in2(in2_mat[i].data(),
                                        n_rows[i], n_cols[i]-col_offset[i],
                                        strides[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == out.layout());
        }

        // permute rows
        conflux::permute_rows<double>(in, out, perm[i]);
        conflux::inverse_permute_rows<double>(out, in2, perm[i]);

        std::cout << "Input matrix:" << std::endl;
        std::cout << in.to_string() << std::endl;
        std::cout << "Output matrix:" << std::endl;
        std::cout << in2.to_string() << std::endl;
        std::cout << "-------------------" << std::endl;
        std::cout << "Target Matrix:" << std::endl;
        std::cout << in.to_string() << std::endl;
        std::cout << "===================" << std::endl;

        EXPECT_TRUE(in.submatrix(0, col_offset[i], n_rows[i], n_cols[i]) == in2);
    }
}

TEST(inverse_permute_rows, col_major) {
    // input matrix
    std::vector<std::vector<double>> in_mat = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            1, 4, 7,
            2, 5, 8,
            3, 6, 9
        },
        // rectangular case
        {
            1, 5, 9,
            2, 6, 10,
            3, 7, 11,
            4, 8, 12
        },
        // rectangular case - strided
        {
            1, 5, 9, -1, -1,
            2, 6, 10, -1, -1,
            3, 7, 11, -1, -1,
            4, 8, 12, -1, -1
        },
        // rectangular case - strided (with col-offset)
        {
            1, 5, 9, -1, -1,
            2, 6, 10, -1, -1,
            3, 7, 11, -1, -1,
            4, 8, 12, -1, -1
        },
        // rectangular case - strided 
        // (with col-offset, and different target stride)
        {
            1, 5, 9, -1, -1,
            2, 6, 10, -1, -1,
            3, 7, 11, -1, -1,
            4, 8, 12, -1, -1
        }
    };
    // this is the target (correct) result
    std::vector<std::vector<double>> result = {
        // weird cases
        {
            1
        },
        // weird case - strided
        {
            1, -1
        },
        // square case
        {
            4, 7, 1,
            5, 8, 2,
            6, 9, 3
        },
        // rectangular case
        {
            5, 9,  1,
            6, 10, 2,
            7, 11, 3,
            8, 12, 4
        },
        // rectangular case - strided
        {
            5, 9, 1, -1, -1,
            6, 10, 2, -1, -1,
            7, 11, 3, -1, -1,
            8, 12, 4, -1, -1
        },
        // rectangular case - strided (with col-offset)
        {
            6, 10, 2, -1, -1,
            7, 11, 3, -1, -1,
            8, 12, 4, -1, -1
        },
        // rectangular case - strided 
        // (with col-offset, and different target stride)
        {
            6, 10, 2, -1, -1, -1, -1, -1,
            7, 11, 3, -1, -1, -1, -1, -1,
            8, 12, 4, -1, -1, -1, -1, -1
        }
    };

    // permutation: i-th row -> perm[i]-th row
    // 0th row -> 2nd row
    // 1st row -> 0th row
    // 2nd row -> 1st row
    // 3rd row -> 3rd row
    std::vector<std::vector<int>> perm = {
        {0},
        {0},
        {2, 0, 1},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3},
        {2, 0, 1, 3}
    };

    std::vector<int> col_offset = {0, 0, 0, 0, 0, 1, 1};
    std::vector<int> n_rows = {1, 1, 3, 3, 3, 3, 3};
    std::vector<int> n_cols = {1, 1, 3, 4, 4, 4, 4};
    std::vector<int> strides = {1, 2, 3, 3, 5, 5, 5};
    std::vector<int> strides_offset = {0, 0, 0, 0, 0, 0, 3};

    // output matrix
    std::vector<std::vector<double>> out_mat(result.size());
    std::vector<std::vector<double>> in2_mat(in_mat.size());

    // row major data layout
    auto layout = conflux::order::col_major;

    for (int i = 0; i < in_mat.size(); ++i) {
        // flush the output
        out_mat[i] = std::vector<double>(result[i].size());
        in2_mat[i] = std::vector<double>(in_mat[i].size());

        // input matrix view
        conflux::matrix_view<double> in(in_mat[i].data(),
                                        n_rows[i], n_cols[i],
                                        strides[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == in.layout());
        }

        // output matrix view
        conflux::matrix_view<double> out(out_mat[i].data(),
                                        n_rows[i], n_cols[i]-col_offset[i],
                                        strides[i] + strides_offset[i],
                                        layout);

        // input matrix view
        conflux::matrix_view<double> in2(in2_mat[i].data(),
                                        n_rows[i], n_cols[i]-col_offset[i],
                                        strides[i],
                                        layout);

        // layout makes sense if n_rows > 1 and n_cols > 1
        if (n_rows[i] > 1 && n_cols[i] > 1) {
            EXPECT_TRUE(layout == out.layout());
        }

        // permute rows
        conflux::permute_rows<double>(in, out, perm[i]);
        conflux::inverse_permute_rows<double>(out, in2, perm[i]);

        std::cout << "Input matrix:" << std::endl;
        std::cout << in.to_string() << std::endl;
        std::cout << "Output matrix:" << std::endl;
        std::cout << in2.to_string() << std::endl;
        std::cout << "-------------------" << std::endl;
        std::cout << "Target Matrix:" << std::endl;
        std::cout << in.to_string() << std::endl;
        std::cout << "===================" << std::endl;

        EXPECT_TRUE(in.submatrix(0, col_offset[i], n_rows[i], n_cols[i]) == in2);
    }
}
