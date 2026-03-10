/**
 * @file test_bfloat16.cpp
 * @brief BFloat16 type support tests for COSTA
 * @author David Sanftenberg
 * @date 2025-10-19
 *
 * Tests BF16 type integration in COSTA's core infrastructure:
 * - MPI type wrapper (MPI_UINT16_T mapping)
 * - Template instantiations for block, local_blocks, message
 * - Template instantiations for transform operations
 * - ADL support for abs() and conjugate_f()
 */

#include "../gtest.h"
#include <costa/bfloat16.hpp>
#include <costa/grid2grid/block.hpp>
#include <costa/grid2grid/mpi_type_wrapper.hpp>
#include <costa/grid2grid/transform.hpp>
#include <costa/layout.hpp>
#include <vector>

using namespace costa;

/**
 * @brief Test BF16 type basic properties
 *
 * Validates that bfloat16 type has correct size (2 bytes) and
 * can be used with COSTA's type system.
 */
TEST(BFloat16COSTA, TypeProperties) {
    EXPECT_EQ(sizeof(bfloat16), 2) << "BF16 should be 2 bytes";

    // Test that we can create vectors of BF16
    std::vector<bfloat16> vec(10);
    EXPECT_EQ(vec.size(), 10);

    // Test basic assignment
    vec[0] = bfloat16(3.14f);
    float result = static_cast<float>(vec[0]);
    EXPECT_NEAR(result, 3.14f, 0.02f) << "BF16 conversion should preserve value (within BF16 precision)";
}

/**
 * @brief Test MPI type wrapper for BF16
 *
 * Validates that BF16 uses MPI_UINT16_T for MPI communication
 * as defined in mpi_type_wrapper.hpp.
 */
TEST(BFloat16COSTA, MPITypeWrapper) {
    auto mpi_type = mpi_type_wrapper<bfloat16>::type();
    EXPECT_EQ(mpi_type, MPI_UINT16_T) << "BF16 should use MPI_UINT16_T";
}

/**
 * @brief Test conjugate_f function for BF16
 *
 * Validates the conjugate_f template instantiation added in block.cpp.
 * For real types like BF16, conjugate should be identity.
 */
TEST(BFloat16COSTA, ConjugateFunction) {
    bfloat16 val(3.14f);
    bfloat16 conj = conjugate_f(val);
    
    EXPECT_FLOAT_EQ(static_cast<float>(conj), static_cast<float>(val))
        << "Conjugate of real BF16 value should be itself";
}

/**
 * @brief Test abs() function via ADL for BF16
 *
 * Validates that abs() can be found via ADL for bfloat16 type,
 * using the abs() function defined in costa namespace.
 */
TEST(BFloat16COSTA, AbsFunction) {
    bfloat16 positive(3.14f);
    bfloat16 negative(-3.14f);
    
    // ADL should find costa::abs()
    bfloat16 abs_pos = abs(positive);
    bfloat16 abs_neg = abs(negative);
    
    EXPECT_GT(static_cast<float>(abs_pos), 0.0f) << "abs(positive) should be positive";
    EXPECT_GT(static_cast<float>(abs_neg), 0.0f) << "abs(negative) should be positive";
    EXPECT_NEAR(static_cast<float>(abs_pos), 3.14f, 0.02f);
    EXPECT_NEAR(static_cast<float>(abs_neg), 3.14f, 0.02f);
}

/**
 * @brief Test block<bfloat16> template instantiation
 *
 * Validates that block template can be instantiated with bfloat16
 * as added in block.cpp template instantiations.
 */
TEST(BFloat16COSTA, BlockInstantiation) {
    // Create a simple interval and block
    interval row_int(0, 10);
    interval col_int(0, 10);
    block_coordinates coords(0, 0);
    
    std::vector<bfloat16> data(100, bfloat16(1.0f));
    
    // This tests that block<bfloat16> template is properly instantiated
    block<bfloat16> blk(row_int, col_int, coords, data.data(), 0);
    
    EXPECT_EQ(blk.total_size(), 100) << "Block should have 100 elements";
}

/**
 * @brief Test local_blocks<bfloat16> template instantiation
 *
 * Validates that local_blocks template can be instantiated with bfloat16.
 */
TEST(BFloat16COSTA, LocalBlocksInstantiation) {
    interval row_int(0, 5);
    interval col_int(0, 5);
    block_coordinates coords(0, 0);
    
    std::vector<bfloat16> data(25, bfloat16(2.0f));
    
    std::vector<block<bfloat16>> blocks;
    blocks.emplace_back(row_int, col_int, coords, data.data(), 0);
    
    // This tests that local_blocks<bfloat16> template is properly instantiated
    local_blocks<bfloat16> local_blks(std::move(blocks));
    
    EXPECT_EQ(local_blks.num_blocks(), 1) << "Should have 1 block";
}

/**
 * @brief Test that block<bfloat16> can transpose
 *
 * Validates that the transpose() method works for bfloat16 blocks.
 */
TEST(BFloat16COSTA, BlockTranspose) {
    interval row_int(0, 4);
    interval col_int(0, 3);
    block_coordinates coords(0, 0);
    
    // 4x3 matrix in row-major: [[1,2,3], [4,5,6], [7,8,9], [10,11,12]]
    std::vector<bfloat16> data = {
        bfloat16(1.0f), bfloat16(2.0f), bfloat16(3.0f),
        bfloat16(4.0f), bfloat16(5.0f), bfloat16(6.0f),
        bfloat16(7.0f), bfloat16(8.0f), bfloat16(9.0f),
        bfloat16(10.0f), bfloat16(11.0f), bfloat16(12.0f)
    };
    
    block<bfloat16> blk(row_int, col_int, coords, data.data(), 0);
    
    // Test that transpose() is properly instantiated
    EXPECT_NO_THROW({
        blk.transpose();
    }) << "block<bfloat16>::transpose() should be properly instantiated";
}

/**
 * @brief Test local_blocks<bfloat16> transpose
 *
 * Validates that local_blocks::transpose() works for bfloat16.
 */
TEST(BFloat16COSTA, LocalBlocksTranspose) {
    interval row_int(0, 3);
    interval col_int(0, 3);
    block_coordinates coords(0, 0);
    
    std::vector<bfloat16> data(9, bfloat16(1.5f));
    
    std::vector<block<bfloat16>> blocks;
    blocks.emplace_back(row_int, col_int, coords, data.data(), 0);
    
    local_blocks<bfloat16> local_blks(std::move(blocks));
    
    // Test that transpose() is properly instantiated
    EXPECT_NO_THROW({
        local_blks.transpose();
    }) << "local_blocks<bfloat16>::transpose() should be properly instantiated";
}
