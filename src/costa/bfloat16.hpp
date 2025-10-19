/**
 * @file bfloat16.hpp
 * @brief BFloat16 (Brain Floating Point) type definition
 * @author David Sanftenberg
 * @date 2025-10-19
 * 
 * Implements the BFloat16 format: 16-bit floating point with 1 sign bit,
 * 8 exponent bits, and 7 mantissa bits. This format is compatible with
 * FP32's exponent range but has reduced precision, making it suitable for
 * deep learning and scientific computing where memory bandwidth is critical.
 * 
 * Memory layout (big-endian bit ordering):
 * [15]: Sign bit
 * [14:7]: Exponent (8 bits, same as FP32)
 * [6:0]: Mantissa (7 bits, truncated from FP32's 23 bits)
 */

#pragma once

#include <cstdint>
#include <cstring>
#include <limits>
#include <ostream>

namespace costa {

/**
 * @brief BFloat16 (Brain Floating Point) 16-bit floating point type
 * 
 * This class provides a compact 16-bit floating point representation with
 * the same exponent range as FP32 but reduced mantissa precision. It is
 * designed for use in neural networks and matrix operations where memory
 * bandwidth is more critical than precision.
 * 
 * Key properties:
 * - Size: 16 bits (2 bytes)
 * - Range: Same as FP32 (~1e-38 to ~3e38)
 * - Precision: ~3 decimal digits (vs ~7 for FP32)
 * - Conversion to/from FP32: Simple bit truncation/extension
 */
class bfloat16 {
private:
    uint16_t data_;  ///< Raw 16-bit storage

public:
    /**
     * @brief Default constructor (initializes to zero)
     */
    constexpr bfloat16() : data_(0) {}

    /**
     * @brief Construct from raw uint16_t bits
     * @param raw Raw 16-bit representation
     */
    explicit constexpr bfloat16(uint16_t raw) : data_(raw) {}

    /**
     * @brief Construct from int (convenience for literals like 0, 1)
     * @param value Integer value to convert
     */
    bfloat16(int value) : bfloat16(static_cast<float>(value)) {}

    /**
     * @brief Construct from float (FP32)
     * @param value FP32 value to convert
     * 
     * Conversion truncates the lower 16 bits of the FP32 mantissa.
     * This is a simple right-shift operation on the bit representation.
     * 
     * Note: Non-explicit to allow implicit conversion from float for convenience
     */
    bfloat16(float value) {
        uint32_t fp32_bits;
        std::memcpy(&fp32_bits, &value, sizeof(float));
        
        // BF16 is the upper 16 bits of FP32
        // Simple truncation (round-to-nearest-even would be more accurate but slower)
        data_ = static_cast<uint16_t>(fp32_bits >> 16);
    }

    /**
     * @brief Convert to float (FP32)
     * @return FP32 representation
     * 
     * Conversion extends the BF16 by appending 16 zero bits to the mantissa.
     * This is a simple left-shift operation.
     */
    explicit operator float() const {
        // Extend BF16 to FP32 by shifting left and zero-padding
        uint32_t fp32_bits = static_cast<uint32_t>(data_) << 16;
        
        float result;
        std::memcpy(&result, &fp32_bits, sizeof(float));
        return result;
    }

    /**
     * @brief Get raw 16-bit representation
     * @return Raw bits as uint16_t
     */
    constexpr uint16_t raw() const { return data_; }

    /**
     * @brief Set raw 16-bit representation
     * @param raw Raw bits to set
     */
    constexpr void set_raw(uint16_t raw) { data_ = raw; }

    // Comparison operators
    bool operator==(const bfloat16& other) const {
        return data_ == other.data_;
    }

    bool operator!=(const bfloat16& other) const {
        return data_ != other.data_;
    }

    bool operator<(const bfloat16& other) const {
        return static_cast<float>(*this) < static_cast<float>(other);
    }

    bool operator>(const bfloat16& other) const {
        return static_cast<float>(*this) > static_cast<float>(other);
    }

    bool operator<=(const bfloat16& other) const {
        return static_cast<float>(*this) <= static_cast<float>(other);
    }

    bool operator>=(const bfloat16& other) const {
        return static_cast<float>(*this) >= static_cast<float>(other);
    }

    // Arithmetic operators (implemented via conversion to FP32)
    bfloat16 operator+(const bfloat16& other) const {
        return bfloat16(static_cast<float>(*this) + static_cast<float>(other));
    }

    bfloat16 operator-(const bfloat16& other) const {
        return bfloat16(static_cast<float>(*this) - static_cast<float>(other));
    }

    bfloat16 operator*(const bfloat16& other) const {
        return bfloat16(static_cast<float>(*this) * static_cast<float>(other));
    }

    bfloat16 operator/(const bfloat16& other) const {
        return bfloat16(static_cast<float>(*this) / static_cast<float>(other));
    }

    bfloat16& operator+=(const bfloat16& other) {
        *this = *this + other;
        return *this;
    }

    bfloat16& operator-=(const bfloat16& other) {
        *this = *this - other;
        return *this;
    }

    bfloat16& operator*=(const bfloat16& other) {
        *this = *this * other;
        return *this;
    }

    bfloat16& operator/=(const bfloat16& other) {
        *this = *this / other;
        return *this;
    }

    // Unary operators
    bfloat16 operator-() const {
        return bfloat16(-static_cast<float>(*this));
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const bfloat16& bf) {
        os << static_cast<float>(bf);
        return os;
    }
};

// Mathematical functions for bfloat16
inline float abs(const bfloat16& x) {
    return std::abs(static_cast<float>(x));
}

} // namespace costa

// Specialization of std::numeric_limits for bfloat16
namespace std {
    template<>
    class numeric_limits<costa::bfloat16> {
    public:
        static constexpr bool is_specialized = true;
        static constexpr bool is_signed = true;
        static constexpr bool is_integer = false;
        static constexpr bool is_exact = false;
        static constexpr bool has_infinity = true;
        static constexpr bool has_quiet_NaN = true;
        static constexpr bool has_signaling_NaN = true;
        static constexpr float_denorm_style has_denorm = denorm_present;
        static constexpr bool has_denorm_loss = false;
        static constexpr float_round_style round_style = round_to_nearest;
        static constexpr bool is_iec559 = false;
        static constexpr bool is_bounded = true;
        static constexpr bool is_modulo = false;
        static constexpr int digits = 8;        // Mantissa bits + 1 (implicit)
        static constexpr int digits10 = 2;      // Decimal digits of precision
        static constexpr int max_digits10 = 4;  // Max decimal digits for round-trip
        static constexpr int radix = 2;
        static constexpr int min_exponent = -125;
        static constexpr int min_exponent10 = -37;
        static constexpr int max_exponent = 128;
        static constexpr int max_exponent10 = 38;
        static constexpr bool traps = false;
        static constexpr bool tinyness_before = false;

        static constexpr costa::bfloat16 min() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x0080)); // Smallest normalized positive value
        }

        static constexpr costa::bfloat16 lowest() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0xFF7F)); // Most negative finite value
        }

        static constexpr costa::bfloat16 max() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x7F7F)); // Largest finite value
        }

        static constexpr costa::bfloat16 epsilon() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x3C00)); // 2^-7 (smallest x where 1+x != 1)
        }

        static constexpr costa::bfloat16 round_error() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x3F00)); // 0.5 in BF16
        }

        static constexpr costa::bfloat16 infinity() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x7F80)); // +Infinity
        }

        static constexpr costa::bfloat16 quiet_NaN() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x7FC0)); // Quiet NaN
        }

        static constexpr costa::bfloat16 signaling_NaN() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x7F80 | 1)); // Signaling NaN
        }

        static constexpr costa::bfloat16 denorm_min() noexcept {
            return costa::bfloat16(static_cast<uint16_t>(0x0001)); // Smallest denormalized positive value
        }
    };
} // namespace std
