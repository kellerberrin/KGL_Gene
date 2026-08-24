//
// Created by kellerberrin on 10/8/26.
//

#ifndef KEL_STRING_CRC64HASH_H
#define KEL_STRING_CRC64HASH_H

#include <array>
#include <string_view>
#include <cstdint>
#include <cstddef>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implement the 64 bit CRC-64-ECMA-182 hash function.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kellerberrin { //  organization level namespace

// Prefer namespace instead of class
namespace StringHash64 {

namespace internal {

// CRC-64-ECMA-182 parameters in normal (MSB-first) representation.
inline constexpr std::uint64_t POLY   = 0x42F0E1EBA9EA3693ULL;
inline constexpr std::uint64_t INIT   = 0x0000000000000000ULL;
inline constexpr std::uint64_t XOR_TABLE_SIZE = 256;


// Normal (MSB-first) lookup table.  Each entry is the remainder of the byte
// placed in the most-significant byte position followed by 56 zero bits.
[[nodiscard]] consteval std::array<std::uint64_t, XOR_TABLE_SIZE> generate_crc64_table() {

    std::array<std::uint64_t, XOR_TABLE_SIZE> table{};

    for (std::uint64_t i = 0; i < XOR_TABLE_SIZE; ++i) {

        std::uint64_t crc = i << 56;

        for (std::size_t j = 0; j < 8; ++j) {

            if (crc & (1ULL << 63)) {

                crc = (crc << 1) ^ POLY;

            } else {

                crc <<= 1;

            }

        }

        table[i] = crc;

    }

    return table;
}


inline constexpr std::array<std::uint64_t, XOR_TABLE_SIZE> crc64_table = generate_crc64_table();

} // namespace internal

/// Computes the CRC-64 checksum of a given string.
[[nodiscard]] constexpr std::uint64_t CRC_64_ECMA_182(std::string_view str_text) noexcept {
    // Initial value according to ECMA-182 standard is 0x0000000000000000
    std::uint64_t crc = internal::INIT;

    for (auto text_char : str_text) {
        // Extract topmost byte of the current CRC and combine it with input byte
        auto index = static_cast<std::uint8_t>((crc >> 56) ^ static_cast<std::uint8_t>(text_char));
        crc = (crc << 8) ^ internal::crc64_table[index];

    }

    // Final XOR value for ECMA-182 is 0x0000000000000000 (no final modification)
    return crc;
}

/// FNV-1a 64-bit compile-time string hashing.
[[nodiscard]] constexpr std::uint64_t FNV_1A(std::string_view str_text) noexcept {

    std::uint64_t hash = 0xcbf29ce484222325ULL; // FNV offset basis

    for (auto text_char : str_text) {

        hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(text_char));
        hash *= 0x00000100000001B3ULL; // FNV prime

    }

    return hash;

}


} // namespace StringHash64


// Assert correct CRC-64-ECMA-182 crc values
static_assert( StringHash64::CRC_64_ECMA_182("") == 0x0000000000000000ULL, "CRC_64_ECMA_182 empty-string check vector mismatch");
static_assert( StringHash64::CRC_64_ECMA_182("123456789") == 0x6C40DF5F0B497347ULL, "CRC_64_ECMA_182 test vector mismatch");
// Assert correct FNV_1A hash values
static_assert( StringHash64::FNV_1A("") == 0xcbf29ce484222325ULL, "FNV_1A empty-string check vector mismatch");
static_assert( StringHash64::FNV_1A("Hello World") == 0x3d58dee72d4e0c27ULL, "FNV_1A test vector mismatch");


} // namespace kellerberrin

#endif //KEL_STRING_CRC64HASH_H
