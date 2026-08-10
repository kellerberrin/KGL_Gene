//
// Created by kellerberrin on 10/8/26.
//

#ifndef KEL_STRING_CRC64HASH_H
#define KEL_STRING_CRC64HASH_H

#include <iostream>
#include <array>
#include <string_view>
#include <cstdint>



namespace kellerberrin {
//  organization level namespace

class StringHash64Table {
public:

  StringHash64Table() = delete;
  ~StringHash64Table() = delete;

  // 1. Generate the CRC-64 lookup table at compile time
  static constexpr std::array<uint64_t, 256> generate_crc64_table() {

    uint64_t polynomial = 0x42F0E1EBA9EA3693ULL;
    std::array<uint64_t, 256> table{};

    for (uint32_t i = 0; i < 256; ++i) {
      uint64_t crc = i;
      for (uint32_t j = 0; j < 8; ++j) {
        if (crc & 1) {
          crc = (crc >> 1) ^ polynomial;
        } else {
          crc >>= 1;
        }
      }
      table[i] = crc;
    }

    return table;

  }

}; // Helper class creates the 64 bit hash table.

class UtilityStringHash64 {


  // 2. Global compile-time table instantiation
  constexpr static std::array<uint64_t, 256> crc64_table = StringHash64Table::generate_crc64_table();

public:

  UtilityStringHash64() = delete;
  ~UtilityStringHash64() = delete;


  // 3. Compile-time capable CRC-64 calculation function
  [[nodiscard]] static constexpr uint64_t crc64StringHash(std::string_view data) {

    uint64_t crc = 0xFFFFFFFFFFFFFFFFULL; // Initial value

    for (char c : data) {
      uint8_t index = static_cast<uint8_t>(crc ^ c);
      crc = (crc >> 8) ^ crc64_table[index];
    }

    return crc ^ 0xFFFFFFFFFFFFFFFFULL; // Final XOR
  }


}; // UtilityStringHash64



}  // kellerberrin namespace

#endif //KEL_STRING_CRC64HASH_H
