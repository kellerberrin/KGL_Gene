//
// Created by kellerberrin on 10/8/26.
//

#ifndef KEL_ENTROPY_SOURCE_H
#define KEL_ENTROPY_SOURCE_H


#include <algorithm>
#include <array>
#include <functional>
#include <random>

namespace kellerberrin {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Random number generators and Pdfs for the Uniform, Discrete Uniform, Gaussian, Binomial and BetaBinomial distributions.
// Note that object copy semantics have been disabled as there seems no reasonable reason to do this,
// and in the case of the random number generators, is probably dangerous.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Use the 64 bit Mersenne Twister as the random number entropy source.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using EntropyGenerator = std::mt19937_64;

// Production entropy source: properly seeded 64-bit Mersenne Twister.
class RandomEntropySource {
public:
    using result_type = EntropyGenerator::result_type;

    RandomEntropySource() : generator_(make_seeded_engine()) {}
    explicit RandomEntropySource(std::seed_seq& seed_sequence) : generator_(seed_sequence) {}

    RandomEntropySource(const RandomEntropySource&) = delete;
    RandomEntropySource(RandomEntropySource&&) = default;
    RandomEntropySource& operator=(const RandomEntropySource&) = delete;
    RandomEntropySource& operator=(RandomEntropySource&&) = default;
    ~RandomEntropySource() = default;

    [[nodiscard]] EntropyGenerator& generator() { return generator_; }
    [[nodiscard]] result_type operator()() { return generator_(); }

    [[nodiscard]] static constexpr result_type min() { return EntropyGenerator::min(); }
    [[nodiscard]] static constexpr result_type max() { return EntropyGenerator::max(); }

private:
    EntropyGenerator generator_;

    [[nodiscard]] static EntropyGenerator make_seeded_engine() {
        std::random_device rd;
        // 64-bit mt19937_64 has 312 64-bit state words; seed it with twice as many
        // 32-bit random-device outputs to fully populate the state.
        std::array<std::random_device::result_type, EntropyGenerator::state_size * 2> seed_data{};
        std::generate(seed_data.begin(), seed_data.end(), std::ref(rd));
        std::seed_seq seq(seed_data.begin(), seed_data.end());
        return EntropyGenerator{seq};
    }
};

// Deterministic source for reproducible runs / debugging. Use explicitly when needed.
class DeterministicEntropySource {
public:
    using result_type = EntropyGenerator::result_type;

    explicit DeterministicEntropySource(std::uint64_t seed = 1111) : generator_(seed) {}

    DeterministicEntropySource(const DeterministicEntropySource&) = delete;
    DeterministicEntropySource(DeterministicEntropySource&&) = default;
    DeterministicEntropySource& operator=(const DeterministicEntropySource&) = delete;
    DeterministicEntropySource& operator=(DeterministicEntropySource&&) = default;
    ~DeterministicEntropySource() = default;

    [[nodiscard]] EntropyGenerator& generator() { return generator_; }
    [[nodiscard]] result_type operator()() { return generator_(); }

    [[nodiscard]] static constexpr result_type min() { return EntropyGenerator::min(); }
    [[nodiscard]] static constexpr result_type max() { return EntropyGenerator::max(); }

private:
    EntropyGenerator generator_;
};

// Default production alias. Pick the deterministic source explicitly when required.
//using EntropySource = DeterministicEntropySource;
using EntropySource = RandomEntropySource;


} // namespace kellerberrin


#endif //KEL_ENTROPY_SOURCE_H
