//
// kgl_sequence.h — owning Sequence<Policy> over std::vector<Symbol>.
//
// Composition, not inheritance: Sequence builds a SequenceView on demand, so there is no
// stored span to invalidate and move operations stay defaulted. Mutation returns a
// bool-compatible MutationResult so both `if (not seq.append(x))` and `bool ok = ...` compile.
//

#ifndef KGL_SEQUENCE_H
#define KGL_SEQUENCE_H


#include <algorithm>
#include <array>
#include <compare>
#include <concepts>
#include <cstddef>
#include <expected>
#include <iterator>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "kgl_genome_types.h"
#include "kgl_alphabet.h"
#include "kgl_sequence_view.h"
#include "kel_interval_unsigned.h"
#include "kel_exec_env.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// Forward declaration (defined in kgl_sequence_dna.h).
class DNA5SequenceCoding;

// The mutation failure reasons.
enum class SeqError { OffsetOutOfBounds, IntervalNotContained, InsertOffsetInvalid, EmptySource };

[[nodiscard]] constexpr std::string_view toString(SeqError error) noexcept {

  switch (error) {
    case SeqError::OffsetOutOfBounds:  return "offset out of bounds";
    case SeqError::IntervalNotContained: return "interval not contained";
    case SeqError::InsertOffsetInvalid: return "insert offset invalid";
    case SeqError::EmptySource: return "empty source";
  }
  std::unreachable();   // SeqError is a closed enum; no trailing-return noise.

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MutationResult — a bool-compatible std::expected<void, SeqError>.
// std::expected's operator bool is explicit, so `bool ok = seq.append(x)` would not compile;
// this thin adapter restores the implicit conversion the reference's bool-returning API had.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MutationResult {

public:

  MutationResult() noexcept = default;
  MutationResult(std::expected<void, SeqError> result) noexcept : result_(result) {}
  MutationResult(std::unexpected<SeqError> error) noexcept : result_(error) {}

  [[nodiscard]] bool has_value() const noexcept { return result_.has_value(); }
  [[nodiscard]] operator bool() const noexcept { return result_.has_value(); }   // implicit (bool-compatible)
  [[nodiscard]] SeqError error() const noexcept { return result_.error(); }

private:

  std::expected<void, SeqError> result_{};

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sequence<Policy> — the owning sequence.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Policy>
class Sequence {

public:

  using Symbol = typename Policy::Alphabet;

  Sequence() = default;
  explicit Sequence(std::string_view ascii_string) { convertFromCharString(ascii_string); }
  explicit Sequence(std::vector<Symbol>&& symbols) noexcept : bases_(std::move(symbols)) {}
  Sequence(Sequence&& moved) noexcept = default;
  Sequence& operator=(Sequence&& moved) noexcept = default;
  Sequence(const Sequence&) = delete;                 // move-only policy retained
  Sequence& operator=(const Sequence&) = delete;
  ~Sequence() = default;

  /// Explicit deep copy (replaces the reference's implicit AlphabetString copy).
  [[nodiscard]] Sequence clone() const { return Sequence(std::vector<Symbol>(bases_)); }

  // ---- observers (on-demand view) ----
  [[nodiscard]] SequenceView<Policy> getView() const noexcept { return SequenceView<Policy>(*this); }
  [[nodiscard]] SequenceView<Policy> view() const noexcept { return getView(); }
  [[nodiscard]] ContigSize_t length() const noexcept { return bases_.size(); }
  [[nodiscard]] bool empty() const noexcept { return bases_.empty(); }
  [[nodiscard]] OpenRightUnsigned interval() const noexcept { return {0, bases_.size()}; }
  [[nodiscard]] Symbol operator[](ContigOffset_t offset) const noexcept { return bases_[offset]; }
  [[nodiscard]] std::optional<Symbol> at(ContigOffset_t offset) const noexcept {
    if (offset >= bases_.size()) {
      return std::nullopt;
    }
    return bases_[offset];
  }
  [[nodiscard]] std::string_view getStringView() const noexcept { return getView().getStringView(); }
  [[nodiscard]] std::string_view ascii() const noexcept { return getStringView(); }
  [[nodiscard]] const Symbol* data() const noexcept { return bases_.data(); }
  [[nodiscard]] std::string str() const { return std::string(getStringView()); }

  void push_back(Symbol symbol) { bases_.push_back(symbol); }
  void reserve(ContigSize_t size) { bases_.reserve(size); }
  void clear() { bases_.clear(); }

  // ---- owning-type observer forwarders (delegate to the view) ----
  [[nodiscard]] std::vector<ContigOffset_t> findAll(SequenceView<Policy> sub_sequence) const noexcept {
    return getView().findAll(sub_sequence);
  }
  [[nodiscard]] std::size_t commonPrefix(SequenceView<Policy> cmp_sequence) const noexcept {
    return getView().commonPrefix(cmp_sequence);
  }
  [[nodiscard]] std::size_t commonSuffix(SequenceView<Policy> cmp_sequence) const noexcept {
    return getView().commonSuffix(cmp_sequence);
  }
  [[nodiscard]] bool compareSubSequence(ContigOffset_t offset, SequenceView<Policy> sub_sequence) const noexcept {

    auto sub_view_opt = getView().getSubView({offset, offset + sub_sequence.length()});
    if (not sub_view_opt) {
      return false;
    }
    return sub_view_opt.value() == sub_sequence;

  }
  [[nodiscard]] bool compareLetter(ContigOffset_t offset, Symbol symbol) const noexcept {

    if (offset < bases_.size()) {
      return bases_[offset] == symbol;
    }
    return false;

  }
  [[nodiscard]] bool verifyString() const {

    for (std::size_t index = 0; index < bases_.size(); ++index) {

      if (not Policy::validAlphabet(bases_[index])) {

        ExecEnv::log().error("AlphabetString::verifyString(), invalid Alphabet value (int): {} found at index: {}", static_cast<std::size_t>(bases_[index]), index);
        return false;

      }

    }
    return true;

  }

  // ---- mutation: one bounds policy (reject, never clamp, never log inside) ----
  // Every call consumes or returns the result; none is [[nodiscard]]-ignored (v4).
  [[nodiscard]] MutationResult append(const Sequence& inserted_sequence) {

    bases_.insert(bases_.end(), inserted_sequence.bases_.begin(), inserted_sequence.bases_.end());
    return {};

  }
  [[nodiscard]] MutationResult append(SequenceView<Policy> inserted_sequence) {

    bases_.insert(bases_.end(), inserted_sequence.span().begin(), inserted_sequence.span().end());
    return {};

  }
  [[nodiscard]] MutationResult insertSubSequence(ContigOffset_t insert_offset, const Sequence& inserted_sequence) {

    if (insert_offset > bases_.size()) {
      return std::unexpected(SeqError::InsertOffsetInvalid);
    }
    bases_.insert(bases_.begin() + insert_offset, inserted_sequence.bases_.begin(), inserted_sequence.bases_.end());
    return {};

  }
  [[nodiscard]] MutationResult insertSubSequence(ContigOffset_t insert_offset, SequenceView<Policy> inserted_sequence) {

    if (insert_offset > bases_.size()) {
      return std::unexpected(SeqError::InsertOffsetInvalid);
    }
    bases_.insert(bases_.begin() + insert_offset, inserted_sequence.span().begin(), inserted_sequence.span().end());
    return {};

  }
  [[nodiscard]] MutationResult deleteSubSequence(const OpenRightUnsigned& delete_interval) {

    if (delete_interval.upper() > bases_.size()) {
      return std::unexpected(SeqError::IntervalNotContained);
    }
    bases_.erase(bases_.begin() + delete_interval.lower(), bases_.begin() + delete_interval.upper());
    return {};

  }
  [[nodiscard]] MutationResult modifyBase(ContigOffset_t base_offset, Symbol nucleotide) {

    if (base_offset >= bases_.size()) {
      return std::unexpected(SeqError::OffsetOutOfBounds);
    }
    bases_[base_offset] = nucleotide;
    return {};

  }

  // ---- subsequence / prefix / suffix ----
  [[nodiscard]] std::optional<Sequence> subSequence(const OpenRightUnsigned& sub_interval) const {

    auto sub_sequence_opt = getSubsequence(sub_interval);
    if (not sub_sequence_opt) {
      return std::nullopt;
    }
    return sub_sequence_opt;

  }
  [[nodiscard]] std::optional<Sequence> getSubsequence(const OpenRightUnsigned& sub_interval) const {

    if (not interval().containsInterval(sub_interval)) {

      ExecEnv::log().warn("Sub interval: {} not contained in interval: {}.", sub_interval.toString(), interval().toString());
      return std::nullopt;

    }
    return Sequence(std::vector<Symbol>(bases_.begin() + sub_interval.lower(), bases_.begin() + sub_interval.upper()));

  }
  [[nodiscard]] Sequence removePrefix(ContigSize_t prefix_size) const noexcept {
    return removePrefixSuffix(prefix_size, 0);
  }
  [[nodiscard]] Sequence removeSuffix(ContigSize_t suffix_size) const noexcept {
    return removePrefixSuffix(0, suffix_size);
  }
  [[nodiscard]] Sequence removePrefixSuffix(ContigSize_t prefix_size, ContigSize_t suffix_size) const noexcept {

    auto from_iter = std::ranges::next(bases_.begin(), prefix_size, bases_.end());
    auto to_iter = std::ranges::prev(bases_.end(), suffix_size, bases_.begin());
    if (std::distance(from_iter, to_iter) > 0) {
      return Sequence(std::vector<Symbol>(from_iter, to_iter));
    }
    return Sequence();

  }

  // ---- statistics ----
  [[nodiscard]] std::vector<std::pair<Symbol, std::size_t>> countSymbols() const {

    const auto& alphabet = Policy::enumerateAlphabet();
    std::vector<std::pair<Symbol, std::size_t>> symbol_count_vector;
    symbol_count_vector.reserve(alphabet.size());
    for (auto const symbol : alphabet) {
      symbol_count_vector.emplace_back(symbol, 0);
    }
    for (auto const symbol : bases_) {
      symbol_count_vector[Policy::symbolToColumn(symbol)].second++;
    }
    return symbol_count_vector;

  }
  [[nodiscard]] std::size_t countTwoSymbols(Symbol first_symbol, Symbol second_symbol) const noexcept {

    std::size_t count{0};
    for (std::size_t index = 1; index < bases_.size(); ++index) {

      if (bases_[index - 1] == first_symbol and bases_[index] == second_symbol) {
        ++count;
        ++index;   // skip to examine the next two symbols
      }

    }
    return count;

  }

  /// Concatenate the sub-views extracted in interval (sorted) order.
  [[nodiscard]] std::optional<Sequence> concatSequences(const IntervalSetLower& interval_set) const {

    std::optional<Sequence> concatenated_sequence;
    for (auto const& sub_interval : interval_set) {

      auto sub_view_opt = getView().getSubView(sub_interval);
      if (not sub_view_opt) {

        ExecEnv::log().warn("Unable to extract sub-sequence: {} for interval: {}", sub_interval.toString(), interval().toString());
        return std::nullopt;

      }

      if (concatenated_sequence) {

        // The reference logged and rejected when the concatenating append failed (v4: the
        // result is consumed — never [[nodiscard]]-ignored).
        if (not concatenated_sequence->append(sub_view_opt.value())) {

          ExecEnv::log().warn("Unable to concatenate modified sequence for interval");
          return std::nullopt;

        }

      } else {

        concatenated_sequence.emplace(std::vector<Symbol>(sub_view_opt->span().begin(), sub_view_opt->span().end()));

      }

    }

    if (not concatenated_sequence) {

      ExecEnv::log().warn("No concat sub-sequences for interval set size: {}", interval_set.size());
      return std::nullopt;

    }

    return concatenated_sequence;

  }

  // ---- comparison ----
  [[nodiscard]] auto operator<=>(const Sequence& rhs) const noexcept { return getStringView() <=> rhs.getStringView(); }
  [[nodiscard]] bool operator==(const Sequence& rhs) const noexcept { return getStringView() == rhs.getStringView(); }

  // ---- DNA5-only retained spellings (constrained members; defined in kgl_sequence_dna.h) ----
  [[nodiscard]] DNA5SequenceCoding codingSequence(StrandSense strand) const
    requires std::same_as<Policy, DNA5>;
  [[nodiscard]] static Sequence downConvertToLinear(const DNA5SequenceCoding& stranded_sequence)
    requires std::same_as<Policy, DNA5>;

private:

  std::vector<Symbol> bases_;

  void convertFromCharString(std::string_view alphabet_str) {

    bases_.reserve(alphabet_str.length());

    // Aggregate diagnostics (improvement over per-character logging): one report per string.
    // The conversion loop tallies as it converts, so the string is scanned once (v4).
    ParseReport report;
    for (const char chr : alphabet_str) {
      bases_.push_back(Policy::convertChar(chr));
      if constexpr (std::same_as<Policy, DNA5> or std::same_as<Policy, CodingDNA5>) {
        detail::tallyInvalidDNA5(report, chr);
      } else if constexpr (std::same_as<Policy, AminoAcid>) {
        detail::tallyInvalidAminoAcid(report, chr);
      }
    }

    if constexpr (std::same_as<Policy, DNA5> or std::same_as<Policy, CodingDNA5>) {

      if (report.extended_chars > 0) {
        ExecEnv::log().warn("Sequence: IUPAC extended nucleotides detected, all converted to the unknown nucleotide 'N'");
      }
      if (report.invalid_chars > 0) {
        ExecEnv::log().error("Sequence: {} unknown nucleotide(s) detected, all converted to 'N'. Input is probably corrupt or not DNA text.",
                             report.invalid_chars);
      }

    } else if constexpr (std::same_as<Policy, AminoAcid>) {

      if (report.invalid_chars > 0) {
        ExecEnv::log().error("Sequence: {} invalid amino acid(s) detected, all converted to 'Z'. Input is probably corrupt or not protein text.",
                             report.invalid_chars);
      }

    }

  }

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Deferred SequenceView definitions (both types now complete).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Policy>
constexpr SequenceView<Policy>::SequenceView(const Sequence<Policy>& sequence) noexcept
  : view_(sequence.data(), sequence.length()) {}

template<class Policy>
Sequence<Policy> SequenceView<Policy>::toSequence() const {
  return Sequence<Policy>(std::vector<Symbol>(view_.begin(), view_.end()));
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Concrete type aliases (both historical spellings denote the same instantiation).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

using StringDNA5 = Sequence<DNA5>;
using StringCodingDNA5 = Sequence<CodingDNA5>;
using StringAminoAcid = Sequence<AminoAcid>;
using DNA5SequenceLinear = Sequence<DNA5>;
using DNA5SequenceLinearView = SequenceView<DNA5>;
using AminoSequence = Sequence<AminoAcid>;
using AminoSequenceView = SequenceView<AminoAcid>;


}   // end namespace


#endif //KGL_SEQUENCE_H
