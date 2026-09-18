//
// kgl_sequence_view.h — non-owning SequenceView<Policy>, concepts and the fasta SequenceRef.
//
// A SequenceView does not own its data; the referenced sequence must outlive the view.
// It is trivially copyable and assignable (the reference's const data member made views
// non-assignable). Comparison is explicitly delegated to the byte view because std::span
// has no relational operators in C++20/23.
//

#ifndef KGL_SEQUENCE_VIEW_H
#define KGL_SEQUENCE_VIEW_H


#include <algorithm>
#include <compare>
#include <concepts>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

#include "kgl_genome_types.h"
#include "kgl_alphabet.h"
#include "kel_interval_unsigned.h"
#include "kel_search.h"
#include "kel_exec_env.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// Forward declarations.
template<class Policy> class Sequence;
class DNA5SequenceCoding;
enum class StrandSense : char;   // defined in kgl_genome_prelim.h


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SequenceView<Policy> — a non-owning view over contiguous symbols.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Policy>
class SequenceView {

public:

  using Symbol = typename Policy::Alphabet;

  constexpr SequenceView() noexcept = default;
  constexpr SequenceView(const Symbol* data, std::size_t size) noexcept : view_(data, size) {}
  constexpr SequenceView(std::span<const Symbol> span) noexcept : view_(span) {}
  constexpr SequenceView(const Sequence<Policy>& sequence) noexcept;

  SequenceView(const SequenceView& copy) noexcept = default;
  SequenceView& operator=(const SequenceView& copy) noexcept = default;   // assignable (was deleted by const member)
  ~SequenceView() = default;

  /// Random access to the view symbols (unchecked; precondition offset < length).
  [[nodiscard]] constexpr Symbol operator[](ContigOffset_t offset) const noexcept { return view_[offset]; }
  /// Random access with genuine bounds checking; empty optional if out of range.
  [[nodiscard]] constexpr std::optional<Symbol> at(ContigOffset_t offset) const noexcept {
    if (offset >= view_.size()) {
      return std::nullopt;
    }
    return view_[offset];
  }

  [[nodiscard]] constexpr ContigSize_t length() const noexcept { return view_.size(); }
  [[nodiscard]] constexpr bool empty() const noexcept { return view_.empty(); }
  [[nodiscard]] constexpr OpenRightUnsigned interval() const noexcept { return {0, view_.size()}; }

  /// The contiguous symbols as bytes (all alphabets are byte sized).
  [[nodiscard]] constexpr std::span<const Symbol> span() const noexcept { return view_; }

  /// Zero-copy std::string_view (name retained — the most-used downstream accessor).
  [[nodiscard]] constexpr std::string_view getStringView() const noexcept {
    return std::string_view{reinterpret_cast<const char*>(view_.data()), view_.size()};
  }
  [[nodiscard]] constexpr std::string_view ascii() const noexcept { return getStringView(); }

  /// Extract a subsequence view, or std::nullopt if the interval is not contained.
  [[nodiscard]] std::optional<SequenceView> getSubView(const OpenRightUnsigned& sub_interval) const noexcept {

    if (not interval().containsInterval(sub_interval)) {

      ExecEnv::log().warn("Sub interval: {} not contained in interval: {}.", sub_interval.toString(), interval().toString());
      return std::nullopt;

    }
    return SequenceView(view_.subspan(sub_interval.lower(), sub_interval.size()));

  }

  /// Returns the intersection of the view and the interval.
  [[nodiscard]] SequenceView getIntersection(const OpenRightUnsigned& sub_interval) const noexcept {

    const auto truncate_interval = interval().intersection(sub_interval);
    return SequenceView(view_.subspan(truncate_interval.lower(), truncate_interval.size()));

  }

  /// Search for all (overlapping) subsequences.
  [[nodiscard]] std::vector<ContigOffset_t> findAll(SequenceView sub_sequence) const noexcept {

    std::vector<ContigOffset_t> offset_vector;
    std::string_view haystack = getStringView();
    const std::string_view needle = sub_sequence.getStringView();
    if (needle.empty()) {
      return offset_vector;
    }

    std::size_t offset = haystack.find(needle);
    while (offset != std::string_view::npos) {

      offset_vector.push_back(static_cast<ContigOffset_t>(offset));
      offset = haystack.find(needle, offset + 1);

    }

    return offset_vector;

  }

  /// Longest common prefix length.
  [[nodiscard]] constexpr std::size_t commonPrefix(SequenceView cmp_view) const noexcept {

    const auto [this_iter, cmp_iter] = std::mismatch(view_.begin(), view_.end(), cmp_view.view_.begin(), cmp_view.view_.end());
    return static_cast<std::size_t>(std::distance(view_.begin(), this_iter));

  }

  /// Longest common suffix length.
  [[nodiscard]] constexpr std::size_t commonSuffix(SequenceView cmp_view) const noexcept {

    const auto [this_iter, cmp_iter] = std::mismatch(view_.rbegin(), view_.rend(), cmp_view.view_.rbegin(), cmp_view.view_.rend());
    return static_cast<std::size_t>(std::distance(view_.rbegin(), this_iter));

  }

  /// Explicit materialisation of the view into an owning sequence.
  [[nodiscard]] Sequence<Policy> toSequence() const;

  /// Retained view-side conversion spelling (DNA5 only; defined in kgl_sequence_dna.h).
  [[nodiscard]] DNA5SequenceCoding codingSequence(StrandSense strand) const
    requires std::same_as<Policy, DNA5>;

  /// Lexical comparison over the byte view (std::span has no relational operators).
  [[nodiscard]] constexpr auto operator<=>(const SequenceView& rhs) const noexcept {
    return getStringView() <=> rhs.getStringView();
  }
  [[nodiscard]] constexpr bool operator==(const SequenceView& rhs) const noexcept {
    return getStringView() == rhs.getStringView();
  }

protected:

  std::span<const Symbol> view_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Concepts.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class S>
concept SequenceLike = requires(const S& s) {
  { s.getStringView() } -> std::convertible_to<std::string_view>;
  { s.length() } -> std::convertible_to<std::size_t>;
};

template<class S, class Policy>
concept SequenceOf = SequenceLike<S> and std::same_as<typename S::Symbol, typename Policy::Alphabet>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Free regex search — replaces VirtualSequence::regexSearch. A cached-regex overload avoids
// rebuilding std::regex for repeated motif searches.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

[[nodiscard]] inline std::vector<OpenRightUnsigned> regexSearch(SequenceLike auto const& sequence,
                                                                std::string_view regex) {
  return Search::searchView(regex, sequence.getStringView());
}

[[nodiscard]] inline std::vector<OpenRightUnsigned> regexSearch(SequenceLike auto const& sequence,
                                                                const std::regex& compiled) {
  return Search::searchView(compiled, sequence.getStringView());
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SequenceRef — type-erased sequence reference for the fasta boundary. The owning overload
// keeps a shared_ptr source alive; the borrowing overload requires the caller to ensure lifetime.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SequenceRef {

public:

  // Owning overload. `chars_` is initialised before `keep_alive_` moves the owner, so the view
  // is taken from the live pointer (avoids a moved-from dereference).
  template<class S> requires SequenceLike<S>
  SequenceRef(std::shared_ptr<S> owner) noexcept
    : chars_(owner ? owner->getStringView() : std::string_view{}),
      keep_alive_(std::move(owner)) {}

  template<class S> requires SequenceLike<S>
  explicit SequenceRef(const S& sequence) noexcept : chars_(sequence.getStringView()) {}

  [[nodiscard]] std::string_view getStringView() const noexcept { return chars_; }

private:

  std::string_view chars_;
  std::shared_ptr<const void> keep_alive_;   // null for non-owning references

};


}   // end namespace


#endif //KGL_SEQUENCE_VIEW_H
