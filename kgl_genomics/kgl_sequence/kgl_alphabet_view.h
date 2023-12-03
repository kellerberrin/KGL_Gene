//
// Created by kellerberrin on 20/11/23.
//

#ifndef KGL_ALPHABET_VIEW_H
#define KGL_ALPHABET_VIEW_H



#include "kgl_sequence_virtual.h"

#include <string_view>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A template Sequence class to specify a view of a sequence.
// As with char strings, the view does not contain any sequence data but is an offset to an Alphabet sequence.
// The creator of the Alphabet view must ensure its lifetime does not exceed the corresponding sequence.
// Sequence views only have const operations on the underlying sequence.
// A sequence can create a view (via sub_sequence() or assignment) and a view can recreate a copy of the sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward declaration of the sequence template.
template<class alphabet> class AlphabetSequence;

template<typename Alphabet>
class AlphabetView : public VirtualSequence {

public:

  explicit AlphabetView(const AlphabetSequence<Alphabet>& sequence) : alphabet_view_(sequence.data(), sequence.length()) {}
  explicit AlphabetView(const std::basic_string_view<typename Alphabet::Alphabet>& view) : alphabet_view_(view) {}
  AlphabetView(const AlphabetView& copy) = default;

  ~AlphabetView() override = default;

  [[nodiscard]] auto operator[] (ContigOffset_t offset) const { return alphabet_view_[offset]; }
  [[nodiscard]] auto at(ContigOffset_t offset) const { return alphabet_view_[offset]; }

  [[nodiscard]] constexpr ContigSize_t length() const { return alphabet_view_.length(); }
  [[nodiscard]] OpenRightUnsigned interval() const { return {0, alphabet_view_.length() }; }

  // Assumes that sequence alphabets are 1 byte (char) and map onto ascii char types. Avoids a sequence byte copy and conversion.
  [[nodiscard]] constexpr std::string_view getStringView() const override { return std::string_view{reinterpret_cast<const char*>(alphabet_view_.data()), alphabet_view_.length()}; }

  // Get the physical sequence from the view.
  [[nodiscard]] AlphabetSequence<Alphabet> getSequence() const { return AlphabetSequence<Alphabet>(std::move(AlphabetString<Alphabet>(alphabet_view_))); }

  // Search for all subsequences.
  [[nodiscard]] std::vector<ContigOffset_t> findAll(const AlphabetView& sub_sequence) const { return alphabet_view_.findAll(sub_sequence.alphabet_view_); }

  // Ptr to the base of the alphabet string. Used to initialize a SequenceView object.
  [[nodiscard]] constexpr const Alphabet::Alphabet* data() const { return  alphabet_view_.data(); }

  // Returns std::nullopt if offset and/or size are out of bounds.
  [[nodiscard]] std::optional<AlphabetView> getSubView(const OpenRightUnsigned& sub_interval) const;

  // Returns std::nullopt if offset and/or size are out of bounds.
  [[nodiscard]] AlphabetView getTruncate(const OpenRightUnsigned& sub_interval) const;

  // Sequence view comparison using the spaceship operator. Sequence views (and thus sequences) are ordered lexically using std::string_view.
  [[nodiscard]] constexpr auto operator<=>(const AlphabetView& rhs) const { return getStringView() <=> rhs.getStringView(); }
  [[nodiscard]] constexpr bool operator==(const AlphabetView& rhs) const { return getStringView() == rhs.getStringView(); }

protected:

  const std::basic_string_view<typename Alphabet::Alphabet> alphabet_view_;

};


template<typename Alphabet>
std::optional<AlphabetView<Alphabet>> AlphabetView<Alphabet>::getSubView(const OpenRightUnsigned& sub_interval) const {

  if (not interval().containsInterval(sub_interval)) {

    ExecEnv::log().warn("Sub interval: {} not contained in interval: {}.", sub_interval.toString(), interval().toString());
    return std::nullopt;

  }

  return AlphabetView<Alphabet>(alphabet_view_.substr(sub_interval.lower(), sub_interval.size()));

}

template<typename Alphabet>
AlphabetView<Alphabet> AlphabetView<Alphabet>::getTruncate(const OpenRightUnsigned& sub_interval) const {

  auto truncate_interval = interval().intersection(sub_interval);
  return AlphabetView<Alphabet>(alphabet_view_.substr(truncate_interval.lower(), truncate_interval.size()));

}


}   // end namespace


#endif //KGL_ALPHABET_VIEW_H
