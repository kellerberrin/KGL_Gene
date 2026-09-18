//
// kgl_sequence_dna.h — strand subclass, DNA conversions and the constrained DNA5 members.
//
// The constrained member definitions are out-of-line here so that both Sequence<Policy>
// and DNA5SequenceCoding are complete. This is what keeps the reference's member spellings
// (`linear.codingSequence(strand)`) and static-member spelling
// (`DNA5SequenceLinear::downConvertToLinear(c)`) compiling with no downstream edits.
//

#ifndef KGL_SEQUENCE_DNA_H
#define KGL_SEQUENCE_DNA_H


#include <optional>
#include <span>
#include <vector>

#include "kgl_sequence.h"
#include "kgl_genome_prelim.h"


namespace kellerberrin::genome {   //  organization::project level namespace


template<class Policy> class Sequence;
template<class Policy> class SequenceView;
class DNA5SequenceCoding;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA5SequenceCodingView — a SequenceView<CodingDNA5> that also carries the strand.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DNA5SequenceCodingView : public SequenceView<CodingDNA5> {

public:

  DNA5SequenceCodingView() = default;
  DNA5SequenceCodingView(std::span<const Nucleotide> symbols, StrandSense strand) noexcept
    : SequenceView<CodingDNA5>(symbols), strand_(strand) {}
  DNA5SequenceCodingView(const DNA5SequenceCodingView& copy) = default;
  ~DNA5SequenceCodingView() = default;

  /// Returns a subsequence view or std::nullopt if the interval is out of bounds.
  [[nodiscard]] std::optional<DNA5SequenceCodingView> subView(const OpenRightUnsigned& sub_interval) const;

  /// Returns the sequence strand, FORWARD '+' or REVERSE '-'.
  [[nodiscard]] StrandSense strand() const noexcept { return strand_; }

private:

  StrandSense strand_{StrandSense::FORWARD};

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA5SequenceCoding — a STRANDED DNA sequence that can be translated to protein.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DNA5SequenceCoding : public Sequence<CodingDNA5> {

public:

  DNA5SequenceCoding() = default;
  DNA5SequenceCoding(Sequence<CodingDNA5>&& bases, StrandSense strand) noexcept
    : Sequence<CodingDNA5>(std::move(bases)), strand_(strand) {}
  DNA5SequenceCoding(const DNA5SequenceCodingView& sequence_view);
  DNA5SequenceCoding(DNA5SequenceCoding&&) noexcept = default;
  DNA5SequenceCoding& operator=(DNA5SequenceCoding&&) noexcept = default;
  DNA5SequenceCoding(const DNA5SequenceCoding&) = delete;   // move-only policy retained
  DNA5SequenceCoding& operator=(const DNA5SequenceCoding&) = delete;
  ~DNA5SequenceCoding() = default;

  /// ReturnType a view that also carries the strand.
  [[nodiscard]] DNA5SequenceCodingView getView() const noexcept {
    return DNA5SequenceCodingView(std::span<const Nucleotide>(this->data(), this->length()), strand_);
  }

  /// Returns the sequence strand, FORWARD '+' or REVERSE '-'.
  [[nodiscard]] StrandSense strand() const noexcept { return strand_; }

private:

  StrandSense strand_{StrandSense::FORWARD};

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Free conversions. The symbol domains are identical so FORWARD and down-conversion are
// re-tags (O(n) move, no per-base cast); REVERSE is the reverse complement.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

[[nodiscard]] inline DNA5SequenceCoding codingSequence(const Sequence<DNA5>& linear, StrandSense strand) {

  std::vector<Nucleotide> output;
  output.reserve(linear.length());

  if (strand == StrandSense::REVERSE) {

    for (const Nucleotide nucleotide : linear.getView().span() | std::views::reverse) {
      output.push_back(DNA5::complementNucleotide(nucleotide));
    }

  } else {

    output.assign(linear.getView().span().begin(), linear.getView().span().end());

  }

  return DNA5SequenceCoding(Sequence<CodingDNA5>(std::move(output)), strand);

}

[[nodiscard]] inline Sequence<DNA5> asLinear(const DNA5SequenceCoding& coding) {

  return Sequence<DNA5>(std::vector<Nucleotide>(coding.getView().span().begin(),
                                                coding.getView().span().end()));

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constrained DNA5-only member definitions (both types complete).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Policy>
DNA5SequenceCoding Sequence<Policy>::codingSequence(StrandSense strand) const
  requires std::same_as<Policy, DNA5>
{
  return ::kellerberrin::genome::codingSequence(*this, strand);
}

template<class Policy>
Sequence<Policy> Sequence<Policy>::downConvertToLinear(const DNA5SequenceCoding& stranded_sequence)
  requires std::same_as<Policy, DNA5>
{
  return asLinear(stranded_sequence);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Deferred view/strand definitions.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline std::optional<DNA5SequenceCodingView> DNA5SequenceCodingView::subView(const OpenRightUnsigned& sub_interval) const {

  auto sub_view_opt = getSubView(sub_interval);
  if (not sub_view_opt) {
    return std::nullopt;
  }
  return DNA5SequenceCodingView(sub_view_opt->span(), strand_);

}

// View-side conversion spelling: DNA5SequenceLinearView::codingSequence(strand).
template<class Policy>
DNA5SequenceCoding SequenceView<Policy>::codingSequence(StrandSense strand) const
  requires std::same_as<Policy, DNA5>
{
  return ::kellerberrin::genome::codingSequence(this->toSequence(), strand);
}

inline DNA5SequenceCoding::DNA5SequenceCoding(const DNA5SequenceCodingView& sequence_view)
  : Sequence<CodingDNA5>(sequence_view.toSequence()), strand_(sequence_view.strand()) {}


}   // end namespace


#endif //KGL_SEQUENCE_DNA_H
