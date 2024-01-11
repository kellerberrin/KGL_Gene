//
// Created by kellerberrin on 22/02/18.
//

#ifndef KGL_SEQUENCE_DISTANCE_IMPL_H
#define KGL_SEQUENCE_DISTANCE_IMPL_H


#include <memory>
#include <string>
#include <map>
#include "kgl_genome_types.h"
#include "kgl_sequence_amino.h"


namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance is conceptually different from comparison.
// 1. Distances are always normalized for sequence length.
// For local distances, this is the size of the match.
// For global distances, this is the length of the sequences.
// 2. Distances are always positive or zero.
// 3. Distances are symmetric, d(x,y) = d(y,x)
// 4. Distances observe the triangle inequality, d(x,y) + d(y,z) >= d(x,z)
// 5. Distances are returned as a CompareDistance_t which is a double.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


CompareDistance_t LevenshteinGlobalImpl(const char* sequenceA,
                                        size_t sequenceA_size,
                                        const char* sequenceB,
                                        size_t sequenceB_size);

template<typename Seq>
CompareDistance_t LevenshteinGlobalImpl(const Seq& sequenceA, const Seq& sequenceB) {

  auto sequenceA_view = sequenceA.getStringView();
  auto sequenceB_view = sequenceB.getStringView();
  return LevenshteinGlobalImpl(sequenceA_view.data(),
                               sequenceA_view.length(),
                               sequenceB_view.data(),
                               sequenceB_view.length());

}

CompareDistance_t LevenshteinLocalImpl(const char* sequenceA,
                                       size_t sequenceA_size,
                                       const char* sequenceB,
                                       size_t sequenceB_size);

template<typename Seq>
CompareDistance_t LevenshteinLocalImpl(const Seq& sequenceA, const Seq& sequenceB) {

  auto sequenceA_view = sequenceA.getStringView();
  auto sequenceB_view = sequenceB.getStringView();
  return LevenshteinLocalImpl(sequenceA_view.data(),
                              sequenceA_view.length(),
                              sequenceB_view.data(),
                              sequenceB_view.length());

}

template<typename Seq>
CompareDistance_t globalblosum80Impl(const Seq&, const Seq&) {

  return 0.0;

}

template<typename Seq>
CompareDistance_t localblosum80Impl(const Seq&, const Seq&) {

  return 0.0;

}

// Define the parentDistance function object.
template<typename Seq>
using SequenceDistanceMetricFn = std::function<CompareDistance_t(const Seq&, const Seq&)>;

// Delete the default constructor.
// The parentDistance function object must be initialized.
template<typename Seq>
class SequenceDistanceMetric : public SequenceDistanceMetricFn<Seq> {

public:
  using SequenceDistanceMetricFn<Seq>::SequenceDistanceMetricFn;
  SequenceDistanceMetric() = delete;

};


using VirtualDistanceMetric = SequenceDistanceMetric<VirtualSequence>;

inline static const VirtualDistanceMetric LevenshteinGlobalVirtual{LevenshteinGlobalImpl<VirtualSequence>};
inline static const VirtualDistanceMetric LevenshteinLocalVirtual{LevenshteinLocalImpl<VirtualSequence>};

using AminoDistanceMetric = SequenceDistanceMetric<AminoSequence>;

inline static const AminoDistanceMetric LevenshteinGlobalAmino{LevenshteinGlobalImpl<AminoSequence>};
inline static const AminoDistanceMetric LevenshteinLocalAmino{LevenshteinLocalImpl<AminoSequence>};
inline static const AminoDistanceMetric globalblosum80Amino{globalblosum80Impl<AminoSequence>};
inline static const AminoDistanceMetric localblosum80Amino{localblosum80Impl<AminoSequence>};


using CodingDistanceMetric = SequenceDistanceMetric<DNA5SequenceCoding>;

inline static const CodingDistanceMetric LevenshteinGlobalCoding{LevenshteinGlobalImpl<DNA5SequenceCoding>};
inline static const CodingDistanceMetric LevenshteinLocalCoding{LevenshteinLocalImpl<DNA5SequenceCoding>};

using CodingDistanceMetricView = SequenceDistanceMetric<DNA5SequenceCodingView>;

inline static const CodingDistanceMetricView LevenshteinGlobalCodingView{LevenshteinGlobalImpl<DNA5SequenceCodingView>};
inline static const CodingDistanceMetricView LevenshteinLocalCodingView{LevenshteinLocalImpl<DNA5SequenceCodingView>};

using LinearDistanceMetric = SequenceDistanceMetric<DNA5SequenceLinear>;

inline static const LinearDistanceMetric LevenshteinGlobalLinear{LevenshteinGlobalImpl<DNA5SequenceLinear>};
inline static const LinearDistanceMetric LevenshteinLocalLinear{LevenshteinLocalImpl<DNA5SequenceLinear>};

}   // end namespace




#endif //KGL_SEQUENCE_DISTANCE_IMPL_H
