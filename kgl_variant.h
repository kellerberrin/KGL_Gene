// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// Created by kellerberrin on 13/10/17.
//

#ifndef KGL_VARIANT_H
#define KGL_VARIANT_H

#include <map>
#include <memory>
#include "kgl_genome_types.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Base class for a genome variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class Variant;
class ReadCountVariant;
// Variant filters use the vistor pattern.
class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  bool applyFilter(Variant& variant) const { return true; }  // Filter is ignored by variant.
  virtual bool applyFilter(ReadCountVariant& variant) const = 0;

private:


};



class Variant {

public:


  Variant(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  virtual ~Variant() = default;
  bool filterVariant(const VariantFilter& filter) { return applyFilter(filter); }

private:

  ContigId_t contig_id_;

  virtual bool applyFilter(const VariantFilter& filter) = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A read count variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ReadCountVariant : public Variant {

public:

  ReadCountVariant(const ContigId_t& contig_id, ContigOffset_t contig_offset, NucleotideReadCount_t read_count)
      : Variant(contig_id), contig_offset_(contig_offset), read_count_(read_count) {}
  ~ReadCountVariant() override = default;

  ContigOffset_t SNPOffset() const { return contig_offset_; }
  NucleotideReadCount_t readCount() const { return read_count_; }

private:

  ContigOffset_t contig_offset_;
  NucleotideReadCount_t read_count_;

  bool applyFilter(const VariantFilter& filter) override { return filter.applyFilter(*this); }

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SNPVariant : public ReadCountVariant {

public:

  SNPVariant(const ContigId_t& contig_id, ContigOffset_t contig_offset, NucleotideReadCount_t read_count)
      : ReadCountVariant(contig_id, contig_offset, read_count) {}
  ~SNPVariant() override = default;

private:

  bool applyFilter(const VariantFilter& filter) final { return true; }

};


class ReadCountFilter : public VariantFilter {

public:

  explicit ReadCountFilter(NucleotideReadCount_t read_count) : read_count_(read_count) {}
  ~ReadCountFilter() override = default;

  bool applyFilter(ReadCountVariant& variant) const override { return variant.readCount() >= read_count_; }

private:

  NucleotideReadCount_t read_count_;

};


class MutantProportionFilter : public VariantFilter {

public:

  explicit MutantProportionFilter(double proportion) : mutant_proportion_(proportion) {}
  ~MutantProportionFilter() override = default;

  bool applyFilter(ReadCountVariant& variant) const override { return true; }

private:

  double mutant_proportion_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<Variant>>;
class ContigVariant {

public:

  ContigVariant(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  ContigVariant(const ContigVariant&) = default;
  ~ContigVariant() = default;

  ContigVariant& operator=(const ContigVariant&) = default;

  void addVariant(ContigOffset_t contig_offset, std::shared_ptr<Variant>& variant_ptr);
  const ContigId_t& contigId() const { return contig_id_; }
  size_t variantCount() const { return offset_variant_map_.size(); }

  size_t filterVariants(const VariantFilter& filter);

private:

  ContigId_t contig_id_;
  OffsetVariantMap offset_variant_map_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenomeVariantMap = std::map<ContigId_t, std::shared_ptr<ContigVariant>>;
class GenomeVariant {

public:

  explicit GenomeVariant() {}
  GenomeVariant(const GenomeVariant&) = default;
  ~GenomeVariant() = default;

  GenomeVariant& operator=(const GenomeVariant&) = default;

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  void filterVariants(const VariantFilter& filter);

private:

  GenomeVariantMap genome_variant_map_;

};




}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
