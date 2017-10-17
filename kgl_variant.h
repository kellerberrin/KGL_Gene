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
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class Variant; // Forward decl.
class ReadCountVariant; // Forward decl.
class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  bool applyFilter(const Variant& variant) const { return true; }  // Catchall for filters not defined for variant type.
  virtual bool applyFilter(const ReadCountVariant& variant) const = 0;
  virtual std::string filterName() const = 0;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Base class for a genome variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Variant {

public:


  Variant(const ContigId_t& contig_id, ContigOffset_t contig_offset) : contig_id_(contig_id),
                                                                       contig_offset_(contig_offset) {}
  virtual ~Variant() = default;
  bool filterVariant(const VariantFilter& filter) const { return applyFilter(filter); }

  const ContigId_t& contigId() const { return contig_id_; }
  ContigOffset_t contigOffset() const { return contig_offset_; }
  bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };

private:

  ContigId_t contig_id_;
  ContigOffset_t contig_offset_;

  virtual bool applyFilter(const VariantFilter& filter) const = 0;
  virtual bool equivalent(const Variant& cmp_var) const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A read count variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ReadCountVariant : public Variant {

public:

  ReadCountVariant(const ContigId_t& contig_id,
                   ContigOffset_t contig_offset,
                   NucleotideReadCount_t read_count,
                   NucleotideReadCount_t mutant_count)
      : Variant(contig_id, contig_offset), read_count_(read_count), mutant_count_(mutant_count) {}
  ~ReadCountVariant() override = default;

  NucleotideReadCount_t readCount() const { return read_count_; }
  NucleotideReadCount_t mutantCount() const { return mutant_count_; }

  double proportion() const { return static_cast<double>(mutant_count_) / static_cast<double>(read_count_); }

private:

  NucleotideReadCount_t read_count_;
  NucleotideReadCount_t mutant_count_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T>
class SNPVariant : public ReadCountVariant {

public:

  SNPVariant(const ContigId_t& contig_id,
             ContigOffset_t contig_offset,
             NucleotideReadCount_t read_count,
             NucleotideReadCount_t mutant_count,
             T reference,
             T mutant)
      : ReadCountVariant(contig_id, contig_offset, read_count, mutant_count), reference_(reference) , mutant_(mutant) {}
  ~SNPVariant() override = default;

  bool equivalent(const Variant& cmp_var) const final;

  const T& reference() const { return reference_; }
  const T& mutant() const { return mutant_; }

private:

  T reference_;
  T mutant_;

  bool applyFilter(const VariantFilter& filter) const final { return filter.applyFilter(*this); }

};

template<class T>
bool SNPVariant<T>::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const SNPVariant<T>*>(&cmp_var);

  if (cmp_snp == nullptr) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<const Variant>>;
class ContigVariant {

public:

  explicit ContigVariant(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  ContigVariant(const ContigVariant&) = default;
  ~ContigVariant() = default;

  ContigVariant& operator=(const ContigVariant&) = default;

  void addVariant(ContigOffset_t contig_offset, std::shared_ptr<const Variant>& variant_ptr);
  const ContigId_t& contigId() const { return contig_id_; }
  size_t variantCount() const { return offset_variant_map_.size(); }

  // Set functions.
  bool isElement(const Variant& variant) const;
  ContigVariant& Union(const ContigVariant& contig_variant) const;
  ContigVariant& Intersection(const ContigVariant& contig_variant) const;
  ContigVariant& Difference(const ContigVariant& contig_variant) const;

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

  GenomeVariant& operator=(const GenomeVariant& genome_variant) = default;

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  void filterVariants(const VariantFilter& filter);

  bool isElement(const Variant& variant) const;
  std::shared_ptr<GenomeVariant> Union(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Intersection(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Difference(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;

private:

  GenomeVariantMap genome_variant_map_;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
