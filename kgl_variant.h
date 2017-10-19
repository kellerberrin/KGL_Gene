///
// Created by kellerberrin on 13/10/17.
//

#ifndef KGL_VARIANT_H
#define KGL_VARIANT_H

#include <map>
#include <memory>
#include <vector>
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

  virtual bool applyFilter(const Variant& variant) const = 0;
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
  friend std::ostream& operator<<(std::ostream &os, const Variant& variant) { return variant.output(os); }

private:

  ContigId_t contig_id_;
  ContigOffset_t contig_offset_;

  virtual bool applyFilter(const VariantFilter& filter) const { return filter.applyFilter(*this); }
  virtual bool equivalent(const Variant& cmp_var) const = 0;
  virtual std::ostream& output(std::ostream& os) const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A read count variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ReadCountVariant : public Variant {

public:

  ReadCountVariant(const ContigId_t& contig_id,
                   ContigOffset_t contig_offset,
                   NucleotideReadCount_t read_count,
                   NucleotideReadCount_t mutant_count,
                   NucleotideReadCount_t const count_array[],
                   ContigSize_t count_array_size) : Variant(contig_id, contig_offset),
                                                            read_count_(read_count),
                                                            mutant_count_(mutant_count) {

    for(ContigOffset_t idx = 0; idx < count_array_size; ++idx) {

      count_array_.push_back(count_array[idx]);

    }

  }
  ~ReadCountVariant() override = default;

  NucleotideReadCount_t readCount() const { return read_count_; }
  NucleotideReadCount_t mutantCount() const { return mutant_count_; }
  const std::vector<NucleotideReadCount_t>& countArray() const { return count_array_; }

  double proportion() const { return static_cast<double>(mutant_count_) / static_cast<double>(read_count_); }

private:

  NucleotideReadCount_t read_count_;
  NucleotideReadCount_t mutant_count_;
  std::vector<NucleotideReadCount_t> count_array_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }


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
             NucleotideReadCount_t const count_array[],
             ContigSize_t  count_array_size,
             typename T::NucleotideType reference,
             typename T::NucleotideType mutant)
      : ReadCountVariant(contig_id, contig_offset, read_count, mutant_count, count_array, count_array_size),
        reference_(reference),
        mutant_(mutant) {}
  ~SNPVariant() override = default;

  bool equivalent(const Variant& cmp_var) const final;

  const typename T::NucleotideType& reference() const { return reference_; }
  const typename T::NucleotideType& mutant() const { return mutant_; }

private:

  typename T::NucleotideType reference_;
  typename T::NucleotideType mutant_;

  std::ostream & output(std::ostream &os) const final;

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

template<class T>
std::ostream& SNPVariant<T>::output(std::ostream& os) const
{
  os << contigId() << " " << " " << mutantCount() << "/" << readCount()
     << " " << reference() << (contigOffset() + 1) << mutant() << " [";
  for (auto count : countArray()) {
    os << count << ",";
  }
  return os;
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
  std::shared_ptr<ContigVariant> Union(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Intersection(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Difference(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;

  std::shared_ptr<ContigVariant> filterVariants(const VariantFilter& filter) const;

  friend std::ostream & operator<<(std::ostream &os, const ContigVariant& contig_variant);

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

  explicit GenomeVariant(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  GenomeVariant(const GenomeVariant&) = default;
  ~GenomeVariant() = default;

  GenomeVariant& operator=(const GenomeVariant& genome_variant) = default;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& contig_id) { genome_id_ = contig_id; }

  size_t contigCount() const { return genome_variant_map_.size(); }

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  std::shared_ptr<GenomeVariant> filterVariants(const VariantFilter& filter) const;

  bool isElement(const Variant& variant) const;
  std::shared_ptr<GenomeVariant> Union(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Intersection(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Difference(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;

  friend std::ostream & operator<<(std::ostream &os, const GenomeVariant& contig_variant);

private:

  GenomeId_t genome_id_;
  GenomeVariantMap genome_variant_map_;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
