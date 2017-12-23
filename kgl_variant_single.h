//
// Created by kellerberrin on 23/12/17.
//

#ifndef KGL_VARIANT_SINGLE_H
#define KGL_VARIANT_SINGLE_H


#include "kgl_variant_subordinate.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPVariant : public SubordinateSNP {

public:

  SNPVariant(const std::string& variant_source,
             std::shared_ptr<const ContigFeatures> contig_ptr,
             ContigOffset_t contig_offset,
             Phred_t quality,
             std::shared_ptr<const VariantEvidence> evidence_ptr,
             DNA5::Alphabet reference,
             ExtendDNA5::Alphabet mutant) : SubordinateSNP(variant_source,
                                                           contig_ptr,
                                                           contig_offset,
                                                           quality,
                                                           evidence_ptr,
                                                           reference,
                                                           mutant) {}

  SNPVariant(const SNPVariant& variant) = default;
  ~SNPVariant() override = default;

  virtual std::shared_ptr<Variant> clone() const override { return std::shared_ptr<SNPVariant>(std::make_shared<SNPVariant>(*this)); }

  VariantType variantType() const override { return VariantType::SNP; }

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;


  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple deletion variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DeleteVariant : public SubordinateSNP {

public:

  DeleteVariant(const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                ContigOffset_t contig_offset,
                Phred_t quality,
                std::shared_ptr<const VariantEvidence> evidence_ptr,
                DNA5::Alphabet reference,
                ExtendDNA5::Alphabet mutant) : SubordinateSNP(variant_source,
                                                              contig_ptr,
                                                              contig_offset,
                                                              quality,
                                                              evidence_ptr,
                                                              reference,
                                                              mutant) {}

  DeleteVariant(const DeleteVariant& variant) = default;
  ~DeleteVariant() override = default;

  // Use this - not the copy contructor.
  virtual std::shared_ptr<Variant> clone() const override { return std::shared_ptr<DeleteVariant>(std::make_shared<DeleteVariant>(*this)); }

  VariantType variantType() const override { return VariantType::DELETE; }

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;


  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple insertion variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InsertVariant : public SubordinateSNP {

public:

  InsertVariant(const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                ContigOffset_t contig_offset,
                Phred_t quality,
                std::shared_ptr<const VariantEvidence> evidence_ptr,
                DNA5::Alphabet reference,
                ExtendDNA5::Alphabet mutant) : SubordinateSNP(variant_source,
                                                              contig_ptr,
                                                              contig_offset,
                                                              quality,
                                                              evidence_ptr,
                                                              reference,
                                                              mutant) {}

  InsertVariant(const InsertVariant& variant) = default;
  ~InsertVariant() override = default;

  virtual std::shared_ptr<Variant> clone() const override { return std::shared_ptr<InsertVariant>(std::make_shared<InsertVariant>(*this)); }

  VariantType variantType() const override { return VariantType::INSERT; }

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;


  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_SINGLE_H
