//
// Created by kellerberrin on 20/11/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_FREQDB_H
#define KGL_ANALYSIS_MUTATION_INBREED_FREQDB_H

#include "kgl_variant.h"

namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Different variant databases such as Gnomad 2.1, Gnomad 3.1 and 1000 Genomes use slightly different field codes
// to access variant data such as super population variant frequencies. This object hides these differences.
// To access super population EAS frequencies "EAS" is passed to this object and the correct database code
// is looked up. So for Gnomad 2.1 "EAS" is translated internally to "AF_eas" and the variant frequency is
// looked up from the variant database which has been previously parsed and is held in memory.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class FrequencyDatabaseSource { GNOMAD2_1, GNOMAD3_1, GENOMES_1000 };

class FrequencyDatabaseRead {

public:

  explicit FrequencyDatabaseRead(FrequencyDatabaseSource source) : source_(source) {}
  ~FrequencyDatabaseRead() = default;

  // Get a scalar floating Info field.
  [[nodiscard]] std::optional<double> processFloatField( const Variant& variant, const std::string& super_population) const;


  // List the super populations supported.
  [[nodiscard]] static const std::vector<std::string>& superPopulations() { return super_populations_; }

  // Valid super population codes.
  constexpr static const char* SUPER_POP_AFR_ {"AFR"} ;  // African
  constexpr static const char* SUPER_POP_AMR_ {"AMR"};  // American
  constexpr static const char* SUPER_POP_EAS_ {"EAS"};  // East Asian
  constexpr static const char* SUPER_POP_EUR_{"EUR"};  // European
  constexpr static const char* SUPER_POP_SAS_{"SAS"};  // South Asian
  constexpr static const char* SUPER_POP_ALL_{"ALL"};  // All populations.

private:

  FrequencyDatabaseSource source_;


  inline static std::vector<std::string> super_populations_ = { SUPER_POP_AFR_,
                                                                SUPER_POP_AMR_,
                                                                SUPER_POP_EAS_,
                                                                SUPER_POP_EUR_,
                                                                SUPER_POP_SAS_,
                                                                SUPER_POP_ALL_ };

  // Lookup a variant super population frequency code. The field code varies with the (FrequencyDatabaseSource) database source.
  [[nodiscard]] std::string lookupVariantSuperPopField(const std::string& super_population) const;

  // Use a super population code to lookup a corresponding AF field.
  [[nodiscard]] static std::string lookupGnomad_2_1_Field(const std::string& super_population);

  using SuperPopPair = std::pair<const char*, const char*>;
  // The Info field identifiers for allele frequency for Gnomad 2.1 and 3.0
  constexpr static const SuperPopPair SUPER_POP_AFR_GNOMAD_2_1 {SUPER_POP_AFR_, "AF_afr"} ;  // African
  constexpr static const SuperPopPair SUPER_POP_AMR_GNOMAD_2_1 {SUPER_POP_AMR_, "AF_amr"};  // American
  constexpr static const SuperPopPair SUPER_POP_EAS_GNOMAD_2_1 {SUPER_POP_EAS_, "AF_eas"};  // East Asian
  constexpr static const SuperPopPair SUPER_POP_EUR_GNOMAD_2_1 {SUPER_POP_EUR_, "AF_nfe"};  // European
  constexpr static const SuperPopPair SUPER_POP_SAS_GNOMAD_2_1 {SUPER_POP_SAS_, "AF"};  // South Asian
  constexpr static const SuperPopPair SUPER_POP_ALL_GNOMAD_2_1 {SUPER_POP_ALL_, "AF"};  // All Super Populations

  // Use a super population code to lookup a corresponding AF field.
  [[nodiscard]] static std::string lookupGnomad_3_1_Field(const std::string& super_population);

  // The Info field identifiers for allele frequency for Gnomad 3.1
  constexpr static const SuperPopPair SUPER_POP_AFR_GNOMAD_3_1 {SUPER_POP_AFR_, "AF-afr"} ;  // African
  constexpr static const SuperPopPair SUPER_POP_AMR_GNOMAD_3_1 {SUPER_POP_AMR_, "AF-amr"};  // American
  constexpr static const SuperPopPair SUPER_POP_EAS_GNOMAD_3_1 {SUPER_POP_EAS_, "AF-eas"};  // East Asian
  constexpr static const SuperPopPair SUPER_POP_EUR_GNOMAD_3_1 {SUPER_POP_EUR_, "AF-nfe"};  // European
  constexpr static const SuperPopPair SUPER_POP_SAS_GNOMAD_3_1 {SUPER_POP_SAS_, "AF-sas"};  // South Asian
  constexpr static const SuperPopPair SUPER_POP_ALL_GNOMAD_3_1 {SUPER_POP_ALL_, "AF"};  // All Super Populations

  // Use a super population code to lookup a corresponding AF field.
  [[nodiscard]] static std::string lookup_1000_Field(const std::string& super_population);

  // The Info field identifiers for allele frequency for the 1000 genomes project.
  constexpr static const SuperPopPair SUPER_POP_AFR_1000 {SUPER_POP_AFR_, "AFR_AF"} ;  // African
  constexpr static const SuperPopPair SUPER_POP_AMR_1000 {SUPER_POP_AMR_, "AMR_AF"};  // American
  constexpr static const SuperPopPair SUPER_POP_EAS_1000 {SUPER_POP_EAS_, "EAS_AF"};  // East Asian
  constexpr static const SuperPopPair SUPER_POP_EUR_1000 {SUPER_POP_EUR_, "EUR_AF"};  // European
  constexpr static const SuperPopPair SUPER_POP_SAS_1000 {SUPER_POP_SAS_, "SAS_AF"};  // South Asian
  constexpr static const SuperPopPair SUPER_POP_ALL_1000 {SUPER_POP_ALL_, "AF"};  // All Super Populations

};


} // namespace

#endif //KGL_ANALYSIS_MUTATION_INBREED_FREQDB_H
