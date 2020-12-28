//
// Created by kellerberrin on 20/11/20.
//

#ifndef KGL_VARIANT_DB_FREQ_H
#define KGL_VARIANT_DB_FREQ_H

#include "kgl_variant.h"
#include "kgl_variant_db_population.h"

namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convenience object.
//
// Different variant databases such as Gnomad 2.1, Gnomad 3.1 and 1000 Genomes use slightly different field codes
// to access variant data such as super population variant frequencies. This object hides these differences.
// To access super population EAS frequencies "EAS" is passed to this object and the correct database code
// is looked up. So for Gnomad 2.1 "EAS" is translated internally to "AF_eas" and the variant frequency is
// looked up from the variant database which has been previously parsed and is held in memory.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FrequencyDatabaseRead {

public:

  FrequencyDatabaseRead() = delete;
  ~FrequencyDatabaseRead() = delete;

  // Read super population frequencies
  [[nodiscard]] static std::optional<double> processSuperPopField(const Variant& variant,
                                                                  DataSourceEnum data_source,
                                                                  const std::string& super_population);
  // List the super populations supported.
  [[nodiscard]] static const std::vector<std::string>& superPopulations() { return super_populations_; }

  // Read any frequency with supplied INFO field code.
  [[nodiscard]] static std::optional<double> processFloatField(const Variant& variant, const std::string& database_field);

  // Valid super population codes.
  constexpr static const char* SUPER_POP_AFR_{"AFR"} ;  // African
  constexpr static const char* SUPER_POP_AMR_{"AMR"};  // American
  constexpr static const char* SUPER_POP_EAS_{"EAS"};  // East Asian
  constexpr static const char* SUPER_POP_EUR_{"EUR"};  // European
  constexpr static const char* SUPER_POP_SAS_{"SAS"};  // South Asian
  constexpr static const char* SUPER_POP_ALL_{"ALL"};  // All populations.

private:

  inline static std::vector<std::string> super_populations_ = { SUPER_POP_AFR_,
                                                                SUPER_POP_AMR_,
                                                                SUPER_POP_EAS_,
                                                                SUPER_POP_EUR_,
                                                                SUPER_POP_SAS_,
                                                                SUPER_POP_ALL_ };

  // Lookup a variant super population frequency code. The field code varies with the (FrequencyDatabaseSource) database source.
  [[nodiscard]] std::string lookupVariantSuperPopField(const std::string& super_population) const;

  // Lookup a variant super population frequency code. The field code varies with the (FrequencyDatabaseSource) database source.
  [[nodiscard]] static std::optional<std::string> lookupVariantSuperPopField( DataSourceEnum data_source,
                                                                              const std::string& super_population);

  // Use a super population code to lookup a corresponding AF field.
  [[nodiscard]] static std::string lookupGnomad_2_1_Field(const std::string& super_population);

  using SuperPopPair = std::pair<const char*, const char*>;
  // The Info field identifiers for allele frequency for Gnomad 2.1.
  constexpr static const SuperPopPair SUPER_POP_AFR_GNOMAD_2_1 {SUPER_POP_AFR_, "AF_afr"} ;  // African
  constexpr static const SuperPopPair SUPER_POP_AMR_GNOMAD_2_1 {SUPER_POP_AMR_, "AF_amr"};  // American
  constexpr static const SuperPopPair SUPER_POP_EAS_GNOMAD_2_1 {SUPER_POP_EAS_, "AF_eas"};  // East Asian
  constexpr static const SuperPopPair SUPER_POP_EUR_GNOMAD_2_1 {SUPER_POP_EUR_, "AF_nfe"};  // European
  constexpr static const SuperPopPair SUPER_POP_SAS_GNOMAD_2_1 {SUPER_POP_SAS_, "AF"};  // South Asian
  constexpr static const SuperPopPair SUPER_POP_ALL_GNOMAD_2_1 {SUPER_POP_ALL_, "AF"};  // All Super Populations

  // Use a super population code to lookup a corresponding AF field.
  [[nodiscard]] static std::string lookupGnomad_ex_2_1_Field(const std::string& super_population);

  // The Info field identifiers for allele frequency for Gnomes 2.1 exomes.
  constexpr static const SuperPopPair SUPER_POP_AFR_GNOMAD_EX_2_1 {SUPER_POP_AFR_, "AF_afr"} ;  // African
  constexpr static const SuperPopPair SUPER_POP_AMR_GNOMAD_EX_2_1 {SUPER_POP_AMR_, "AF_amr"};  // American
  constexpr static const SuperPopPair SUPER_POP_EAS_GNOMAD_EX_2_1 {SUPER_POP_EAS_, "AF_eas"};  // East Asian
  constexpr static const SuperPopPair SUPER_POP_EUR_GNOMAD_EX_2_1 {SUPER_POP_EUR_, "AF_nfe"};  // European
  constexpr static const SuperPopPair SUPER_POP_SAS_GNOMAD_EX_2_1 {SUPER_POP_SAS_, "AF_sas"};  // South Asian
  constexpr static const SuperPopPair SUPER_POP_ALL_GNOMAD_EX_2_1 {SUPER_POP_ALL_, "AF"};  // All Super Populations

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
  [[nodiscard]] static std::string lookupGnomadGenome_3_1_Field(const std::string& super_population);

  // The Info field identifiers for allele frequency for Gnomad 3.1
  constexpr static const SuperPopPair SUPER_POP_AFR_GNOMADGENOME_3_1 {SUPER_POP_AFR_, "gnomad-AF-afr"} ;  // African
  constexpr static const SuperPopPair SUPER_POP_AMR_GNOMADGENOME_3_1 {SUPER_POP_AMR_, "gnomad-AF-amr"};  // American
  constexpr static const SuperPopPair SUPER_POP_EAS_GNOMADGENOME_3_1 {SUPER_POP_EAS_, "gnomad-AF-eas"};  // East Asian
  constexpr static const SuperPopPair SUPER_POP_EUR_GNOMADGENOME_3_1 {SUPER_POP_EUR_, "gnomad-AF-nfe"};  // European
  constexpr static const SuperPopPair SUPER_POP_SAS_GNOMADGENOME_3_1 {SUPER_POP_SAS_, "gnomad-AF-sas"};  // South Asian
  constexpr static const SuperPopPair SUPER_POP_ALL_GNOMADGENOME_3_1 {SUPER_POP_ALL_, "gnomad-AF"};  // All Super Populations


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

#endif //KGL_VARIANT_DB_FREQ_H
