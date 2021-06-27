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
// looked up from the variant database which has been previously parsed and is held in memory attached to a variant object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FrequencyDatabaseRead {

public:

  FrequencyDatabaseRead() = delete;
  ~FrequencyDatabaseRead() = delete;

  // Read super population allele frequencies "AF"
  [[nodiscard]] static std::optional<double> superPopFrequency(const Variant& variant,
                                                               const std::string& super_population);

  // Read super population total number of alleles in called genotypes "AN" (integer).
  [[nodiscard]] static std::optional<int64_t> superPopTotalAlleles( const Variant& variant,
                                                                    const std::string& super_population);

  // Read super population alternate allele count for samples "AC" (integer).
  [[nodiscard]] static std::optional<int64_t> superPopAltAlleles( const Variant& variant,
                                                                  const std::string& super_population);

  // List the super populations supported.
  [[nodiscard]] static const std::vector<std::string>& superPopulations() { return super_populations_; }

  // Read any float field with supplied INFO field code.
  [[nodiscard]] static std::optional<double> infoFloatField(const Variant& variant, const std::string& database_field);

  // Read any integer field with supplied INFO field code.
  [[nodiscard]] static std::optional<int64_t> infoIntegerField(const Variant& variant, const std::string& database_field);

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



  struct SuperPopFieldText {

    std::string gnomad_2_1;
    std::string gnomad_ex_2_1;
    std::string gnomad_3_1;
    std::string gnomadgenome_3_1;
    std::string genome_1000;

  };

  // The Info field identifiers for alternative allele frequency (float).
  inline static const std::map<std::string, SuperPopFieldText> field_text_map_AF_ {

      {SUPER_POP_AFR_, {"AF_afr", "AF_afr", "AF_afr", "gnomad-AF-afr", "AFR_AF"}},  // African
      {SUPER_POP_AMR_, {"AF_amr", "AF_amr", "AF_amr", "gnomad-AF-amr", "AMR_AF"}},   // American
      {SUPER_POP_EAS_, {"AF_eas", "AF_eas", "AF_eas", "gnomad-AF-eas", "EAS_AF"}},   // East Asian
      {SUPER_POP_EUR_, {"AF_nfe", "AF_nfe", "AF_nfe", "gnomad-AF-nfe", "EUR_AF"}},   // European
      {SUPER_POP_SAS_, {"AF", "AF_sas", "AF_sas", "gnomad-AF-sas", "SAS_AF"}},   // South Asian
      {SUPER_POP_ALL_, {"AF", "AF", "AF", "gnomad-AF", "AF"}}   // All Super Populations

  };

  // The Info field identifiers for total number of alleles in called genotypes (integer).
  inline static const std::map<std::string, SuperPopFieldText> field_text_map_AN_ {

      {SUPER_POP_AFR_, {"AN_afr", "AN_afr", "AN_afr", "gnomad-AN-afr", "AFR_AN"}},  // African
      {SUPER_POP_AMR_, {"AN_amr", "AN_amr", "AN_amr", "gnomad-AN-amr", "AMR_AN"}},   // American
      {SUPER_POP_EAS_, {"AN_eas", "AN_eas", "AN_eas", "gnomad-AN-eas", "EAS_AN"}},   // East Asian
      {SUPER_POP_EUR_, {"AN_nfe", "AN_nfe", "AN_nfe", "gnomad-AN-nfe", "EUR_AN"}},   // European
      {SUPER_POP_SAS_, {"AN", "AN_sas", "AN_sas", "gnomad-AN-sas", "SAS_AN"}},   // South Asian
      {SUPER_POP_ALL_, {"AN", "AN", "AN", "gnomad-AN", "AN"}}   // All Super Populations

  };


  // The Info field identifiers for alternate allele count for samples (integer).
  inline static const std::map<std::string, SuperPopFieldText> field_text_map_AC_ {

      {SUPER_POP_AFR_, {"AC_afr", "AC_afr", "AC_afr", "gnomad-AC-afr", "AFR_AC"}},  // African
      {SUPER_POP_AMR_, {"AC_amr", "AC_amr", "AC_amr", "gnomad-AC-amr", "AMR_AC"}},   // American
      {SUPER_POP_EAS_, {"AC_eas", "AC_eas", "AC_eas", "gnomad-AC-eas", "EAS_AC"}},   // East Asian
      {SUPER_POP_EUR_, {"AC_nfe", "AC_nfe", "AC_nfe", "gnomad-AC-nfe", "EUR_AC"}},   // European
      {SUPER_POP_SAS_, {"AC", "AC_sas", "AC_sas", "gnomad-AC-sas", "SAS_AC"}},   // South Asian
      {SUPER_POP_ALL_, {"AC", "AC", "AC", "gnomad-AC", "AC"}}   // All Super Populations

  };

  // Lookup a variant super population frequency code. The field code varies with the (FrequencyDatabaseSource) database source.
  [[nodiscard]] static std::optional<std::string> lookupVariantSuperPopField( const std::map<std::string, SuperPopFieldText>& lookup_map,
                                                                              DataSourceEnum data_source,
                                                                              const std::string& super_population);

};


} // namespace

#endif //KGL_VARIANT_DB_FREQ_H
