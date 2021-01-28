//
// Created by kellerberrin on 28/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_CLINVAR_H
#define KGL_ANALYSIS_MUTATION_CLINVAR_H

#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct ClinvarInfo {

  std::string rsid;          // The variant rs identifier from dbSNP
  std::string clnsig;        // Textual clinical significance.
  std::string clndn;         // Clinical description.
  std::string clnisdb;       // Database information about the condition
  std::shared_ptr<const Variant> variant_ptr;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AnalyzeClinvar {

public:

  AnalyzeClinvar() = delete;
  ~AnalyzeClinvar() = delete;

  static std::shared_ptr<const ContigDB> getClinvarContig(const ContigId_t& contig_id,
                                                          const std::shared_ptr<const PopulationDB>& clinvar_population_ptr);

  static std::shared_ptr<const ContigDB> FilterPathogenic(std::shared_ptr<const ContigDB> clinvar_contig);

  static std::shared_ptr<const ContigDB> findClinvar(const std::shared_ptr<const ContigDB>& subject_contig,
                                                     const std::shared_ptr<const ContigDB>& clinvar_contig_ptr);

  static std::vector<ClinvarInfo> clinvarInfo(const std::shared_ptr<const ContigDB>& clinvar_contig_ptr);

  static std::vector<std::string> clinvarVectorDesc(const std::vector<ClinvarInfo>& clinvar_vector);


private:

  // Clinvar fields.
  constexpr static const char* CLINVAR_RS_FIELD = "RS";
  constexpr static const char* CLINVAR_CLNDN_FIELD = "CLNDN";
  constexpr static const char* CLINVAR_CLNSIG_FIELD = "CLNSIG";
  constexpr static const char* CLINVAR_CLNDISDB_FIELD = "CLNDISDB";
  constexpr static const char* CLINVAR_PATH_SIGNIF = "PATHOGENIC";
  constexpr static const char* CLINVAR_RISK_SIGNIF = "RISK";
  constexpr static const char* CLINVAR_DESC_CONCAT = "&";


};




} // namespace



#endif // KGL_ANALYSIS_MUTATION_CLINVAR_H
