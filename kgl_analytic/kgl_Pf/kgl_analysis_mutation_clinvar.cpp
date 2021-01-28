//
// Created by kellerberrin on 28/1/21.
//

#include "kgl_analysis_mutation_clinvar.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"



namespace kgl = kellerberrin::genome;




std::shared_ptr<const kgl::ContigDB> kgl::AnalyzeClinvar::getClinvarContig( const ContigId_t& contig_id,
                                                                            const std::shared_ptr<const PopulationDB>& clinvar_population_ptr) {

  std::shared_ptr<ContigDB> null_contig_ptr(std::make_shared<ContigDB>(contig_id));

  if (clinvar_population_ptr->getMap().size() != 1) {

    ExecEnv::log().error("AnalyzeClinvar::getClinvarContig; expected clinvar population to have 1 genome, actual size: {}",
                         clinvar_population_ptr->getMap().size());

    return null_contig_ptr;

  }

  auto [genomne_id, genome_ptr] = *(clinvar_population_ptr->getMap().begin());

  auto contig_opt = genome_ptr->getContig(contig_id);

  if (not contig_opt) {

    return null_contig_ptr;

  }

  auto clinvar_contig_ptr = contig_opt.value();

  return clinvar_contig_ptr;

}



std::shared_ptr<const kgl::ContigDB> kgl::AnalyzeClinvar::FilterPathogenic(std::shared_ptr<const ContigDB> clinvar_contig) {

  return clinvar_contig->filterVariants(InfoSubStringFilter(CLINVAR_CLNSIG_FIELD, CLINVAR_PATH_SIGNIF));

}


std::shared_ptr<const kgl::ContigDB> kgl::AnalyzeClinvar::findClinvar( const std::shared_ptr<const ContigDB>& subject_contig,
                                                                       const std::shared_ptr<const ContigDB>& clinvar_contig_ptr) {

  auto found_contig = clinvar_contig_ptr->findContig(subject_contig);

  return found_contig;


}


std::vector<kgl::ClinvarInfo> kgl::AnalyzeClinvar::clinvarInfo(const std::shared_ptr<const ContigDB>& clinvar_contig_ptr) {

  std::vector<ClinvarInfo> clinvarVector;

  for (auto const& [offset, offset_ptr] : clinvar_contig_ptr->getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      ClinvarInfo clinvar_record;

      auto field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, CLINVAR_CLNSIG_FIELD);
      if (field_opt) {

        std::vector<std::string> signif_vector = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());
        if (not signif_vector.empty()) {

          clinvar_record.clnsig = signif_vector.front();

        }

      }

      field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, CLINVAR_CLNDN_FIELD);
      if (field_opt) {

        std::vector<std::string> desc_vector = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());
        if (not desc_vector.empty()) {

          clinvar_record.clndn = desc_vector.front();

        }

      }

      field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, CLINVAR_CLNDISDB_FIELD);
      if (field_opt) {

        std::vector<std::string> db_vector = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());
        if (not db_vector.empty()) {

          clinvar_record.clnisdb = db_vector.front();

        }

      }

      field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, CLINVAR_RS_FIELD);
      if (field_opt) {

        std::vector<std::string> rs_vector = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());
        if (not rs_vector.empty()) {

          clinvar_record.clnisdb = CLINVAR_RS_FIELD + rs_vector.front();

        }

      }

      clinvar_record.variant_ptr = variant_ptr;

      clinvarVector.push_back(clinvar_record);

    }

  }

  return clinvarVector;

}



std::vector<std::string> kgl::AnalyzeClinvar::clinvarVectorDesc(const std::vector<ClinvarInfo>& clinvar_vector) {

  std::set<std::string> desc_set;
  std::vector<std::string> desc_vector;

  for (auto const& clinvar_record: clinvar_vector) {

    desc_set.insert(clinvar_record.clndn);

  }

  for (auto const& desc: desc_set) {

    desc_vector.push_back(desc);

  }

  return desc_vector;

}