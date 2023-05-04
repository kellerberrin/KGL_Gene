//
// Created by kellerberrin on 31/1/21.
//

#include "kgl_analysis_mutation_gene_clinvar.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_variant_filter_info.h"
#include "kgl_variant_filter_db.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::GeneClinvar::writeOutput(  const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                     std::ostream& out_file,
                                     char output_delimiter) const {


  std::string concat_desc{"\""};
  for (auto const& desc : getClinvarDesc()) {

    if (desc != *getClinvarDesc().begin()) {

      concat_desc += CONCAT_TOKEN_;

    }

    concat_desc += desc;

  }
  concat_desc += "\"";
  out_file << concat_desc << output_delimiter;

  out_file << genome_count_ << output_delimiter
           << hom_genome_ << output_delimiter;

  getEthnicity().writeOutput(genome_aux_data, out_file, output_delimiter);

}


void kgl::GeneClinvar::writeHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                    std::ostream& out_file,
                                    char output_delimiter) const {


  out_file << "CLV_Desc" << output_delimiter
           << "CLV_All"  << output_delimiter
           << "CLV_Hom"  << output_delimiter;

  getEthnicity().writeHeader(genome_aux_data, out_file, output_delimiter);

}


void kgl::GeneClinvar:: processClinvar( const GenomeId_t& genome_id,
                                        const ContigId_t& contig_id,
                                        const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                        const std::shared_ptr<const ContigDB>& gene_variants,
                                        const std::shared_ptr<const HsGenomeAux>& genome_aux_data) {

  if (clinvar_contig_->contigId() != contig_id) {

    clinvar_contig_ = GeneClinvar::getClinvarContig(contig_id, clinvar_population_ptr);
    clinvar_contig_ = GeneClinvar::FilterPathogenic(clinvar_contig_);

  }

  processClinvar( genome_id, gene_variants, genome_aux_data);

}


void kgl::GeneClinvar::processClinvar( const GenomeId_t& genome_id,
                                       const std::shared_ptr<const ContigDB>& subject_variants,
                                       const std::shared_ptr<const HsGenomeAux>& genome_aux_data) {

  auto subject_clinvar = clinvar_contig_->findContig(subject_variants);
  auto info_vector = clinvarInfo(subject_clinvar);
  if (subject_clinvar->variantCount() > 0) {

    ++genome_count_;
    updateEthnicity().genomeAnalysis(genome_id, 1, genome_aux_data);

  }

  for (auto const& record : info_vector) {

    // Save any unique descriptions.
    clinvar_desc_.insert(record.clndn);

  }

  auto hom_clinvar = subject_clinvar->viewFilter(HomozygousFilter());
  if (hom_clinvar->variantCount() > 0) {

    ++hom_genome_;

  }

}


std::shared_ptr<const kgl::ContigDB> kgl::GeneClinvar::getClinvarContig( const ContigId_t& contig_id,
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



std::shared_ptr<const kgl::ContigDB> kgl::GeneClinvar::FilterPathogenic(std::shared_ptr<const ContigDB> clinvar_contig) {

  auto substringLambda = [](const std::string& field_str)->bool {

    return Utility::toupper(field_str).find(CLINVAR_PATH_SIGNIF) != std::string::npos;

  };

  return clinvar_contig->viewFilter(InfoFilter<std::string, false>(CLINVAR_CLNSIG_FIELD, substringLambda));

}


std::vector<kgl::ClinvarInfo> kgl::GeneClinvar::clinvarInfo(const std::shared_ptr<const ContigDB>& clinvar_contig_ptr) {

  std::vector<ClinvarInfo> clinvarVector;

  for (auto const& [offset, offset_ptr] : clinvar_contig_ptr->getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      ClinvarInfo clinvar_record;

      auto field_opt = InfoEvidenceAnalysis::getInfoData(*variant_ptr, CLINVAR_CLNDN_FIELD);
      if (field_opt) {

        std::vector<std::string> desc_vector = InfoEvidenceAnalysis::varianttoStrings(field_opt.value());
        if (not desc_vector.empty()) {

          clinvar_record.clndn = desc_vector.front();

        }

      }

      clinvar_record.variant_ptr = variant_ptr;

      clinvarVector.push_back(clinvar_record);

    }

  }

  return clinvarVector;

}
