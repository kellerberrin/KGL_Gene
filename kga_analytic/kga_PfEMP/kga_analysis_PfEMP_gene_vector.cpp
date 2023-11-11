//
// Created by kellerberrin on 6/04/23.
//


#include "kga_analysis_PfEMP.h"


namespace kgl = kellerberrin::genome;


kgl::GeneVector kgl::PfEMPAnalysis::getGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr
    , const std::string& desc_uc_text) const {

  GeneVector gene_vector;

  for (auto const& [contig_ident, contig_ptr] : genome_ptr->getMap()) {

    for (auto const &[gene_offset, gene_ptr]: contig_ptr->getGeneMap()) {

      auto description_vector = gene_ptr->getAttributes().getDescription();
      for (auto const& description : description_vector) {

        if (Utility::toupper(description).find(desc_uc_text) != std::string::npos) {

          gene_vector.push_back(gene_ptr);

        }

      } // Gene

    } // Contig

  } // Genome

  return gene_vector;

}


kgl::GeneVector kgl::PfEMPAnalysis::getncRNAGeneVector( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                                        const std::string& desc_uc_text,
                                                        size_t max_size) const {

  GeneVector gene_vector;
  size_t contig_count{0};

  for (auto const& [contig_ident, contig_ptr] : genome_ptr->getMap()) {

    ++contig_count;

    for (auto const &[gene_offset, gene_ptr]: contig_ptr->getGeneMap()) {

      auto description = gene_ptr->descriptionText();
      if (GeneFeature::ncRNACoding(gene_ptr)
          and  (gene_ptr->sequence().length() <= max_size or max_size == 0)
          and (Utility::toupper(description).find(desc_uc_text) != std::string::npos or desc_uc_text.empty())) {

        gene_vector.push_back(gene_ptr);

      }

    } // Gene

  } // Contig

  return gene_vector;

}

