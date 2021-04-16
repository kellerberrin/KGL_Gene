//
// Created by kellerberrin on 14/4/21.
//

#include "kel_exec_env.h"
#include "kgl_analysis_mutation_gene_ontology.h"


namespace kgl = kellerberrin::genome;


void kgl::OntologyStats::processOntologyStats(const GeneCharacteristic& gene_info, const std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr) {

//  ExecEnv::log().info("OntologyStats::processOntologyStats; ********** Annotation: GO terms : {}, Genes: {}",
//                      ontology_db_ptr->annotation()->getNumGoTerms(), ontology_db_ptr->annotation()->getNumGenes());
//  ExecEnv::log().info("OntologyStats::processOntologyStats; ***********Go Graph Go Vertices: {}, Edges: {}",
//                      ontology_db_ptr->goGraph()->getNumVertices(), ontology_db_ptr->goGraph()->getNumEdges());

  if (ontology_db_ptr->annotation()->hasGene(gene_info.gafId())) {

    ++go_id_count_;
    std::set<std::string> go_set;
    for (auto const& go_term : ontology_db_ptr->annotation()->getGoTermsForGene(gene_info.gafId())) {

      go_set.insert(go_term);

    }

    if (go_set == gene_info.goSet()) {

      ++go_term_count_;

    }
  }

  auto result = malaria_genes_.find(gene_info.gafId());
  if (result != malaria_genes_.end()) {

    auto const& [gaf_id, gene_id] = *result;
    if (gene_info.geneId() != gene_id) {

      ExecEnv::log().error("OntologyStats::processOntologyStats; mismatch MAP GafId: {}, GeneId: {}; INFO GafId: {}, GeneId: {}",
                           gaf_id, gene_id, gene_info.gafId(), gene_info.geneId());

    }

  }


}


void kgl::OntologyStats::writeOntology(std::ostream& out_file, char output_delimiter) const {

  out_file << go_term_count_ << output_delimiter
           << go_id_count_;

}

void kgl::OntologyStats::writeOntologyHeader(std::ostream& out_file, char output_delimiter) const {

  out_file << "GeneGOTerms" << output_delimiter
           << "GOIdGO";

}
