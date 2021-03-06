//
// Created by kellerberrin on 4/3/21.
//

#include "kgl_analysis_mutation_gene_stats.h"


namespace kgl = kellerberrin::genome;



// Perform the genetic analysis per iteration.
bool kgl::GeneCharacteristic::geneDefinition( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                              const GenomeId_t& genome_id,
                                              const std::string& name,
                                              const std::string& gaf_ident)
{

  std::vector<std::string> description_vec;
  gene_ptr->getAttributes().getDescription(description_vec);
  std::string description_str = description_vec.empty() ? "" : description_vec.front();

  std::vector<std::string> gene_biotype_vec;
  gene_ptr->getAttributes().getGeneBioType(gene_biotype_vec);
  std::string biotype_str = gene_biotype_vec.empty() ? "" : gene_biotype_vec.front();

  genome_ = genome_id;
  contig_ = gene_ptr->contig()->contigId();
  gene_ptr_ = gene_ptr;
  gene_id_ = gene_ptr->id();
  gene_name_ = name;
  description_ = description_str;
  biotype_ = biotype_str;
  valid_protein_ = ContigReference::verifyGene(gene_ptr);
  gaf_id_ = gaf_ident;
  gene_begin_ = gene_ptr->sequence().begin();
  gene_end_ = gene_ptr->sequence().end();
  gene_span_ = gene_ptr->sequence().length();
  strand_ = gene_ptr->sequence().strandText();

  const std::shared_ptr<const CodingSequenceArray> sequence_array = GeneFeature::getCodingSequences(gene_ptr);

  sequences_ = sequence_array->size();
  if (not sequence_array->getMap().empty()) {

    auto [sequence_name, seq_ptr] = *(sequence_array->getMap().begin());
    seq_name_ = sequence_name;
    nucleotides_ = seq_ptr->codingNucleotides();
    exons_ = seq_ptr->exons();

  }
  auto const& attribute_map = gene_ptr->getAttributes().getMap();
  attribute_size_ = attribute_map.size();

  std::string concat_attributes{"\""};
  for (auto const& [key, attrib] : attribute_map) {

    if ( key != attribute_map.begin()->first
         or attrib != attribute_map.begin()->second) {

      concat_attributes += CONCAT_TOKEN_;

    }

    concat_attributes += key;
    concat_attributes += ":";
    concat_attributes += attrib;

  }

  concat_attributes += "\"";
  attributes_ = concat_attributes;

  return true;

}


