//
// Created by kellerberrin on 4/3/21.
//

#include "kga_analysis_mutation_gene_stats.h"


namespace kga = kellerberrin::genome::analysis;



// Perform the genetic analysis per iteration.
bool kga::GeneCharacteristic::geneDefinition( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                              const GenomeId_t& genome_id,
                                              const std::string& name,
                                              const std::string& hgnc_id,
                                              const std::vector<std::string>& ensembl_ids,
                                              const std::string& uniprot_id,
                                              const std::string& entrez_id)
{

  std::vector<std::string> description_vec;
  gene_ptr->getAttributes().getDescription(description_vec);
  std::string description_str = description_vec.empty() ? "" : description_vec.front();

  std::vector<std::string> gene_biotype_vec;
  gene_ptr->getAttributes().getGeneBioType(gene_biotype_vec);
  std::string biotype_str = gene_biotype_vec.empty() ? "" : gene_biotype_vec.front();

  genome_ = genome_id;
  contig_ = gene_ptr->contig_ref_ptr()->contigId();
  gene_ptr_ = gene_ptr;
  gene_id_ = gene_ptr->id();
  symbol_id_ = name;
  description_ = description_str;
  biotype_ = biotype_str;
  valid_protein_ = ContigReference::verifyGene(gene_ptr);
  gene_begin_ = gene_ptr->sequence().begin();
  gene_end_ = gene_ptr->sequence().end();
  gene_span_ = gene_ptr->sequence().length();
  strand_ = gene_ptr->sequence().strandText();

  auto transcript_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);

  sequences_ = transcript_array_ptr->size();
  if (not transcript_array_ptr->getMap().empty()) {

    auto& [transcript_id, transcript_ptr] = *(transcript_array_ptr->getMap().begin());
    seq_name_ = transcript_id;
    nucleotides_ = transcript_ptr->codingNucleotides();
    exons_ = transcript_ptr->getFeatureMap().size();

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
    concat_attributes += "-";
    concat_attributes += attrib;

  }

  concat_attributes += "\"";
  attributes_ = concat_attributes;

  uniprotKB_id_ = uniprot_id;
  ensembl_ids_ = ensembl_ids;
  HGNC_id_ = hgnc_id;
  entrez_id_ = entrez_id;

  return true;

}


