//
// Created by kellerberrin on 4/3/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_STATS_H
#define KGL_ANALYSIS_MUTATION_GENE_STATS_H

#include "kgl_genome_genome.h"
#include "kgl_hsgenealogy_parser.h"
#include "kgl_variant_db_population.h"
#include "kgl_pubmed_resource.h"



namespace kellerberrin::genome {   //  organization::project level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class GeneCharacteristic {

public:

  GeneCharacteristic() = default;
  ~GeneCharacteristic() = default;

  GeneCharacteristic(const GeneCharacteristic&) = default;
  GeneCharacteristic& operator=(const GeneCharacteristic&) = default;

  void writeGene(std::ostream& out_file, char output_delimiter) const;
  void writeGeneHeader(std::ostream& out_file, char output_delimiter) const;
  bool geneDefinition( const std::shared_ptr<const GeneFeature>& gene_ptr,
                       const GenomeId_t& genome_id,
                       const std::string& name,
                       const std::string& hgnc_id,
                       const std::vector<std::string>& ensembl_ids,
                       const std::string& uniprot_id,
                       const std::string& entrez_id);

  [[nodiscard]] const ContigId_t& contigId() const { return contig_; }
  [[nodiscard]] const std::shared_ptr<const GeneFeature>& genePtr() const { return gene_ptr_; }
  [[nodiscard]] ContigOffset_t geneBegin() const { return gene_begin_; }
  [[nodiscard]] ContigOffset_t geneEnd() const { return gene_end_; }
  [[nodiscard]] ContigSize_t nucleotides() const { return nucleotides_; }
  [[nodiscard]] ContigSize_t geneSpan() const { return gene_span_; }
  [[nodiscard]] const std::vector<std::string>& ensemblIds() const { return ensembl_ids_; }
  [[nodiscard]] const std::string& gafId() const { return uniprotKB_id_; }
  [[nodiscard]] const std::string& symbolId() const { return symbol_id_; }
  [[nodiscard]] const std::string& entrezId() const { return entrez_id_; }
  [[nodiscard]] const std::set<std::string>& goSet() const { return GO_set_; }
  [[nodiscard]] const std::set<std::string>& diseasePublications() const { return disease_cites_; }


  void  writeGenePublications( std::ostream& out_file,
                               const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr) const;
  void update_pmid(size_t all_pmid, std::set<std::string> disease_pmids) { citations_ = all_pmid; disease_cites_ = std::move(disease_pmids); }

private:

  constexpr static const char* CONCAT_TOKEN_ = "&";

  GenomeId_t genome_;
  ContigId_t contig_;
  FeatureIdent_t gene_id_;
  std::string symbol_id_;
  std::string description_;
  std::string biotype_;
  bool valid_protein_{false};
  std::string uniprotKB_id_;
  std::set<std::string> GO_set_;
  std::string entrez_id_;
  std::string HGNC_id_;
  std::vector<std::string> ensembl_ids_;
  size_t citations_{0};
  std::set<std::string> disease_cites_;
  ContigOffset_t gene_begin_{0};
  ContigOffset_t gene_end_{0};
  ContigSize_t gene_span_{0};
  ContigSize_t exons_{0};
  ContigSize_t nucleotides_{0};
  std::string strand_;
  size_t sequences_{0};
  std::string seq_name_;
  size_t attribute_size_{0};
  std::string attributes_;
  std::shared_ptr<const GeneFeature> gene_ptr_;

};






} // namespace.

#endif //KGL_ANALYSIS_MUTATION_GENE_STATS_H
