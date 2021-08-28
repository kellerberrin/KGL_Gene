//
// Created by kellerberrin on 27/8/21.
//

#ifndef KGL_ANALYSIS_LITERATURE_GENE_H
#define KGL_ANALYSIS_LITERATURE_GENE_H

#include "kgl_uniprot_parser.h"
#include "kgl_citation_parser.h"
#include "kgl_entrez_parser.h"
#include "kgl_bio_pmid_parser.h"
#include "kgl_pubmed_resource.h"
#include "kgl_mutation/kgl_analysis_mutation_gene_stats.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GeneLiterature {

public:

  GeneLiterature() = default;
  ~GeneLiterature() = default;

  // Create the gene vector.
  void defineGenes( const std::shared_ptr<const GenomeReference>& genome_ptr,
                    const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                    const std::shared_ptr<const EntrezResource>& entrez_nomenclature_ptr);

  void updatePMIDStatistics( const std::set<std::string>& disease_pmid_set,
                             const std::shared_ptr<const BioPMIDFileData>& bio_pmid_ptr);

  void outputGenePmid(const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr,
                      const std::string& literature_directory,
                      size_t pmid_count) const;

  void outputPmidGene(const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr,
                      const std::string& literature_directory,
                      size_t max_genes,
                      size_t min_genes,
                      size_t min_citations) const;

  static bool filterPublication(const PubMedPublicationSummary& publication);

  static void writeGenePublications( std::ostream& out_file,
                                     const GeneCharacteristic& gene,
                                     const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr);
private:

//  std::shared_ptr<const PubmedRequester> pubmed_requestor_ptr_;

  // All the genes found in the genome.
  std::vector<GeneCharacteristic> gene_vector_;



};



} // namespace.


#endif // KGL_ANALYSIS_LITERATURE_GENE_H
