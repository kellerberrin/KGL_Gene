//
// Created by kellerberrin on 24/8/21.
//

#ifndef KGL_ANALYSIS_LITERATURE_H
#define KGL_ANALYSIS_LITERATURE_H


#include "kgl_package_analysis_virtual.h"
#include "kgl_hsgenealogy_parser.h"
#include "kgl_uniprot_parser.h"
#include "kgl_variant_sort_analysis.h"
#include "kga_analysis_mutation_gene.h"
#include "kga_analysis_mutation_gene_allele.h"
#include "kgl_citation_parser.h"
#include "kgl_entrez_parser.h"
#include "kgl_bio_pmid_parser.h"
#include "kgl_pubmed_resource.h"

#include "kga_analysis_literature_gene.h"


namespace kellerberrin::genome::analysis {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LiteratureAnalysis : public VirtualAnalysis {

public:

  LiteratureAnalysis() = default;
  ~LiteratureAnalysis() override = default;

  // The ident must match the ident used in the package XML.
  constexpr static std::string IDENT {"LITERATURE"};
  // Need a polymorphic version to interrogate VirtualAnalysis pointers.
  [[nodiscard]] std::string ident() const override { return IDENT; }
  // Simple creation factory function.
  [[nodiscard]] static std::unique_ptr<VirtualAnalysis> factory() { return std::make_unique<LiteratureAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const ActiveParameterList& named_parameters,
                                         const std::shared_ptr<const AnalysisResources>& resource_ptr) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  GeneLiterature gene_literature_;

  // Various (low memory usage) requested resources
  std::shared_ptr<const GenomeReference> ref_genome_ptr_;
  std::shared_ptr<const UniprotResource> uniprot_nomenclature_ptr_;
  std::shared_ptr<const EnsemblHGNCResource> ensembl_nomenclature_ptr_;
  std::shared_ptr<const EntrezResource> entrez_nomenclature_ptr_;
  std::shared_ptr<const CitationResource> allele_citation_ptr_;
  std::shared_ptr<const PubmedRequester> pubmed_requestor_ptr_;

  // Parameters and output files.
  inline static const std::string GENE_DIRECTORY_NAME_{"Gene"};
  inline static const std::string ANALYSIS_DIRECTORY_NAME_{"Analysis"};
  std::string ident_work_directory_;
  std::string gene_work_directory_;
  std::string analysis_work_directory_;


};


} // namespace

#endif //KGL_ANALYSIS_LITERATURE_H
