//
// Created by kellerberrin on 30/11/20.
//

#ifndef KGL_ANALYSIS_INBREED_OUTPUT_H
#define KGL_ANALYSIS_INBREED_OUTPUT_H


#include "kgl_hsgenealogy_parser.h"
#include "kga_analysis_inbreed_args.h"


namespace kellerberrin::genome::analysis {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures to hold threadpool calculation results.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct LocusResults {

  GenomeId_t genome;
  size_t major_hetero_count{0};   // 1 minor allele only.
  double major_hetero_freq{0.0};
  size_t minor_hetero_count{0};   // 2 different minor alleles.
  double minor_hetero_freq{0.0};
  size_t minor_homo_count{0};     // 2 identical minor alleles.
  double minor_homo_freq{0.0};
  size_t major_homo_count{0};     // 2 identical major alleles (generally not recorded).
  double major_homo_freq{0.0};
  size_t total_allele_count{0};  // All alleles.
  double inbred_allele_sum{0.0}; // inbreeding coefficient

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Indexes the above results by analyzed genome (output rows).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ResultsMap = std::map<GenomeId_t, LocusResults>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A named column of results, rows are indexed by GenomeId_t in the ResultsMap.
// Typically this will be the results for a contig or part of a contig_ref_ptr.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InbreedingResultColumn {

public:

  explicit InbreedingResultColumn(std::string column_ident, ResultsMap results) : column_identifier_(std::move(column_ident)),
                                                                                  results_(std::move(results)) {}
  ~InbreedingResultColumn() = default;

  [[nodiscard]] const std::string& columnIdent() const { return column_identifier_; }
  [[nodiscard]] const ResultsMap& results() const { return results_; }

  // Helper function generates a column ident based on a contig_ref_ptr id and the region analyzed.
  static std::string generateIdent(const ContigId_t& contig_id, ContigOffset_t lower, ContigOffset_t upper);

private:

  std::string column_identifier_;
  ResultsMap results_;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The final layer of results. Each object corresponds to an InbreedingParameters record defined in the XML parameter file.
// The object contains multiple output result columns.
// Each object produces a separate .csv output file with the named columns defined above.
// Note that the output file is defined in the InbreedingParameters record.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InbreedParamOutput {

public:

  explicit InbreedParamOutput(const InbreedingParameters& parameters) : parameters_(parameters) {}
  ~InbreedParamOutput() = default;

  void addColumn(const InbreedingResultColumn& column) { column_results_.push_back(column); }
  [[nodiscard]] const std::vector<InbreedingResultColumn>& getColumns() const { return column_results_; }
  [[nodiscard]] const InbreedingParameters& getParameters() const { return parameters_; }
  // Non-const version
  [[nodiscard]] InbreedingParameters& getParameters() { return parameters_; }
  // Check that each column has the same genome structure.
  [[nodiscard]] bool verifyResults () const;

private:

  InbreedingParameters parameters_;
  std::vector<InbreedingResultColumn> column_results_;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Format and output results to a csv file.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////




class InbreedingOutput  {

public:

  InbreedingOutput() = delete;
  ~InbreedingOutput() = delete;



  // Write the analysis results to a CSV file.
  static bool writePedResults( const InbreedParamOutput& output_results,
                               const HsGenomeGenealogyData& ped_data,
                               const std::string& file_path);

  static bool writeNoPedResults(const InbreedParamOutput& output_results, const std::string& file_path);
  static bool writeSynthetic(const InbreedParamOutput& output_results, const std::string& file_path);


private:

  constexpr static const char DELIMITER_ = ',';
  constexpr static const char* FILE_EXT_ = ".csv";


};



} // namespace


#endif //KGL_ANALYSIS_MUTATION_INBREED_OUTPUT_H
