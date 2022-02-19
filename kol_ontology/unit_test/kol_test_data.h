//
// Created by kellerberrin on 9/4/21.
//

#ifndef KGL_KOL_TEST_DATA_H
#define KGL_KOL_TEST_DATA_H


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the necessary file and directory definitions.
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class UnitTestDefinitions {

public:

  // Just static definitions.
  UnitTestDefinitions() = delete;

  // GoGraphImpl file names.
  [[nodiscard]] static std::string oboFileName() { return std::string(GRAPH_DIRECTORY_) + std::string(GO_OBO_); }
  [[nodiscard]] static std::string newOboFileName() { return std::string(GRAPH_DIRECTORY_) + std::string("new_") + std::string(GO_OBO_); }
  [[nodiscard]] static std::string xmlFileName() { return std::string(GRAPH_DIRECTORY_) + std::string(GO_XML_); }
  // Association file names.
  [[nodiscard]] static std::string gafFileName() { return std::string(ANNOTATION_DIRECTORY_) + std::string(HUMAN_GAF_ANNOTATION_); }
  [[nodiscard]] static std::string newGafFileName() { return std::string(ANNOTATION_DIRECTORY_) + std::string("new_") + std::string(HUMAN_GAF_ANNOTATION_); }
  [[nodiscard]] static std::string geneFileName() { return std::string(ANNOTATION_DIRECTORY_) + std::string(GENE_ASSOCIATION_); }
  [[nodiscard]] static std::string entrezFileName() { return std::string(ANNOTATION_DIRECTORY_) + std::string(HUMAN_ENTREZ_ANNOTATION_); }
  // Matrix File Names for Different Ontologies.
  [[nodiscard]] static std::string matrixFileNameBP() { return std::string(MATRIX_DIRECTORY_) + std::string(MATRIX_FILE_BP_); }
  [[nodiscard]] static std::string matrixFileNameMF() { return std::string(MATRIX_DIRECTORY_) + std::string(MATRIX_FILE_MF_); }
  [[nodiscard]] static std::string matrixFileNameCC() { return std::string(MATRIX_DIRECTORY_) + std::string(MATRIX_FILE_CC_); }

private:

  // Modify these directory constants to change the relative (or absolute) location of the data files.
  static const constexpr char* GRAPH_DIRECTORY_ = "Additional/ggtk/example_graphs/";
  static const constexpr char* ANNOTATION_DIRECTORY_ = "Additional/ggtk/example_annotations/";
  static const constexpr char* MATRIX_DIRECTORY_ = "Additional/ggtk/python/tests/matrix_files/";


  static const constexpr char* GO_OBO_ = "go-basic.obo";
  static const constexpr char* GO_XML_ = "go_daily-termdb.obo-xml";

  static const constexpr char* HUMAN_GAF_ANNOTATION_ = "goa_human.gaf";
  static const constexpr char* GENE_ASSOCIATION_ = "gene_association.mgi";
  static const constexpr char* HUMAN_ENTREZ_ANNOTATION_ = "human_gene2go";

  static const constexpr char* MATRIX_FILE_BP_ = "test_bp_mat.txt";
  static const constexpr char* MATRIX_FILE_MF_ = "test_mf_mat.txt";
  static const constexpr char* MATRIX_FILE_CC_ = "test_cc_mat.txt";

};


#endif //KGL_KOL_TEST_DATA_H
