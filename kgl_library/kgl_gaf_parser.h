//
// Created by kellerberrin on 26/01/18.
//

#ifndef KGL_GAF_PARSER_H
#define KGL_GAF_PARSER_H

#include <string>
#include <vector>

#include "kgl_genome_types.h"


namespace kellerberrin::genome {   //  organization level namespace


class GAFRecord {

public:

  GAFRecord() = default;
  GAFRecord(const GAFRecord&) = default;
  ~GAFRecord() = default;

  std::string DB_ID;              // required
  FeatureIdent_t gene_id;    // required
  std::string symbolic_ref;	  // required
  std::string Qualifier;    	    // optional
  OntologyIdent_t ontolotgy_id;	        // required
  std::string DB_Reference;       // required
  std::string	Evidence_Code;	    // required
  std::string With_From;	        // optional
  std::string Aspect;             // required
  std::string description;	      // optional
  std::string	alt_symbolic_ref;  // required
  std::string DB_Object_Type;	    // required
  std::string	Taxon;              // required
  std::string Date;               // required
  std::string	Assigned_By;        // required
  std::string Annotation;         // optional
  std::string	Gene_Product;       // optional

};


using GafGoMap = std::multimap<OntologyIdent_t, GAFRecord>;
class OntologyRecord {

public:

  OntologyRecord(FeatureIdent_t gene_id,
                std::string symbolic_ref,
                std::string alt_symbolic_ref,
                std::string description) : gene_id_(std::move(gene_id)),
                                           symbolic_ref_(std::move(symbolic_ref)),
                                           alt_symbolic_ref_(std::move(alt_symbolic_ref)),
                                           description_(std::move(description)) {}

  OntologyRecord(const OntologyRecord&) = default;
  ~OntologyRecord() = default;

  void addGafRecord(GAFRecord gaf_record);

  [[nodiscard]] const FeatureIdent_t& gene_id() const { return gene_id_; }
  [[nodiscard]] const std::string& symbolicReference() const { return symbolic_ref_;	}
  [[nodiscard]] const std::string& description() const { return description_; }
  [[nodiscard]] const std::string& altSymbolicReference() const { return alt_symbolic_ref_; }
  [[nodiscard]] const GafGoMap& goRecords() const { return go_records_; }

private:

  FeatureIdent_t gene_id_;        // required (1) key
  std::string symbolic_ref_;	    // required (1) (gene family).
  std::string alt_symbolic_ref_;  // required (1) (alt gene family)
  std::string description_;	      // optional (1)
  GafGoMap go_records_;    // vector of GO: ontology records for the gene.

};




/// \class GeneOntology
///	\brief A class to parse and store a gene ontology gaf file.
///
///	This class will read a gaf file and add those annotations to the GafRecordMap container.
/// Indexed by gene id (const FeatureIdent_t&).
///
using GafRecordMap = std::map<FeatureIdent_t, std::shared_ptr<OntologyRecord>>;
class GeneOntology {

public:

  GeneOntology() = default;
  ~GeneOntology() = default;
  GeneOntology(const GeneOntology&) = default;

  /// Reads and parses a gaf 2.0 or 2.1 file.
  ///
  ///  \param const std::string& filename.
  ///  \return bool indicating if the read was successful.
  ///
  bool readGafFile(const std::string& filename);

  /// Returns a map with all gene indexed gaf annotations.
  ///
  /// \return const GafRecordMap&
  ///
  const GafRecordMap& getMap() const { return gaf_record_map_; }


  /// Returns the gaf annotations of a gene.
  ///
  ///  \param const FeatureIdent_t& gene_feature. This parameter specifies the gene.
  ///  \param std::shared_ptr<const OntologyRecord>. The gaf records are returned with this obj pointer.
  ///  Important: the std::shared_ptr == null_ptr if there are no ontology records for this gene.
  ///  \return bool false if no annotations for gene
  bool getGafFeatureVector(const FeatureIdent_t& gene_feature, std::shared_ptr<const OntologyRecord>& ontology_record) const;

private:

  constexpr static const size_t REPORT_INCREMENT_ = 50000;
  constexpr static const size_t GO_FIELD_COUNT_ = 17;

  void parseGafRecord(const std::string& record_str);
  void semanticGafParse(const std::vector<std::string> &field_vec);

  GafRecordMap gaf_record_map_;

};


}   // end namespace




#endif //KGL_KGL_GAF_PARSER_H
