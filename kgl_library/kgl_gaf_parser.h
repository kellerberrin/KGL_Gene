//
// Created by kellerberrin on 26/01/18.
//

#ifndef KGL_GAF_PARSER_H
#define KGL_GAF_PARSER_H

#include <string>
#include <vector>

#include "kgl_genome_types.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class GAFRecord {

public:

  std::string DB_ID;              // required
  FeatureIdent_t DB_Object_ID;    // required
  std::string DB_Object_Symbol;	  // required
  std::string Qualifier;    	    // optional
  std::string GO_ID;	            // required
  std::string DB_Reference;       // required
  std::string	Evidence_Code;	    // required
  std::string With_From;	        // optional
  std::string Aspect;             // required
  std::string Object_Name;	      // optional
  std::string	DB_Object_Synonym;  // required
  std::string DB_Object_Type;	    // required
  std::string	Taxon;              // required
  std::string Date;               // required
  std::string	Assigned_By;        // required
  std::string Annotation;         // optional
  std::string	Gene_Product;       // optional

};

using GafRecordVector = std::vector<const GAFRecord>;
using GafRecordMap = std::map<FeatureIdent_t, std::shared_ptr<GafRecordVector>>;

class GeneOntology {

public:

  GeneOntology() = default;
  ~GeneOntology() = default;
  GeneOntology(const GeneOntology&) = default;

  bool readGafFile(const std::string& filename);

  const GafRecordMap& getMap() const { return gaf_record_map_; }

  // Note that the requested Gaf record may not exist
  // always check the return flag or pointer = null_ptr
  bool getGafFeatureVector(const FeatureIdent_t& gene_feature, std::shared_ptr<const GafRecordVector>& gaf_vector_ptr);

private:


  constexpr static size_t REPORT_INCREMENT_ = 50000;

  void parseGafRecord(const std::string& record_str);

  GafRecordMap gaf_record_map_;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_KGL_GAF_PARSER_H
