//
// Created by kellerberrin on 26/01/18.
//

#ifndef KGL_GAF_PARSER_H
#define KGL_GAF_PARSER_H

#include <string>
#include <vector>
#include <map>

#include "kgl_genome_types.h"
#include "kgl_ensembl_id_parser.h"
#include "kol_ParserGafRecord.h"

namespace kol = kellerberrin::ontology;
namespace kellerberrin::genome {   //  organization level namespace


/// \class GeneOntology
///	\brief A class to parse and store a gene ontology gaf file.
///
///	This class will read a gaf file and add those annotations to the GafRecordMap container.
/// Indexed by gene id (const FeatureIdent_t&).
///

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
  bool readIdFile(const std::string& filename);

  /// Returns a map with all gene indexed gaf annotations.
  ///
  /// \return const GafRecordMap&
  ///
  [[nodiscard]] const GeneSynonymVector& getSynonymVector() const { return synonym_vector_; }
  [[nodiscard]] const std::vector<std::shared_ptr<const kol::GAFRecord>>& getGafRecordVector() const { return gaf_record_vector_; }


private:

  void semanticGafParse(const std::shared_ptr<const kol::GAFRecord>& gaf_record_ptr);

  GeneSynonymVector  synonym_vector_;
  std::vector<std::shared_ptr<const kol::GAFRecord>> gaf_record_vector_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Resort the gene identifiers
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GeneIdentMap = std::map<std::string, std::string>;
class ResortIds {

public:

  ResortIds() = default;
  ~ResortIds() = default;

  void sortByHGNC(const GeneSynonymVector& ident_vector);
  void sortByEnsembl(const GeneSynonymVector& ident_vector);

  [[nodiscard]] const GeneIdentMap& getMap() const { return gene_id_map_; }

private:

  GeneIdentMap gene_id_map_;

};




}   // end namespace




#endif //KGL_KGL_GAF_PARSER_H
