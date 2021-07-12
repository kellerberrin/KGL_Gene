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
#include "kol_TermAnnotation.h"

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
  bool readGafFile( const std::string& filename);

  /// Returns a map with all gene indexed gaf annotations.
  ///
  /// \return const GafRecordMap&
  ///
  [[nodiscard]] const std::vector<std::shared_ptr<const kol::GAFRecord>>& getGafRecordVector() const { return gaf_record_vector_; }


private:


  std::vector<std::shared_ptr<const kol::GAFRecord>> gaf_record_vector_;

};



}   // end namespace




#endif //KGL_KGL_GAF_PARSER_H
