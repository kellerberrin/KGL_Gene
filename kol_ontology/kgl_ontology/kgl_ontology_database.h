//
// Created by kellerberrin on 12/4/21.
//

#ifndef KGL_ONTOLOGY_DATABASE_H
#define KGL_ONTOLOGY_DATABASE_H


#include <string>
#include <memory>

#include "kol_TermAnnotation.h"

namespace kellerberrin::ontology {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Forward declaration of graph and annotation objects hides boost implementation detail from consumers of
// ontology functionality.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GoGraph;
class TermAnnotation;


class OntologyDatabase {

public:

  OntologyDatabase(const std::string& ontology_ident,
                   const std::string& go_graph_file,
                   const std::string& annotation_file);

  ~OntologyDatabase() = default;

  [[nodiscard]] const std::string& ontologyIdent() const { return ontology_ident_; }
  [[nodiscard]] const std::shared_ptr<const GoGraph>& goGraph() const { return go_graph_ptr_; }
  [[nodiscard]] const std::shared_ptr<const TermAnnotation>& annotation() const { return annotation_ptr_; }

private:

  const std::string ontology_ident_;
  std::shared_ptr<const GoGraph> go_graph_ptr_;
  std::shared_ptr<const TermAnnotation> annotation_ptr_;

  [[nodiscard]] std::shared_ptr<const TermAnnotation> getAnnotation(const std::string& annotation_file);
  [[nodiscard]] std::shared_ptr<const GoGraph> getGoGraph(const std::string& go_graph_file);

};

} // namespace

#endif //KGL_ONTOLOGY_DATABASE_H
