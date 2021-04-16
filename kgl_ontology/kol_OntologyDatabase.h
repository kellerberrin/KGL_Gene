//
// Created by kellerberrin on 12/4/21.
//

#ifndef KGL_KOL_ONTOLOGYDATABASE_H
#define KGL_KOL_ONTOLOGYDATABASE_H

#include "kol_GoParserFactory.h"
#include "kol_GoGraph.h"
#include "kol_AnnotationParserFactory.h"
#include "kol_AnnotationData.h"


namespace kellerberrin::ontology {


class OntologyDatabase {

public:

  OntologyDatabase(const std::string& ontology_ident,
                   const std::string& go_graph_file,
                   const std::string& annotation_file) : ontology_ident_(ontology_ident) {

    go_graph_ptr_ = getGoGraph(go_graph_file);
    annotation_ptr_ = getAnnotation(annotation_file);

  }
  ~OntologyDatabase() = default;

  [[nodiscard]] const std::string& ontologyIdent() const { return ontology_ident_; }
  [[nodiscard]] const std::shared_ptr<const GoGraph>& goGraph() const { return go_graph_ptr_; }
  [[nodiscard]] const std::shared_ptr<const AnnotationData>& annotation() const { return annotation_ptr_; }

private:

  const std::string ontology_ident_;
  std::shared_ptr<const GoGraph> go_graph_ptr_;
  std::shared_ptr<const AnnotationData> annotation_ptr_;

  [[nodiscard]] std::unique_ptr<const AnnotationData> getAnnotation(const std::string& annotation_file) {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    return anno_parser_ptr->parseAnnotationFile(annotation_file);

  }

  [[nodiscard]] std::unique_ptr<const GoGraph> getGoGraph(const std::string& go_graph_file) {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    return go_parser_ptr->parseGoFile(go_graph_file);

  }

};

} // namespace

#endif //KGL_KOL_ONTOLOGYDATABASE_H
