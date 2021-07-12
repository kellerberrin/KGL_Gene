//
// Created by kellerberrin on 24/4/21.
//


#include "kgl_ontology_database.h"

#include "kol_ParserGoFactory.h"
#include "kol_GoGraph.h"
#include "kol_ParserAnnotationGaf.h"
#include "kol_TermAnnotation.h"



namespace kol = kellerberrin::ontology;
namespace kgl = kellerberrin::genome;



kol::OntologyDatabase::OntologyDatabase( const std::string& ontology_ident,
                                         const std::string& go_graph_file,
                                         const std::string& annotation_file) : kgl::ResourceBase(ontology_ident) {

  go_graph_ptr_ = getGoGraph(go_graph_file);
  annotation_ptr_ = getAnnotation(annotation_file);

}


std::shared_ptr<const kol::TermAnnotation> kol::OntologyDatabase::getAnnotation(const std::string& annotation_file) {

  PolicyEvidence default_evidence;
  return ParserAnnotationGaf::parseAnnotationFile(default_evidence, annotation_file);

}

std::shared_ptr<const kol::GoGraph> kol::OntologyDatabase::getGoGraph(const std::string& go_graph_file) {

  auto relationship_policy = PolicyRelationship(GO::allRelationships());
  auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO, relationship_policy);
  return go_parser_ptr->parseGoFile(go_graph_file);

}
