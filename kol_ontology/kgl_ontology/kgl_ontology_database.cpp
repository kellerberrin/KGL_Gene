//
// Created by kellerberrin on 24/4/21.
//


#include "kgl_ontology_database.h"

#include "kol_ParserGoFactory.h"
#include "kol_GoGraph.h"
#include "kol_ParserAnnotationFactory.h"
#include "kol_AnnotationData.h"



namespace kol = kellerberrin::ontology;



kol::OntologyDatabase::OntologyDatabase( const std::string& ontology_ident,
                                         const std::string& go_graph_file,
                                         const std::string& annotation_file) : ontology_ident_(ontology_ident) {

    go_graph_ptr_ = getGoGraph(go_graph_file);
    annotation_ptr_ = getAnnotation(annotation_file);

}


kol::OntologyDatabase::~OntologyDatabase() {

    go_graph_ptr_ = nullptr;
    annotation_ptr_ = nullptr;

}


std::unique_ptr<const kol::AnnotationData> kol::OntologyDatabase::getAnnotation(const std::string& annotation_file) {

    auto anno_parser_ptr = ParserAnnotationFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           PolicyAllowedEvidence());
    return anno_parser_ptr->parseAnnotationFile(annotation_file);

}

std::unique_ptr<const kol::GoGraph> kol::OntologyDatabase::getGoGraph(const std::string& go_graph_file) {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::OBO_GO_STANDARD);
    return go_parser_ptr->parseGoFile(go_graph_file);

}
