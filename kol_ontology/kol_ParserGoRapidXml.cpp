//
// Created by kellerberrin on 19/4/21.
//



#include "kol_OntologyTypes.h"
#include "kol_ParserGoRapidXml.h"

namespace kol = kellerberrin::ontology;


//! Method to parse the go file, should be an XML file.
/*!
  This method will read a Gene Ontology XML file and add all relationships
   to the graph.

*/

std::shared_ptr<kol::GoGraph> kol::GoParserRapidXml::parseGoFile(const std::string &filename) const {

  //graph object to be returned
  std::shared_ptr<GoGraph> graph(std::make_shared<GoGraph>());

  //open xmlfile
  rapidxml::file<> xmlFile(filename.c_str());

  //initialize document
  rapidxml::xml_document<> doc;
  doc.parse<0>(xmlFile.data());

  //create an xml node
  rapidxml::xml_node<> *node;

  for (node = doc.first_node("obo")->first_node(); node; node = node->next_sibling()) {

    //get node name
    std::string type = node->name();

    //if node is a GO term
    if (type == "term") {

      std::string term, name, description, ontology;

      //vector to hold go accession ids of related terms.
      std::vector<std::string> relatedTerms;
      std::vector<std::string> relatedTermRelationship;

      bool isObsolete = false;

      rapidxml::xml_node<> *childNode;
      for (childNode = node->first_node(); childNode; childNode = childNode->next_sibling()) {

        std::string attr = childNode->name();
        //attribute specific action
        if (attr == "id") {
          //go accession id
          term = childNode->value();

        } else if (attr == "name") {
          //go natural language name/title
          name = childNode->value();

        } else if (attr == "namespace") {
          //go ontology, BIOLOGICAL_PROCESS,MOLECULAR_FUNCTION,CELLULAR_COMPONENT
          ontology = childNode->value();

        } else if (attr == "def") {
          //go term definition, long form description
          description = childNode->first_node("defstr")->value();

        } else if (attr == "is_a") {
          //is_a relationship
          std::string relatedTerm = childNode->value();

          //add related term to vector for adding later
          relatedTerms.push_back(relatedTerm);

          //add attr as relationship, attr must = "is_a" here
          relatedTermRelationship.push_back(attr);


        } else if (attr == "is_obsolete") {
          //if obsolete set flag.
          isObsolete = true;
          //early exit from loop
          break;

        } else if (attr == "relationship") {
          //extract relationship
          std::string type = childNode->first_node("type")->value();
          //extract term to which related
          std::string toTerm = childNode->first_node("to")->value();

          //add all types of relationships
          relatedTerms.push_back(toTerm);

          relatedTermRelationship.push_back(type);

        }//END if-else blocks, term attributes


      }//for each attribute node

      //if term is obsolete continue, do not add to graph
      if (isObsolete) { continue; }

      ///////////////////////////////////////////////////////
      // Main work, add nodes, add related nodes, add edges
      ///////////////////////////////////////////////////////


      //add the current term to the graph
      graph->insertTerm(term, name, description, ontology);


      //loop over all related terms
      for (std::size_t i = 0; i < relatedTerms.size(); ++i) {

        //create a temp variable for readability
        std::string relatedTerm = relatedTerms.at(i);
        //std::cout << relatedTerm << "\t";

        //create a temp variable for readability
        std::string relationship = relatedTermRelationship.at(i);
        //std::cout << relationship << std::endl;

        //insert related terms, these are just stubs to be overwritten later on
        graph->insertTerm(relatedTerm, "name", "description", "ontology");

        //insert edge
        graph->insertRelationship(term, relatedTerm, relationship);

      }//end for, each related term

    }//end if, is term?

  }//end for, for each node in obo, source, header, term


  //call to initialize the graph's vertex to index maps
  graph->initMaps();

  //return the graph pointer
  return graph;

}//end method parseGoFile

//! A method to test if a file fits the accepted format
/*!
Returns true if the file matches accepted format, false otherwise
*/
bool kol::GoParserRapidXml::isFileGood(const std::string &filename) const {

  std::ifstream in(filename.c_str());
  if (!in.good()) {
    return false;
  } else {
    in.close();
  }


  std::size_t count = 0;

  //open xmlfile
  rapidxml::file<> xmlFile(filename.c_str());

  //initialize document
  rapidxml::xml_document<> doc;
  try {
    doc.parse<0>(xmlFile.data());
  } catch (rapidxml::parse_error &pe) {
    return false;
  }


  //create an xml node
  rapidxml::xml_node<> *node;
  for (node = doc.first_node("obo")->first_node(); node; node = node->next_sibling()) {

    std::string type = node->name();
    if (type == "term") {
      count += 1;
      if (count == 3) {
        break;
      }
    }
  }

  if (count < 3) {
    return false;
  } else {
    return true;
  }
}
