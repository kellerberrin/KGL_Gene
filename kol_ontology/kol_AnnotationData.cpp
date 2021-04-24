//
// Created by kellerberrin on 16/4/21.
//

#include "kol_AnnotationData.h"

namespace kol = kellerberrin::ontology;


//! A Method to add annotations to the dataset.
/*!
  This method adds annotations to the database. It takes a gene, a goTerm, and an evidence code.
    This method checks existence and indexing to remove the burden from parser implementations.
*/
void kol::AnnotationData::addAssociation(const std::string &gene, const std::string &goTerm, const std::string &evidenceCode) {

  //Index variables
  std::size_t geneIndex;
  std::size_t goIndex;

  //try{

  //first time finding gene
  if (_stringToGene.find(gene) == _stringToGene.end()) {
    //set the index for gene (it will be inserted next, reason for _genes.size())
    _stringToGene[gene] = _genes.size();
    //add gene to the list
    try {
      _genes.push_back(gene);
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      std::cout << "push_back _gene" << std::endl;
      std::cin.get();
    }

    //initialize the vectors for go and evidence
    try {
      _geneToGos.emplace_back(std::vector<std::size_t>());
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      std::cout << "push_back _geneToGos" << std::endl;
      std::cout << "push_back _geneToGos" << _geneToGos.size() << std::endl;
      std::cin.get();
    }

    try {
      _geneToGosEvidence.emplace_back(std::vector<GO::EvidenceCode>());
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      std::cout << "push_back _geneToGosEvidence" << std::endl;
      std::cin.get();
    }


  }//end if, first finding gene
  //get the index in constant time
  geneIndex = _stringToGene[gene];


  /*
  }catch(std::exception e){
    std::cout << e.what() << std::endl;
    std::cout << "add gene" << std::endl;
    std::cout << gene << " " << goTerm  << " " << evidenceCode << std::endl;
    std::cout << "size " << _genes.size() << std::endl;
  }
  */

  try {

    //first time finding term
    if (_stringToGo.find(goTerm) == _stringToGo.end()) {
      //get the index for go term, next to be inserted, use vec.size()
      _stringToGo[goTerm] = _goTerms.size();
      //add goTerm to list
      _goTerms.push_back(goTerm);


      //initialize the vectors for genes and evidence
      _goToGenes.emplace_back(std::vector<std::size_t>());
      _goToGenesEvidence.emplace_back(std::vector<GO::EvidenceCode>());

    }//end if, first finding goTerm
    //get the index in constant time
    goIndex = _stringToGo[goTerm];

  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    std::cout << "add term" << std::endl;
    std::cout << gene << " " << goTerm << " " << evidenceCode << std::endl;
    std::cout << "size " << _goTerms.size() << std::endl;
  }


  try {

    //get the EvidenceCode enum from the string
    GO::EvidenceCode eCode = GO::evidenceStringToCode(evidenceCode);

    try {
      //add go term (index) to gene's list
      _geneToGos.at(geneIndex).push_back(goIndex);
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      std::cout << "evidence _geneToGos" << std::endl;
      std::cin.get();
    }

    try {
      //add evidence code to gene's list, in parallel with go term
      _geneToGosEvidence.at(geneIndex).push_back(eCode);
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      std::cout << "evidence _geneToGosEvidence" << std::endl;
      std::cin.get();
    }


    try {
      //add gene (index) to go's list
      _goToGenes.at(goIndex).push_back(geneIndex);
    } catch (std::exception &e) {
      std::cout << e.what() << std::endl;
      std::cout << "evidence _goToGenes" << std::endl;
      std::cout << "size      " << _goToGenes.size() << std::endl;
      std::cout << "size sub  " << _goToGenes.at(goIndex).size() << std::endl;
      std::cin.get();
    }

    try {
      //add evidence code to go's list, in parallel with gene
      _goToGenesEvidence.at(goIndex).push_back(eCode);

    } catch (std::exception &e) {

      std::cout << e.what() << std::endl;
      std::cout << "evidence _goToGenesEvidence" << std::endl;
      std::cout << "size      " << _goToGenesEvidence.size() << std::endl;
      std::cout << "size sub  " << _goToGenesEvidence.at(goIndex).size() << std::endl;
      std::cin.get();

    }

  } catch (std::exception &e) {

    std::cout << e.what() << std::endl;
    std::cout << "add evidence" << std::endl;
    std::cout << gene << " " << goTerm << " " << evidenceCode << std::endl;
    std::cout << "size " << _goTerms.size() << std::endl;
    std::cin.get();

  }

}//end method addAssciation

//! This method gets the go terms for a gene.
/*!
  A helper method to return, for a gene, a list of go terms as a vector of strings.
*/
std::vector<std::string> kol::AnnotationData::getGoTermsForGene(const std::string &gene) const {

  std::vector<std::string> goTerms;

  //test if gene exists,prevent key error
  if (hasGene(gene)) {

    auto const&[key, index] = *(_stringToGene.find(gene));

    //move other the indices placing the correct go in the list
    for (auto const &element : _geneToGos.at(index)) {

      goTerms.push_back(_goTerms.at(element));

    }

  }

  return goTerms;

}//end method, getGoTermsForGene


//! This method gets the go terms for a gene within the specified onotlogy.
/*!
  A helper method to return a list of go terms as a set of strings for a gene
  given the sub ontology BIOLOGICAL_PROCESS, MOLECULAR_FUNCTION, or CELLULAR_COMPONENT.
*/
kol::OntologySetType<std::string> kol::AnnotationData::getGoTermsForGeneByOntology(const std::string &gene, GO::Ontology filterOntology, const GoGraph &G) const {
  //temparary storage variable
  OntologySetType<std::string> goTerms;

  //test if gene exists,prevent key error
  if (hasGene(gene)) {

    auto const&[key, index] = *(_stringToGene.find(gene));

    //move other the indices placeing the correct go in the list
    for (auto const &element : _geneToGos.at(index)) {

      std::string term = _goTerms.at(element);
      if (G.getTermOntology(term) == filterOntology) {

        goTerms.insert(term);

      }

    }

  }

  return goTerms;

}


//! A method to get the evidence codes for a list of go terms.
/*!
  This method returns the evidence codes for a list of go terms. It parallels the
    getGoTermsForGene method and is used for printing and testing.
*/
std::vector<std::string> kol::AnnotationData::getGoTermsEvidenceForGene(const std::string &gene) const {

  std::vector<std::string> eCodes;
  if (hasGene(gene)) {

    auto const&[key, index] = *(_stringToGene.find(gene));
    for (auto const &go_element : _geneToGosEvidence.at(index)) {

      std::string code = GO::evidenceToString(go_element);
      eCodes.push_back(code);

    }

  }

  return eCodes;

}

//! This method gets the genes for a go term.
/*!
  A helper method to return, for a go term, a list of genes as a vector of strings.
*/
std::vector<std::string> kol::AnnotationData::getGenesForGoTerm(const std::string &goTerm) const {
  //temparary storage variable
  std::vector<std::string> genes;

  //test if gene exists,prevent key error
  if (hasGoTerm(goTerm)) {

    auto const&[key, index] = *(_stringToGo.find(goTerm));

    //move other the indices placeing the correct go in the list
    for (auto const &element : _goToGenes.at(index)) {

      genes.push_back(_genes.at(element));

    }

  }

  return genes;

}//end method, getGenesForGoTerm

//! This method adds the genes for a go term to a set.
/*!
  A helper method to add genes associated to a term to a set of genes.
    Used in enrichment calculation
*/
void kol::AnnotationData::addGenesForGoTerm(const std::string &goTerm, OntologySetType<std::string> &geneSet) const {
  //test if gene exists,prevent key error
  if (hasGoTerm(goTerm)) {

    auto const&[key, index] = *(_stringToGo.find(goTerm));

    //move other the indices placeing the correct go in the list
    for (auto const &element : _goToGenes.at(index)) {

      geneSet.insert(_genes.at(element));

    }

  }

}

//! A method to get the evidence codes for a list of genes.
/*!
  This method returns the evidence codes for a list of genes. It parallels the
    getGenesForGoTerm method and is used for printing and testing.
*/
std::vector<std::string> kol::AnnotationData::getGenesEvidenceForGoTerm(const std::string &goTerm) const {

  std::vector<std::string> eCodes;
  if (hasGoTerm(goTerm)) {

    auto const&[key, index] = *(_stringToGo.find(goTerm));
    for (auto const &go_element : _goToGenesEvidence.at(index)) {

      std::string code = GO::evidenceToString(go_element);
      eCodes.push_back(code);

    }
  }

  return eCodes;

}


//! A method to get the number of annotations of a particular go term.
/*!
  This method returns the number of annotations for a go term. Queries the data base
    rather than extracting a vector. Used to calculate information content.
*/
size_t kol::AnnotationData::getNumAnnotationsForGoTerm(const std::string &goTerm) const {

  if (not hasGoTerm(goTerm)) {

    return 0;

  }

  auto const&[key, index] = *(_stringToGo.find(goTerm));

  return _goToGenes.at(index).size();

}

//! A method to get the number of annotations of a particular gene.
/*!
  This method returns the number of annotations for a go term. Queries the data base
    rather than extracting a vector.
*/
size_t kol::AnnotationData::getNumAnnotationsForGene(const std::string &gene) const {

  if (not hasGene(gene)) {

    return 0;

  }

  auto const&[key, index] = *(_stringToGene.find(gene));

  return _geneToGos.at(index).size();

}

std::vector<std::string> kol::AnnotationData::getOntologyTerms(const GoGraph &graph, GO::Ontology ontology) const {

  std::vector<std::string> ontologyTerms;
  //Use only terms in the annotation database, this will save on space and computation time.
  for (auto const &term : _goTerms) {

    if (graph.getTermOntology(term) == ontology) {

      ontologyTerms.push_back(term);

    }

  }

  return ontologyTerms;

}

