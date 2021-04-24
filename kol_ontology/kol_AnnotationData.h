/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ANNOTATION_DATA
#define KGL_ANNOTATION_DATA

#include "kol_SetUtilities.h"
#include "kol_GoEnums.h"
#include "kol_GoGraph.h"

#include <fstream>
#include <vector>
#include <iostream> 
#include <map>
#include <string>


namespace kellerberrin::ontology {


//! A class for storing information about genes annotated with go terms.
/*!
	This class hold all information about a set of annotations for genes annotated with go terms.
	 It holds a list of genes, a list of go terms, as well as mappings from a gene to a list of go terms,
	 and mappings from a go term to a list of annotated genes. This class allows querying go annotations
	 and their evidence codes.
*/
class AnnotationData {

public:


  //! Class constructor.
  /*!
    This constructor initialized each vector as an empty vector of the correct type.
  */
  AnnotationData() = default;
  //! class destructor
  /*!
    This destructor clears all maps and vectors.
  */
  ~AnnotationData() = default;


  //! A Method to add annotations to the dataset.
  /*!
    This method adds annotations to the database. It takes a gene, a goTerm, and an evidence code.
      This method checks existence and indexing to remove the burden from parser implementations.
  */
  void addAssociation(const std::string &gene, const std::string &goTerm, const std::string &evidenceCode);


  //! This method tests the existence of a term in the database.
  /*!
    A helper method to check if a term exists in the database.
  */
  [[nodiscard]] bool hasGoTerm(const std::string &goTerm) const { return _stringToGo.find(goTerm) != _stringToGo.end(); }


  //! This method tests the existence of a gene in the database.
  /*!
    A helper method to check if a gene exists in the database.
  */
  [[nodiscard]] bool hasGene(const std::string &gene) const { return _stringToGene.find(gene) != _stringToGene.end(); }

  //! This method returns all the go terms in the database
  /*!
    A helper method to return all the GO terms in the database
  */
  [[nodiscard]] const std::vector<std::string> &getAllGoTerms() const { return _goTerms; }

  //! This method returns all genes in the database
  /*!
    A helper method to return all the genes in the databse
  */
  [[nodiscard]] const std::vector<std::string> &getAllGenes() const { return _genes; }

  //! This method gets the go terms for a gene.
  /*!
    A helper method to return, for a gene, a list of go terms as a vector of strings.
  */
  [[nodiscard]] std::vector<std::string> getGoTermsForGene(const std::string &gene) const;

  //! This method gets the go terms for a gene within the specified onotlogy.
  /*!
    A helper method to return a list of go terms as a set of strings for a gene
    given the sub ontology BIOLOGICAL_PROCESS, MOLECULAR_FUNCTION, or CELLULAR_COMPONENT.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneByOntology(const std::string &gene, GO::Ontology filterOntology, const GoGraph &G) const;
  //! This method gets the biological process go terms for a gene.
  /*!
    A helper method to return a list of BIOLOGICAL_PROCESS go terms for a gene.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneBP(const std::string &gene, const GoGraph &graph) const {

    return getGoTermsForGeneByOntology(gene, GO::Ontology::BIOLOGICAL_PROCESS, graph);

  }

  //! This method gets the molecular function go terms for a gene.
  /*!
    A helper method to return a list of MOLECULAR_FUNCTION go terms for a gene.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneMF(const std::string &gene, const GoGraph &graph) const {

    return getGoTermsForGeneByOntology(gene, GO::Ontology::MOLECULAR_FUNCTION, graph);

  }

  //! This method gets the cellular component go terms for a gene.
  /*!
  A helper method to return a list of CELLULAR_COMPONENT go terms for a gene.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneCC(const std::string &gene, const GoGraph &graph) const {

    return getGoTermsForGeneByOntology(gene, GO::Ontology::CELLULAR_COMPONENT, graph);

  }

  //! A method to get the evidence codes for a list of go terms.
  /*!
    This method returns the evidence codes for a list of go terms. It parallels the
      getGoTermsForGene method and is used for printing and testing.
  */
  [[nodiscard]] std::vector<std::string> getGoTermsEvidenceForGene(const std::string &gene) const;

  //! This method gets the genes for a go term.
  /*!
    A helper method to return, for a go term, a list of genes as a vector of strings.
  */
  [[nodiscard]] std::vector<std::string> getGenesForGoTerm(const std::string &goTerm) const;

  //! This method adds the genes for a go term to a set.
  /*!
    A helper method to add genes associated to a term to a set of genes.
      Used in enrichment calculation
  */
  void addGenesForGoTerm(const std::string &goTerm, OntologySetType<std::string> &geneSet) const;
  //! A method to get the evidence codes for a list of genes.
  /*!
    This method returns the evidence codes for a list of genes. It parallels the
      getGenesForGoTerm method and is used for printing and testing.
  */
  [[nodiscard]] std::vector<std::string> getGenesEvidenceForGoTerm(const std::string &goTerm) const;

  //! A method to get the number of annotations of a particular go term.
  /*!
    This method returns the number of annotations for a go term. Queries the data base
      rather than extracting a vector. Used to calculate information content.
  */
  [[nodiscard]] size_t getNumAnnotationsForGoTerm(const std::string &goTerm) const;

  //! A method to get the number of annotations of a particular gene.
  /*!
    This method returns the number of annotations for a go term. Queries the data base
      rather than extracting a vector.
  */
  [[nodiscard]] size_t getNumAnnotationsForGene(const std::string &gene) const;
  //! A helper method to get the number of genes in the db
  /*!
    This method reutrns the size of the _genes vector.
  */
  [[nodiscard]] size_t getNumGenes() const { return _genes.size(); }


  //! A helper method to get the number of go terms in the db
  /*!
    This method reutrns the size of the _goTerms vector.
  */
  [[nodiscard]] size_t getNumGoTerms() const { return _goTerms.size(); }

  //!	A helper method to return only the terms of the give ontology.
  /*!
    Returns only those terms used that occur for the given ontology.
  */
  [[nodiscard]] std::vector<std::string> getOntologyTerms(const GoGraph &graph, GO::Ontology ontology) const;

private:

  /////////////////////////////////////////////////////////
  //  Lists of gene names and go terms used as mapping keys
  /////////////////////////////////////////////////////////
  //! A list of genes stored by the annotation data object.
  /*!
    This storage variable stores the gene names.
  */
  std::vector<std::string> _genes;

  //! A list of go terms stored by the annotation data object.
  /*!
    This storage variable stores the go terms.
  */
  std::vector<std::string> _goTerms;



  /////////////////////////////////
  // A map of keys to index in list
  /////////////////////////////////
  //! A map from a gene strings to a gene index.
  /*!
    This map accespts gene strings and returns gene indices.
      boost unordered_map ensures O(1) constant time find/has_key queries (hash table).
  */
  OntologyMapType<std::string, std::size_t> _stringToGene;

  //! A map from a go term strings to a go term index.
  /*!
    This map accespts go term strings and returns go term indices.
      boost unordered_map ensures O(1) constant time find/has_key queries (hash table).
  */
  OntologyMapType<std::string, std::size_t> _stringToGo;
  ////////////////////////////////////////////////////////////////////////////////
  // Main data storage objects 2d vectors storing gos for genes and genes for gos.
  ////////////////////////////////////////////////////////////////////////////////
  //! A list of lists of genes, one for each go term.
  /*!
    This vector holds one entry for each go term. Each entry holds a list of genes
      annotated to that go term.
  */
  std::vector<std::vector<std::size_t> > _goToGenes;
  //! A list of lists of evidence codes, one for each go term. Parallel to _goToGenes.
  /*!
    This vector holds one entry for each go term. Each entry holds a list of evidence
      codes for each gene annotated to that go term. It parallels the _goToGenes vectors
      having the same size and dimensions for each element.
  */
  std::vector<std::vector<GO::EvidenceCode> > _goToGenesEvidence;

  //! A list of lists of go terms, one for each gene.
  /*!
    This vector holds one entry for each gene. Each entry holds a list of go terms
      annotated to that gene.
  */
  std::vector<std::vector<std::size_t> > _geneToGos;
  //! A list of lists of evidence codes, one for each gene. Parallel to _geneToGos.
  /*!
    This vector holds one entry for each gene. Each entry holds a list of evidence
      codes for each go term annotated to that gene. It parallels the _geneToGos vectors
      having the same size and dimensions for each element.
  */
  std::vector<std::vector<GO::EvidenceCode> > _geneToGosEvidence;


};

}


#endif

