//
// Created by kellerberrin on 27/5/21.
//

#ifndef KGL_KOL_NEWANNOTATIONDATA_H
#define KGL_KOL_NEWANNOTATIONDATA_H



#include "kol_SetUtilities.h"
#include "kol_GoEnums.h"
#include "kol_GoGraph.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <map>
#include <string>


namespace kellerberrin::ontology {

// A map keyed by gene name, value is a map of GO terms and associated ontologies and evidence codes.
using GeneAnnotationMap = OntologyMapType<std::string, OntologyMapType<std::string, std::pair<GO::Ontology, GO::EvidenceCode>>>;
// A map keyed by go term, value is the GO ontology and a map of annotated genes with associated evidence codes.
using GOAnnotationMap = OntologyMapType<std::string, std::pair<GO::Ontology, OntologyMapType<std::string, GO::EvidenceCode>>>;

//! A class for storing information about genes annotated with go terms.
/*!
	This class hold all information about a set of annotations for genes annotated with go terms.
	 It holds a list of genes, a list of go terms, as well as mappings from a gene to a list of go terms,
	 and mappings from a go term to a list of annotated genes. This class allows querying go annotations
	 and their evidence codes.
*/
class AnnotationDataNew {

public:


  //! Class constructor.
  /*!
    This constructor initialized each vector as an empty vector of the correct type.
  */
  AnnotationDataNew() = default;
  //! class destructor
  /*!
    This destructor clears all maps and vectors.
  */
  ~AnnotationDataNew() = default;


  //! A Method to add annotations to the dataset.
  /*!
    This method adds annotations to the database. It takes a gene, a goTerm, and an evidence code.
      This method checks existence and indexing to remove the burden from parser implementations.
  */
  bool addAssociation(const std::string &gene_id, const std::string &go_term, GO::Ontology go_ontology, GO::EvidenceCode evidence_code);


  //! This method tests the existence of a term in the database.
  /*!
    A helper method to check if a term exists in the database.
  */
  [[nodiscard]] bool hasGoTerm(const std::string &goTerm) const { return go_annotation_map_.contains(goTerm); }


  //! This method tests the existence of a gene in the database.
  /*!
    A helper method to check if a gene exists in the database.
  */
  [[nodiscard]] bool hasGene(const std::string &gene) const { return gene_annotation_map_.contains(gene); }

  //! This method returns all the go terms in the database
  /*!
    A helper method to return all the GO terms in the database
  */
  [[nodiscard]] const GOAnnotationMap& getAllGoTerms() const { return go_annotation_map_; }

  //! This method returns all genes in the database
  /*!
    A helper method to return all the genes in the databse
  */
  [[nodiscard]] const GeneAnnotationMap& getAllGenes() const { return gene_annotation_map_; }

  //! This method gets the go terms for a gene.
  /*!
    A helper method to return, for a gene, a list of go terms as a vector of strings.
  */
  [[nodiscard]] const std::vector<std::string>& getGoTermsForGene(const std::string &gene) const;

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
  [[nodiscard]] size_t getNumGenes() const { return gene_annotation_map_.size(); }


  //! A helper method to get the number of go terms in the db
  /*!
    This method reutrns the size of the _goTerms vector.
  */
  [[nodiscard]] size_t getNumGoTerms() const { return go_annotation_map_.size(); }

  //!	A helper method to return only the terms of the give ontology.
  /*!
    Returns only those terms used that occur for the given ontology.
  */
  [[nodiscard]] std::vector<std::string> getOntologyTerms(GO::Ontology ontology) const;

private:

  GeneAnnotationMap gene_annotation_map_;
  GOAnnotationMap go_annotation_map_;

};

}


#endif //KGL_KOL_NEWANNOTATIONDATA_H
