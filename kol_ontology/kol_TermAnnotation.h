//
// Created by kellerberrin on 27/5/21.
//

#ifndef KOL_TERMANNOTATION_H
#define KOL_TERMANNOTATION_H


#include "kol_GoEnums.h"
#include "kol_OntologyTypes.h"
#include "kol_ParserGafRecord.h"
#include "kol_PolicyEvidence.h"

#include <vector>
#include <map>
#include <string>
#include <memory>


namespace kellerberrin::ontology {

// A map keyed by gene name, value is a map of GO terms and GAFRecords.
using GeneGOMap = OntologyMapType<std::string, std::vector<std::shared_ptr<const GAFRecord>>>;
using GeneAnnotationMap = OntologyMapType<std::string, GeneGOMap>;
// A map keyed by go term, value is a map of annotated genes and GAFRecords.
using GOGeneMap = OntologyMapType<std::string, std::vector<std::shared_ptr<const GAFRecord>>>;
using GOAnnotationMap = OntologyMapType<std::string, GOGeneMap>;


enum class AnnotationGeneName { UNIPROT_GENE_ID, SYMBOLIC_GENE_ID};
//! A class for storing information about genes annotated with go terms.
/*!
	This class hold all information about a set of annotations for genes annotated with go terms.
	 It holds a list of genes, a list of go terms, as well as mappings from a gene to a list of go terms,
	 and mappings from a go term to a list of annotated genes. This class allows querying go annotations
	 and their evidence codes.
*/
class TermAnnotation {

public:


  explicit TermAnnotation( const std::vector<std::shared_ptr<const GAFRecord>>& gaf_records,
                           AnnotationGeneName gene_id_type = AnnotationGeneName::UNIPROT_GENE_ID) : gene_id_type_(gene_id_type) {

    createAnnotationMap(gaf_records, gene_id_type);

  }

  TermAnnotation( const PolicyEvidence& evidence_policy,
                  const std::vector<std::shared_ptr<const GAFRecord>>& gaf_records,
                  AnnotationGeneName gene_id_type = AnnotationGeneName::UNIPROT_GENE_ID) : gene_id_type_(gene_id_type) {

    createAnnotationMap(filterGAFRecords(evidence_policy, gaf_records), gene_id_type);

  }
  ~TermAnnotation() = default;

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
  [[nodiscard]] std::vector<std::string> getGoTermsForGene(const std::string &gene) const;

  //! This method gets the go terms for a gene within the specified onotlogy.
  /*!
    A helper method to return a list of go terms as a set of strings for a gene
    given the sub ontology BIOLOGICAL_PROCESS, MOLECULAR_FUNCTION, or CELLULAR_COMPONENT.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneByOntology(const std::string &gene, GO::Ontology filterOntology) const;
  //! This method gets the biological process go terms for a gene.
  /*!
    A helper method to return a list of BIOLOGICAL_PROCESS go terms for a gene.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneBP(const std::string &gene) const {

    return getGoTermsForGeneByOntology(gene, GO::Ontology::BIOLOGICAL_PROCESS);

  }

  //! This method gets the molecular function go terms for a gene.
  /*!
    A helper method to return a list of MOLECULAR_FUNCTION go terms for a gene.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneMF(const std::string &gene) const {

    return getGoTermsForGeneByOntology(gene, GO::Ontology::MOLECULAR_FUNCTION);

  }

  //! This method gets the cellular component go terms for a gene.
  /*!
  A helper method to return a list of CELLULAR_COMPONENT go terms for a gene.
  */
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneCC(const std::string &gene) const {

    return getGoTermsForGeneByOntology(gene, GO::Ontology::CELLULAR_COMPONENT);

  }


  [[nodiscard]] const GeneGOMap& getGeneGOMap(const std::string &gene) const;
  [[nodiscard]] const GOGeneMap& getGOGeneMap(const std::string &go_ident) const;
  [[nodiscard]] std::vector<std::shared_ptr<const GAFRecord>> getAllGAFRecords() const;

  //! This method gets the genes for a go term.
  /*!
    A helper method to return, for a go term, a list of genes as a vector of strings.
  */
  [[nodiscard]] std::vector<std::string> getGenesForGoTerm(const std::string &goTerm) const;

  //! A method to get the number of annotations of a particular go term.
  /*!
    This method returns the number of annotations for a go term. Queries the data base
      rather than extracting a vector. Used to calculate information content.
  */
  [[nodiscard]] size_t getNumAnnotationsForGoTerm(const std::string &go_term) const { return getGOGeneMap(go_term).size(); }

  //! A method to get the number of annotations of a particular gene.
  /*!
    This method returns the number of annotations for a go term. Queries the data base
      rather than extracting a vector.
  */
  [[nodiscard]] size_t getNumAnnotationsForGene(const std::string &gene) const { return getGeneGOMap(gene).size(); }
  //! A helper method to get the number of genes in the db
  /*!
    This method reutrns the size of the _genes vector.
  */
  [[nodiscard]] size_t getNumGenes() const { return gene_annotation_map_.size(); }

  void addGenesForGoTerm(const std::string &goTerm, OntologySetType<std::string> &geneSet) const;

  //! A helper method to get the number of go terms in the db
  /*!
    This method reutrns the size of the _goTerms vector.
  */
  [[nodiscard]] size_t getNumGoTerms() const { return go_annotation_map_.size(); }

  //!	A helper method to return only the terms of the given ontology.
  /*!
    Returns only those terms used that occur for the given ontology.
  */
  [[nodiscard]] std::vector<std::string> getOntologyTerms(GO::Ontology ontology) const;

  [[nodiscard]] static std::vector<std::shared_ptr<const GAFRecord>> filterGAFRecords(const PolicyEvidence& evidence_policy,
                                                                                      const std::vector<std::shared_ptr<const GAFRecord>>& gaf_records);

  [[nodiscard]] AnnotationGeneName geneIdType() const { return gene_id_type_; }

private:

  GeneAnnotationMap gene_annotation_map_;
  GOAnnotationMap go_annotation_map_;
  const AnnotationGeneName gene_id_type_;

  bool addGAFRecord(const std::shared_ptr<const GAFRecord>& gaf_record_ptr, AnnotationGeneName gene_name_type);
  void createAnnotationMap(const std::vector<std::shared_ptr<const GAFRecord>>& gaf_records, AnnotationGeneName gene_name_type);
  const std::string& geneName(const std::shared_ptr<const GAFRecord>& gaf_record_ptr, AnnotationGeneName gene_name_type);

};



} // Namespace.


#endif //KOL_ANNOTATIONDATA_H

