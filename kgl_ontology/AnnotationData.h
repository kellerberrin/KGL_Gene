/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef ANNOTATION_DATA
#define ANNOTATION_DATA

#include <SetUtilities.h>
#include <GoEnums.h>
#include <GoGraph.h>

#include <fstream>
#include <vector>
#include <iostream> 
#include <map>
#include <string>


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
	void addAssociation(const std::string &gene, const std::string &goTerm, const std::string &evidenceCode) {

		//Index variables
		std::size_t geneIndex;
		std::size_t goIndex;

		//try{

		//first time finding gene
		if(_stringToGene.find(gene) == _stringToGene.end()){
			//set the index for gene (it will be inserted next, reason for _genes.size())
			_stringToGene[gene] = _genes.size();
			//add gene to the list
			try{
				_genes.push_back(gene);
			}catch(std::exception& e){
				std::cout << e.what() << std::endl;
				std::cout << "push_back _gene" << std::endl;
				std::cin.get();
			}

			//initialize the vectors for go and evidence
			try{
				_geneToGos.emplace_back(std::vector<std::size_t>());
			}catch(std::exception& e){
				std::cout << e.what() << std::endl;
				std::cout << "push_back _geneToGos" << std::endl;
				std::cout << "push_back _geneToGos" << _geneToGos.size() << std::endl;
				std::cin.get();
			}

			try{
				_geneToGosEvidence.emplace_back(std::vector<GO::EvidenceCode>());
			}catch(std::exception& e){
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
		
		try{

		//first time finding term
		if(_stringToGo.find(goTerm) == _stringToGo.end()){
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

		}catch(std::exception& e){
			std::cout << e.what() << std::endl;
			std::cout << "add term" << std::endl;
			std::cout << gene << " " << goTerm  << " " << evidenceCode << std::endl;
			std::cout << "size " << _goTerms.size() << std::endl;
		}


		try{

		//get the EvidenceCode enum from the string
		GO::EvidenceCode eCode = GO::evidenceStringToCode(evidenceCode);

		try{
		//add go term (index) to gene's list
		_geneToGos.at(geneIndex).push_back(goIndex);
		}catch(std::exception& e){
			std::cout << e.what() << std::endl;
			std::cout << "evidence _geneToGos" << std::endl;
			std::cin.get();
		}

		try{
		//add evidence code to gene's list, in parallel with go term
		_geneToGosEvidence.at(geneIndex).push_back(eCode);
		}catch(std::exception& e){
			std::cout << e.what() << std::endl;
			std::cout << "evidence _geneToGosEvidence" << std::endl;
			std::cin.get();
		}



		try{
		//add gene (index) to go's list
		_goToGenes.at(goIndex).push_back(geneIndex);
		}catch(std::exception& e){
			std::cout << e.what() << std::endl;
			std::cout << "evidence _goToGenes" << std::endl;
			std::cout << "size      " << _goToGenes.size() << std::endl;
			std::cout << "size sub  " << _goToGenes.at(goIndex).size() << std::endl;
			std::cin.get();
		}

		try{
		//add evidence code to go's list, in parallel with gene
		_goToGenesEvidence.at(goIndex).push_back(eCode);

		} catch(std::exception& e) {

			std::cout << e.what() << std::endl;
			std::cout << "evidence _goToGenesEvidence" << std::endl;
			std::cout << "size      " <<  _goToGenesEvidence.size() << std::endl;
			std::cout << "size sub  " <<  _goToGenesEvidence.at(goIndex).size() << std::endl;
			std::cin.get();

		}

  } catch(std::exception& e) {

    std::cout << e.what() << std::endl;
    std::cout << "add evidence" << std::endl;
    std::cout << gene << " " << goTerm  << " " << evidenceCode << std::endl;
    std::cout << "size " << _goTerms.size() << std::endl;
    std::cin.get();

  }
	
	}//end method addAssciation


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
  [[nodiscard]] const std::vector<std::string>& getAllGoTerms() const { return _goTerms; }

	//! This method returns all genes in the database
	/*!
		A helper method to return all the genes in the databse
	*/
	[[nodiscard]] const std::vector<std::string>& getAllGenes() const { return _genes; }

	//! This method gets the go terms for a gene.
	/*!
		A helper method to return, for a gene, a list of go terms as a vector of strings.
	*/
	[[nodiscard]] std::vector<std::string> getGoTermsForGene(const std::string &gene) const {

		std::vector<std::string> goTerms;

		//test if gene exists,prevent key error
		if (hasGene(gene)) {

			auto const& [key, index] = *(_stringToGene.find(gene));

			//move other the indices placing the correct go in the list
			for(auto const& element : _geneToGos.at(index)) {

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
	[[nodiscard]] OntologySetType<std::string> getGoTermsForGeneByOntology(const std::string &gene, GO::Ontology filterOntology, const GoGraph& G) const {
		//temparary storage variable
    OntologySetType<std::string> goTerms;

		//test if gene exists,prevent key error
		if (hasGene(gene)) {

			auto const& [key, index] = *(_stringToGene.find(gene));

			//move other the indices placeing the correct go in the list
			for(auto const& element : _geneToGos.at(index)) {

				std::string term = _goTerms.at(element);
				if (G.getTermOntology(term) == filterOntology) {

          goTerms.insert(term);

				}

			}

		}

		return goTerms;
		
	}

	//! This method gets the biological process go terms for a gene.
	/*!
		A helper method to return a list of BIOLOGICAL_PROCESS go terms for a gene.
	*/
  [[nodiscard]] OntologySetType<std::string> getGoTermsForGeneBP(const std::string &gene, const GoGraph& graph) const {

		return getGoTermsForGeneByOntology(gene, GO::Ontology::BIOLOGICAL_PROCESS, graph);

	}

	//! This method gets the molecular function go terms for a gene.
	/*!
		A helper method to return a list of MOLECULAR_FUNCTION go terms for a gene.
	*/
	[[nodiscard]] OntologySetType<std::string> getGoTermsForGeneMF(const std::string &gene, const GoGraph& graph) const {

		return getGoTermsForGeneByOntology(gene, GO::Ontology::MOLECULAR_FUNCTION, graph);

	}

	//! This method gets the cellular component go terms for a gene.
	/*!
	A helper method to return a list of CELLULAR_COMPONENT go terms for a gene.
	*/
	[[nodiscard]] OntologySetType<std::string> getGoTermsForGeneCC(const std::string &gene, const GoGraph& graph) const {

		return getGoTermsForGeneByOntology(gene, GO::Ontology::CELLULAR_COMPONENT, graph);

	}

	//! A method to get the evidence codes for a list of go terms.
	/*!
		This method returns the evidence codes for a list of go terms. It parallels the 
		  getGoTermsForGene method and is used for printing and testing.
	*/
	[[nodiscard]] std::vector<std::string> getGoTermsEvidenceForGene(const std::string &gene) const {

		std::vector<std::string> eCodes;
		if( hasGene(gene) ) {

			auto const& [key, index] = *(_stringToGene.find(gene));
			for(auto const& go_element : _geneToGosEvidence.at(index)) {

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
	[[nodiscard]] std::vector<std::string> getGenesForGoTerm(const std::string &goTerm) const {
		//temparary storage variable
		std::vector<std::string> genes;

		//test if gene exists,prevent key error
		if (hasGoTerm(goTerm)) {

			auto const& [key, index] = *(_stringToGo.find(goTerm));

			//move other the indices placeing the correct go in the list
			for(auto const& element : _goToGenes.at(index)) {

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
	void addGenesForGoTerm(const std::string &goTerm, OntologySetType<std::string> &geneSet) const {
		//test if gene exists,prevent key error
		if(hasGoTerm(goTerm)) {

			auto const& [key, index] = *(_stringToGo.find(goTerm));

			//move other the indices placeing the correct go in the list
			for(auto const& element : _goToGenes.at(index)) {

				geneSet.insert(_genes.at(element));

			}

		}

	}

	//! A method to get the evidence codes for a list of genes.
	/*!
		This method returns the evidence codes for a list of genes. It parallels the 
		  getGenesForGoTerm method and is used for printing and testing.
	*/
	[[nodiscard]] std::vector<std::string> getGenesEvidenceForGoTerm(const std::string &goTerm) const {

		std::vector<std::string> eCodes;
		if (hasGoTerm(goTerm)) {

			auto const& [key, index] = *(_stringToGo.find(goTerm));
			for(auto const& go_element : _goToGenesEvidence.at(index)) {

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
	[[nodiscard]] size_t getNumAnnotationsForGoTerm(const std::string &goTerm) const {

		if (not hasGoTerm(goTerm)) {

			return 0;

		}

		auto const& [key, index] = *(_stringToGo.find(goTerm));

		return _goToGenes.at(index).size();

	}

	//! A method to get the number of annotations of a particular gene.
	/*!
		This method returns the number of annotations for a go term. Queries the data base
		  rather than extracting a vector.
	*/
	[[nodiscard]] size_t getNumAnnotationsForGene(const std::string &gene) const {

		if (not hasGene(gene)){

			return 0;

		}

		auto const& [key, index] = *(_stringToGene.find(gene));

		return _geneToGos.at(index).size();

	}

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
	[[nodiscard]] std::vector<std::string> getOntologyTerms(const GoGraph& graph, GO::Ontology ontology) const {

		std::vector<std::string> ontologyTerms;
		//Use only terms in the annotation database, this will save on space and computation time.
		for (auto const& term : _goTerms) {

			if (graph.getTermOntology(term) == ontology){

				ontologyTerms.push_back(term);

			}

		}

		return ontologyTerms;

	}

	// Accessor routines to check database integrity.
	// Number of genes sorted by gene index.
  [[nodiscard]] const std::vector<std::string>& genes() const { return _genes; }
  // Number of go terms sorted by go index.
  [[nodiscard]] const std::vector<std::string>& goTerms() const { return _goTerms; }
  // Number of genes, value = gene index.
  [[nodiscard]] const OntologyMapType<std::string,std::size_t>& geneIndex() const { return _stringToGene; }
  // Number of go terms value = go index.
  [[nodiscard]] const OntologyMapType<std::string,std::size_t>& goIndex() const { return _stringToGo; }
  // Number of go terms X the gene indices for each go term.
  [[nodiscard]] const std::vector<std::vector<std::size_t> >& goIndexToGeneIndex() const { return _goToGenes; }
  // Number of go terms X the gene evidence code for each each. Same dimensions as above.
  [[nodiscard]] const std::vector<std::vector<GO::EvidenceCode> >& goIndexToGeneEvidence() const { return _goToGenesEvidence; }
  // Number of genes X the index of the go terms for each gene.
  [[nodiscard]] const std::vector<std::vector<std::size_t> >& geneIndexToGoIndex() const { return _geneToGos; }
  // Number of genes X the go evidence code for each gene. Same dimensions as above.
  [[nodiscard]] const std::vector<std::vector<GO::EvidenceCode> >& geneIndexToGoEvidence() const { return _geneToGosEvidence; }

  [[nodiscard]] bool databaseIntegrityCheck() const {

	  // Ensure the uniqueness of the gene and go vectors.
	  auto const gene_set = SetUtilities::convert_vector(genes());
	  if (gene_set.size() != genes().size()) {

	    return false;

	  }
    auto const go_set = SetUtilities::convert_vector(genes());
    if (go_set.size() != goTerms().size()) {

      return false;

    }
    // Check the size of the index lookups
    if (genes().size() != geneIndex().size()) {

      return false;

    }

    // Check the size of the index lookups
    if (goTerms().size() != goIndex().size()) {

      return false;

    }

    if (goTerms().size() != goIndexToGeneIndex().size()) {

      return false;

    }

    if (goTerms().size() != goIndexToGeneEvidence().size()) {

      return false;

    }


    for (size_t idx = 0; idx < goTerms().size(); ++idx) {

      if (goIndexToGeneIndex()[idx].size() != goIndexToGeneEvidence()[idx].size()) {

        return false;

      }

    }

    if (genes().size() != geneIndexToGoIndex().size()) {

      return false;

    }

    if (genes().size() != geneIndexToGoEvidence().size()) {

      return false;

    }

    for (size_t idx = 0; idx < genes().size(); ++idx) {

      if (geneIndexToGoIndex()[idx].size() != geneIndexToGoEvidence()[idx].size()) {

        return false;

      }

    }

    return true;

  }



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
  OntologyMapType<std::string,std::size_t> _stringToGene;

  //! A map from a go term strings to a go term index.
  /*!
    This map accespts go term strings and returns go term indices.
      boost unordered_map ensures O(1) constant time find/has_key queries (hash table).
  */
  OntologyMapType<std::string,std::size_t> _stringToGo;
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
#endif

