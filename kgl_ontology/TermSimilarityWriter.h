/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/


#ifndef TERM_SIMILARITY_WRITER
#define TERM_SIMILARITY_WRITER

#include <vector>
#include <string>
#include <iostream>

#include <GoEnums.h>
#include <GoGraph.h>
#include <AnnotationData.h>
#include <TermSimilarityInterface.h>

//! A class write a term similarity matrix to file. Companion to PrecomputedMatrixTermSimilarity
/*! \class TermSimilarityWriter
	This class writes a term similarity matrix to file by first calculating all similarity values
	  between pairs of terms. This class will be use to write matrices 
	  read by PrecomputedMatrixTermSimilarity.
*/
class TermSimilarityWriter {


public:
	//! Parameterized constructor
	/*!
		A simple parameterized constructor.
		This class take an instance of the GO Graph.
	*/
	TermSimilarityWriter(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const AnnotationData> annoData)
	: _goGraph(std::move(goGraph)), _annoData(std::move(annoData)) {}
  ~TermSimilarityWriter() = default;

	//! A method to write a term similarity matrix
	/*!
		Calculates and writes the similarity matrix to file.
		Calculates the similarity between all pairs of terms.
		The complexity is O(N^2) * O(Term Pair Calculation Cost).
		O(Term Pair Calculation Cost) is usually near constant time,
		 but some method will be extremely slow and require other methods
	*/
	void writeSimilarityMatrix(const std::shared_ptr<const TermSimilarityInterface>& termSim, const std::string& fileName, long ontology_code) const {

		auto ontology = static_cast<GO::Ontology>(ontology_code);
		std::vector<std::string> ontologyTerms = _annoData->getOntologyTerms(_goGraph, ontology);

		// Initialze a matrix
		std::size_t nTerms = ontologyTerms.size();
		std::vector<std::vector<double> > matrix;
		matrix.reserve(nTerms);

		for (std::size_t i = 0; i < nTerms; ++i){

			matrix.push_back(std::vector<double>());

			for (std::size_t j = 0; j < nTerms; ++j){
				matrix[i].push_back(0.0);
			}

		}

		// --- Main Calculation ---
		for (std::size_t i = 0; i < nTerms; ++i){

			std::string termA = ontologyTerms.at(i);
			matrix[i][i] = termSim->calculateNormalizedTermSimilarity(termA, termA);

			for (std::size_t j = i + 1; j < nTerms; ++j){

				std::string termB = ontologyTerms.at(j);
				double value = termSim->calculateNormalizedTermSimilarity(termA, termB);
				//cout << termA << "\t" << termB << "\t" << value << endl;
				matrix[i][j] = value;
				matrix[j][i] = value;

			}

		}

		// write the matrix
		writeMatrix(fileName, matrix, ontologyTerms);

	}

	//! A method to write a term similarity matrix for Biological Process
	/*!
		Calcualtes and writes the biological process term similarity matrix
	*/
	void writeSimilarityMatrixBP(const std::shared_ptr<const TermSimilarityInterface>& termSim, const std::string& fileName) const {

		writeSimilarityMatrix(termSim, fileName, GO::Ontology::BIOLOGICAL_PROCESS);

	}

	//! A method to write a term similarity matrix for Molecular Function
	/*!
		Calcualtes and writes the molecular function term similarity matrix
	*/
	void writeSimilarityMatrixMF(const std::shared_ptr<const TermSimilarityInterface>& termSim, const std::string& fileName) const {

		writeSimilarityMatrix(termSim, fileName, GO::Ontology::MOLECULAR_FUNCTION);

	}

	//! A method to write a term similarity matrix for Cellular Component
	/*!
	Calcualtes and writes the cellular component term similarity matrix
	*/
	void writeSimilarityMatrixCC(const std::shared_ptr<const TermSimilarityInterface>& termSim, const std::string& fileName) const {

		writeSimilarityMatrix(termSim, fileName, GO::Ontology::CELLULAR_COMPONENT);

	}

private:

	//! A method to write a matrix to file
	/*!
		writes the matrix to file.
	*/
	static void writeMatrix(const std::string& fname, const std::vector<std::vector<double> > &matrix, const std::vector<std::string> &terms){

		std::ofstream out(fname.c_str());

		std::size_t matrix_size = matrix.size();
		for (std::size_t i = 0; i < matrix_size; ++i){

			std::string termA = terms.at(i);
			out << termA;
			for (std::size_t j = 0; j < matrix_size; ++j){

				out << '\t' << matrix[i][j];

			}

			out << '\n';

		}

	}

	//! A pointer to the gene ontology graph
	std::shared_ptr<const GoGraph> _goGraph;
	//! A pointer to the annotation database
	std::shared_ptr<const AnnotationData> _annoData;

};



#endif