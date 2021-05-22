/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/


#ifndef KGL_TERM_SIMILARITY_WRITER
#define KGL_TERM_SIMILARITY_WRITER

#include <vector>
#include <string>
#include <iostream>

#include "kol_GoEnums.h"
#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_SimilarityInterface.h"

namespace kellerberrin::ontology {


//! A class write a term similarity matrix to file. Companion to SimilarityMatrix
/*! \class SimilarityWriter
	This class writes a term similarity matrix to file by first calculating all similarity values
	  between pairs of terms. This class will be use to write matrices 
	  read by SimilarityMatrix.
*/
class SimilarityWriter {


public:
  //! Parameterized constructor
  /*!
    A simple parameterized constructor.
    This class take an instance of the GO Graph.
  */
  SimilarityWriter(const std::shared_ptr<const GoGraph> &graph_ptr,
                   const std::shared_ptr<const AnnotationData> &annotation_ptr)
      : graph_ptr_(graph_ptr), annotation_ptr_(annotation_ptr) {}

  ~SimilarityWriter() = default;

  //! A method to write a term similarity matrix
  /*!
    Calculates and writes the similarity matrix to file.
    Calculates the similarity between all pairs of terms.
    The complexity is O(N^2) * O(Term Pair Calculation Cost).
    O(Term Pair Calculation Cost) is usually near constant time,
     but some method will be extremely slow and require other methods
  */
  void writeSimilarityMatrix(const std::shared_ptr<const SimilarityInterface> &termSim,
                             const std::string &fileName,
                             GO::Ontology ontology) const;

  //! A method to write a term similarity matrix for Biological Process
  /*!
    Calcualtes and writes the biological process term similarity matrix
  */
  void writeSimilarityMatrixBP(const std::shared_ptr<const SimilarityInterface> &termSim, const std::string &fileName) const {

    writeSimilarityMatrix(termSim, fileName, GO::Ontology::BIOLOGICAL_PROCESS);

  }

  //! A method to write a term similarity matrix for Molecular Function
  /*!
    Calcualtes and writes the molecular function term similarity matrix
  */
  void writeSimilarityMatrixMF(const std::shared_ptr<const SimilarityInterface> &termSim, const std::string &fileName) const {

    writeSimilarityMatrix(termSim, fileName, GO::Ontology::MOLECULAR_FUNCTION);

  }

  //! A method to write a term similarity matrix for Cellular Component
  /*!
  Calcualtes and writes the cellular component term similarity matrix
  */
  void writeSimilarityMatrixCC(const std::shared_ptr<const SimilarityInterface> &termSim, const std::string &fileName) const {

    writeSimilarityMatrix(termSim, fileName, GO::Ontology::CELLULAR_COMPONENT);

  }

private:

  //! A pointer to the gene ontology graph
  std::shared_ptr<const GoGraph> graph_ptr_;
  //! A pointer to the annotation database
  std::shared_ptr<const AnnotationData> annotation_ptr_;

  //! A method to write a matrix to file
  /*!
    writes the matrix to file.
  */
  static void writeMatrix( const std::string &fname,
                           const std::vector<std::vector<double> > &matrix,
                           const std::vector<std::string> &terms);


};

} // namespace

#endif