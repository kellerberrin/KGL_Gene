//
// Created by kellerberrin on 11/12/19.
//

#ifndef KPL_GENETICCODE_H
#define KPL_GENETICCODE_H


#include "kpl_xstrom.h"

#include <boost/algorithm/string.hpp>

#include <vector>
#include <map>
#include <iostream>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class Data;
class Model;
class QMatrix;

class GeneticCode {

//  friend class Data;
//  friend class Model;
//  friend class QMatrix;

public:

  typedef std::map<int, int> genetic_code_map_t;
  typedef std::map<char, unsigned> amino_acid_map_t;
  typedef std::vector<unsigned> amino_acid_vect_t;
  typedef std::vector<std::string> codon_vect_t;
  typedef std::vector<char> amino_acid_symbol_vect_t;
  typedef std::map<std::string, std::string> genetic_code_definitions_t;
  typedef std::vector<std::string> genetic_code_names_t;


  GeneticCode();

  GeneticCode(std::string name);

  ~GeneticCode();

  std::string getGeneticCodeName() const;

  void useGeneticCode(std::string name);

  unsigned getNumNonStopCodons() const;

  int getStateCode(int triplet_index) const;

  char getAminoAcidAbbrev(unsigned aa_index) const;

  void copyCodons(std::vector<std::string> &codon_vect) const;

  [[nodiscard]] const codon_vect_t& getCodons() const { return _codons; }

  void copyAminoAcids(std::vector<unsigned> &aa_vect) const;

  static genetic_code_names_t getRecognizedGeneticCodeNames();

  static bool isRecognizedGeneticCodeName(const std::string &name);

  static void ensureGeneticCodeNameIsValid(const std::string &name);

private:

  void buildGeneticCodeTranslators();

  std::string _genetic_code_name;

  genetic_code_map_t _genetic_code_map;
  amino_acid_map_t _amino_acid_map;

  amino_acid_vect_t _amino_acids;
  codon_vect_t _codons;

  const amino_acid_symbol_vect_t _all_amino_acids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
                                                     'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
  const std::vector<std::string> _all_codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC",
                                                "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT",
                                                "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC",
                                                "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                                                "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC",
                                                "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
                                                "TTA", "TTC", "TTG", "TTT"};

  static genetic_code_definitions_t _definitions;

public:

  typedef std::shared_ptr<GeneticCode> SharedPtr;
};



}  // phylogenetic
}  // kellerberrin



#endif //KPL_GENETICCODE_H
