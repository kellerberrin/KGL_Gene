//
// Created by kellerberrin on 12/12/19.
//

#include "kpl_geneticcode.h"


namespace kpl = kellerberrin::phylogenetic;


// The Codon Amino Acid conversion tables.

kpl::GeneticCode::genetic_code_definitions_t kpl::GeneticCode::_definitions = {
    // codon order is alphabetical: i.e. AAA, AAC, AAG, AAT, ACA, ..., TTT
    {"standard",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"vertmito",             "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"yeastmito",            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"moldmito",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"invertmito",           "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"ciliate",              "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"},
    {"echinomito",           "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"euplotid",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"},
    {"plantplastid",         "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"altyeast",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"ascidianmito",         "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"altflatwormmito",      "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"},
    {"blepharismamacro",     "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"},
    {"chlorophyceanmito",    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"},
    {"trematodemito",        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"scenedesmusmito",      "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"},
    {"thraustochytriummito", "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"}
};


// member function bodies go here

kpl::GeneticCode::GeneticCode() {
  //std::cout << "Constructing a standard GeneticCode" << std::endl;
  useGeneticCode("standard");
}

kpl::GeneticCode::GeneticCode(std::string name) {
  //std::cout << "Constructing a " << name << " GeneticCode" << std::endl;
  useGeneticCode(name);
}

kpl::GeneticCode::~GeneticCode() {
  //std::cout << "Destroying a GeneticCode" << std::endl;
}


std::string kpl::GeneticCode::getGeneticCodeName() const {
  return _genetic_code_name;
}


void kpl::GeneticCode::useGeneticCode(std::string name) {
  _genetic_code_name = name;
  buildGeneticCodeTranslators();
}


unsigned kpl::GeneticCode::getNumNonStopCodons() const {
  return (unsigned) _codons.size();
}


int kpl::GeneticCode::getStateCode(int triplet_index) const {
  return _genetic_code_map.at(triplet_index);
}


char kpl::GeneticCode::getAminoAcidAbbrev(unsigned aa_index) const {
  return _all_amino_acids[aa_index];
}


void kpl::GeneticCode::copyCodons(std::vector<std::string> &codon_vect) const {
  codon_vect.resize(_codons.size());
  std::copy(_codons.begin(), _codons.end(), codon_vect.begin());
}


void kpl::GeneticCode::copyAminoAcids(std::vector<unsigned> &aa_vect) const {
  aa_vect.resize(_amino_acids.size());
  std::copy(_amino_acids.begin(), _amino_acids.end(), aa_vect.begin());
}



void kpl::GeneticCode::buildGeneticCodeTranslators() {
  _amino_acid_map.clear();
  for (unsigned i = 0; i < 20; ++i) {
    char aa = _all_amino_acids[i];
    _amino_acid_map[aa] = i;
  }

  ensureGeneticCodeNameIsValid(_genetic_code_name);
  std::string gcode_aa = _definitions[_genetic_code_name];  // e.g. "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"

  int k = 0;
  int state_code = 0;
  _codons.clear();
  _amino_acids.clear();
  _genetic_code_map.clear();
  for (char ch : gcode_aa) {
    if (ch != '*') {
      _genetic_code_map[k] = state_code++;
      _codons.push_back(_all_codons[k]);
      _amino_acids.push_back(_amino_acid_map[ch]);
    }
    ++k;
  }
}


kpl::GeneticCode::genetic_code_names_t kpl::GeneticCode::getRecognizedGeneticCodeNames() {

  genetic_code_names_t names;

  for (auto it = _definitions.begin(); it != _definitions.end(); ++it) {

    names.push_back(it->first);

  }

  return names;

}


bool kpl::GeneticCode::isRecognizedGeneticCodeName(const std::string &name) {
  std::string lcname = name;
  boost::to_lower(lcname);
  return (_definitions.find(lcname) != _definitions.end());
}


void kpl::GeneticCode::ensureGeneticCodeNameIsValid(const std::string &name) {
  if (!isRecognizedGeneticCodeName(name)) {
    auto valid_genetic_code_names = getRecognizedGeneticCodeNames();
    std::cout << "Recognized genetic codes:\n";
    for (std::string name : valid_genetic_code_names) {
      std::cout << "  " << name << "\n";
    }
    std::cout << std::endl;
    throw XStrom(boost::format("%s is not a recognized genetic code") % name);
  }
}
