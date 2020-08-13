//
// Created by kellerberrin on 12/12/19.
//

#include "kpl_genetictype.h"


namespace kpl = kellerberrin::phylogenetic;


// member function bodies go here

kpl::DataType::DataType() : _datatype(0), _num_states(0) {
  //std::cout << "Creating a DataType object" << std::endl;
  setNucleotide();
}


kpl::DataType::~DataType() {
  //std::cout << "Destroying a DataType object" << std::endl;
}


void kpl::DataType::setNucleotide() {
  _datatype = 1;
  _num_states = 4;
  _genetic_code = nullptr;
}


void kpl::DataType::setCodon() {
  _datatype = 2;
  _genetic_code = GeneticCode::SharedPtr(new GeneticCode("standard"));
  _num_states = _genetic_code->getNumNonStopCodons();
}


void kpl::DataType::setProtein() {
  _datatype = 3;
  _num_states = 20;
  _genetic_code = nullptr;
}


void kpl::DataType::setStandard() {
  _datatype = 4;
  _num_states = 2;
  _genetic_code = nullptr;
}


bool kpl::DataType::isNucleotide() const {

  return (_datatype == 1);

}


bool kpl::DataType::isCodon() const {

  return (_datatype == 2);

}


bool kpl::DataType::isProtein() const {

  return (_datatype == 3);

}


bool kpl::DataType::isStandard() const {

  return (_datatype == 4);

}


void kpl::DataType::setGeneticCodeFromName(std::string genetic_code_name) {

  assert(isCodon());
  _genetic_code = GeneticCode::SharedPtr(new GeneticCode(genetic_code_name));

}


void kpl::DataType::setGeneticCode(GeneticCode::SharedPtr gcode) {

  assert(isCodon());
  assert(gcode);

  _genetic_code = gcode;

}


void kpl::DataType::setStandardNumStates(unsigned nstates) {

  _datatype = 4;
  _num_states = nstates;
  _genetic_code = nullptr;

}


unsigned kpl::DataType::getDataType() const {

  return _datatype;

}


unsigned kpl::DataType::getNumStates() const {

  return _num_states;

}


const kpl::GeneticCode::SharedPtr kpl::DataType::getGeneticCode() const {

  assert(isCodon());
  return _genetic_code;

}


std::string kpl::DataType::getDataTypeAsString() const {

  std::string s = translateDataTypeToString(_datatype);

  if (isCodon()) {

    s += boost::str(boost::format(",%s") % _genetic_code->getGeneticCodeName());

  }

  return s;

}


std::string kpl::DataType::translateDataTypeToString(unsigned datatype) {

  assert(datatype == 1 || datatype == 2 || datatype == 3 || datatype == 4);

  switch(datatype) {

    case 1: return std::string("nucleotide");

    case 2:  return std::string("codon");

    case 3:  return std::string("protein");

    default:  return std::string("standard");

  }

}
