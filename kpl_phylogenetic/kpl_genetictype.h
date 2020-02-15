//
// Created by kellerberrin on 11/12/19.
//

#ifndef KPL_GENETICTYPE_H
#define KPL_GENETICTYPE_H


#include "kpl_geneticcode.h"

#include <boost/format.hpp>


namespace kellerberrin::phylogenetic {   //  organization level namespace

class DataType {
public:
  DataType();

  ~DataType();

  void setNucleotide();

  void setCodon();

  void setProtein();

  void setStandard();

  bool isNucleotide() const;

  bool isCodon() const;

  bool isProtein() const;

  bool isStandard() const;

  void setStandardNumStates(unsigned nstates);

  void setGeneticCodeFromName(std::string genetic_code_name);

  void setGeneticCode(GeneticCode::SharedPtr gcode);

  unsigned getDataType() const;

  unsigned getNumStates() const;

  std::string getDataTypeAsString() const;

  const GeneticCode::SharedPtr getGeneticCode() const;

  static std::string translateDataTypeToString(unsigned datatype);

private:

  unsigned _datatype;
  unsigned _num_states;
  GeneticCode::SharedPtr _genetic_code;
};



} // end Namespace


#endif //KPL_GENETICTYPE_H
