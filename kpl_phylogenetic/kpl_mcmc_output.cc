//
// Created by kellerberrin on 16/12/19.
//

#include "kpl_mcmc_output.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::OutputManager::OutputManager() {
  //std::cout << "Constructing an OutputManager" << std::endl;
  _tree_file_name = "trees.t";
  _param_file_name = "params.p";

}


kpl::OutputManager::~OutputManager() {
  //std::cout << "Destroying an OutputManager" << std::endl;
}


void kpl::OutputManager::openTreeFile(std::string filename, Data::SharedPtr data) {

  assert(!_treefile.is_open());
  _tree_file_name = filename;
  _treefile.open(_tree_file_name.c_str());

  if (!_treefile.is_open()) {

    throw XStrom(boost::str(boost::format("Could not open tree file \"%s\"") % _tree_file_name));

  }

  _treefile << "#nexus\n\n";
  _treefile << data->createTaxaBlock() << std::endl;

  _treefile << "begin trees;\n";
  _treefile << data->createTranslateStatement() << std::endl;

}


void kpl::OutputManager::closeTreeFile() {

  assert(_treefile.is_open());
  _treefile << "end;\n";
  _treefile.close();

}


void kpl::OutputManager::openParameterFile(std::string filename, Model::SharedPtr model) {

  assert(model);
  assert(!_parameterfile.is_open());

  _param_file_name = filename;
  _parameterfile.open(_param_file_name.c_str());

  if (!_parameterfile.is_open()) {

    throw XStrom(boost::str(boost::format("Could not open parameter file \"%s\"") % _param_file_name));

  }

  _parameterfile << boost::str(boost::format("%s\t%s\t%s\t%s\t%s\t%s") % "iter" % "lnL" % "lnPr" % "TL" % "m" % model->paramNamesAsString("\t")) << std::endl;

}


void kpl::OutputManager::closeParameterFile() {

  assert(_parameterfile.is_open());
  _parameterfile.close();

}


void kpl::OutputManager::outputConsole(std::string s) {

  std::cout << s << std::endl;

}


void kpl::OutputManager::outputTree(unsigned iter, TreeManip::SharedPtr tm) {

  assert(_treefile.is_open());
  assert(tm);
  _treefile << boost::str(boost::format("  tree iter_%d = %s;") % iter % tm->makeNewick(5)) << std::endl;

}


void kpl::OutputManager::outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, Model::SharedPtr model) {

  assert(model);
  assert(_parameterfile.is_open());

  _parameterfile << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%d\t%s") % iter % lnL % lnP % TL % m % model->paramValuesAsString("\t")) << std::endl;

}
