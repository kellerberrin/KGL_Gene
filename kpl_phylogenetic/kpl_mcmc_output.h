//
// Created by kellerberrin on 16/12/19.
//

#ifndef KPL_MCMC_OUTPUT_H
#define KPL_MCMC_OUTPUT_H


#include "kpl_geneticdata.h"
#include "kpl_treemanip.h"
#include "kpl_model.h"
#include "kpl_xstrom.h"


#include <fstream>



namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace



class OutputManager {
public:
  typedef std::shared_ptr< OutputManager > SharedPtr;

  OutputManager();
  ~OutputManager();

  void                       setModel(Model::SharedPtr model) {_model = model;}
  void                       setTreeManip(TreeManip::SharedPtr tm) {_tree_manip = tm;}

  void                       openTreeFile(std::string filename, Data::SharedPtr data);
  void                       openParameterFile(std::string filename, Model::SharedPtr model);

  void                       closeTreeFile();
  void                       closeParameterFile();

  void                       outputConsole(std::string s);
  void                       outputTree(unsigned iter, TreeManip::SharedPtr tm);
  void                       outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, Model::SharedPtr model);


private:

  TreeManip::SharedPtr       _tree_manip;
  Model::SharedPtr           _model;
  std::ofstream              _treefile;
  std::ofstream              _parameterfile;
  std::string                _tree_file_name;
  std::string                _param_file_name;
};



} // phylogenetic
} // kellerberrin



#endif //KPL_MCMC_OUTPUT_H
