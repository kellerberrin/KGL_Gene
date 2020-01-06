//
// Created by kellerberrin on 10/12/19.
//

#include <iostream>

#include <iostream>
#include "kpl_node.h"
#include "kpl_tree.h"
#include "kpl_treemanip.h"
#include "kpl_treesummary.h"
#include "kpl_strom.h"


namespace kpl = kellerberrin::phylogenetic;


int main(int argc, const char * argv[]) {
  std::cout << "Main Starting..." << std::endl;

  kpl::Strom strom;

  try {

    strom.processCommandLineOptions(argc, argv);
    strom.run();

  }
  catch(std::exception & x) {

    std::cerr << "Exception: " << x.what() << std::endl;
    std::cerr << "Aborted." << std::endl;

  }
  catch(...) {

    std::cerr << "Exception of unknown type!\n";

  }


  return 0;
}


