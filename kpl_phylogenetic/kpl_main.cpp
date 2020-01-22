//
// Created by kellerberrin on 10/12/19.
//

#include <iostream>

#include "kgl_exec_env.h"
#include "kpl_strom.h"


namespace kpl = kellerberrin::phylogenetic;

#define EXECENV 1
#ifndef EXECENV


int main(int argc, const char * argv[]) {
  std::cout << "Main Starting..." << std::endl;

  kpl::Strom strom;

  try {

    strom.parseCommandLine(argc, argv);
    strom.executeApp();

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


#else

/// The mainline.
int main(int argc, char const ** argv)
{

  return kpl::ExecEnv::runApplication<kpl::Strom>(argc, argv);

}


#endif