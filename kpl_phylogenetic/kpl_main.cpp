//
// Created by kellerberrin on 10/12/19.
//

#include <iostream>

#include "kel_exec_env.h"
#include "kel_exec_env_app.h"
#include "kpl_strom.h"

/// The mainline.
int main(int argc, char const ** argv)
{
  namespace kpl = kellerberrin::phylogenetic;
  namespace kel = kellerberrin;

  return kel::ExecEnv::runApplication<kpl::Strom>(argc, argv);

}
