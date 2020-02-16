//
// Created by kellerberrin on 30/09/17.
//
#include "kgl_phylogenetic_app.h"

/// The mainline.
int main(int argc, char const ** argv)
{

  namespace kgl = kellerberrin::genome;
  namespace kel = kellerberrin;

  return kel::ExecEnv::runApplication<kgl::PhylogeneticExecEnv>(argc, argv);

}

