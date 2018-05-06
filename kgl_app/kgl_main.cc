//
// Created by kellerberrin on 30/09/17.
//
#include "kgl_phylogenetic_app.h"


namespace kgl = kellerberrin::genome;


/// The mainline.
int main(int argc, char const ** argv)
{

  return kgl::ExecEnv::runApplication<kgl::PhylogeneticExecEnv>(argc, argv);

}

