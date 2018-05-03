//
// Created by kellerberrin on 30/09/17.
//
#include "kgl_application.h"
#include "kgl_phylogenetic_env.h"
#include "kgl_phylogenetic_app.h"


namespace kgl = kellerberrin::genome;


/// The mainline.
int main(int argc, char const ** argv)
{

  return kgl::application<kgl::PhylogeneticExecEnv, kgl::PhylogeneticApp>(argc, argv);

}

