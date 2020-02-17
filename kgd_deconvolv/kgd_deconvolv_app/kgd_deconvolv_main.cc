
#include "kgd_deconvolv_app.h"



/// The mainline.
int main(int argc, char const ** argv)
{

  namespace kgd = kellerberrin::deconvolv;
  namespace kel = kellerberrin;

  return kel::ExecEnv::runApplication<kgd::Deconvolv>(argc, argv);

}

