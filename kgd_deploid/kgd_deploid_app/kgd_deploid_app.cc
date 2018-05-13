//
// Created by kellerberrin on 4/05/18.
//

#include <iostream>
#include "kgd_mt_random.h"
#include "kgd_mcmc.h"
#include "kgd_deploid_io.h"

#include "kgd_deploid_app.h"
#include "kgl_utility.h"


namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;


// Static private member declarations.
kgd::DeploidArgs kgd::DeploidExecEnv::args_;

// Public static member functions.
const kgd::DeploidArgs& kgd::DeploidExecEnv::getArgs() { return args_; }

// Constants for the executable.
constexpr const char* kgd::DeploidExecEnv::MODULE_NAME;
constexpr const char* kgd::DeploidExecEnv::VERSION;


void kgd::DeploidExecEnv::executeApp() {

  try {

    std::shared_ptr<DEploidIO> dEploidIO_ptr(std::make_shared<DEploidIO>());

    std::shared_ptr<MersenneTwister> random_generator(std::make_shared<MersenneTwister>(dEploidIO_ptr->randomSeed()));

    if ( dEploidIO_ptr->doComputeLLK() ) {

      dEploidIO_ptr->computeLLKfromInitialHap();

    } else if ( dEploidIO_ptr->doLsPainting() ) {

      dEploidIO_ptr->chromPainting(random_generator);

    } else if ( dEploidIO_ptr->doIbdPainting() ) {

      dEploidIO_ptr->paintIBD(random_generator);

    }else{

      if (dEploidIO_ptr->useIBD()) { // ibd

        std::shared_ptr<McmcSample> ibdMcmcSample(std::make_shared<McmcSample>());

        McmcMachinery ibdMcmcMachinery(dEploidIO_ptr, ibdMcmcSample, random_generator, true /* use IBD */);

        ibdMcmcMachinery.runMcmcChain(true /* show progress */,  true /* use IBD */);

      }

      std::shared_ptr<McmcSample> mcmcSample(std::make_shared<McmcSample>());

      McmcMachinery mcmcMachinery(dEploidIO_ptr, mcmcSample, random_generator, false /* don't use IBD */);

      mcmcMachinery.runMcmcChain(true /* show progress */, false /* don't use IBD */);

      dEploidIO_ptr->paintIBD(random_generator);

    }
    // Finishing, write log
    dEploidIO_ptr->wrapUp();

  }
  catch (const std::exception &e) {

    kgl::ExecEnv::log().critical("Caught Runtime Error: {}", e.what());

  }

}
