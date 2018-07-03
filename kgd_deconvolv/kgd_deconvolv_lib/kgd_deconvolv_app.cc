//
// Created by kellerberrin on 4/05/18.
//

#include <iostream>
#include <kgd_mcmc_ibd.h>
#include <kgd_mcmc_hap.h>
#include "kgd_mt_random.h"
#include "kgd_deploid_io.h"

#include "kgd_deconvolv_app.h"
#include "kgl_utility.h"


namespace kgd = kellerberrin::deconvolv;


// Static private member declarations.
kgd::DeconvolvArgs kgd::ExecEnv::args_;

// Public static member functions.
const kgd::DeconvolvArgs& kgd::ExecEnv::getArgs() { return args_; }

// Constants for the executable.
constexpr const char* kgd::ExecEnv::MODULE_NAME;
constexpr const char* kgd::ExecEnv::VERSION;


void kgd::ExecEnv::executeApp() {

  try {

    std::shared_ptr<DEploidIO> dEploidIO_ptr(std::make_shared<DEploidIO>());

    std::shared_ptr<MersenneTwister> random_generator(std::make_shared<MersenneTwister>(dEploidIO_ptr->getRandomSeed()));

    if ( dEploidIO_ptr->getMixtureControl().doComputeLLK() ) {

      dEploidIO_ptr->computeLLKfromInitialHap();

    } else if ( dEploidIO_ptr->getMixtureControl().doLsPainting() ) {

      dEploidIO_ptr->chromPainting(random_generator);

    } else if ( dEploidIO_ptr->getMixtureControl().doIbdPainting() ) {

      dEploidIO_ptr->paintIBD(random_generator);

    } else {

      if (dEploidIO_ptr->useIBD()) { // ibd

        std::shared_ptr<McmcSample> ibdMcmcSample(std::make_shared<McmcSample>());

        MCMCIBD ibdMcmcMachinery(dEploidIO_ptr, ibdMcmcSample, random_generator);

        ibdMcmcMachinery.runMcmcChain(true /* show progress */);

        dEploidIO_ptr->writeMcmcRelated(ibdMcmcSample, true /* ibd */);

      }

      std::shared_ptr<McmcSample> mcmcSample(std::make_shared<McmcSample>());

      MCMCHAP hapMcmc(dEploidIO_ptr, mcmcSample, random_generator);

      hapMcmc.runMcmcChain(true /* show progress */);

      dEploidIO_ptr->writeMcmcRelated(mcmcSample, false /* mcmchap */);

      dEploidIO_ptr->paintIBD(random_generator);

    }
    // Finishing, write log
    dEploidIO_ptr->wrapUp();

  }
  catch (const std::exception &e) {

    ExecEnv::log().critical("Caught Runtime Error: {}", e.what());

  }

}
