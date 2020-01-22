//
// Created by kellerberrin on 4/05/18.
//

#include <iostream>
#include <kgd_mcmc_ibd.h>
#include <kgd_mcmc_hap.h>
#include "kgd_deploid_io.h"

#include "kgd_deconvolv_app.h"
#include "kgl_utility.h"


namespace kgd = kellerberrin::deconvolv;


// Static private member declarations.
kgd::DeconvolvArgs kgd::Deconvolv::args_;

// Public static member functions.
const kgd::DeconvolvArgs& kgd::Deconvolv::getArgs() { return args_; }

// Constants for the executable.
constexpr const char* kgd::Deconvolv::MODULE_NAME;
constexpr const char* kgd::Deconvolv::VERSION;


void kgd::Deconvolv::executeApp() {

  try {

    std::shared_ptr<DEploidIO> dEploidIO_ptr(std::make_shared<DEploidIO>());

    if ( dEploidIO_ptr->getMixtureControl().doComputeLLK() ) {

      dEploidIO_ptr->computeLLKfromInitialHap();

    } else if ( dEploidIO_ptr->getMixtureControl().doLsPainting() ) {

      dEploidIO_ptr->chromPainting();

    } else if ( dEploidIO_ptr->getMixtureControl().doIbdPainting() ) {

      dEploidIO_ptr->paintIBD();

    } else {

      if (dEploidIO_ptr->useIBD()) { // ibd

        std::shared_ptr<McmcSample> ibdMcmcSample(std::make_shared<McmcSample>());

        MCMCIBD ibdMcmcMachinery(dEploidIO_ptr, ibdMcmcSample);

        ibdMcmcMachinery.runMcmcChain(true /* show progress */);

        dEploidIO_ptr->writeMcmcRelated(ibdMcmcSample, true /* ibd */);

      }

      std::shared_ptr<McmcSample> mcmcSample(std::make_shared<McmcSample>());

      MCMCHAP hapMcmc(dEploidIO_ptr, mcmcSample);

      hapMcmc.runMcmcChain(true /* show progress */);

      dEploidIO_ptr->writeMcmcRelated(mcmcSample, false /* mcmchap */);

      dEploidIO_ptr->paintIBD();

    }
    // Finishing, write log
    dEploidIO_ptr->wrapUp();

  }
  catch (const std::exception &e) {

    DeconvolvApp::log().critical("Caught Runtime Error: {}", e.what());

  }

}


void kgd::Deconvolv::executeLib(const MixtureDataObj& mixture_data) {

  try {

    std::shared_ptr<DEploidIO> dEploidIO_ptr(std::make_shared<DEploidIO>(mixture_data));

    dEploidIO_ptr->hapParameters().setMcmcSample(500);
    dEploidIO_ptr->hapParameters().setMcmcMachineryRate(2);
    dEploidIO_ptr->ibdParameters().setMcmcSample(50);
    dEploidIO_ptr->ibdParameters().setMcmcMachineryRate(2);
    dEploidIO_ptr->ibdParameters().setProposalScaling(50.0);

    dEploidIO_ptr->getMixtureControl().setUseIBD(true);

    if ( dEploidIO_ptr->getMixtureControl().doComputeLLK() ) {

      dEploidIO_ptr->computeLLKfromInitialHap();

    } else if ( dEploidIO_ptr->getMixtureControl().doLsPainting() ) {

      dEploidIO_ptr->chromPainting();

    } else if ( dEploidIO_ptr->getMixtureControl().doIbdPainting() ) {

      dEploidIO_ptr->paintIBD();

    } else {

      if (dEploidIO_ptr->useIBD()) { // ibd

        std::shared_ptr<McmcSample> ibdMcmcSample(std::make_shared<McmcSample>());

        MCMCIBD ibdMcmcMachinery(dEploidIO_ptr, ibdMcmcSample);

        ibdMcmcMachinery.runMcmcChain(true /* show progress */);

        dEploidIO_ptr->writeMcmcRelated(ibdMcmcSample, true /* ibd */);

      }

      std::shared_ptr<McmcSample> mcmcSample(std::make_shared<McmcSample>());

      MCMCHAP hapMcmc(dEploidIO_ptr, mcmcSample);

      hapMcmc.runMcmcChain(true /* show progress */);

      dEploidIO_ptr->writeMcmcRelated(mcmcSample, false /* mcmchap */);

      dEploidIO_ptr->paintIBD();

    }
    // Finishing, write log
    dEploidIO_ptr->wrapUp();

  }
  catch (const std::exception &e) {

    DeconvolvApp::log().critical("Caught Runtime Error: {}", e.what());

  }

}

