/*
 * kgl_deploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgl_deploid.
 *
 * kgl_deploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream> // std::cout
#include "mcmc.hpp"
#include "dEploidIO.hpp"
#include "kgl_deploid_app.h"

namespace kgl = kellerberrin::genome;



/// The mainline.
int main(int argc, char const ** argv)
{

    return kgl::ExecEnv::runApplication<kgl::DeploidExecEnv>(argc, argv);

}

/// Original mainline.
int kgl::deploidMain( int argc, const char **argv) {

    try {

        DEploidIO dEploidIO(argc, argv);
        std::ostream *output = &std::cout;

        if ( dEploidIO.version() ) {

            dEploidIO.printVersion(*output);
            return EXIT_SUCCESS;

        }

        if ( dEploidIO.help()) {

            dEploidIO.printHelp(*output);
            return EXIT_SUCCESS;

        }

        if ( dEploidIO.doComputeLLK() ) {

            dEploidIO.computeLLKfromInitialHap();

        } else if ( dEploidIO.doLsPainting() ) {

            dEploidIO.chromPainting();

        } else if ( dEploidIO.doIbdPainting() ) {

            dEploidIO.paintIBD();

        }else{

            if (dEploidIO.useIBD()) { // ibd

                std::shared_ptr<McmcSample> ibdMcmcSample(std::make_shared<McmcSample>());

                MersenneTwister ibdRg(dEploidIO.randomSeed());

                McmcMachinery ibdMcmcMachinery(&dEploidIO, ibdMcmcSample, &ibdRg, true);

                ibdMcmcMachinery.runMcmcChain(true, // show progress
                                              true);  // use IBD

            }

            std::shared_ptr<McmcSample> mcmcSample(std::make_shared<McmcSample>());

            MersenneTwister rg(dEploidIO.randomSeed());

            McmcMachinery mcmcMachinery(&dEploidIO,
                                        mcmcSample,
                                        &rg,
                                        false); // use IBD

            mcmcMachinery.runMcmcChain(true, // show progress
                                       false); // use IBD

            dEploidIO.paintIBD();

        }
        // Finishing, write log
        dEploidIO.wrapUp();

    }
    catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}
