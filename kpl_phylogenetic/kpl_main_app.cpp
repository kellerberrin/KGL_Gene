//
// Created by kellerberrin on 8/12/22.
//

#include "kpl_main_app.h"


namespace kpl = kellerberrin::phylogenetic;
namespace kel = kellerberrin;



void kpl::PhyloExecEnv::executeApp() {

  strom_app_.executeApp();

}


std::unique_ptr<kel::ExecEnvLogger> kpl::PhyloExecEnv::createLogger() {

  return ExecEnv::createLogger(MODULE_NAME, args_.logFile, args_.max_error_count, args_.max_warn_count);

}


bool kpl::PhyloExecEnv::parseCommandLine(int argc, char const ** argv) {

  return strom_app_.parseCommandLine(argc, argv);

}
