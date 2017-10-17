//
// Created by kellerberrin on 10/10/17.
//


#include <iostream>
#include "kgl_exec_env.h"
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/timer/timer.hpp>


// Define namespace alias
namespace bt = boost::timer;
namespace kgl = kellerberrin::genome;



// Hide the boost cpu timer in an anonymous namespace.
namespace {  bt::cpu_timer cpu_timer; }

void kgl::ExecEnv::getElpasedTime(double& Clock, double& User, double& System) {

  Clock = 0; User = 0; System = 0;
  bt::cpu_times elapsedtime = cpu_timer.elapsed();
  Clock = elapsedtime.wall / 1e09; // Convert from nanoseconds to seconds
  User = elapsedtime.user / 1e09;
  System = elapsedtime.system / 1e09;

}

