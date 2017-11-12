//
// Created by kellerberrin on 10/10/17.
//


#include <iostream>
#include "kgl_exec_env.h"
#include <boost/timer/timer.hpp>
#define BOOST_FILESYSTEM_NO_DEPRECATED // Recommended by boost filesystem documentation.
#include <boost/filesystem.hpp>


// Define namespace alias
namespace bt = boost::timer;
namespace fs = boost::filesystem;
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


// Returns the filename with the path directory appended to it "path/file".
std::string kgl::ExecEnv::filePath(const std::string& file_name, const std::string& path) {

  fs::path directory_path = fs::path(path);
  fs::path file_path = directory_path / fs::path(file_name);
  return file_path.string();

}

