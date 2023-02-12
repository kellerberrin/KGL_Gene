//
// Created by kellerberrin on 1/4/21.
//

#define BOOST_TEST_NO_MAIN
//#define BOOST_TEST_MODULE "Kellerberrin Ontology Library Unit Test"
//#define BOOST_TEST_DYN_LINK
#include "kel_exec_env_app.h"

#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;
namespace kel = kellerberrin;


test_suite*
init_unit_test_suite( int argc, char ** argv )
{

  framework::master_test_suite().p_name.value = "Kellerberrin Ontology Library Unit Test";

  return nullptr;
}

// The Runtime environment.
class TestExecEnv {

public:

  TestExecEnv()=delete;
  ~TestExecEnv()=delete;

  // The following 4 static members are required for all applications.
  inline static constexpr const char* VERSION = "0.9";
  inline static constexpr const char* MODULE_NAME = "kolTest";
  inline static constexpr const char* LOG_FILE = "kol_test.log";
  inline static constexpr const size_t MAX_ERROR_MESSAGES = 1000;
  inline static constexpr const size_t MAX_WARNING_MESSAGES = 1000;

  static void executeApp() {

    boost::unit_test::unit_test_main(&init_unit_test_suite, argc_, argv_);

  }
  [[nodiscard]] static bool parseCommandLine(int argc, char const ** argv) {

    argc_ = argc;
    argv_ = const_cast<char **>(argv);
    return true;

  }

  // Create application logger.
  [[nodiscard]] static std::unique_ptr<kel::Logger> createLogger() {

    return kel::ExecEnv::createLogger(MODULE_NAME, LOG_FILE, MAX_ERROR_MESSAGES, MAX_WARNING_MESSAGES);

  }


private:

  inline static int argc_;
  inline static char ** argv_;

};


int main(int argc, const char* argv[]) {

  return kel::ExecEnv::runApplication<TestExecEnv>(argc, argv);

}

