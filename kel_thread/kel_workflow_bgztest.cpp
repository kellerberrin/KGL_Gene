//
// Created by kellerberrin on 21/01/23.
//

#include "kel_exec_env.h"
#include "kel_exec_env_app.h"
#include "kel_bzip.h"
#include "kel_bzip_workflow.h"
#include "kel_thread_pool.h"



namespace kellerberrin {   //  organization level namespace


// Holds the Commandline Arguments.
struct CmdLineArgs {

  std::string workDirectory{"./"};
  std::string logFile{"test_bgz.log"};
  int max_error_count{1000};
  int max_warn_count{1000};
  std::string bgz_file_name{"/media/kellerberrin/DataGenome/Genetics/Genome/HomoSapien/GRCh38/Gnomad3_1/gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz"};
};

// The Runtime environment.
class ExecEnvBGZ {

public:

  ExecEnvBGZ() = delete;

  ~ExecEnvBGZ() = delete;

  [[nodiscard]] inline static const CmdLineArgs &getArgs() { return args_; }

  // The following 5 static members are required for all applications.
  inline static constexpr const char *VERSION = "0.1";
  inline static constexpr const char *MODULE_NAME = "Test BGZ";

  static void executeApp(); // Application mainline.
  // Parse command line arguments.
  [[nodiscard]] static bool parseCommandLine(int /* argc */, char const ** /* argv */) { return true; }

  // Create application logger.
  [[nodiscard]] static std::unique_ptr<Logger> createLogger() { return ExecEnv::createLogger(MODULE_NAME, args_.logFile, args_.max_error_count, args_.max_warn_count); }

private:

  inline static CmdLineArgs args_;

};


void testWorkflowReader() {

  size_t count{0};
  std::chrono::time_point<std::chrono::system_clock> prior = std::chrono::system_clock::now();
  ExecEnv::log().info("Begin Workflow Parsing BGZ file: {}", ExecEnvBGZ::getArgs().bgz_file_name);

  BGZStream test_stream(20);
  test_stream.open(ExecEnvBGZ::getArgs().bgz_file_name);
  while (true) {

    auto line = test_stream.readLine();
    if (not line) break;
    ++count;
    if (count % 1000000 == 0) {

      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - prior);
      prior = now;
      ExecEnv::log().info("Workflow, Lines Read: {}, WF Input Size: {}, WF Output Size: {}, Line Queue Size: {}, Elapsed Time (ms): {}", count, test_stream.workFlow().inputQueue().size(), test_stream.workFlow().outputQueue().size(), test_stream.lineQueue().size(), milliseconds.count());

    }

  }

  ExecEnv::log().info("End Workflow, Lines Read: {}", count);

  ExecEnv::log().info("Begin Reader Parsing BGZ file: {}", ExecEnvBGZ::getArgs().bgz_file_name);
  BGZReader test_reader(15);
  test_reader.open(ExecEnvBGZ::getArgs().bgz_file_name);
  count = 0;
  while (true) {

    auto line = test_reader.readLine();
    if (not line) break;
    ++count;
    if (count % 1000000 == 0) {

      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - prior);
      prior = now;
      ExecEnv::log().info("Reader, Lines Read: {}, Decompress Queue Size: {}, Line Queue Size: {}, Elapsed Time (ms): {}", count, test_reader.decompressQueue().size(), test_reader.lineQueue().size(), milliseconds.count());

    }

  }
  ExecEnv::log().info("End Reader, Lines Read: {}", count);

}


template <typename QueueType> void testBoundedTidal(std::string queue_name) {

  constexpr static const size_t QUEUE_LOW_TIDE_{2000};
  constexpr static const size_t QUEUE_HIGH_TIDE_{4000};
  constexpr static const size_t QUEUE_SAMPLE_FREQ_{100};

  std::atomic<size_t> producer_count{0};
  std::atomic<size_t> consumer_count{0};
  const size_t iterations = 10000000;
  WorkflowThreads consumers(10);
  WorkflowThreads producers(40);
  std::vector<std::future<bool>> thread_pool_futures;
  QueueType tidal_queue(QUEUE_HIGH_TIDE_, QUEUE_LOW_TIDE_, queue_name, QUEUE_SAMPLE_FREQ_);
  auto push_lambda = [&]()->bool{

    while (producer_count.fetch_add(1) < iterations) {

      tidal_queue.push(std::make_unique<size_t>(0));

    }
    return true;

  };
  auto pop_lambda = [&]()->bool{

    while (consumer_count.fetch_add(1) < iterations) {

      tidal_queue.waitAndPop();

    }
    return true;

  };

  ExecEnv::log().info("Queue: {} Begins; iterations: {}", queue_name, iterations);
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  // queue the thread pools
  for (size_t i = 0; i < consumers.threadCount(); ++i) {

    thread_pool_futures.push_back(consumers.enqueueTask(pop_lambda));

  }
  for (size_t i = 0; i < producers.threadCount(); ++i) {

    thread_pool_futures.push_back(producers.enqueueTask(push_lambda));

  }

// Do the work and finish.

  for (size_t i = 0; i < thread_pool_futures.size(); ++i) {

    thread_pool_futures[i].get();

  }

  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  ExecEnv::log().info("Queue: {} Size: {}, Elapsed Time (ms): {}", queue_name, tidal_queue.size(), elapsed_ms.count());

}


void ExecEnvBGZ::executeApp() {

//  testWorkflowReader();
  while (true) {

    testBoundedTidal<BoundedMtQueue<std::unique_ptr<size_t>>>("BoundedMtQueue<size_t>");
    testBoundedTidal<TidalQueue<std::unique_ptr<size_t>>>("TidalQueue<size_t>");

  }

}


} // namespace


/// The mainline.
int main(int argc, char const ** argv)
{
  namespace kel = kellerberrin;

  return kel::ExecEnv::runApplication<kel::ExecEnvBGZ>(argc, argv);

}
