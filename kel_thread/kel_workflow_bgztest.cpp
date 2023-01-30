//
// Created by kellerberrin on 21/01/23.
//

#include "kel_exec_env.h"
#include "kel_exec_env_app.h"
#include "kel_bzip.h"
#include "kel_bzip_workflow.h"
#include "kel_workflow_threads.h"



namespace kellerberrin {   //  organization level namespace


// Holds the Commandline Arguments.
struct CmdLineArgs {

  std::string workDirectory{"./"};
  std::string logFile{"test_bgz.log"};
  int max_error_count{100};
  int max_warn_count{100};
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


template <typename BGZdecoder> void testWorkflowReader(std::string decoder_name, size_t thread_count) {

  size_t count{0};
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  ExecEnv::log().info("{} begins Parsing BGZ file: {}", decoder_name, ExecEnvBGZ::getArgs().bgz_file_name);

  BGZdecoder test_stream(thread_count);
  test_stream.open(ExecEnvBGZ::getArgs().bgz_file_name);
  std::chrono::time_point<std::chrono::system_clock> prior = std::chrono::system_clock::now();
  size_t total_char_size{0};
  while (true) {

    auto line = test_stream.readLine();
    if (not line) break;
    ++count;
    total_char_size += line.value().second->size();
    if (count % 1000000 == 0) {

      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - prior);
      prior = now;
      ExecEnv::log().info("{}, Lines Read: {}, Line queue size: {}, Elapsed Time (ms): {}"
                          , decoder_name, count, test_stream.lineQueue().size(), milliseconds.count());

    }

//    std::this_thread::sleep_for(std::chrono::milliseconds{1});

  }

  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  ExecEnv::log().info("{} ends, Elapsed time (sec): {}, Lines Read: {}, Total char size: {}, Thread count:{}"
                      , decoder_name, elapsed.count(), count, total_char_size, thread_count);

}


template <template<typename> typename Queue, typename Type> void testBoundedTidal( const std::string queue_name
                                                                                 , const size_t producers
                                                                                 , const size_t consumers
                                                                                 , const size_t iterations) {

  const size_t QUEUE_LOW_TIDE_{2000};
  const size_t QUEUE_HIGH_TIDE_{4000};
  const size_t QUEUE_SAMPLE_FREQ_{100};
  std::vector<std::thread> thread_vector;
//  auto queue_ptr = std::make_unique<Queue<Type>>(QUEUE_HIGH_TIDE_, QUEUE_LOW_TIDE_, queue_name, QUEUE_SAMPLE_FREQ_);
  auto queue_ptr = std::make_unique<Queue<Type>>();

  class PushPull {

  public:

    PushPull(size_t iterations, std::unique_ptr<Queue<Type>> queue_ptr) : iterations_(iterations), queue_ptr_(std::move(queue_ptr)) {}

    void push() {

      while (producer_count.fetch_add(1) < iterations_) {

        queue_ptr_->push(Type());

      }

    };

    void pull() {

      static size_t next_count{0};

      while (consumer_count.fetch_add(1) < iterations_) {

        auto volatile dequeued = queue_ptr_->waitAndPop();
        //      if (*dequeued != next_count) {

        //        ExecEnv::log().error("Out of Sequence; expected: {}, dequeued : {}", next_count, *dequeued);

        //      }
        ++next_count;

      }

    }

  private:

    const size_t iterations_;
    std::atomic<size_t> producer_count{0};
    std::atomic<size_t> consumer_count{0};
    std::unique_ptr<Queue<Type>> queue_ptr_;

  };   // ~PushPull
  PushPull push_pull(iterations, std::move(queue_ptr));

  ExecEnv::log().info("Queue: {} Begins; Producers: {}, Consumers: {}, iterations: {}", queue_name, producers, consumers, iterations);
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  // queue the thread pools
  for (size_t i = 0; i < producers; ++i) {

    thread_vector.emplace_back(&PushPull::push, &push_pull);

  }
  for (size_t i = 0; i < consumers; ++i) {

    thread_vector.emplace_back(&PushPull::pull, &push_pull);

  }

// Do the work and finish.
  for (size_t i = 0; i < thread_vector.size(); ++i) {

    thread_vector[i].join();

  }

  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  ExecEnv::log().info("Queue: {}, Elapsed Time (ms): {}", queue_name, elapsed_ms.count());

}


void ExecEnvBGZ::executeApp() {

//  const size_t thread_count{25};
//  testWorkflowReader<BGZStream>("Workflow", thread_count);
//  testWorkflowReader<BGZReader>("Reader", thread_count);

//  return;

  constexpr const size_t producers{20};
  constexpr const size_t consumers{20};
  constexpr const size_t iterations{100000000};

  while (true) {

    testBoundedTidal<MtQueue, size_t>("BoundedMtQueue<size_t>", producers, consumers, iterations);
//      testBoundedTidal<TidalQueue<std::unique_ptr<size_t>>>("TidalQueue<size_t>");

  }

}


} // namespace


/// The mainline.
int main(int argc, char const ** argv)
{
  namespace kel = kellerberrin;

  return kel::ExecEnv::runApplication<kel::ExecEnvBGZ>(argc, argv);

}
