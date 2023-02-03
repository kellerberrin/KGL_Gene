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


template <template<typename> typename Queue, typename Type> class PushPull {

public:

  PushPull(size_t iterations) : iterations_(iterations), queue_ptr_(std::make_unique<Queue<Type>>(5000, 2000, "tidal queue", 100)) {}
  ~PushPull() = default;

  void push() {

    while (producer_count_.fetch_add(1) < iterations_) {

      queue_ptr_->push(Type());

      ++push_count_;

    }

  };


  void pullNoWait() {


    while (consumer_count_ < iterations_) {

      auto dequeued = queue_ptr_->pop();
      if (not dequeued.has_value()) {

        ++no_value_count_;

      } else {

        ++consumer_count_;

      }

      ++pull_count_;

    }

  }

  void pullWait() {

    while (consumer_count_.fetch_add(1) < iterations_) {

      auto dequeued = queue_ptr_->waitAndPop();
      ++pull_count_;

    }

  }



  [[nodiscard]] size_t noValueCount() const { return no_value_count_; }
  [[nodiscard]] size_t pushCount() const { return push_count_; }
  [[nodiscard]] size_t pullCount() const { return pull_count_; }
  [[nodiscard]] size_t queueSize() const { return queue_ptr_->size(); }


private:

  const size_t iterations_;
  std::atomic<size_t> no_value_count_{0};
  std::atomic<size_t> pull_count_{0};
  std::atomic<size_t> push_count_{0};
  std::atomic<size_t> producer_count_{0};
  std::atomic<size_t> consumer_count_{0};
  std::unique_ptr<Queue<Type>> queue_ptr_;

};   // ~PushPull

template <template<typename> typename Queue, typename Type> void testBoundedTidal( const std::string queue_name
                                                                                 , const size_t producers
                                                                                 , const size_t consumers
                                                                                 , const size_t iterations) {

  const size_t QUEUE_LOW_TIDE_{2000};
  const size_t QUEUE_HIGH_TIDE_{4000};
  const size_t QUEUE_SAMPLE_FREQ_{100};
  std::vector<std::thread> thread_vector;

  PushPull<Queue, Type> push_pull(iterations);

  ExecEnv::log().info("Queue: {} Begins; Producers: {}, Consumers: {}, iterations: {}", queue_name, producers, consumers, iterations);
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  WorkflowThreads producer_threads(producers);
  producer_threads.enqueueVoid(producers, &PushPull<Queue, Type>::push, &push_pull);

  WorkflowThreads consumer_threads(consumers);
  consumer_threads.enqueueVoid(consumers, &PushPull<Queue, Type>::pullWait, &push_pull);

  producer_threads.joinThreads();
  consumer_threads.joinThreads();

  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  ExecEnv::log().info("Queue: {}, Elapsed Time (ms): {}, Push Count: {}, Pull Count: {}, No Value Pull: {}, Queue Size: {}"
                      , queue_name, elapsed_ms.count(), push_pull.pushCount(), push_pull.pullCount(), push_pull.noValueCount(), push_pull.queueSize());

}


void ExecEnvBGZ::executeApp() {

//  const size_t thread_count{25};
//  testWorkflowReader<BGZStream>("Workflow", thread_count);
//  testWorkflowReader<BGZReader>("Reader", thread_count);

//  return;

  constexpr const size_t producers{5};
  constexpr const size_t consumers{7};
  constexpr const size_t iterations{100000000};

  while (true) {

    testBoundedTidal<QueueTidal, size_t>("QueueTidal<size_t>", producers, consumers, iterations);
    testBoundedTidal<QueueTidal, size_t>("QueueTidal<size_t>", consumers, producers, iterations);

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
