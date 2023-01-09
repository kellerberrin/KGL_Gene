// Copyright 2023 Kellerberrin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//

#ifndef KEL_WORKFLOW_UNITTEST_H
#define KEL_WORKFLOW_UNITTEST_H


#include "kel_workflow_async.h"
#include "kel_workflow_sync.h"
#include "kel_exec_env_app.h"
#include "kel_logging.h"

#include <map>


namespace kellerberrin {  //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Unit test object tests simple asynchronous work queues.
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//template<typename T> using AsynchQueueType = MtQueue<T>;
template <typename T> using AsyncQueueType = BoundedMtQueue<T>;

using AsynchInQueue = WorkflowAsyncQueue<InputType, AsyncQueueType>;
using AsynchMedQueue = WorkflowAsyncQueue<IntermediateType, AsyncQueueType>;
using AsynchOutQueue = WorkflowAsyncQueue<OutputType, AsyncQueueType>;


class AsynchQueueUnitTest {

public:

  AsynchQueueUnitTest() = default;
  ~AsynchQueueUnitTest() = default;

  static void asynchMoveable();

  void asynchInputWork( std::shared_ptr<AsynchMedQueue> med_queue
                       , InputType input_item);

  void asynchIntermediateWork( std::shared_ptr<const AsynchMedQueue> med_queue
                              , std::shared_ptr<AsynchOutQueue> output_queue
                              , IntermediateType intermediate_item);

  void asynchOutputWork(std::shared_ptr<const AsynchOutQueue> output_queue, OutputType output_item);

  void pushInput(std::shared_ptr<AsynchInQueue> input_queue);


  inline static const size_t iterations{1000000};
  inline static const size_t report_iterations{100000};
  inline static const size_t work_iterations{1000000}; // work for each thread

  inline static const size_t high_tide{5000};
  inline static const size_t low_tide{1000};
  inline static const size_t mon_freq_ms{100};

  inline static const size_t input_thread_count{50};
  inline static const size_t intermediate_thread_count{50};
  inline static const size_t output_thread_count{50};

private:

  std::atomic<size_t> input_count_{0};
  std::atomic<size_t> intermediate_count_{0};
  std::atomic<size_t> output_count_{0};
  size_t check_sum_{0};


};




} // end namespace.


#endif //KEL_WORKFLOW_UNITTEST_H
