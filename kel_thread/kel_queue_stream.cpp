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
//

#include <iostream>

// Replace the spdlog based messaging system.
enum class MessageType { INFO, WARNING, ERROR};

#define WORKFLOW_STAND_ALONE 1
#ifdef WORKFLOW_STAND_ALONE

void workflowStreamOut(MessageType type, const std::string& message) {

  switch (type) {

    case MessageType::INFO:
      std::clog << "INFO - ";
      break;

    case MessageType::WARNING:
      std::clog << "WARNING - ";
      break;

    case MessageType::ERROR:
      std::clog << "ERROR - ";
      break;

  }

  std::clog << message << std::endl;

}

#else

#include "kel_exec_env.h"

void workflowStreamOut(MessageType type, const std::string& message) {

  switch (type) {

    case MessageType::INFO:
      kellerberrin::ExecEnv::log().info(message);
      break;

    case MessageType::WARNING:
      kellerberrin::ExecEnv::log().warn(message);
      break;

    case MessageType::ERROR:
      kellerberrin::ExecEnv::log().error(message);
      break;

  }

}

#endif


