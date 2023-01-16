//
// Created by kellerberrin on 16/01/23.
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


