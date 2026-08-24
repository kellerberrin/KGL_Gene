//
// Created by kellerberrin on 9/8/21.
//

#include "kel_rest_api.h"
#include "kel_exec_env.h"


#include <curl/curl.h>

#include <mutex>
#include <string>
#include <vector>


namespace kellerberrin { // namespace



class CurlRestAPI {

public:

  CurlRestAPI();
  ~CurlRestAPI();

  [[nodiscard]] std::pair<bool, std::string> synchronousRequest(const std::string& url, const std::string& raw_request);


private:

  CURL* lib_curl_handle_;
  // Serialize concurrent requests on the same CURL handle. A single handle is not thread-safe.
  std::mutex request_mutex_;

  // Response callback.
  static size_t append(char* buffer, size_t item_size, size_t num_items, void* user_data);

};

// curl_global_init()/curl_global_cleanup() must be called exactly once and are not thread-safe.
// Guard the global library initialisation with a once_flag and schedule cleanup at process exit.
namespace {

void ensureCurlGlobalInit() {

  static std::once_flag init_flag;
  std::call_once(init_flag, [] {
    curl_global_init(CURL_GLOBAL_DEFAULT);
    std::atexit([] { curl_global_cleanup(); });
  });

}

} // anonymous namespace


CurlRestAPI::CurlRestAPI() {

  ensureCurlGlobalInit();
  lib_curl_handle_ = curl_easy_init();
  if (lib_curl_handle_ == nullptr) {

    ExecEnv::log().error("CurlRestAPI::CurlRestAPI; Unable to initialise Curl Library");

  }

}


CurlRestAPI::~CurlRestAPI() {

  if (lib_curl_handle_ != nullptr) {

    curl_easy_cleanup(lib_curl_handle_);

  }

}

std::pair<bool, std::string> CurlRestAPI::synchronousRequest(const std::string& url, const std::string& raw_request) {

  std::lock_guard<std::mutex> lock(request_mutex_);

  if (lib_curl_handle_ == nullptr) {

    return {false, "CurlRestAPI::synchronousRequest; Invalid Curl Handle"};

  }

  // Use a per-request response buffer passed through CURLOPT_WRITEDATA. This makes the
  // callback instance-safe and thread-safe (no shared state).
  std::string response_buffer;

  auto setopt = [&](auto option, auto value) -> bool {

    const auto result = curl_easy_setopt(lib_curl_handle_, option, value);
    if (result != CURLE_OK) {

      ExecEnv::log().error("CurlRestAPI::synchronousRequest; curl_easy_setopt failed: {}", curl_easy_strerror(result));
      return false;

    }
    return true;

  };

  if (not setopt(CURLOPT_POSTFIELDS, raw_request.c_str())) return {false, response_buffer};
  if (not setopt(CURLOPT_URL, url.c_str())) return {false, response_buffer};
  // Register the callback.
  if (not setopt(CURLOPT_WRITEFUNCTION, &CurlRestAPI::append)) return {false, response_buffer};
  if (not setopt(CURLOPT_WRITEDATA, &response_buffer)) return {false, response_buffer};

  // The call is synchronous (blocking).
  const auto result = curl_easy_perform(lib_curl_handle_);

  // Error check.
  if (result != CURLE_OK) {

    ExecEnv::log().error("CurlRestAPI::synchronousRequest; error making curl request: {}", curl_easy_strerror(result));
    return {false, response_buffer};

  }

  return {true, response_buffer};

}


// Static callback function.
size_t CurlRestAPI::append(char* data, size_t buffer_size, size_t num_items, void* user_data) {

  auto* response_buffer = static_cast<std::string*>(user_data);
  if (response_buffer == nullptr) {

    return 0;

  }

  try {

    response_buffer->append(data, buffer_size * num_items);

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("CurlRestAPI::append; exception appending response data: {}", e.what());
    return 0;  // Signal error to curl.

  }
  return buffer_size * num_items;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Just re-directed PIMPL calls.

RestAPI::RestAPI() {

  restapi_ptr_ = std::make_unique<CurlRestAPI>();

}

// Must be explicitly declared.
RestAPI::~RestAPI() {}


std::pair<bool, std::string> RestAPI::synchronousRequest(const std::string& url, const std::string& raw_request) {

  return restapi_ptr_->synchronousRequest(url, raw_request);

}


} // namespace
