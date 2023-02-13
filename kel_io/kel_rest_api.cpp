//
// Created by kellerberrin on 9/8/21.
//

#include "kel_rest_api.h"
#include "kel_exec_env.h"


#include <curl/curl.h>


namespace kellerberrin { // namespace



class CurlRestAPI {

public:

  CurlRestAPI();
  ~CurlRestAPI();

  [[nodiscard]] std::pair<bool, std::string> synchronousRequest(const std::string& url, const std::string& raw_request);


private:

  CURL* lib_curl_handle_;
  inline static std::vector<std::string> request_responses_;

  // Response callback
  static size_t requestCallback(char *buffer, size_t item_size, size_t num_items, void *ignored);

};

CurlRestAPI::CurlRestAPI() {

  curl_global_init(CURL_GLOBAL_DEFAULT);
  lib_curl_handle_ = curl_easy_init();
  if (lib_curl_handle_ == nullptr) {

    ExecEnv::log().error("CurlRestAPI::CurlRestAPI; Unable to initialise Curl Library");

  }

}


CurlRestAPI::~CurlRestAPI() {

  if (lib_curl_handle_ != nullptr) {

    curl_easy_cleanup(lib_curl_handle_);

  }
  curl_global_cleanup();

}

std::pair<bool, std::string> CurlRestAPI::synchronousRequest(const std::string& url, const std::string& raw_request) {

  if (lib_curl_handle_ == nullptr) {

    return {false, "CurlRestAPI::synchronousRequest; Invalid Curl Handle"};

  }

  request_responses_.clear();
  curl_easy_setopt(lib_curl_handle_, CURLOPT_POSTFIELDS, raw_request.c_str());
  curl_easy_setopt(lib_curl_handle_, CURLOPT_URL, url.c_str());
  // Register the callback.
  curl_easy_setopt(lib_curl_handle_, CURLOPT_WRITEFUNCTION, &CurlRestAPI::requestCallback);

  // The call is synchronous (blocking).
  auto result = curl_easy_perform(lib_curl_handle_);

  // Error check
  if (result != CURLE_OK) {

    ExecEnv::log().error("urlRestAPI::synchronousRequest; error making curl request: {}", curl_easy_strerror(result));
    // Unpack the response for any error message.
    std::string error_response;
    for (auto const& response : request_responses_) {

      error_response += response;

    }

    return { false, error_response};

  }

  std::string success_response;
  for (auto const& response : request_responses_) {

    success_response += response;

  }

  return {true, success_response};

}


// Static callback function
size_t CurlRestAPI::requestCallback(char *buffer, size_t item_size, size_t num_items, void * /*ignored*/ )
{

  size_t bytes_returned = item_size * num_items;
  request_responses_.emplace_back(std::string(buffer, bytes_returned));

  return bytes_returned;

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

