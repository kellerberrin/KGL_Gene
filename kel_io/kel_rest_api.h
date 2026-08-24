//
// Created by kellerberrin on 9/8/21.
//

#ifndef KEL_REST_API_H
#define KEL_REST_API_H

#include <string>
#include <memory>
#include <utility>


namespace kellerberrin {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Object to query Restful HTTP APIs (Curl implementation hidden using the PIMPL idiom).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Forward decl of implementation class.
class CurlRestAPI;

/// Synchronous REST API client using libcurl (implementation hidden via PIMPL idiom).
/// Each instance owns its own CURL handle. The global curl library is initialised once per process.
class RestAPI {

public:

  RestAPI();
  ~RestAPI();

  /// Performs a synchronous POST request to the given URL with the given raw request body.
  /// Returns {true, response_body} on success, or {false, error_or_partial_response} on failure.
  [[nodiscard]] std::pair<bool, std::string> synchronousRequest(const std::string& url, const std::string& raw_request);


private:

  // PIMPL idiom implementation
  std::unique_ptr<CurlRestAPI> restapi_ptr_;

};



} // namespace

#endif //KEL_REST_API_H
