//
// Created by kellerberrin on 9/8/21.
//

#ifndef KEL_REST_API_H
#define KEL_REST_API_H

#include <string>
#include <memory>


namespace kellerberrin {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Object to query Restful HTTP APIs (Curl implementation hidden using the PIMPL idiom).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Forward decl of implementation class.
class CurlRestAPI;

class RestAPI {

public:

  RestAPI();
  ~RestAPI();

  [[nodiscard]] std::pair<bool, std::string> synchronousRequest(const std::string& url, const std::string& raw_request);


private:

  // PIMPL idiom implementation
  std::unique_ptr<CurlRestAPI> restapi_ptr_;

};



} // namespace

#endif //KEL_REST_API_H
