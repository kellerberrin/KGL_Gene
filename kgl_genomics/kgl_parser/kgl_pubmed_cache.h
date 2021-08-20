//
// Created by kellerberrin on 20/8/21.
//

#ifndef KGL_PUBMED_CACHE_H
#define KGL_PUBMED_CACHE_H


#include "kgl_pubmed.h"
#include "kel_rest_api.h"
#include "kgl_resource_db.h"

#include <string>
#include <memory>
#include <map>
#include <set>
#include <vector>
#include <chrono>



namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////?////////////////////////////////////////////////////
//
// Pubmed API stores pubmed returned XML records in a file.
// Loads and parses the file and passes back any requested pubmed records
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class pubMedCachedAPI {

public:

  pubMedCachedAPI() = default;
  ~pubMedCachedAPI() = default;


private:





};



} // namespace.


#endif //KGL_PUBMED_CACHE_H
