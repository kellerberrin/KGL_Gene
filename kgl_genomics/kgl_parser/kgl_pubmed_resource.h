//
// Created by kellerberrin on 11/8/21.
//

#ifndef KGL_PUBMED_RESOURCE_H
#define KGL_PUBMED_RESOURCE_H

#include "kgl_pubmed_api.h"
#include "kgl_resource_db.h"

#include <string>
#include <memory>
#include <map>
#include <set>
#include <vector>
#include <chrono>


namespace kellerberrin::genome {   //  organization level namespace



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Request the pubmed literature. The high-level resource used by packages requiring this functionality.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PubmedRequester : public ResourceBase {

public:

  explicit PubmedRequester(const std::string& identifier, const std::string& cache_file_name)
  : ResourceBase(identifier) , cache_file_name_(cache_file_name) {}
  ~PubmedRequester() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::PUBMED_API;  }

  // Only the data for unique pmids is returned, for convenience the requesting array of pmids can contain non-unique pmids.
  // Large requests are automatically throttled to Pubmed requirements (max 10 api calls a second).
  // Citation and reference information.
  [[nodiscard]] LitPublicationMap getPublicationDetails(const std::vector<std::string>& pmid_vector) const;

private:

  const PubmedAPIRequester pubmed_rest_api_;
  const std::string cache_file_name_;

};





} // namespace


#endif //KGL_PUBMED_RESOURCE_H
