//
// Created by kellerberrin on 11/8/21.
//

#ifndef KGL_PUBMED_RESOURCE_H
#define KGL_PUBMED_RESOURCE_H

#include "kgl_pubmed_api.h"
#include "kgl_resource_db.h"

#include "kel_utility.h"

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

  PubmedRequester(const std::string& identifier, const std::string& cache_file_spec)
  : ResourceBase(identifier) , cache_file_spec_(cache_file_spec) {}
  ~PubmedRequester() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::PUBMED_API; }

  // Only the data for unique pmids is returned, for convenience the requesting array of pmids can contain non-unique pmids.
  // Large requests are automatically throttled to Pubmed requirements (max 10 api calls a second).
  // Citation and reference information. No cache records are written or read.
  // Pubmed publication caching is not used.
  [[nodiscard]] LitPublicationMap getPublications(const std::vector<std::string>& pmid_vector) const { return pubmed_rest_api_.getPublications(pmid_vector, false); }

  // Same functionality as above but checks if the publications are held on a disk cache before sending API requests
  // to Pubmed. Any records not found in the cache are requested using the Pubmed API and then written to the cache.
  [[nodiscard]] LitPublicationMap getCachedPublications(const std::vector<std::string>& pmid_vector) const { return pubmed_rest_api_.getCachedPublications(pmid_vector); }

  // Empty the cache files.
  [[nodiscard]] bool flushCache() const { return pubmed_rest_api_.flushCache(); }

  // Return the resource cache file spec;
  [[nodiscard]] const std::string& cacheFileSpec() const { return cache_file_spec_; }

  // Set the partial file spec for the cache files.
  void setWorkDirectory(const std::string& directory) const { pubmed_rest_api_.setCacheFile(Utility::filePath(cache_file_spec_, directory)); }

private:

  const PubmedAPIRequester pubmed_rest_api_;
  const std::string cache_file_spec_;

};





} // namespace


#endif //KGL_PUBMED_RESOURCE_H
