//
// Created by kellerberrin on 11/8/21.
//

#ifndef KGL_PUBMED_RESOURCE_H
#define KGL_PUBMED_RESOURCE_H

#include "kgl_pubmed_api.h"
#include "kgl_runtime_resource.h"

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

  PubmedRequester(const std::string& identifier, const std::string& publication_cache_file, const std::string& citation_cache_file)
  : ResourceBase(identifier), pubmed_rest_api_(publication_cache_file, citation_cache_file) {}
  ~PubmedRequester() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::PUBMED_API; }

  // Returns all cached publications. Note that this assumes that all required publications have already been downloaded and are cached.
  [[nodiscard]] const LitPublicationMap& getAllCachedPublications() const { return pubmed_rest_api_.getAllCachedPublications(); }
  // Same functionality as above but checks if the publications are held on a disk/memory cache before sending API requests to Pubmed.
  // Any records not found in the cache are requested using the Pubmed API and then written to the disk/memory cache.
  [[nodiscard]] LitPublicationMap getCachedPublications(const std::vector<std::string>& pmid_vector) const { return pubmed_rest_api_.getCachedPublications(pmid_vector); }
  // Same as above but with a set of pmids.
  [[nodiscard]] LitPublicationMap getCachedPublications(const std::set<std::string>& pmid_set) const { return pubmed_rest_api_.getCachedPublications(pmid_set); }

private:

  const PubmedAPIRequester pubmed_rest_api_;

};





} // namespace


#endif //KGL_PUBMED_RESOURCE_H
