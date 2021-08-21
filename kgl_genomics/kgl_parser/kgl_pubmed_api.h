//
// Created by kellerberrin on 20/8/21.
//

#ifndef KGL_PUBMED_API_H
#define KGL_PUBMED_API_H

#include "kgl_pubmed.h"
#include "kel_rest_api.h"
#include "kgl_pubmed_cache.h"

#include <string>
#include <memory>
#include <map>
#include <set>
#include <vector>
#include <chrono>


namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////?////////////////////////////////////////////////////
//
// Pubmed API requestor. This is a resource that can be requested by analysis packages at runtime.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PubmedAPIRequester {

public:

  PubmedAPIRequester() = default;
  ~PubmedAPIRequester() = default;

  // Returns cited publications after making calls for publication details and publication citations.
  // If write cache flag is set to true then XML records are written to the cache files, no cached entries are returned.
  [[nodiscard]] LitPublicationMap getPublications(const std::vector<std::string>& pmid_vector, bool write_cache) const;

  // Request Pubmed publications, use cached results when available, else request from Pubmed using the API.
  // Always writes correctly parsed API results to cache.
  [[nodiscard]] LitPublicationMap getCachedPublications(const std::vector<std::string>& pmid_vector) const;

  // Empty the Pubmed cache files.
  [[nodiscard]] bool flushCache() const { return pubmed_cache_.flushCache(); }

  // Set the cache file.
  void setCacheFile(const std::string& cache_file) const { pubmed_cache_.setCacheFile(cache_file); }

  // Make separate api calls.
  [[nodiscard]] LitPublicationMap getPublicationDetails(const std::vector<std::string>& pmid_vector, bool write_cache) const;
  [[nodiscard]] LitCitationMap getCitations(const std::vector<std::string>& pmid_vector, bool write_cache) const;
  [[nodiscard]] LitCitationMap getReferences(const std::vector<std::string>& pmid_vector, bool write_cache) const;

private:

  mutable RestAPI pubmed_rest_api_;
  const pubMedAPICache pubmed_cache_;

  // Kellerberrin Pubmed API key.
  const std::string API_KEY {"&api_key=8cd3dde4cbf1eeb71b5ae469ae8a99247609"};
  // Pubmed resource constraints.
  constexpr const static size_t REQUESTS_PER_SECOND_{10};
  constexpr const static size_t PMID_PER_REQUEST_{10};
  constexpr const static std::chrono::milliseconds BATCH_TIME_INTERVAL_{1000};

  // The citation/reference constants.
  const std::string PUBMED_ELINK_URL_{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"};
  const std::string PUBMED_ARTICLE_CITEDBY_ARGS_{"dbfrom=pubmed&linkname=pubmed_pubmed_citedin"};
  const std::string PUBMED_ARTICLE_REFERENCES_ARGS_{"dbfrom=pubmed&linkname=pubmed_pubmed_refs"};

  [[nodiscard]] LitCitationMap getCitationReference(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args, bool write_cache) const;
  [[nodiscard]] LitCitationMap citationBatch(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args, bool write_cache) const;

  // The publication detail constants.
  const std::string PUBMED_EFETCH_URL_{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"};
  const std::string PUBMED_ARTICLE_DETAIL_ARGS_{"db=pubmed&retmode=xml"};

  [[nodiscard]] LitPublicationMap publicationBatch(const std::vector<std::string>& pmid_vector, bool write_cache) const;

};




} // namespace



#endif // KGL_PUBMED_API_H
