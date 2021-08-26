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
// Stores Pubmed returned XML records in cache files.
// On request, reads and parses the cache files and passes back any requested pubmed record found in the caches.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PubmedAPICache {

public:

  PubmedAPICache() = default;
  ~PubmedAPICache() = default;

  // Empty both cache files.
  [[nodiscard]] bool flushCache() const;

  // Write publication detail XML records to cache.
  [[nodiscard]] bool writeDetailCache(const std::string& xml_cache_records) const;

  // Write publication citation XML records to cache.
  [[nodiscard]] bool writeCitationCache(const std::string& xml_cache_records) const;

  // Unconditionally read all cached publications.
  [[nodiscard]] LitPublicationMap readCachedPublications() const;

  // Sets the cache file location. This is a partial file spec prefixed to the cite and publication files specified below.
  void setCacheFilePrefix(const std::string& cache_file_prefix) const { cache_file_prefix_ = cache_file_prefix; }

private:

  mutable std::string cache_file_prefix_;

  // Literature Cache files.
  const std::string PUBLICATION_CACHE_{"pubmed_pub_cache.xml"};
  const std::string CITATION_CACHE_{"pubmed_cite_cache.xml"};

  // Pseudo XML to delimit the Pubmed XML text.
  const std::string START_CACHE_NODE_{"<CacheBlock Size="};
  constexpr static const char START_CACHE_NODE_END_{'>'};
  const std::string END_CACHE_NODE_{"</CacheBlock>"};

  constexpr static const size_t MAX_CACHE_SIZE_{1000000};
  constexpr static const size_t MIN_CACHE_SIZE_{0};

  // Read and parse publication details
  [[nodiscard]] LitPublicationMap readPublicationCache() const;
  // Read and parse literature citations.
  [[nodiscard]] LitCitationMap readCitationCache() const;
  // Read cached XML text.
  [[nodiscard]] bool readCacheRecord(std::istream& input, std::string& record_string) const;

};



} // namespace.


#endif //KGL_PUBMED_CACHE_H
