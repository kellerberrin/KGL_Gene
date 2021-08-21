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


class pubMedAPICache {

public:

  pubMedAPICache() = default;
  ~pubMedAPICache() = default;

  // Empty both cache files.
  [[nodiscard]] bool flushCache() const;

  // Write publication detail XML records to cache.
  [[nodiscard]] bool writeDetailCache(const std::string& xml_cache_records) const;

  // Write publication citation XML records to cache.
  [[nodiscard]] bool writeCitationCache(const std::string& xml_cache_records) const;

  // Return any requested publications found in the caches.
  [[nodiscard]] LitPublicationMap getCachedPublications(const std::vector<std::string>& pmid_vector) const;

  // Sets the cache file name, This is a partial file spec.
  void setCacheFile(const std::string& cache_file) const { cache_file_name_ = cache_file; }

private:

  mutable std::string cache_file_name_;

  constexpr static const char* PUBLICATION_CACHE_{"pub_cache.xml"};
  constexpr static const char* CITATION_CACHE_{"cite_cache.xml"};

  const std::string START_CACHE_NODE_{"<CacheBlock Size="};
  constexpr static const char START_CACHE_NODE_END_{'>'};
  const std::string END_CACHE_NODE_{"</CacheBlock>"};

  constexpr static const size_t MAX_CACHE_SIZE_{1000000};
  constexpr static const size_t MIN_CACHE_SIZE{0};

  [[nodiscard]] LitPublicationMap readParseCachedPublications() const;
  [[nodiscard]] LitCitationMap readParseCachedCitations() const;
  [[nodiscard]] bool readCacheRecord(std::istream& input, std::string& record_string) const;

};



} // namespace.


#endif //KGL_PUBMED_CACHE_H
