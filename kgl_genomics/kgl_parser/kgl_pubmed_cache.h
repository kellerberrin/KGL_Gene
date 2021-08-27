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

  PubmedAPICache(const std::string& publication_cache_file, const std::string& citation_cache_file)
  : publication_cache_file_(publication_cache_file), citation_cache_file_(citation_cache_file) {}
  ~PubmedAPICache() = default;

  // Write publication detail XML records to cache.
  [[nodiscard]] bool writePublicationCache(const std::string& xml_cache_record) const;

  // Write publication citation XML records to cache.
  [[nodiscard]] bool writeCitationCache(const std::string& xml_cache_records) const;

  // Unconditionally read all cached publications.
  [[nodiscard]] LitPublicationMap readCachedPublications() const;


private:

  // Literature Cache files.
  const std::string publication_cache_file_;
  const std::string citation_cache_file_;

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
