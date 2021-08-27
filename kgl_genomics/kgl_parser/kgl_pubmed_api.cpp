//
// Created by kellerberrin on 20/8/21.
//

#include "kgl_pubmed_api.h"
#include "kel_exec_env.h"
#include "kgl_pubmed_xml_parser.h"  // Do not include in header files, contains 3rd party (rapidxml) references.


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Publication API details.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::LitPublicationMap kgl::PubmedAPIRequester::getCachedPublications(const std::vector<std::string>& pmid_vector) const {


  // For efficiency ensure that all pmids are unique.
  std::set<std::string> unique_pmids;
  for (auto const& pmid : pmid_vector) {

    unique_pmids.insert(pmid);

  }
  std::vector<std::string> unique_pmid_vector;
  for (auto const& pmid : unique_pmids) {

    unique_pmid_vector.push_back(pmid);

  }

  // Only retrieve from the Pubmed API publications not found in the cache.
  std::vector<std::string> unique_uncached_vector;
  for (auto const& pmid : unique_pmid_vector) {

    if (not cached_publications_.contains(pmid)) {

      unique_uncached_vector.push_back(pmid);

    }

  }

  auto api_publications = getAPIPublications(unique_uncached_vector, true); // Write to cache flag is set.

//  ExecEnv::log().info("PubmedAPIRequester::requestCachedPublications; unique publications: {}, found cached: {}, requested Pubmed api: {}",
//                      unique_pmid_vector.size(), (unique_pmid_vector.size()-api_publications.size()), api_publications.size());

  // Combine with the cache.
  cached_publications_.merge(api_publications);

  // Cache should now contain all the required publications.
  LitPublicationMap found_pub_map;
  for (auto const& pmid :  unique_pmid_vector) {

    auto result = cached_publications_.find(pmid);
    if (result != cached_publications_.end()) {

      const auto& [pmid_key, publication] = *result;
      found_pub_map.emplace(pmid, publication);

    } else {

      ExecEnv::log().error("PubmedAPIRequester::requestCachedPublications; requested pmid: {} not found in cache", pmid);

    }

  }

  return found_pub_map;

}


kgl::LitPublicationMap kgl::PubmedAPIRequester::getAPIPublications(const std::vector<std::string>& pmid_vector, bool write_cache) const {

  LitPublicationMap publication_map;

  // Allow trivial requests.
  if (pmid_vector.empty()) {

    return publication_map;

  }

  // For efficiency ensure that all pmids are unique.
  std::set<std::string> unique_pmids;
  for (auto const& pmid : pmid_vector) {

    unique_pmids.insert(pmid);

  }
  std::vector<std::string> unique_pmid_vector;
  for (auto const& pmid : unique_pmids) {

    unique_pmid_vector.push_back(pmid);

  }

  publication_map = getPublicationDetails(unique_pmid_vector, write_cache);

  // Get the citations
  std::vector<std::string> unique_vector;
  for (auto const& pmid : unique_pmids) {

    unique_vector.push_back(pmid);

  }

  // Merge the citations.
  auto cite_map = getCitations(unique_vector, write_cache);
  for (auto& [pmid, publication] : publication_map) {

    auto result = cite_map.find(pmid);
    if (result == cite_map.end()) {

      ExecEnv::log().warn("PubmedRequester::getAPIPublications; no citations found for publication (pmid): {}", pmid);

    } else {

      auto const& [cite_pmid, cite_set] = *result;

      publication.citations(cite_set);

    }

  }

  return publication_map;

}



kgl::LitPublicationMap kgl::PubmedAPIRequester::getPublicationDetails(const std::vector<std::string>& pmid_vector, bool write_cache) const {

  LitPublicationMap publication_map;

  std::chrono::steady_clock::time_point begin_processing_time = std::chrono::steady_clock::now();
  auto begin_batch_time = begin_processing_time;

  std::vector<std::string> batch_pmids;
  size_t api_requests{0};
  size_t batch_requests{0};
  size_t pmid_count{0};

  for (auto const& pmid : pmid_vector) {

    batch_pmids.push_back(pmid);
    if (batch_pmids.size() >= PMID_PER_REQUEST_) {

      std::chrono::steady_clock::time_point begin_api_time = std::chrono::steady_clock::now();
      // Get the results of the api call and add to the citation map.
      auto publication_batch_map = publicationBatch(batch_pmids, write_cache);
      publication_map.merge(publication_batch_map);

      std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
      auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
      auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
      pmid_count += batch_pmids.size();
      ++api_requests;
      ++batch_requests;

      auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
      ExecEnv::log().info("PubmedAPIRequester::getAPIPublications; api call {}ms;  total api calls/elapsed : {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
                          api_duration.count(), api_requests, total_duration.count(), av_api_duration, pmid_count);

      batch_pmids.clear();

      if (batch_requests >= REQUESTS_PER_SECOND_) {

        auto batch_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_batch_time);

        if (batch_milliseconds < BATCH_TIME_INTERVAL_) {

          auto sleep_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(BATCH_TIME_INTERVAL_ - batch_milliseconds);
          std::this_thread::sleep_for(sleep_milliseconds);

        }

        begin_batch_time = end_api_time;
        batch_requests = 0;

      }

    }

  }

  if (not batch_pmids.empty()) {

    std::chrono::steady_clock::time_point begin_api_time = std::chrono::steady_clock::now();
    // Get the results of the api call and add to the citation map.
    auto publication_batch_map = publicationBatch(batch_pmids, write_cache);
    publication_map.merge(publication_batch_map);

    std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
    auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
    pmid_count += batch_pmids.size();
    ++api_requests;

    auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
    ExecEnv::log().info("PubmedAPIRequester::getAPIPublications; api call {}ms;  total api calls/elapsed: {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
                        api_duration.count(), api_requests, total_duration.count(), av_api_duration, pmid_count);

    batch_pmids.clear();

  }

  return publication_map;

}


kgl::LitPublicationMap kgl::PubmedAPIRequester::publicationBatch(const std::vector<std::string>& pmid_vector, bool write_cache) const {

  LitPublicationMap publication_map;
  std::string request_string;
  const std::string id_prefix{"&id="};
  for (auto const& pmid : pmid_vector) {

    request_string += id_prefix + pmid;

  }

  std::string citation_request = PUBMED_ARTICLE_DETAIL_ARGS_ + request_string + API_KEY;

  auto [publication_result, publication_text] = pubmed_rest_api_.synchronousRequest(PUBMED_EFETCH_URL_, citation_request);

  if (not publication_result) {

    ExecEnv::log().error("PubmedAPIRequester::publicationBatch; problem with Pubmed API: {}", publication_text);

  } else {

    auto [parse_result, parsed_map] = ParsePublicationXMLImpl::parsePublicationXML(publication_text);
    if (parse_result and write_cache) {

      if (not pubmed_cache_.writePublicationCache(publication_text)) {

        ExecEnv::log().error("PubmedAPIRequester::publicationBatch; problem writing to publication detail cache");

      }

    }

    publication_map = std::move(parsed_map);

  }

  return publication_map;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::LitCitationMap kgl::PubmedAPIRequester::getCitations(const std::vector<std::string>& pmid_vector, bool write_cache) const {

  return getCitationReference(pmid_vector, PUBMED_ARTICLE_CITEDBY_ARGS_, write_cache);

}

kgl::LitCitationMap kgl::PubmedAPIRequester::getReferences(const std::vector<std::string>& pmid_vector, bool write_cache) const {

  return getCitationReference(pmid_vector, PUBMED_ARTICLE_REFERENCES_ARGS_, write_cache);

}


kgl::LitCitationMap kgl::PubmedAPIRequester::getCitationReference(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args, bool write_cache) const {

  LitCitationMap citation_map;
  std::set<std::string> unique_pmids;

  // Allow trivial requests.
  if (pmid_vector.empty()) {

    return citation_map;

  }

  // For efficiency ensure that all pmids are unique.
  for (auto const& pmid : pmid_vector) {

    unique_pmids.insert(pmid);

  }

  std::chrono::steady_clock::time_point begin_processing_time = std::chrono::steady_clock::now();
  auto begin_batch_time = begin_processing_time;

  std::vector<std::string> batch_pmids;
  size_t api_requests{0};
  size_t batch_requests{0};
  size_t pmid_count{0};

  for (auto const& pmid : unique_pmids) {

    batch_pmids.push_back(pmid);
    if (batch_pmids.size() >= PMID_PER_REQUEST_) {

      std::chrono::steady_clock::time_point begin_api_time = std::chrono::steady_clock::now();
      // Get the results of the api call and add to the citation map.
      auto citation_batch_map = citationBatch(batch_pmids, cite_type_args, write_cache);
      citation_map.merge(citation_batch_map);

      std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
      auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
      auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
      pmid_count += batch_pmids.size();
      ++api_requests;
      ++batch_requests;

      auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
      ExecEnv::log().info("PubmedAPIRequester::getCitationReference; api call {}ms;  total api calls/elapsed: {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
                          api_duration.count(), api_requests, total_duration.count(), av_api_duration, pmid_count);

      batch_pmids.clear();

      if (batch_requests >= REQUESTS_PER_SECOND_) {

        auto batch_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_batch_time);

        if (batch_milliseconds < BATCH_TIME_INTERVAL_) {

          auto sleep_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(BATCH_TIME_INTERVAL_ - batch_milliseconds);
          std::this_thread::sleep_for(sleep_milliseconds);

        }

        begin_batch_time = end_api_time;
        batch_requests = 0;

      }

    }

  }

  if (not batch_pmids.empty()) {

    std::chrono::steady_clock::time_point begin_api_time = std::chrono::steady_clock::now();
    // Get the results of the api call and add to the citation map.
    auto citation_batch_map = citationBatch(batch_pmids, cite_type_args, write_cache);
    citation_map.merge(citation_batch_map);

    std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
    auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
    pmid_count += batch_pmids.size();
    ++api_requests;

    auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
    ExecEnv::log().info("PubmedAPIRequester::getCitationReference; api call {}ms;  total api calls/elapsed: {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
                        api_duration.count(), api_requests, total_duration.count(), av_api_duration, pmid_count);

    batch_pmids.clear();

  }

  return citation_map;

}

kgl::LitCitationMap kgl::PubmedAPIRequester::citationBatch(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args, bool write_cache) const {

  LitCitationMap citation_map;
  std::string request_string;
  const std::string id_prefix{"&id="};
  for (auto const& pmid : pmid_vector) {

    request_string += id_prefix + pmid;

  }

  std::string citation_request = cite_type_args + request_string + API_KEY;

  auto [citation_result, citation_text] = pubmed_rest_api_.synchronousRequest(PUBMED_ELINK_URL_, citation_request);

  if (not citation_result) {

    ExecEnv::log().error("PubmedAPIRequester::getCitationReference; problem with Pubmed API: {}", citation_text);

  } else {

     auto [parse_result, parsed_map] = ParseCitationXMLImpl::parseCitationXML(citation_text);
     if (parse_result and write_cache) {

       if (not pubmed_cache_.writeCitationCache(citation_text)) {

         ExecEnv::log().error("PubmedAPIRequester::getCitationReference; problem writing to citation cache");

       }

     }

     citation_map = std::move(parsed_map);

  }

  return citation_map;

}


