//
// Created by kellerberrin on 11/8/21.
//

#include "kgl_pubmed_api.h"
#include "kgl_pubmed_xml_parser.h"
#include "kel_exec_env.h"

#include <set>

namespace kgl = kellerberrin::genome;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Citations/References for each publication.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::LitCitationMap kgl::PubmedRequester::getCitations(const std::vector<std::string>& pmid_vector) const {

  return getCitationReference(pmid_vector, PUBMED_ARTICLE_CITEDBY_ARGS_);

}

kgl::LitCitationMap kgl::PubmedRequester::getReferences(const std::vector<std::string>& pmid_vector) const {

  return getCitationReference(pmid_vector, PUBMED_ARTICLE_REFERENCES_ARGS_);

}


kgl::LitCitationMap kgl::PubmedRequester::getCitationReference(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) const {

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
      auto citation_batch_map = citationBatch(batch_pmids, cite_type_args);
      citation_map.merge(citation_batch_map);

      std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
      auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
      auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
      pmid_count += batch_pmids.size();
      ++api_requests;
      ++batch_requests;

      auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
      ExecEnv::log().info("PubmedRequester::getCitationReference; api call {}ms;  total api calls/elapsed: {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
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
    auto citation_batch_map = citationBatch(batch_pmids, cite_type_args);
    citation_map.merge(citation_batch_map);

    std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
    auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
    pmid_count += batch_pmids.size();
    ++api_requests;

    auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
    ExecEnv::log().info("PubmedRequester::getCitationReference; api call {}ms;  total api calls/elapsed: {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
                        api_duration.count(), api_requests, total_duration.count(), av_api_duration, pmid_count);

    batch_pmids.clear();

  }

  return citation_map;

}

kgl::LitCitationMap kgl::PubmedRequester::citationBatch(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) const {

  LitCitationMap citation_map;
  std::string request_string;
  const std::string id_prefix{"&id="};
  for (auto const& pmid : pmid_vector) {

    request_string += id_prefix + pmid;

  }

  std::string citation_request = cite_type_args + request_string + API_KEY;

  auto [citation_result, citation_text] = pubmed_rest_api_.synchronousRequest(PUBMED_ELINK_URL_, citation_request);

  if (not citation_result) {

    ExecEnv::log().error("PubmedRequester::getCitationReference; problem with Pubmed API: {}", citation_text);

  } else {

    citation_map = ParseCitationXMLImpl::parseCitationXML(citation_text);

  }

  return citation_map;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Publication details.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::LitPublicationMap kgl::PubmedRequester::getPublicationDetails(const std::vector<std::string>& pmid_vector) const {

  LitPublicationMap publication_map;
  std::set<std::string> unique_pmids;

  // Allow trivial requests.
  if (pmid_vector.empty()) {

    return publication_map;

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
      auto publication_batch_map = publicationBatch(batch_pmids);
      publication_map.merge(publication_batch_map);

      std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
      auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
      auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
      pmid_count += batch_pmids.size();
      ++api_requests;
      ++batch_requests;

      auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
      ExecEnv::log().info("PubmedRequester::getPublicationDetails; api call {}ms;  total api calls/elapsed : {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
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
    auto publication_batch_map = publicationBatch(batch_pmids);
    publication_map.merge(publication_batch_map);

    std::chrono::steady_clock::time_point end_api_time = std::chrono::steady_clock::now();
    auto api_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_api_time);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_api_time - begin_processing_time);
    pmid_count += batch_pmids.size();
    ++api_requests;

    auto av_api_duration = static_cast<double>(total_duration.count()) / static_cast<double>(api_requests);
    ExecEnv::log().info("PubmedRequester::getPublicationDetails; api call {}ms;  total api calls/elapsed: {}/{}ms, av. api call {:.1f}ms, pmids retrieved: {}",
                        api_duration.count(), api_requests, total_duration.count(), av_api_duration, pmid_count);

    batch_pmids.clear();

  }

  // Get the citations
  std::vector<std::string> unique_vector;
  for (auto const& pmid : unique_pmids) {

    unique_vector.push_back(pmid);

  }

  // Merge the citations.
  auto cite_map = getCitations(unique_vector);
  for (auto& [pmid, publication] : publication_map) {

    auto result = cite_map.find(pmid);
    if (result == cite_map.end()) {

      ExecEnv::log().warn("PubmedRequester::getPublicationDetails; no citations found for publication (pmid): {}", pmid);

    } else {

      auto const& [cite_pmid, cite_set] = *result;

      publication.citations(cite_set);

    }

  }

  return publication_map;

}


kgl::LitPublicationMap kgl::PubmedRequester::publicationBatch(const std::vector<std::string>& pmid_vector) const {

  LitPublicationMap publication_map;
  std::string request_string;
  const std::string id_prefix{"&id="};
  for (auto const& pmid : pmid_vector) {

    request_string += id_prefix + pmid;

  }

  std::string citation_request = PUBMED_ARTICLE_DETAIL_ARGS_ + request_string + API_KEY;

  auto [publication_result, publication_text] = pubmed_rest_api_.synchronousRequest(PUBMED_EFETCH_URL_, citation_request);

  if (not publication_result) {

    ExecEnv::log().error("PubmedRequester::publicationBatch; problem with Pubmed API: {}", publication_text);

  } else {

    publication_map = ParsePublicationXMLImpl::parsePublicationXML(publication_text);

  }

  return publication_map;

}


