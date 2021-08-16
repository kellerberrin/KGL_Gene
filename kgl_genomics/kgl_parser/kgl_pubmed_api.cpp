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

  return getCitationArgs(pmid_vector, PUBMED_ARTICLE_CITEDBY_ARGS_);

}

kgl::LitCitationMap kgl::PubmedRequester::getReferences(const std::vector<std::string>& pmid_vector) const {

  return getCitationArgs(pmid_vector, PUBMED_ARTICLE_REFERENCES_ARGS_);

}


kgl::LitCitationMap kgl::PubmedRequester::getCitationArgs(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) const {

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

  std::vector<std::string> batch_pmids;
  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();
  size_t time_period_requests{0};
  for (auto const& pmid : unique_pmids) {

    batch_pmids.push_back(pmid);
    if (batch_pmids.size() >= PMID_PER_REQUEST_) {

      // Get api results.
      citation_map.merge(citationBatch(batch_pmids, cite_type_args));
      batch_pmids.clear();
      ++time_period_requests;
      if (time_period_requests >= REQUESTS_PER_SECOND_) {

        std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();
        auto batch_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch-begin_batch);

        if (batch_milliseconds < BATCH_TIME_INTERVAL_) {

          auto sleep_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(BATCH_TIME_INTERVAL_ - batch_milliseconds);
          std::this_thread::sleep_for(sleep_milliseconds);

        }

        time_period_requests = 0;
        begin_batch = std::chrono::steady_clock::now();

      }

    }

  }

  if (not batch_pmids.empty()) {

    citation_map.merge(citationBatch(batch_pmids, cite_type_args));
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

  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();

  auto [citation_result, citation_text] = pubmed_rest_api_.synchronousRequest(PUBMED_ELINK_URL_, citation_request);

  std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();

  if (not citation_result) {

    ExecEnv::log().error("PubmedRequester::getCitationArgs; problem with Pubmed API: {}", citation_text);

  } else {

    citation_map = ParseCitationXMLImpl::parseCitationXML(citation_text);

  }

  auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch - begin_batch);

  ExecEnv::log().info("Batch Execution time: {} (ms)", batch_duration.count());

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

  std::vector<std::string> batch_pmids;
  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();
  size_t time_period_requests{0};
  for (auto const& pmid : unique_pmids) {

    batch_pmids.push_back(pmid);
    if (batch_pmids.size() >= PMID_PER_REQUEST_) {

      // Get api results.
      publication_map.merge(publicationBatch(batch_pmids));
      batch_pmids.clear();
      ++time_period_requests;
      if (time_period_requests >= REQUESTS_PER_SECOND_) {

        std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();
        auto batch_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch-begin_batch);

        if (batch_milliseconds < BATCH_TIME_INTERVAL_) {

          auto sleep_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(BATCH_TIME_INTERVAL_ - batch_milliseconds);
          std::this_thread::sleep_for(sleep_milliseconds);

        }

        time_period_requests = 0;
        begin_batch = std::chrono::steady_clock::now();

      }

    }

  }

  if (not batch_pmids.empty()) {

    publication_map.merge(publicationBatch(batch_pmids));
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

  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();

  auto [publication_result, publication_text] = pubmed_rest_api_.synchronousRequest(PUBMED_EFETCH_URL_, citation_request);

  std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();

  if (not publication_result) {

    ExecEnv::log().error("PubmedRequester::publicationBatch; problem with Pubmed API: {}", publication_text);

  } else {

    publication_map = ParsePublicationXMLImpl::parsePublicationXML(publication_text);

  }

  auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch - begin_batch);

  ExecEnv::log().info("Batch Execution time: {} (ms)", batch_duration.count());

  return publication_map;

}


