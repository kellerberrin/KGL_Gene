//
// Created by kellerberrin on 11/8/21.
//

#include "kgl_pubmed_api.h"
#include "kel_exec_env.h"

#include <rapidxml.h>

#include <set>

namespace kgl = kellerberrin::genome;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Citations/References for each publication.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::LitCitationMap kgl::PubmedRequester::getCitations(const std::vector<std::string>& pmid_vector) {

  return getCitationArgs(pmid_vector, PUBMED_ARTICLE_CITEDBY_ARGS_);

}

kgl::LitCitationMap kgl::PubmedRequester::getReferences(const std::vector<std::string>& pmid_vector) {

  return getCitationArgs(pmid_vector, PUBMED_ARTICLE_REFERENCES_ARGS_);

}


kgl::LitCitationMap kgl::PubmedRequester::getCitationArgs(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) {

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

kgl::LitCitationMap kgl::PubmedRequester::citationBatch(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) {

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

    citation_map = parseCitationXML(citation_text);

  }

  auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch - begin_batch);

  ExecEnv::log().info("Batch Execution time: {} (ms)", batch_duration.count());

  return citation_map;

}


kgl::LitCitationMap kgl::PubmedRequester::parseCitationXML(const std::string& citation_text) {

  LitCitationMap citation_map;

  // All this raw pointer stuff is very nasty, is rapidxml the right library?
  rapidxml::xml_document<> doc;
  std::unique_ptr<char[]> text_buffer_ptr(std::make_unique<char[]>(citation_text.size() + 1));
  std::memcpy(text_buffer_ptr.get(), citation_text.c_str(), citation_text.size() + 1);

  // Parse the buffer using the xml file parsing library into doc
  doc.parse<0>(text_buffer_ptr.get());

  rapidxml::xml_node<>* root_node = doc.first_node(CITATION_ROOT_NODE_.c_str());

  if (root_node != nullptr) {

    rapidxml::xml_node<> * citation_node = root_node->first_node(CITATION_RECORD_.c_str());
    while (citation_node != nullptr)
    {

      auto citation_pmid = citation_node->first_node(CITATION_PMID_.c_str());
      if (citation_pmid != nullptr) {

        auto pmid_node = citation_pmid->first_node();
        if (pmid_node != nullptr) {

          std::vector<std::string> citation_pmid_vector;
          auto link_db = citation_node->first_node(CITATION_LINK_DB_.c_str());
          if (link_db != nullptr) {

            auto link_node = link_db->first_node(CITATION_LINK_SET_.c_str());
            while (link_node != nullptr) {

              auto cite_pmid_node = link_node->first_node();
              if (cite_pmid_node != nullptr) {

                citation_pmid_vector.emplace_back(cite_pmid_node->value());

              } else {

                ExecEnv::log().error("PubmedRequester::parseCitationXML; No citation node found");

              }

              link_node = link_node->next_sibling();

            }

          }

          auto [iter, result] = citation_map.try_emplace(pmid_node->value(), citation_pmid_vector);
          if (not result) {

            ExecEnv::log().error("PubmedRequester::parseCitationXML; expected duplicate for pmid: {}", pmid_node->value());

          }

        } else {

          ExecEnv::log().error("PubmedRequester::parseCitationXML; Pmid Id Attribute: {} does not exist", CITATION_PMID_);

        }

      } else {

        ExecEnv::log().error("PubmedRequester::parseCitationXML; Pmid Id Attribute: {} does not exist", CITATION_PMID_);

      }

      citation_node = citation_node->next_sibling();

    }

  } else  {

    ExecEnv::log().error("PubmedRequester::parseCitationXML; Citation root node: {} is null", CITATION_ROOT_NODE_);

  }

  return citation_map;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Publication details.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::LitPublicationMap kgl::PubmedRequester::getPublicationDetails(const std::vector<std::string>& pmid_vector) {

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

  return publication_map;

}


kgl::LitPublicationMap kgl::PubmedRequester::publicationBatch(const std::vector<std::string>& pmid_vector) {

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

    ExecEnv::log().info("Publication text: {}", publication_text);
    publication_map = parsePublicationXML(publication_text);

  }

  auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch - begin_batch);

  ExecEnv::log().info("Batch Execution time: {} (ms)", batch_duration.count());

  return publication_map;

}


kgl::LitPublicationMap kgl::PubmedRequester::parsePublicationXML(const std::string& publication_xml_text) {

  LitPublicationMap publication_map;

  // Normally I avoid try/catch blocks but these XML documents are just too complex for if/then error handling.
  try {

    rapidxml::xml_document<> doc;
    std::unique_ptr<char[]> text_buffer_ptr(std::make_unique<char[]>(publication_xml_text.size() + 1));
    std::memcpy(text_buffer_ptr.get(), publication_xml_text.c_str(), publication_xml_text.size() + 1);

    // Parse the buffer using the xml file parsing library into doc
    doc.parse<0>(text_buffer_ptr.get());

    rapidxml::xml_node<>* root_node = doc.first_node(PUBLICATION_ROOT_NODE_.c_str());
    if (root_node == nullptr) {

      std::string error_message = "Publication Root Node: " + PUBLICATION_ROOT_NODE_ + " not found";
      throw std::runtime_error(error_message);

    }

    // Loop through the articles.
    size_t article_count{0};
    rapidxml::xml_node<> * article_node = root_node->first_node(PUBLICATION_NODE_.c_str());
    while (article_node != nullptr) {

      ++article_count;

      article_node = article_node->next_sibling();

    }

    ExecEnv::log().info("PubmedRequester::parsePublicationXML; article count {}", article_count);

  } catch(std::exception& e) {

    ExecEnv::log().error("PubmedRequester::parsePublicationXML; error parsing Pubmed publication XML, error: {}", e.what());
    publication_map.clear();

  }

  return publication_map;

}

