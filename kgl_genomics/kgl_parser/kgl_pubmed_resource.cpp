//
// Created by kellerberrin on 11/8/21.
//

#include "kgl_pubmed_resource.h"
#include "kel_exec_env.h"


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Publication details.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::LitPublicationMap kgl::PubmedRequester::getPublicationDetails(const std::vector<std::string>& pmid_vector) const {

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

  publication_map = pubmed_rest_api_.getPublicationDetails(unique_pmid_vector);

  // Get the citations
  std::vector<std::string> unique_vector;
  for (auto const& pmid : unique_pmids) {

    unique_vector.push_back(pmid);

  }

  // Merge the citations.
  auto cite_map = pubmed_rest_api_.getCitations(unique_vector);
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



