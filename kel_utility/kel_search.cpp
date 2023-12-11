//
// Created by kellerberrin on 11/12/23.
//

#include "kel_search.h"



namespace kel = kellerberrin;


std::vector<kel::OpenRightUnsigned> kel::Search::searchView(const std::regex& search_spec, const std::string_view& sequence_view) {

  std::vector<OpenRightUnsigned> search_matches;

  auto iter_begin = std::regex_iterator(sequence_view.begin(), sequence_view.end(), search_spec);
  auto iter_end = std::cregex_iterator();
  for (auto iter = iter_begin; iter != iter_end; ++iter) {

    auto const& match = *iter;
    if (match.empty()) {

      ExecEnv::log().warn("Unexpected empty search results for regex");
      break;

    }

    search_matches.emplace_back(match.position(), match.position() + match.length());

  }


  return search_matches;

}


std::vector<kel::OpenRightUnsigned> kel::Search::searchView(const std::string_view& search_spec, const std::string_view& sequence_view) {

  return searchView(std::regex(std::string(search_spec), std::regex::icase), sequence_view);

}
