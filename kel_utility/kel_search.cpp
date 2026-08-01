//
// Created by kellerberrin on 11/12/23.
//

#include "kel_search.h"



namespace kel = kellerberrin;


std::vector<kel::OpenRightUnsigned> kel::Search::searchView(const std::regex& search_spec, const std::string_view& sequence_view) {

  std::vector<OpenRightUnsigned> search_matches;

  using RegexIter = std::regex_iterator<std::string_view::const_iterator>;
  RegexIter iter_begin(sequence_view.begin(), sequence_view.end(), search_spec);
  RegexIter iter_end; // Default constructor creates an end() guard.
  for (auto iter = iter_begin; iter != iter_end; ++iter) {

    auto const& match = *iter;
    if (match.empty()) {

      ExecEnv::log().warn("Unexpected empty search results for regex");
      continue;

    }
    if (match.position() < 0) {

      ExecEnv::log().error("Found -ve position in regex sequence search");
      continue;

    }

    search_matches.emplace_back(match.position(), match.position() + match.length());

  }


  return search_matches;

}


std::vector<kel::OpenRightUnsigned> kel::Search::searchView(const std::string_view& search_spec, const std::string_view& sequence_view) {

  try {

    auto regex_query = std::regex(std::string(search_spec), std::regex::icase);
    auto search_results = searchView(regex_query, sequence_view);
    return search_results;

  } catch (const std::regex_error& e) {

  ExecEnv::log().error("Invalid regex text: {}, error: {}", search_spec, e.what());
  return {};

  } catch (const std::exception& e) {

  ExecEnv::log().error("Unexpected error searching sequence for regex text: {}, error: {}", search_spec, e.what());
  return {};

  }

}
