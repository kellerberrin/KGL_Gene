//
// Created by kellerberrin on 11/12/23.
//

#include "kel_search.h"



namespace kel = kellerberrin;


#include <cstddef>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

namespace kel = kellerberrin;

std::vector<kel::OpenRightUnsigned>
kel::Search::searchView(const std::regex& search_spec, std::string_view sequence_view) {

  std::vector<OpenRightUnsigned> search_matches;

  const char* const first = sequence_view.data();
  const char* const last  = first + sequence_view.size();

  using RegexIter = std::cregex_iterator;
  RegexIter iter_begin(first, last, search_spec);
  RegexIter iter_end; // Default-constructed end sentinel.

  for (auto iter = iter_begin; iter != iter_end; ++iter) {
    const auto& match = *iter;

    if (match.empty()) {

      ExecEnv::log().warn("Unexpected empty regex match encountered; skipping.");
      continue;

    }

    const auto pos = match.position();
    const auto len = static_cast<std::size_t>(match.length());

    if (pos < 0) {

      ExecEnv::log().error("Invalid negative regex match position ({}); skipping.", pos);
      continue;

    }

    const std::size_t start = static_cast<std::size_t>(pos);
    search_matches.emplace_back(start, start + len);

  }

  return search_matches;

}

std::vector<kel::OpenRightUnsigned>
kel::Search::searchView(std::string_view search_spec, std::string_view sequence_view) {

  try {
    // Case-insensitive search by default.
    const std::regex regex_query{std::string(search_spec), std::regex::icase};
    return searchView(regex_query, sequence_view);

  } catch (const std::regex_error& e) {

    ExecEnv::log().error("Invalid regex text '{}': {}", search_spec, e.what());
    return {};

  } catch (const std::exception& e) {

    ExecEnv::log().error("Unexpected error searching with regex '{}': {}", search_spec, e.what());
    return {};

  }

}
