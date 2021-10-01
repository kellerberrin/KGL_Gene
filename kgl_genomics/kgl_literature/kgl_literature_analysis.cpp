//
// Created by kellerberrin on 20/9/21.
//

#include "kel_utility.h"
#include "kgl_literature_filter.h"
#include "kgl_literature_analysis.h"

#include "kel_utility.h"

#include <numeric>

namespace kgl = kellerberrin::genome;


kgl::LiteratureAnalysis::LiteratureAnalysis(const LitPublicationMap& publication_map) {

  static const PlasmodiumFilter pf_filter;

  for (auto const& [pmid, publication_ptr] : publication_map) {

    if (not pf_filter.applyFilter(*publication_ptr)) {

      publication_map_.emplace(pmid, publication_ptr);

    }

  }

}


kgl::LitAuthorMap kgl::LiteratureAnalysis::analyseAuthors() const {

  LitAuthorMap author_map;

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    for (auto const& [surname, initials] : publication_ptr->authors()) {

      // Construct Author key

      std::string author_key = Utility::findAndReplaceAll(surname, ",", "_");

      if (not initials.empty()) {

        author_key += '_';
        std::string upper_initials = Utility::toupper(initials);
        author_key += upper_initials[0];

      }

      // Check the Author map and update.
      auto result = author_map.find(author_key);
      if (result == author_map.end()) {

        std::multimap<size_t, std::shared_ptr<const PublicationSummary>> pub_map{{publication_ptr->citedBy().size(), publication_ptr}};
        author_map.emplace(author_key, pub_map);

      } else {

        auto& [key, pub_cite_map] = *result;
        pub_cite_map.emplace(publication_ptr->citedBy().size(), publication_ptr);

      }

    } // For author

  } // For publication.

  return author_map;

}


kgl::LitYearMap kgl::LiteratureAnalysis::analyseYears() const {

  LitYearMap year_map;

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    // Check the Author map and update.
    auto result = year_map.find(publication_ptr->publicationYear());
    if (result == year_map.end()) {

      std::multimap<size_t, std::shared_ptr<const PublicationSummary>> pub_map{{publication_ptr->citedBy().size(), publication_ptr}};
      year_map.emplace(publication_ptr->publicationYear(), pub_map);

    } else {

      auto& [key, pub_cite_map] = *result;
      pub_cite_map.emplace(publication_ptr->citedBy().size(), publication_ptr);

    }

  } // For publication.

  return year_map;

}



kgl::LitJournalMap kgl::LiteratureAnalysis::analyseJournal() const {

  LitJournalMap journal_map;

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    std::string journal_key = Utility::findAndReplaceAll(publication_ptr->journal(), ",", "_");

    // Check the Author map and update.
    auto result = journal_map.find(journal_key);
    if (result == journal_map.end()) {

      std::multimap<size_t, std::shared_ptr<const PublicationSummary>> pub_map{{publication_ptr->citedBy().size(), publication_ptr}};
      journal_map.emplace(journal_key, pub_map);

    } else {

      auto& [key, pub_cite_map] = *result;
      pub_cite_map.emplace(publication_ptr->citedBy().size(), publication_ptr);

    }

  } // For publication.

  return journal_map;

}


kgl::LitCitationPeriodMap kgl::LiteratureAnalysis::analyseCitationPeriod() const {

  LitCitationPeriodMap citation_period_map;
  const size_t minimum_citation_count = 0;
  const size_t maximum_citation_count = 1000000;

  size_t found_cite{0};
  size_t cite_not_found{0};
  for (auto const& [pmid, publication_ptr] : publication_map_) {

    if (publication_ptr->citedBy().size() < minimum_citation_count or publication_ptr->citedBy().size() > maximum_citation_count) {

      continue;

    }

    for (auto const& cite_pmid : publication_ptr->citedBy()) {

      auto result = publication_map_.find(cite_pmid);
      if (result == publication_map_.end()) {

        ++cite_not_found;

      } else {

        ++found_cite;
        auto const& [_pmid, cite_publication_ptr] = *result;
        size_t cite_months = DateGP::monthsDifference(publication_ptr->publicationDate(), cite_publication_ptr->publicationDate());

        auto period_result = citation_period_map.find(cite_months);
        if (period_result == citation_period_map.end()) {

          citation_period_map.emplace(cite_months, 1);

        } else {

          auto& [_months, cite_count] = *period_result;
          ++cite_count;

        }

      }

    }

  }

  ExecEnv::log().info("LiteratureAnalysis::analyseCitationPeriod; Found citations: {}, Citations Not Found: {}", found_cite, cite_not_found);

  return citation_period_map;

}


kgl::LitDateMap kgl::LiteratureAnalysis::sortPublicationDate() const {

  LitDateMap date_sorted_publications;

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    date_sorted_publications.emplace(publication_ptr->publicationDate(), publication_ptr);

  }

  return date_sorted_publications;

}




kgl::LitCitationVarianceMap kgl::LiteratureAnalysis::analyseCitationPercent() const {

  const size_t minimum_citation_count = 0;
  const size_t maximum_citation_count = 1000000;
  const size_t max_period_months = 120;

  size_t found_cite{0};
  size_t cite_not_found{0};
  std::vector<std::vector<double>> cite_vector{max_period_months};
  for (auto const& [pmid, publication_ptr] : publication_map_) {

    if (publication_ptr->citedBy().size() < minimum_citation_count or publication_ptr->citedBy().size() > maximum_citation_count) {

      continue;

    }

    LitCitationPeriodMap citation_period_map;
    for (auto const& cite_pmid : publication_ptr->citedBy()) {

      auto result = publication_map_.find(cite_pmid);
      if (result == publication_map_.end()) {

        ++cite_not_found;

      } else {

        ++found_cite;
        auto const& [_pmid, cite_publication_ptr] = *result;
        size_t cite_months = DateGP::monthsDifference(publication_ptr->publicationDate(), cite_publication_ptr->publicationDate());

        if (cite_months < max_period_months) {

          auto period_result = citation_period_map.find(cite_months);
          if (period_result == citation_period_map.end()) {

            citation_period_map.emplace(cite_months, 1);

          } else {

            auto& [_months, cite_count] = *period_result;
            ++cite_count;

          }

        }

      }

    }

    if (not citation_period_map.empty()) {

      size_t total_cite_count{0};
      for (auto const& [period, cites] : citation_period_map) {

        total_cite_count += cites;

      }

      size_t cumulative_cite_count{0};
      for (size_t period = 0; period < max_period_months; ++period) {

        auto result = citation_period_map.find(period);
        if (result != citation_period_map.end()) {

          auto const& [_period, cite_count] = *result;
          cumulative_cite_count += cite_count;

        }

        double cite_proportion = static_cast<double>(cumulative_cite_count) / static_cast<double>(total_cite_count);

        cite_vector[period].push_back(cite_proportion);

      }

    }

  }

  LitCitationVarianceMap citation_variance_map;
  for (size_t period = 0; period < max_period_months; ++period) {

    auto [mean, stdev] = Utility::stddev(cite_vector[period]);
    citation_variance_map.emplace(period, std::pair<double, double>{mean, stdev});

  }

  ExecEnv::log().info("LiteratureAnalysis::analyseCitationPercent; Found citations: {}, Citations Not Found: {}", found_cite, cite_not_found);

  return citation_variance_map;

}

kgl::CitationQuantile kgl::LiteratureAnalysis::analyseCitationQuartiles() const {

  CitationQuantile citation_quantile;
  const size_t MONTHS_ELAPSED{120};

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    auto months = DateGP::monthsDifference(publication_ptr->downloadDate(), publication_ptr->publicationDate());
    if (months >= MONTHS_ELAPSED) {

      citation_quantile.addElement(publication_ptr->citedBy().size(), publication_ptr);

    }

  }

  return citation_quantile;

}
