//
// Created by kellerberrin on 20/9/21.
//

#include "kel_utility.h"
#include "kgl_literature_filter.h"
#include "kgl_literature_analysis.h"

namespace kgl = kellerberrin::genome;


kgl::LitAuthorMap kgl::LiteratureAnalysis::AnalyseAuthors() const {

  LitAuthorMap author_map;
  static const PlasmodiumFilter pf_filter;

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    if (not pf_filter.applyFilter(*publication_ptr)) {

      continue;

    }

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


kgl::LitYearMap kgl::LiteratureAnalysis::AnalyseYears() const {

  LitYearMap year_map;
  static const PlasmodiumFilter pf_filter;

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    if (not pf_filter.applyFilter(*publication_ptr)) {

      continue;

    }

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



kgl::LitJournalMap kgl::LiteratureAnalysis::AnalyseJournal() const {

  LitJournalMap journal_map;
  static const PlasmodiumFilter pf_filter;

  for (auto const& [pmid, publication_ptr] : publication_map_) {

    if (not pf_filter.applyFilter(*publication_ptr)) {

      continue;

    }

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