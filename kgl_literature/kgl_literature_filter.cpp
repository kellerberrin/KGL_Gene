//
// Created by kellerberrin on 29/8/21.
//


#include "kgl_literature.h"
#include "kgl_literature_filter.h"

namespace kgl = kellerberrin::genome;



// Filters.

bool kgl::ChemicalFilter::applyFilter(const PublicationSummary& publication) const {

  for (auto const& [chem_code, description] : publication.chemicals()) {

    for (auto const& find_chem : chemical_list_) {

      if (chem_code == find_chem) {

        return true;

      }

    }

  }

  return false;

}


bool kgl::MeSHFilter::applyFilter(const PublicationSummary& publication) const {

  for (auto const& [MeSH_code, description] : publication.chemicals()) {

    for (auto const& find_MeSH : MeSH_list_) {

      if (MeSH_code == find_MeSH) {

        return true;

      }

    }

  }

  return false;

}


bool kgl::AbstractTextFilter::applyFilter(const PublicationSummary& publication) const {

  for (auto const& text : text_list_) {

    if (publication.abstract().find(text) != std::string::npos) {

      return true;

    }

  }

  return false;

}


bool kgl::TitleTextFilter::applyFilter(const PublicationSummary& publication) const {

  for (auto const& text : text_list_) {

    if (publication.abstract().find(text) != std::string::npos) {

      return true;

    }

  }

  return false;

}

