//
// Created by kellerberrin on 29/8/21.
//


#include "kgl_literature.h"

namespace kgl = kellerberrin::genome;



// Filters.

bool kgl::PublicationSummary::hasChemical(const std::vector<std::string>& chemical_list) const {

  for (auto const& [MeSH_code, description] : chemicals()) {

    for (auto const& find_MeSH : chemical_list) {

      if (MeSH_code == find_MeSH) {

        return true;

      }

    }

  }

  return false;

}


bool kgl::PublicationSummary::hasMeSHCode(const std::vector<std::string>& MeSH_list) const {

  for (auto const& [MeSH_code, description] : MeshCodes()) {

    for (auto const& find_MeSH : MeSH_list) {

      if (MeSH_code == find_MeSH) {

        return true;

      }

    }

  }

  return false;

}


bool kgl::PublicationSummary::hasTitleText(const std::vector<std::string>& search_text) const {

  for (auto const& text : search_text) {

    if (title().find(text) != std::string::npos) {

      return true;

    }

  }

  return false;

}


bool kgl::PublicationSummary::hasAbstractText(const std::vector<std::string>& search_text) const {

  for (auto const& text : search_text) {

    if (abstract().find(text) != std::string::npos) {

      return true;

    }

  }

  return false;

}

bool kgl::PublicationSummary::hasAuthor(const std::vector<std::pair<std::string, std::string>>& author_list) const {

  for (auto const& [surname, initials] : author_list) {

    for (auto const& [pub_surname, pub_initials] : authors()) {

      if (pub_surname == surname and (initials.empty() or initials == pub_initials)) {

        return true;

      }

    }

  }

  return false;

}



bool kgl::PublicationSummary::PfalciparumFilter(const PublicationSummary& publication) {

  static const std::vector<std::string> Mesh_codes{ "D010963" /* Plasmodium falciparum */
                                                    ,"D016778" /*  Malaria, Falciparum */
                                                    , "D008288" /* Malaria */ };

  static const std::vector<std::string> search_text { "Plasmodium", "plasmodium", "Falciparum",  "falciparum", "Malaria", "malaria" };

  if (publication.hasMeSHCode(Mesh_codes)) {

    return true;

  }

  if (publication.hasTitleText(search_text)) {

    return true;

  }

  if (publication.hasAbstractText(search_text)) {

    return true;

  }

  return false;

}
