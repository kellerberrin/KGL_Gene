//
// Created by kellerberrin on 16/8/21.
//

#include "kgl_pubmed.h"

namespace kgl = kellerberrin::genome;

// DD-MM-YYYY or MM-YYYY or YYYY
std::string kgl::PubMedPublicationSummary::publicationDate() const {

  std::string date_string =  publicationYear();

  if (not publicationMonth().empty()) {

    date_string += "-";
    date_string += publicationMonth();

  }

  if (not publicationDay().empty()) {

    date_string += "-";
    date_string += publicationDay();

  }

  return date_string;

}


std::ostream& kgl::PubMedPublicationSummary::output(std::ostream& out_stream, bool detail, char delimiter) const {

  out_stream << "PMID" << delimiter << pmid() << '\n';
  out_stream << "PUB_DATE" << delimiter << publicationDate() << '\n'; // DD-MM-YYYY or MM-YYYY
  out_stream << "JOURNAL" << delimiter << journal() << '\n';
  out_stream << "ISSN" << delimiter << journalISSN() << '\n';
  out_stream << "ISSUE" << delimiter << journalIssue() << '\n'; // may be empty
  out_stream << "VOLUME" << delimiter << journalVolume() << '\n';  // may be empty
  out_stream << "DOI" << delimiter << doi() << '\n';
  out_stream << "TITLE" << delimiter << title() << '\n';
  out_stream << "ABSTRACT" << delimiter << abstract() << '\n';
  out_stream << "AUTHORS" << delimiter;
  for (auto const& author : authors()) {

    out_stream << author.first << "&" << author.second;
    if (author != authors().back()) {

      out_stream << delimiter;

    }

  }
  out_stream << '\n';

  out_stream << "CHEMICALS" << delimiter;
  for (auto const& [MeSH, desc] : chemicals()) {

    out_stream << MeSH << "&" << desc;
    if (MeSH != chemicals().back().first) {

      out_stream << delimiter << '\n';

    }

  }
  out_stream << '\n';


  out_stream << "MESH" << delimiter;
  for (auto const& [MeSH, desc] : MeshCodes()) {

    out_stream << MeSH << "&" << desc;
    if (MeSH != MeshCodes().back().first) {

      out_stream << delimiter << '\n';

    }

  }
  out_stream << '\n';

  out_stream << "CITED_BY(" << citedBy().size() << ")" << delimiter;
  if (detail) {

    for (auto const& cites : citedBy()) {

      out_stream << cites;
      if (cites != *citedBy().rbegin()) {

        out_stream << delimiter;

      }

    }

  }
  out_stream << '\n';

  out_stream << "REFERENCES(" << references().size() << ")" << delimiter;
  if (detail) {

    for (auto const& ref : references()) {

      if (ref.first.empty()) {

        out_stream << "&" << ref.second << '\n';

      } else {

        out_stream << ref.first;

      }
      if (ref != *references().rbegin()) {

        out_stream << delimiter;

      }

    }

  }
  out_stream << '\n';

  return out_stream;

}

std::ostream& kgl::PubMedPublicationSummary::latexBiblio(std::ostream& out_stream) const {

  out_stream << "@Article{pmid" << pmid() << ",\n";

  // Authors
  std::string author_string{"author={"};
  for (auto const& author : authors()) {

    author_string += author.first;
    if (not author.second.empty()) {

      author_string += ", ";
      author_string += author.second;

    }
    if (author != authors().back()) {

      author_string += " and \n";

    }

  }

  author_string += "}";
  out_stream << author_string;
  out_stream << ",\ntitle={{" << title() << "}}";

  out_stream << ",\njournal={" << journal() << "}";

  out_stream << ",\nyear={" << publicationYear() << "}";

  if (not publicationMonth().empty()) {

    out_stream << ",\nmonth={" << publicationMonth() << "}";

  }

  if (not publicationDay().empty()) {

    out_stream << ",\nday={" << publicationDay() << "}";

  }

  if (not journalVolume().empty()) {

    out_stream << ",\nvolume={" << journalVolume() << "}";

  }

  if (not journalIssue().empty()) {

    out_stream << ",\nnumber={" << journalIssue() << "}";

  }


  if (not abstract().empty()) {

    out_stream << ",\nabstract={" << abstract() << "}";

  }

  if (not journalISSN().empty()) {

    out_stream << ",\nissn={" << journalISSN() << "}";

  }

  if (not doi().empty()) {

    out_stream << ",\ndoi={" << doi() << "}";

  }

  out_stream << "\n}\n";

  return out_stream;

}


std::ostream& kgl::PubMedPublicationSummary::extendedBiblio(std::ostream& out_stream, bool detail, char delimiter) const {

  latexBiblio(out_stream);

  out_stream << '\n';

  out_stream << "CHEMICALS" << delimiter;
  for (auto const& [MeSH, desc] : chemicals()) {

    out_stream << MeSH << "&" << desc;
    if (MeSH != chemicals().back().first) {

      out_stream << delimiter << '\n';

    }

  }
  out_stream << '\n';


  out_stream << "MESH" << delimiter;
  for (auto const& [MeSH, desc] : MeshCodes()) {

    out_stream << MeSH << "&" << desc;
    if (MeSH != MeshCodes().back().first) {

      out_stream << delimiter << '\n';

    }

  }
  out_stream << '\n';

  out_stream << "CITED_BY(" << citedBy().size() << ")" << delimiter;
  if (detail) {

    for (auto const& cites : citedBy()) {

      out_stream << cites;
      if (cites != *citedBy().rbegin()) {

        out_stream << delimiter;

      }

    }

  }
  out_stream << '\n';

  out_stream << "REFERENCES(" << references().size() << ")" << delimiter;
  if (detail) {

    for (auto const& ref : references()) {

      if (ref.first.empty()) {

        out_stream << "&" << ref.second << '\n';

      } else {

        out_stream << ref.first;

      }
      if (ref != *references().rbegin()) {

        out_stream << delimiter;

      }

    }

  }
  out_stream << '\n';

  return out_stream;

}


// Filters.

bool kgl::PubMedPublicationSummary::hasChemical(const std::vector<std::string>& chemical_list) const {

  for (auto const& [MeSH_code, description] : chemicals()) {

    for (auto const& find_MeSH : chemical_list) {

      if (MeSH_code == find_MeSH) {

        return true;

      }

    }

  }

  return false;

}


bool kgl::PubMedPublicationSummary::hasMeSHCode(const std::vector<std::string>& MeSH_list) const {

  for (auto const& [MeSH_code, description] : MeshCodes()) {

    for (auto const& find_MeSH : MeSH_list) {

      if (MeSH_code == find_MeSH) {

        return true;

      }

    }

  }

  return false;

}


bool kgl::PubMedPublicationSummary::hasTitleText(const std::vector<std::string>& search_text) const {

  for (auto const& text : search_text) {

    if (title().find(text) != std::string::npos) {

      return true;

    }

  }

  return false;

}


bool kgl::PubMedPublicationSummary::hasAbstractText(const std::vector<std::string>& search_text) const {

  for (auto const& text : search_text) {

    if (abstract().find(text) != std::string::npos) {

      return true;

    }

  }

  return false;

}

bool kgl::PubMedPublicationSummary::hasAuthor(const std::vector<std::pair<std::string, std::string>>& author_list) const {

  for (auto const& [surname, initials] : author_list) {

    for (auto const& [pub_surname, pub_initials] : authors()) {

      if (pub_surname == surname and (initials.empty() or initials == pub_initials)) {

        return true;

      }

    }

  }

  return false;

}