//
// Created by kellerberrin on 16/8/21.
//

#include "kgl_literature/kgl_literature.h"

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


// If not initialized then today's date is assumed.
kel::DateGP kgl::PublicationSummary::downloadDate() const {

  if (download_date_.notInitialized()) {

    DateGP todays_date;
    todays_date.setToday();
    return todays_date;

  }

  return download_date_;

}


std::ostream& kgl::PublicationSummary::output(std::ostream& out_stream, bool detail, char delimiter) const {

  out_stream << "PMID" << delimiter << pmid() << '\n';
  out_stream << "PUB_DATE" << delimiter << publicationDate().text() << '\n'; // YYYY-MMM-DD
  out_stream << "DOWNLOAD_DATE" << delimiter << downloadDate().text() << '\n'; // YYYY-MMM-DD
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

std::ostream& kgl::PublicationSummary::latexBiblio(std::ostream& out_stream) const {

  out_stream << "@Article{pmid" << pmid() << ",\n";

  // Authors
  size_t author_count{0};
  const size_t authors_per_line{10};
  std::string author_string{"author={"};
  for (auto const& author : authors()) {

    ++author_count;
    author_string += '{' + author.first + '}';
    if (not author.second.empty()) {

      author_string += ", ";
      author_string += author.second;

    }
    if (author != authors().back()) {

      author_string += " and ";

    }
    if (author_count % authors_per_line == 0) {

      author_string += '\n';

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


std::ostream& kgl::PublicationSummary::extendedBiblio(std::ostream& out_stream, bool detail, char delimiter) const {

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

