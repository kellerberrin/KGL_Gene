//
// Created by kellerberrin on 16/8/21.
//

#include "kgl_pubmed.h"

namespace kgl = kellerberrin::genome;



std::ostream& kgl::PubMedPublicationSummary::output(std::ostream& out_stream, char delimiter) const {

  out_stream << "PMID" << delimiter << pmid() << '\n';
  out_stream << "PUB_DATE" << delimiter << publicationDate() << '\n'; // DD-MM-YYYY or MM-YYYY
  out_stream << "JOURNAL" << delimiter << journal() << '\n';
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
  for (auto const& cites : citedBy()) {

    out_stream << cites;
    if (cites != *citedBy().rbegin()) {

      out_stream << delimiter;

    }

  }
  out_stream << '\n';

  out_stream << "REFERENCES(" << references().size() << ")" << delimiter;
  for (auto const& ref : references()) {

    if (ref.first.empty()) {

      out_stream << "&" << ref.second << '\n';

    } else {

      out_stream << ref.first;

    }
    if (ref != *references().rbegin()) {

      out_stream << delimiter;

    } else {

      out_stream << '\n';

    }

  }

  return out_stream;

}
