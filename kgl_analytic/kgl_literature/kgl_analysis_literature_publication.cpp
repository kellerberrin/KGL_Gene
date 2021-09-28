//
// Created by kellerberrin on 21/9/21.
//

#include "kgl_analysis_literature_publication.h"
#include "kgl_literature_analysis.h"

#include <fstream>

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


void kgl::PublicationLiterature::writeAuthorAnalysis(const std::string& literature_directory) {

  std::string author_file{"author_analysis.csv"};
  std::string author_file_path = Utility::filePath(author_file, literature_directory);
  std::ofstream out_file(author_file_path);

  if (not out_file.good()) {

    ExecEnv::log().error("GeneLiterature::writeAuthorAnalysis; cannot open file: {} for output", author_file_path);
    return;

  }

  out_file << "Author" << ',' << "PubCount" << ',' << "Citations" << ',' << "TopCitedPubs" << '\n';

  LiteratureAnalysis literature_analysis(publication_map_);

  auto author_map = literature_analysis.AnalyseAuthors();

  for (auto const& [author, pub_cite_map] : author_map) {

    size_t citations{0};
    size_t publications{0};
    for (auto const& [cites, publication_ptr] : pub_cite_map) {

      ++publications;
      citations += cites;

    }

    out_file << author << ',' << publications << ',' << citations;

    // top publications.
    std::string pmids;
    const size_t max_pmid_count{10};
    size_t pmid_count{0};
    for (auto iter = pub_cite_map.rbegin(); pmid_count < max_pmid_count and iter != pub_cite_map.rend(); ++iter) {

      auto const& [cites, publication_ptr] = *iter;
      pmids += publication_ptr->pmid();
      pmids += '&';
      ++pmid_count;

    }

    out_file << ',' << pmids << '\n';

  }

}



void kgl::PublicationLiterature::writeYearAnalysis(const std::string& literature_directory) {

  std::string year_file{"year_analysis.csv"};
  std::string year_file_path = Utility::filePath(year_file, literature_directory);
  std::ofstream out_file(year_file_path);

  if (not out_file.good()) {

    ExecEnv::log().error("PublicationLiterature::writeYearAnalysis; cannot open file: {} for output", year_file_path);
    return;

  }

  out_file << "Year" << ',' << "PubCount" << ',' << "Citations" << ',' << "TopCitedPubs" << '\n';

  LiteratureAnalysis literature_analysis(publication_map_);

  auto year_map = literature_analysis.AnalyseYears();

  for (auto const& [year, pub_cite_map] : year_map) {

    size_t citations{0};
    size_t publications{0};
    for (auto const& [cites, publication_ptr] : pub_cite_map) {

      ++publications;
      citations += cites;

    }

    out_file << year << ',' << publications << ',' << citations;

    // top publications.
    std::string pmids;
    const size_t max_pmid_count{10};
    size_t pmid_count{0};
    for (auto iter = pub_cite_map.rbegin(); pmid_count < max_pmid_count and iter != pub_cite_map.rend(); ++iter) {

      auto const& [cites, publication_ptr] = *iter;
      pmids += publication_ptr->pmid();
      pmids += '&';
      ++pmid_count;

    }

    out_file << ',' << pmids << '\n';

  }

}




void kgl::PublicationLiterature::writeJournalAnalysis(const std::string& literature_directory) {

  std::string journal_file{"journal_analysis.csv"};
  std::string journal_file_path = Utility::filePath(journal_file, literature_directory);
  std::ofstream out_file(journal_file_path);

  if (not out_file.good()) {

    ExecEnv::log().error("PublicationLiterature::writeJournalAnalysis; cannot open file: {} for output", journal_file_path);
    return;

  }

  out_file  << "Journal" << ',' << "PubCount" << ',' << "Citations" << ',' << "TopCitedPubs" << '\n';

  LiteratureAnalysis literature_analysis(publication_map_);

  auto journal_map = literature_analysis.AnalyseJournal();

  for (auto const& [journal, pub_cite_map] : journal_map) {

    size_t citations{0};
    size_t publications{0};
    for (auto const& [cites, publication_ptr] : pub_cite_map) {

      ++publications;
      citations += cites;

    }

    out_file << journal << ',' << publications << ',' << citations;

    // top publications.
    std::string pmids;
    const size_t max_pmid_count{10};
    size_t pmid_count{0};
    for (auto iter = pub_cite_map.rbegin(); pmid_count < max_pmid_count and iter != pub_cite_map.rend(); ++iter) {

      auto const& [cites, publication_ptr] = *iter;
      pmids += publication_ptr->pmid();
      pmids += '&';
      ++pmid_count;

    }

    out_file << ',' << pmids << '\n';

  }

}


void kgl::PublicationLiterature::writeCitationPeriod(const std::string& literature_directory) {

  std::string citation_period_file{"citation_period_analysis.csv"};
  std::string citation_file_path = Utility::filePath(citation_period_file, literature_directory);
  std::ofstream out_file(citation_file_path);

  if (not out_file.good()) {

    ExecEnv::log().error("PublicationLiterature::writeCitationPeriod; cannot open file: {} for output", citation_file_path);
    return;

  }

  out_file  << "Months" << ',' << "CitationCount" << '\n';

  LiteratureAnalysis literature_analysis(publication_map_);

  auto citation_period_map = literature_analysis.AnalyseCitationPeriod();

  for (auto const& [month, cite_count] : citation_period_map) {


    out_file << month << ',' << cite_count << '\n';

  }


}


kel::DateGP kgl::PublicationLiterature::mostRecentPublication() {

  LiteratureAnalysis literature_analysis(publication_map_);

  auto date_publication_map = literature_analysis.SortPublicationDate();

  if (not date_publication_map.empty()) {

    auto recent_iter = date_publication_map.rbegin();
    auto [pub_date, publication_ptr] = *recent_iter;

    return pub_date;

  }

  return DateGP();

}