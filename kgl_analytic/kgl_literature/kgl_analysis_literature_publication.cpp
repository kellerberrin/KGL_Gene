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

  auto author_map = literature_analysis.analyseAuthors();

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

  auto year_map = literature_analysis.analyseYears();

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

  auto journal_map = literature_analysis.analyseJournal();

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

  auto citation_period_map = literature_analysis.analyseCitationPeriod();

  for (auto const& [month, cite_count] : citation_period_map) {


    out_file << month << ',' << cite_count << '\n';

  }


}



void kgl::PublicationLiterature::writeCitationVariance(const std::string& literature_directory) {

  std::string citation_variance_file{"citation_variance_analysis.csv"};
  std::string citation_file_path = Utility::filePath(citation_variance_file, literature_directory);
  std::ofstream out_file(citation_file_path);

  if (not out_file.good()) {

    ExecEnv::log().error("PublicationLiterature::writeCitationVeriance; cannot open file: {} for output", citation_file_path);
    return;

  }

  out_file  << "Months" << ',' << "CitationMean" << ',' << "CitationStdev" << '\n';

  LiteratureAnalysis literature_analysis(publication_map_);

  auto citation_variance_map = literature_analysis.analyseCitationPercent();

  for (auto const& [month, cite_pair] : citation_variance_map) {

    auto const& [mean, stdev] = cite_pair;

    out_file << month << ',' << mean << ',' << stdev << '\n';

  }

}


std::shared_ptr<const kgl::PublicationSummary> kgl::PublicationLiterature::mostRecentPublication() {

  LiteratureAnalysis literature_analysis(publication_map_);

  auto date_publication_map = literature_analysis.sortPublicationDate();


  auto pub_iter = date_publication_map.rbegin();
  while(pub_iter != date_publication_map.rend()) {

    auto const& [pub_date, publication_ptr] = *pub_iter;

    if (pub_date > publication_ptr->downloadDate()) {

      ExecEnv::log().info("Publication: {}, publication date: {}, greater than download date: {}",
                          publication_ptr->pmid(), pub_date.text(), publication_ptr->downloadDate().text());

    } else {

      break;

    }

    ++pub_iter;

  }


  if (not date_publication_map.empty()) {

    auto recent_iter = date_publication_map.rbegin();

    auto [pub_date, publication_ptr] = *recent_iter;

    return publication_ptr;

  }

  // Should not happen.
  ExecEnv::log().critical("PublicationLiterature::mostRecentPublication; no publications in the publication map");

  return std::shared_ptr<const kgl::PublicationSummary>();

}

void kgl::PublicationLiterature::writeCitationQuantiles(const std::string& literature_directory) {


  std::string citation_quantile_file{"citation_quantile_analysis.csv"};
  std::string citation_file_path = Utility::filePath(citation_quantile_file, literature_directory);
  std::ofstream out_file(citation_file_path);

  if (not out_file.good()) {

    ExecEnv::log().error("PublicationLiterature::writeCitationQuantiles; cannot open file: {} for output", citation_file_path);
    return;

  }

  out_file  << "Quantile" << ',' << "Citations" << '\n';

  LiteratureAnalysis literature_analysis(publication_map_);

  auto citation_quantiles = literature_analysis.analyseCitationQuartiles();

  for (double quantile = 0.01; quantile <= 1.0; quantile += 0.01) {

    out_file << quantile << ',';

    auto quant_opt = citation_quantiles.percentile(quantile);
    if (quant_opt) {

      auto const& [cite_count, publication_ptr] = quant_opt.value();
      out_file << cite_count << '\n';

    } else {

      out_file << "N/A\n";

    }

  }

}

