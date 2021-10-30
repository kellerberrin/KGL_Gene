//
// Created by kellerberrin on 21/9/21.
//

#ifndef KGL_ANALYSIS_LITERATURE_PUBLICATION_H
#define KGL_ANALYSIS_LITERATURE_PUBLICATION_H


#include "kgl_pubmed_resource.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PublicationLiterature {

public:

  explicit PublicationLiterature(const LitPublicationMap& publication_map) : publication_map_(publication_map) {}
  ~PublicationLiterature() = default;

  void writeAuthorAnalysis(const std::string& literature_directory);
  void writeYearAnalysis(const std::string& literature_directory);
  void writeJournalAnalysis(const std::string& literature_directory);
  void writeCitationPeriod(const std::string& literature_directory);
  void writeCitationVariance(const std::string& literature_directory);
  std::shared_ptr<const PublicationSummary> mostRecentPublication();
  void writeCitationQuantiles(const std::string& literature_directory);
  void writeCitationHistogram(const std::string& literature_directory);
  void writeCitationData(const std::string& literature_directory);
  void writePublicationCitations(const std::string& literature_directory, const std::string& publication_pmid);

private:

  const LitPublicationMap publication_map_;

};



} // namespace


#endif //KGL_ANALYSIS_LITERATURE_PUBLICATION_H
