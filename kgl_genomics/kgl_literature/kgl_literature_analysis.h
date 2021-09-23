//
// Created by kellerberrin on 20/9/21.
//


#include "kgl_pubmed_resource.h"

#include <memory>

#ifndef KGL_LITERATURE_ANALYSIS_H
#define KGL_LITERATURE_ANALYSIS_H


namespace kellerberrin::genome {   //  organization level namespace

// The key is "Surname_FirstInitial" for example "Jones_B", the value map is .key=citations and .value=publication_ptr.
using LitAuthorMap = std::map<std::string, std::multimap<std::size_t, std::shared_ptr<const PublicationSummary>>>;
// The key is Year as integer, the value map is .key=citations and .value=publication_ptr.
using LitYearMap = std::map<std::string, std::multimap<std::size_t, std::shared_ptr<const PublicationSummary>>>;
// The key is Journal, the value map is .key=citations and .value=publication_ptr.
using LitJournalMap = std::map<std::string, std::multimap<std::size_t, std::shared_ptr<const PublicationSummary>>>;

class LiteratureAnalysis {

public:

  explicit LiteratureAnalysis(const LitPublicationMap& publication_map) : publication_map_(publication_map) {}
  ~LiteratureAnalysis() = default;

  LitAuthorMap AnalyseAuthors() const;
  LitYearMap AnalyseYears() const;
  LitJournalMap AnalyseJournal() const;

private:

  LitPublicationMap publication_map_;

};



} // namespace


#endif // KGL_LITERATURE_ANALYSIS_H
