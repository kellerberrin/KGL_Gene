//
// Created by kellerberrin on 11/8/21.
//

#ifndef KGL_PUBMED_API_H
#define KGL_PUBMED_API_H

#include "kel_rest_api.h"

#include <string>
#include <memory>
#include <map>
#include <vector>
#include <chrono>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Pubmed Specific Objects. A vector of article Pubmed pmids (article idents) can be used to query article references, citations, and details.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct PubMedPublicationDetails {

  std::string pmid;
  std::string publication_date;
  std::string journal;
  std::string journal_issue;
  std::string journal_volume;
  std::string doi;
  std::string title;
  std::string abstract;
  std::vector<std::string> authors;
  std::vector<std::string> MeSHcodes;
  std::vector<std::string> cited_by_articles;
  std::vector<std::string> references;

};
// key = pmid, value = publication details.
using LitPublicationMap = std::map<std::string, PubMedPublicationDetails>;

// key = pmid, value = vector of pmids that cite, or are referenced by, the key pmid.
using LitCitationMap = std::map<std::string, std::vector<std::string>>;

class PubmedRequester {

public:

  PubmedRequester() = default;
  ~PubmedRequester() = default;

  [[nodiscard]] LitCitationMap getCitations(const std::vector<std::string>& pmid_vector);
  [[nodiscard]] LitCitationMap getReferences(const std::vector<std::string>& pmid_vector);
  [[nodiscard]] LitPublicationMap getPublicationDetails(const std::vector<std::string>& pmid_vector);

private:

  RestAPI pubmed_rest_api_;

  // Pubmed API key.
  const std::string API_KEY {"&api_key=8cd3dde4cbf1eeb71b5ae469ae8a99247609"};
  // Pubmed resource constraints.
  constexpr const static size_t REQUESTS_PER_SECOND_{10};
  constexpr const static size_t PMID_PER_REQUEST_{10};
  constexpr const static std::chrono::milliseconds BATCH_TIME_INTERVAL_{1000};

  // The publication citation/reference constants.
  const std::string PUBMED_ELINK_URL_{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"};
  const std::string PUBMED_ARTICLE_CITEDBY_ARGS_{"dbfrom=pubmed&linkname=pubmed_pubmed_citedin"};
  const std::string PUBMED_ARTICLE_REFERENCES_ARGS_{"dbfrom=pubmed&linkname=pubmed_pubmed_refs"};

  // Returned citation/reference XML nodes
  const std::string CITATION_ROOT_NODE_{"eLinkResult"};
  const std::string CITATION_RECORD_{"LinkSet"};
  const std::string CITATION_PMID_{"IdList"};
  const std::string CITATION_LINK_DB_{"LinkSetDb"};
  const std::string CITATION_LINK_SET_{"Link"};

  [[nodiscard]] LitCitationMap getCitationArgs(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args);
  [[nodiscard]] LitCitationMap citationBatch(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args);
  [[nodiscard]] LitCitationMap parseCitationXML(const std::string& citation_text);

  // The publication detail constants.
  const std::string PUBMED_EFETCH_URL_{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"};
//  const std::string PUBMED_ARTICLE_DETAIL_ARGS_{"db=pubmed&rettype=xml&retmode=xml"};
  const std::string PUBMED_ARTICLE_DETAIL_ARGS_{"db=pubmed&retmode=xml"};

  [[nodiscard]] LitPublicationMap publicationBatch(const std::vector<std::string>& pmid_vector);
  [[nodiscard]] LitPublicationMap parsePublicationXML(const std::string& publication_xml_text);

  // Returned publication XML nodes
  const std::string PUBLICATION_ROOT_NODE_{"PubmedArticleSet"};
  const std::string PUBLICATION_NODE_{"PubmedArticle"};


};




} // namespace


#endif //KGL_PUBMED_API_H
