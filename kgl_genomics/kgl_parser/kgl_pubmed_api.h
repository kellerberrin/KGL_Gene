//
// Created by kellerberrin on 11/8/21.
//

#ifndef KGL_PUBMED_API_H
#define KGL_PUBMED_API_H

#include "kel_rest_api.h"

#include <string>
#include <memory>
#include <map>
#include <set>
#include <vector>
#include <chrono>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Pubmed Specific Objects. A vector of article Pubmed pmids (article idents) can be used to query article references, citations, and details.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PubMedPublicationSummary {

public:

  PubMedPublicationSummary() = default;
  ~PubMedPublicationSummary() = default;

  // Getters
  [[nodiscard]] const std::string& pmid() const { return pmid_; }
  [[nodiscard]] const std::string& publicationDate() const { return publication_date_; } // DD-MM-YYYY or MM-YYYY
  [[nodiscard]] const std::string& journal() const { return journal_; }
  [[nodiscard]] const std::string& journalIssue() const { return journal_issue_;  } // may be empty
  [[nodiscard]] const std::string& journalVolume() const { return journal_volume_; }  // may be empty
  [[nodiscard]] const std::string& doi() const { return doi_; }
  [[nodiscard]] const std::string& title() const { return title_; }
  [[nodiscard]] const std::string& abstract() const { return abstract_; }
  [[nodiscard]] const std::vector<std::string>& authors() const { return authors_; }  // surname, initials
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& chemicals() const { return chemicals_; }  // .first is MeshCode .second is description
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& MeshCodes() const { return MeSHcodes_; }  // .first is MeshCode .second is description
  [[nodiscard]] const std::set<std::string>& citedBy() const { return cited_by_articles_; } // pmids only.
  [[nodiscard]] const std::set<std::string>& reference() const { return references_;  } // pmids only.

  // Setters
  void pmid(const std::string& id) { pmid_ = id; }
  void publicationDate(const std::string& pub_date) { publication_date_ = pub_date; } // DD-MM-YYYY or MM-YYYY
  void journal(const std::string& name) { journal_ = name; }
  void journalIssue(const std::string& issue) { journal_issue_ = issue;  } // may be empty
  void journalVolume(const std::string& volume) { journal_volume_ = volume; }  // may be empty
  void doi(const std::string& id) { doi_ = id; }
  void title(const std::string& text) { title_ = text; }
  void abstract(const std::string& text) { abstract_ = text; }
  void authors(const std::vector<std::string>& author_list) { authors_ = author_list; }  // surname, initials
  void  chemicals(const std::vector<std::pair<std::string, std::string>>& chem_list) { chemicals_ = chem_list; }  // .first is MeshCode .second is description
  void MeshCodes(const std::vector<std::pair<std::string, std::string>>& mesh_list) { MeSHcodes_ = mesh_list; }  // .first is MeshCode .second is description
  void citedBy(const std::set<std::string>& cite_set) { cited_by_articles_ = cite_set; } // pmids only.
  void reference(const std::set<std::string>& ref_set) { references_ = ref_set;  } // pmids only.

  std::string pmid_;
  std::string publication_date_;  // DD-MM-YYYY or MM-YYYY
  std::string journal_;
  std::string journal_issue_;  // may be empty
  std::string journal_volume_;  // may be empty
  std::string doi_;
  std::string title_;
  std::string abstract_;
  std::vector<std::string> authors_;  // surname, initials
  std::vector<std::pair<std::string, std::string>> chemicals_;   // .first is MeshCode .second is description
  std::vector<std::pair<std::string, std::string>> MeSHcodes_;   // .first is MeshCode .second is description
  std::set<std::string> cited_by_articles_;  // pmids only.
  std::set<std::string> references_;   // pmids only.

};
// key = pmid_, value = publication details.
using LitPublicationMap = std::map<std::string, PubMedPublicationSummary>;

// key = pmid_, value = vector of pmids that cite, or are referenced by, the key pmid_.
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
