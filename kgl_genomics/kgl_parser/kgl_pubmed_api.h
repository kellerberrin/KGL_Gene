//
// Created by kellerberrin on 11/8/21.
//

#ifndef KGL_PUBMED_API_H
#define KGL_PUBMED_API_H

#include "kgl_pubmed.h"
#include "kel_rest_api.h"
#include "kgl_resource_db.h"

#include <string>
#include <memory>
#include <map>
#include <set>
#include <vector>
#include <chrono>


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////?////////////////////////////////////////////////////
//
// Pubmed API requestor. This is a resource that can be requested by analysis packages at runtime.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// key = pmid_, value = publication details record.
using LitPublicationMap = std::map<std::string, PubMedPublicationSummary>;

// key = pmid_, value = unique set of pmids that cite, or are referenced by, the key pmid_.
using LitCitationMap = std::map<std::string, std::set<std::string>>;


class PubmedRequester : public ResourceBase {

public:

  explicit PubmedRequester(const std::string& identifier) : ResourceBase(identifier) {}
  ~PubmedRequester() override = default;

  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::PUBMED_API;  }

  // Only the data for unique pmids is returned, for convenience the requesting array of pmids can contain non-unique pmids.
  // Large requests are automatically throttled to Pubmed requirements (max 10 api calls a second).
  // Citation and reference information.
  [[nodiscard]] LitCitationMap getCitations(const std::vector<std::string>& pmid_vector) const;
  [[nodiscard]] LitCitationMap getReferences(const std::vector<std::string>& pmid_vector) const;
  // Get all publication summary data including citations and references.
  [[nodiscard]] LitPublicationMap getPublicationDetails(const std::vector<std::string>& pmid_vector) const;

private:

  mutable RestAPI pubmed_rest_api_;

  // Kellerberrin Pubmed API key.
  const std::string API_KEY {"&api_key=8cd3dde4cbf1eeb71b5ae469ae8a99247609"};
  // Pubmed resource constraints.
  constexpr const static size_t REQUESTS_PER_SECOND_{10};
  constexpr const static size_t PMID_PER_REQUEST_{10};
  constexpr const static std::chrono::milliseconds BATCH_TIME_INTERVAL_{1000};

  // The citation/reference constants.
  const std::string PUBMED_ELINK_URL_{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"};
  const std::string PUBMED_ARTICLE_CITEDBY_ARGS_{"dbfrom=pubmed&linkname=pubmed_pubmed_citedin"};
  const std::string PUBMED_ARTICLE_REFERENCES_ARGS_{"dbfrom=pubmed&linkname=pubmed_pubmed_refs"};

  [[nodiscard]] LitCitationMap getCitationReference(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) const;
  [[nodiscard]] LitCitationMap citationBatch(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) const;

  // The publication detail constants.
  const std::string PUBMED_EFETCH_URL_{"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"};
  const std::string PUBMED_ARTICLE_DETAIL_ARGS_{"db=pubmed&retmode=xml"};

  [[nodiscard]] LitPublicationMap publicationBatch(const std::vector<std::string>& pmid_vector) const;


};




} // namespace


#endif //KGL_PUBMED_API_H
