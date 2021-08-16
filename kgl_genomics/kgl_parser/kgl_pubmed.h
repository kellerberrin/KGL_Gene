//
// Created by kellerberrin on 16/8/21.
//

#ifndef KGL_PUBMED_H
#define KGL_PUBMED_H


#include <string>
#include <memory>
#include <map>
#include <set>
#include <vector>
#include <chrono>
#include <ostream>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Pubmed Specific Objects. A vector of article Pubmed pmids (article idents) can be used to query article references, citations, and details.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PubMedPublicationSummary {

public:

  explicit PubMedPublicationSummary(std::string id) : pmid_(std::move(id)) {}
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
  [[nodiscard]] const std::set<std::string>& references() const { return references_;  } // pmids only.

  // Setters
  void publicationDate(const std::string& pub_date) { publication_date_ = pub_date; } // DD-MM-YYYY or MM-YYYY
  void journal(const std::string& name) { journal_ = name; }
  void journalIssue(const std::string& issue) { journal_issue_ = issue;  } // may be empty
  void journalVolume(const std::string& volume) { journal_volume_ = volume; }  // may be empty
  void doi(const std::string& id) { doi_ = id; }
  void title(const std::string& text) { title_ = text; }
  void abstract(const std::string& text) { abstract_ = text; }
  void addAuthor(const std::string& author) { authors_.push_back(author); }  // surname, initials
  void addChemical(const std::string& MeSHCode, const std::string& description) { chemicals_.emplace_back(MeSHCode, description); }  // .first is MeshCode .second is description
  void addMeSHCode(const std::string& MeSHCode, const std::string& description) { MeSHcodes_.emplace_back(MeSHCode, description); }  // .first is MeshCode .second is description
  void citations(const std::set<std::string>& cite_pmid_set) { cited_by_articles_ = cite_pmid_set; } // pmids only.
  void reference(const std::set<std::string>& ref_pmid_set) { references_ = ref_pmid_set; } // pmids only.
  [[nodiscard]] bool addCitation(const std::string& cite_pmid) { auto [iter, result ] = cited_by_articles_.insert(cite_pmid); return result; } // pmids only.
  [[nodiscard]] bool addReference(const std::string& ref_pmid) { auto [iter, result ] = references_.insert(ref_pmid); return result;  } // pmids only.

  // delimites stream output
  std::ostream& output(std::ostream& out_stream, char delimiter = SV_DELIMITER_) const;

private:

  const std::string pmid_;
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

  constexpr const static char SV_DELIMITER_{','};

};


} // namespace


#endif // KGL_PUBMED_H
