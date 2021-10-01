//
// Created by kellerberrin on 16/8/21.
//

#ifndef KGL_LITERATURE_H
#define KGL_LITERATURE_H

#include "kel_date_time.h"

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



class PublicationSummary;
// key = pmid_, value = publication details record.
using LitPublicationMap = std::map<std::string, std::shared_ptr<const PublicationSummary>>;

using APIPublicationMap = std::map<std::string, std::shared_ptr<PublicationSummary>>;

// key = pmid_, value is a pair of .first = download-date and .second = a unique unique set of pmids that cite, or are referenced by, the key pmid_.
using LitCitationMap = std::map<std::string, std::pair<DateGP, std::set<std::string>>>;


class PublicationSummary {

public:

  explicit PublicationSummary(std::string id) : pmid_(std::move(id)) {}
  ~PublicationSummary() = default;

  // Getters
  [[nodiscard]] const std::string& pmid() const { return pmid_; }
  [[nodiscard]] DateGP publicationDate() const { return publication_date_; }
  [[nodiscard]] DateGP downloadDate() const; // If not initialized then today's date is assumed.
  [[nodiscard]] std::string publicationYear() const { return publication_date_.year_text(); }
  [[nodiscard]] std::string publicationMonth() const { return publication_date_.month_text(); }
  [[nodiscard]] std::string publicationDay() const { return publication_date_.day_text(); }
  [[nodiscard]] const std::string& journal() const { return journal_; }
  [[nodiscard]] const std::string& journalISSN() const { return journal_issn_;  } // may be empty
  [[nodiscard]] const std::string& journalIssue() const { return journal_issue_;  } // may be empty
  [[nodiscard]] const std::string& journalVolume() const { return journal_volume_; }  // may be empty
  [[nodiscard]] const std::string& doi() const { return doi_; }
  [[nodiscard]] const std::string& title() const { return title_; }
  [[nodiscard]] const std::string& abstract() const { return abstract_; }
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& authors() const { return authors_; }  // surname, initials
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& chemicals() const { return chemicals_; }  // .first is MeshCode .second is description
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& MeshCodes() const { return MeSHcodes_; }  // .first is MeshCode .second is description
  [[nodiscard]] const std::set<std::string>& citedBy() const { return cited_by_articles_; } // pmids only.
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& references() const { return references_;  } // pmids only.

  // Setters
  void publicationDate(const DateGP& publication_date) { publication_date_ = publication_date; }
  void downloadDate(const DateGP& download_date) { download_date_ = download_date; }
  void journal(const std::string& name) { journal_ = name; }
  void journalISSN(const std::string& issn) { journal_issn_ = issn;  } // may be empty
  void journalIssue(const std::string& issue) { journal_issue_ = issue;  } // may be empty
  void journalVolume(const std::string& volume) { journal_volume_ = volume; }  // may be empty
  void doi(const std::string& id) { doi_ = id; }
  void title(const std::string& text) { title_ = text; }
  void abstract(const std::string& text) { abstract_ = text; }
  void addAuthor(const std::string& author, const std::string& initials) { authors_.emplace_back(author, initials); }  // surname, initials
  void addChemical(const std::string& MeSHCode, const std::string& description) { chemicals_.emplace_back(MeSHCode, description); }  // .first is MeshCode .second is description
  void addMeSHCode(const std::string& MeSHCode, const std::string& description) { MeSHcodes_.emplace_back(MeSHCode, description); }  // .first is MeshCode .second is description
  void citations(const std::set<std::string>& cite_pmid_set) { cited_by_articles_ = cite_pmid_set; } // pmids only.
  void reference(const std::vector<std::pair<std::string, std::string>>& ref_pmid_set) { references_ = ref_pmid_set; } // pmids only.
  [[nodiscard]] bool addCitation(const std::string& cite_pmid) { auto [iter, result ] = cited_by_articles_.insert(cite_pmid); return result; } // pmids only.
  void addReference(const std::string& ref_pmid, const std::string& citation_text) { references_.emplace_back(ref_pmid, citation_text);  } // pmids only.

  // Output the record to an ostream, one field per line.
  std::ostream& output(std::ostream& out_stream, bool detail = false, char delimiter = SV_DELIMITER_) const;
  // Output in latex biblographic (.bib) format.
  std::ostream& latexBiblio(std::ostream& out_stream) const;
  // Adds reference, citation and Mash data to the biblio.
  std::ostream& extendedBiblio(std::ostream& out_stream, bool detail = false, char delimiter = SV_DELIMITER_) const;

private:

  const std::string pmid_;
  DateGP publication_date_;
  DateGP download_date_;   // The date that this record was downloaded from Pubmed, if not initialized then today's date is assumed.
  std::string journal_;
  std::string journal_issn_;
  std::string journal_issue_;  // may be empty
  std::string journal_volume_;  // may be empty
  std::string doi_;
  std::string title_;
  std::string abstract_;
  std::vector<std::pair<std::string, std::string>> authors_;  // .first = surname, .second = initials (can be empty).
  std::vector<std::pair<std::string, std::string>> chemicals_;   // .first is MeshCode .second is description
  std::vector<std::pair<std::string, std::string>> MeSHcodes_;   // .first is MeshCode .second is description
  std::set<std::string> cited_by_articles_;  // pmids only.
  std::vector<std::pair<std::string, std::string>> references_;   // .first=pmid (can be empty), .second is citation_text

  constexpr const static char SV_DELIMITER_{'|'};

};


} // namespace


#endif // KGL_LITERATURE_H
