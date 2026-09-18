//
// Created by kellerberrin on 22/9/21.
//

#ifndef KGL_LITERATURE_FILTER_H
#define KGL_LITERATURE_FILTER_H


#include "kgl_literature.h"
#include "kgl_literature_filter_virtual.h"
#include "kel_utility.h"


namespace kellerberrin::genome {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter returns true if any of the text list items are found in the publication MeSH list.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MeSHFilter : public LiteratureFilter {

public:

  explicit MeSHFilter(std::vector<std::string> MeSH_list) : MeSH_list_(std::move(MeSH_list)) {}
  MeSHFilter(const MeSHFilter&) = default;
  ~MeSHFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override;

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<MeSHFilter>(*this); }

private:

  std::vector<std::string> MeSH_list_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter returns true if any of the text list items are found in the publication chemical list.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ChemicalFilter : public LiteratureFilter {

public:

  explicit ChemicalFilter(std::vector<std::string> chemical_list) : chemical_list_(std::move(chemical_list)) {}
  ChemicalFilter(const ChemicalFilter&) = default;
  ~ChemicalFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override;

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<ChemicalFilter>(*this); }

private:

  std::vector<std::string> chemical_list_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter returns true if any of the text list items are found in the publication abstract.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AbstractTextFilter : public LiteratureFilter {

public:

  explicit AbstractTextFilter(std::vector<std::string> text_list) : text_list_(std::move(text_list)) {}
  AbstractTextFilter(const AbstractTextFilter&) = default;
  ~AbstractTextFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override;

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<AbstractTextFilter>(*this); }

private:

  std::vector<std::string> text_list_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter returns true if any of the text list items are found in the publication title.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TitleTextFilter : public LiteratureFilter {

public:

  explicit TitleTextFilter(std::vector<std::string> text_list) : text_list_(std::move(text_list)) {}
  TitleTextFilter(const TitleTextFilter&) = default;
  ~TitleTextFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override;

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<TitleTextFilter>(*this); }

private:

  std::vector<std::string> text_list_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Logical filters, used to compose a compound filter.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// True Filter - performs no filtering. If combined with the NotFilter below, would filter every publication.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TrueLitFilter : public LiteratureFilter {

public:

  TrueLitFilter() = default;
  TrueLitFilter(const TrueLitFilter&) = default;
  ~TrueLitFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary&) const override { return true; }

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<TrueLitFilter>(); }

private:

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// False Filter - unconditionally filters all variants. Filters out all publications.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FalseLitFilter : public LiteratureFilter {

public:

  FalseLitFilter() = default;
  FalseLitFilter(const FalseLitFilter&) = default;
  ~FalseLitFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary&) const override { return false; }

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<FalseLitFilter>(); }


private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Negation Filter, the logical negation of a supplied filter.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NotLitFilter : public LiteratureFilter {

public:

  explicit NotLitFilter(const LiteratureFilter& filter) : filter_ptr_(filter.clone()) {}
  NotLitFilter(const NotLitFilter& copy) : filter_ptr_(copy.filter_ptr_->clone()) {}
  ~NotLitFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override { return not filter_ptr_->applyFilter(publication); }

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<NotLitFilter>(*this); }

private:

  std::unique_ptr<LiteratureFilter> filter_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// And Filter, logical and of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AndLitFilter : public LiteratureFilter {

public:

  AndLitFilter(const LiteratureFilter& filter1, const LiteratureFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {}
  AndLitFilter(const AndLitFilter& copy) : filter1_ptr_(copy.filter1_ptr_->clone()), filter2_ptr_(copy.filter2_ptr_->clone()) {}
  ~AndLitFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override { return filter1_ptr_->applyFilter(publication) and filter2_ptr_->applyFilter(publication); }

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<AndLitFilter>(*this); }

private:

  std::unique_ptr<LiteratureFilter> filter1_ptr_;
  std::unique_ptr<LiteratureFilter> filter2_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Or Filter, logical or of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class OrLitFilter : public LiteratureFilter {

public:

  OrLitFilter(const LiteratureFilter& filter1, const LiteratureFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {}
  OrLitFilter(const OrLitFilter& copy) { filter1_ptr_ = copy.filter1_ptr_->clone(); filter2_ptr_ = copy.filter2_ptr_->clone(); }
  ~OrLitFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override { return filter1_ptr_->applyFilter(publication) or filter2_ptr_->applyFilter(publication); }

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<OrLitFilter>(*this); }

private:

  std::unique_ptr<LiteratureFilter> filter1_ptr_;
  std::unique_ptr<LiteratureFilter> filter2_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter returns true if Plasmodium falciparum found in title, abstract or Mesh codes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PlasmodiumFilter : public LiteratureFilter {

public:

  PlasmodiumFilter() : pf_filter_(OrLitFilter(OrLitFilter(TitleTextFilter(search_text_), AbstractTextFilter(search_text_)), MeSHFilter(Mesh_codes_)).clone()) {}
  PlasmodiumFilter(const PlasmodiumFilter& copy) : pf_filter_(copy.pf_filter_->clone()) {}
  ~PlasmodiumFilter() override = default;

  [[nodiscard]] bool applyFilter(const PublicationSummary& publication) const override { return pf_filter_->applyFilter(publication); }

  [[nodiscard]] std::unique_ptr<LiteratureFilter> clone() const override { return std::make_unique<PlasmodiumFilter>(*this); }

private:

  std::unique_ptr<LiteratureFilter> pf_filter_;

  inline static const std::vector<std::string> Mesh_codes_{  "D010963" /* Plasmodium falciparum */
                                                             ,"D016778" /*  Malaria, Falciparum */
                                                             , "D008288" /* Malaria */ };

  inline static const std::vector<std::string> search_text_ { "Plasmodium",
                                                              "plasmodium",
                                                              "Falciparum",
                                                              "falciparum",
                                                              "Malaria",
                                                              "malaria" };

};



}   // end namespace


#endif //KGL_LITERATURE_FILTER_H
