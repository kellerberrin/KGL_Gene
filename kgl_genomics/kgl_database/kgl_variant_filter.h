///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant.h"
#include "kgl_runtime_resource.h"

#include "kel_utility.h"

#include <unordered_set>


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// These objects filter on Info field criteria specified in the supplied filter lambda.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using InfoFilterLambda = std::function<bool(const InfoDataVariant&)>;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper class finds the Info field (if it exists).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoFilterImpl  {

public:

  explicit InfoFilterImpl(std::string field_name) : field_name_(std::move(field_name)) {}
  InfoFilterImpl(const InfoFilterImpl&) = default;
  ~InfoFilterImpl() = default;

  [[nodiscard]] bool applyFilter(const InfoFilterLambda& filter_lambda, const Variant& variant) const;

  [[nodiscard]] std::string fieldName() const { return field_name_; }

private:

  const std::string field_name_;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// General Info filter class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoFilter : public VariantFilter {

public:

  InfoFilter(const std::string& field_name, InfoFilterLambda filter_lambda)
  : info_filter_(field_name), filter_lambda_(std::move(filter_lambda)) {

    filterName(field_name);

  }
  InfoFilter(const InfoFilter&) = default;
  ~InfoFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::InfoFilter; }

private:

  const InfoFilterImpl info_filter_;
  const InfoFilterLambda filter_lambda_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter on a Vep subfield found in the Gnomad Homosapien data. Will silently return false for all other data.
//
// If the vep field contains the specified sub-string then the filter returns 'true'.
// If the the empty string "" is specified then the corresponding vep field must be empty to return true.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VepSubStringFilter : public VariantFilter {

public:

  VepSubStringFilter(std::string vep_field_name, std::string sub_string)  // How the filter responds if the data is missing.
  : vep_field_name_(std::move(vep_field_name)),  sub_string_(std::move(sub_string)) {

    std::stringstream ss;
    ss << "Vep Info SubField: " << vep_field_name_ << " contains sub string \"" << sub_string_ <<"\"";
    filterName(ss.str());

  }
  VepSubStringFilter(const VepSubStringFilter&) = default;
  ~VepSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<VepSubStringFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::VepSubStringFilter; }

private:

  const std::string vep_field_name_;
  const std::string sub_string_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Integer Info Filter. Returns true if a scalar integer and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQIntegerFilter : public VariantFilter {

public:

  InfoGEQIntegerFilter(const std::string& field_name, int64_t comparison_value)  // How the filter responds if the data is missing.
  : info_filter_(field_name), comparison_value_(comparison_value) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " >= " << comparison_value;
    filterName(ss.str());

  }
  InfoGEQIntegerFilter(const InfoGEQIntegerFilter&) = default;
  ~InfoGEQIntegerFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoGEQIntegerFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::InfoGEQIntegerFilter; }


private:

  const InfoFilterImpl info_filter_;
  const int64_t comparison_value_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_integer_vector = std::get_if<std::vector<int64_t>>(&data_variant);
    if (p_integer_vector != nullptr) {

      if (p_integer_vector->size() == 1) {

        return (p_integer_vector->front() >= this->comparison_value_);

      } else {

        return false;

      }

    } else {

      return false;

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Float Info Filter. Returns true if a scalar float and greater or equal to the comparison value.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoGEQFloatFilter : public VariantFilter {

public:

  InfoGEQFloatFilter(const std::string& field_name, double comparison_value)  // How the filter responds if the data is missing.
  : info_filter_(field_name), comparison_value_(comparison_value) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " >= " << comparison_value;
    filterName(ss.str());

  }
  InfoGEQFloatFilter(const InfoGEQFloatFilter&) = default;
  ~InfoGEQFloatFilter() override = default;


  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoGEQFloatFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::InfoGEQFloatFilter; }

private:

  const InfoFilterImpl info_filter_;
  const double comparison_value_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_float_vector = std::get_if<std::vector<double>>(&data_variant);
    if (p_float_vector != nullptr) {

      if (p_float_vector->size() == 1) {

        return (p_float_vector->front() >= this->comparison_value_);

      } else {

        return false;

      }

    } else {

      return false;

    }

  };

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// String Info Filter. Returns true if a string Info field contains a specified substring (including itself).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoSubStringFilter : public VariantFilter {

public:

  InfoSubStringFilter(const std::string& field_name, std::string sub_string, bool case_sensitive = false)  // How the filter responds if the data is missing.
  : info_filter_(field_name), sub_string_(std::move(sub_string)), case_sensitive_(case_sensitive) {

    std::stringstream ss;
    ss << "VCF Info Field: " << field_name << " contains sub string \"" << sub_string <<"\"";
    filterName(ss.str());
    if (not case_sensitive_) {

      sub_string_ = Utility::toupper(sub_string_);

    }

  }
  InfoSubStringFilter(const InfoSubStringFilter&) = default;
  ~InfoSubStringFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoSubStringFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::InfoSubStringFilter; }

private:

  const InfoFilterImpl info_filter_;
  std::string sub_string_;
  const bool case_sensitive_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_string_vector = std::get_if<std::vector<std::string>>(&data_variant);
    if (p_string_vector != nullptr) {

      if (not p_string_vector->empty()) {

        for (auto const& value : *p_string_vector) {

          bool result;
          if (case_sensitive_) {

            result = value.find(sub_string_) != std::string::npos;

          } else {

            result = Utility::toupper(value).find(sub_string_) != std::string::npos;

          }

          if (result) {

            return true;

          }

        }

        return false;

      } else {

        return false;

      }

    } else {

      return false;

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Boolean Filter. Returns true if a boolean info field is true. Can be used with the Negation filter below.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InfoBooleanFilter : public VariantFilter {

public:

  explicit InfoBooleanFilter(const std::string& field_name)  // How the filter responds if the data is missing.
  : info_filter_(field_name) {

    std::stringstream ss;
    ss << "VCF Boolean Info Field: " << field_name;
    filterName(ss.str());

  }
  InfoBooleanFilter(const InfoBooleanFilter&) = default;
  ~InfoBooleanFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return info_filter_.applyFilter(filter_lambda_, variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<InfoBooleanFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::InfoBooleanFilter; }

private:

  const InfoFilterImpl info_filter_;
  const InfoFilterLambda filter_lambda_ = [this] (const InfoDataVariant& data_variant) -> bool {

    auto p_bool = std::get_if<bool>(&data_variant);
    if (p_bool != nullptr) {

      return *p_bool;

    } else {

      return false;

    }

  };

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RefAltCountFilter : public VariantFilter {

public:

  explicit RefAltCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    std::stringstream ss;
    ss << "Variants with minimum Ref+Alt base count:" << minimum_count_;
    filterName(ss.str());

  }
  RefAltCountFilter(const RefAltCountFilter&) = default;
  ~RefAltCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<RefAltCountFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::RefAltCountFilter; }


private:

  size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a specified minimum DP counts.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DPCountFilter : public VariantFilter {

public:

  explicit DPCountFilter(size_t minimum_count) : minimum_count_(minimum_count) {

    std::stringstream ss;
    ss << "Variants with minimum DP base count:" << minimum_count_;
    filterName(ss.str());

  }
  DPCountFilter(const DPCountFilter&) = default;
  ~DPCountFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<DPCountFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::DPCountFilter; }

private:

  size_t minimum_count_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unique variants disregarding phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UniqueUnphasedFilter: public VariantFilter {

public:

  UniqueUnphasedFilter() { filterName("UniqueUnphasedFilter"); }
  UniqueUnphasedFilter(const UniqueUnphasedFilter&) = default;
  ~UniqueUnphasedFilter() override = default;

  UniqueUnphasedFilter& operator=(const UniqueUnphasedFilter&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<UniqueUnphasedFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::UniqueUnphasedFilter; }


private:

  mutable std::unordered_set<std::string> hashed_variants_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unique variants including phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UniquePhasedFilter: public VariantFilter {

public:

  UniquePhasedFilter() { filterName("UniquePhasedFilter"); }
  UniquePhasedFilter(const UniquePhasedFilter&) = default;
  ~UniquePhasedFilter() override = default;

  UniquePhasedFilter& operator=(const UniquePhasedFilter&) = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override;

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<UniquePhasedFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::UniquePhasedFilter; }


private:

  mutable std::unordered_set<std::string> hashed_variants_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ReturnType both variants if and only there are 2 variants at the location that are identical disregarding phase.
// Note If the the variants are phased then the filter will pass back 2 identical variants which may need to be
// further filtered for uniqueness. The inSituFilter version deletes all variant offsets that are not homozygous.
// The homozygous filter must be called first before any other filter and cannot be used in a compound filter.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class HomozygousFilter: public VariantFilter {

public:

  HomozygousFilter() { filterName("Homozygous"); }
  HomozygousFilter(const HomozygousFilter&) = default;
  ~HomozygousFilter() override = default;

  HomozygousFilter& operator=(const HomozygousFilter&) = default;

  // Dummy implementation. Implemented at offset level
  [[nodiscard]] bool applyFilter(const Variant&) const override { return true; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<HomozygousFilter>(); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::HomozygousFilter; }


private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Ensure max 2 variants per offset and properly phased (or unphased).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DiploidFilter: public VariantFilter {

public:

  DiploidFilter() { filterName("DiploidFilter"); }
  DiploidFilter(const DiploidFilter&) = default;
  ~DiploidFilter() override = default;

  DiploidFilter& operator=(const DiploidFilter&) = default;

  // Dummy implementation. Implemented at offset level
  [[nodiscard]] bool applyFilter(const Variant&) const override { return true; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<DiploidFilter>(); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::DiploidFilter; }


private:

};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter a population for genomes.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GenomeFilter: public VariantFilter {

public:

  explicit GenomeFilter(const std::vector<GenomeId_t>& filtered_genomes) {

    filterName("GenomeFilter");

    for (auto const& genome_id : filtered_genomes) {

      filter_genomes_.insert(genome_id);

    }

  }
  GenomeFilter(const GenomeFilter&) = default;
  ~GenomeFilter() override = default;

  GenomeFilter& operator=(const GenomeFilter&) = default;

  // Dummy implementation. Implemented at genome level
  [[nodiscard]] bool applyFilter(const Variant&) const override { return true; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<GenomeFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::GenomeFilter; }

  [[nodiscard]] bool filterGenome(const GenomeId_t& genome_id) const { return filter_genomes_.find(genome_id) != filter_genomes_.end(); }

private:

  std::unordered_set<GenomeId_t> filter_genomes_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants on phasing.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PhaseFilter : public VariantFilter {

public:

  explicit PhaseFilter(VariantPhase phase) : phase_(phase) {

    filterName("PhaseFilter");

  }
  PhaseFilter(const PhaseFilter&) = default;
  ~PhaseFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.phaseId() == phase_; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<PhaseFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::PhaseFilter; }


private:

  const VariantPhase phase_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PassFilter : public VariantFilter {

public:

  explicit PassFilter() { filterName("Variant marked 'Pass' for filters in VCF"); }
  PassFilter(const PassFilter&) = default;
  ~PassFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.evidence().passFilter(); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<PassFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::PassFilter; }


private:


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPFilter : public VariantFilter {

public:

  explicit SNPFilter() { filterName("SNP and MNP Variants"); }
  SNPFilter(const SNPFilter&) = default;
  ~SNPFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.isSNP(); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<SNPFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::SNPFilter; }


private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident) : contig_ident_(contig_ident) {

    filterName(std::string("Variants in Contig: ") + contig_ident_);

  }
  ContigFilter(const ContigFilter&) = default;
  ~ContigFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return implementFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<ContigFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::ContigFilter; }


private:

  const ContigId_t contig_ident_;

  [[nodiscard]] bool implementFilter(const Variant& variant) const;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Region filter - generally used as AndFilter(ContigFilter("ContigName"), RegionFilter(Start, End))
// Uses the half open interval convention [Begin, End) with zero offset (the first contig element has zero offset).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RegionFilter : public VariantFilter {

public:

  explicit RegionFilter(ContigOffset_t start, ContigOffset_t end) : start_(start), end_(end) {

    std::stringstream ss;
    ss << "Variant in the half-interval [" << start_ << ", " << end_ << ")";
    filterName(ss.str());

  }
  RegionFilter(const RegionFilter&) = default;
  ~RegionFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return variant.offset() >= start_ and variant.offset() < end_; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<RegionFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::RegionFilter; }


private:

  ContigOffset_t start_;
  ContigOffset_t end_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// True Filter - performs no filtering. If combined with the NotFilter below, would filter every variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TrueFilter : public VariantFilter {

public:

  explicit TrueFilter() {

    filterName("TrueFilter");

  }
  TrueFilter(const TrueFilter&) = default;
  ~TrueFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant&) const override { return true; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<TrueFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::TrueFilter; }


private:

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// False Filter - unconditionally filters all variants. Useful for deleting large populations.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FalseFilter : public VariantFilter {

public:

  explicit FalseFilter() {

    filterName("FalseFilter");

  }
  FalseFilter(const FalseFilter&) = default;
  ~FalseFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant&) const override { return false; }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<FalseFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::FalseFilter; }


private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Negation Filter, the logical negation of a supplied filter.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NotFilter : public VariantFilter {

public:

  explicit NotFilter(const VariantFilter& filter) : filter_ptr_(filter.clone()) {

    filterName(std::string("NOT(") + filter_ptr_->filterName() + std::string(")"));

  }
  NotFilter(const NotFilter& copy) { filter_ptr_ = copy.filter_ptr_->clone(); }
  ~NotFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return not filter_ptr_->applyFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<NotFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::NotFilter; }


private:


  std::shared_ptr<VariantFilter> filter_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// And Filter, logical and of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AndFilter : public VariantFilter {

public:

  explicit AndFilter(const VariantFilter& filter1, const VariantFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {

    filterName("AND(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")");

  }
  AndFilter(const AndFilter& copy) { filter1_ptr_ = copy.filter1_ptr_->clone(); filter2_ptr_ = copy.filter2_ptr_->clone();}
  ~AndFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) and filter2_ptr_->applyFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<AndFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::AndFilter; }


private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Or Filter, logical or of two supplied filters.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class OrFilter : public VariantFilter {

public:

  explicit OrFilter(const VariantFilter& filter1, const VariantFilter& filter2) : filter1_ptr_(filter1.clone()), filter2_ptr_(filter2.clone()) {

    filterName("OR(" + filter1_ptr_->filterName() + ", " + filter2_ptr_->filterName() + ")");

  }
  OrFilter(const OrFilter& copy ) { filter1_ptr_ = copy.filter1_ptr_->clone(); filter2_ptr_ = copy.filter2_ptr_->clone(); }
  ~OrFilter() override = default;

  [[nodiscard]] bool applyFilter(const Variant& variant) const override { return filter1_ptr_->applyFilter(variant) or filter2_ptr_->applyFilter(variant); }

  [[nodiscard]] std::shared_ptr<VariantFilter> clone() const override { return std::make_shared<OrFilter>(*this); }

  [[nodiscard]] FilterType filterType() const override { return FilterType::OrFilter; }

private:

  std::shared_ptr<VariantFilter> filter1_ptr_;
  std::shared_ptr<VariantFilter> filter2_ptr_;

};



}   // end namespace


#endif //KGL_FILTER_H
