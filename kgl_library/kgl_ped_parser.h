//
// Created by kellerberrin on 1/8/20.
//

#ifndef KGL_PED_PARSER_H
#define KGL_PED_PARSER_H


#include "kgl_data_base.h"
#include "kgl_data_file_impl.h"

namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file data structure.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PEDRecord {

public:

  PEDRecord(std::string family_id,
            std::string individual_id,
            std::string paternal_id,
            std::string maternal_id,
            std::string sex,
            std::string pheno_type,
            std::string population,
            std::string relationship,
            std::string siblings,
            std::string second_order,
            std::string third_order,
            std::string comments) : family_id_(std::move(family_id)),
                                    individual_id_(std::move(individual_id)),
                                    paternal_id_(std::move(paternal_id)),
                                    maternal_id_(std::move(maternal_id)),
                                    sex_(std::move(sex)),
                                    pheno_type_(std::move(pheno_type)),
                                    population_(std::move(population)),
                                    relationship_(std::move(relationship)),
                                    siblings_(std::move(siblings)),
                                    second_order_(std::move(second_order)),
                                    third_order_(std::move(third_order)),
                                    comments_(std::move(comments)) {}

  PEDRecord(const PEDRecord&) = default;
  ~PEDRecord() = default;

  [[nodiscard]] const std::string& familyId() const { return family_id_; }
  [[nodiscard]] const std::string& individualId() const { return individual_id_; }
  [[nodiscard]] const std::string& paternalId() const { return paternal_id_; }
  [[nodiscard]] const std::string& maternalId() const { return maternal_id_; }
  [[nodiscard]] const std::string& sex() const { return sex_; }
  [[nodiscard]] const std::string& phenoType() const { return pheno_type_; }
  [[nodiscard]] const std::string& population() const { return population_; }
  [[nodiscard]] const std::string& relationship() const { return relationship_; }
  [[nodiscard]] const std::string& siblings() const { return siblings_; }
  [[nodiscard]] const std::string& secondOrder() const { return second_order_; }
  [[nodiscard]] const std::string& thirdOrder() const { return third_order_; }
  [[nodiscard]] const std::string& comments() const { return comments_; }

  [[nodiscard]] static size_t PEDFieldCount() { return PED_FIELD_COUNT_; }

private:

  std::string family_id_;
  std::string individual_id_;
  std::string paternal_id_;
  std::string maternal_id_;
  std::string sex_;
  std::string pheno_type_;
  std::string population_;
  std::string relationship_;
  std::string siblings_;
  std::string second_order_;
  std::string third_order_;
  std::string comments_;

  static constexpr const size_t PED_FIELD_COUNT_{12};         // Threads parsing PED records

};


using PEDRecordMap = std::map<std::string, const PEDRecord>;

class GenomePEDData : public DataObjectBase {

public:

  explicit GenomePEDData(std::string ped_ident) : DataObjectBase(std::move(ped_ident)) {}
  ~GenomePEDData() override = default;

  bool addPEDRecord(const PEDRecord& record) {

    auto [iter, result] = PED_record_map_.try_emplace(record.individualId(), record);

    if (not result) {

      ExecEnv::log().error("GenomePEDData::addPEDRecord, could add PED for sample: {} (duplicate)", record.individualId());
      return false;

    }

    return true;

  }

  [[nodiscard]] const PEDRecordMap& getMap() const { return PED_record_map_; }

  [[nodiscard]] DataTypeEnum dataType() const override { return DataTypeEnum::PedAncestor; }

private:

  PEDRecordMap PED_record_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file parser.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ParsePedFile : private FileDataIO {

public:

  explicit ParsePedFile(std::shared_ptr<GenomePEDData> ped_data,
                        std::string file_name) : FileDataIO(std::move(file_name)), ped_data_(ped_data) {}
  ~ParsePedFile() override = default;

  void readParsePEDImpl();

private:

  std::shared_ptr<GenomePEDData> ped_data_;
  static constexpr const long PARSER_THREADS_{1};         // Threads parsing PED records
  static constexpr const char PED_FIELD_DELIMITER_CHAR_{'\t'};   // PED Field separator (char).

  bool moveToPEDRecord(std::string&& line_record);

};



} // namespace

#endif //KGL_KGL_PED_PARSER_H
