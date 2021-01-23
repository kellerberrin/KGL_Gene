//
// Created by kellerberrin on 1/8/20.
//

#ifndef KGL_PED_PARSER_H
#define KGL_PED_PARSER_H


#include "kgl_variant_db_type.h"
#include "kgl_square_parser.h"

namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file data structure.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class PedSexType { MALE, FEMALE };

class PEDRecord {

public:

  PEDRecord(std::string family_id,
            std::string individual_id,
            std::string paternal_id,
            std::string maternal_id,
            std::string sex,
            std::string pheno_type,
            std::string population,
            std::string population_description,
            std::string super_population,
            std::string super_description,
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
                                    population_description_(std::move(population_description)),
                                    super_population_(std::move(super_population)),
                                    super_description_(std::move(super_description)),
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
  [[nodiscard]] PedSexType sexType() const { return (sex_ == MALE_ ? PedSexType::MALE : PedSexType::FEMALE); }
  [[nodiscard]] const std::string& phenoType() const { return pheno_type_; }
  [[nodiscard]] const std::string& population() const { return population_; }
  [[nodiscard]] const std::string& populationDescription() const { return population_description_; }
  [[nodiscard]] const std::string& superPopulation() const { return super_population_; }
  [[nodiscard]] const std::string& superDescription() const { return super_description_; }
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
  std::string population_description_;
  std::string super_population_;
  std::string super_description_;
  std::string relationship_;
  std::string siblings_;
  std::string second_order_;
  std::string third_order_;
  std::string comments_;

  static constexpr const size_t PED_FIELD_COUNT_{15};         // Threads parsing PED records

  constexpr static const char* MALE_ = "1";
  constexpr static const char* FEMALE_ = "2";


};


using PEDRecordMap = std::map<std::string, const PEDRecord>;

class GenomePEDData : public DataDB {

public:

  explicit GenomePEDData(std::string ped_ident, DataSourceEnum data_source) : DataDB(data_source),
                                                                              ped_ident_(std::move(ped_ident)) {}
  ~GenomePEDData() override = default;

  [[nodiscard]] const std::string& fileId() const override { return ped_ident_; }
  [[nodiscard]] const PEDRecordMap& getMap() const { return PED_record_map_; }

  bool addPEDRecord(const PEDRecord& record);

private:

  PEDRecordMap PED_record_map_;
  std::string ped_ident_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file parser.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ParsePedFile  {

public:

  explicit ParsePedFile(std::shared_ptr<GenomePEDData> ped_data) : ped_data_(ped_data) {}
  ~ParsePedFile() = default;

  void readParsePEDImpl(const std::string& file_name);

private:

  std::shared_ptr<GenomePEDData> ped_data_;
  SquareTextParser text_parser_;

  static constexpr const char PED_FIELD_DELIMITER_CHAR_{'\t'};   // PED Field separator (char).


  bool moveToPEDRecord(const std::vector<std::string>& record_fields);

};



} // namespace

#endif //KGL_KGL_PED_PARSER_H
