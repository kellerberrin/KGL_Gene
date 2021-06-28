//
// Created by kellerberrin on 1/8/20.
//

#ifndef KGL_PED_PARSER_H
#define KGL_PED_PARSER_H


#include "kgl_variant_db_type.h"
#include "kgl_square_parser.h"
#include "kgl_resource_db.h"

namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information data structure.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GenealogyRecord;

enum class AuxSexType { MALE, FEMALE };

class GenomeAuxRecord {

public:

  GenomeAuxRecord(const GenealogyRecord& ped_record);
  ~GenomeAuxRecord() = default;

  [[nodiscard]] const std::string& individualId() const { return individual_id_; }
  [[nodiscard]] AuxSexType sexType() const { return sex_; }
  [[nodiscard]] const std::string& population() const { return population_; }
  [[nodiscard]] const std::string& populationDescription() const { return population_description_; }
  [[nodiscard]] const std::string& superPopulation() const { return super_population_; }
  [[nodiscard]] const std::string& superDescription() const { return super_description_; }

private:

  std::string individual_id_;
  AuxSexType sex_;
  std::string population_;
  std::string population_description_;
  std::string super_population_;
  std::string super_description_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information pure virtual object.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using AuxPopulationMap = std::map<std::string, std::string>; // Population (key) and description (value).

class GenomeAuxInfo : public ResourceBase {

public:

  GenomeAuxInfo() = default;
  ~GenomeAuxInfo() override = default;

  [[nodiscard]] virtual const AuxPopulationMap& populationList() const = 0;
  [[nodiscard]] virtual const AuxPopulationMap& superPopulationList() const = 0;
  [[nodiscard]] virtual std::vector<GenomeId_t> getGenomeList() const = 0;
  [[nodiscard]] virtual std::optional<GenomeAuxRecord> getGenome(const std::string& genome) const = 0;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file data structure.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class GenealogyRecord {

public:

  GenealogyRecord(std::string family_id,
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

  GenealogyRecord(const GenealogyRecord&) = default;
  ~GenealogyRecord() = default;

  [[nodiscard]] const std::string& familyId() const { return family_id_; }
  [[nodiscard]] const std::string& individualId() const { return individual_id_; }
  [[nodiscard]] const std::string& paternalId() const { return paternal_id_; }
  [[nodiscard]] const std::string& maternalId() const { return maternal_id_; }
  [[nodiscard]] const std::string& sex() const { return sex_; }
  [[nodiscard]] AuxSexType sexType() const { return (sex_ == MALE_ ? AuxSexType::MALE : AuxSexType::FEMALE); }
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


using GenealogyRecordMap = std::map<std::string, const GenealogyRecord>;

class GenomeGenealogyData : public GenomeAuxInfo {

public:

  explicit GenomeGenealogyData(std::string genealogy_ident) : genealogy_ident_(std::move(genealogy_ident)) {}
  ~GenomeGenealogyData() override = default;

  [[nodiscard]] const std::string& genealogyId() const { return genealogy_ident_; }
  [[nodiscard]] std::optional<GenealogyRecord> getGenomePedRecord(const std::string& genome) const;

  bool addPEDRecord(const GenealogyRecord& record);

  // Re-create the population lists any time after construction.
  void refreshPopulationLists();

  // Virtual functions.
  [[nodiscard]] const AuxPopulationMap& populationList() const override { return population_list_; }
  [[nodiscard]] const AuxPopulationMap& superPopulationList() const override { return super_population_list_; }
  [[nodiscard]] std::vector<GenomeId_t> getGenomeList() const override;
  [[nodiscard]] std::optional<GenomeAuxRecord> getGenome(const std::string& genome) const override;
  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::GENOME_GENEALOGY; }

private:

  GenealogyRecordMap genealogy_record_map_;
  std::string genealogy_ident_;
  AuxPopulationMap population_list_;
  AuxPopulationMap super_population_list_;

  [[nodiscard]] const GenealogyRecordMap& getMap() const { return genealogy_record_map_; }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The PED file parser.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ParseGenomeGenealogyFile  {

public:

  explicit ParseGenomeGenealogyFile(std::shared_ptr<GenomeGenealogyData> genealogy_data) : genealogy_data_(genealogy_data) {}
  ~ParseGenomeGenealogyFile() = default;

  void readParseImpl(const std::string& file_name);

private:

  std::shared_ptr<GenomeGenealogyData> genealogy_data_;
  SquareTextParser text_parser_;

  static constexpr const char FIELD_DELIMITER_CHAR_{'\t'};   // PED Field separator (char).


  bool moveToRecord(const std::vector<std::string>& record_fields);

};



} // namespace

#endif //KGL_KGL_PED_PARSER_H
