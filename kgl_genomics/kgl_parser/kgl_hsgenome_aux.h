//
// Created by kellerberrin on 29/6/21.
//

#ifndef KGL_HSGENOME_AUX_H
#define KGL_HSGENOME_AUX_H


#include "kgl_variant_db_type.h"
#include "kgl_square_parser.h"
#include "kgl_resource_db.h"

namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information data structure.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class HsGenealogyRecord;

enum class AuxSexType { MALE, FEMALE };

class HsGenomeAuxRecord {

public:

  explicit HsGenomeAuxRecord(const HsGenealogyRecord& ped_record);
  HsGenomeAuxRecord(std::string individual_id,
                    const std::string& sex_text,
                    std::string super_population,
                    std::string super_description,
                    std::string population,
                    std::string population_description) :
  individual_id_(std::move(individual_id)),
  population_(std::move(population)),
  population_description_(std::move(population_description)),
  super_population_(std::move(super_population)),
  super_description_(std::move(super_description)) {

    sex_ = (sex_text == MALE_ ? AuxSexType::MALE : AuxSexType::FEMALE);

  }
  ~HsGenomeAuxRecord() = default;

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


  constexpr static const char* MALE_ = "XY";
  constexpr static const char* FEMALE_ = "XX";


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information pure virtual object.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using HsAuxPopulationMap = std::map<std::string, std::string>; // Population (key) and description (value).

class HsGenomeAux : public ResourceBase {

public:

  explicit HsGenomeAux(std::string identifier) : ResourceBase(std::move(identifier)) {}
  ~HsGenomeAux() override = default;

  [[nodiscard]] virtual const HsAuxPopulationMap& populationList() const = 0;
  [[nodiscard]] virtual const HsAuxPopulationMap& superPopulationList() const = 0;
  [[nodiscard]] virtual std::vector<GenomeId_t> getGenomeList() const = 0;
  [[nodiscard]] virtual std::optional<HsGenomeAuxRecord> getGenome(const std::string& genome) const = 0;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information implementation object.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



using HsGenomeAuxRecordMap = std::map<std::string, const HsGenomeAuxRecord>;

class HsGenomeAuxData : public HsGenomeAux {

public:

  explicit HsGenomeAuxData(std::string genome_aux_ident) : HsGenomeAux(std::move(genome_aux_ident)) {}
  ~HsGenomeAuxData() override = default;


  // Re-create the population lists any time after construction.
  void refreshPopulationLists();
  bool addGenomeAuxRecord(const HsGenomeAuxRecord& record);

  // Virtual functions.
  [[nodiscard]] const HsAuxPopulationMap& populationList() const override { return population_list_; }
  [[nodiscard]] const HsAuxPopulationMap& superPopulationList() const override { return super_population_list_; }
  [[nodiscard]] std::vector<GenomeId_t> getGenomeList() const override;
  [[nodiscard]] std::optional<HsGenomeAuxRecord> getGenome(const std::string& genome) const override;
  [[nodiscard]] RuntimeResourceType getResourceType() const override { return RuntimeResourceType::GENOME_AUX_INFO; }

private:

  HsGenomeAuxRecordMap genome_aux_record_map_;
  HsAuxPopulationMap population_list_;
  HsAuxPopulationMap super_population_list_;

  [[nodiscard]] const HsGenomeAuxRecordMap& getMap() const { return genome_aux_record_map_; }

};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information parser object.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class ParseHsGenomeAuxFile  {

public:

  explicit ParseHsGenomeAuxFile(const std::shared_ptr<HsGenomeAuxData>& genealogy_data) : genome_aux_data_(genealogy_data) {}
  ~ParseHsGenomeAuxFile() = default;

  void readParseImpl(const std::string& file_name);

private:

  std::shared_ptr<HsGenomeAuxData> genome_aux_data_;
  SquareTextParser text_parser_;

  static constexpr const char FIELD_DELIMITER_CHAR_{'\t'};   // File field separator (char).
  static constexpr const size_t GENOME_AUX_FIELD_COUNT_{7};  // File contains 7 tab delimited fields.


  bool moveToRecord(const std::vector<std::string>& record_fields);

};




} // namespace


#endif //KGL_HSGENOME_AUX_H
