//
// Created by kellerberrin on 31/1/21.
//

#include "kgl_analysis_mutation_gene_ethnic.h"


namespace kgl = kellerberrin::genome;




bool kgl::GeneEthnicitySex::pedAnalysis(const GenomeId_t& genome_id,
                                      size_t count,
                                      const std::shared_ptr<const GenomePEDData>& ped_data) {


  if (not ped_data) {

    ExecEnv::log().critical("GenomeMutation::pedAnalysis; PED file pointer net defined");
    return false;

  }

  total_ += count;

  auto result = ped_data->getMap().find(genome_id);

  if (result == ped_data->getMap().end()) {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; Genome sample: {} does not have a PED record", genome_id);
    return false;

  }

  auto const& [sample_id, ped_record] = *result;

  if (ped_record.sexType() == PedSexType::MALE) {

    male_ += count;

  } else {

    female_ += count;

  }

  std::string super_pop = ped_record.superPopulation();
  auto find_result = super_population_.find(super_pop);
  if (find_result != super_population_.end()) {

    auto& [super_pop, super_count] = *find_result;
    super_count += count;

  } else {

    ExecEnv::log().error("GenomeMutation::pedAnalysis; Unable to find super population: {}", super_pop);

  }

  std::string pop = ped_record.population();
  auto pop_result = population_.find(pop);
  if (pop_result != population_.end()) {

    auto& [pop, pop_count] = *pop_result;
    pop_count += count;

  } else {

    ExecEnv::log().error("GenomeMutation::pedAnalysis; Unable to find population: {}", pop);

  }

  return true;

}


void kgl::GeneEthnicitySex::updatePopulations(const std::shared_ptr<const GenomePEDData>& ped_data) {

  // Generate the template populations.
  std::map<std::string, size_t> population_map;
  for (auto const& [population, description] : ped_data->populationList()) {

    population_map[population] = 0;

  }

  // Generate the template super populations.
  std::map<std::string, size_t> super_population_map;
  for (auto const& [super_population, description] : ped_data->superPopulationList()) {

    super_population_map[super_population] = 0;

  }


  super_population_ = super_population_map;
  population_ = population_map;

}


void kgl::GeneEthnicitySex::writeHeader( const std::shared_ptr<const GenomePEDData>& ped_data,
                                         std::ostream& out_file,
                                         char output_delimiter) const {

  if (display_flags_ & DISPLAY_SEX_FLAG) {

    writeSexHeader(out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_SUPER_POP_FLAG) {

    writeSuperPopHeader( ped_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_POPULATION_FLAG) {

    writePopHeader( ped_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }

}



void kgl::GeneEthnicitySex::writeSexHeader( std::ostream& out_file, char output_delimiter) const {

  out_file << header_prefix_ << "Male" << output_delimiter
           << header_prefix_ << "Female";

}

void kgl::GeneEthnicitySex::writeSuperPopHeader( const std::shared_ptr<const GenomePEDData>& ped_data,
                                                 std::ostream& out_file,
                                                 char output_delimiter) const {

  for (auto const& [population, description] : ped_data->superPopulationList()) {

    if (population != ped_data->superPopulationList().begin()->first) {

      out_file << output_delimiter;

    }

    out_file << header_prefix_ << population;


  }

}

void kgl::GeneEthnicitySex::writePopHeader( const std::shared_ptr<const GenomePEDData>& ped_data,
                                            std::ostream& out_file,
                                            char output_delimiter) const {

  for (auto const& [population, description] : ped_data->populationList()) {

    if (population != ped_data->populationList().begin()->first) {

      out_file << output_delimiter;

    }

    out_file << header_prefix_ << population;


  }

}



void kgl::GeneEthnicitySex::writeOutput(const std::shared_ptr<const GenomePEDData>& ped_data,
                                        std::ostream& out_file,
                                        char output_delimiter) const {

  if (display_flags_ & DISPLAY_SEX_FLAG) {

    kgl::GeneEthnicitySex::writeSex(out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_SUPER_POP_FLAG) {

    writeSuperPop( ped_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }

  if (display_flags_ & DISPLAY_POPULATION_FLAG) {

    writePop( ped_data, out_file, output_delimiter);
    out_file << output_delimiter;

  }


}


void kgl::GeneEthnicitySex::writeSex( std::ostream& out_file,
                                      char output_delimiter) const {

  out_file << male_ << output_delimiter
           << female_;

}


void kgl::GeneEthnicitySex::writeSuperPop( const std::shared_ptr<const GenomePEDData>& ped_data,
                                           std::ostream& out_file,
                                           char output_delimiter) const {

  if (super_population_.size() != ped_data->superPopulationList().size()) {

    ExecEnv::log().error( "GeneEthnicitySex::writeSuperPop; Mismatch between data super population size: {}, and Ped size: {}",
                          super_population_.size(), ped_data->superPopulationList().size());

  }

  for (auto const& [super_pop, super_pop_count] : super_population_) {

    if (super_pop != super_population_.begin()->first) {

      out_file << output_delimiter;

    }

    out_file << super_pop_count;

  }

}


void kgl::GeneEthnicitySex::writePop( const std::shared_ptr<const GenomePEDData>& ped_data,
                                      std::ostream& out_file,
                                      char output_delimiter) const {

  if (population_.size() != ped_data->populationList().size()) {

    ExecEnv::log().error("GeneEthnicitySex::writePop; Mismatch between data super population size: {}, and Ped size: {}",
                         population_.size(), ped_data->populationList().size());

  }

  for (auto const& [pop, pop_count] : population_) {

    if (pop != population_.begin()->first) {

      out_file << output_delimiter;

    }

    out_file << pop_count;

  }

}


