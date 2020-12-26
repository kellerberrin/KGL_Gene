//
// Created by kellerberrin on 31/7/20.
//

#ifndef KGL_DATA_BASE_H
#define KGL_DATA_BASE_H

#include "kel_exec_env.h"

#include <string>
#include <vector>
#include <optional>

namespace kellerberrin::genome {   //  organization::project


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object links the data descriptor in the XML files to the required data parser.
// If any data file types are added then this object must be updated (and possibly a new parser written).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Particular genetic data sources.
enum class DataSourceEnum { Genome1000,
                            GnomadGenome3_1,
                            Falciparum,
                            GnomadExomes3_1,
                            GnomadExomes2_1,
                            Gnomad3_1,
                            Gnomad3_0,
                            Gnomad2_1,
                            PedGnome1000,
                            Clinvar,
                            dbSNP,
                            NotImplemented};   // Error condition if the data source is not found.

// Parsers available for genetic sources.
enum class ParserTypeEnum { DiploidPhased,
                            DiploidFalciparum,
                            DiploidGnomad,
                            MonoGenomeUnphased,
                            PedGenome1000};

// The conceptual structure of the genetic information.
enum class DataStructureEnum { DiploidPhased,   // Phased Diploid Genome1000 only (PopulationVariant)
                               DiploidUnphased,  // Unphased Diploid GnomadGenome3_1 (PopulationVariant)
                               UnphasedMonoGenome, // Genomic data that contains allele information (PopulationVariant)
                               PedGenome1000};  // Additional data to complement the Genome1000 data. (GenomePEDData)

// The actual C++ implementation of the data type. Used for casting from the DataObjectBase class.
enum class DataImplEnum { PopulationVariant,  GenomePEDData };


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A structure to hold a vector of data source characteristics, these are typically indexed by the source_text in the
// XML data definition files.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct DataCharacteristic {

  const std::string source_text;
  const DataSourceEnum data_source;
  const ParserTypeEnum parser_type;
  const DataStructureEnum data_structure;
  const DataImplEnum data_implementation;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Information about data files (populations) returned from the data parsers.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DataObjectBase {

public:

  DataObjectBase(std::string data_identifier, DataSourceEnum data_source) : data_identifier_(std::move(data_identifier)),
                                                                            data_source_(data_source) {}
  virtual ~DataObjectBase() = default;

  [[nodiscard]] const std::string& fileId() const { return data_identifier_; }
  [[nodiscard]] DataSourceEnum dataSource() const { return data_source_; }
  [[nodiscard]] DataCharacteristic dataCharacteristic() const {

    for (auto const& type : data_characteristics_) {

      if (data_source_ == type.data_source) return type;

    }

    // Should never happen.
    ExecEnv::log().critical("DataObjectBase::dataCharacteristic; critical error unknown population type, program terminates.");
    return data_characteristics_[0];

  }

  void setPopulationId(const std::string& data_id) { data_identifier_ = data_id; }

  static std::optional<DataCharacteristic> findCharacteristic(const std::string& file_type) {

    for (auto const& type : data_characteristics_) {

      if (file_type == type.source_text) return type;

    }

    return std::nullopt;

  }


private:

// The data file name.
  std::string data_identifier_;
  DataSourceEnum data_source_;

// DataCharacteristic vector.
  static const std::vector<DataCharacteristic> data_characteristics_;

};


// Data type characteristics vector.
// Additional data file types must be added here.
inline const std::vector<DataCharacteristic>  DataObjectBase::data_characteristics_ = {

    { "Genome1000", DataSourceEnum::Genome1000, ParserTypeEnum::DiploidPhased, DataStructureEnum::DiploidPhased, DataImplEnum::PopulationVariant},
    { "GnomadGenome3_1", DataSourceEnum::GnomadGenome3_1, ParserTypeEnum::DiploidGnomad, DataStructureEnum::DiploidUnphased, DataImplEnum::PopulationVariant },
    { "Falciparum", DataSourceEnum::Falciparum, ParserTypeEnum::DiploidFalciparum, DataStructureEnum::DiploidUnphased, DataImplEnum::PopulationVariant },
    { "GnomadExomes3_1", DataSourceEnum::GnomadExomes3_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant },
    { "GnomadExomes2_1", DataSourceEnum::GnomadExomes2_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant },
    { "Gnomad3_1", DataSourceEnum::Gnomad3_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant },
    { "Gnomad3_0", DataSourceEnum::Gnomad3_0, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant },
    { "Gnomad2_1", DataSourceEnum::Gnomad2_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant },
    { "Clinvar", DataSourceEnum::Clinvar, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant },
    { "dbSNP", DataSourceEnum::dbSNP, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant },
    { "PedGnome1000", DataSourceEnum::PedGnome1000, ParserTypeEnum::PedGenome1000, DataStructureEnum::PedGenome1000, DataImplEnum::GenomePEDData },

};



} // namespace


#endif //KGL_KGL_DATA_FILE_H
