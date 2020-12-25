//
// Created by kellerberrin on 31/7/20.
//

#ifndef KGL_DATA_BASE_H
#define KGL_DATA_BASE_H

#include <string>
#include <vector>

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
                            Gnomad3_1,
                            Gnomad2_1,
                            PedGnome1000,
                            Clinvar,
                            dbSNP,
                            NotImplemented};   // Error condition if the data source is not found.

// Parsers available for genetic sources.
enum class ParserTypeEnum { DiploidPopulation,
                            DiploidFalciparum,
                            MonoGenomeUnphased,
                            PedGenome1000};

// The structure of the genetic information.
enum class DataStructureEnum { DiploidPhased,   // Phased Diploid Genome1000 only
                               DiploidUnphased,  // Unphased Diploid GnomadGenome3_1
                               DiploidFalciparum,  // All the Falciparum data is unphased diploid to identify complexity of Infection (CoI).
                               UnphasedMonoGenome, // Genomic data that contains allele information (GnomadExomes3_1, Gnomad3_1, Gnomad2_1, etc)
                               PedGenome1000};  // Additional data to complement the Genome1000 data.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A structure to hold a vector of data source characteristics, these are typically indexed by the source_text in the
// XML data definition files.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct DataCharacteristic {

  const char* source_text;
  const DataSourceEnum data_source;
  const ParserTypeEnum parser_type;
  const DataStructureEnum data_structure;

};

// Legacy, to be removed.
enum class DataTypeEnum { DiploidPopulation, HaploidPopulation, UnphasedPopulation, PedAncestor, NotImplemented};


class DataObjectBase {

public:

  DataObjectBase(std::string identifier) : data_identifier_(std::move(identifier)) {}
  virtual ~DataObjectBase() = default;

  [[nodiscard]] const std::string& Id() const { return data_identifier_; }
  void setId(const std::string& data_id) { data_identifier_ = data_id; }

  virtual DataTypeEnum dataType() const = 0;



private:

// The data file name.
   std::string data_identifier_;

// DataCharacteristic vector.
  static const std::vector<DataCharacteristic> data_characteristics_;

};


// Data type characteristics vector.
// Additional data file types must be added here.
inline const std::vector<DataCharacteristic>  DataObjectBase::data_characteristics_ = {

    { "Genome1000", DataSourceEnum::Genome1000, ParserTypeEnum::DiploidPopulation, DataStructureEnum::DiploidPhased },
    { "GnomadGenome3_1", DataSourceEnum::GnomadGenome3_1, ParserTypeEnum::DiploidPopulation, DataStructureEnum::DiploidUnphased },
    { "Falciparum", DataSourceEnum::Falciparum, ParserTypeEnum::DiploidFalciparum, DataStructureEnum::DiploidUnphased },
    { "GnomadExomes3_1", DataSourceEnum::GnomadExomes3_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome },
    { "Gnomad3_1", DataSourceEnum::Gnomad3_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome },
    { "Gnomad2_1", DataSourceEnum::Gnomad2_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome },
    { "Clinvar", DataSourceEnum::Clinvar, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome },
    { "dbSNP", DataSourceEnum::dbSNP, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome },
    { "PedGnome1000", DataSourceEnum::PedGnome1000, ParserTypeEnum::PedGenome1000, DataStructureEnum::PedGenome1000 },

};



} // namespace


#endif //KGL_KGL_DATA_FILE_H
