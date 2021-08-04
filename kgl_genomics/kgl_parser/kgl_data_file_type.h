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


// The organism that a data file refers.
enum class DataOrganism { HomoSapien,
                          PlasmodiumFalciparum,
                          PlasmodiumVivax,
                          NoOrganism};   // Data file is not organism specific or organism not defined.

// Particular genetic data sources.
enum class DataSourceEnum { Genome1000,
                            GnomadGenome3_1,
                            Falciparum,
                            GnomadExomes3_1,
                            GnomadExomes2_1,
                            Gnomad3_1,
                            Gnomad3_0,
                            Gnomad2_1,
                            Clinvar,
                            dbSNP,
                            JSONdbSNP,
                            BioPMID,
                            NotImplemented};   // Data source is not specified.

// Parsers available for genetic sources.
enum class ParserTypeEnum { DiploidPhased,
                            DiploidFalciparum,
                            DiploidGnomad,
                            MonoGenomeUnphased,
                            MonoDBSNPUnphased,
                            MonoJSONdbSNPUnphased,
                            ParseBioPMID,
                            FilenameOnly};  // File name is passed to the requesting package - no parsing performed.

// The conceptual structure of the genetic information.
enum class DataStructureEnum { DiploidPhased,   // Phased Diploid Genome1000 only (PopulationDB)
                               DiploidUnphased,  // Unphased Diploid GnomadGenome3_1 (PopulationDB)
                               UnphasedMonoGenome,
                               CitationMap,
                               BioPMIDMap,
                               NoStructure}; // No data passed to the requesting package (only the file name)

// The actual C++ implementation of the data type. Used for casting from the DataDB class.
enum class DataImplEnum { PopulationVariant, PMIDCitationMap, PMIDBioMap, FileName };


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
  const DataOrganism data_organism;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Virtual file DB Base class, information about data files (populations) returned from the data parsers.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DataDB {

public:

  DataDB(DataSourceEnum data_source) : data_source_(data_source) {}
  virtual ~DataDB() = default;

  [[nodiscard]] virtual const std::string& fileId() const = 0;
  [[nodiscard]] DataSourceEnum dataSource() const { return data_source_; }
  [[nodiscard]] DataCharacteristic dataCharacteristic() const;

  [[nodiscard]] static std::optional<DataCharacteristic> findCharacteristic(const std::string& source_text);
  [[nodiscard]] static std::optional<DataCharacteristic> findCharacteristic(DataSourceEnum data_source);
  [[nodiscard]] static std::string dataSource(DataSourceEnum data_source);

private:

// The data file name.
  DataSourceEnum data_source_;

// DataCharacteristic vector.
  static const std::vector<DataCharacteristic> data_characteristics_;

};


// Data type characteristics vector.
// Additional data file types must be added here.
inline const std::vector<DataCharacteristic>  DataDB::data_characteristics_ = {

    { "Genome1000", DataSourceEnum::Genome1000, ParserTypeEnum::DiploidPhased, DataStructureEnum::DiploidPhased, DataImplEnum::PopulationVariant, DataOrganism::HomoSapien},
    { "GnomadGenome3_1", DataSourceEnum::GnomadGenome3_1, ParserTypeEnum::DiploidGnomad, DataStructureEnum::DiploidUnphased, DataImplEnum::PopulationVariant,         DataOrganism::HomoSapien },
    { "Falciparum", DataSourceEnum::Falciparum, ParserTypeEnum::DiploidFalciparum, DataStructureEnum::DiploidUnphased, DataImplEnum::PopulationVariant,               DataOrganism::PlasmodiumFalciparum },
    { "GnomadExomes3_1", DataSourceEnum::GnomadExomes3_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant, DataOrganism::HomoSapien },
    { "GnomadExomes2_1", DataSourceEnum::GnomadExomes2_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant, DataOrganism::HomoSapien },
    { "Gnomad3_1", DataSourceEnum::Gnomad3_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant,             DataOrganism::HomoSapien },
    { "Gnomad3_0", DataSourceEnum::Gnomad3_0, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant,             DataOrganism::HomoSapien },
    { "Gnomad2_1", DataSourceEnum::Gnomad2_1, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant,             DataOrganism::HomoSapien },
    { "Clinvar", DataSourceEnum::Clinvar, ParserTypeEnum::MonoGenomeUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant,                 DataOrganism::HomoSapien },
    { "dbSNP", DataSourceEnum::dbSNP, ParserTypeEnum::MonoDBSNPUnphased, DataStructureEnum::UnphasedMonoGenome, DataImplEnum::PopulationVariant,                      DataOrganism::HomoSapien },
    { "JSONdbSNP", DataSourceEnum::JSONdbSNP, ParserTypeEnum::MonoJSONdbSNPUnphased, DataStructureEnum::CitationMap, DataImplEnum::PMIDCitationMap,                   DataOrganism::HomoSapien },
    { "FileNameOnly", DataSourceEnum::NotImplemented, ParserTypeEnum::FilenameOnly, DataStructureEnum::NoStructure, DataImplEnum::FileName,                           DataOrganism::NoOrganism},
    { "BioPMID", DataSourceEnum::BioPMID, ParserTypeEnum::ParseBioPMID, DataStructureEnum::BioPMIDMap, DataImplEnum::PMIDBioMap,                                      DataOrganism::NoOrganism},

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Data file processing handled within the requesting package.
// Just pass the data source and file name back to the requesting package.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FilenameDataDB : public DataDB {

public:

  FilenameDataDB(DataSourceEnum data_source, const std::string& file_name) : DataDB(data_source) , file_name_(file_name) {}
  ~FilenameDataDB() override = default;

  [[nodiscard]] const std::string& fileId() const override { return file_name_; }


private:

  const std::string file_name_;

};







} // namespace


#endif //KGL_KGL_DATA_FILE_H
