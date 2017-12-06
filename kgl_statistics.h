//
// Created by kellerberrin on 2/12/17.
//

#ifndef KGL_STATISTICS_H
#define KGL_STATISTICS_H


#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_attributes.h"
#include "kgl_variant_db.h"




namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Auxillary container objects to hold variants for statistical and phylogenetic analysis.
// Similar to the variant_db.h objects.
// Implemented separately so that statistics can be added and deleted as required, without
// changing the underlying the variant_db.h objects.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for per coding feature for statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using FeatureStatiticsMap = std::multimap<ContigOffset_t , std::shared_ptr<const Variant>>;
class FeatureStatistics {

public:

  explicit FeatureStatistics(const FeatureIdent_t& feature_id,
                             const std::shared_ptr<const CodingSequence>& coding_sequence) : feature_id_(feature_id),
                                                                                             coding_sequence_(coding_sequence) {}
  FeatureStatistics(const FeatureStatistics&) = default;
  ~FeatureStatistics() = default;

  void addVariant(const std::shared_ptr<const Variant>& variant_ptr);

  const FeatureIdent_t& featureId() const { return feature_id_; }
  const FeatureStatiticsMap& getMap() const { return offset_variant_map_; }
  std::shared_ptr<const CodingSequence> codingSequence() const { return coding_sequence_; }

  size_t size() const { return offset_variant_map_.size(); }


private:

  FeatureIdent_t feature_id_;
  std::shared_ptr<const CodingSequence> coding_sequence_;
  FeatureStatiticsMap offset_variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for per contig statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ContigStatiticsMap = std::map<FeatureIdent_t, std::shared_ptr<FeatureStatistics>>;
class ContigStatistics {

public:

  explicit ContigStatistics(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  ContigStatistics(const ContigStatistics&) = default;
  ~ContigStatistics() = default;


  bool addVariant(std::shared_ptr<const Variant>& variant_ptr);
  bool addFeatureStatistics(const FeatureIdent_t& feature_id,
                            const std::shared_ptr<const CodingSequence>& coding_sequence);
  bool getFeatureStatistics(const FeatureIdent_t& feature_id,
                            std::shared_ptr<FeatureStatistics>& feature_statistics) const;

  const ContigId_t& contigId() const { return contig_id_; }
  const ContigStatiticsMap& getMap() const { return feature_map_; }


private:

  ContigId_t contig_id_;
  ContigStatiticsMap feature_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for genome statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenomeStatisticsMap = std::map<ContigId_t, std::shared_ptr<ContigStatistics>>;
class GenomeStatistics {

public:

  GenomeStatistics(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                   const std::shared_ptr<const GenomeVariant>& genome_variant_ptr);
  GenomeStatistics(const GenomeStatistics&) = default;
  ~GenomeStatistics() = default;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& contig_id) { genome_id_ = contig_id; }

  bool addContigStatistics(std::shared_ptr<ContigStatistics>& contig_statistics);
  bool getContigStatistics(const ContigId_t& contig_id, std::shared_ptr<ContigStatistics>& contig_statistics) const;

  bool addVariant(std::shared_ptr<const Variant> variant);

  const GenomeStatisticsMap& getMap() const { return genome_statistics_map_; }

  std::string output(char field_delimiter, VariantOutputIndex output_index) const;
  bool outputCSV(const std::string& file_name, VariantOutputIndex output_index) const;

private:

  GenomeId_t genome_id_;
  GenomeStatisticsMap genome_statistics_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for population statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using PopulationStatisticsMap = std::map<GenomeId_t, std::shared_ptr<const GenomeStatistics>>;
class PopulationStatistics {

public:

  explicit PopulationStatistics() = default;
  PopulationStatistics(const PopulationStatistics&) = default;
  ~PopulationStatistics() = default;

  bool addGenomeStatistics(std::shared_ptr<const GenomeStatistics> genome_statistics);

  const PopulationStatisticsMap& getMap() const { return population_statistics_map_; }

private:

  PopulationStatisticsMap population_statistics_map_;

};



}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_STATISTICS_H
