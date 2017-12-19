//
// Created by kellerberrin on 2/12/17.
//

#ifndef KGL_STATISTICS_H
#define KGL_STATISTICS_H


#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include <iostream>
#include "kgl_attributes.h"
#include "kgl_variant_db.h"
#include "kgl_statistics_upgma.h"




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

  ContigSize_t size() const { return offset_variant_map_.size(); }

private:

  FeatureIdent_t feature_id_;
  std::shared_ptr<const CodingSequence> coding_sequence_;
  FeatureStatiticsMap offset_variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for per contig statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ContigFeatureMap = std::map<FeatureIdent_t, std::shared_ptr<FeatureStatistics>>;
using ContigVariantMap = std::multimap<ContigOffset_t , std::shared_ptr<const Variant>>;

class ContigStatistics {

public:

  explicit ContigStatistics(std::shared_ptr<const ContigFeatures> contig_ptr) : contig_ptr_(contig_ptr) {}
  ContigStatistics(const ContigStatistics&) = default;
  ~ContigStatistics() = default;


  bool addVariant(std::shared_ptr<const Variant>& variant_ptr);
  bool addFeatureStatistics(const FeatureIdent_t& feature_id,
                            const std::shared_ptr<const CodingSequence>& coding_sequence);
  bool getFeatureStatistics(const FeatureIdent_t& feature_id,
                            std::shared_ptr<FeatureStatistics>& feature_statistics) const;

  const ContigId_t& contigId() const { return contig_ptr_->contigId(); }
  std::shared_ptr<const ContigFeatures> contig() const { return contig_ptr_; }

  const ContigFeatureMap& getFeatureMap() const { return feature_map_; }
  const ContigVariantMap& getVariantMap() const { return variant_map_; }

  void getVariant(ContigOffset_t begin,
                  ContigSize_t size,
                  std::vector<std::shared_ptr<const Variant>>& variant_vector) const;
  void prime_5(FeatureIdent_t featureId, ContigSize_t size, std::vector<std::shared_ptr<const Variant>>& variant_vector) const;
  void prime_3(FeatureIdent_t featureId, ContigSize_t size, std::vector<std::shared_ptr<const Variant>>& variant_vector) const;
  void prime_5(ContigSize_t size, std::vector<std::shared_ptr<const Variant>>& variant_vector) const;
  void prime_3(ContigSize_t size, std::vector<std::shared_ptr<const Variant>>& variant_vector) const;

  ContigSize_t size() const { return variant_map_.size(); }

private:

  std::shared_ptr<const ContigFeatures> contig_ptr_;
  ContigFeatureMap feature_map_;
  ContigVariantMap variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for genome statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenomeStatisticsMap = std::map<ContigId_t, std::shared_ptr<ContigStatistics>>;
class GenomeStatistics {

public:

  explicit GenomeStatistics(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                            const std::shared_ptr<const GenomeVariant>& genome_variant_ptr);
  explicit GenomeStatistics(const GenomeId_t& genome_id) { genomeId(genome_id); }
  explicit GenomeStatistics(const GenomeStatistics& copy) = delete;
  ~GenomeStatistics() = default;

  // Use this instead of the copy constructor.
  std::shared_ptr<GenomeStatistics> deepCopy() const;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& genome_id) { genome_id_ = genome_id; }

  bool addContigStatistics(std::shared_ptr<ContigStatistics>& contig_statistics);
  bool getContigStatistics(const ContigId_t& contig_id, std::shared_ptr<ContigStatistics>& contig_statistics) const;

  bool addVariant(std::shared_ptr<const Variant> variant);

  const GenomeStatisticsMap& getMap() const { return genome_statistics_map_; }

  ContigSize_t size() const;

  bool isElement(std::shared_ptr<const Variant> variant) const;
  void getVariants(std::vector<std::shared_ptr<const Variant>>& variant_vector) const;

  void prime_5(ContigSize_t size, std::vector<std::shared_ptr<const Variant>>& variant_vector) const;
  void prime_3(ContigSize_t size, std::vector<std::shared_ptr<const Variant>>& variant_vector) const;

  static std::string outputFeatureHeader(char delimiter);
  std::string outputFeature(char field_delimiter, VariantOutputIndex output_index) const;
  bool outputFeatureCSV(const std::string &file_name, VariantOutputIndex output_index) const;


  // UPGMA Classification functions
  void write_node(std::ofstream& outfile) const { outfile << genomeId(); }
  DistanceType_t distance(std::shared_ptr<const GenomeStatistics> GenomeStatistics) const;


private:

  GenomeId_t genome_id_;
  GenomeStatisticsMap genome_statistics_map_;

  static constexpr const char* DESCRIPTION_KEY_{"DESCRIPTION"};
  static constexpr const ContigSize_t PRIME_5_SIZE_{1000};
  static constexpr const ContigSize_t PRIME_3_SIZE_{1000};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for population statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using PopulationStatisticsMap = std::map<GenomeId_t, std::shared_ptr<const GenomeStatistics>>;
class PopulationStatistics {

public:

  explicit PopulationStatistics(const std::string& population_id) : population_id_(population_id) {}
  PopulationStatistics(const PopulationStatistics&) = default;
  ~PopulationStatistics() = default;

  bool addGenomeStatistics(std::shared_ptr<const GenomeStatistics> genome_statistics);

  const PopulationStatisticsMap& getMap() const { return population_statistics_map_; }

  // Initialize the UPGMA analysis.
  std::shared_ptr<NodeVector<const GenomeStatistics>> initUPGMA() const;

private:

  PopulationStatisticsMap population_statistics_map_;
  std::string population_id_;

};



}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_STATISTICS_H
