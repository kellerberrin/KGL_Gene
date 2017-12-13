//
// Created by kellerberrin on 2/12/17.
//


#include <memory>
#include <fstream>
#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_statistics.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for per coding feature for statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This function will insert multiple variants for a contig offset in a std::multimap
// Only adds a variant if it does not already exist (not equivalent, but can have the same offset).
void kgl::FeatureStatistics::addVariant(const std::shared_ptr<const kgl::Variant>& variant_ptr) {

  auto result = offset_variant_map_.equal_range(variant_ptr->contigOffset());

  for (auto it = result.first; it != result.second; ++it) {

    if (it->second->equivalent(*variant_ptr)) return;

  }

  offset_variant_map_.insert(std::make_pair(variant_ptr->contigOffset(), variant_ptr));

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for per contig statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ContigStatistics::addVariant(std::shared_ptr<const Variant>& variant_ptr) {

  // if not coding then ignore for now

  if (variant_ptr->type() != VariantSequenceType::CDS_CODING) {

    return true;

  }

  std::shared_ptr<FeatureStatistics> feature_statistics;
  if (not getFeatureStatistics(variant_ptr->codingSequenceId(), feature_statistics)) {

    ExecEnv::log().error("Feature: {} not found, variant: {}",
                         variant_ptr->codingSequenceId(), variant_ptr->output(' ', VariantOutputIndex::START_0_BASED));
    return false;
  }

  feature_statistics->addVariant(variant_ptr);

  return true;

}


bool kgl::ContigStatistics::getFeatureStatistics(const FeatureIdent_t& feature_id,
                                                 std::shared_ptr<FeatureStatistics>& feature_statistics) const {

  bool result;

  auto feature = feature_map_.find(feature_id);

  if (feature != feature_map_.end()) {

    feature_statistics = feature->second;
    result = true;

  } else {

    feature_statistics = nullptr;
    result = false;

  }

  return result;

}


bool kgl::ContigStatistics::addFeatureStatistics(const FeatureIdent_t& feature_id,
                                                 const std::shared_ptr<const CodingSequence>& coding_sequence) {


  std::shared_ptr<FeatureStatistics> feature_statistics(std::make_shared<FeatureStatistics>(feature_id, coding_sequence));
  std::pair<FeatureIdent_t, std::shared_ptr<FeatureStatistics>> insert_pair(feature_id, feature_statistics);
  auto result = feature_map_.insert(insert_pair);

  if (not result.second) {

    ExecEnv::log().error("Probable duplicate, could add statistical feature: {}", feature_id);
    return false;

  }

  return true;

}


kgl::ContigSize_t kgl::ContigStatistics::size() const {

  ContigSize_t size_count = 0;
  for (auto feature : feature_map_) {

    size_count += feature.second->size();

  }

  return size_count;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds variants for genome statistical analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Creates an empty genome variant with the same contig structure as the genome database.
kgl::GenomeStatistics::GenomeStatistics(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                        const std::shared_ptr<const GenomeVariant>& genome_variant_ptr) {


  genome_id_ = genome_variant_ptr->genomeId();

  for (auto contig_db : genome_db_ptr->getMap()) {

    std::shared_ptr<ContigStatistics> contig_statistics(std::make_shared<ContigStatistics>(contig_db.first));
    if (addContigStatistics(contig_statistics)) {

      for (auto gene : contig_db.second->getGeneMap()) {

        const std::shared_ptr<const CodingSequenceArray> coding_array = kgl::GeneFeature::getCodingSequences(gene.second);

        for (auto coding_seq : coding_array->getMap()) {

          contig_statistics->addFeatureStatistics(coding_seq.first, coding_seq.second);

        } // for all sequences.

      } // for all genes

    } else {

      ExecEnv::log().error("GenomeStatistics(), could not add contig statistics: {}", contig_db.first);

    }

  } // for all contigs.

  // Add in all the variants to the newly constructed statistics object.

  std::vector<std::shared_ptr<const Variant>> variant_vector;
  genome_variant_ptr->getVariants(variant_vector);

  for (auto variant : variant_vector) {

    if (not addVariant(variant)) {

      ExecEnv::log().error("GenomeStatistics(), could not add variant: {}",
                           variant->output(' ', VariantOutputIndex::START_0_BASED));

    }

  } // for variants.

}


bool kgl::GenomeStatistics::addContigStatistics(std::shared_ptr<kgl::ContigStatistics>& contig_statistics) {

  auto result = genome_statistics_map_.insert(std::make_pair(contig_statistics->contigId(), contig_statistics));

  return result.second;

}


bool kgl::GenomeStatistics::getContigStatistics(const ContigId_t& contig_id,
                                                std::shared_ptr<ContigStatistics>& contig_statistics) const {
  bool result;

  auto contig = genome_statistics_map_.find(contig_id);

  if (contig != genome_statistics_map_.end()) {

    contig_statistics = contig->second;
    result = true;

  } else {

    contig_statistics = nullptr;
    result = false;

  }

  return result;

}


bool kgl::GenomeStatistics::addVariant(std::shared_ptr<const Variant> variant) {

  std::shared_ptr<ContigStatistics> contig_variant;
  if (not getContigStatistics(variant->contigId(), contig_variant)) {

    ExecEnv::log().error("Contig: {} not found, variant: {}",
                         variant->contigId(), variant->output(' ', VariantOutputIndex::START_0_BASED));
    return false;
  }

  contig_variant->addVariant(variant);

  return true;

}


std::string kgl::GenomeStatistics::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  for (const auto& contig_statistics : genome_statistics_map_) {

    for (const auto& feature : contig_statistics.second->getMap()) {

      std::vector<std::pair<std::string,std::string>> description_vec;
      feature.second->codingSequence()->getGene()->getAttributes().getAllAttributes(description_vec);

      ss << genomeId() << delimiter << contig_statistics.first;
      ss << delimiter << feature.first << delimiter;
      ss << feature.second->codingSequence()->codingNucleotides() << delimiter;
      ss << feature.second->size();
      for (const auto& description : description_vec) {

        if (description.first == DESCRIPTION_KEY_) {

          ss << delimiter << description.second;

        }

      }

      ss << '\n';

    }

  }

  return ss.str();

}

bool kgl::GenomeStatistics::outputCSV(const std::string& file_name, VariantOutputIndex output_index) const {

  // open the file.
  std::fstream out_file(file_name, std::fstream::out | std::fstream::app);
  if (!out_file) {

    ExecEnv::log().error("Cannot open output CSV file (--outCSVFile): {}", file_name);
    return false;

  }

  out_file << output(',', output_index);

  return out_file.good();

}


kgl::ContigSize_t kgl::GenomeStatistics::size() const {

  ContigSize_t size_count = 0;
  for (auto contig : genome_statistics_map_) {

    size_count += contig.second->size();

  }

  return size_count;

}


bool kgl::GenomeStatistics::isElement(std::shared_ptr<const Variant> variant_ptr) const {

  // get contig
  std::shared_ptr<ContigStatistics> contig_variant;
  if (not getContigStatistics(variant_ptr->contigId(), contig_variant)) {

    return false;

  }

  // check if not coding
  if (variant_ptr->type() != VariantSequenceType::CDS_CODING) {

    return false;

  }

  // get feature.
  std::shared_ptr<FeatureStatistics> feature_statistics;
  if (not contig_variant->getFeatureStatistics(variant_ptr->codingSequenceId(), feature_statistics)) {

    return false;

  }

  // check if variant is present in the feature.
  auto result = feature_statistics->getMap().equal_range(variant_ptr->contigOffset());

  for (auto it = result.first; it != result.second; ++it) {

    if (variant_ptr->equivalent(*(it->second))) return true;

  }

  return false;

}


void kgl::GenomeStatistics::getVariants(std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  variant_vector.clear();

  for (auto contig : getMap()) {

    for(auto feature : contig.second->getMap()) {

      for (auto contig : feature.second->getMap()) {

        variant_vector.push_back(contig.second);

      }

    }

  }

}




kgl::DistanceType_t kgl::GenomeStatistics::distance(std::shared_ptr<const GenomeStatistics> genome_stats_ptr) const {

  DistanceType_t distance = 0;
  std::vector<std::shared_ptr<const Variant>> variant_vector;

  // find variants in genome_stats_ptr but not in this genome.
  genome_stats_ptr->getVariants(variant_vector);

  for (auto variant : variant_vector) {

    if (not isElement(variant)) {

      distance++;

    }

  }

  // find variants in this genome but not in this genome_stats_ptr
  getVariants(variant_vector);

  for (auto variant : variant_vector) {

    if (not genome_stats_ptr->isElement(variant)) {

      distance++;

    }

  }

  return distance;

}


std::shared_ptr<const kgl::GenomeStatistics> kgl::GenomeStatistics::merge(std::shared_ptr<const GenomeStatistics> merge_genome) const {

  std::shared_ptr<GenomeStatistics> merged = merge_genome->deepCopy();

  std::vector<std::shared_ptr<const Variant>> variant_vector;
  getVariants(variant_vector);

  for (const auto& variant : variant_vector) {

    merged->addVariant(variant);

  }

  return merged;

}


std::shared_ptr<kgl::GenomeStatistics> kgl::GenomeStatistics::deepCopy() const {

  std::shared_ptr<GenomeStatistics> copy = std::shared_ptr<GenomeStatistics>(std::make_shared<GenomeStatistics>(genomeId()));

  // Copy the index structure.
  for (auto contig : getMap()) {

    std::shared_ptr<ContigStatistics> contig_stats(std::make_shared<ContigStatistics>(contig.second->contigId()));
    copy->addContigStatistics(contig_stats);

    for (auto feature : contig.second->getMap()) {

      contig_stats->addFeatureStatistics(feature.second->featureId(), feature.second->codingSequence());

    }

  }

  // Add in the variants.
  std::vector<std::shared_ptr<const Variant>> variant_vector;
  getVariants(variant_vector);

  for (const auto& variant : variant_vector) {

    copy->addVariant(variant);

  }

  return copy;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Holds population genome statistical objects. For Phylogenetic analysis.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::PopulationStatistics::addGenomeStatistics(std::shared_ptr<const kgl::GenomeStatistics> genome_stats_ptr) {

  std::pair<GenomeId_t, std::shared_ptr<const kgl::GenomeStatistics>> insert_pair(genome_stats_ptr->genomeId(), genome_stats_ptr);
  auto result = population_statistics_map_.insert(insert_pair);

  if (not result.second) {

    ExecEnv::log().error("Could not add genome: {} to population statistics object - probable duplicate",
                         genome_stats_ptr->genomeId());

  }

  return result.second;

}


void kgl::PopulationStatistics::initUPGMA(NodeVector<const GenomeStatistics>& node_vector) const {

  node_vector.clear();
  for (auto genome : population_statistics_map_) {

    std::shared_ptr<PhyloNode<const GenomeStatistics>>
    phylo_node_ptr(std::make_shared<PhyloNode<const GenomeStatistics>>(genome.second, genome.second->size()));

    node_vector.push_back(phylo_node_ptr);

  }

}
