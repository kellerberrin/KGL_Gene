//
// Created by kellerberrin on 27/6/20.
//

#include "kgl_variant_factory_1000_impl.h"



namespace kgl = kellerberrin::genome;


// Process VCF header information.
void kgl::Genome1000VCFImpl::processVCFHeader(const VCFHeaderInfo& header_info) {

  // Investigate header.
  VCFContigMap vcf_contig_map;
  VCFInfoRecordMap vcf_info_map;
  if (not VCFParseHeader::parseVcfHeader( header_info, vcf_contig_map, vcf_info_map)) {

    ExecEnv::log().error("PfVCFImpl::processVCFHeader, Problem parsing header information in VCF file. No variants processed.");

  }

  // Store all the available VCF info fields.
  evidence_factory_.availableInfoFields(vcf_info_map);

  if (VCFParseHeader::checkVCFReferenceContigs(vcf_contig_map, genome_db_ptr_)) {

    ExecEnv::log().info("VCF File and Reference Genome Contig Sizes All Match");

  } else {

    ExecEnv::log().info("VCF File and Reference Genome Contig Size/Name Mis-Match");

  }

}

void kgl::Genome1000VCFImpl::readParseVCFImpl(const std::string &vcf_file_name) {

  readVCFFile(vcf_file_name);  // parsing.

}

// This is multi-threaded code called from the reader defined above.
void kgl::Genome1000VCFImpl::ProcessVCFRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) {

  try {

    ParseRecord(std::move(vcf_record_ptr));

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("Genome1000VCFImpl::ProcessVCFRecord, Unexpected Exception: {} thrown record ignored", e.what());

  }

}


// This is multithreaded code called from the reader defined above.
void kgl::Genome1000VCFImpl::ParseRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) {

  // Parse the info fields into a map.
  // For performance reasons the info field is std::moved - don't reference again.
  auto mutable_info = const_cast<std::string&>(vcf_record_ptr->info);
  std::shared_ptr<const DataMemoryBlock> info_evidence_ptr = evidence_factory_.createVariantEvidence(std::move(mutable_info));  // Each vcf record.

  // Look at the filter field for "Pass"
  bool passed_filter = Utility::toupper(vcf_record_ptr->filter) == PASSED_FILTERS_;

  // Convert VCF contig to genome contig.
  std::string contig = contig_alias_map_.lookupAlias(vcf_record_ptr->contig_id);

  if (getGenomeNames().size() != vcf_record_ptr->genotypeInfos.size()) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), vcf_record_ptr->genotypeInfos.size());

  }

  std::vector<std::string> alt_vector = Utility::charTokenizer(vcf_record_ptr->alt, MULTIPLE_ALT_SEPARATOR_);

  if (alt_vector.empty()) {

    ExecEnv::log().error("Genome1000VCFImpl::ProcessVCFRecord, Zero sized alt vector, alt: {}", vcf_record_ptr->alt);

  }

  // For all genotypes.

  std::map<size_t, std::vector<GenomeId_t>> phase_A_map, phase_B_map;
  for (size_t genotype_count = 0; genotype_count< vcf_record_ptr->genotypeInfos.size(); ++genotype_count)
  {

    auto const& genotype = vcf_record_ptr->genotypeInfos[genotype_count];
    auto indices = alternateIndex(contig, genotype, alt_vector);

    if (indices.first != REFERENCE_VARIANT_INDEX_) {

      phase_A_map[indices.first-1].push_back(getGenomeNames()[genotype_count]);

    }

    if (indices.second != REFERENCE_VARIANT_INDEX_) {

      phase_B_map[indices.second-1].push_back(getGenomeNames()[genotype_count]);

    }

  }

  addVariants(phase_A_map,
              contig,
              VariantPhase::DIPLOID_PHASE_A,
              vcf_record_ptr->offset,
              passed_filter,
              info_evidence_ptr,
              vcf_record_ptr->ref,
              vcf_record_ptr->id,
              alt_vector,
              vcf_record_ptr->line_number);

  addVariants(phase_B_map,
              contig,
              VariantPhase::DIPLOID_PHASE_B,
              vcf_record_ptr->offset,
              passed_filter,
              info_evidence_ptr,
              vcf_record_ptr->ref,
              vcf_record_ptr->id,
              alt_vector,
              vcf_record_ptr->line_number);

  if (vcf_record_ptr->line_number % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("Processed :{} records, Virtual variants processed: {}", vcf_record_ptr->line_number, variant_count_);
    ExecEnv::log().info("Offset: {}, Actual variants: {}", vcf_record_ptr->offset, actual_variant_count_);

  }

}


std::pair<size_t, size_t> kgl::Genome1000VCFImpl::alternateIndex( const std::string& contig,
                                                                  const std::string& genotype,
                                                                  const std::vector<std::string>& alt_vector) const {

  // Trim any whitespace.
  auto trim_genotype = Utility::trimEndWhiteSpace(genotype);

  // Check for info.
  if (trim_genotype.empty()) {

    ExecEnv::log().warn("Genome1000VCFImpl::alternateIndex, No Phase Vector information");
    return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

  }

  // Size of the GT Block.
  size_t GT_size{0};
  auto GT_index = trim_genotype.find(GT_SEPARATOR_);
  if (GT_index != std::string::npos) {

    GT_size = GT_index;

  } else {

    GT_size = trim_genotype.size();

  }


  // The phenotype size or the first ":" if that exists.
  std::string_view unphased_view(trim_genotype.c_str(), GT_size);
  // Look for the "|" phase separator.
  std::vector<std::string_view> phase_vector = Utility::viewTokenizer(unphased_view, PHASE_MARKER_);

  size_t phase_A_alt{REFERENCE_VARIANT_INDEX_};
  size_t phase_B_alt{REFERENCE_VARIANT_INDEX_};


  try {

    if (phase_vector.size() == 1) {

      // If the phase vector size is 1 then we assume that we are processing the X/Y chromosomes of a male_.
      if (unphased_view != REFERENCE_VARIANT_INDICATOR_ and unphased_view != ALT_REFERENCE_VARIANT_INDICATOR_) {

        switch(contig_alias_map_.lookupType(contig)) {

          case ChromosomeType::ALLOSOME_X:
            phase_A_alt = std::stoul(std::string(unphased_view));
            phase_B_alt = REFERENCE_VARIANT_INDEX_;
            break;

          case ChromosomeType::ALLOSOME_Y:
            phase_A_alt = REFERENCE_VARIANT_INDEX_;
            phase_B_alt = std::stoul(std::string(unphased_view));
            break;

          default:
            ExecEnv::log().warn("Genome1000VCFImpl::alternateIndex, Expected Autosomal Chromosome: {} to have 2 phases", contig);
            break;

        }

      }

      if (phase_A_alt > alt_vector.size() or phase_B_alt > alt_vector.size()) {

        ExecEnv::log().warn("Genome1000VCFImpl::alternateIndex, phase A index: {}, phase B index: {}, exceed alternate vector size: {}, genotype: {}",
                            phase_A_alt, phase_B_alt, phase_vector.size(), genotype);
        return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

      }

      return {phase_A_alt, phase_B_alt};

    }

    // Size of the vector is 2
    if (phase_vector[0].find(ABSTRACT_ALT_BRACKET_) == std::string::npos) {

      if (phase_vector[0] != REFERENCE_VARIANT_INDICATOR_ and phase_vector[0] != ALT_REFERENCE_VARIANT_INDICATOR_) {

        phase_A_alt = std::stoul(std::string(phase_vector[0]));

      }

    } else {

      ++abstract_variant_count_;

    }

    if (phase_vector[1].find(ABSTRACT_ALT_BRACKET_) == std::string::npos) {

      if (phase_vector[1] != REFERENCE_VARIANT_INDICATOR_ and phase_vector[0] != ALT_REFERENCE_VARIANT_INDICATOR_) {

        phase_B_alt = std::stoul(std::string(phase_vector[1]));

      }

    } else {

      ++abstract_variant_count_;

    }

  }
  catch(...) {

    ExecEnv::log().warn("Genome1000VCFImpl::alternateIndex, Problem converting phase indexes to unsigned longs, phase text: {}", genotype);
    return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

  }

  if (phase_A_alt > alt_vector.size() or phase_B_alt > alt_vector.size()) {

    ExecEnv::log().warn("Genome1000VCFImpl::alternateIndex, phase A index: {}, phase B index: {}, exceed alternate vector size: {}, genotype: {}",
                        phase_A_alt, phase_B_alt, phase_vector.size(), genotype);
    return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

  }

  return {phase_A_alt, phase_B_alt};

}

void kgl::Genome1000VCFImpl::addVariants( const std::map<size_t, std::vector<GenomeId_t>>& phase_map,
                                          const ContigId_t& contig,
                                          VariantPhase phase,
                                          ContigOffset_t offset,
                                          bool passed_filters,
                                          const std::shared_ptr<const DataMemoryBlock>& info_evidence_ptr,
                                          const std::string& reference,
                                          const std::string& identifier,
                                          const std::vector<std::string>& alt_vector,
                                          size_t vcf_record_count) {

  for (const auto& [alt_allele, genome_vector] : phase_map) {

    // Setup the evidence object.
    VariantEvidence evidence(vcf_record_count,
                             diploid_population_ptr_->dataSource(),
                             passed_filters,
                             info_evidence_ptr,
                             nullptr,
                             alt_allele,
                             alt_vector.size());
    // Add the variant.
    StringDNA5 reference_str(reference);
    StringDNA5 alternate_str(alt_vector[alt_allele]);

    std::shared_ptr<const Variant> variant_ptr(std::make_shared<const Variant>( contig,
                                                                                offset,
                                                                                phase,
                                                                                identifier,
                                                                                std::move(reference_str),
                                                                                std::move(alternate_str),
                                                                                evidence));

    if (addThreadSafeVariant(variant_ptr, genome_vector)) {

      ++actual_variant_count_;
      variant_count_ += genome_vector.size();

    } else {

      ExecEnv::log().error("Genome1000VCFImpl::addVariants, problem adding: {} variants", genome_vector.size());

    }

  }

}


bool kgl::Genome1000VCFImpl::addThreadSafeVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                                   const std::vector<GenomeId_t>& genome_vector) const {

  // The population structure can be updated concurrently (embedded mutexes).
  return diploid_population_ptr_->addVariant(variant_ptr, genome_vector);

}

