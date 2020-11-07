//
// Created by kellerberrin on 27/6/20.
//

#include "kgl_variant_factory_1000_impl.h"



namespace kgl = kellerberrin::genome;


// Process VCF header information.
void kgl::Genome1000VCFImpl::processVCFHeader(const VcfHeaderInfo& header_info) {

  // Investigate header.
  VCFContigMap vcf_contig_map;
  VCFInfoRecordMap vcf_info_map;
  if (not VCFParseHeader::parseVcfHeader( header_info, vcf_contig_map, vcf_info_map)) {

    ExecEnv::log().error("Pf3kVCFImpl::processVCFHeader, Problem parsing header information in VCF file. No variants processed.");

  }

  // Store all the available VCF info fields.
  evidence_factory_.availableInfoFields(vcf_info_map);

  if (VCFParseHeader::checkVCFReferenceContigs(vcf_contig_map, genome_db_ptr_)) {

    ExecEnv::log().info("VCF File and Reference Genome Contig Sizes All Match");

  } else {

    ExecEnv::log().info("VCF File and Reference Genome Contig Size/Name Mis-Match");

  }

}

void kgl::Genome1000VCFImpl::readParseVCFImpl() {

  readVCFFile();  // parsing.

}

// This is multi-threaded code called from the reader defined above.
void kgl::Genome1000VCFImpl::ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) {

  try {

    ParseRecord(vcf_record_count, vcf_record);

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("Genome1000VCFImpl::ProcessVCFRecord, Unexpected Exception: {} thrown record ignored", e.what());

  }

}


// This is multithreaded code called from the reader defined above.
void kgl::Genome1000VCFImpl::ParseRecord(size_t vcf_record_count, const VcfRecord& record) {

  // Parse the info fields into a map.
  // For performance reasons the info field is std::moved - don't reference again.
  auto mutable_info = const_cast<std::string&>(record.info);
  InfoDataEvidence info_evidence_opt = evidence_factory_.createVariantEvidence(std::move(mutable_info));  // Each vcf record.

  // Look at the filter field for "Pass"
  bool passed_filter = Utility::toupper(record.filter) == PASSED_FILTERS_;

  // Convert VCF contig to genome contig.
  std::string contig = contig_alias_map_.lookupAlias(record.contig_id);

  if (getGenomeNames().size() != record.genotypeInfos.size()) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), record.genotypeInfos.size());

  }

  std::vector<std::string> alt_vector = Utility::char_tokenizer(record.alt, MULTIPLE_ALT_SEPARATOR_);

  if (alt_vector.empty()) {

    ExecEnv::log().error("Genome1000VCFImpl::ProcessVCFRecord, Zero sized alt vector, alt: {}", record.alt);

  }

  // For all genotypes.
  size_t genotype_count = 0;
  std::map<size_t, std::vector<GenomeId_t>> phase_A_map, phase_B_map;
  for (auto const& genotype : record.genotypeInfos)
  {

    auto indices = alternateIndex(genotype, alt_vector);

    if (indices.first != REFERENCE_VARIANT_INDEX_) {

      phase_A_map[indices.first-1].push_back(getGenomeNames()[genotype_count]);

    }

    if (indices.second != REFERENCE_VARIANT_INDEX_) {

      phase_B_map[indices.second-1].push_back(getGenomeNames()[genotype_count]);

    }

    // Next genome name.
    ++genotype_count;

  }

  addVariants(phase_A_map, contig, VariantSequence::DIPLOID_PHASE_A, record.offset, passed_filter, info_evidence_opt, record.ref, alt_vector, vcf_record_count);
  addVariants(phase_B_map, contig, VariantSequence::DIPLOID_PHASE_B, record.offset, passed_filter, info_evidence_opt, record.ref, alt_vector, vcf_record_count);

  if (vcf_record_count % VARIANT_REPORT_INTERVAL_ == 0) {

    ExecEnv::log().info("Processed :{} records, Virtual variants processed: {}", vcf_record_count, variant_count_);
    ExecEnv::log().info("Offset: {}, Actual variants: {}", record.offset, actual_variant_count_);

  }

}


std::pair<size_t, size_t> kgl::Genome1000VCFImpl::alternateIndex(const std::string& genotype, const std::vector<std::string>& alt_vector) const {

  std::vector<std::string> phase_vector = Utility::char_tokenizer(genotype, PHASE_MARKER_);

  if (phase_vector.size() != 2) {

    ExecEnv::log().warn("Genome1000VCFImpl::alternateIndex, Unexpected Phase Vector Size: {} Genotype Info {}", phase_vector.size(), genotype);
    return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

  }

  size_t phase_A_alt{REFERENCE_VARIANT_INDEX_};
  size_t phase_B_alt{REFERENCE_VARIANT_INDEX_};

  try {

    if (phase_vector[0].find(ABSTRACT_ALT_BRACKET_) == std::string::npos) {

      if (phase_vector[0] != REFERENCE_VARIANT_INDICATOR_) {

        phase_A_alt = std::stoul(phase_vector[0]);

      }

    } else {

      ++abstract_variant_count_;

    }

    if (phase_vector[1].find(ABSTRACT_ALT_BRACKET_) == std::string::npos) {

      if (phase_vector[1] != REFERENCE_VARIANT_INDICATOR_) {

        phase_B_alt = std::stoul(phase_vector[1]);

      }

    } else {

      ++abstract_variant_count_;

    }

  }
  catch(...) {

    ExecEnv::log().warn("Genome1000VCFImpl::ParseRecord, Problem converting phase indexes to unsigned longs, phase A: {}, phase B: {} genotype: {}",
                          phase_vector[0], phase_vector[1], genotype);
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
                                          PhaseId_t phase,
                                          ContigOffset_t offset,
                                          bool passedFilters,
                                          const InfoDataEvidence info_evidence_opt,
                                          const std::string& reference,
                                          const std::vector<std::string>& alt_vector,
                                          size_t vcf_record_count) {

  for (const auto& [alt_allele, genome_vector] : phase_map) {

    std::optional<std::shared_ptr<FormatData>> null_format_data = std::nullopt;
// Setup the evidence object.
    VariantEvidence evidence(vcf_record_count, info_evidence_opt, null_format_data, alt_allele, alt_vector.size());
    // Add the variant.
    StringDNA5 reference_str(reference);
    StringDNA5 alternate_str(alt_vector[alt_allele]);

    std::unique_ptr<const Variant> variant_ptr(std::make_unique<Variant>( contig,
                                                                          phase,
                                                                          offset,
                                                                          passedFilters,
                                                                          evidence,
                                                                          std::move(reference_str),
                                                                          std::move(alternate_str)));


    if (addThreadSafeVariant(std::move(variant_ptr), genome_vector)) {

      ++actual_variant_count_;
      variant_count_ += genome_vector.size();

    } else {

      ExecEnv::log().error("Genome1000VCFImpl::addVariants, problem adding: {} variants", genome_vector.size());

    }

  }

}


bool kgl::Genome1000VCFImpl::addThreadSafeVariant( std::unique_ptr<const Variant>&& variant_ptr,
                                                   const std::vector<GenomeId_t>& genome_vector) const {

  // The population structure can be updated concurrently (embedded mutexes).
  return diploid_population_ptr_->addVariant(std::move(variant_ptr), genome_vector);

}

