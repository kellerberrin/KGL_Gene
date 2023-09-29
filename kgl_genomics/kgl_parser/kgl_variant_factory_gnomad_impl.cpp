//
// Created by kellerberrin on 20/12/20.
//

#include "kgl_variant_factory_gnomad_impl.h"



namespace kgl = kellerberrin::genome;


// Process VCF header information.
void kgl::GenomeGnomadVCFImpl::processVCFHeader(const VCFHeaderInfo& header_info) {

  // Investigate header.
  VCFContigMap vcf_contig_map;
  VCFInfoRecordMap vcf_info_map;
  if (not VCFParseHeader::parseVcfHeader( header_info, vcf_contig_map, vcf_info_map)) {

    ExecEnv::log().error("GenomeGnomadVCFImpl::processVCFHeader; Problem parsing header information in VCF file. No variants processed.");

  }

  // Store all the available VCF info fields.
  evidence_factory_.availableInfoFields(vcf_info_map);

  if (VCFParseHeader::checkVCFReferenceContigs(vcf_contig_map, genome_db_ptr_)) {

    ExecEnv::log().info("VCF File and Reference Genome Contig Sizes All Match");

  } else {

    ExecEnv::log().info("VCF File and Reference Genome Contig Size/Name Mis-Match");

  }

}

void kgl::GenomeGnomadVCFImpl::readParseVCFImpl(const std::string &vcf_file_name) {

  readVCFFile(vcf_file_name);  // parsing.

}

// This is multi-threaded code called from the reader defined above.
void kgl::GenomeGnomadVCFImpl::ProcessVCFRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) {

  try {

    ParseRecord(std::move(vcf_record_ptr));

  }
  catch(const std::exception& e) {

    ExecEnv::log().error("GenomeGnomadVCFImpl::processVCFRecord; Unexpected Exception: '{}' thrown, VCF record ignored", e.what());

  }

}


// This is multithreaded code called from the reader defined above.
void kgl::GenomeGnomadVCFImpl::ParseRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) {

  // Parse the info fields into a map.
  // For performance reasons the info field is std::moved - don't reference again.
  auto mutable_info = const_cast<std::string&>(vcf_record_ptr->info);
  std::shared_ptr<const DataMemoryBlock> info_evidence_ptr = evidence_factory_.createVariantEvidence(std::move(mutable_info));  // Each vcf record.

  // Look at the filter field for "Pass"
  bool passed_filter = Utility::toupper(vcf_record_ptr->filter) == PASSED_FILTERS_;

  // Convert VCF contig to genome contig_ref_ptr.
  std::string contig = contig_alias_map_.lookupAlias(vcf_record_ptr->contig_id);

  if (getGenomeNames().size() != vcf_record_ptr->genotypeInfos.size()) {

    ExecEnv::log().warn("Genome Name Size: {}, Genotype count: {}",
                        getGenomeNames().size(), vcf_record_ptr->genotypeInfos.size());

  }

  std::vector<std::string> alt_vector = Utility::charTokenizer(vcf_record_ptr->alt, MULTIPLE_ALT_SEPARATOR_);

  if (alt_vector.empty()) {

    ExecEnv::log().error("GenomeGnomadVCFImpl::processVCFHeader; Zero sized alt vector, alt: {}", vcf_record_ptr->alt);

  }

  // For all genotypes.

  std::map<size_t, std::vector<GenomeId_t>> phase_A_map, phase_B_map;
  for (size_t genotype_count = 0; genotype_count < vcf_record_ptr->genotypeInfos.size(); ++genotype_count)
  {

    const std::string& genotype = vcf_record_ptr->genotypeInfos[genotype_count];
    const std::string& genome = getGenomeNames()[genotype_count];

    auto const [A_index, B_index] = alternateIndex(genotype, alt_vector);

    if (A_index != REFERENCE_VARIANT_INDEX_) {

      phase_A_map[A_index-1].push_back(genome);

    }

    if (B_index != REFERENCE_VARIANT_INDEX_) {

      phase_B_map[B_index-1].push_back(genome);

    }

  }

  addVariants(phase_A_map,
              contig,
              VariantPhase::UNPHASED,
              vcf_record_ptr->offset,
              passed_filter,
              info_evidence_ptr,
              vcf_record_ptr->ref,
              vcf_record_ptr->id,
              alt_vector,
              vcf_record_ptr->line_number);
  addVariants(phase_B_map,
              contig,
              VariantPhase::UNPHASED,
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


std::pair<size_t, size_t> kgl::GenomeGnomadVCFImpl::alternateIndex(const std::string& genotype, const std::vector<std::string>& alt_vector) const {

  if (genotype.size() < MINIMUM_GENOTYPE_SIZE_) {

    ExecEnv::log().warn("GenomeGnomadVCFImpl::alternateIndex; Genotype: {} Size: {} < {} characters", genotype, genotype.size(), MINIMUM_GENOTYPE_SIZE_);
    return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

  }

  // Size of the GT Block.
  size_t GT_size{0};
  auto GT_index = genotype.find(GT_SEPARATOR_);
  if (GT_index != std::string::npos) {

    GT_size = GT_index;

  } else {

    GT_size = genotype.size();

  }

  // The first 3 characters or the first ":" if that exists.
  std::string_view unphased_view(genotype.c_str(), GT_size);
  // Look for the "/" phase separator.
  std::vector<std::string_view> phase_vector = Utility::viewTokenizer(unphased_view, PHASE_MARKER_);

  // Init to Ref allele (which is a no-op).
  size_t phase_A_alt{REFERENCE_VARIANT_INDEX_};
  size_t phase_B_alt{REFERENCE_VARIANT_INDEX_};

  if (phase_vector.size() == 2) {
  // Phase separator found.

    try {

      if (phase_vector[0] != REFERENCE_VARIANT_INDICATOR_) {

        phase_A_alt = std::stoul(std::string(phase_vector[0]));

      }

      if (phase_vector[1] != REFERENCE_VARIANT_INDICATOR_) {

        phase_B_alt = std::stoul(std::string(phase_vector[1]));

      }

    }
    catch(...) {

      ExecEnv::log().warn("GenomeGnomadVCFImpl::alternateIndex; Problem converting phase indexes to unsigned longs, phase A: {}, phase B: {} genotype: {}",
                          std::string(phase_vector[0]), std::string(phase_vector[1]), genotype);
      return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

    }

  } else {
  // Phase separator not found, assume an X or Y chromosome with single Genotype indicator for males.

    try {

      if (unphased_view != REFERENCE_VARIANT_INDICATOR_) {

        phase_A_alt = std::stoul(std::string(unphased_view));

      }

    }
    catch(...) {

      ExecEnv::log().warn("GenomeGnomadVCFImpl::alternateIndex; Problem converting phase index to unsigned long, phase A: {} (only), genotype: {}",
                          std::string(unphased_view), genotype);
      return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

    }

  }

  if (phase_A_alt > alt_vector.size() or phase_B_alt > alt_vector.size()) {

    ExecEnv::log().warn("GenomeGnomadVCFImpl::alternateIndex; phase A index: {}, phase B index: {}, exceed alternate vector size: {}, genotype: {}",
                        phase_A_alt, phase_B_alt, phase_vector.size(), genotype);
    return {REFERENCE_VARIANT_INDEX_, REFERENCE_VARIANT_INDEX_};

  }

  return {phase_A_alt, phase_B_alt};

}

void kgl::GenomeGnomadVCFImpl::addVariants( const std::map<size_t, std::vector<GenomeId_t>>& phase_map,
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
                             population_ptr_->dataSource(),
                             passed_filters,
                             info_evidence_ptr,
                             nullptr,
                             alt_allele,
                             alt_vector.size());

    // Add the variant.
    std::shared_ptr<const Variant> variant_ptr(std::make_shared<const Variant>( contig,
                                                                                offset,
                                                                                phase,
                                                                                identifier,
                                                                                DNA5SequenceLinear(StringDNA5(reference)),
                                                                                DNA5SequenceLinear(StringDNA5(alt_vector[alt_allele])),
                                                                                evidence));

    if (addThreadSafeVariant(variant_ptr, genome_vector)) {

      ++actual_variant_count_;
      variant_count_ += genome_vector.size();

    } else {

      ExecEnv::log().error("GenomeGnomadVCFImpl::addVariants; problem adding: {} variants", genome_vector.size());

    }

  }

}


bool kgl::GenomeGnomadVCFImpl::addThreadSafeVariant(const std::shared_ptr<const Variant>& variant_ptr,
                                                    const std::vector<GenomeId_t>& genome_vector) const {

  // The population structure can be updated concurrently (embedded mutexes).
  return population_ptr_->addVariant(variant_ptr, genome_vector);

}

