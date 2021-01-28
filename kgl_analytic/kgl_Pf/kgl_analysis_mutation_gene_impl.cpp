//
// Created by kellerberrin on 27/1/21.
//


#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_mutation.h"

#include <fstream>


namespace kgl = kellerberrin::genome;



std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneContig( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                         const GeneMutation& gene_mutation) {

  std::shared_ptr<ContigDB> gene_contig(std::make_shared<ContigDB>(gene_mutation.seq_name));
  OffsetVariantMap variant_map;

  contig_ptr->getSortedVariants( VariantSequence::UNPHASED,
                                 gene_mutation.gene_begin,
                                 gene_mutation.gene_end,
                                 variant_map);

  for (auto const& [offset, variant_ptr] : variant_map) {

    if (not gene_contig->addVariant(variant_ptr)) {

      ExecEnv::log().error("GenomeMutation::getGeneContig; contig: {} cannot add variant: {}",
                           gene_contig->contigId(), variant_ptr->output(',',VariantOutputIndex::START_0_BASED, false));

    }

  }


  return gene_contig;

}

std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneSpan(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneMutation& gene_mutation) {

  return contig_ptr->subset(gene_mutation.gene_begin, gene_mutation.gene_end);

}


std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneExon(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneMutation& gene_mutation) {

  std::shared_ptr<ContigDB> gene_contig(std::make_shared<ContigDB>(gene_mutation.seq_name));

  const std::shared_ptr<const CodingSequenceArray> sequence_array = GeneFeature::getCodingSequences(gene_mutation.gene_ptr);

  if (not sequence_array->getMap().empty()) {

    auto [seq_name, seq_ptr] = *(sequence_array->getMap().begin());

    SortedCDS sorted_cds = seq_ptr->getSortedCDS();
    for (auto const& [offset, cds_ptr] : sorted_cds) {

      // Get the exon dimensions
      size_t begin = cds_ptr->sequence().begin();
      size_t end = cds_ptr->sequence().end();

      // Get the variants in the exon.
      auto contig_sub_set_ptr = contig_ptr->subset(begin, end);
      // Add exon variants to the gene_contig object.
      contig_sub_set_ptr->processAll(*gene_contig, &ContigDB::addVariant);

    }

  }

  return gene_contig;

}



bool kgl::GenomeMutation::pedAnalysis( GeneMutation& gene_mutation,
                                       const GenomeId_t& genome_id,
                                       size_t data_count,
                                       const std::shared_ptr<const GenomePEDData>& ped_data) {


  auto result = ped_data->getMap().find(genome_id);

  if (result == ped_data->getMap().end()) {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; Genome sample: {} does not have a PED record", genome_id);
    return false;

  }

  auto const& [sample_id, ped_record] = *result;

  if (data_count > 0) {

    if (ped_record.sexType() == PedSexType::MALE) {

      ++gene_mutation.male_value;

    } else {

      ++gene_mutation.female_value;

    }

  }

  std::string super_pop = ped_record.superPopulation();

  if (super_pop == "EAS") {

    gene_mutation.EAS += data_count;

  } else if (super_pop == "EUR") {

    gene_mutation.EUR += data_count;

  } else if (super_pop == "SAS") {

    gene_mutation.SAS += data_count;

  } else if (super_pop == "AFR") {

    gene_mutation.AFR += data_count;

  } else if (super_pop == "AMR") {

    gene_mutation.AMR += data_count;

  } else {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; Unknown super population: {}", super_pop);

  }

  return true;

}


kgl::VepInfo kgl::GenomeMutation::geneSpanVep( const std::shared_ptr<const ContigDB>& span_contig,
                                               const std::shared_ptr<const PopulationDB>& unphased_population_ptr) {

  VepInfo vep_info;

  if (unphased_population_ptr->getMap().size() != 1) {

    ExecEnv::log().error("GenomeMutation::geneSpanVep; expected unphased population to have 1 genome, size if: {}", unphased_population_ptr->getMap().size());
    return vep_info;

  }

  auto [genomne_id, genome_ptr] = *(unphased_population_ptr->getMap().begin());

  auto contig_opt = genome_ptr->getContig(span_contig->contigId());

  if (not contig_opt) {

    return vep_info;

  }

  auto unphased_contig = contig_opt.value();

  auto phase_A_variants = span_contig->filterVariants(PhaseFilter(VariantSequence::DIPLOID_PHASE_A));
  auto found_unphased_A = unphased_contig->findContig(phase_A_variants);


  auto phase_B_variants = span_contig->filterVariants(PhaseFilter(VariantSequence::DIPLOID_PHASE_B));
  auto found_unphased_B = unphased_contig->findContig(phase_B_variants);

  vep_info.female_lof = VepCount(found_unphased_A, LOF_VEP_FIELD, LOF_HC_VALUE);
  vep_info.male_lof = VepCount(found_unphased_B, LOF_VEP_FIELD, LOF_HC_VALUE);

  VepSubFieldValues vep_field(IMPACT_VEP_FIELD);

  vep_field.getContigValues(phase_A_variants);

  if (vep_field.getMap().size() > 0) ExecEnv::log().info("GenomeMutation::geneSpanVep; IMPACT entries: {}", vep_field.getMap().size());

  auto result = vep_field.getMap().find(IMPACT_HIGH_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.female_high_effect = count;

  }

  result = vep_field.getMap().find(IMPACT_MODERATE_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.female_moderate_effect = count;

  }

  vep_field.getContigValues(phase_B_variants);

  result = vep_field.getMap().find(IMPACT_HIGH_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.male_high_effect = count;

  }

  result = vep_field.getMap().find(IMPACT_MODERATE_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.male_moderate_effect = count;

  }

  return vep_info;

}


size_t kgl::GenomeMutation::VepCount( const std::shared_ptr<const ContigDB>& vep_contig,
                                      const std::string& vep_field_ident,
                                      const std::string& vep_field_value) {

  VepSubFieldValues vep_field(vep_field_ident);
  vep_field.getContigValues(vep_contig);

  auto result = vep_field.getMap().find(vep_field_value);
  if (result != vep_field.getMap().end()) {

    auto [field_value, count] = *result;
    return count;

  }

  return 0;

}

