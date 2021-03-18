//
// Created by kellerberrin on 15/1/21.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_analysis_mutation_clinvar.h"
#include "kgl_variant_mutation.h"



#include <fstream>
#include <memory_resource>


namespace kgl = kellerberrin::genome;



// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                          const std::shared_ptr<const GenomePEDData>& ped_data)
{

  const GafRecordMap& ont_map = genome_ptr->geneOntology().getMap();
  ResortGaf symbolic_gaf; // re-sort by the symbolic reference.
  symbolic_gaf.sortBySymbolic(ont_map);
  ResortGaf gene_id_gaf; // re-sort by the gene id field.
  gene_id_gaf.sortByGeneId(ont_map);

  gene_vector_.clear();

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [offset, gene_ptr] : contig_ptr->getGeneMap()) {

      std::vector<std::string> name_vec;
      gene_ptr->getAttributes().getName(name_vec);
      std::string name;
      std::string gaf_id;
      std::set<std::string> GO_set;

      if (not name_vec.empty()) {

        name = name_vec.front();
        auto result = symbolic_gaf.getMap().find(name);
        if (result != symbolic_gaf.getMap().end()) {

          const auto& [symbolic, ontology_ptr] = *result;
          gaf_id = ontology_ptr->gene_id();
          for (const auto& GO_record : ontology_ptr->goRecords()) {

            GO_set.insert(GO_record.second.ontolotgy_id);

          }

        }

      }

      // If gaf_id_ is empty then try a lookup with the gene id.
      if (gaf_id.empty()) {

        auto result = gene_id_gaf.getMap().find(gene_ptr->id());
        if (result != gene_id_gaf.getMap().end()) {

          const auto& [gene_id, ontology_ptr] = *result;
          gaf_id = ontology_ptr->gene_id();
          for (const auto& GO_record : ontology_ptr->goRecords()) {

            GO_set.insert(GO_record.second.ontolotgy_id);

          }

        }

      }

      GeneCharacteristic gene_characteristic;
      gene_characteristic.geneDefinition(gene_ptr, genome_ptr->genomeId(), name, gaf_id, GO_set);
      GeneMutation mutation;
      mutation.gene_characteristic = gene_characteristic;
      mutation.clinvar.updateEthnicity().updatePopulations(ped_data);
      mutation.gene_variants.updateLofEthnicity().updatePopulations(ped_data);
      mutation.gene_variants.updateHighEthnicity().updatePopulations(ped_data);
      mutation.gene_variants.updateModerateEthnicity().updatePopulations(ped_data);
      gene_vector_.push_back(mutation);

    } // Gene.

  } // Contig.

  return true;

}


bool kgl::GenomeMutation::variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr,
                                          const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                          const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                          const std::shared_ptr<const GenomePEDData>& ped_data) {



  ThreadPool thread_pool(ThreadPool::hardwareThreads());
  // A vector for futures.
  std::vector<std::future<GeneMutation>> future_vector;
  // Index by Ensembl gene code.
  EnsemblIndexMap ensembl_index_map = ensemblIndex(unphased_population_ptr);
  // Queue a thread for each gene.
  for (auto& gene_mutation : gene_vector_) {

    std::future<GeneMutation> future = thread_pool.enqueueTask( &GenomeMutation::geneSpanAnalysis,
                                                                this,
                                                                population_ptr,
                                                                unphased_population_ptr,
                                                                clinvar_population_ptr,
                                                                ped_data,
                                                                ensembl_index_map,
                                                                gene_mutation);
    future_vector.push_back(std::move(future));

  } // for genes

  gene_vector_.clear();

  // Wait on completed threads
  for (auto& future : future_vector) {

    gene_vector_.push_back(future.get());

  }

  return true;

}


kgl::GeneMutation kgl::GenomeMutation::geneSpanAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                                                         const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                                         const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                                         const std::shared_ptr<const GenomePEDData>& ped_data,
                                                         const EnsemblIndexMap& ensembl_index_map,
                                                         GeneMutation gene_mutation) {


  std::shared_ptr<const ContigDB> clinvar_contig;
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

    auto contig_opt = genome_ptr->getContig(gene_mutation.gene_characteristic.contigId());
    if (contig_opt) {

      std::shared_ptr<const ContigDB> contig_ptr = contig_opt.value();

      if (not contig_ptr->getMap().empty()) {

        // Determine which method to determine gene membership of variants.
        std::shared_ptr<const ContigDB> span_variant_ptr;
        switch(gene_membership_) {

          case VariantGeneMembership::BY_SPAN:
            span_variant_ptr = getGeneSpan(contig_ptr, gene_mutation.gene_characteristic);
            break;

          case VariantGeneMembership::BY_EXON:
            span_variant_ptr = getGeneExon(contig_ptr, gene_mutation.gene_characteristic);
            break;

          default:
          case VariantGeneMembership::BY_ENSEMBL:
            span_variant_ptr = getGeneEnsembl( ensembl_index_map, gene_mutation.gene_characteristic);
            break;

        }

        if (not clinvar_contig) {

          clinvar_contig = AnalyzeClinvar::getClinvarContig(gene_mutation.gene_characteristic.contigId(), clinvar_population_ptr);
          clinvar_contig = AnalyzeClinvar::FilterPathogenic(clinvar_contig);

        } else {

          if (clinvar_contig->contigId() != gene_mutation.gene_characteristic.contigId()) {

            clinvar_contig = AnalyzeClinvar::getClinvarContig(gene_mutation.gene_characteristic.contigId(), clinvar_population_ptr);
            clinvar_contig = AnalyzeClinvar::FilterPathogenic(clinvar_contig);

          }

        }

        gene_mutation.clinvar.processClinvar(genome_id, span_variant_ptr, clinvar_contig, ped_data);
        gene_mutation.gene_variants.ProcessVariantStats(genome_id, span_variant_ptr, unphased_population_ptr, ped_data);

      } // contig not empty

    } // if contig

  } // for genome

  return gene_mutation;

}


// Gets variants over the whole gene span.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneSpan(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneCharacteristic& gene_char) {

  auto subset_contig = contig_ptr->subset(gene_char.geneBegin(), gene_char.geneEnd());

  // Remove duplicates.
  subset_contig->inSituFilter(UniquePhasedFilter());

  return subset_contig;

}

// Get variants only occurring within exons for all mRNA sequences.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneExon(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneCharacteristic& gene_char) {

  std::shared_ptr<ContigDB> gene_contig(std::make_shared<ContigDB>(contig_ptr->contigId()));
  std::shared_ptr<const CodingSequenceArray> coding_sequence_array = GeneFeature::getCodingSequences(gene_char.genePtr());

  for (auto const& [sequence_id, sequence_ptr] : coding_sequence_array->getMap()) {

    for (const auto& [cds_id, cds_ptr] : sequence_ptr->getSortedCDS()) {

      gene_contig->merge(contig_ptr->subset(cds_ptr->sequence().begin(), cds_ptr->sequence().end()));

    }

  }

  // Remove duplicates.
  gene_contig->inSituFilter(UniquePhasedFilter());

  return gene_contig;

}

// Get variants matching the ensembl.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneEnsembl( const EnsemblIndexMap& ensembl_index_map,
                                                                          const GeneCharacteristic& gene_char) {

  std::shared_ptr<ContigDB> gene_contig(std::make_shared<ContigDB>(gene_char.contigId()));

  auto lower_iterator = ensembl_index_map.lower_bound(gene_char.ensemblId());
  auto const upper_iterator = ensembl_index_map.upper_bound(gene_char.ensemblId());

  while(lower_iterator != upper_iterator) {

    auto const& [ensembl_id, variant_ptr] = *lower_iterator;
    if (not gene_contig->addVariant(variant_ptr)) {

      ExecEnv::log().error( "GenomeMutation::getGeneEnsembl, unable to add variant: {}",
                            variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

    }

    ++lower_iterator;

  }

  // Remove duplicates.
  gene_contig->inSituFilter(UniquePhasedFilter());

  return gene_contig;

}


// Index variants by the Ensembl gene code in the vep field.
kgl::EnsemblIndexMap kgl::GenomeMutation::ensemblIndex(const std::shared_ptr<const PopulationDB>& unphased_population_ptr) {

  // Local object performs the indexing.
  struct IndexMap {

    IndexMap() : pmr_memory_(pmr_byte_array_, sizeof(pmr_byte_array_)), unique_ident_(&pmr_memory_) {}
    ~IndexMap() = default;

    EnsemblIndexMap ensembl_indexed_variants_;
    std::byte pmr_byte_array_[4096];
    std::pmr::monotonic_buffer_resource pmr_memory_;
    std::pmr::set<std::string> unique_ident_;
    VepIndexVector field_index_;
    bool initialized_{false};

    bool ensemblIndex(const std::shared_ptr<const Variant>& variant) {

      if (not initialized_) {

        field_index_ = InfoEvidenceAnalysis::getVepIndexes(*variant, std::vector<std::string>{VEP_ENSEMBL_FIELD_});
        initialized_ = true;

      }

      auto field_vector = InfoEvidenceAnalysis::getVepData(*variant, field_index_);

      // Only unique gene idents.
      for (auto const& field : field_vector) {

        // Only 1 field in the map.
        if (not field.empty()) {

          const auto& [field_ident, field_value] = *field.begin();
          unique_ident_.insert(field_value);

        }

      }

      for (auto const& ident : unique_ident_) {

        if (not ident.empty()) {

          ensembl_indexed_variants_.emplace(ident, variant);

        }

      }

      unique_ident_.clear();

      return true;

    }

  }; // IndexMap object.

  IndexMap index_map;

  unphased_population_ptr->processAll(index_map, &IndexMap::ensemblIndex);

  return index_map.ensembl_indexed_variants_;

}
