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

void kgl::GenomeMutation::analysisType() {

  std::string member_text;
  switch(gene_membership_) {

    case VariantGeneMembership::BY_SPAN:
      member_text = "gene variant membership by gene span";
      break;

    case VariantGeneMembership::BY_EXON:
      member_text = "gene variant membership by defined exons";
      break;

    default:
    case VariantGeneMembership::BY_ENSEMBL:
      member_text = "gene variant membership by Ensembl gene code (vep field)";
      break;

  }

  ExecEnv::log().info("Gene variant analysis begins with: {}", member_text);

}

// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_ptr,
                                          const std::shared_ptr<const GenomePEDData>& ped_data,
                                          const std::shared_ptr<const OntologyDatabase>& ontology_db_ptr)
{

  // Only execute this function once.

  if (analysis_initialized_) {

    return  true;

  }

  analysis_initialized_ = true;
  const GafRecordMap& ont_map = genome_ptr->geneOntology().getMap();
  ResortGaf symbolic_gaf; // re-sort by the symbolic reference.
  symbolic_gaf.sortBySymbolic(ont_map);
  ResortGaf gene_id_gaf; // re-sort by the gene id field.
  gene_id_gaf.sortByGeneId(ont_map);
  const GeneSynonymVector synonym_vector = genome_ptr->geneOntology().getSynonymVector();
  ResortIds resort_ids;
  resort_ids.sortByHGNC(synonym_vector);
  ExecEnv::log().info("GenomeMutation::genomeAnalysis; HGNC sorted Gene Ids: {}", resort_ids.getMap().size());

  gene_vector_.clear();
  ExecEnv::log().info("Creating Ontology Cache ...");
  OntologyCache ontology_cache(ontology_db_ptr);
  ExecEnv::log().info("Ontology Cache Created");

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
      gene_characteristic.geneDefinition(gene_ptr, genome_ptr->genomeId(), name, gaf_id, GO_set, resort_ids.getMap());
      GeneMutation mutation;
      mutation.gene_characteristic = gene_characteristic;
      mutation.clinvar.updateEthnicity().updatePopulations(ped_data);
      mutation.gene_variants.updateLofEthnicity().updatePopulations(ped_data);
      mutation.gene_variants.updateHighEthnicity().updatePopulations(ped_data);
      mutation.gene_variants.updateModerateEthnicity().updatePopulations(ped_data);
      mutation.ontology.processOntologyStats(mutation.gene_characteristic, ontology_db_ptr, ontology_cache);
      gene_vector_.push_back(mutation);

    } // Gene.

  } // Contig.

  return true;

}


bool kgl::GenomeMutation::variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr,
                                          const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                          const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                          const std::shared_ptr<const GenomePEDData>& ped_data,
                                          const std::shared_ptr<const OntologyDatabase>& ontology_db_ptr) {

  ThreadPool thread_pool(ThreadPool::hardwareThreads());
  // A vector for futures.
  std::vector<std::future<GeneMutation>> future_vector;
  // Index by Ensembl gene code.
  std::shared_ptr<const EnsemblIndexMap> ensembl_index_map_ptr = ensemblIndex(unphased_population_ptr);
  ExecEnv::log().info("Unphased variants sorted by Ensembl Gene code: {}, Total Unphased Variants: {}",
                      ensembl_index_map_ptr->size(), unphased_population_ptr->variantCount());
  // Queue a thread for each gene.
  for (auto& gene_mutation : gene_vector_) {

    std::future<GeneMutation> future = thread_pool.enqueueTask( &GenomeMutation::geneSpanAnalysis,
                                                                this,
                                                                population_ptr,
                                                                unphased_population_ptr,
                                                                clinvar_population_ptr,
                                                                ped_data,
                                                                ontology_db_ptr,
                                                                ensembl_index_map_ptr,
                                                                gene_mutation);
    future_vector.push_back(std::move(future));

  } // for genes


  std::vector<GeneMutation> gene_vector;
  // Wait on completed threads
  gene_vector.reserve(future_vector.size());
  for (auto& future : future_vector) {

    gene_vector.emplace_back(future.get());

  }

  // todo: This logic is inefficient, the entire gene vector is copied for each VCF file (24 times). Re-design and Re-code.
  gene_vector_ = std::move(gene_vector);

  return true;

}


kgl::GeneMutation kgl::GenomeMutation::geneSpanAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                                                         const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                                         const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                                         const std::shared_ptr<const GenomePEDData>& ped_data,
                                                         const std::shared_ptr<const OntologyDatabase>& ontology_db_ptr,
                                                         const std::shared_ptr<const EnsemblIndexMap>& ensembl_index_map_ptr,
                                                         GeneMutation gene_mutation) {


  std::shared_ptr<const ContigDB> clinvar_contig;
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

    auto contig_opt = genome_ptr->getContig(gene_mutation.gene_characteristic.contigId());
    if (contig_opt) {

      std::shared_ptr<const ContigDB> contig_ptr = contig_opt.value();

      if (not contig_ptr->getMap().empty()) {

        // Which method gene membership of variants is determined.
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
            span_variant_ptr = getGeneEnsemblSpan( contig_ptr, *ensembl_index_map_ptr, gene_mutation.gene_characteristic);
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
//        gene_mutation.ontology.processOntologyStats(gene_mutation.gene_characteristic, ontology_db_ptr);

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
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneEnsemblSpan( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                              const EnsemblIndexMap& ensembl_index_map,
                                                                              const GeneCharacteristic& gene_char) {
  if (not gene_char.ensemblId().empty()) {

    return getGeneEnsembl(contig_ptr, ensembl_index_map, gene_char);

  } else {

    return getGeneSpan(contig_ptr, gene_char);

  }

}


// Get variants matching the ensembl.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneEnsembl( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                          const EnsemblIndexMap& ensembl_index_map,
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

  // Select all variants with the ensembl identifier.
  auto ensembl_variants = contig_ptr->findContig(gene_contig);

  // Remove any duplicates.
  ensembl_variants->inSituFilter(UniquePhasedFilter());

  return ensembl_variants;

}


// Index variants by the Ensembl gene code in the vep field.
std::shared_ptr<const kgl::EnsemblIndexMap> kgl::GenomeMutation::ensemblIndex(const std::shared_ptr<const PopulationDB>& unphased_population_ptr) {

  // Local object performs the indexing.
  class IndexMap {

  public:

    IndexMap() : ensembl_indexed_variants_ptr_(std::make_shared<EnsemblIndexMap>()),
                 pmr_memory_(pmr_byte_array_, sizeof(pmr_byte_array_)),
                 unique_ident_(&pmr_memory_) {}
    ~IndexMap() = default;

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
          if (not field_value.empty()) {

            unique_ident_.insert(field_value);

          }

        }

      }

      for (auto const& ident : unique_ident_) {

        if (not ident.empty()) {

          ensembl_indexed_variants_ptr_->emplace(ident, variant);

        }

      }

      unique_ident_.clear();

      return true;

    }

    [[nodiscard]] std::shared_ptr<const EnsemblIndexMap> ensemblIndexedMap() const { return ensembl_indexed_variants_ptr_; }

  private:

    std::shared_ptr<EnsemblIndexMap> ensembl_indexed_variants_ptr_;
    std::byte pmr_byte_array_[PMR_BUFFER_SIZE_];
    std::pmr::monotonic_buffer_resource pmr_memory_;
    std::pmr::set<std::string> unique_ident_;
    VepIndexVector field_index_;
    bool initialized_{false};


  }; // IndexMap object.

  IndexMap index_map;

  unphased_population_ptr->processAll(index_map, &IndexMap::ensemblIndex);

  return index_map.ensemblIndexedMap();

}
