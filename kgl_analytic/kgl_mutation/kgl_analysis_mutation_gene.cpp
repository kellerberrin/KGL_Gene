//
// Created by kellerberrin on 15/1/21.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_analysis_mutation_clinvar.h"
#include "kgl_variant_mutation.h"
#include "kgl_variant_sort.h"



#include <fstream>
#include <memory_resource>


namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;

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
bool kgl::GenomeMutation::genomeAnalysis( const std::vector<std::string>& target_genes,
                                          const std::shared_ptr<const GenomeReference>& genome_ptr,
                                          const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                          const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                                          const std::shared_ptr<const EnsemblHGNCResource>& nomenclature_ptr)
{

  // Only execute this function once.

  if (analysis_initialized_) {

    return  true;

  }

  analysis_initialized_ = true;



  //  kol::PolicyEvidence evidence_policy(kol::GO::getEvidenceType(kol::GO::EvidenceType::EXPERIMENTAL));

  kol::PolicyEvidence evidence_policy;  // all evidence codes.

  const std::vector<kol::GO::EvidenceCode> automated_electronic_evidence{kol::GO::EvidenceCode::IEA};
  evidence_policy.excludeEvidence(automated_electronic_evidence); // exclude IEA

  std::shared_ptr<const kol::TermAnnotation> term_annotation_ptr(std::make_shared<const kol::TermAnnotation>(evidence_policy,
                                                                                                             genome_ptr->geneOntology().getGafRecordVector(),
                                                                                                             kol::AnnotationGeneName::SYMBOLIC_GENE_ID));
  const GeneSynonymVector synonym_vector = nomenclature_ptr->getGeneSynonym();
  ResortIds resort_ids;
  resort_ids.sortByHGNC(synonym_vector);
  ExecEnv::log().info("GenomeMutation::genomeAnalysis; HGNC sorted Gene Ids: {}", resort_ids.getMap().size());

  gene_vector_.clear();
  ExecEnv::log().info("Creating Ontology Cache ...");

  OntologyCache ontology_cache(target_genes, term_annotation_ptr, ontology_db_ptr->goGraph());

  ExecEnv::log().info("Ontology Cache Created");

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [offset, gene_ptr] : contig_ptr->getGeneMap()) {

      std::vector<std::string> name_vec;
      gene_ptr->getAttributes().getName(name_vec);
      std::string name;
      std::string gaf_id;

      if (not name_vec.empty()) {

        name = name_vec.front();
        auto [uniprot_id, symbolic_id] = term_annotation_ptr->getGeneIdentifiers(name);
        gaf_id = uniprot_id;

      }

      // If gaf_id_ is empty then try a lookup with the gene id.
      if (gaf_id.empty()) {

        auto [uniprot_id, symbolic_id] = term_annotation_ptr->getGeneIdentifiers(gene_ptr->id());
        gaf_id = uniprot_id;

      }


      GeneCharacteristic gene_characteristic;
      gene_characteristic.geneDefinition(gene_ptr, genome_ptr->genomeId(), name, gaf_id, resort_ids.getMap());
      GeneMutation mutation;
      mutation.gene_characteristic = gene_characteristic;
      mutation.clinvar.updateEthnicity().updatePopulations(genome_aux_data);
      mutation.gene_variants.updateLofEthnicity().updatePopulations(genome_aux_data);
      mutation.gene_variants.updateHighEthnicity().updatePopulations(genome_aux_data);
      mutation.gene_variants.updateModerateEthnicity().updatePopulations(genome_aux_data);
      mutation.ontology.processOntologyStats(mutation.gene_characteristic.geneId(), ontology_cache);
      gene_vector_.push_back(mutation);

    } // Gene.

  } // Contig.

  return true;

}


bool kgl::GenomeMutation::variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr,
                                          const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                          const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                          const std::shared_ptr<const HsGenomeAux>& genome_aux_data) {

  // Count the ethnic samples in the populations.
  ethnic_statistics_.updatePopulations(genome_aux_data);
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    ethnic_statistics_.pedAnalysis(genome_id, 1, genome_aux_data);

  }
  if (not ethnic_statistics_.auditTotals()) {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; problem with ethnic statistics totals");

  }

  ThreadPool thread_pool(ThreadPool::hardwareThreads());
  // A vector for futures.
  std::vector<std::future<GeneMutation>> future_vector;
  // Index by Ensembl gene code.
  std::shared_ptr<const EnsemblIndexMap> ensembl_index_map_ptr = VariantSort::ensemblIndex(unphased_population_ptr);
  ExecEnv::log().info("Unphased variants sorted by Ensembl Gene code: {}, Total Unphased Variants: {}",
                      ensembl_index_map_ptr->size(), unphased_population_ptr->variantCount());
  // Queue a thread for each gene.
  for (auto& gene_mutation : gene_vector_) {

    std::future<GeneMutation> future = thread_pool.enqueueTask(&GenomeMutation::geneSpanAnalysis,
                                                               this,
                                                               population_ptr,
                                                               unphased_population_ptr,
                                                               clinvar_population_ptr,
                                                               genome_aux_data,
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
                                                         const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                                         const std::shared_ptr<const EnsemblIndexMap>& ensembl_index_map_ptr,
                                                         GeneMutation gene_mutation) {



  std::shared_ptr<const ContigDB> clinvar_contig;
  bool contig_data{false};
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

    auto contig_opt = genome_ptr->getContig(gene_mutation.gene_characteristic.contigId());
    if (contig_opt) {

      contig_data = true;
      std::shared_ptr<const ContigDB> contig_ptr = contig_opt.value();

      if (not contig_ptr->getMap().empty()) {

        // Which method gene membership of variants is determined.
        std::shared_ptr<const ContigDB> gene_variant_ptr;
        switch(gene_membership_) {

          case VariantGeneMembership::BY_SPAN:
            gene_variant_ptr = getGeneSpan(contig_ptr, gene_mutation.gene_characteristic);
            break;

          case VariantGeneMembership::BY_EXON:
            gene_variant_ptr = getGeneExon(contig_ptr, gene_mutation.gene_characteristic);
            break;

          default:
          case VariantGeneMembership::BY_ENSEMBL:
            gene_variant_ptr = getGeneEnsemblSpan(contig_ptr, *ensembl_index_map_ptr, gene_mutation.gene_characteristic);
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

        gene_mutation.clinvar.processClinvar(genome_id, gene_variant_ptr, clinvar_contig, genome_aux_data);
        gene_mutation.gene_variants.processVariantStats(genome_id, gene_variant_ptr, unphased_population_ptr, genome_aux_data);

      } // contig not empty

    } // if contig

  } // for genome

  // Generate aggregate variant statistics for each gene.
  if (contig_data) {

    if (not gene_mutation.gene_variants.processSummaryStatistics( population_ptr,
                                                                  ethnic_statistics_,
                                                                  gene_mutation.gene_characteristic.geneId())) {

      ExecEnv::log().warn("GenomeMutation::geneSpanAnalysis; problem with processSummaryStatistics, gene: {}",
                          gene_mutation.gene_characteristic.geneId());

    }

  } else {


    gene_mutation.gene_variants.initializeSummaryStatistics( ethnic_statistics_);

  }

  return gene_mutation;

}


// Gets variants over the whole gene span.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneSpan(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneCharacteristic& gene_char) {

  auto subset_contig = contig_ptr->subset(gene_char.geneBegin(), gene_char.geneEnd());

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

  return gene_contig;

}

// Get variants matching the ensembl.
std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneEnsemblSpan( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                              const EnsemblIndexMap& ensembl_index_map,
                                                                              const GeneCharacteristic& gene_char) {
  if (not gene_char.ensemblId().empty()) {

    return getGeneEnsembl(contig_ptr, ensembl_index_map, gene_char);

  } else {

    return getGeneExon(contig_ptr, gene_char);

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


  return ensembl_variants;

}


