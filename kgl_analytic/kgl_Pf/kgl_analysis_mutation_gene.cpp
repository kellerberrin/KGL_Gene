//
// Created by kellerberrin on 15/1/21.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_analysis_mutation_clinvar.h"
#include "kgl_variant_mutation.h"



#include <fstream>


namespace kgl = kellerberrin::genome;






// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_ptr)
{

  const GafRecordMap& ont_map = genome_ptr->geneOntology().getMap();
  ResortGaf symbolic_gaf; // re-sort by the symbolic reference.
  symbolic_gaf.sortBySymbolic(ont_map);
  ResortGaf gene_id_gaf; // re-sort by the gene id field.
  gene_id_gaf.sortByGeneId(ont_map);

  gene_vector_.clear();

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [offset, gene_ptr] : contig_ptr->getGeneMap()) {

      const std::shared_ptr<const CodingSequenceArray> sequence_array = GeneFeature::getCodingSequences(gene_ptr);

      std::vector<std::string> name_vec;
      gene_ptr->getAttributes().getName(name_vec);
      std::string name;
      std::string gaf_id;
      if (not name_vec.empty()) {

        name = name_vec.front();
        auto result = symbolic_gaf.getMap().find(name);
        if (result != symbolic_gaf.getMap().end()) {

          gaf_id = result->second->gene_id();

        }

      }

      // If gaf_id is empty then try a lookup with the gene id.
      if (gaf_id.empty()) {

        auto result = gene_id_gaf.getMap().find(gene_ptr->id());
        if (result != gene_id_gaf.getMap().end()) {

          gaf_id = result->second->gene_id();

        }

      }

      std::vector<std::string> description_vec;
      gene_ptr->getAttributes().getDescription(description_vec);
      std::string description = description_vec.empty() ? "" : description_vec.front();

      std::vector<std::string> gene_biotype_vec;
      gene_ptr->getAttributes().getGeneBioType(gene_biotype_vec);
      std::string biotype = gene_biotype_vec.empty() ? "" : gene_biotype_vec.front();

      GeneCharacteristic gene_characteristic;

      gene_characteristic.genome = genome_ptr->genomeId();
      gene_characteristic.contig  = contig_id;
      gene_characteristic.gene_id = gene_ptr->id();
      gene_characteristic.gene_name = name;
      gene_characteristic.description = description;
      gene_characteristic.biotype = biotype;
      gene_characteristic.valid_protein = ContigReference::verifyGene(gene_ptr);
      gene_characteristic.gaf_id = gaf_id;
      gene_characteristic.gene_begin = gene_ptr->sequence().begin();
      gene_characteristic.gene_end = gene_ptr->sequence().end();
      gene_characteristic.gene_span = gene_ptr->sequence().length();
      gene_characteristic.strand = gene_ptr->sequence().strandText();
      gene_characteristic.sequences = sequence_array->size();
      if (not sequence_array->getMap().empty()) {

        auto [seq_name, seq_ptr] = *(sequence_array->getMap().begin());
        gene_characteristic.seq_name = seq_name;
        gene_characteristic.nucleotides = seq_ptr->codingNucleotides();
        gene_characteristic.exons = seq_ptr->exons();

      }
      auto const& attribute_map = gene_ptr->getAttributes().getMap();
      gene_characteristic.attribute_size = attribute_map.size();

      std::string concat_attributes{"\""};
      for (auto const& [key, attrib] : attribute_map) {

        if ( key != attribute_map.begin()->first
            or attrib != attribute_map.begin()->second) {

          concat_attributes += CONCAT_TOKEN;

        }

        concat_attributes += key;
        concat_attributes += ":";
        concat_attributes += attrib;

      }

      concat_attributes += "\"";

      gene_characteristic.attributes = concat_attributes;


      GeneMutation mutation;
      mutation.gene_characteristic = gene_characteristic;
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

  // Queue a thread for each gene.
  for (auto& gene_mutation : gene_vector_) {

    std::future<GeneMutation> future = thread_pool.enqueueTask( &GenomeMutation::geneSpanAnalysis,
                                                                this,
                                                                population_ptr,
                                                                unphased_population_ptr,
                                                                clinvar_population_ptr,
                                                                ped_data,
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
                                                         GeneMutation gene_mutation) {


  std::map<std::string, size_t> variant_distribution;
  std::shared_ptr<const ContigDB> clinvar_contig;
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

    auto contig_opt = genome_ptr->getContig(gene_mutation.gene_characteristic.contig);
    if (contig_opt) {

      std::shared_ptr<const ContigDB> contig_ptr = contig_opt.value();

      if (not contig_ptr->getMap().empty()) {


        ++gene_mutation.genome_count;


        std::shared_ptr<const ContigDB> span_variant_ptr = getGeneSpan(contig_ptr, gene_mutation.gene_characteristic);

        if (not clinvar_contig) {

          clinvar_contig = AnalyzeClinvar::getClinvarContig(gene_mutation.gene_characteristic.contig, clinvar_population_ptr);
          clinvar_contig = AnalyzeClinvar::FilterPathogenic(clinvar_contig);

        } else {

          if (clinvar_contig->contigId() != gene_mutation.gene_characteristic.contig) {

            clinvar_contig = AnalyzeClinvar::getClinvarContig(gene_mutation.gene_characteristic.contig, clinvar_population_ptr);
            clinvar_contig = AnalyzeClinvar::FilterPathogenic(clinvar_contig);

          }

        }

        processClinvar(genome_id, span_variant_ptr, clinvar_contig, ped_data, gene_mutation.clinvar);

        // Get phased VEP info.
        VepInfo lof_het_hom = geneSpanVep( span_variant_ptr, unphased_population_ptr);

        if (lof_het_hom.female_lof > 0) {

          ++gene_mutation.female_lof;

        }

        if (lof_het_hom.male_lof > 0) {

          ++gene_mutation.male_lof;

        }

        // Loss of Function in both chromosomes, Mendelian genetics.
        size_t hom_lof{0};
        if (lof_het_hom.female_lof > 0 and lof_het_hom.male_lof > 0) {

          hom_lof = 1;

        }

        gene_mutation.hom_lof += hom_lof;

        if (lof_het_hom.female_high_effect > 0) {

          ++gene_mutation.female_high_effect;

        }

        if (lof_het_hom.male_high_effect > 0) {

          ++gene_mutation.male_high_effect;

        }

        // High Impact in both chromosomes, Mendelian genetics.
        if (lof_het_hom.male_high_effect > 0 and lof_het_hom.female_high_effect > 0) {

          ++gene_mutation.hom_high_effect;

        }



        size_t span_variant_count = span_variant_ptr->variantCount();

        if (span_variant_count > 0) {

          ++gene_mutation.genome_variant;

        }

        gene_mutation.variant_count += span_variant_count;

        gene_mutation.span_variant_count += span_variant_count;

        double indel{0.0}, transition{0.0}, transversion{0.0};

        for (auto const& [offset, offset_ptr] : span_variant_ptr->getMap()) {

          OffsetDBArray variant_array = offset_ptr->getVariantArray();

          for (auto const& variant_ptr : variant_array) {


            if (variant_ptr->phaseId() == VariantPhase::DIPLOID_PHASE_A) {

              ++gene_mutation.female_phase;

            } else {

              ++gene_mutation.male_phase;

            }

            switch (variant_ptr->variantType()) {

              case VariantType::INDEL:
                indel += 1.0;
                break;

              case VariantType::TRANSITION:
                transition += 1.0;
                break;

              case VariantType::TRANSVERSION:
                transversion += 1.0;
                break;

            }

            auto find_result = variant_distribution.find(variant_ptr->identifier());
            if (find_result == variant_distribution.end()) {

              auto[it, insert_result] = variant_distribution.try_emplace(variant_ptr->identifier(), 1);
              if (not insert_result) {

                ExecEnv::log().error("GenomeMutation::variantAnalysis; cannot insert (duplicate) variant with identifier: {}",
                                     variant_ptr->identifier());

              }

            } else {

              auto&[variant_ident, count] = *find_result;
              ++count;

            }

          } //for variant

        } //for offset

        double sum = indel + transition + transversion;
        if (sum >= 1.0) {
          gene_mutation.indel = indel / sum;
          gene_mutation.transition = transition / sum;
          gene_mutation.transversion = transversion / sum;

        } else {

          gene_mutation.indel = 0.0;
          gene_mutation.transition = 0.0;
          gene_mutation.transversion = 0.0;

        }

        gene_mutation.unique_variants = variant_distribution.size();

        for (auto const& [offset, offset_ptr] : span_variant_ptr->getMap()) {

          if (offset_ptr->getVariantArray().size() == 1) {

            ++gene_mutation.heterozygous;

          } else if (offset_ptr->getVariantArray().size() == 2) {

            auto offset_array = offset_ptr->getVariantArray();
            if (offset_array.front()->homozygous(*offset_array.back())) {

              ++gene_mutation.homozygous;

            }

          }

        }

      } // contig not empty

    } // if contig

  } // for genome

  return gene_mutation;

}


void kgl::GenomeMutation::processClinvar( const GenomeId_t& genome_id,
                                          const std::shared_ptr<const ContigDB>& subject_variants,
                                          const std::shared_ptr<const ContigDB>& clinvar_contig,
                                          const std::shared_ptr<const GenomePEDData>& ped_data,
                                          GeneClinvar& results) {


  auto subject_clinvar = AnalyzeClinvar::findClinvar(subject_variants, clinvar_contig);
  auto info_vector = AnalyzeClinvar::clinvarInfo(subject_clinvar);
  size_t clinvar_variants = subject_clinvar->variantCount();
  if (clinvar_variants > 0) {

    ++results.phase.phase_either_;

  }

  for (auto const& record : info_vector) {

    // Save any unique descriptions.
    results.clinvar_desc.insert(record.clndn);

  }

  auto male_variants = subject_variants->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_B));
  auto male_clinvar = AnalyzeClinvar::findClinvar(male_variants, clinvar_contig);
  size_t male_phase_variants = male_clinvar->variantCount();
  if (male_phase_variants > 0) {

    ++results.phase.phase_male_;

  }

  auto female_variants = subject_variants->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_A));
  auto female_clinvar = AnalyzeClinvar::findClinvar(female_variants, clinvar_contig);
  size_t female_phase_variants = female_clinvar->variantCount();
  if (female_phase_variants > 0) {

    ++results.phase.phase_female_;

  }

  if (male_phase_variants > 0 and female_phase_variants > 0) {

    ++results.phase.phase_hom_;

  }

  size_t count = clinvar_variants > 0 ? 1 : 0;
  results.results.pedAnalysis(genome_id, count, ped_data);

}


