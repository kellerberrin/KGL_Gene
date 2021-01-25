//
// Created by kellerberrin on 15/1/21.
//

#include <kgl_variant_factory_vcf_evidence_analysis.h>
#include "kgl_analysis_mutation_gene.h"
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


      GeneMutation gene_mutation;

      gene_mutation.gene_ptr = gene_ptr;
      gene_mutation.genome = genome_ptr->genomeId();
      gene_mutation.contig  = contig_id;
      gene_mutation.gene_id = gene_ptr->id();
      gene_mutation.gene_name = name;
      gene_mutation.description = description;
      gene_mutation.biotype = biotype;
      gene_mutation.valid_protein = ContigReference::verifyGene(gene_ptr);
      gene_mutation.gaf_id = gaf_id;
      gene_mutation.gene_begin = gene_ptr->sequence().begin();
      gene_mutation.gene_end = gene_ptr->sequence().end();
      gene_mutation.gene_span = gene_ptr->sequence().length();
      gene_mutation.strand = gene_ptr->sequence().strandText();
      gene_mutation.sequences = sequence_array->size();
      if (not sequence_array->getMap().empty()) {

        auto [seq_name, seq_ptr] = *(sequence_array->getMap().begin());
        gene_mutation.seq_name = seq_name;
        gene_mutation.nucleotides = seq_ptr->codingNucleotides();
        gene_mutation.exons = seq_ptr->exons();

      }
      auto const& attribute_map = gene_ptr->getAttributes().getMap();
      gene_mutation.attribute_size = attribute_map.size();

      gene_vector_.push_back(gene_mutation);

    } // Gene.

  } // Contig.

  return true;

}


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


bool kgl::GenomeMutation::variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr,
                                          const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
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
                                                         const std::shared_ptr<const GenomePEDData>& ped_data,
                                                         GeneMutation gene_mutation) {

  std::map<std::string, size_t> variant_distribution;
  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

    auto contig_opt = genome_ptr->getContig(gene_mutation.contig);
    if (contig_opt) {

      std::shared_ptr<const ContigDB> contig_ptr = contig_opt.value();

      if (not contig_ptr->getMap().empty()) {

        ++gene_mutation.genome_count;

        std::shared_ptr<const ContigDB> exome_variant_ptr = getGeneExon(contig_ptr, gene_mutation);

        size_t variant_count = exome_variant_ptr->variantCount();

        std::shared_ptr<const ContigDB> span_variant_ptr = getGeneSpan(contig_ptr, gene_mutation);

        // .first is the count of phase A lof, .second is the count of phase B lof.
        std::pair<size_t, size_t> lof_het_hom = geneSpanVep( span_variant_ptr, unphased_population_ptr);

        if (lof_het_hom.first > 0 and lof_het_hom.second > 0) {

          ++gene_mutation.hom_lof;

        } else if ((lof_het_hom.first + lof_het_hom.second) > 0) {

          ++gene_mutation.het_lof;

        }

        size_t span_variant_count = span_variant_ptr->variantCount();

        if (variant_count > 0) {

          ++gene_mutation.genome_variant;

        }

        gene_mutation.variant_count += variant_count;

        gene_mutation.span_variant_count += span_variant_count;

        if (ped_data) {

          pedAnalysis( gene_mutation, genome_id, variant_count, ped_data);

        }

        double indel{0.0}, transition{0.0}, transversion{0.0};

        for (auto const& [offset, offset_ptr] : exome_variant_ptr->getMap()) {

          OffsetDBArray variant_array = offset_ptr->getVariantArray();

          for (auto const& variant_ptr : variant_array) {


            if (variant_ptr->phaseId() == VariantSequence::DIPLOID_PHASE_A) {

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

        for (auto const& [offset, offset_ptr] : exome_variant_ptr->getMap()) {

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



bool kgl::GenomeMutation::pedAnalysis( GeneMutation& gene_mutation,
                                       const GenomeId_t& genome_id,
                                       size_t variant_count,
                                       const std::shared_ptr<const GenomePEDData>& ped_data) {


  auto result = ped_data->getMap().find(genome_id);

  if (result == ped_data->getMap().end()) {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; Genome sample: {} does not have a PED record", genome_id);
    return false;

  }

  auto const& [sample_id, ped_record] = *result;

  if (variant_count > 0) {

    if (ped_record.sexType() == PedSexType::MALE) {

      ++gene_mutation.male_variant;

    } else {

      ++gene_mutation.female_variant;

    }

  }

  std::string super_pop = ped_record.superPopulation();

  if (super_pop == "EAS") {

    gene_mutation.EAS_variant_count += variant_count;

  } else if (super_pop == "EUR") {

    gene_mutation.EUR_variant_count += variant_count;

  } else if (super_pop == "SAS") {

    gene_mutation.SAS_variant_count += variant_count;

  } else if (super_pop == "AFR") {

    gene_mutation.AFR_variant_count += variant_count;

  } else if (super_pop == "AMR") {

    gene_mutation.AMR_variant_count += variant_count;

  } else {

    ExecEnv::log().error("GenomeMutation::variantAnalysis; Unknown super population: {}", super_pop);

  }

  return true;

}


std::pair<size_t, size_t> kgl::GenomeMutation::geneSpanVep( const std::shared_ptr<const ContigDB>& span_contig,
                                                            const std::shared_ptr<const PopulationDB>& unphased_population_ptr) {

  if (unphased_population_ptr->getMap().size() != 1) {

    ExecEnv::log().error("GenomeMutation::geneSpanVep; expected unphased population to have 1 genome, size if: {}", unphased_population_ptr->getMap().size());
    return {0, 0};

  }

  auto [genomne_id, genome_ptr] = *(unphased_population_ptr->getMap().begin());

  auto contig_opt = genome_ptr->getContig(span_contig->contigId());

  if (not contig_opt) {

    return {0, 0};

  }

  auto unphased_contig = contig_opt.value();

  auto phase_A_variants = span_contig->filterVariants(PhaseFilter(VariantSequence::DIPLOID_PHASE_A));
  auto found_unphased_A = unphased_contig->findContig(phase_A_variants);


  auto phase_B_variants = span_contig->filterVariants(PhaseFilter(VariantSequence::DIPLOID_PHASE_B));
  auto found_unphased_B = unphased_contig->findContig(phase_B_variants);

  return { VepCount(found_unphased_A), VepCount(found_unphased_B) };

}


size_t kgl::GenomeMutation::VepCount( const std::shared_ptr<const ContigDB>& vep_contig) {

  VepSubFieldValues vep_field(LOF_VEP_FIELD);
  vep_field.getContigValues(vep_contig);

  auto result = vep_field.getMap().find(LOF_HC_FIELD_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_value, count] = *result;
    return count;

  }

  return 0;

}


// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::writeOutput(const std::string& output_file_name, char output_delimiter) const {

  const double homozygous_bias{0.1};
  const double sex_bias{0.1};

  std::ofstream out_file(output_file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("GenomeMutation::writeOutput; could not open file: {} for output", output_file_name);
    return false;

  } else {

    ExecEnv::log().info("GenomeMutation writing output to file: {}", output_file_name);

  }

  writeHeader(out_file, output_delimiter);

  for (auto const& gene : gene_vector_) {
    
    out_file << gene.genome << output_delimiter
             << gene.contig << output_delimiter
             << gene.gene_id << output_delimiter
             << gene.gene_name << output_delimiter
             << gene.description << output_delimiter
             << gene.biotype << output_delimiter
             << (gene.valid_protein ? "Valid" : "Invalid") << output_delimiter
             << gene.gaf_id << output_delimiter
             << gene.gene_begin << output_delimiter
             << gene.gene_end << output_delimiter
             << gene.gene_span << output_delimiter
             << gene.strand << output_delimiter
             << gene.sequences << output_delimiter
             << gene.seq_name << output_delimiter
             << gene.nucleotides << output_delimiter
             << gene.exons << output_delimiter
             << gene.attribute_size << output_delimiter
             << gene.unique_variants << output_delimiter
             << gene.variant_count << output_delimiter
             << gene.span_variant_count << output_delimiter
             << gene.het_lof << output_delimiter
             << gene.hom_lof << output_delimiter
             << (static_cast<double>(gene.female_phase) / static_cast<double>(gene.male_phase + sex_bias)) << output_delimiter
             << (static_cast<double>(gene.female_variant) / static_cast<double>(gene.male_variant + sex_bias)) << output_delimiter
             << gene.genome_count << output_delimiter
             << gene.genome_variant << output_delimiter
             << (static_cast<double>(gene.span_variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter;

    if (gene.nucleotides > 0) {

      out_file << (static_cast<double>(gene.variant_count) / static_cast<double>(gene.nucleotides)) << output_delimiter
               << (static_cast<double>(gene.EAS_variant_count) / static_cast<double>(gene.nucleotides)) << output_delimiter
               << (static_cast<double>(gene.EUR_variant_count) / static_cast<double>(gene.nucleotides)) << output_delimiter
               << (static_cast<double>(gene.SAS_variant_count) / static_cast<double>(gene.nucleotides)) << output_delimiter
               << (static_cast<double>(gene.AMR_variant_count) / static_cast<double>(gene.nucleotides)) << output_delimiter
               << (static_cast<double>(gene.AFR_variant_count) / static_cast<double>(gene.nucleotides)) << output_delimiter;

    } else {

      out_file << 0 << output_delimiter
               << 0 << output_delimiter
               << 0 << output_delimiter
               << 0 << output_delimiter
               << 0 << output_delimiter
               << 0 << output_delimiter;

    }

    out_file << gene.heterozygous << output_delimiter
             << gene.homozygous << output_delimiter;
    double ratio = static_cast<double>(gene.heterozygous) / (static_cast<double>(gene.homozygous) + homozygous_bias);
    out_file << ratio << output_delimiter
             << (gene.indel * 100) << output_delimiter
             << (gene.transition * 100) << output_delimiter
             << (gene.transversion * 100) << '\n';

  } // Gene

  return true;

}


void kgl::GenomeMutation::writeHeader(std::ostream& out_file, char output_delimiter) const {

  out_file << "Genome" << output_delimiter
           << "Contig" << output_delimiter
           << "Gene" << output_delimiter
           << "Name" << output_delimiter
           << "Description" << output_delimiter
           << "BioType" << output_delimiter
           << "ValidProtein" << output_delimiter
           << "GafId" <<  output_delimiter
           << "Begin" << output_delimiter
           << "End" << output_delimiter
           << "Span" << output_delimiter
           << "Strand" << output_delimiter
           << "Sequences" << output_delimiter
           << "SeqName" << output_delimiter
           << "Nucleotides" << output_delimiter
           << "Exons" << output_delimiter
           << "Attributes" << output_delimiter
           << "UniqueVariants" << output_delimiter
           << "VariantCount" << output_delimiter
           << "SpanVariantCount" << output_delimiter
           << "HetLoF" << output_delimiter
           << "HomLoF" << output_delimiter
           << "F/MPhase" << output_delimiter
           << "F/MVariant" << output_delimiter
           << "GenomeCount" << output_delimiter
           << "GenomeVariant" << output_delimiter
           << "SpanDensity" << output_delimiter
           << "VariantDensity" << output_delimiter
           << "EASDensity" << output_delimiter
           << "EURDensity" << output_delimiter
           << "SASDensity" << output_delimiter
           << "AMRDensity" << output_delimiter
           << "AFRDensity" << output_delimiter
           << "Heterozygous" << output_delimiter
           << "Homozygous" << output_delimiter
           << "Het/Hom" << output_delimiter
           << "Indel%" << output_delimiter
           << "Transition%" << output_delimiter
           << "Transversion%" << '\n';

}

