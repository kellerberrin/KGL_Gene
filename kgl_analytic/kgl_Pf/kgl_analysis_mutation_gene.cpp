//
// Created by kellerberrin on 15/1/21.
//

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

        const std::shared_ptr<const CodingSequence> coding_sequence = sequence_array->getFirst();
        gene_mutation.nucleotides = coding_sequence->codingNucleotides();
        gene_mutation.exons = coding_sequence->exons();

      }
      auto const& attribute_map = gene_ptr->getAttributes().getMap();
      gene_mutation.attribute_size = attribute_map.size();

      gene_vector_.push_back(gene_mutation);

    } // Gene.

  } // Contig.

  return true;

}


bool kgl::GenomeMutation::variantAnalysis100(const std::shared_ptr<const PopulationDB>& population_ptr,
                                             const std::shared_ptr<const GenomePEDData>& ped_data) {

  for (auto& gene_mutation : gene_vector_) {

    for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

      auto contig_opt = genome_ptr->getContig(gene_mutation.contig);
      if (contig_opt) {

        auto contig_ptr = contig_opt.value();

        if (not contig_ptr->getMap().empty()) {

          OffsetVariantMap variant_map;
          contig_ptr->getSortedVariants( VariantSequence::UNPHASED,
                                         gene_mutation.gene_begin,
                                         gene_mutation.gene_end,
                                         variant_map);

          gene_mutation.variant_count += variant_map.size();

          auto result = ped_data->getMap().find(genome_id);

          if (result == ped_data->getMap().end()) {

            ExecEnv::log().error("GenomeMutation::variantAnalysis100; Genome sample: {} does not have a PED record", genome_id);
            continue;

          }

          auto const& [sample_id, ped_record] = *result;

          std::string super_pop = ped_record.superPopulation();

          if (super_pop == "EAS") {

            gene_mutation.EAS_variant_count += variant_map.size();

          } else if (super_pop == "EUR") {

            gene_mutation.EUR_variant_count += variant_map.size();

          } else if (super_pop == "SAS") {

            gene_mutation.SAS_variant_count += variant_map.size();

          } else if (super_pop == "AFR") {

            gene_mutation.AFR_variant_count += variant_map.size();

          } else if (super_pop == "AMR") {

            gene_mutation.AMR_variant_count += variant_map.size();

          } else {

            ExecEnv::log().error("GenomeMutation::variantAnalysis100; Unknown super population: {}", super_pop);

          }

          auto subset_ptr = contig_opt.value()->subset(gene_mutation.gene_begin, gene_mutation.gene_end);
          for (auto const& [offset, offset_ptr] : subset_ptr->getMap()) {

            if (offset_ptr->getVariantArray().size() == 1) {

              ++gene_mutation.heterozygous;

            } else if (offset_ptr->getVariantArray().size() == 2) {

              auto offset_array = offset_ptr->getVariantArray();
              if (offset_array.front()->homozygous(*offset_array.back())) {

                ++gene_mutation.homozygous;

              }

            }

          }



          ++gene_mutation.genome_count;
          if (not variant_map.empty()) {

            ++gene_mutation.genome_variant;

            double indel{0.0}, transition{0.0}, transversion{0.0};
            for (auto const& [offset, variant] : variant_map) {

              switch(variant->variantType()) {

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

            } // for variant

            double sum = indel + transition + transversion;
            gene_mutation.indel = indel / sum;
            gene_mutation.transition = transition / sum;
            gene_mutation.transversion = transversion / sum;

          } else {

            gene_mutation.indel = 0.0;
            gene_mutation.transition = 0.0;
            gene_mutation.transversion = 0.0;

          }

        } // contig not empty

      } // if contig

    } // for genome

  } // for genes

  return true;

}


bool kgl::GenomeMutation::variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr) {

  for (auto& gene_mutation : gene_vector_) {

    for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

      auto contig_opt = genome_ptr->getContig(gene_mutation.contig);
      if (contig_opt) {

        auto contig_ptr = contig_opt.value();

        if (not contig_ptr->getMap().empty()) {

          OffsetVariantMap variant_map;
          contig_ptr->getSortedVariants( VariantSequence::UNPHASED,
                                         gene_mutation.gene_begin,
                                         gene_mutation.gene_end,
                                         variant_map);

          gene_mutation.variant_count += variant_map.size();

          auto subset_ptr = contig_ptr->subset(gene_mutation.gene_begin, gene_mutation.gene_end);
          for (auto const& [offset, offset_ptr] : subset_ptr->getMap()) {

            if (offset_ptr->getVariantArray().size() == 1) {

              ++gene_mutation.heterozygous;

            } else if (offset_ptr->getVariantArray().size() == 2) {

              auto offset_array = offset_ptr->getVariantArray();
              if (offset_array.front()->homozygous(*offset_array.back())) {

                ++gene_mutation.homozygous;

              }

            }

          }



          ++gene_mutation.genome_count;
          if (not variant_map.empty()) {

            ++gene_mutation.genome_variant;

            double indel{0.0}, transition{0.0}, transversion{0.0};
            for (auto const& [offset, variant] : variant_map) {

              switch(variant->variantType()) {

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

            } // for variant

            double sum = indel + transition + transversion;
            gene_mutation.indel = indel / sum;
            gene_mutation.transition = transition / sum;
            gene_mutation.transversion = transversion / sum;

          } else {

            gene_mutation.indel = 0.0;
            gene_mutation.transition = 0.0;
            gene_mutation.transversion = 0.0;

          }

        } // contig not empty

      } // if contig

    } // for genome

  } // for genes

  return true;

}



// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::writeOutput100(const std::string& output_file_name, char output_delimiter) const {

  const double homozygous_bias{0.1};

  std::ofstream out_file(output_file_name);

  if (not out_file.good()) {

    ExecEnv::log().error("GenomeMutation::writeOutput; could not open file: {} for output", output_file_name);
    return false;

  } else {

    ExecEnv::log().info("GenomeMutation writing output to file: {}", output_file_name);

  }

  writeHeader100(out_file, output_delimiter);

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
             << gene.nucleotides << output_delimiter
             << gene.exons << output_delimiter
             << gene.attribute_size << output_delimiter
             << gene.variant_count << output_delimiter
             << gene.genome_count << output_delimiter
             << gene.genome_variant << output_delimiter
             << (static_cast<double>(gene.variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter
             << (static_cast<double>(gene.EAS_variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter
             << (static_cast<double>(gene.EUR_variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter
             << (static_cast<double>(gene.SAS_variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter
             << (static_cast<double>(gene.AMR_variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter
             << (static_cast<double>(gene.AFR_variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter
             << gene.heterozygous << output_delimiter
             << gene.homozygous << output_delimiter;
    double ratio = static_cast<double>(gene.heterozygous) / (static_cast<double>(gene.homozygous) + homozygous_bias);
    out_file << ratio << output_delimiter
             << (gene.indel * 100) << output_delimiter
             << (gene.transition * 100) << output_delimiter
             << (gene.transversion * 100) << '\n';

  } // Gene

  return true;

}


void kgl::GenomeMutation::writeHeader100(std::ostream& out_file, char output_delimiter) const {

  out_file << "Genome" << output_delimiter << "Contig" << output_delimiter
           << "Gene" << output_delimiter << "Name" << output_delimiter
           << "Description" << output_delimiter << "BioType" << output_delimiter
           << "ValidProtein" << output_delimiter << "GafId" <<  output_delimiter
           << "Begin" << output_delimiter << "End" << output_delimiter
           << "Span" << output_delimiter << "Strand" << output_delimiter
           << "Sequences" << output_delimiter << "Nucleotides" << output_delimiter
           << "Exons" << output_delimiter << "Attributes" << output_delimiter
           << "VariantCount" << output_delimiter << "GenomeCount" << output_delimiter
           << "GenomeVariant" << output_delimiter << "VariantDensity" <<  output_delimiter
           << "EASDensity" << output_delimiter << "EURDensity" << output_delimiter
           << "SASDensity" << output_delimiter << "AMRDensity" << output_delimiter
           << "AFRDensity" << output_delimiter
           << "Heterozygous" << output_delimiter << "Homozygous" << output_delimiter
           << "Het/Hom" << output_delimiter << "Indel%" << output_delimiter
           << "Transition%" << output_delimiter << "Transversion%" << '\n';

}


// Perform the genetic analysis per iteration.
bool kgl::GenomeMutation::writeOutput(const std::string& output_file_name, char output_delimiter) const {

  const double homozygous_bias{0.1};

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
             << gene.nucleotides << output_delimiter
             << gene.exons << output_delimiter
             << gene.attribute_size << output_delimiter
             << gene.variant_count << output_delimiter
             << gene.genome_count << output_delimiter
             << gene.genome_variant << output_delimiter
             << (static_cast<double>(gene.variant_count) / static_cast<double>(gene.gene_span)) << output_delimiter
             << gene.heterozygous << output_delimiter
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

  out_file << "Genome" << output_delimiter << "Contig" << output_delimiter
           << "Gene" << output_delimiter << "Name" << output_delimiter
           << "Description" << output_delimiter << "BioType" << output_delimiter
           << "ValidProtein" << output_delimiter << "GafId" <<  output_delimiter
           << "Begin" << output_delimiter << "End" << output_delimiter
           << "Span" << output_delimiter << "Strand" << output_delimiter
           << "Sequences" << output_delimiter << "Nucleotides" << output_delimiter
           << "Exons" << output_delimiter <<  "Attributes" << output_delimiter
           << "VariantCount" << output_delimiter << "GenomeCount" << output_delimiter
           << "GenomeVariant" << output_delimiter << "VariantDensity" <<  output_delimiter
           << "Heterozygous" << output_delimiter << "Homozygous" << output_delimiter
           << "Het/Hom" << output_delimiter << "Indel%" << output_delimiter
           << "Transition%" << output_delimiter << "Transversion%" << '\n';

}
