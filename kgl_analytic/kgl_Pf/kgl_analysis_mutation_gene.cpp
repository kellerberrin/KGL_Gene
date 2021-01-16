//
// Created by kellerberrin on 15/1/21.
//

#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_mutation.h"


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

      std::shared_ptr<const CodingSequenceArray> sequence_array = GeneFeature::getCodingSequences(gene_ptr);

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
      gene_mutation.gene_size = gene_ptr->sequence().length();
      gene_mutation.strand = gene_ptr->sequence().strandText();
      gene_mutation.exons = sequence_array->size();
      auto const& attribute_map = gene_ptr->getAttributes().getMap();
      gene_mutation.attribute_size = attribute_map.size();

      gene_vector_.push_back(gene_mutation);

    } // Gene.

  } // Contig.

  return true;

}


bool kgl::GenomeMutation::variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr) {

  for (auto& gene_mutation : gene_vector_) {

    for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()){

      auto contig_opt = genome_ptr->getContig(gene_mutation.contig);
      if (contig_opt) {

        if (not contig_opt.value()->getMap().empty()) {

          OffsetVariantMap variant_map;
          contig_opt.value()->getSortedVariants( VariantSequence::UNPHASED,
                                                 gene_mutation.gene_begin,
                                                 gene_mutation.gene_end,
                                                 variant_map);

          gene_mutation.variant_count += variant_map.size();

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

          }

        } // contig not empty

      } // if contig

    } // for genome

  } // for genes

  return true;

}
