//
// Created by kellerberrin on 23/12/17.
//


#include "kgl_variant_single.h"
#include "kgl_variant_factory.h"
#include "kgl_alphabet_extend.h"


namespace kgl = kellerberrin::genome;


// Generate SNP variants.
std::shared_ptr<const kgl::GenomeVariant>
kgl::SingleFactory::createSingleVariants(const std::string &genome_name,
                                         const std::shared_ptr<const ContigCountData> &count_data,
                                         const std::shared_ptr<const GenomeDatabase> &genome_db,
                                         NucleotideReadCount_t minimum_read_count,
                                         double minimum_proportion,
                                         Phred_t read_quality) const {

  std::shared_ptr<GenomeVariant> genome_single_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_name, genome_db);
  size_t snp_count = 0;

  for (auto& contig_block : count_data->getMap()) {   // For each contig block.

    // Get the sequence.
    std::shared_ptr<ContigFeatures> contig_ptr;
    if (not genome_db->getContigSequence(contig_block.first, contig_ptr)) {

      ExecEnv::log().error("Contig: {} not found in SingleFactory.create()", contig_block.first);
      continue;

    } else {

      const DNA5SequenceContig& contig_sequence = contig_ptr->sequence();

      const auto &nucleotide_array = contig_block.second->getNucleotideArray();

      for (ContigOffset_t contig_offset = 0; contig_offset < nucleotide_array.contigSize(); ++contig_offset) {

        // The nucleotide_count_ptr[] is an array of 11 counts A, C, G, T, N, -, +A, +C, +G, +T, +N
        const NucleotideReadCount_t* nucleotide_count_ptr = nucleotide_array.readCount(contig_offset);

        DNA5::Alphabet reference_nucleotide = contig_sequence[contig_offset];

        // The nucleotide_count_ptr[] is an array of 11 counts A, C, G, T, N, -, +A, +C, +G, +T, +N
        // The first six A, C, G, T, N, - are used to generate SNP and delete variants.
        // The remaining five +A, +C, +G, +T, +N are used to generate insert variants.

        snp_count += GenerateSNPDelete(genome_name,
                                       contig_ptr,
                                       contig_offset,
                                       reference_nucleotide,
                                       nucleotide_count_ptr,
                                       minimum_read_count,
                                       minimum_proportion,
                                       read_quality,
                                       genome_single_variants);

        // Warning; raw pointer arithmetic - very nasty.
        // Increment the array pointer by six to access the remaining five elements of the count array.
        const NucleotideReadCount_t* nucleotide_insert_count_ptr = nucleotide_count_ptr + ExtendCountColumns::NUCLEOTIDE_COLUMNS;

        snp_count += GenerateInsert(genome_name,
                                    contig_ptr,
                                    contig_offset,
                                    reference_nucleotide,
                                    nucleotide_insert_count_ptr,
                                    minimum_read_count,
                                    minimum_proportion,
                                    read_quality,
                                    genome_single_variants);

      }  // for all sequence elements

      ExecEnv::log().info("Contig: {} has: {} raw SNPs", contig_ptr->contigId(), snp_count);
      snp_count = 0;

    } // found contig.

  }  // for all contigs.

  return genome_single_variants;

}


size_t kgl::SingleFactory::GenerateSNPDelete(const std::string &genome_name,
                                           std::shared_ptr<ContigFeatures> contig_ptr,
                                           ContigOffset_t contig_offset,
                                           DNA5::Alphabet reference_nucleotide,
                                           const NucleotideReadCount_t nucleotide_count_array[],
                                           NucleotideReadCount_t minimum_read_count,
                                           double minimum_proportion,
                                           Phred_t read_quality,
                                           std::shared_ptr<GenomeVariant> genome_single_variants) const {

  size_t variant_count = 0;
  NucleotideReadCount_t read_count = 0;
  for(ContigOffset_t idx = 0; idx <  ExtendCountColumns::NUCLEOTIDE_COLUMNS; ++idx) {

    read_count += nucleotide_count_array[idx];

  }

  for(size_t idx = 0; idx <  ExtendCountColumns::NUCLEOTIDE_COLUMNS; ++idx) {

    double proportion = static_cast<double>(nucleotide_count_array[idx]) / static_cast<double>(read_count);

    if (ExtendCountColumns::offsetToNucleotide(idx) != ExtendCountColumns::baseToExtend(reference_nucleotide)
        and nucleotide_count_array[idx] > 0
        and read_count >= minimum_read_count
        and proportion >= minimum_proportion) {

      ExtendCountColumns::Alphabet mutant_nucleotide = ExtendCountColumns::offsetToNucleotide(idx);

      std::shared_ptr<const VariantEvidence>
      evidence_ptr(std::make_shared<const ReadCountEvidence<ExtendCountColumns>>(idx, nucleotide_count_array));

      Phred_t quality = evidence_ptr->calculateQuality();

      if (ExtendCountColumns::isBaseCode(mutant_nucleotide)) {

        SNPVariant snp_variant(genome_name,
                               contig_ptr,
                               contig_offset,
                               quality,
                               evidence_ptr,
                               reference_nucleotide,
                               ExtendCountColumns::extendToBase(mutant_nucleotide));

        variant_count += addSingleVariant(genome_single_variants, snp_variant); // Annotate with genome information

      } else if (ExtendCountColumns::isDeletion(mutant_nucleotide)) {

        DeleteVariant delete_variant(genome_name,
                                     contig_ptr,
                                     contig_offset,
                                     quality,
                                     evidence_ptr,
                                     reference_nucleotide);

        variant_count += addSingleVariant(genome_single_variants, delete_variant); // Annotate with genome information

      } else {

        ExecEnv::log().error("ExtendCountColumns Unknown letter type: {}", ExtendCountColumns::convertToChar(mutant_nucleotide));

      }

    }

  }

  return variant_count;

}


size_t kgl::SingleFactory::GenerateInsert(const std::string &genome_name,
                                          std::shared_ptr<ContigFeatures> contig_ptr,
                                          ContigOffset_t contig_offset,
                                          DNA5::Alphabet reference_nucleotide,
                                          const NucleotideReadCount_t nucleotide_count_array[],
                                          NucleotideReadCount_t minimum_read_count,
                                          double minimum_proportion,
                                          Phred_t read_quality,
                                          std::shared_ptr<GenomeVariant> genome_single_variants) const {

  size_t variant_count = 0;
  NucleotideReadCount_t read_count = 0;
  for(ContigOffset_t idx = 0; idx <  DNA5::NUCLEOTIDE_COLUMNS; ++idx) {

    read_count += nucleotide_count_array[idx];

  }

  for(size_t idx = 0; idx <  DNA5::NUCLEOTIDE_COLUMNS; ++idx) {

    double proportion = static_cast<double>(nucleotide_count_array[idx]) / static_cast<double>(read_count);

    if (DNA5::offsetToNucleotide(idx) != reference_nucleotide
        and nucleotide_count_array[idx] > 0
        and read_count >= minimum_read_count
        and proportion >= minimum_proportion) {

      DNA5::Alphabet mutant_nucleotide = DNA5::offsetToNucleotide(idx);

      std::shared_ptr<const VariantEvidence>
      evidence_ptr(std::make_shared<const ReadCountEvidence<DNA5>>(idx, nucleotide_count_array));

      Phred_t quality = evidence_ptr->calculateQuality();

      InsertVariant insert_variant(genome_name,
                                   contig_ptr,
                                   contig_offset,
                                   quality,
                                   evidence_ptr,
                                   reference_nucleotide,
                                   mutant_nucleotide);

      variant_count += addSingleVariant(genome_single_variants, insert_variant); // Annotate with genome information

    }

  }

  return variant_count;

}


// This function will insert multiple variants for each CDS sequence within each gene.
size_t kgl::SingleFactory::addSingleVariant(std::shared_ptr<GenomeVariant> genome_single_variants,
                                            const Variant &variant) const {

  // Annotate the variant with genome information.
  size_t variant_count = 0;
  GeneVector gene_vector;
  ContigOffset_t variant_offset = variant.contigOffset();
  if (variant.contig()->findGenes(variant_offset, gene_vector)) {

    for (const auto& gene_ptr : gene_vector) {

      std::shared_ptr<const CodingSequenceArray> sequence_array = kgl::GeneFeature::getCodingSequences(gene_ptr);
      if (sequence_array->empty()) {

        // create a variant copy and annotate with a gene.
        std::shared_ptr<Variant> intron_single_ptr = variant.clone();
        intron_single_ptr->defineIntron(gene_ptr); // intron
        genome_single_variants->addVariant(intron_single_ptr);
        ++variant_count;

      } else {

        for (const auto& sequence : sequence_array->getMap()) {

          if (sequence.second->isWithinCoding(variant_offset)) {

            // create a variant copy and annotate with a coding sequence.
            std::shared_ptr<Variant> coding_single_ptr = variant.clone();
            coding_single_ptr->defineCoding(sequence.second); // coding
            genome_single_variants->addVariant(coding_single_ptr);
            ++variant_count;

          } else {  // an intron for this sequence

            // create a variant copy and annotate with a gene.
            std::shared_ptr<Variant> intron_single_ptr = variant.clone();
            intron_single_ptr->defineIntron(gene_ptr); // intron
            genome_single_variants->addVariant(intron_single_ptr);
            ++variant_count;

          } // if valid sequence for offset

        } // for all sequences within a gene

      } // if gene has a valid sequence.

    } // for all genes.

  } else {

    // create a variant copy and tag as non-coding.
    std::shared_ptr<Variant> noncoding_single_ptr  = variant.clone();
    noncoding_single_ptr->defineNonCoding(); // non coding
    genome_single_variants->addVariant(noncoding_single_ptr);
    ++variant_count;

  }

  return variant_count;

}


