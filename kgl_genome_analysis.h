// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// Created by kellerberrin on 13/10/17.
//

#ifndef KGL_GENOME_ANALYSIS_H
#define KGL_GENOME_ANALYSIS_H

#include "kgl_variant.h"
#include "kgl_genome_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class GenomeAnalysis {

public:

  explicit GenomeAnalysis() = default;
  ~GenomeAnalysis() = default;

  template <typename T>
  std::shared_ptr<GenomeVariant> simpleSNPVariants(std::shared_ptr<ContigCountData>& count_data,
                                                   std::shared_ptr<GenomeDatabase>& genome_db);

private:


};


template<typename T>
std::shared_ptr<GenomeVariant> GenomeAnalysis::simpleSNPVariants(std::shared_ptr<ContigCountData>& count_data,
                                                                 std::shared_ptr<GenomeDatabase>& genome_db) {

  std::shared_ptr<GenomeVariant> snp_variants(std::make_shared<GenomeVariant>());

  for (auto& contig_block : count_data->getMap()) {   // For each contig block.

    // Get the sequence.
    std::shared_ptr<ContigRecord> contig_ptr;
    if (not genome_db->getContigSequence(contig_block.first, contig_ptr)) {

      ExecEnv::log().error("Contig: {} not found in simpleSNPVariants()", contig_block.first);
      continue;

    } else {

      const Sequence_t& contig_sequence = contig_ptr->sequence();
      std::shared_ptr<ContigVariant> contig_variant_ptr(std::make_shared<ContigVariant>(contig_ptr->contigId()));

      const auto &nucleotide_array = contig_block.second->getNucleotideArray();
      for (ContigOffset_t contig_offset = 0; contig_offset < nucleotide_array.contigSize(); ++contig_offset) {

        const NucleotideReadCount_t* nucleotide_count_ptr = nucleotide_array.readCount(contig_offset);

        ContigOffset_t max_count_offset = 0;
        NucleotideReadCount_t read_count = 0;
        for(ContigOffset_t idx = 0; idx <  T::NUCLEOTIDE_COLUMNS; ++idx) {

          if (idx > 0 and nucleotide_count_ptr[idx-1] < nucleotide_count_ptr[idx]) {

            max_count_offset = idx;

          }

          read_count += nucleotide_count_ptr[idx];

        }

        typename T::NucleotideType max_count_nucleotide = T::offsetToNucleotide(max_count_offset);

        if (max_count_nucleotide != contig_sequence[contig_offset] and read_count > 0) {

          std::shared_ptr<const Variant>
          snp_variant(std::make_shared<const SNPVariant<typename T::NucleotideType>>(contig_block.first,
                                                      contig_offset,
                                                      read_count,
                                                      nucleotide_count_ptr[max_count_offset],
                                                      max_count_nucleotide,
                                                      contig_sequence[contig_offset]));

          contig_variant_ptr->addVariant(contig_offset, snp_variant);

        }

      }  // for all sequence elements

      ExecEnv::log().info("Contig: {} has: {} raw SNPs",
                          contig_variant_ptr->contigId(),
                          contig_variant_ptr->variantCount());


      if (not snp_variants->addContigVariant(contig_variant_ptr)) {

        ExecEnv::log().info("Duplicate Contig: {} variant. SNP variants not added to genome variants",
                            contig_variant_ptr->contigId());

      }
    } // found contig.

  }  // for all contigs.

  return snp_variants;

}


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_GENOME_ANALYSIS_H
