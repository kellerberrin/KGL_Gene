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

#include "kgl_genome_analysis.h"


namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeAnalysis::simpleSNPVariants(std::shared_ptr<kgl::ContigCountData>& count_data,
                                       std::shared_ptr<GenomeDatabase>& genome_db) {

  std::shared_ptr<kgl::GenomeVariant> snp_variant(std::make_shared<kgl::GenomeVariant>());

  for (auto& contig_block : count_data->getMap()) {   // For each contig block.

    const auto& nucleotide_array = contig_block.second->getNucleotideArray();
    for (ContigOffset_t contig_offset = 0; contig_offset < nucleotide_array.contigSize(); ++contig_offset) {

      NucleotideReadCount_t A_count = nucleotide_array.readCount(contig_offset, 'A');

      NucleotideReadCount_t C_count = nucleotide_array.readCount(contig_offset, 'C');
      NucleotideReadCount_t G_count = nucleotide_array.readCount(contig_offset, 'G');
      NucleotideReadCount_t T_count = nucleotide_array.readCount(contig_offset, 'T');
      NucleotideReadCount_t N_count = nucleotide_array.readCount(contig_offset, 'N');
      NucleotideReadCount_t Delete_count = nucleotide_array.readCount(contig_offset,
                                                                      StandardNucleotideColumn::DELETE_NUCLEOTIDE);
      NucleotideReadCount_t Insert_count = nucleotide_array.readCount(contig_offset,
                                                                      StandardNucleotideColumn::INSERT_SEQUENCE);

    }

  }

  return snp_variant;

}
