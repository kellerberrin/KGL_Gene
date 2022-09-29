///
// Created by kellerberrin on 3/10/17.
//

#ifndef KGL_GFF_FASTA_H
#define KGL_GFF_FASTA_H



#include "kgl_genome_genome.h"
#include "kgl_gff3.h"
#include "kgl_fasta.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Read and write Fasta files.
// Static object to provide data hiding and namespace.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseGffFasta {

public:

  ParseGffFasta() =delete;
  ~ParseGffFasta() = delete;

  [[nodiscard]] static std::shared_ptr<GenomeReference> readFastaGffFile( const std::string& organism,
                                                                          const std::string& fasta_file_name,
                                                                          const std::string& gff_file_name);

private:


};





}   // end namespace


#endif