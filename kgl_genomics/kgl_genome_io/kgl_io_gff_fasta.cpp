//
// Created by kellerberrin on 3/10/17.
//


#include "kgl_io_gff_fasta.h"



namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::GenomeReference> kgl::ParseGffFasta::readFastaGffFile( const std::string& organism,
                                                                            const std::string& fasta_file_name,
                                                                            const std::string& gff_file_name ) {

  std::shared_ptr<GenomeReference> genome_db_ptr = ParseFasta::readFastaFile(organism, fasta_file_name);
  ParseGff3::readGffFile(gff_file_name, *genome_db_ptr);
  return genome_db_ptr;


}




