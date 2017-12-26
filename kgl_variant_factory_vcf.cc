//
// Created by kellerberrin on 26/12/17.
//


#include "kgl_utility.h"
#include "kgl_sam_process.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_single.h"
#include "kgl_variant_factory_single.h"

#include <seqan/vcf_io.h>

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory::VcfFileImpl does all the heavy lifting using 3rd a party library. In this case; Seqan.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class kgl::VcfFactory::VcfFileImpl {

public:

  VcfFileImpl() = default;
  ~VcfFileImpl() = default;

  bool readVcfFile(const std::string &vcf_file_name);

private:


};


bool kgl::VcfFactory::VcfFileImpl::readVcfFile(const std::string &vcf_file_name) {

  // Open input file.
  seqan::VcfFileIn vcfIn(seqan::toCString(vcf_file_name));

  // Attach to standard output.
  //  VcfFileOut vcfOut(vcfIn);
  //  open(vcfOut, std::cout, Vcf());

  // Copy over header.
  seqan::VcfHeader header;
  readHeader(header, vcfIn);
  // writeHeader(vcfOut, header);

  // Copy the file record by record.

  size_t vcf_record_count = 0;

  seqan::VcfRecord record;
  while (!seqan::atEnd(vcfIn))
  {
    readRecord(record, vcfIn);
    ++vcf_record_count;
    // writeRecord(vcfOut, record);
  }

  ExecEnv::log().info("Read: {} records from VCF file: {}", vcf_record_count, vcf_file_name);

  return true;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::VcfFileImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::VcfFactory::VcfFactory() : vcf_file_impl_ptr_(std::make_unique<kgl::VcfFactory::VcfFileImpl>()) {}
kgl::VcfFactory::~VcfFactory() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


std::shared_ptr<kgl::GenomeVariant> kgl::VcfFactory::readParseVcf(const std::string& genome_name,
                                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                  const std::string& vcf_file_name,
                                                                  Phred_t read_quality,
                                                                  NucleotideReadCount_t min_read_count,
                                                                  double min_proportion) {

  std::shared_ptr<GenomeVariant> genome_single_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_name, genome_db_ptr);

  // To be implemented.
  ExecEnv::log().warn("VCF file factory is not yet implmented, 0 variants will be returned for VCF file :{}", vcf_file_name);

  readVcfFile(vcf_file_name);

  ExecEnv::log().info("VCF file {} has: {} raw variants", vcf_file_name, genome_single_variants->size());
  return genome_single_variants;

}

bool kgl::VcfFactory::readVcfFile(const std::string& vcf_file_name) {

  return vcf_file_impl_ptr_->readVcfFile(vcf_file_name);

}
