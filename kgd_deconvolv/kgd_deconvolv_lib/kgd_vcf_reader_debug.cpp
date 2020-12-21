
#include "kgd_global.h"
#include <iostream>      // std::cout
#include "kgd_vcf_reader.h"


namespace kgd = kellerberrin::deconvolv;


bool kgd::VcfReader::printSampleName() {

  dout << "Sample name is " << this->sampleName_ << std::endl;
  return true;

}
