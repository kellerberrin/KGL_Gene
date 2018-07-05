//
// Created by kellerberrin on 30/06/18.
//

#include "kgd_ctl_data.h"
#include "kgd_vcf_reader.h"
#include "kgd_txt_reader.h"

namespace kgd = kellerberrin::deconvolv;


bool kgd::MixtureDataObj::readVCFPlafExclude(const std::string& vcf_filename,
                                          const std::string& plaf_filename,
                                          const std::string& exclude_filename) {

  std::shared_ptr<ExcludeMarker> excluded_reader_ptr(std::make_shared<ExcludeMarker>());
  excluded_reader_ptr->readFromFile(exclude_filename.c_str());

  VcfReader vcf_reader(vcf_filename);
  vcf_reader.findAndKeepMarkers(excluded_reader_ptr);
  vcf_reader.finalize(); // Finalize after remove variantlines

  refCount_ = vcf_reader.getRefCount();
  altCount_ = vcf_reader.getAltCount();

  TxtReader plaf;
  plaf.readFromFile(plaf_filename.c_str());
  plaf.findAndKeepMarkers(excluded_reader_ptr);

  plaf_ = plaf.getInfo();

  chrom_ = plaf.getChrom();
  position_ = plaf.getPosition();
  indexOfChromStarts_ = plaf.getIndexChromStarts();

  // nLoci() can only be used after plaf is set.
  if (nLoci() != refCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading VCF file: {}, Exclude file: {}, PLAF file: {}",
                         vcf_filename, exclude_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  if (nLoci() != altCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading VCF file: {}, Exclude file: {}, PLAF file: {}",
                         vcf_filename, exclude_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  return true;

}


bool kgd::MixtureDataObj::readVCFPlaf(const std::string& vcf_filename,
                                   const std::string& plaf_filename) {

  std::shared_ptr<ExcludeMarker> excluded_reader_ptr = std::make_shared<ExcludeMarker>();

  VcfReader vcf_reader(vcf_filename);
  vcf_reader.finalize(); // Finalize after remove variantlines

  refCount_ = vcf_reader.getRefCount();
  altCount_ = vcf_reader.getAltCount();

  TxtReader plaf;
  plaf.readFromFile(plaf_filename.c_str());

  plaf_ = plaf.getInfo();

  chrom_ = plaf.getChrom();
  position_ = plaf.getPosition();
  indexOfChromStarts_ = plaf.getIndexChromStarts();

  // nLoci() can only be used after plaf is set.
  if (nLoci() != refCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading VCF file: {}, Exclude file: {}, PLAF file: {}",
                         vcf_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  if (nLoci() != altCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading VCF file: {}, Exclude file: {}, PLAF file: {}",
                         vcf_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  return true;

}


bool kgd::MixtureDataObj::readRefAltPlafExclude(const std::string& ref_filename,
                                             const std::string& alt_filename,
                                             const std::string& plaf_filename,
                                             const std::string& exclude_filename) {

  std::shared_ptr<ExcludeMarker> excluded_reader_ptr(std::make_shared<ExcludeMarker>());
  excluded_reader_ptr->readFromFile(exclude_filename.c_str());

  TxtReader ref;
  ref.readFromFile(ref_filename.c_str());
  ref.findAndKeepMarkers(excluded_reader_ptr);

  refCount_ = ref.getInfo();

  TxtReader alt;
  alt.readFromFile(alt_filename.c_str());
  alt.findAndKeepMarkers(excluded_reader_ptr);

  altCount_ = alt.getInfo();

  TxtReader plaf;
  plaf.readFromFile(plaf_filename.c_str());
  plaf.findAndKeepMarkers(excluded_reader_ptr);

  plaf_ = plaf.getInfo();

  chrom_ = plaf.getChrom();
  position_ = plaf.getPosition();
  indexOfChromStarts_ = plaf.getIndexChromStarts();

  // nLoci() can only be used after plaf is set.
  if (nLoci() != refCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading Ref file: {}, Alt file: {}, Exclude file: {}, PLAF file: {}",
                         ref_filename, alt_filename, exclude_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  if (nLoci() != altCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading Ref file: {}, Alt file: {}, Exclude file: {}, PLAF file: {}",
                         ref_filename, alt_filename, exclude_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  return true;

}



bool kgd::MixtureDataObj::readRefAltPlaf(const std::string& ref_filename,
                                      const std::string& alt_filename,
                                      const std::string& plaf_filename) {

  TxtReader ref;
  ref.readFromFile(ref_filename.c_str());

  refCount_ = ref.getInfo();

  TxtReader alt;
  alt.readFromFile(alt_filename.c_str());

  altCount_ = alt.getInfo();

  TxtReader plaf;
  plaf.readFromFile(plaf_filename.c_str());

  plaf_ = plaf.getInfo();

  chrom_ = plaf.getChrom();
  position_ = plaf.getPosition();
  indexOfChromStarts_ = plaf.getIndexChromStarts();

  // nLoci() can only be used after plaf is set.
  if (nLoci() != refCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading Ref file: {}, Alt file: {}, PLAF file: {}",
                         ref_filename, alt_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  if (nLoci() != altCount_.size()) {

    // Fatal error program terminates.
    ExecEnv::log().error("Error reading Ref file: {}, Alt file: {}, PLAF file: {}",
                         ref_filename, alt_filename, plaf_filename);
    ExecEnv::log().critical("Mismatch between ref count size: {}, alt count size: {}, plaf count size:{}",
                            refCount_.size(), altCount_.size(), plaf_.size());

  }

  return true;

}