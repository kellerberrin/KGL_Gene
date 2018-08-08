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

  return verifyPrint(true);

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

  return verifyPrint(true);


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

  return verifyPrint(true);

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

  return verifyPrint(true);

}

bool kgd::MixtureDataObj::verifyPrint(bool print) const {

  // Perform sanity checks.
  // The number contigs should equal the size of the position vector.
  if (getChrom().size() != getPosition().size()) {

    ExecEnv::log().error("MixtureDataObj::verifyPrint(); Contig count: {}, does not equal Position size: {}", getChrom().size(), getPosition().size());
    return false;

  }

  size_t total_positions = 0;
  for (const auto& contig_positions : getPosition()) {

    total_positions += contig_positions.size();

  }

  if (getPlaf().size() != total_positions) {

    ExecEnv::log().error("MixtureDataObj::verifyPrint(); Plaf size: {}, does not equal Variant count: {}", getPlaf().size(), total_positions);
    return false;

  }

  if (getRefCount().size() != total_positions) {

    ExecEnv::log().error("MixtureDataObj::verifyPrint(); Refcount size: {}, does not equal Variant count: {}", getRefCount().size(), total_positions);
    return false;

  }

  if (getAltCount().size() != total_positions) {

    ExecEnv::log().error("MixtureDataObj::verifyPrint(); Altcount size: {}, does not equal Variant count: {}", getAltCount().size(), total_positions);
    return false;

  }

  // Check that the offsets correspond to the positions sizes.
  if (indexOfChromStarts().size() != getPosition().size()) {

    ExecEnv::log().error("MixtureDataObj::verifyPrint(); indexOfChromStarts size: {}, does not equal does not equal Position size: {}",
                         indexOfChromStarts().size(), getPosition().size());
    return false;

  }

  size_t position_offset = 0;
  for (size_t idx = 1; idx < indexOfChromStarts().size(); ++idx) {

    position_offset += getPosition()[idx-1].size();
    if (indexOfChromStarts()[idx] != position_offset) {

      ExecEnv::log().error("MixtureDataObj::verifyPrint(); indexOfChromStarts[{}]: {}, does not equal does not summed Position[{}] size: {}",
                           idx, indexOfChromStarts()[idx], idx-1, position_offset);
      return false;

    }

  }

  if (print) {

    for (size_t idx = 0; idx < indexOfChromStarts().size(); ++idx) {

      size_t variant_count = 0;
      if (idx < (indexOfChromStarts().size() - 1)) {

        for (size_t count_idx = indexOfChromStarts()[idx]; count_idx < indexOfChromStarts()[idx+1]; ++count_idx) {

          if ((getRefCount()[count_idx] + getAltCount()[count_idx]) > 0) {

            variant_count++;

          }

        }

      } else {

        for (size_t count_idx = indexOfChromStarts()[idx]; count_idx < getRefCount().size(); ++count_idx) {

          if (getRefCount()[count_idx] + getAltCount()[count_idx] > 0) {

            variant_count++;

          }

        }

      }

      ExecEnv::log().info("MixtureDataObj::verifyPrint(); contig: {}, total variants: {}, genome variants:{}",
                          getChrom()[idx], getPosition()[idx].size(), variant_count);

    }

  }

  return true;

}