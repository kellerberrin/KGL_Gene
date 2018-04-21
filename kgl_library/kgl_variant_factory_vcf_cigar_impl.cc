//
// Created by kellerberrin on 11/03/18.
//

#include "kgl_variant_factory_vcf_cigar_impl.h"
#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


bool kgl::ParseCigarImpl::parseCigarItems(const std::string& genome_name,
                                            std::shared_ptr<const ContigFeatures> contig_ptr,
                                            const std::vector<CigarEditItem>& parsed_cigar,
                                            ContigOffset_t contig_offset,
                                            const std::string& reference,
                                            const std::string& alternate,
                                            Phred_t quality,
                                            const std::string& info,
                                            size_t& record_variants)  {


  size_t reference_index = 0;
  size_t alternate_index = 0;
  bool result = false;

  // Generate variants.
  for (auto cigar_item : parsed_cigar) {

    switch(cigar_item.second) {


      case CigarEditType::UNCHANGED: // 'M'
        result = parseCheck(cigar_item.first,
                            contig_ptr,
                            reference,
                            alternate,
                            reference_index,
                            alternate_index,
                            contig_offset);
        break;

      case CigarEditType::CHANGED: // 'X'
        result = parseSNP(cigar_item.first,
                          genome_name,
                          contig_ptr,
                          quality,
                          info,
                          reference,
                          alternate,
                          reference_index,
                          alternate_index,
                          contig_offset,
                          record_variants);
        break;

      case CigarEditType::INSERT: // 'I'
        result = parseInsert(cigar_item.first,
                             genome_name,
                             contig_ptr,
                             quality,
                             info,
                             alternate,
                             contig_offset,
                             alternate_index,
                             record_variants);
        break;

      case CigarEditType::DELETE: // 'D'
        result = parseDelete(cigar_item.first,
                             genome_name,
                             contig_ptr,
                             quality,
                             info,
                             reference,
                             reference_index,
                             contig_offset,
                             record_variants);
        break;

    }

    if (not result) {

      ExecEnv::log().error("VCF file, problem parsing cigar element {}:{}", cigar_item.first, static_cast<char>(cigar_item.second));
      return false;

    }

  }

  return true;

}


bool kgl::ParseCigarImpl::parseCheck(size_t cigar_count,
                                       std::shared_ptr<const ContigFeatures> contig_ptr,
                                       const std::string& reference,
                                       const std::string& alternate,
                                       size_t& reference_index,
                                       size_t& alternate_index,
                                       ContigOffset_t& contig_offset) const {

  for (size_t idx = 0; idx < cigar_count; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_index, reference[reference_index], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
      return false;

    }

    if (alternate[alternate_index] != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match VCF record alternate[{}] = base: {}",
                           reference_index, reference[reference_index], alternate_index, alternate[alternate_index]);
      return false;

    }

    ++reference_index;
    ++alternate_index;
    ++contig_offset;

  }

  return true;

}


bool kgl::ParseCigarImpl::parseSNP(size_t cigar_count,
                                     const std::string& variant_source,
                                     std::shared_ptr<const ContigFeatures> contig_ptr,
                                     Phred_t quality,
                                     const std::string&, // info,
                                     const std::string& reference,
                                     const std::string& alternate,
                                     size_t& reference_index,
                                     size_t& alternate_index,
                                     ContigOffset_t& contig_offset,
                                     size_t& variant_count) {

  for (size_t idx = 0; idx < cigar_count; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_index, reference[reference_index], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
      return false;

    }

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));

    std::shared_ptr<SNPVariant> snp_variant_ptr(std::make_shared<SNPVariant>(variant_source,
                                                                             contig_ptr->contigId(),
                                                                             contig_offset,
                                                                             quality,
                                                                             evidence_ptr,
                                                                             DNA5::convertChar(reference[reference_index]),
                                                                             DNA5::convertChar(alternate[alternate_index])));

    variant_count += addThreadSafeGenomeVariant(snp_variant_ptr); // Annotate with genome information

    ++reference_index;
    ++alternate_index;
    ++contig_offset;

  }

  return true;

}


bool kgl::ParseCigarImpl::parseInsert(size_t cigar_count,
                                        const std::string& variant_source,
                                        std::shared_ptr<const ContigFeatures> contig_ptr,
                                        Phred_t quality,
                                        const std::string&, // info,
                                        const std::string& alternate,
                                        ContigOffset_t contig_offset,
                                        size_t& alternate_index,
                                        size_t& variant_count)  {

  CompoundVariantMap compound_variant_map;

  for (size_t idx = 0; idx < cigar_count; ++idx) {

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));

    std::shared_ptr<InsertVariant> insert_variant_ptr(std::make_shared<InsertVariant>(variant_source,
                                                                                      contig_ptr->contigId(),
                                                                                      contig_offset,
                                                                                      quality,
                                                                                      evidence_ptr,
                                                                                      contig_ptr->sequence().at(contig_offset),
                                                                                      DNA5::convertChar(alternate[alternate_index])));


    std::pair<ContigOffset_t, std::shared_ptr<const SingleVariant>> insert_pair(contig_offset, insert_variant_ptr);
    auto result = compound_variant_map.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("parseInsert(), unable to insert variant with duplicate offset: {}", contig_offset);

    }

    ++alternate_index;
    ++contig_offset;

  }

  if (compound_variant_map.size() > 1) {

    variant_count += addThreadSafeGenomeVariant(CompoundInsertFactory().createCompoundVariant(compound_variant_map));

  } else if (compound_variant_map.size() == 1) {

    std::shared_ptr<Variant> single_variant = std::const_pointer_cast<SingleVariant>(compound_variant_map.begin()->second);
    variant_count += addThreadSafeGenomeVariant(single_variant); // Annotate with genome information

  } else {

    ExecEnv::log().error("parseInsert(), cigar size: {} but zero insert variant generated", cigar_count);

  }


  return true;

}

bool kgl::ParseCigarImpl::parseDelete(size_t cigar_count,
                                        const std::string& variant_source,
                                        std::shared_ptr<const ContigFeatures> contig_ptr,
                                        Phred_t quality,
                                        const std::string&, // info,
                                        const std::string& reference,
                                        size_t& reference_index,
                                        ContigOffset_t& contig_offset,
                                        size_t& variant_count) {

  CompoundVariantMap compound_variant_map;

  for (size_t idx = 0; idx < cigar_count; ++idx) {

    if (DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)) != reference[reference_index]) {

      ExecEnv::log().error("VCF record reference[{}] = base: {} does not match contig: {}[{}] = base: {}",
                           reference_index, reference[reference_index], contig_ptr->contigId(),
                           contig_offset, DNA5::convertToChar(contig_ptr->sequence().at(contig_offset)));
      return false;

    }

//    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>(info, quality));
    std::shared_ptr<VCFEvidence> evidence_ptr(std::make_shared<VCFEvidence>("", quality));


    std::shared_ptr<DeleteVariant> delete_variant_ptr(std::make_shared<DeleteVariant>(variant_source,
                                                                                      contig_ptr->contigId(),
                                                                                      contig_offset,
                                                                                      quality,
                                                                                      evidence_ptr,
                                                                                      contig_ptr->sequence().at(contig_offset)));

    std::pair<ContigOffset_t, std::shared_ptr<const SingleVariant>> insert_pair(contig_offset, delete_variant_ptr);
    auto result = compound_variant_map.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("parseDelete(), unable to insert variant with duplicate offset: {}", contig_offset);

    }

    ++reference_index;
    ++contig_offset;

  }

  if (compound_variant_map.size() > 1) {

    variant_count += addThreadSafeGenomeVariant(CompoundDeleteFactory().createCompoundVariant(compound_variant_map));

  } else if (compound_variant_map.size() == 1) {

    std::shared_ptr<Variant> single_variant = std::const_pointer_cast<SingleVariant>(compound_variant_map.begin()->second);
    variant_count += addThreadSafeGenomeVariant(single_variant); // Annotate with genome information

  } else {

    ExecEnv::log().error("parseDelete(), cigar size: {} but zero insert variant generated", cigar_count);

  }

  return true;

}
