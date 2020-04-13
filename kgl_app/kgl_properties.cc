//
// Created by kellerberrin on 11/11/18.
//

#include "kgl_properties.h"
#include "kgl_table_impl.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// High level application specific property retrieval
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::RuntimeProperties::readProperties(const std::string& properties_file) {

  std::string properties_path = Utility::filePath(properties_file, work_directory_);
  return property_tree_.readProperties(properties_path);

}



void kgl::RuntimeProperties::getGenomeDBFiles(const std::string& organism,
                                              std::string& fasta_file,
                                              std::string& gff_file,
                                              std::string& gaf_file,
                                              std::string& translationTable) const {



  // The GAF file is optional.
  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + organism + std::string(DOT_) + std::string(GAF_FILE_) + std::string(DOT_) + std::string(VALUE_);

  if (not property_tree_.getOptionalFileProperty(key, work_directory_, gaf_file)) {

    ExecEnv::log().info("Optional runtime XML property not present: {}", key);

  }

  key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + organism + std::string(DOT_) + std::string(FASTA_FILE_) + std::string(DOT_) + std::string(VALUE_);
  if (not property_tree_.getFileProperty(key, work_directory_, fasta_file)) {

    ExecEnv::log().critical("Required Fasta File: {} does not exist.", fasta_file);

  }

  key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + organism + std::string(DOT_) + std::string(GFF_FILE_) + std::string(DOT_) + std::string(VALUE_);
  if (not property_tree_.getFileProperty(key, work_directory_, gff_file)) {

    ExecEnv::log().critical("Required GFF3 File: {} does not exist.", gff_file);

  }

  key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + organism + std::string(DOT_) + std::string(TRANSLATION_TABLE_) + std::string(DOT_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, translationTable)) {

    translationTable = Tables::STANDARDTABLE->table_name; // The standard amino acid translation table.
    ExecEnv::log().warn("Amino Acid translation table not specified, using standard table: {}", translationTable);

  }

}



bool kgl::RuntimeProperties::getMixtureFile(std::string& mixture_file) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(MIXTURE_FILE_) + std::string(DOT_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, mixture_file)) {

    return false;

  }

  mixture_file = Utility::filePath(mixture_file, work_directory_);

  return Utility::fileExists(mixture_file);

}


bool kgl::RuntimeProperties::getActiveGenomes(std::vector<std::string>& genome_list) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(GENOME_LIST_) + std::string(DOT_) + std::string(ACTIVE_);
  if (not property_tree_.getPropertyVector(key, genome_list)) {

    return false;

  }

  for (auto genome : genome_list) {

    ExecEnv::log().info("Loading Genome: {}", genome);

  }

  return true;

}


bool kgl::RuntimeProperties::getGenomeAuxFiles(const std::string& organism, std::vector<AuxFileInfo>& auxfiles) const {

  auxfiles.clear();

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + organism;

//  std::string key = std::string(RUNTIME_ROOT_) + organism + dot + AUX_GENOME_FILE_;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    return false;

  }

  for (auto sub_tree : property_tree_vector) {

    if (sub_tree.first == AUX_GENOME_FILE_) {

      std::string aux_file;
      if (not sub_tree.second.getFileProperty(AuxFileInfo::AUX_FILE_NAME_, work_directory_, aux_file)) {

        ExecEnv::log().error("RuntimeProperties::getGenomeAuxFiles, problem retrieving auxiliary file name");
        sub_tree.second.treeTraversal();
        continue;

      }

      std::string aux_type;
      if (not sub_tree.second.getProperty(AuxFileInfo::AUX_FILE_TYPE_, aux_type)) {

        ExecEnv::log().error("RuntimeProperties::getGenomeAuxFiles, problem retrieving auxiliary file type");
        sub_tree.second.treeTraversal();
        continue;

      }

      auxfiles.emplace_back(AuxFileInfo(aux_file, aux_type));

    }

  }

  return true;

}



std::vector<kgl::VCFFileInfo> kgl::RuntimeProperties::getVCFFileVector() const {


  std::vector<kgl::VCFFileInfo> vcf_files;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + VCF_LIST_ + std::string(DOT_) + std::string(ACTIVE_);

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getVCFFileVector, no VCF files specified");
    return vcf_files;

  }

  ExecEnv::log().info("RuntimeProperties::getVCFFileVector size: {}", property_tree_vector.size());

  for (const auto& sub_tree : property_tree_vector) {

    key = std::string(VCF_FILE_NAME_) + std::string(DOT_) + std::string(VALUE_);
    std::string vcf_file_name;
    if (not sub_tree.second.getFileProperty( key, workDirectory(), vcf_file_name)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFileVector, No VCF file name information");
      continue;

    }

    key = std::string(VCF_FILE_GENOME_) + std::string(DOT_) + std::string(VALUE_);
    std::string vcf_reference_genome;

    if (not sub_tree.second.getProperty( key, vcf_reference_genome)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFileVector, No reference genome information for VCF file: {}", vcf_file_name);
      continue;

    }

    key = std::string(VCF_FILE_PLOIDY_) + std::string(DOT_) + std::string(VALUE_);
    size_t vcf_ploidy;
    if (not sub_tree.second.getProperty( key, vcf_ploidy)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFileVector, No ploidy information for VCF file: {}", vcf_file_name);
      continue;

    }

    vcf_files.emplace_back(VCFFileInfo(vcf_file_name, vcf_reference_genome, vcf_ploidy));

  }

  return vcf_files;

}


bool kgl::RuntimeProperties::getPropertiesAuxFile(std::string &aux_file) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(AUX_FILE_) + std::string(DOT_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, aux_file)) {

    return false;

  }

  aux_file = Utility::filePath(aux_file, work_directory_);

  return Utility::fileExists(aux_file);

}


