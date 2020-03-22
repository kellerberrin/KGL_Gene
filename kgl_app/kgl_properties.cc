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


  const std::string dot = ".";

  std::string key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(GAF_FILE_) + std::string(VALUE_);
  property_tree_.getOptionalFileProperty(key, work_directory_, gaf_file);

  key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(FASTA_FILE_) + std::string(VALUE_);
  if (not property_tree_.getFileProperty(key, work_directory_, fasta_file)) {

    ExecEnv::log().critical("Required Fasta File: {} does not exist.", fasta_file);

  }

  key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(GFF_FILE_) + std::string(VALUE_);
  if (not property_tree_.getFileProperty(key, work_directory_, gff_file)) {

    ExecEnv::log().critical("Required GFF3 File: {} does not exist.", gff_file);

  }

  key = std::string(RUNTIME_ROOT_) + organism + dot + std::string(TRANSLATION_TABLE_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, translationTable)) {

    translationTable = Tables::STANDARDTABLE->table_name; // The standard amino acid translation table.
    ExecEnv::log().warn("Amino Acid translation table not specified, using standard table: {}", translationTable);

  }

}


size_t kgl::RuntimeProperties::getVCFPloidy() const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(VCF_PLOIDY_) + std::string(VALUE_);
  size_t ploidy;

  if (not property_tree_.getProperty(key, ploidy)) {

    ExecEnv::log().warn("VCF ploidy is not defined (need not be organism ploidy)");
    ploidy = DEFAULT_PLOIDY_;

  }

  return ploidy;

}


bool kgl::RuntimeProperties::getMixtureFile(std::string& mixture_file) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(MIXTURE_FILE_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, mixture_file)) {

    return false;

  }

  mixture_file = Utility::filePath(mixture_file, work_directory_);

  return Utility::fileExists(mixture_file);

}


bool kgl::RuntimeProperties::getVCFFiles(std::vector<std::string>& vcf_files) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(VCF_LIST_) + std::string(FILE_LIST_);
  std::vector<std::string> truncated_vcf_vector;
  if (not property_tree_.getPropertyVector(key, truncated_vcf_vector)) {

    return false;

  }

  for (auto file : truncated_vcf_vector) {

    vcf_files.push_back(Utility::filePath(file, work_directory_));

  }

  return true;

}


bool kgl::RuntimeProperties::getActiveGenomes(std::vector<std::string>& genome_list) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(GENOME_LIST_) + std::string(ACTIVE_);
  if (not property_tree_.getPropertyVector(key, genome_list)) {

    return false;

  }

  for (auto genome : genome_list) {

    ExecEnv::log().info("Loading Genome: {}", genome);

  }

  return true;

}


bool kgl::RuntimeProperties::getGenomeAuxFiles(const std::string& organism, std::vector<AuxFileProperty>& auxfiles) const {

  auxfiles.clear();

  std::string key = std::string(RUNTIME_ROOT_) + organism;

//  std::string key = std::string(RUNTIME_ROOT_) + organism + dot + AUX_GENOME_FILE_;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    return false;

  }

  for (auto sub_tree : property_tree_vector) {

    if (sub_tree.first == AUX_GENOME_FILE_) {

      std::string aux_file;
      if (not sub_tree.second.getFileProperty(AuxFileProperty::AUX_FILE_NAME_, work_directory_, aux_file)) {

        ExecEnv::log().error("RuntimeProperties::getGenomeAuxFiles, problem retrieving auxiliary file name");
        sub_tree.second.treeTraversal();
        continue;

      }

      std::string aux_type;
      if (not sub_tree.second.getProperty(AuxFileProperty::AUX_FILE_TYPE_, aux_type)) {

        ExecEnv::log().error("RuntimeProperties::getGenomeAuxFiles, problem retrieving auxiliary file type");
        sub_tree.second.treeTraversal();
        continue;

      }

      auxfiles.emplace_back(AuxFileProperty(aux_file, aux_type));

    }

  }

  return true;

}


bool kgl::RuntimeProperties::getPropertiesAuxFile(std::string &aux_file) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(AUX_FILE_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, aux_file)) {

    return false;

  }

  aux_file = Utility::filePath(aux_file, work_directory_);

  return Utility::fileExists(aux_file);

}


