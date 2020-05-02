//
// Created by kellerberrin on 11/11/18.
//

#include "kgl_properties.h"
#include "kgl_table_impl.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Boilerplate code to extract structured analytic and file information from "runtime.xml"
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Runtime xml retrieval
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// A map of analysis
kgl::RuntimePackageMap kgl::RuntimeProperties::getPackageMap() const {

  RuntimePackageMap package_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + PACKAGE_LIST_;
  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().error("RuntimeProperties::getPackageMap, no Analysis Packages Specified (must be a least one)");
    return package_map; // return empty map.

  }

  // For each package
  for (auto const& sub_tree : property_tree_vector) {

    // Only process Package records.
    if (sub_tree.first != PACKAGE_) continue;

    // Get the Package Ident.
    key = std::string(PACKAGE_IDENT_) + std::string(DOT_) + std::string(VALUE_);
    std::string package_ident;
    if (not sub_tree.second.getProperty(key,  package_ident)) {

      ExecEnv::log().error("RuntimeProperties::getPackageMap, No  Package Identifier");
      continue;

    }

    // Get a Vector of Analysis functions to perform. These must be defined in code.
    // This means we can perform more than one analysis for each time-expensive data file read.
    std::vector<SubPropertyTree> analysis_sub_trees;
    if (not sub_tree.second.getPropertyTreeVector( PACKAGE_ANALYSIS_LIST_, analysis_sub_trees))  {

      ExecEnv::log().warn("RuntimeProperties::getPackageMap, No Analysis specified for Package: {} (use 'NULL' analysis)", package_ident);

    }
    std::vector<std::string> analysis_vector;
    for (auto const& analysis_sub_tree : analysis_sub_trees) {

      if (analysis_sub_tree.first == VALUE_) {

        analysis_vector.push_back(analysis_sub_tree.second.getValue());

      }

    }

    // Get a Vector of Genome Databases to read.
    // This is done on initialization.
    std::vector<SubPropertyTree> genome_vector;
    if (not sub_tree.second.getPropertyTreeVector( PACKAGE_GENOME_LIST_, genome_vector))  {

      ExecEnv::log().error("RuntimeProperties::getPackageMap, No Organism Reference Genomes specified for Package: {}", package_ident);

    }
    std::vector<std::string> genomes;
    for (auto const& genome_sub_tree : genome_vector) {

      if (genome_sub_tree.first == VALUE_) {

        genomes.push_back(genome_sub_tree.second.getValue());

      }

    }

    // Get a Vector of items loaded on initialization.
    // These are generally smaller data (VCF) files that can be loaded into memory
    std::vector<SubPropertyTree> load_vector;
    if (not sub_tree.second.getPropertyTreeVector( PACKAGE_LOAD_LIST_, load_vector))  {

      ExecEnv::log().warn("RuntimeProperties::getPackageMap, No Initialization items specified for Package: {}", package_ident);

    }
    std::vector<std::string> load_files;
    for (auto const& load_sub_tree : load_vector) {

      if (load_sub_tree.first == VALUE_) {

        load_files.push_back(load_sub_tree.second.getValue());

      }

    }

    // Get a Vector of iteration items (VCF files) to load iteratively.
    // Each of these larger data (VCF) files are iteratively loaded into memory and then presented to the analysis functions.
    // The data memory is recovered and the next file is loaded and analysed.
    std::vector<SubPropertyTree> iteration_vector;
    if (not sub_tree.second.getPropertyTreeVector(PACKAGE_ITERATION_LIST_, iteration_vector))  {

      ExecEnv::log().warn("RuntimeProperties::getPackageMap, No Iteration items (VCF Files) specified for Package: {}", package_ident);

    }

    std::vector<std::string> iteration_files;
    for (auto const& iteration_sub_tree : iteration_vector) {

      if (iteration_sub_tree.first == VALUE_) {

        iteration_files.push_back(iteration_sub_tree.second.getValue());

      }

    }

    std::pair<std::string, RuntimePackage> new_package(package_ident, RuntimePackage(package_ident, analysis_vector, genomes, load_files, iteration_files));

    auto result = package_map.insert(new_package);

    if (not result.second) {

      ExecEnv::log().error("RuntimeProperties::getPackageMap, Could not add Package Ident: {} to map (duplicate)", package_ident);

    }

  }

  return package_map;

}


// A map of analysis
kgl::RuntimeAnalysisMap kgl::RuntimeProperties::getAnalysisMap() const {

  RuntimeAnalysisMap analysis_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + ANALYSIS_LIST_;
  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getAnalysisMap, no Analysis specified (NULL analysis available)");
    return analysis_map; // return empty map.

  }

  for (const auto& sub_tree : property_tree_vector) {

    // Only process Analysis records.
    if (sub_tree.first != ANALYSIS_) continue;

    key = std::string(ANALYSIS_IDENT_) + std::string(DOT_) + std::string(VALUE_);
    std::string analysis_ident;
    if (not sub_tree.second.getProperty( key, analysis_ident)) {

      ExecEnv::log().error("RuntimeProperties::getAnalysisMap, No  Analysis Identifier");
      continue;

    }

    key = PARAMETER_LIST_;
    std::vector<SubPropertyTree> parameter_tree_vector;
    if (not property_tree_.getPropertyTreeVector(key, parameter_tree_vector)) {

      ExecEnv::log().info("RuntimeProperties::getAnalysisMap, no Parameters Specified for Analysis: {}", analysis_ident);
      return analysis_map; // return empty map.

    }

    RuntimeParameterMap parameter_map;
    for (const auto& parameter_sub_tree : parameter_tree_vector) {

      // Only process parameter records.
      if (sub_tree.first != PARAMETER_) continue;

      key = std::string(PARAMETER_IDENT_) + std::string(DOT_) + std::string(VALUE_);
      std::string parameter_ident;
      if (not parameter_sub_tree.second.getProperty( key, parameter_ident)) {

        ExecEnv::log().error("RuntimeProperties::getAnalysisMap, No Parameter Identifier Found for Analysis: {}", analysis_ident);
        continue;

      }

      key = std::string(PARAMETER_VALUE_) + std::string(DOT_) + std::string(VALUE_);
      std::string parameter_value;
      if (not parameter_sub_tree.second.getProperty( key, parameter_value)) {

        ExecEnv::log().error("RuntimeProperties::getAnalysisMap, No Parameter Value Found for Parameter: {}, Analysis: {}", parameter_ident, analysis_ident);
        continue;

      }

      std::pair<std::string, std::string> new_parameter(parameter_ident, parameter_value);
      auto result = parameter_map.insert(new_parameter);

      if (not result.second) {

        ExecEnv::log().error("RuntimeProperties::getAnalysisMap, Duplicate found for Parameter: {}, Analysis: {}", parameter_ident, analysis_ident);

      }

    }  // for parameter

    std::pair<std::string, RuntimeAnalysis> new_analysis(analysis_ident, RuntimeAnalysis(analysis_ident, parameter_map));

    auto result = analysis_map.insert(new_analysis);

    if (not result.second) {

      ExecEnv::log().error("RuntimeProperties::getAnalysisMap, Could not add Analysis Ident: {} to map (duplicate)", analysis_ident);

    }

  }

  return analysis_map;

}

// A map of genomes.
kgl::RuntimeGenomeDatabaseMap kgl::RuntimeProperties::getGenomeReferenceMap() const {

  RuntimeGenomeDatabaseMap genome_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + GENOME_LIST_;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getGenomeReferenceMap, no VCF files specified");
    return genome_map; // return empty map.

  }

  for (const auto& sub_tree : property_tree_vector) {

    // Only process Genome records.
    if (sub_tree.first != GENOME_DATABASE_) continue;

    key = std::string(GENOME_IDENT_) + std::string(DOT_) + std::string(VALUE_);
    std::string genome_ident;
    if (not sub_tree.second.getProperty( key, genome_ident)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No Genome Database Identifier");
      continue;

    }

    key = std::string(FASTA_FILE_) + std::string(DOT_) + std::string(VALUE_);
    std::string fasta_file_name;
    if (not sub_tree.second.getFileProperty( key, workDirectory(), fasta_file_name)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No Fasta file name information");
      continue;

    }

    key = std::string(GFF_FILE_) + std::string(DOT_) + std::string(VALUE_);
    std::string gff_file_name;
    if (not sub_tree.second.getFileProperty( key, workDirectory(), gff_file_name)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No Gff file name information");
      continue;

    }

    key = std::string(TRANSLATION_TABLE_) + std::string(DOT_) + std::string(VALUE_);
    std::string translation_table;
    if (not sub_tree.second.getProperty( key, translation_table)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No DNA/Amino translationb table specified");
      continue;

    }

    key = std::string(GAF_FILE_) + std::string(DOT_) + std::string(VALUE_);
    std::string gaf_file_name;
    if (not sub_tree.second.getOptionalFileProperty(key, workDirectory(), gaf_file_name)) {

      gaf_file_name.clear();

    }

    std::pair<std::string, RuntimeGenomeProperty> new_genome(genome_ident, RuntimeGenomeProperty(genome_ident, fasta_file_name, gff_file_name, translation_table));

    new_genome.second.setGafFileName(gaf_file_name);

    auto result = genome_map.insert(new_genome);

    if (not result.second) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, Could not add Genome Database ident: {} to map (duplicate)", genome_ident);

    }

  }

  return genome_map;

}

// A map of VCF files.
kgl::RuntimeVCFFileMap kgl::RuntimeProperties::getVCFFiles() const {

  RuntimeVCFFileMap vcf_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + VCF_LIST_;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getVCFFiles, no VCF files specified");
    return vcf_map; // return empty map.

  }

  for (const auto& sub_tree : property_tree_vector) {

    // Only process VCF file records.
    if (sub_tree.first != VCF_FILE_) continue;

    key = std::string(VCF_IDENT_) + std::string(DOT_) + std::string(VALUE_);
    std::string vcf_ident;
    if (not sub_tree.second.getProperty( key, vcf_ident)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFiles, No VCF Identifier");
      continue;

    }

    key = std::string(VCF_FILE_NAME_) + std::string(DOT_) + std::string(VALUE_);
    std::string vcf_file_name;
    if (not sub_tree.second.getFileProperty( key, workDirectory(), vcf_file_name)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFiles, No VCF file name information");
      continue;

    }

    key = std::string(VCF_PARSER_TYPE_) + std::string(DOT_) + std::string(VALUE_);
    std::string vcf_parser_type;
    if (not sub_tree.second.getProperty( key, vcf_parser_type)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFiles, No VCF file parser type information");
      continue;

    }

    key = std::string(VCF_FILE_GENOME_) + std::string(DOT_) + std::string(VALUE_);
    std::string vcf_reference_genome;

    if (not sub_tree.second.getProperty( key, vcf_reference_genome)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFiles, No reference genome information for VCF file: {}", vcf_file_name);
      continue;

    }

    key = std::string(VCF_FILE_PLOIDY_) + std::string(DOT_) + std::string(VALUE_);
    size_t vcf_ploidy;
    if (not sub_tree.second.getProperty( key, vcf_ploidy)) {

      ExecEnv::log().error("RuntimeProperties::getVCFFiles, No ploidy information for VCF file: {}", vcf_file_name);
      continue;

    }

    std::pair<std::string, RuntimeVCFFileInfo> new_vcf(vcf_ident, RuntimeVCFFileInfo(vcf_ident, vcf_file_name, vcf_parser_type, vcf_reference_genome, vcf_ploidy));
    auto result = vcf_map.insert(new_vcf);

    if (not result.second) {

      ExecEnv::log().error("RuntimeProperties::getVCFFiles, Could not add VCF file ident: {} to map (duplicate)", vcf_ident);

    }

  }

  return vcf_map;

}

// Parse the list of VCF Alias for Chromosome/Contigs in the Genome Database.
kgl::ContigAliasMap kgl::RuntimeProperties::getContigAlias() const {


  ContigAliasMap contig_alias_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(ALIAS_LIST_);

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getContigAlias, No Contig Alias Specified");
    return contig_alias_map;

  }

  for (const auto& sub_tree : property_tree_vector) {

    std::string contig_ident;
    if (not sub_tree.second.getProperty(ALIAS_IDENT_, contig_ident)) {

      ExecEnv::log().error("RuntimeProperties::getContigAlias, No Chromosome/Contig Identifier specified for Alias.");
      continue;

    }

    // Alias is idempotent.
    contig_alias_map.setAlias(contig_ident, contig_ident);

    // Get a vector of alias
    std::vector<std::string> alias_vector;
    if (not sub_tree.second.getNodeVector(ALIAS_ENTRY_, alias_vector))  {

      ExecEnv::log().warn("RuntimeProperties::getContigAlias, No Alias Specified for Contig: {}", sub_tree.first);

    }

    for (auto const& alias : alias_vector) {

      contig_alias_map.setAlias(alias, contig_ident);

    }

  }

  return contig_alias_map;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Legacy Code (To be demolished)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::RuntimeProperties::getPropertiesAuxFile(std::string &aux_file) const {

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(AUX_FILE_) + std::string(DOT_) + std::string(VALUE_);
  if (not property_tree_.getProperty(key, aux_file)) {

    return false;

  }

  aux_file = Utility::filePath(aux_file, work_directory_);

  return Utility::fileExists(aux_file);

}

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





