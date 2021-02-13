//
// Created by kellerberrin on 11/11/18.
//

#include "kgl_properties.h"
#include "kgl_table_impl.h"
#include "kel_utility.h"

#include <set>

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Boilerplate code to extract structured analytic and file information from "runtime.xml"
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Runtime xml retrieval
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// A vector of active packages.



// Parse the list of VCF Alias for Chromosome/Contigs in the Genome Database.
kgl::ActivePackageVector kgl::RuntimeProperties::getActivePackages() const {


  ActivePackageVector active_packages;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(EXECUTE_LIST_) + std::string(DOT_) + std::string(ACTIVE_);

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getActivePackages, No Active Packages Specified");
    return active_packages;

  }

  for (const auto& sub_tree : property_tree_vector) {

    if (sub_tree.first == PACKAGE_) {

      std::string active_package = sub_tree.second.getValue();
      active_packages.emplace_back(active_package);

    }

  }

  return active_packages;

}


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
    key = std::string(PACKAGE_IDENT_);
    std::string package_ident;
    if (not sub_tree.second.getProperty(key,  package_ident)) {

      ExecEnv::log().error("RuntimeProperties::getPackageMap, No  Package Identifier");
      continue;

    }

    // Get a Vector of Analysis functions to perform. These must be defined in code.
    // This means we can perform more than one analysis for each time-expensive data file read.
    std::vector<SubPropertyTree> analysis_sub_trees;
    std::string key = PACKAGE_ANALYSIS_LIST_;
    if (not sub_tree.second.getPropertyTreeVector(key , analysis_sub_trees))  {

      ExecEnv::log().warn("RuntimeProperties::getPackageMap, No Analysis specified for Package: {} (use 'NULL' analysis)", package_ident);

    }
    std::vector<std::string> analysis_vector;
    for (auto const& analysis_sub_tree : analysis_sub_trees) {

      if (analysis_sub_tree.first == ANALYSIS_) {

        std::string analysis_ident = analysis_sub_tree.second.getValue();
        if (analysis_ident.empty()) {

          ExecEnv::log().critical("RuntimeProperties::getPackageMap, Package: {}, No Analysis identifier specified", package_ident);

        } else {

          analysis_vector.push_back(analysis_ident);

        }

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

      if (genome_sub_tree.first == GENOME_DATABASE_) {

        std::string genome_database = genome_sub_tree.second.getValue();
        if (genome_database.empty()) {

          ExecEnv::log().critical("RuntimeProperties::getPackageMap, Package: {}, No Genome Database identifier specified", package_ident);

        } else {

          genomes.push_back(genome_database);

        }

      }

    }

    // Get a Vector of iteration items (VCF files) to load iteratively.
    // Each of these larger data (VCF) files are iteratively loaded into memory and then presented to the analysis functions.
    // The data memory is recovered and the next file is loaded and analysed.
    std::vector<SubPropertyTree> iteration_vector;
    if (not sub_tree.second.getPropertyTreeVector(PACKAGE_ITERATION_LIST_, iteration_vector))  {

      ExecEnv::log().warn("RuntimeProperties::getPackageMap, No Iteration items (VCF Files) specified for Package: {}", package_ident);

    }

    std::vector<std::vector<std::string>> vector_iteration_files;

    for (auto const& iteration_sub_tree : iteration_vector) {

      std::vector<std::string> iteration_files;
      if (iteration_sub_tree.first == PACKAGE_ITERATION_) {

        if (not iteration_sub_tree.second.getNodeVector(DATA_FILE_IDENT_, iteration_files))  {

          ExecEnv::log().critical("RuntimeProperties::getPackageMap, No Iteration items (VCF Files) specified for Package: {}", package_ident);

        }

      }

      if (not iteration_files.empty()) {

        vector_iteration_files.push_back(iteration_files);

      }

    }

    std::pair<std::string, RuntimePackage> new_package(package_ident, RuntimePackage(package_ident, analysis_vector, genomes, vector_iteration_files));

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

    key = std::string(ANALYSIS_IDENT_);
    std::string analysis_ident;
    if (not sub_tree.second.getProperty( key, analysis_ident)) {

      ExecEnv::log().error("RuntimeProperties::getAnalysisMap, No  Analysis Identifier");
      continue;

    }

    key = std::string(PARAMETER_RUNTIME_) +  std::string(DOT_) + std::string(ACTIVE_);
    std::vector<SubPropertyTree> parameter_tree_vector;
    if (not sub_tree.second.getPropertyTreeVector(key, parameter_tree_vector)) {

      ExecEnv::log().info("RuntimeProperties::getAnalysisMap, no Parameters Specified for Analysis: {}", analysis_ident);
      return analysis_map; // return empty map.

    }

    RuntimeParameterMap parameter_map;
    for (const auto& parameter_sub_tree : parameter_tree_vector) {

      // Only process parameter records.
      if (parameter_sub_tree.first == PARAMETER_BLOCK_) {

        std::string parameter_ident = parameter_sub_tree.second.getValue();
        if (parameter_ident.empty()) {

          ExecEnv::log().error("RuntimeProperties::getAnalysisMap, No Parameter Identifier Found for Analysis: {}", analysis_ident);
          continue;

        } else {

          parameter_map.push_back(parameter_ident);

        }

      } // if parameter block.

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

    key = std::string(GENOME_IDENT_);
    std::string genome_ident;
    if (not sub_tree.second.getProperty( key, genome_ident)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No Genome Database Identifier");
      continue;

    }

    key = std::string(FASTA_FILE_);
    std::string fasta_file_name;
    if (not sub_tree.second.getFileProperty( key, workDirectory(), fasta_file_name)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No Fasta file name information");
      continue;

    }

    key = std::string(GFF_FILE_);
    std::string gff_file_name;
    if (not sub_tree.second.getFileProperty( key, workDirectory(), gff_file_name)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No Gff file name information");
      continue;

    }

    key = std::string(TRANSLATION_TABLE_);
    std::string translation_table;
    if (not sub_tree.second.getProperty( key, translation_table)) {

      ExecEnv::log().error("RuntimeProperties::getGenomeReferenceMap, No DNA/Amino translationb table specified");
      continue;

    }

    key = std::string(GAF_FILE_);
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
kgl::RuntimeDataFileMap kgl::RuntimeProperties::getDataFiles() const {

  RuntimeDataFileMap data_file_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + DATA_FILE_LIST_;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getDataFiles, no Data files specified");
    return data_file_map; // return empty map.

  }

  for (const auto& sub_tree : property_tree_vector) {

    // Process VCF file record.
    if (sub_tree.first == VCF_DATA_FILE_TYPE_) {

      key = std::string(DATA_FILE_IDENT_);
      std::string vcf_ident;
      if (not sub_tree.second.getProperty( key, vcf_ident)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No VCF Identifier");
        continue;

      }

      key = std::string(DATA_FILE_NAME_);
      std::string vcf_file_name;
      if (not sub_tree.second.getFileProperty( key, workDirectory(), vcf_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No VCF file name information");
        continue;

      }

      key = std::string(DATA_PARSER_TYPE_);
      std::string vcf_parser_type;
      if (not sub_tree.second.getProperty( key, vcf_parser_type)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No VCF file parser type information");
        continue;

      }

      key = std::string(VCF_FILE_GENOME_);
      std::string vcf_reference_genome;

      if (not sub_tree.second.getProperty( key, vcf_reference_genome)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No reference genome information for VCF file: {}", vcf_file_name);
        continue;

      }

      key = std::string(VCF_INFO_EVIDENCE_);;
      std::string evidence_ident;
      if (not sub_tree.second.getProperty(key, evidence_ident)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No VCF Info evidence specified for VCF file: {}", vcf_file_name);

      }

      std::shared_ptr<BaseFileInfo> file_info_ptr = std::make_shared<RuntimeVCFFileInfo>( vcf_ident,
                                                                                          vcf_file_name,
                                                                                          vcf_parser_type,
                                                                                          vcf_reference_genome,
                                                                                          evidence_ident);

      auto result = data_file_map.try_emplace(vcf_ident, file_info_ptr);

      if (not result.second) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; Could not add VCF file ident: {} to map (duplicate)", vcf_ident);

      }

    } else if (sub_tree.first == GENERAL_DATA_FILE_TYPE_) { // Process General file record.

      key = std::string(DATA_FILE_IDENT_);
      std::string ped_ident;
      if (not sub_tree.second.getProperty(key, ped_ident)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No File Identifier");
        continue;

      }

      key = std::string(DATA_FILE_NAME_);
      std::string ped_file_name;
      if (not sub_tree.second.getFileProperty(key, workDirectory(), ped_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No Ped file name information");
        continue;

      }

      key = std::string(DATA_PARSER_TYPE_);
      std::string ped_parser_type;
      if (not sub_tree.second.getProperty(key, ped_parser_type)) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; No Ped file parser type information");
        continue;

      }

      std::shared_ptr<BaseFileInfo> file_info_ptr = std::make_shared<BaseFileInfo>( ped_ident,
                                                                                    ped_file_name,
                                                                                    ped_parser_type);

      auto result = data_file_map.try_emplace(ped_ident, file_info_ptr);

      if (not result.second) {

        ExecEnv::log().error("RuntimeProperties::getDataFiles; Could not add PED file ident: {} to map (duplicate)", ped_ident);

      }

    }

  } // for

  return data_file_map;

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


kgl::VariantEvidenceMap kgl::RuntimeProperties::getEvidenceMap() const {

  VariantEvidenceMap variant_evidence_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(EVIDENCE_LIST_);

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().info("RuntimeProperties::getEvidenceMap, No Variant Evidence items specified.");
    return variant_evidence_map;

  }

  for (const auto& sub_tree : property_tree_vector) {

    std::string evidence_ident;
    key = std::string(EVIDENCE_IDENT_);
    if (not sub_tree.second.getProperty(key, evidence_ident)) {

      ExecEnv::log().error("RuntimeProperties::getEvidenceMap, No Evidence Identifier specified for evidence list.");
      continue;

    }

    key = EVIDENCE_INFO_LIST_;
    std::vector<SubPropertyTree> info_tree_vector;
    if (not sub_tree.second.getPropertyTreeVector(key, info_tree_vector)) {

      ExecEnv::log().info("RuntimeProperties::getEvidenceMap, no info items specified for evidence ident: {}", evidence_ident);

    }

    std::set<std::string> evidence_list;
    for (const auto& info_sub_tree : info_tree_vector) {

      // Only process info records.
      if (info_sub_tree.first != EVIDENCE_INFO_ITEM_) continue;

      std::string info = info_sub_tree.second.getValue();

      auto result = evidence_list.insert(info);
      if (not result.second) {

        ExecEnv::log().warn("RuntimeProperties::getEvidenceMap, Duplicate Info item: {} specified for evidence ident: {}", info, evidence_ident);

      }

    }  // for info item

    variant_evidence_map.setEvidence(evidence_ident, evidence_list);

  }

  return variant_evidence_map;

}



kgl::ActiveParameterList kgl::RuntimeProperties::getParameterMap() const {

  ActiveParameterList defined_named_parameters;

  ContigAliasMap contig_alias_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + std::string(PARAMETER_LIST_);

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().warn("RuntimeProperties::getParameterMap; No Parameter Blocks specified");

  } else {

    for (auto const& subtree : property_tree_vector) {

      auto const& [sub_tree_tag, sub_tree] = subtree;

      // Ignore any comments.
      if (sub_tree_tag == PARAMETER_BLOCK_) {

        NamedParameterVector named_parameter_vector;
        // Get the block name.
        std::string block_name;
        if (not sub_tree.getProperty(PARAMETER_NAME_, block_name)) {

          ExecEnv::log().error("RuntimeProperties::getParameterMap; No block <parameterName> specified, block skipped");
          continue;

        }

        // Store the block name.
        named_parameter_vector.first = Utility::trimEndWhiteSpace(block_name);

        std::vector<SubPropertyTree> property_block_vector;
        if (not sub_tree.getPropertySubTreeVector(property_block_vector)) {

          ExecEnv::log().warn("RuntimeProperties::getParameterMap; No Parameter Vectors specified");

        }

        ParameterVector parameter_vector;
        for (auto const &block_tree : property_block_vector) {

          auto& [block_tree_tag, block_sub_tree] = block_tree;
          if (block_tree_tag == PARAMETER_VECTOR_) {

            std::vector<SubPropertyTree> property_vector;
            if (not block_sub_tree.getPropertySubTreeVector(property_vector)) {

              continue;

            }

            // Create a parameter map for each parameter vector.
            ParameterMap parameter_map;
            // Unpack the parameter vector.
            for (auto const& vector_item : property_vector) {

              auto const& [item_tag, item_tree] = vector_item;

              // Ignore help and info tags.
              if (item_tag == PARAMETER_) {

                std::string parameter_ident;
                if (not item_tree.getProperty(PARAMETER_IDENT_, parameter_ident)) {

                  ExecEnv::log().error("RuntimeProperties::getParameterMap; No block <parameterIdent> specified, parameter skipped");
                  continue;

                }
                std::string parameter_value;
                if (not item_tree.getProperty(PARAMETER_VALUE_, parameter_value)) {

                  ExecEnv::log().error("RuntimeProperties::getParameterMap; No block <parameterValue> specified, parameter skipped");
                  continue;

                }

                parameter_ident = Utility::trimEndWhiteSpace(parameter_ident);
                parameter_value = Utility::trimEndWhiteSpace(parameter_value);
                parameter_map.insert(parameter_ident, parameter_value);

              } // if Parameter

            } // for all vector tags.

            parameter_vector.push_back(parameter_map);

          } // if vector

        } // for all block tags

        named_parameter_vector.second = parameter_vector;
        defined_named_parameters.addNamedParameterVector(named_parameter_vector);

      } // if block.

    } // for all blocks.

  } // if any blocks.

  return defined_named_parameters;

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





