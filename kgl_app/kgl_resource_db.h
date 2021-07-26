//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_RESOURCE_DB_H
#define KGL_RESOURCE_DB_H

#include "kel_exec_env.h"

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ResourceBase inherited by all resources.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



enum class RuntimeResourceType { GENOME_DATABASE,       // Fasta and GFF resources
                                 ONTOLOGY_DATABASE,     // GO ontology data files.
                                 GENE_NOMENCLATURE,     // Homo Sapien gene equivalent naming codes (symbol, ensembl etc).
                                 GENOME_GENEALOGY,      // Genome aux info including genealogy (1000Genomes only).
                                 GENOME_AUX_INFO,       // Genome aux info, sex, population and super-population data (all individual genomes).
                                 ALLELE_CITATION };     // PMID citation identifiers, indexed by allele rsid ('rsXXXXXXXXX').

class ResourceBase {

public:

  explicit ResourceBase(std::string identifier) : identifier_(std::move(identifier)) {}
  virtual ~ResourceBase() = default;

  [[nodiscard]] const std::string& identifier() const { return identifier_; }
  [[nodiscard]] virtual RuntimeResourceType getResourceType() const = 0;

  // Gene Nomenclature identifiers.
  static const constexpr char* NOMENCLATURE_UNIPROTID{"UniprotID"};   // The uniprot nomenclature file, class id 'UniprotResource'
  static const constexpr char* NOMENCLATURE_ENSEMBL{"EnsemblHGNC"};   // The ensembl nomenclature file, class id 'EnsemblHGNCResource'

private:

  const std::string identifier_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container object to hold resources as they are passed to analysis packages.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ResourceMap = std::multimap<RuntimeResourceType, std::shared_ptr<const ResourceBase>>;


class AnalysisResources {

public:

  AnalysisResources() = default;
  ~AnalysisResources() = default;

  void addResource(const std::shared_ptr<const ResourceBase>& resource_ptr) {

    resource_map_.emplace(resource_ptr->getResourceType(), resource_ptr);

  }


  [[nodiscard]] std::vector<std::shared_ptr<const ResourceBase>> getResources( RuntimeResourceType resource,
                                                                               const std::string& resource_ident = "") const;

  [[nodiscard]] const ResourceMap& getMap() const { return resource_map_; }

  template <class ResourceClass>
  [[nodiscard]] std::shared_ptr<const ResourceClass> getSingleResource(RuntimeResourceType resource, std::string resource_ident = "") const {

    auto resource_vector = getResources(resource, resource_ident);
    if (resource_vector.size() != 1) {

      ExecEnv::log().critical( "Request Resource: {}, Ident: {} expected 1 resource, found: {} resources - unrecoverable error",
                               resourceDescription(resource), resource_ident, resource_vector.size());

    }

    auto resource_ptr = std::dynamic_pointer_cast<const ResourceClass>(resource_vector.front());

    if (not resource_ptr) {

      ExecEnv::log().critical( "Request Resource: {}, Ident: {} invalid resource type found - unrecoverable error",
                               resourceDescription(resource), resource_ident);

    }

    return resource_ptr;

  }

private:

  ResourceMap resource_map_;

  inline static std::vector<std::pair<RuntimeResourceType, std::string>> resource_description = {

      { RuntimeResourceType::GENOME_DATABASE, "RuntimeResourceType::GENOME_DATABASE"} ,      // Fasta and GFF resources
      { RuntimeResourceType::ONTOLOGY_DATABASE,"RuntimeResourceType::ONTOLOGY_DATABASE"},     // GO ontology data files.
      { RuntimeResourceType::GENE_NOMENCLATURE, "RuntimeResourceType::GENE_NOMENCLATURE"},    // Homo Sapien gene equivalent naming codes (symbol, ensembl etc).
      { RuntimeResourceType::GENOME_GENEALOGY, "RuntimeResourceType::GENOME_GENEALOGY"},      // Genome aux info including genealogy (1000Genomes only).
      { RuntimeResourceType::GENOME_AUX_INFO, "RuntimeResourceType::GENOME_AUX_INFO"},     // Genome aux info, sex, population and super-population data (all individual genomes).
      { RuntimeResourceType::ALLELE_CITATION, "RuntimeResourceType::ALLELE_CITATION"}     // PMID citation identifiers, indexed by allele rsid ('rsXXXXXXXXX').

  };

  static std::string resourceDescription(RuntimeResourceType resource);

};



}   // end namespace

#endif //KGL_RESOURCE_DB_H
