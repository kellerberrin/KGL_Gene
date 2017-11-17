//
// Created by kellerberrin on 17/11/17.
//

#ifndef KGL_PHYLOGENETIC_ANALYSIS_H
#define KGL_PHYLOGENETIC_ANALYSIS_H


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_filter.h"
#include "kgl_gff_fasta.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object implements a series of high level functions for the kgl_phylogenetic application.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class ApplicationAnalysis {

public:

  ApplicationAnalysis() = default;
  virtual ~ApplicationAnalysis() = default;

  static bool writeMutantProtein(const std::string& fastaFile,
                                 const std::string& sequenceName,
                                 const ContigId_t& contig_id,
                                 const FeatureIdent_t& gene_id,
                                 const FeatureIdent_t& sequence_id,
                                 const std::shared_ptr<const GenomeDatabase>& genome_db,
                                 const std::shared_ptr<const GenomeVariant>& genome_variant);

  static bool readMutantProtein(const std::string& fastaFile,
                                const std::string& sequenceName,
                                const ContigId_t& contig_id,
                                const FeatureIdent_t& gene_id,
                                const FeatureIdent_t& sequence_id,
                                const std::shared_ptr<const GenomeDatabase>& genome_db,
                                const std::shared_ptr<const GenomeVariant>& genome_variant,
                                std::string& comparison_string);

private:


};


}   // namespace genome
}   // namespace kellerberrin





#endif //KGL_PHYLOGENETIC_ANALYSIS_H
