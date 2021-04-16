//
// Created by kellerberrin on 14/4/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H
#define KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H

#include <fstream>
#include "kol_OntologyDatabase.h"
#include "kgl_analysis_mutation_gene_stats.h"

namespace kol = kellerberrin::ontology;
namespace kellerberrin::genome {   //  organization::project level namespace



class OntologyStats {

public:

  OntologyStats() = default;
  OntologyStats(const OntologyStats &) = default;
  ~OntologyStats() = default;

  OntologyStats &operator=(const OntologyStats &) = default;

  void writeOntology(std::ostream& out_file, char output_delimiter) const;
  void writeOntologyHeader(std::ostream& out_file, char output_delimiter) const;

  void processOntologyStats(const GeneCharacteristic& gene_info, const std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr);


private:

  size_t go_term_count_{0};
  size_t go_id_count_{0};
  inline static const std::map<std::string, std::string> malaria_genes_ {
      { "P16671", "CD36"}, { "P06028","GYPB"}, { "P12318", "FCGR2A"}, { "P31994", "FCGR2B"}, { "P05362", "ICAM1"},
      {  "O14931", "NCR3"}, { "P68871", "HBB"}, { "P35228", "NOS2"}, { "P01375", "TNF"}, { "O14931", "NCR3"},
      { "Q9UNN8", "PROCR"}, { "P02730", "SLC4A1"}, { "Q9NSE2", "CISH"}, { "Q96A59", "MARVELD3"}, { "Q9Y231", "FUT9"},
      { "P19320", "VCAM1"}, { "P58753", "TIRAP"}, { "P04921", "GYPC"}, { "P0C091", "FREM3"}, { "P02724", "GYPA"},
      { "P11413", "G6PD"}, { "Q8N126", "CADM3"},  { "Q16570", "ACKR1"}, { "P23634", "ATP2B4"}, { "P17927", "CR1"},
      { "P16442", "ABO"}, {"P69905", "HBA1"}, {"P35613", "BSG"}, {"P08174", "CD55"}, {"Q8NHL6",  "LILRB1"},
      {"Q6GTX8", "LAIR1"} };
};


} // namespace


#endif //KGL_KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H
