//
// Created by kellerberrin on 19/05/18.
//

#ifndef KGD_ALLELEFREQ_H
#define KGD_ALLELEFREQ_H

#include <vector>
#include <map>
#include <algorithm>
#include "kgd_deconvolv_app.h"
#include "kgd_deconvolv_types.h"



namespace kellerberrin::deconvolv {    // organization level namespace


/// Single allele frequency paired with the contig offset.
class AlleleOffset {

public:

  AlleleOffset(ContigOffset_t offset, AlleleFreq_t allelefreq) : offset_(offset), allelefreq_(allelefreq) {}
  AlleleOffset(const AlleleOffset&) = default;
  ~AlleleOffset() = default;

  ContigOffset_t offset() const { return offset_; }
  AlleleFreq_t alleleFreq() const { return allelefreq_; }

private:

  ContigOffset_t offset_;
  AlleleFreq_t allelefreq_;

};


/// A vector of allele frequencies for a contig
using AlleleVector = std::vector<AlleleOffset>;
class ContigAlleles {

public:


  ContigAlleles(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  ContigAlleles(const ContigAlleles&) = default;
  ~ContigAlleles() = default;

  // Access routines
  const ContigId_t& contigId() const { return contig_id_; }
  void addAlleleFreq(const AlleleOffset& allele_offset) { allele_vector_.push_back(allele_offset); }
  const AlleleVector& getAlleleVector() const { return allele_vector_; }

  // Sort the alleles in increasing offset order.
  void sortAlleleFreq() {  std::sort(allele_vector_.begin(), allele_vector_.end(), [](AlleleOffset& a, AlleleOffset& b) { return a.offset() > b.offset(); }); }

private:

  ContigId_t contig_id_;
  AlleleVector allele_vector_;

};


/// A map of contig allele frequencies for a genome.
using ContigAlleleMap = std::map<ContigId_t, ContigAlleles>;
class GenomeAlleles {

public:

  explicit GenomeAlleles(const GenomeId_t& genome_id) : genome_id_(genome_id) {}
  GenomeAlleles() : genome_id_(GENOME_NOT_SPECIFIED_) {}
  GenomeAlleles(const GenomeAlleles&) = default;
  ~GenomeAlleles() = default;

  // Access routines
  [[nodiscard]] const GenomeId_t& genomeId() const { return genome_id_; }
  [[nodiscard]] const ContigAlleleMap& getMap() const { return contig_allele_map_; }

  void addAlleleFreq(const ContigId_t& contig_id, ContigOffset_t offset, AlleleFreq_t frequency);

  // Methods
  void sortAlleleFreq() { for(auto contig : contig_allele_map_)  contig.second.sortAlleleFreq(); }

private:

  GenomeId_t genome_id_;
  ContigAlleleMap contig_allele_map_;

  constexpr static const char* GENOME_NOT_SPECIFIED_ = "Genome Id Not Specified";

  // Get or Create ContigAlleles
  ContigAlleles& getCreateContigAllele(const ContigId_t& contig_id);

};



}   // end amespace


#endif //KGD_ALLELEFREQ_H
