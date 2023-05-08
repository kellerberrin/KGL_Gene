//
// Created by kellerberrin on 28/02/18.
//

#ifndef KGL_VARIANT_FACTORY_VCF_CIGAR_H
#define KGL_VARIANT_FACTORY_VCF_CIGAR_H



#include <map>

#include "kel_utility.h"
#include "kgl_genome_types.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF parser. Miscellaneous parser functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returned from the cigar functions.
enum class CigarEditType : char { UNCHANGED = 'M', INSERT = 'I', DELETE = 'D', CHANGED = 'X'};
using CigarEditItem = std::pair<size_t, CigarEditType>; // Used to specify edit as vector '1M,1X,3D,3I'.
using CigarVector = std::vector<CigarEditItem>;


class ParseVCFCigar {

public:

  ParseVCFCigar() = default;
  ~ParseVCFCigar() = default;


  [[nodiscard]] static bool parseCigar( const std::string& cigar,
                                        size_t& check_reference_size,
                                        size_t& check_alternate_size,
                                        std::vector<CigarEditItem>& parsed_cigar);


  // The CigarVector contains (n x 'M') + (1 x 'X') and is an SNP.
  [[nodiscard]] static bool isSNP(const std::string& reference, const std::string& alternate);

  // Generate a CIGAR from two sequences.
  [[nodiscard]] static std::string generateCigar(const std::string& reference, const std::string& alternate);

  // Generate a CIGAR from a cigar vector
  [[nodiscard]] static std::string generateCigar(const CigarVector& cigar_vector);

// Use edlib to generate a cigar vector.
  [[nodiscard]] static CigarVector generateEditVector(const std::string& reference, const std::string& alternate);

// Calculate the allele offset, the non-zero value of the first CigarVector item if it is 'UNCHANGED'
  [[nodiscard]] static ContigOffset_t alleleEditOffset(const std::string& reference, const std::string& alternate);
  // Calculate the allele offset, the non-zero value of the first CigarVector item if it is 'UNCHANGED'
  [[nodiscard]] static ContigOffset_t alleleOffset(const std::string& reference, const std::string& alternate);


  // Given a reference count and a cigar vector compute a number that calculates the equivalent
  // size of the alternate string.
  // For UNCHANGED = 'M' and CHANGED = 'X' cigar items the reference_count and alternate count are incremented.
  // For INSERT = 'I' the alternate is incremented and the reference_count is not.
  // For DELETE = 'D' the reference count is incremented and the alternate is not.
  [[nodiscard]] static size_t alternateCount(size_t reference_count, const CigarVector& cigar_vector);


private:


// Use edlib to generate a cigar string.
  static std::vector<CigarEditType> generateEditString(const std::string& reference, const std::string& alternate);

};



}   // end namespace








#endif //KGL_VARIANT_FACTORY_VCF_CIGAR_H
