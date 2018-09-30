//
// Created by kellerberrin on 29/09/18.
//

#ifndef KGL_FINESTRUCTURE_ANALYSIS_H
#define KGL_FINESTRUCTURE_ANALYSIS_H

#include <string>
#include <kgl_variant_db.h>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



// This class generates the input files for the FineStructure population analysis
// software. "Inference of population structure using dense haplotype data, Daniel Lawson, Garrett Hellenthal,
// Simon Myers, and Daniel Falush, 2012. PLoS Genetics, Vol. 8(1): e1002453,"

// Input is a pointer to a phased population "std::shared_ptr<const PhasedPopulation> population_ptr"

// Three files are generated:
// 1. The phase file "<Filename>.phase"
// 2. The recombination file "<Filename>.recomb"
// 3. The ID file "<Filename>.ids"


//##########################
//IDFILE FORMAT:
//This specifies the names of the individuals in the data, as well as (optionally)
//which population they are from and
//whether they are included.
//Format: N lines, one per individual, containing the following columns:
//<NAME> <POPULATION> <INCLUSION> <ignored extra info>
//Where <NAME> and <POPULATION> are strings and <INCLUSION> is 1 to include an
//individual and 0 to exclude them. The
//second and third columns can be omitted (but the second must be present if the
//third is). Currently <POPULATION> is not
//used by this version of fs.
//EXAMPLE IDFILE:
//Ind1 Pop1 1
//Ind2 Pop1 1
//Ind3 Pop2 0
//Ind4 Pop2 1
//Ind5 Pop2 1

// THe phase file (one for each chromosome)
//##########################
//CHROMOPAINTER’S v2 ’PHASE’ FORMAT:
//This is heavily based on ’FastPhase’ output.
//* The first line contains the number of *haplotypes* (i.e. for diploids, 2* the
//number of individuals).
//The
//second line contains the number of SNPs.
//*
//The
//third line contains the letter P, followed by the basepair location of each
//*
//SNP (space separated). These must
//match the recombination file. Within each chromosome, basepairs must be in order.
//* Each additional line contains a haplotype, in the order specified in the IDFILE.
//Diploids have two contiguous rows.
//Each character (allowing no spaces!) represents a *biallelic* SNP. Accepted
//characters are 0,1,A,C,G,T, with NO missing
//values!
//EXAMPLE PHASEFILE:
//610
//6
//P 100 200 300 400 500 600
//010101
//011101
//111101
//001101
//011000
//001100
//001001
//001011
//001001
//001111


//##########################
//CHROMOPAINTERS RECOMBINATION FILE FORMAT:
//Required only if running in unlinked mode.
//This specifies the distance between SNPs in ’recombination rate’ units. There
//should be a header line followed by one
//line for each SNP in haplotype infile. Each line should contain two columns, with
//the first column denoting the
//basepair position values given in haplotype infile, in the same order. The second
//column should give the genetic
//distance per basepair between the SNP at the position in the first column of the
//same row and the SNP at the position
//in the first column of the subsequent row. The last row should have a ’0’ in the
//second column (though this is not
//required
//this value is simply ignored by the program). Genetic distance should
// be given in Morgans, or at least the
// relevant output files assume this value is in Morgans.
// If you are including genetic information from multiple chromosomes, put a ’-9’ (or
// any value < 0) next to the last
// basepair position of the preceeding chromosome.
// EXAMPLE RECOMBFILE:
// start.pos recom.rate.perbp
// 100 0.01
// 200 0.02
// 300 -9
// 400 0.02
// 500 0.05
// 600 0


class FineStructureAnalysis {

public:

  FineStructureAnalysis() = default;
  ~FineStructureAnalysis() = default;

  static bool generateFiles(const std::string& Filename,
                            std::shared_ptr<const PhasedPopulation> population_ptr,
                            double bases_per_centimorgan);


private:


  constexpr static const char DELIMITER_ = ' ';

  static bool generateIDFile(const std::string& Filename,
                             std::shared_ptr<const PhasedPopulation> population_ptr);

  static bool generatePhaseFile(const std::string& Filename,
                                std::shared_ptr<const PhasedPopulation> population_ptr,
                                double bases_per_centimorgan);

};





}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_FINESTRUCTURE_ANALYSIS_H
