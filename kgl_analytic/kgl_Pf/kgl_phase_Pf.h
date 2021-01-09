//
// Created by kellerberrin on 9/1/21.
//

#ifndef KGL_PHASE_PF_H
#define KGL_PHASE_PF_H

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Accepts a Diploid Falciparum Variant database and attempts to return a Haploid variant database by examining
// the Complexity of Infection file (not available for all samples), if the variant is heterozygous or homozygous
// (the blood phase of Pf is Haploid) and VCF FORMAT statistics.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PhaseFalciparum {

public:

  PhaseFalciparum() = default;
  ~PhaseFalciparum() = default;



private:



};


#endif //KGL_PHASE_PF_H
