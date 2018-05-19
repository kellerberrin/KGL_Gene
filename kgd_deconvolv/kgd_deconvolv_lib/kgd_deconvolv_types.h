//
// Created by kellerberrin on 18/05/18.
//

#ifndef KGD_DECONVOLV_TYPES_H
#define KGD_DECONVOLV_TYPES_H

#include <string>



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


using GenomeId_t = std::string;
using ContigId_t = std::string;
using AlleleFreq_t = double;
using ContigOffset_t = uint64_t;         // Paris Japonica has 150 billion base pairs; use 64 bit integers.


}   // organization level namespace
}   // project level namespace


#endif //KGD_DECONVOLV_TYPES_H
