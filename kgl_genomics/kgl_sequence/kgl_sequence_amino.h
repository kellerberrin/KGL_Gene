//
// kgl_sequence_amino.h — legacy include shim.  [DEPRECATED — scheduled for removal]
//
// In the refactored module the amino sequence type (`AminoSequence`) lives in `kgl_sequence.h`
// and the translation facade (`TranslateToAmino`) in `kgl_genetic_code.h`. This header keeps
// the reference include path working for one release. Migrate with:
//     #include "kgl_sequence_amino.h"  ->  #include "kgl_sequence.h" + #include "kgl_genetic_code.h"
// Delete this shim after the downstream migration lands.
//

#ifndef KGL_SEQUENCE_AMINO_H
#define KGL_SEQUENCE_AMINO_H


#include "kgl_sequence.h"
#include "kgl_genetic_code.h"


#endif //KGL_SEQUENCE_AMINO_H
