//
// Created by kellerberrin on 17/07/23.
//

#ifndef KGL_VARIANT_FILTER_UNIQUE_H
#define KGL_VARIANT_FILTER_UNIQUE_H

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// For a haploid organism, such as the blood stage of P.Falciparum, only one variant can occur at a particular offset.
// These filters attempt to resolve the situation where there is more than 1 valid variant specified at a location.
// The filters assume that all variants have been converted to "canonical" format where the cigar format of an SNP is '1X',
// Deletes are '1MnD' and Inserts are '1MnI'.
//
// Note that an SNP and Indel specifying the same offset is not a problem, since  by convention, canonical
// Indels actually occur at the next (+1) offset..
////////////////////////////////////////////////////////////////////////////////////////////////////////////////




#endif //KGL_VARIANT_FILTER_UNIQUE_H
