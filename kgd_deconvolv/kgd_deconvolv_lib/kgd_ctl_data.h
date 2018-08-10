//
// Created by kellerberrin on 30/06/18.
//

#ifndef KGD_CTL_DATA_H
#define KGD_CTL_DATA_H



#include <memory>
#include <vector>


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


class MixtureDataObj {

public:

  MixtureDataObj() = default;
  MixtureDataObj(  std::vector<std::string>& chrom,
                   std::vector<size_t>& indexOfChromStarts,
                   std::vector<std::vector<size_t> >& position,
                   std::vector<double>& plaf,
                   std::vector<double>& refCount,
                   std::vector<double>& altCount) : chrom_(chrom), indexOfChromStarts_(indexOfChromStarts), position_(position),
                                                    plaf_(plaf), refCount_(refCount), altCount_(altCount) {}
  ~MixtureDataObj() = default;

  MixtureDataObj& operator=(const MixtureDataObj& copy) = default;

  // Access functions
  const std::vector<size_t>& indexOfChromStarts() const { return indexOfChromStarts_; }
  const std::vector<std::vector<size_t>>& getPosition() const { return position_; }
  const std::vector<double>& getPlaf() const { return plaf_; }
  const std::vector<double>& getRefCount() const { return refCount_; }
  const std::vector<double>& getAltCount() const { return altCount_; }
  std::vector<std::string> getChrom() const { return chrom_; }
  size_t nLoci() const { return plaf_.size(); }

  // Read functions
  bool readVCFPlaf(const std::string& filename,
                   const std::string& plaf_filename);
  bool readVCFPlafExclude(const std::string& vcf_filename,
                          const std::string& plaf_filename,
                          const std::string& exclude_filename);
  bool readRefAltPlaf(const std::string& ref_filename,
                      const std::string& alt_filename,
                      const std::string& plaf_filename);
  bool readRefAltPlafExclude(const std::string& ref_filename,
                             const std::string& alt_filename,
                             const std::string& plaf_filename,
                             const std::string& exclude_filename);

// Verifies the data structure and (optional) prints out contigs and variant stats.
  bool verifyPrint(bool print = false /* check only */) const;

private:


  std::vector<std::string> chrom_;
  std::vector<size_t> indexOfChromStarts_;
  std::vector<std::vector<size_t> > position_;
  std::vector<double> plaf_;
  std::vector<double> refCount_;
  std::vector<double> altCount_;

  double calcMedianCount(std::vector<double>& countVector) const;

};



}   // organization level namespace
}   // project level namespace


#endif //KGD_CTL_DATA_H
