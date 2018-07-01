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
  ~MixtureDataObj() = default;

  MixtureDataObj& operator=(const MixtureDataObj& copy) = default;

  // Access functions
  const std::vector<size_t>& indexOfChromStarts() const { return indexOfChromStarts_; }
  const std::vector<std::vector<int>>& getPosition() const { return position_; }
  const std::vector<double>& getPlaf() const { return plaf_; }
  const std::vector<double>& getRefCount() const { return refCount_; }
  const std::vector<double>& getAltCount() const { return altCount_; }
  std::vector<std::string> getChrom() const { return chrom_; }
  size_t nLoci() const { return nLoci_; }

  // Read functions
  bool readVCFPlaf(const std::string& filename,
                   const std::string& plaf_filename);  // Either read this
  bool readVCFPlafExclude(const std::string& vcf_filename,
                          const std::string& plaf_filename,
                          const std::string& exclude_filename);  // Either read this
  bool readRefAltPlaf(const std::string& ref_filename,
                      const std::string& alt_filename,
                      const std::string& plaf_filename);
  bool readRefAltPlafExclude(const std::string& ref_filename,
                             const std::string& alt_filename,
                             const std::string& plaf_filename,
                             const std::string& exclude_filename);

private:


  std::vector<std::string> chrom_;
  std::vector<size_t> indexOfChromStarts_;
  std::vector<std::vector<int> > position_;
  std::vector<double> plaf_;
  std::vector<double> refCount_;
  std::vector<double> altCount_;
  size_t nLoci_;

};



}   // organization level namespace
}   // project level namespace


#endif //KGD_CTL_DATA_H
