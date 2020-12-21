

#ifndef KGD_TXTREADER_H
#define KGD_TXTREADER_H


#include "kgd_alleleread.h"
#include "kgd_variant_index.h"
#include "kgd_exceptions.h"


namespace kellerberrin::deconvolv {          // project level namespace


class TxtReader : public VariantIndex {

#ifdef UNITTEST
  friend class TestPanel;
  friend class TestTxtReader;
  friend class TestInitialHaplotypes;
#endif


public: // move the following to private

  TxtReader() = default;
  ~TxtReader() override = default;

  // Methods
  virtual void readFromFile(const char inchar[]) { this->readFromFileBase(inchar); };
  void readFromFileBase(const char inchar[]);
  void removeMarkers();

  // Access routines.
  double getContentIndex(size_t loci, size_t strain) const { return content_[loci][strain]; }
  const std::vector<std::vector<double> >& getContent() const { return content_; }
  size_t getInfoLines() const { return nInfoLines_; }
  const std::vector<double>& getInfo() const { return info_; }

  void setInfoLines(size_t lines) { nInfoLines_ = lines; }
  void addContent(const std::vector<double>& content) { content_.push_back(content); }
  std::vector<std::vector<double> >& setContent() { return content_; }

private:

  // content is a matrix of n.loci by n.strains, i.e. content length is n.loci
  std::vector<std::vector<double> > content_;
  std::vector<std::vector<double> > keptContent_;
  // info_ only refers to the first column of the content
  std::vector<double> info_;

  std::vector<size_t> tmpPosition_;

  size_t nInfoLines_;
  int tmpChromInex_;

  std::string fileName_;

  AlleleReader allele_reader_;

  // Methods
  void extractChrom(std::string &tmp_str);
  void extractPOS(std::string &tmp_str);
  void reshapeContentToInfo();


};


class ExcludeMarker : public TxtReader {

  // sorting

public:

  ExcludeMarker() = default;
  ~ExcludeMarker() override = default;

};


}   // organization level namespace


#endif
