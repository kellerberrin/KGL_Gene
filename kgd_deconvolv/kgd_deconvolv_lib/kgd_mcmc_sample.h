//
// Created by kellerberrin on 17/05/18.
//

#ifndef KGD_MCMC_SAMPLE_H
#define KGD_MCMC_SAMPLE_H

#include <vector>



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class McmcSample {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif

public:

  McmcSample() = default;
  ~McmcSample() = default;

  void clear() {

    proportion.clear();
    sumLLKs.clear();
    moves.clear();

  }

  const std::vector<double>& getSiteOfTwoSwitchOne() const { return siteOfTwoSwitchOne; }
  const std::vector<double>& getSiteOfTwoMissCopyOne() const { return siteOfTwoMissCopyOne; }
  const std::vector<double>& getSiteOfTwoSwitchTwo() const { return siteOfTwoSwitchTwo; }
  const std::vector<double>& getSiteOfTwoMissCopyTwo() const { return siteOfTwoMissCopyTwo; }
  const std::vector<double>& getSiteOfOneSwitchOne() const { return siteOfOneSwitchOne; }
  const std::vector<double>& getSiteOfOneMissCopyOne() const { return siteOfOneMissCopyOne; }

  const std::vector<double>& getCurrentsiteOfTwoSwitchOne() const { return currentsiteOfTwoSwitchOne; }
  const std::vector<double>& getCurrentsiteOfTwoMissCopyOne() const { return currentsiteOfTwoMissCopyOne; }
  const std::vector<double>& getCurrentsiteOfTwoSwitchTwo() const { return currentsiteOfTwoSwitchTwo; }
  const std::vector<double>& getCurrentsiteOfTwoMissCopyTwo() const { return currentsiteOfTwoMissCopyTwo; }
  const std::vector<double>& getCurrentsiteOfOneSwitchOne() const { return currentsiteOfOneSwitchOne; }
  const std::vector<double>& getCurrentsiteOfOneMissCopyOne() const { return currentsiteOfOneMissCopyOne; }

  const std::vector<std::vector<double> >& getProportion() const { return proportion; }
  const std::vector<double>& getProportionIndex(size_t index) const { return proportion[index]; }
  double getProportionIndex(size_t i, size_t j) const { return proportion[i][j]; }

  const std::vector<double>& getSumLLKs() const { return sumLLKs; }
  double getSumLLKsIndex(size_t index) const { return sumLLKs[index]; }

  int getMovesIndex(size_t index) const { return moves[index]; }

  const std::vector<std::vector<double> >& getHap() const { return hap; }
  const std::vector<double>& getHapIndex(size_t index) const { return hap[index]; }
  double getHapIndex(size_t i, size_t j) const { return hap[i][j]; }


  void sumSiteOfTwoSwitchOne(size_t index, double value) { siteOfTwoSwitchOne[index] += value; }
  void sumSiteOfTwoMissCopyOne(size_t index, double value) { siteOfTwoMissCopyOne[index] += value; }
  void sumSiteOfTwoSwitchTwo(size_t index, double value) { siteOfTwoSwitchTwo[index] += value; }
  void sumSiteOfTwoMissCopyTwo(size_t index, double value) { siteOfTwoMissCopyTwo[index] += value; }
  void sumSiteOfOneSwitchOne(size_t index, double value) { siteOfOneSwitchOne[index] += value; }
  void sumSiteOfOneMissCopyOne(size_t index, double value) { siteOfOneMissCopyOne[index] += value; }

  void setCurrentsiteOfTwoSwitchOne(size_t index, double value) { currentsiteOfTwoSwitchOne[index] = value; }
  void setCurrentsiteOfTwoMissCopyOne(size_t index, double value) { currentsiteOfTwoMissCopyOne[index] = value; }
  void setCurrentsiteOfTwoSwitchTwo(size_t index, double value) { currentsiteOfTwoSwitchTwo[index] = value; }
  void setCurrentsiteOfTwoMissCopyTwo(size_t index, double value) { currentsiteOfTwoMissCopyTwo[index] = value; }
  void setCurrentsiteOfOneSwitchOne(size_t index, double value) { currentsiteOfOneSwitchOne[index] = value; }
  void setCurrentsiteOfOneMissCopyOne(size_t index, double value) { currentsiteOfOneMissCopyOne[index] = value; }

  void setHap(const std::vector<std::vector<double> >& newhap) { hap = newhap; }

  void setVectorSize(size_t vector_size) {

    siteOfTwoSwitchOne.resize(vector_size, 0);
    siteOfTwoMissCopyOne.resize(vector_size, 0);
    siteOfTwoSwitchTwo.resize(vector_size, 0);
    siteOfTwoMissCopyTwo.resize(vector_size, 0);
    siteOfOneSwitchOne.resize(vector_size, 0);
    siteOfOneMissCopyOne.resize(vector_size, 0);

    currentsiteOfTwoSwitchOne.resize(vector_size, 0);
    currentsiteOfTwoMissCopyOne.resize(vector_size, 0);
    currentsiteOfTwoSwitchTwo.resize(vector_size, 0);
    currentsiteOfTwoMissCopyTwo.resize(vector_size, 0);
    currentsiteOfOneSwitchOne.resize(vector_size, 0);
    currentsiteOfOneMissCopyOne.resize(vector_size, 0);

  }

  void addMove(int value) { moves.push_back(value); }
  void addProportion(const std::vector<double>& prop_vec) { proportion.push_back(prop_vec); }
  void addSumLLKs(double llk_value) { sumLLKs.push_back(llk_value); }

  void divideSiteVectors(double divisor) {

    for (size_t atSiteI = 0; atSiteI < siteOfTwoSwitchOne.size(); ++atSiteI) {

      siteOfTwoSwitchOne[atSiteI] /= divisor;
      siteOfTwoMissCopyOne[atSiteI] /= divisor;
      siteOfTwoSwitchTwo[atSiteI] /= divisor;
      siteOfTwoMissCopyTwo[atSiteI] /= divisor;
      siteOfOneSwitchOne[atSiteI] /= divisor;
      siteOfOneMissCopyOne[atSiteI] /= divisor;

    }

  }


private:

  std::vector<double> siteOfTwoSwitchOne;
  std::vector<double> siteOfTwoMissCopyOne;
  std::vector<double> siteOfTwoSwitchTwo;
  std::vector<double> siteOfTwoMissCopyTwo;
  std::vector<double> siteOfOneSwitchOne;
  std::vector<double> siteOfOneMissCopyOne;

  std::vector<double> currentsiteOfTwoSwitchOne;
  std::vector<double> currentsiteOfTwoMissCopyOne;
  std::vector<double> currentsiteOfTwoSwitchTwo;
  std::vector<double> currentsiteOfTwoMissCopyTwo;
  std::vector<double> currentsiteOfOneSwitchOne;
  std::vector<double> currentsiteOfOneMissCopyOne;

  std::vector<std::vector<double> > proportion;
  std::vector<std::vector<double> > hap;
  std::vector<double> sumLLKs;
  std::vector<int> moves;

};



}   // organization level namespace
}   // project level namespace





#endif //KGD_MCMC_SAMPLE_H
