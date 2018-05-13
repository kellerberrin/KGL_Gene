//
// Created by kellerberrin on 13/05/18.
//


#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <kgl_exec_env.h>
#include "kgd_hprior.h"



namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;


void kgd::Hprior::buildHprior(size_t kStrain, std::vector<double> &plaf) {

  ibdConfig.buildIBDconfiguration(kStrain);

  effectiveK = ibdConfig.effectiveK;
  nState_ = 0;
  plaf_ = plaf;
  setKstrain(kStrain);
  setnLoci(plaf_.size());

  std::vector<std::vector<int> > hSetBase = IBDconfiguration::enumerateBinaryMatrixOfK(getkStrain());
  size_t stateI = 0;

  for (std::vector<int> state : ibdConfig.states) {

    std::set<int> stateUnique(state.begin(), state.end());

    assert(stateUnique.size() == effectiveK[stateI]);

    std::vector<std::vector<int> > hSetBaseTmp = hSetBase;

    for (size_t j = 0; j < (size_t) getkStrain(); j++) {

      for (size_t i = 0; i < hSetBase.size(); i++) {

        hSetBaseTmp[i][j] = hSetBase[i][state[j]];

      }

    }

    std::vector<std::vector<int> > hSetBaseTmpUnique = IBDconfiguration::unique(hSetBaseTmp); // uu

    size_t sizeOfhSetBaseTmpUnique = hSetBaseTmpUnique.size();

    stateIdxFreq.push_back(sizeOfhSetBaseTmpUnique);
    //h.prior.i<-array(0, c(size.h.set.i, n.loci));

    for (size_t i = 0; i < sizeOfhSetBaseTmpUnique; i++) {
      //vector<int> hSetBaseTmpUniqueSubSet(); // uu[i,a.u,drop=F]
      int tmpSum = 0;

      for (int uniqSt : stateUnique) {

        tmpSum += hSetBaseTmpUnique[i][uniqSt];

      }
      // sumOfVec(hSetBaseTmpUnique[i]);
      int tmpDiff = stateUnique.size() - tmpSum;
      //cout << stateUnique.size() << " " << tmpSum << " " << tmpDiff<< endl;
      std::vector<double> hPriorTmp(nLoci());

      for (size_t site = 0; site < nLoci(); site++) {
        //cout << (1.0-plaf_[site]) << " " << tmpDiff << " " <<pow((1.0-plaf_[site]),(double)tmpDiff) << endl;

        hPriorTmp[site] = pow(plaf_[site], (double) tmpSum) * pow((1.0 - plaf_[site]), (double) tmpDiff);

      }

      priorProb.push_back(hPriorTmp);
      hSet.push_back(hSetBaseTmpUnique[i]);

      nState_++;
      stateIdx.push_back(stateI);

    }

    stateI++;

  }

}


void kgd::Hprior::transposePriorProbs() {

  assert(priorProbTrans.size() == 0);

  for (size_t i = 0; i < nLoci(); i++) {

    std::vector<double> priorProbTransTmp(nState());

    for (size_t j = 0; j < nState(); j++) {

      priorProbTransTmp[j] = priorProb[j][i];
      //cout << priorProb[j][i] << " ";
    }

    priorProbTrans.push_back(priorProbTransTmp);
    //cout << endl;
  }

}


std::vector<std::string> kgd::Hprior::getIBDconfigureHeader() {

  return ibdConfig.getIBDconfigureHeader();

}

