/*
 * kgd_deploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgd_deploid.
 *
 * kgd_deploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <kgl_exec_env.h>
#include "kgd_ibdconfig.h"


namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;



void kgd::IBDconfiguration::buildIBDconfiguration(size_t k) {

  setKstrain(k);
  enumerateOp();
  makePairList();
  makePairToEmission();
  findUniqueState();
  findEffectiveK();

}



void kgd::IBDconfiguration::enumerateOp() {
  //#For each configuration, identify which pairs are IBD

  op = enumerateBinaryMatrixOfK(nchoose2(kStrain()));

}


void kgd::IBDconfiguration::makePairList() {
  //#Make a map of pairs to pair value

  assert(pairList.size() == 0);

  for (size_t i = 0; i < kStrain(); i++) { // 0-indexed

    for (size_t j = i + 1; j < kStrain(); j++) {

      pairList.push_back(std::vector<size_t>({i, j}));

    }

  }

  assert(pairList.size() == (size_t) nchoose2(kStrain()));

}


void kgd::IBDconfiguration::makePairToEmission() {

  assert(pairToEmission.size() == 0);
  //prs2ems<-array(0, c(nrow(op), k));
  //for ( size_t i = 0; i < op.size(); i++ ){

  for (std::vector<int> tmpOp : op) {

    std::vector<int> tmpRow = makeTmpRow();
    //for ( size_t i = 0 ; i <(*opIt).size(); i++){
    //cout << (*opIt)[i] << " ";
    //}
    //cout<<endl;

    std::vector<size_t> ii = findWhichIsOne(tmpOp);
    //ii <- which(op[rowI,]==1);
    //cout << ii.size()<<endl;
    //cout << "##############" <<endl;

    if (ii.size() > 0) {

      std::vector<std::vector<size_t> > tmpIBDPairs;

      for (size_t j = 0; j < ii.size(); j++) {
        //cout << "j = " << j <<" ii[j] = "<<ii[j]<<endl;
        tmpIBDPairs.push_back(pairList[ii[j]]);

        //cout << tmpIBDPairs.back()[0]<< " "<<tmpIBDPairs.back()[1]<<endl;
      }

      int tmpIndex = (tmpIBDPairs.size() - 1);

      while (tmpIndex >= 0) {

        //cout << "replacing element "<< tmpIBDPairs[tmpIndex][0] << " by element " << tmpIBDPairs[tmpIndex][1] <<endl;
        tmpRow[tmpIBDPairs[tmpIndex][0]] = tmpRow[tmpIBDPairs[tmpIndex][1]];
        tmpIndex--;

      }

    }

    pairToEmission.push_back(tmpRow);
    //for (size_t i = 0 ; i < pairToEmission.back().size(); i++){
    //cout << pairToEmission.back()[i]<<" ";
    //}
    //cout <<endl;
  }
  //for (size_t i = 0; i < pairToEmission.size(); i++){
  //for (size_t ii = 0 ; ii < pairToEmission.back().size(); ii++){
  //cout << pairToEmission[i][ii]<<" ";
  //}
  //cout <<endl;
  //}

}


std::vector<int> kgd::IBDconfiguration::makeTmpRow() {

  std::vector<int> ret(kStrain());

  for (size_t i = 0; i < ret.size(); i++) {
    ret[i] = (int) i;
  }

  return ret;

}


std::vector<size_t> kgd::IBDconfiguration::findWhichIsOne(std::vector<int> tmpOp) {

  std::vector<size_t> ret;

  for (size_t i = 0; i < tmpOp.size(); i++) {

    if (tmpOp[i] == 1) {

      ret.push_back(i);

    }

  }

  return ret;

}


void kgd::IBDconfiguration::findUniqueState() {

  kgl::ExecEnv::log().info("pairToEmission size:{}, states size: {}", pairToEmission.size(), states.size());
  assert (states.size() == 0);
  //states.push_back(pairToEmission[0]);
  //for (size_t i = 1; i < pairToEmission.size(); i++){
  //bool aNewState = true;
  //for ( vector<int> state : states){
  //if ( twoVectorsAreSame(state, pairToEmission[i]) ){
  //aNewState = false;
  //break;
  //}
  //}
  //if ( aNewState ){
  //states.push_back(pairToEmission[i]);
  //}
  //}

  states = unique(pairToEmission);

  //for ( vector<int> state : states){
  //for (int i : state){
  //cout << i <<" ";
  //}
  //cout <<endl;
  //}

}


void kgd::IBDconfiguration::findEffectiveK() {

  assert(effectiveK.size() == 0);

  for (std::vector<int> state : states) {

    std::set<int> tmpSet(state.begin(), state.end());
    //cout << tmpSet.size() <<endl;
    effectiveK.push_back(tmpSet.size());

  }

  assert(effectiveK.size() == states.size());

}


std::vector<std::string> kgd::IBDconfiguration::getIBDconfigureHeader() {

  std::vector<std::string> ret;

  for (size_t i = 0; i < states.size(); i++) {

    std::string tmp;

    for (size_t j = 0; j < states[i].size(); j++) {

      std::stringstream tmp_ss;
      tmp_ss << states[i][j];
      tmp += tmp_ss.str() + ((j < (states[i].size() - 1)) ? "-" : "");

    }

    ret.push_back(tmp);

  }

  return ret;

}


std::vector<std::vector<int> > kgd::IBDconfiguration::unique(std::vector<std::vector<int> > &mat) {

  std::vector<std::vector<int> > ret;

  ret.push_back(mat[0]);

  for (size_t i = 1; i < mat.size(); i++) {

    bool aNewState = true;

    for (std::vector<int> state : ret) {

      if (twoVectorsAreSame(state, mat[i])) {
        aNewState = false;
        break;
      }

    }

    if (aNewState) {
      ret.push_back(mat[i]);
    }

  }

  return ret;

}


std::vector<std::vector<int> > kgd::IBDconfiguration::enumerateBinaryMatrixOfK(size_t k) {
  // This function enumerate all possible binary combinations of k elements

  int ksq = pow(2, k);

  std::vector<std::vector<int> > ret;

  for (int i = 0; i < ksq; i++) {

    ret.push_back(convertIntToBinary(i, k));

  }

  return ret;

}


std::vector<int> kgd::IBDconfiguration::convertIntToBinary(int x, size_t len) {

  std::vector<int> ret(len);

  size_t idx = 0;

  while (x) {

    ret[idx] = (x & 1) ? 1 : 0;
    idx++;
    //cout << "x " <<x<< " idx "<<idx<<" len "<< len<<endl;

    if (idx > len) {

      throw OutOfVectorSize();
    }

    x >>= 1;

  }

  std::reverse(ret.begin(), ret.end());
  //for (size_t i = 0; i < ret.size(); i++){
  //cout << ret[i] << " ";
  //}
  //cout<<endl;
  return ret;
}


int kgd::IBDconfiguration::nchoose2(int n) {

  if (n < 2) {

    throw InvalidInput("Input must be at least 2!");

  }

  int ret = n * (n - 1) / 2;
  return ret;

}


bool kgd::IBDconfiguration::twoVectorsAreSame(std::vector<int> vec1, std::vector<int> vec2) {

  if (vec1.size() != vec2.size()) {

    throw InvalidInput("Input vectors have different length!");

  }

  bool ret = true;

  for (size_t i = 0; i < vec1.size(); i++) {

    if (vec1[i] != vec2[i]) {

      ret = false;
      break;

    }

  }

  return ret;

}

