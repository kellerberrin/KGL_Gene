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

#include "kgd_dEploidIO.h"
#include "kgd_updateHap.h"
#include "kgd_mcmc.h"


namespace kgd = kellerberrin::deploid;


void kgd::McmcMachinery::writeLastFwdProb(bool useIBD) {

  if (dEploidIO_->doExportPostProb() != true) {

    return;

  }

  for (size_t tmpk = 0; tmpk < kStrain_; tmpk++) {

    if (dEploidIO_->doAllowInbreeding() == true) {

      updateReferencePanel(panel_->truePanelSize() + kStrain_ - 1, tmpk);

    }

    for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts_.size(); chromi++) {

      size_t start = dEploidIO_->indexOfChromStarts_[chromi];
      size_t length = dEploidIO_->position_[chromi].size();

      UpdateSingleHap updatingSingle(dEploidIO_->refCount_,
                                     dEploidIO_->altCount_,
                                     dEploidIO_->plaf_,
                                     currentExpectedWsaf_,
                                     currentProp_,
                                     currentHap_,
                                     hapRg_,
                                     start,
                                     length,
                                     panel_,
                                     dEploidIO_->missCopyProb_,
                                     dEploidIO_->scalingFactor(),
                                     tmpk);

      if (dEploidIO_->doAllowInbreeding() == true) {

        updatingSingle.setPanelSize(panel_->inbreedingPanelSize());

      }

      updatingSingle.core(dEploidIO_->refCount_,
                          dEploidIO_->altCount_,
                          dEploidIO_->plaf_,
                          currentExpectedWsaf_,
                          currentProp_,
                          currentHap_);

      dEploidIO_->writeLastSingleFwdProb(updatingSingle.fwdProbs_, chromi, tmpk, useIBD);

    }
    //UpdatePairHap updating( dEploidIO_->refCount_,
    //dEploidIO_->altCount_,
    //dEploidIO_->plaf_,
    //currentExpectedWsaf_,
    //currentProp_, currentHap_, hapRg_,
    //start, length,
    //panel_, dEploidIO_->missCopyProb_, dEploidIO_->forbidCopyFromSame(),
    //(size_t)0,
    //(size_t)1);
    //updating.core ( dEploidIO_->refCount_, dEploidIO_->altCount_, dEploidIO_->plaf_, currentExpectedWsaf_, currentProp_, currentHap_);
    //dEploidIO_->writeLastPairFwdProb( updating, chromi );
  }
}


void kgd::DEploidIO::writeLastSingleFwdProb(std::vector<std::vector<double> > &probabilities,
                                            size_t chromIndex,
                                            size_t strainIndex,
                                            bool useIBD) {

  if (probabilities.size() == 0) {

    return;

  }

  size_t panelSize = probabilities[0].size();

  std::string strExportFwdProb = ((useIBD == true) ? strIbdExportSingleFwdProbPrefix_ : strExportSingleFwdProbPrefix_) + std::to_string(strainIndex);

  ofstreamExportFwdProb_.open(strExportFwdProb.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  if (chromIndex == 0) { // Print header

    ofstreamExportFwdProb_ << "CHROM" << "\t" << "POS" << "\t";

    for (size_t ii = 0; ii < probabilities[0].size(); ii++) {

      if (doAllowInbreeding() == true) {

        if (ii <= (panelSize - kStrain())) {

          ofstreamExportFwdProb_ << "P" << (ii + 1);

        } else {

          ofstreamExportFwdProb_ << "I" << (ii) - (panelSize - kStrain());

        }

      } else {

        ofstreamExportFwdProb_ << (ii + 1);

      }

      ofstreamExportFwdProb_ << ((ii < (panelSize - 1)) ? "\t" : "\n");

    }

  }

  size_t siteIndex = 0;

  for (size_t posI = 0; posI < position_[chromIndex].size(); posI++) {

    ofstreamExportFwdProb_ << chrom_[chromIndex] << "\t" << (int) position_[chromIndex][posI] << "\t";

    for (size_t ii = 0; ii < probabilities[siteIndex].size(); ii++) {

      ofstreamExportFwdProb_ << probabilities[siteIndex][ii];
      ofstreamExportFwdProb_ << ((ii < (probabilities[siteIndex].size() - 1)) ? "\t" : "\n");

    }

    siteIndex++;

  }

  ofstreamExportFwdProb_.close();

}


//void DEploidIO::writeLastPairFwdProb( UpdatePairHap & updatePair, size_t chromIndex ){
//ofstreamExportFwdProb.open( strExportPairFwdProb.c_str(), ios::out | ios::app | ios::binary );
//if ( chromIndex == 0 ){ // Print header
//ofstreamExportFwdProb << "CHROM" << "\t" << "POS" << "\t";;
//for ( size_t ii = 0; ii < updatePair.fwdProbs_[0].size(); ii++){
//for ( size_t ij = 0; ij < updatePair.fwdProbs_[0][ii].size(); ij++){
//ofstreamExportFwdProb << (ii+1) << "X" << (ij+1);
//ofstreamExportFwdProb << ((((ii+1) * (ij+1)) < (updatePair.fwdProbs_[0].size()*updatePair.fwdProbs_[0][ii].size()))  ? "\t" : "\n") ;
//}
//}
//}

//size_t siteIndex = 0;
//for ( size_t posI = 0; posI < position_[chromIndex].size(); posI++){
//ofstreamExportFwdProb << chrom_[chromIndex] << "\t" << (int)position_[chromIndex][posI] << "\t";
//for ( size_t ii = 0; ii < updatePair.fwdProbs_[siteIndex].size(); ii++){
//for ( size_t ij = 0; ij < updatePair.fwdProbs_[siteIndex][ii].size(); ij++){
//ofstreamExportFwdProb << updatePair.fwdProbs_[siteIndex][ii][ij];
//ofstreamExportFwdProb << ((((ii+1) * (ij+1)) < (updatePair.fwdProbs_[0].size()*updatePair.fwdProbs_[0][ii].size()))  ? "\t" : "\n") ;
//}
//}
//siteIndex++;
//}

//ofstreamExportFwdProb.close();
//}

