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
*/

#include <ctime>
#include <iterator>
#include <cassert>        // assert
#include <iomanip>        // std::setw
#include "kgd_utility.h"    // normailize by sum
#include "kgd_updateHap.h"  // chromPainting
#include "kgd_dEploidIO.h"
#include "kgd_ibd.h"
#include "kgd_deploid_app.h"

namespace kgl = kellerberrin::genome;
namespace kgd = kellerberrin::deploid;


DEploidIO::DEploidIO() {

  init(); // Reset to default values before parsing
  parse();
  checkInput();
  finalize();

}


void DEploidIO::init() {

  excludedMarkers = std::make_shared<ExcludeMarker>();
  panel = nullptr;
  vcfReaderPtr_ = std::make_shared<VcfReader>(kgd::DeploidExecEnv::getArgs().vcfFile);

  setDoExportRecombProb(false);
  setrandomSeedWasGiven(false);
  setCompressVcf(false);
  setInitialPropWasGiven(false);
  setInitialHapWasGiven(false);
  initialProp.clear();
  setPleaseCheckInitialP(true);
  setExcludeSites( false );
  set_seed( (unsigned)0 );
  setUsePanel(true);
  precision_ = 8;
  prefix_ = "pf3k-kgd_deploid";
  setKStrainWasManuallySet(false);
  setKStrainWasSetByHap(false);
  setKStrainWasSetByProp(false);
  setKstrain(5);
  nMcmcSample_ = 800;
  setDoUpdateProp( true );
  setDoUpdatePair( true );
  setDoUpdateSingle( true );
  setDoExportPostProb( false );
  setDoLsPainting( false );
  setDoIbdPainting( false );
  setUseIBD( false );
  setDoExportSwitchMissCopy ( true );
  setDoAllowInbreeding( false );
  mcmcBurn_ = 0.5;
  mcmcMachineryRate_ = 5;
  missCopyProb_ = 0.01;
  useConstRecomb_ = false;
  setForbidCopyFromSame( false );
  constRecombProb_ = 1.0;
  averageCentimorganDistance_ = 15000.0;
  setScalingFactor(100.0);
  setParameterG(20.0);
  setParameterSigma(5.0);
  setIBDSigma(20.0);
  setUseVcf(false);
  setDoExportVcf(false);
  setDoComputeLLK(false);
  refFileName_.clear();
  altFileName_.clear();
  plafFileName_.clear();
  panelFileName_.clear();
  excludeFileName_.clear();
  getTime(true);
  compileTime_ = "";
  dEploidGitVersion_ = kgd::DeploidExecEnv::VERSION;

}


void DEploidIO::getTime( bool isStartingTime ){
    time_t now = time(0);
    // convert now to string form
    char* dt = ctime(&now);
    if ( isStartingTime ){
        startingTime_ = string(dt);
    } else {
        endTime_ = string(dt);
    }
}


void DEploidIO::finalize(){
    if ( doIbdPainting() | doComputeLLK() ){
        if (!initialPropWasGiven()){
            throw InitialPropUngiven("");
        }
    }

    if ( useIBD() && kStrain() == 1){
        throw InvalidK();
    }

    if ( compressVcf() && !doExportVcf() ){
        throw VcfOutUnSpecified("");
    }

    if ( !randomSeedWasGiven_ ){
        set_seed( (unsigned)(time(0)) );
    }

    if ( excludeSites() ){
        excludedMarkers->readFromFile(excludeFileName_.c_str());
    }

    if ( useVcf() ){ // read vcf files, and parse it to refCount and altCount
        if ( excludeSites() ){
            vcfReaderPtr_->findAndKeepMarkers (excludedMarkers);
        }

        vcfReaderPtr_->finalize(); // Finalize after remove variantlines
        refCount_ = vcfReaderPtr_->refCount;
        altCount_ = vcfReaderPtr_->altCount;
    } else {
        TxtReader ref;
        ref.readFromFile(refFileName_.c_str());
        if ( excludeSites() ){
            ref.findAndKeepMarkers( excludedMarkers );
        }
        refCount_ = ref.info_;

        TxtReader alt;
        alt.readFromFile(altFileName_.c_str());
        if ( excludeSites() ){
            alt.findAndKeepMarkers( excludedMarkers );
        }
        altCount_ = alt.info_;
    }

    TxtReader plaf;
    plaf.readFromFile(plafFileName_.c_str());
    if ( excludeSites() ){
        plaf.findAndKeepMarkers( excludedMarkers );
    }
    plaf_ = plaf.info_;
    chrom_ = plaf.chrom_;
    position_ = plaf.position_;
    indexOfChromStarts_ = plaf.indexOfChromStarts_;

    nLoci_ = refCount_.size();

    if ( nLoci_ != plaf_.size() ){
        throw LociNumberUnequal( plafFileName_ );
    }

    if ( nLoci_ != altCount_.size() ){
        throw LociNumberUnequal( altFileName_ );
    }
    (void)removeFilesWithSameName();

    readPanel();
    IBDpathChangeAt = vector <double>(nLoci());
    finalIBDpathChangeAt = vector <double>(nLoci());

    siteOfTwoSwitchOne = vector <double>(nLoci());
    siteOfTwoMissCopyOne = vector <double>(nLoci());
    siteOfTwoSwitchTwo = vector <double>(nLoci());
    siteOfTwoMissCopyTwo = vector <double>(nLoci());
    siteOfOneSwitchOne = vector <double>(nLoci());
    siteOfOneMissCopyOne = vector <double>(nLoci());

    finalSiteOfTwoSwitchOne = vector <double>(nLoci());
    finalSiteOfTwoMissCopyOne = vector <double>(nLoci());
    finalSiteOfTwoSwitchTwo = vector <double>(nLoci());
    finalSiteOfTwoMissCopyTwo = vector <double>(nLoci());
    finalSiteOfOneSwitchOne = vector <double>(nLoci());
    finalSiteOfOneMissCopyOne = vector <double>(nLoci());
}


void DEploidIO::removeFilesWithSameName(){
    strExportProp = prefix_ + ".prop";
    strExportLLK = prefix_ + ".llk";
    strExportHap = prefix_ + ".hap";

    strIbdExportProp = prefix_ + ".ibd.prop";
    strIbdExportLLK = prefix_ + ".ibd.llk";
    strIbdExportHap = prefix_ + ".ibd.hap";
    strIbdExportProbs = prefix_ + ".ibd.probs";

    strExportVcf = prefix_ + ".vcf";
    if ( compressVcf() ){
        strExportVcf += ".gz";
    }
    strExportLog =  prefix_ + ((doLsPainting()) ? ".painting":"") + ".log";
    strExportRecombProb = prefix_ + ".recomb";

    strExportExtra = prefix_ + ".extra";

    if ( doLsPainting() == false ){
        if (useIBD()){
            remove(strIbdExportProp.c_str());
            remove(strIbdExportLLK.c_str());
            remove(strIbdExportHap.c_str());
        }
        remove(strExportLLK.c_str());
        remove(strExportHap.c_str());
        remove(strExportVcf.c_str());
        remove(strExportProp.c_str());
        remove(strExportExtra.c_str());
        remove(strIbdExportProbs.c_str());
    }

    if (doLsPainting() || doExportPostProb() ){
        if (useIBD()){
            strIbdExportSingleFwdProbPrefix = prefix_ + ".ibd.single";
            for ( size_t i = 0; i < kStrain_ ; i++ ){
                string tmpStrExportSingleFwdProb = strIbdExportSingleFwdProbPrefix + to_string(i);
                remove(tmpStrExportSingleFwdProb.c_str());
            }
            strIbdExportPairFwdProb = prefix_ + ".ibd.pair";
            remove(strIbdExportPairFwdProb.c_str());
        }
        strExportSingleFwdProbPrefix = prefix_ + ".single";
        for ( size_t i = 0; i < kStrain_ ; i++ ){
            string tmpStrExportSingleFwdProb = strExportSingleFwdProbPrefix + to_string(i);
            remove(tmpStrExportSingleFwdProb.c_str());
        }
        strExportPairFwdProb = prefix_ + ".pair";
        remove(strExportPairFwdProb.c_str());
    }
    remove(strExportLog.c_str());
    remove(strExportRecombProb.c_str());

}


void DEploidIO::parse () {

  const kgd::DeploidArgs &args = kgd::DeploidExecEnv::getArgs();

  vcfFileName_ = args.vcfFile; // -vcf
  setUseVcf(true);

  if (args.vcfOutFile != kgd::DeploidArgs::NOT_SPECIFIED) {  // -vcfOut

    setDoExportVcf(true);

  }

  plafFileName_ = args.plafFile;  // -plaf

  if (args.panelFile != kgd::DeploidArgs::NOT_SPECIFIED) { // -panel

    panelFileName_ = args.panelFile;
    setUsePanel(true);

  } else {

    setUsePanel(false);

  }

  setUsePanel(not args.noPanelFlag);  // The -noPanel flag is the compliment of UsePanel.
  if (not usePanel()) {

    setDoExportPostProb(false);
    setDoExportSwitchMissCopy(false);
    setDoAllowInbreeding(false);

  }

  if (args.excludeFile != kgd::DeploidArgs::NOT_SPECIFIED) {  // -exclude

    setExcludeSites(true);
    excludeFileName_ = args.excludeFile;

  }

  prefix_ = args.outputTemplate; // -o
  kStrain_ = args.maxStrains; // -k
  setKStrainWasManuallySet(true);
  nMcmcSample_ = args.MCMCSamples; // -nSamples
  mcmcBurn_ = args.MCMCBurnRate;  // -burn
  mcmcMachineryRate_ = args.MCMCSampleRate;  // -rate
  if (args.inbreedingProbabilitiesFlag) { // -inbreeding

    setDoAllowInbreeding(true);
    setDoExportPostProb(true);

  }
  if (args.identityByDescentPainting) setDoIbdPainting(true); // -idbPainting
  initialProp = args.initialStrainProportions; // -initialP
  setUseIBD(args.identityByDescentFlag); // -idb

  set_seed(args.MCMCRandomSeed); // -seed
  setrandomSeedWasGiven(true);

}

void DEploidIO::checkInput(){

    if ( refFileName_.size() == 0 && useVcf() == false ){
        throw FileNameMissing ( "Ref count" );}
    if ( altFileName_.size() == 0 && useVcf() == false ){
        throw FileNameMissing ( "Alt count" );}
    if ( plafFileName_.size() == 0 ){
        throw FileNameMissing ( "PLAF" );}
    if ( usePanel() && panelFileName_.size() == 0 && !doIbdPainting() && !doComputeLLK() ){
        throw FileNameMissing ( "Reference panel" );}
    if ( initialPropWasGiven() && ( abs(sumOfVec(initialProp) - 1.0) > 0.00001 ) && pleaseCheckInitialP() ){
        throw SumOfPropNotOne ( to_string(sumOfVec(initialProp)) );}
    if ( initialPropWasGiven() ){
        if ( kStrainWasManuallySet() == true ){
        } else {
            // set k strain by proportion length
        }
    }
}





vector <double> DEploidIO::computeExpectedWsafFromInitialHap(){
    // Make this a separate function
    // calculate expected wsaf
    vector <double> expectedWsaf (initialHap.size(), 0.0);
    for ( size_t i = 0; i < initialHap.size(); i++ ){
        assert( kStrain_ == initialHap[i].size() );
        for ( size_t k = 0; k < kStrain_; k++){
            expectedWsaf[i] += initialHap[i][k] * finalProp[k];
        }
        assert ( expectedWsaf[i] >= 0 );
        //assert ( expectedWsaf[i] <= 1.0 );
    }
    return expectedWsaf;
}


void DEploidIO::computeLLKfromInitialHap(){
    for ( auto const& value: initialProp ){
        finalProp.push_back(value);
    }

    vector <double> expectedWsaf = computeExpectedWsafFromInitialHap();
    if (expectedWsaf.size() != refCount_.size()){
        throw LociNumberUnequal("Hap length differs from data!");
    }
    vector <double> llk = calcLLKs ( refCount_, altCount_, expectedWsaf, 0, expectedWsaf.size(), scalingFactor());
    llkFromInitialHap_ = sumOfVec(llk);
}


void DEploidIO::chromPainting(){
    dout << "Painting haplotypes in" << initialHapFileName_ <<endl;

    if ( initialPropWasGiven() == false ){
        clog << "Initial proportion was not specified. Set even proportions" << endl;
        double evenProp = 1.0 / (double)kStrain();
        for ( size_t i = 0; i < kStrain(); i++){
            initialProp.push_back(evenProp);
        }
    }

    for ( auto const& value: initialProp ){
        finalProp.push_back(value);
    }

    // Painting posterior probabilities

    // Export the p'
    // Make this a separate class
    //vector < vector <double> > hap;
    //for ( size_t siteI = 0; siteI < decovolutedStrainsToBeRead.content_.size(); siteI++ ){
        //vector <double> tmpHap;
        //for ( size_t tmpk = 0; tmpk < kStrain_; tmpk++ ){
            //tmpHap.push_back(decovolutedStrainsToBeRead.content_[siteI][tmpk]);
        //}
        //hap.push_back(tmpHap);
    //}

    //vector < vector <double>> hap = decovolutedStrainsToBeRead.content_;

    vector <double> expectedWsaf = computeExpectedWsafFromInitialHap();

    MersenneTwister tmpRg(randomSeed());

    if ( doAllowInbreeding() == true ){
        panel->initializeUpdatePanel(panel->truePanelSize()+kStrain_-1);
    }

    for ( size_t tmpk = 0; tmpk < kStrain_; tmpk++ ){
        if ( doAllowInbreeding() == true ){
            panel->updatePanelWithHaps(panel->truePanelSize()+kStrain_-1, tmpk, initialHap);
        }

        for ( size_t chromi = 0 ; chromi < indexOfChromStarts_.size(); chromi++ ){
            size_t start = indexOfChromStarts_[chromi];
            size_t length = position_[chromi].size();
            dout << "Painting Chrom "<< chromi << " from site "<< start << " to " << start+length << endl;

            UpdateSingleHap updatingSingle( refCount_,
                                      altCount_,
                                      plaf_,
                                      expectedWsaf,
                                      finalProp, initialHap, &tmpRg,
                                      start, length,
                                      panel, missCopyProb_, scalingFactor(),
                                      tmpk);

            if ( doAllowInbreeding() == true ){
                updatingSingle.setPanelSize(panel->inbreedingPanelSize());
            }
            updatingSingle.painting( refCount_, altCount_, expectedWsaf, finalProp, initialHap);
            //writeLastSingleFwdProb( updatingSingle.fwdProbs_, chromi, tmpk, false ); // false as not using ibd
            writeLastSingleFwdProb( updatingSingle.fwdBwdProbs_, chromi, tmpk, false ); // false as not using ibd
        }
    }
}


void DEploidIO::readPanel(){
    if ( usePanel() == false ){
        return;
    }
    if ( doIbdPainting() | doComputeLLK() ){
        return;
    }

    panel = std::make_shared<Panel>();
    panel->readFromFile(panelFileName_.c_str());
    if ( excludeSites() ){
        panel->findAndKeepMarkers( excludedMarkers );
    }

    panel->computeRecombProbs( averageCentimorganDistance(), parameterG(), useConstRecomb(), constRecombProb(), forbidCopyFromSame() );
    panel->checkForExceptions( nLoci(), panelFileName_ );
}



void DEploidIO::getIBDprobsIntegrated(vector < vector <double> > &prob){
    if (prob.size() !=  nLoci()){
        throw InvalidInput("Invalid probabilities! Check size!");
    }

    assert(ibdProbsIntegrated.size() == 0);

    for (size_t i = 0; i < prob[0].size(); i++){
        ibdProbsIntegrated.push_back(0.0);
    }

    for ( size_t siteIndex = 0; siteIndex < nLoci(); siteIndex++ ){
        for (size_t i = 0; i < prob[siteIndex].size(); i++){
            ibdProbsIntegrated[i] += prob[siteIndex][i];
        }
    }
    normalizeBySum(ibdProbsIntegrated);
}


void DEploidIO::computeEffectiveKstrain(vector <double> proportion){
    double tmpSumSq = 0.0;
    for (double p : proportion){
        tmpSumSq += p * p;
    }
    effectiveKstrain_ = 1.0 / tmpSumSq;
}


void DEploidIO::computeInferredKstrain(vector <double> proportion){
    inferredKstrain_ = 0;
    for (double p : proportion){
        if ( p > 0.01 ){
            inferredKstrain_ += 1;
        }
    }
}


void DEploidIO::computeAdjustedEffectiveKstrain(){
    adjustedEffectiveKstrain_ = effectiveKstrain_;
    if ( (inferredKstrain_ == 2) & (ibdProbsIntegrated.size() == 2)){
        if ( ibdProbsIntegrated[1] > 0.95 ){
            adjustedEffectiveKstrain_ = 1;
        }
    }
}



void DEploidIO::paintIBD(){
    vector <double> goodProp;
    vector <size_t> goodStrainIdx;

    if ( doIbdPainting() ){
        finalProp = initialProp;
    }

    for ( size_t i = 0; i < finalProp.size(); i++){
        if (finalProp[i] > 0.01){
            goodProp.push_back(finalProp[i]);
            goodStrainIdx.push_back(i);
        }
    }

    if (goodProp.size() == 1){
        return;
    }

    DEploidIO tmpDEploidIO; // (*this);
    tmpDEploidIO.setKstrain(goodProp.size());
    tmpDEploidIO.setInitialPropWasGiven(true);
    tmpDEploidIO.initialProp = goodProp;
    tmpDEploidIO.finalProp = goodProp;
    tmpDEploidIO.refCount_ = refCount_;
    tmpDEploidIO.altCount_ = altCount_;
    tmpDEploidIO.plaf_ = plaf_;
    tmpDEploidIO.nLoci_= nLoci();
    tmpDEploidIO.position_ = position_;
    tmpDEploidIO.chrom_ = chrom_;
    //tmpDEploidIO.useConstRecomb_ = true;
    //tmpDEploidIO.constRecombProb_ = 0.000001;

    //tmpDEploidIO.writeLog (&std::cout);

    MersenneTwister tmpRg(randomSeed());
    IBDpath tmpIBDpath;
    tmpIBDpath.init(tmpDEploidIO, &tmpRg);
    tmpIBDpath.buildPathProbabilityForPainting(goodProp);
    ibdLLK_ = tmpIBDpath.bestPath(goodProp);
    ibdProbsHeader = tmpIBDpath.getIBDprobsHeader();
    getIBDprobsIntegrated(tmpIBDpath.fwdbwd);
    writeIBDpostProb(tmpIBDpath.fwdbwd, ibdProbsHeader);
}

