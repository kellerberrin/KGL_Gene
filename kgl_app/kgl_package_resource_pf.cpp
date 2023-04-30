//
// Created by kellerberrin on 5/04/23.
//

#include "kgl_package.h"
#include "kgl_pf7_sample_parser.h"
#include "kgl_pf7_fws_parser.h"
#include "kgl_pf7_genetic_distance_parser.h"


namespace kgl = kellerberrin::genome;



void kgl::ExecutePackage::loadPf7SampleResource(const std::string& resource_type,
                                                const std::string& Pf7_sample_ident,
                                                const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, Pf7_sample_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource, Pf7 Sample Data: {}, not defined", Pf7_sample_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(ResourceProperties::PF7SAMPLE_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource, Ident: {} Pf7 Sample Data file not defined", Pf7_sample_ident);

  }

  ParsePf7Sample Pf7_sample_parser;
  if (not Pf7_sample_parser.parsePf7SampleFile(file_name_opt.value())) {

    ExecEnv::log().critical("ExecutePackage::loadPf7SampleResource; failed to create Pf7 Sample Data resource from file: {}", file_name_opt.value());

  }

  std::shared_ptr<Pf7SampleResource> Pf7Sample_ptr(std::make_shared<Pf7SampleResource>(Pf7_sample_ident, Pf7_sample_parser.getPf7SampleVector()));

  resource_ptr->addResource(Pf7Sample_ptr);

}

void kgl::ExecutePackage::loadPf7FwsResource(const std::string& resource_type,
                                             const std::string& Pf7_Fws_ident,
                                             const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, Pf7_Fws_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7FwsResource, Pf7 FWS Data: {}, not defined", Pf7_Fws_ident);

  }

  auto const& params = params_opt.value();
  auto file_name_opt = params.getParameter(ResourceProperties::PF7FWS_FILE_);
  if (not file_name_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7FwsResource, Ident: {} Pf7 FWS Data file not defined", Pf7_Fws_ident);

  }

  ParsePf7Fws Pf7_Fws_parser;
  if (not Pf7_Fws_parser.parsePf7FwsFile(file_name_opt.value())) {

    ExecEnv::log().critical("ExecutePackage::loadPf7FwsResource; failed to create Pf7 FWS Data resource from file: {}", file_name_opt.value());

  }

  std::shared_ptr<Pf7FwsResource> Pf7Fws_ptr(std::make_shared<Pf7FwsResource>(Pf7_Fws_ident, Pf7_Fws_parser.getPf7FwsVector()));

  resource_ptr->addResource(Pf7Fws_ptr);

}



void kgl::ExecutePackage::loadPf7DistanceResource(const std::string& resource_type,
                                                  const std::string& Pf7_Distance_ident,
                                                  const std::shared_ptr<AnalysisResources>& resource_ptr) const {

  auto params_opt = runtime_config_.resourceDefMap().retrieve(resource_type, Pf7_Distance_ident);
  if (not params_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7DistanceResource, Pf7 Distance Data: {}, not defined", Pf7_Distance_ident);

  }

  auto const& params = params_opt.value();
  auto matrix_file_opt = params.getParameter(ResourceProperties::PF7DISTANCE_MATRIXFILE_);
  if (not matrix_file_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7DistanceResource, Ident: {} Pf7 Distance matrix file not defined", Pf7_Distance_ident);

  }

  auto sampleid_file_opt = params.getParameter(ResourceProperties::PF7DISTANCE_IDFILE_);
  if (not sampleid_file_opt) {

    ExecEnv::log().critical("ExecutePackage::loadPf7DistanceResource, Ident: {} Pf7 Distance sample id file not defined", Pf7_Distance_ident);

  }

  ParsePf7Distance Pf7_Distance_parser;
  if (not Pf7_Distance_parser.parsePf7Distance(matrix_file_opt.value(), sampleid_file_opt.value())) {

    ExecEnv::log().critical("ExecutePackage::loadPf7FwsResource; failed to create Pf7 Distance Matrix resource from matrix file: {}, sample id file: {}",
                            matrix_file_opt.value(), sampleid_file_opt.value());

  }

  std::shared_ptr<Pf7GeneticDistanceResource> Pf7Distance_ptr(std::make_shared<Pf7GeneticDistanceResource>(Pf7_Distance_ident,
                                                                                                           Pf7_Distance_parser.getSampleMap(),
                                                                                                           Pf7_Distance_parser.getDistanceMatrix()));

  resource_ptr->addResource(Pf7Distance_ptr);

}

