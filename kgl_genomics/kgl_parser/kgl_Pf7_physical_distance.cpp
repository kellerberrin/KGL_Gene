//
// Created by kellerberrin on 30/04/23.
//

#include "kgl_Pf7_physical_distance.h"


#include <cmath>


namespace kgl = kellerberrin::genome;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::LocationCoordinates::LocationCoordinates(const std::string& latitude_text, const std::string& longitude_text) {

  try {

    // In degrees
    latitude_ = std::stod(latitude_text);
    // Convert to radians
    latitude_ = (latitude_ / 360.0) * 2 * std::numbers::pi;

  } catch(std::exception& e) {

    // Allow blank
    if (not latitude_text.empty()) {

      ExecEnv::log().error("LocationCoordinates::LocationCoordinates; Unable to convert text: {} to a latitude; reason: {}", latitude_text, e.what());

    }
    latitude_ = 0.0;

  }

  try {

    // In degrees
    longitude_ = std::stod(longitude_text);
    // Convert to radians
    longitude_ = (longitude_ / 360.0) * 2 * std::numbers::pi;

  } catch(std::exception& e) {

    // Allow blank
    if (not longitude_text.empty()) {

      ExecEnv::log().error("LocationCoordinates::LocationCoordinates; Unable to convert text: {} to a longitude; reason: {}", longitude_text, e.what());
      longitude_ = 0.0;

    }

  }

}


double kgl::LocationCoordinates::distance_km(const LocationCoordinates& other_location) const {

  double spherical_offset = std::sin(latitude_) * std::sin(other_location.latitude_);
  spherical_offset += std::cos(latitude_) * std::cos(other_location.latitude_) * std::cos((other_location.longitude_ - longitude_));

  double distance = std::acos(spherical_offset) * EARTH_RADIUS_KM_;

  return distance;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::Pf7SampleLocation::initializeLocationMaps() {

  for (auto const& [sample_id, sample_record] : sample_resource_ptr_->getMap()) {

    if (not admin_location_.contains(sample_record.location1_)) {

      auto [iter, result] = admin_location_.try_emplace(sample_record.location1_,
                                                                         sample_record.location1_latitude_,
                                                                         sample_record.location1_longitude_);

      if (not result) {

        ExecEnv::log().error("Pf7SampleLocation::initializeLocationMaps, unexpected; could not store location: {}", sample_record.location1_);

      }

    }

    if (not country_location_.contains(sample_record.country_)) {

      auto [iter, result] = country_location_.try_emplace(sample_record.country_,
                                                        sample_record.country_latitude_,
                                                        sample_record.country_longitude_);

      if (not result) {

        ExecEnv::log().error("Pf7SampleLocation::initializeLocationMaps, unexpected; could not store country: {}", sample_record.location1_);

      }

    }

  } // For loop.

}

double kgl::Pf7SampleLocation::locationDistance(const std::string& location1, const std::string& location2) const {

  auto loc1_iter = admin_location_.find(location1);
  if (loc1_iter == admin_location_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::locationDistance; Could not find location: {}", location1);
    return 0.0;

  }

  auto loc2_iter = admin_location_.find(location2);
  if (loc2_iter == admin_location_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::locationDistance; Could not find location: {}", location2);
    return 0.0;

  }

  auto const& [loc1, loc1_latlong] = *loc1_iter;
  auto const& [loc2, loc2_latlong] = *loc2_iter;

  return loc1_latlong.distance_km(loc2_latlong);

}

double kgl::Pf7SampleLocation::countryDistance(const std::string& country1, const std::string& country2) const {


  auto country1_iter = country_location_.find(country1);
  if (country1_iter == country_location_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::contryDistance; Could not find Country: {}", country1);
    return 0.0;

  }

  auto country2_iter = country_location_.find(country2);
  if (country2_iter == country_location_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::countryDistance; Could not find Country: {}", country2);
    return 0.0;

  }

  auto const& [loc1, country1_latlong] = *country1_iter;
  auto const& [loc2, country2_latlong] = *country2_iter;

  return country1_latlong.distance_km(country2_latlong);

}


double kgl::Pf7SampleLocation::sampleLocationDistance(const std::string& sample1, const std::string& sample2) const {

  auto const& sampleMap = sample_resource_ptr_->getMap();

  auto sample1_iter = sampleMap.find(sample1);
  if (sample1_iter == sampleMap.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::sampleLocationDistance; Could not find Sample: {}", sample1);
    return 0.0;

  }

  auto sample2_iter = sampleMap.find(sample1);
  if (sample2_iter == sampleMap.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::sampleLocationDistance; Could not find Sample: {}", sample2);
    return 0.0;

  }

  auto const& [sample1_id, sample1_record] = *sample1_iter;
  auto const& [sample2_id, sample2_record] = *sample2_iter;

  return locationDistance(sample1_record.location1_, sample2_record.location1_);

}


double kgl::Pf7SampleLocation::sampleCountryDistance(const std::string& sample1, const std::string& sample2) const {


  auto const& sampleMap = sample_resource_ptr_->getMap();

  auto sample1_iter = sampleMap.find(sample1);
  if (sample1_iter == sampleMap.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::sampleCountryDistance; Could not find Sample: {}", sample1);
    return 0.0;

  }

  auto sample2_iter = sampleMap.find(sample1);
  if (sample2_iter == sampleMap.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::sampleCountryDistance; Could not find Sample: {}", sample2);
    return 0.0;

  }

  auto const& [sample1_id, sample1_record] = *sample1_iter;
  auto const& [sample2_id, sample2_record] = *sample2_iter;

  return countryDistance(sample1_record.country_, sample2_record.country_);

}
