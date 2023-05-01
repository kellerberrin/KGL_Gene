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

void kgl::LocationCoordinates::convertLatLong(const std::string& latitude_text, const std::string& longitude_text) {

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

  auto const& [location_text, location_type] = location_;
  auto const& [other_location_text, other_location_type] = other_location.location_;

  if (location_text == other_location_text) {

    return 0.0;

  }

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


void kgl::Pf7SampleLocation::initializeLocationMaps(const std::shared_ptr<const Pf7SampleResource>& sample_resource_ptr) {

  for (auto const& [sample_id, sample_record] : sample_resource_ptr->getMap()) {

    if (not sample_record.location1_.empty()) {

      auto location_iter = location_map_.find(sample_record.location1_);
      if (location_iter == location_map_.end()) {

        auto [iter, result] = location_map_.try_emplace(sample_record.location1_,
                                                        sample_record.location1_latitude_,
                                                        sample_record.location1_longitude_,
                                                        sample_record.location1_,
                                                        LocationType::City);

        if (not result) {

          ExecEnv::log().error("Pf7SampleLocation::initializeLocationMaps, unexpected; could not store location: {}", sample_record.location1_);

        }

      } else {

        auto& [location, location_record] = *location_iter;
        location_record.addSample(sample_id);

      }

    } // empty()

    if (not sample_record.country_.empty()) {

      auto location_iter = location_map_.find(sample_record.country_);
      if (location_iter == location_map_.end()) {

        auto [iter, result] = location_map_.try_emplace(sample_record.country_,
                                                        sample_record.country_latitude_,
                                                        sample_record.country_longitude_,
                                                        sample_record.country_,
                                                        LocationType::Country);

        if (not result) {

          ExecEnv::log().error("Pf7SampleLocation::initializeLocationMaps, unexpected; could not store country: {}", sample_record.location1_);

        }

      } else {

        auto& [location, location_record] = *location_iter;
        location_record.addSample(sample_id);

      }

    } // Empty

  } // For loop.

  // Create the distance cache.
  createDistanceCache();

}

double kgl::Pf7SampleLocation::calculateDistance(const std::string& location1, const std::string& location2) const {

  if (location1 == location2) {

    return 0.0;

  }

  auto loc1_iter = location_map_.find(location1);
  if (loc1_iter == location_map_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::calculateDistance; Could not find location: {}", location1);
    return 0.0;

  }

  auto loc2_iter = location_map_.find(location2);
  if (loc2_iter == location_map_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::calculateDistance; Could not find location: {}", location2);
    return 0.0;

  }

  auto const& [loc1, loc1_latlong] = *loc1_iter;
  auto const& [loc2, loc2_latlong] = *loc2_iter;

  return loc1_latlong.distance_km(loc2_latlong);

}


void kgl::Pf7SampleLocation::createDistanceCache() {

  for (auto const& [location1, location_record1] : location_map_) {

    for (auto const& [location2, location_record2] : location_map_) {

      auto& distance_map = distance_cache_[location1];
      distance_map[location2] = calculateDistance(location1, location2);

    }

  }

}

double kgl::Pf7SampleLocation::distance(const std::string& location1, const std::string& location2) const {

  auto location1_iter = distance_cache_.find(location1);
  if (location1_iter == distance_cache_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::distance; Location: {} not found in distance cache", location1);
    return 0.0;

  }

  auto const& [location1_text, location_map] = *location1_iter;

  auto location2_iter = location_map.find(location2);
  if (location2_iter == location_map.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::distance; Location: {} not found in distance cache", location2);
    return 0.0;

  }

  auto const& [location2_text, distance] = *location2_iter;

  return distance;

}


std::vector<std::string> kgl::Pf7SampleLocation::locationRadius(const std::string& location, double radius) const {

  std::vector<std::string> locations;

  auto location_iter = distance_cache_.find(location);
  if (location_iter == distance_cache_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::locationRadius; Location: {} not found in distance cache", location);
    return locations;

  }

  auto const& [location_text, location_map] = *location_iter;

  for (auto const& [radius_location, distance] : location_map) {

    if (distance <= radius) {

      locations.push_back(radius_location);

    }

  }

  return locations;

}


std::vector<std::string> kgl::Pf7SampleLocation::sampleRadius(const std::string& location, double radius) const {

  std::vector<std::string> sample_vec;

  auto location_vector = locationRadius(location, radius);

  for (auto const& location : location_vector) {

    auto location_iter = location_map_.find(location);
    if (location_iter == location_map_.end()) {

      ExecEnv::log().warn("Pf7SampleLocation::sampleRadius; Location: {} not found in distance cache", location);
      sample_vec.clear();
      return sample_vec;

    }

    auto const& [location_text, location_record] = *location_iter;
    for (auto const& sample_text : location_record.locationSamples()) {

      sample_vec.push_back(sample_text);

    }

  }

  return sample_vec;

}
