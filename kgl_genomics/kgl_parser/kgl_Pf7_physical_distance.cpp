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

void kgl::LocationCoordinates::initSampleFields(const Pf7SampleRecord& sample_record) {

  auto const& [location_str, location_type] = location_;

  if (location_type == LocationType::City) {

    city_ = sample_record.location1_;
    convertLatLong(sample_record.location1_latitude_, sample_record.location1_longitude_);

  } else {

    city_.clear();
    convertLatLong(sample_record.country_latitude_, sample_record.country_longitude_);

  }

  country_ = sample_record.country_;
  region_ = sample_record.population_;

}


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


void kgl::LocationCoordinates::addSample(std::string sample_id, const std::string& study, const std::string& year_str) {

  sample_id_vec_.push_back(std::move(sample_id));
  size_t year = std::stoll(year_str);
  studies_[study] = year;

}


void kgl::LocationCoordinates::addSample(const Pf7SampleRecord& sample_record) {

  sample_id_vec_.push_back(sample_record.Pf7Sample_id);
  size_t year = std::stoll(sample_record.year_);
  studies_[sample_record.study_] = year;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::Pf7SampleLocation::initializeLocationMaps(const Pf7SampleResource& sample_resource) {

  for (auto const& [sample_id, sample_record] : sample_resource.getMap()) {

    if (not sample_record.location1_.empty()) {

      auto location_iter = location_map_.find(sample_record.location1_);
      if (location_iter == location_map_.end()) {

        auto [iter, result] = location_map_.try_emplace( sample_record.location1_,
                                                                          sample_record.location1_,
                                                                          LocationType::City,
                                                                          sample_record);

        if (not result) {

          ExecEnv::log().error("Pf7SampleLocation::initializeLocationMaps, unexpected; could not store location: {}", sample_record.location1_);

        }

      } else {

        auto& [location, location_record] = *location_iter;
        location_record.addSample(sample_record);

      }

    } // empty()

    if (not sample_record.country_.empty()) {

      auto location_iter = location_map_.find(sample_record.country_);
      if (location_iter == location_map_.end()) {

        auto [iter, result] = location_map_.try_emplace(sample_record.country_,
                                                                         sample_record.country_,
                                                                         LocationType::Country,
                                                                         sample_record);

        if (not result) {

          ExecEnv::log().error("Pf7SampleLocation::initializeLocationMaps, unexpected; could not store country: {}", sample_record.location1_);

        }

      } else {

        auto& [location, location_record] = *location_iter;
        location_record.addSample(sample_record);

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


// Get all locations within a radius from the specified location.
// The bool 'all' parameter (default false) determines if all locations within the great circle radius are used to generate
// the location vector. Otherwise, if the location is a city then only city locations are returned and if a country
// then only country locations are returned.
std::vector<std::string> kgl::Pf7SampleLocation::locationRadius(const std::string& location, double radius, bool all) const {

  std::vector<std::string> locations;

  auto cache_iter = distance_cache_.find(location);
  if (cache_iter == distance_cache_.end()) {

    ExecEnv::log().warn("Pf7SampleLocation::locationRadius; Location: {} not found in distance cache", location);
    return locations;

  }
  auto const& [cache_text, cache_map] = *cache_iter;

  auto location_iter = location_map_.find(location);
  if (location_iter == location_map_.end()) {

    ExecEnv::log().error("Pf7SampleLocation::locationRadius; Unexpected, Location: {} found in cache but not found in location map.", location);
    return locations;

  }
  auto const& [location_text, location_record] = *location_iter;
  auto [id, locationtype] = location_record.location();


  for (auto const& [proximity_location, distance] : cache_map) {

    if (distance <= radius) {

      if (all) {

        locations.push_back(proximity_location);

      } else {

        // Determine the proximity location type (city or country) and only add if the same as the location type.
        auto proximity_iter = location_map_.find(proximity_location);
        if (proximity_iter == location_map_.end()) {

          ExecEnv::log().error("Pf7SampleLocation::locationRadius; Unexpected, Location: {} found in cache but not found in location map.", proximity_location);
          locations.clear();
          return locations;

        }
        auto const& [proximity_text, proximity_record] = *proximity_iter;
        auto [prox_id, proximitytype] = proximity_record.location();

        if (locationtype == proximitytype) {

          locations.push_back(proximity_location);

        }

      }

    }

  }

  return locations;

}

// Get all sample ids. within a radius from the specified location.
// The bool 'all' parameter (default false) determines if all locations within the great circle radius are used to generate
// the sample vector. Otherwise, if the location is a city then only samples (genomes) located in cities are returned and if a country
// then only samples (genomes) in countries are returned.
std::vector<std::string> kgl::Pf7SampleLocation::sampleRadius(const std::string& location, double radius, bool all) const {

  std::vector<std::string> sample_vec;

  auto location_vector = locationRadius(location, radius, all);

  for (auto const& proximity_location : location_vector) {

    auto proximity_iter = location_map_.find(proximity_location);
    if (proximity_iter == location_map_.end()) {

      ExecEnv::log().warn("Pf7SampleLocation::sampleRadius; Location: {} not found in distance cache", proximity_location);
      sample_vec.clear();
      return sample_vec;

    }

    auto const& [proximity_text, proximity_record] = *proximity_iter;
    for (auto const& genome_sample : proximity_record.locationSamples()) {

      sample_vec.push_back(genome_sample);

    }

  }

  return sample_vec;

}
