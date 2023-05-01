//
// Created by kellerberrin on 30/04/23.
//

#ifndef KGL_PF7_SAMPLE_DISTANCE_H
#define KGL_PF7_SAMPLE_DISTANCE_H


#include "kgl_pf7_sample_parser.h"

#include <numbers>

namespace kellerberrin::genome {   //  organization level namespace

enum class LocationType { City, Country };
struct LocationCoordinates {

public:

  // Assumes floating point text from the Pf7SampleResource object.
  LocationCoordinates(const std::string& latitude_text,
                      const std::string& longitude_text,
                      std::string location,
                      LocationType location_type) : location_({std::move(location), location_type}) {

    convertLatLong(latitude_text, longitude_text);

  }
  ~LocationCoordinates() = default;

  [[nodiscard]] double latitudeRadians() const { return latitude_; }
  [[nodiscard]] double longitudeRadians() const { return longitude_; }
  [[nodiscard]] double latitudeDegrees() const { return (latitude_ / (std::numbers::pi * 2.0)) * 360.0; }
  [[nodiscard]] double longitudeDegrees() const { return (longitude_ / (std::numbers::pi * 2.0)) * 360.0; }
  [[nodiscard]] const std::pair<std::string, LocationType>& location() const { return location_; }
  [[nodiscard]] const std::vector<std::string>& locationSamples() const { return sample_id_vec_; }

  void addSample(std::string sample_id) { sample_id_vec_.push_back(std::move(sample_id)); }

  // Great circle calculateDistance between 2 locations.
  [[nodiscard]] double distance_km(const LocationCoordinates& other_location) const;

private:

  double latitude_{0.0};     // In radians, -ve is South of the equator.
  double longitude_{0.0};    // In radians, -ve is West of Greenwich (the Prime Meridian).
  std::pair<std::string, LocationType> location_;  // City or Country.
  std::vector<std::string> sample_id_vec_; // All samples at this location.

  // Used to calculate great circle calculateDistance between 2 locations.
  constexpr static const double EARTH_RADIUS_KM_{6371.0};

  void convertLatLong(const std::string& latitude_text, const std::string& longitude_text);

};

using DistanceCache = std::map<std::string, std::map<std::string, double>>;
using SampleLocationMap = std::map<std::string, LocationCoordinates>;
class Pf7SampleLocation {

public:

  explicit Pf7SampleLocation(const std::shared_ptr<const Pf7SampleResource>& sample_resource_ptr) {

    initializeLocationMaps(sample_resource_ptr);

  }
  ~Pf7SampleLocation() = default;

  [[nodiscard]] const SampleLocationMap& locationMap() const { return location_map_; }
  // Calculate great circle distances in kilometers
  [[nodiscard]] double calculateDistance(const std::string& location1, const std::string& location2) const;
  // Lookup the cache great circle distance.
  [[nodiscard]] double distance(const std::string& location1, const std::string& location2) const;
  // Get all locations within a radius from the specified location.
  [[nodiscard]] std::vector<std::string> locationRadius(const std::string& location, double radius) const;
  // Get all sample ids. within a radius from the specified location.
  [[nodiscard]] std::vector<std::string> sampleRadius(const std::string& location, double radius) const;

private:

  SampleLocationMap location_map_;
  DistanceCache distance_cache_;
  void initializeLocationMaps(const std::shared_ptr<const Pf7SampleResource>& sample_resource_ptr);
  void createDistanceCache();

};



} // Namespace.

#endif //KGL_PF7_SAMPLE_DISTANCE_H
