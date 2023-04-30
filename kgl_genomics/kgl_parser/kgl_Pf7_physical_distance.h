//
// Created by kellerberrin on 30/04/23.
//

#ifndef KGL_PF7_SAMPLE_DISTANCE_H
#define KGL_PF7_SAMPLE_DISTANCE_H


#include "kgl_pf7_sample_parser.h"

#include <numbers>

namespace kellerberrin::genome {   //  organization level namespace


struct LocationCoordinates {

public:

  // Assumes floating point text from the Pf7SampleResource object.
  LocationCoordinates(const std::string& latitude_text, const std::string& longitude_text);
  ~LocationCoordinates() = default;

  [[nodiscard]] double latitudeRadians() const { return latitude_; }
  [[nodiscard]] double longitudeRadians() const { return longitude_; }
  [[nodiscard]] double latitudeDegrees() const { return (latitude_ / (std::numbers::pi * 2.0)) * 360.0; }
  [[nodiscard]] double longitudeDegrees() const { return (longitude_ / (std::numbers::pi * 2.0)) * 360.0; }

  // Great circle distance between 2 locations.
  [[nodiscard]] double distance_km(const LocationCoordinates& other_location) const;

private:

  double latitude_{0.0};     // In radians, -ve is South of the equator.
  double longitude_{0.0};    // In radians, -ve is West of Greenwich (the Prime Meridian).

  // Used to calculate great circle distance between 2 locations.
  constexpr static const double EARTH_RADIUS_KM_{6371.0};


};
using SampleLocationMap = std::map<std::string, LocationCoordinates>;


class Pf7SampleLocation {

public:

  explicit Pf7SampleLocation(std::shared_ptr<const Pf7SampleResource> sample_resource_ptr)
  : sample_resource_ptr_(std::move(sample_resource_ptr)) {

    initializeLocationMaps();

  }
  ~Pf7SampleLocation() = default;

  [[nodiscard]] const SampleLocationMap& adminLocations() const { return admin_location_; }
  [[nodiscard]] const SampleLocationMap& countryLocations() const { return country_location_; }

  // Great circle distances in kilometers
  [[nodiscard]] double locationDistance(const std::string& location1, const std::string& location2) const;
  [[nodiscard]] double countryDistance(const std::string& country1, const std::string& country2) const;

  [[nodiscard]] double sampleLocationDistance(const std::string& sample1, const std::string& sample2) const;
  [[nodiscard]] double sampleCountryDistance(const std::string& sample1, const std::string& sample2) const;

private:

  SampleLocationMap admin_location_;
  SampleLocationMap country_location_;
  std::shared_ptr<const Pf7SampleResource> sample_resource_ptr_;

  void initializeLocationMaps();

};



} // Namespace.

#endif //KGL_PF7_SAMPLE_DISTANCE_H
