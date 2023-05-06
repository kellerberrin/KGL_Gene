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
                      LocationType location_type,
                      std::string city,
                      std::string country,
                      std::string region) : location_({std::move(location), location_type}),
                                                    city_(std::move(city)),
                                                    country_(std::move(country)),
                                                    region_(std::move(region)) {

    convertLatLong(latitude_text, longitude_text);

  }
  ~LocationCoordinates() = default;

  [[nodiscard]] double latitudeRadians() const { return latitude_; }
  [[nodiscard]] double longitudeRadians() const { return longitude_; }
  [[nodiscard]] double latitudeDegrees() const { return (latitude_ / (std::numbers::pi * 2.0)) * 360.0; }
  [[nodiscard]] double longitudeDegrees() const { return (longitude_ / (std::numbers::pi * 2.0)) * 360.0; }
  [[nodiscard]] const std::pair<std::string, LocationType>& location() const { return location_; }
  [[nodiscard]] const std::vector<std::string>& locationSamples() const { return sample_id_vec_; }
  // Blank if the location is a country.
  [[nodiscard]] const std::string& city() const { return city_; }
  [[nodiscard]] const std::string& country() const { return country_; }
  [[nodiscard]] const std::string& region() const { return region_; }
  // Study/year pairs.
  [[nodiscard]] const std::map<std::string, size_t>& locationStudies() const { return studies_; }


  // Great circle calculateDistance between 2 locations.
  [[nodiscard]] double distance_km(const LocationCoordinates& other_location) const;

  void addSample(std::string sample_id, const std::string& study, const std::string& year_str);

private:

  double latitude_{0.0};     // In radians, -ve is South of the equator.
  double longitude_{0.0};    // In radians, -ve is West of Greenwich (the Prime Meridian).
  std::pair<std::string, LocationType> location_;  // City or Country. This is the lat/long location.
  std::string city_;  // This is blank if the location is a country.
  std::string country_;
  std::string region_;
  std::vector<std::string> sample_id_vec_; // All samples/genomes at this location.
  std::map<std::string, size_t> studies_; // All study and study year pairs

  // Used to calculate great circle distance between 2 locations.
  constexpr static const double EARTH_RADIUS_KM_{6371.0};

  void convertLatLong(const std::string& latitude_text, const std::string& longitude_text);

};

// Use to 2 maps to cache the distance between any two locations (cities or countries).
using DistanceCache = std::map<std::string, std::map<std::string, double>>;
// Map to store location/country latitude and longitude and all samples located there.
using SampleLocationMap = std::map<std::string, LocationCoordinates>;
class Pf7SampleLocation {

public:

  explicit Pf7SampleLocation(const Pf7SampleResource& sample_resource) {

    initializeLocationMaps(sample_resource);

  }
  ~Pf7SampleLocation() = default;

  // Expose all locations.
  [[nodiscard]] const SampleLocationMap& locationMap() const { return location_map_; }
  // Expose the distance cache.
  [[nodiscard]] const DistanceCache& distanceCache() const { return distance_cache_; }
  // Calculate great circle distances in kilometers
  [[nodiscard]] double calculateDistance(const std::string& location1, const std::string& location2) const;
  // Lookup the cache great circle distance.
  [[nodiscard]] double distance(const std::string& location1, const std::string& location2) const;
  // Get all locations within a radius from the specified location.
  // The bool 'all' parameter (default false) determines if all locations within the great circle radius are used to generate
  // the location vector. Otherwise, if the location is a city then only city locations are returned and if a country
  // then only country locations are returned.
  [[nodiscard]] std::vector<std::string> locationRadius(const std::string& location, double radius, bool all = false) const;
  // Get all sample ids. within a radius from the specified location.
  // The bool 'all' parameter (default false) determines if all locations within the great circle radius are used to generate
  // the sample vector. Otherwise, if the location is a city then only samples (genomes) located in cities are returned and if a country
  // then only samples (genomes) in countries are returned.
  [[nodiscard]] std::vector<std::string> sampleRadius(const std::string& location, double radius, bool all = false) const;

private:

  SampleLocationMap location_map_;
  DistanceCache distance_cache_;
  void initializeLocationMaps(const Pf7SampleResource& sample_resource_ptr);
  void createDistanceCache();

};



} // Namespace.

#endif //KGL_PF7_SAMPLE_DISTANCE_H
