//
// Created by kellerberrin on 15/06/23.
//

#ifndef KEL_INTERVAL_H
#define KEL_INTERVAL_H


#include <set>
#include <map>
#include <vector>
#include <string>


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Defines a simple right open interval [lower_, upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class OpenRightInterval {

public:

  OpenRightInterval(size_t lower, size_t upper);
  ~OpenRightInterval() = default;
  OpenRightInterval(const OpenRightInterval &copy) = default;

  OpenRightInterval &operator=(const OpenRightInterval &copy) = default;

  [[nodiscard]] size_t lower() const { return lower_; }
  [[nodiscard]] size_t upper() const { return upper_; }
  [[nodiscard]] size_t size() const { return upper_ - lower_; }

  [[nodiscard]] bool containsOffset(size_t offset) const { return offset >= lower_ and offset < upper_; }
  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const { return interval.lower_ >= lower_ and (interval.lower_ + interval.size()) <= upper_; }
  // Returns the {lower, upper} bounds of the intersection interval or {0, 0} indicating no intersection.
  [[nodiscard]] std::pair<size_t, size_t> intersection(const OpenRightInterval &interval) const;
  [[nodiscard]] bool intersects(const OpenRightInterval &interval) const { auto const [lower, upper] = intersection(interval); return upper != 0; }


private:

  size_t upper_{0};
  size_t lower_{0};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Adapters for std::set and std::map to use OpenRightInterval.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Comparison operator used to order intervals within indexed containers std::set or std::map.
struct CompareInterval {

  bool operator()(const OpenRightInterval &lhs, const OpenRightInterval &rhs) const { return lhs.lower() < rhs.lower(); }

};



} // namespace



#endif //KEL_INTERVAL_H
