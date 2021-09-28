//
// Created by kellerberrin on 23/9/21.
//

#ifndef KEL_DATE_TIME_H
#define KEL_DATE_TIME_H

#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// General purpose date time class. Implemented using boost::gregorian.
// How many times in the history of C++ programming has this been recreated? Unnecessary and annoying.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization::project level namespace


class DateGP {

public:

  DateGP() = default;
  DateGP(const std::string& date); // Must be a valid delimited date e.g. "2020/1/1" or "2001-Feb-28"
  DateGP(size_t year, size_t month, size_t day); // Must be a valid date, for example DateGP(2001, 2, 29) is an error.
  DateGP(const DateGP&) = default;
  ~DateGP() = default;

  DateGP& operator=(const DateGP&) = default;

  [[nodiscard]] size_t year() const { return year_; }
  [[nodiscard]] size_t month() const { return month_; }
  [[nodiscard]] size_t day() const { return day_; }

  [[nodiscard]] std::string text() const; // returned as YYYY-MMM-DD, e.g. "2020-Jan-01"

  [[nodiscard]] std::string year_text() const { return std::to_string(year_); }
  [[nodiscard]] std::string month_text() const { return std::to_string(month_); }
  [[nodiscard]] std::string day_text() const { return std::to_string(day_); }

  [[nodiscard]] bool operator<(const DateGP& cmp) const;
  [[nodiscard]] bool operator==(const DateGP& cmp) const;
  [[nodiscard]] bool notInitialized() const { return this->operator==(DateGP()); }

  [[nodiscard]] static size_t daysDifference(const DateGP& date1, const DateGP& date2);
  [[nodiscard]] static size_t monthsDifference(const DateGP& date1, const DateGP& date2);

private:

  size_t year_{1901};
  size_t month_{1};
  size_t day_{1};

};


} // namespace


#endif // KEL_DATE_TIME_H
