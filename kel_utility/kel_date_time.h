//
// Created by kellerberrin on 23/9/21.
//

#ifndef KEL_DATE_TIME_H
#define KEL_DATE_TIME_H


#include <chrono>
#include <compare>
#include <optional>
#include <string>
#include <string_view>

namespace kellerberrin {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// General purpose date time class.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class DateGP {
public:
  DateGP() = default;

  explicit DateGP(std::string_view date);
  DateGP(size_t year, size_t month, size_t day);

  DateGP(const DateGP&) = default;
  DateGP& operator=(const DateGP&) = default;
  ~DateGP() = default;

  /// Set the date to the current local date.
  void setToday();
  /// Set the date to the current UTC date.
  void setUTCDate();

  /// True only if the object actually holds a valid date.
  [[nodiscard]] bool ok() const noexcept;
  [[nodiscard]] bool notInitialized() const noexcept { return !ok(); }

  /// Return the year, or 0 if the date is not initialized.
  [[nodiscard]] size_t year() const;
  /// Return the month, or 0 if the date is not initialized.
  [[nodiscard]] size_t month() const;
  /// Return the day, or 0 if the date is not initialized.
  [[nodiscard]] size_t day() const;

  /// Return the date as "YYYY-bbb-DD", e.g. "2020-Jan-01", or empty if not initialized.
  [[nodiscard]] std::string text() const;
  static constexpr size_t TEXTSIZE_{11};

  [[nodiscard]] std::string year_text() const { return std::to_string(year()); }
  [[nodiscard]] std::string month_text() const { return std::to_string(month()); }
  [[nodiscard]] std::string day_text() const { return std::to_string(day()); }

  [[nodiscard]] auto operator<=>(const DateGP&) const = default;
  [[nodiscard]] bool operator==(const DateGP&) const = default;

  /// Return the absolute difference in days between two dates, or 0 if either is uninitialized.
  [[nodiscard]] static size_t daysDifference(const DateGP& date1, const DateGP& date2);
  /// Return the absolute difference in months between two dates, or 0 if either is uninitialized.
  [[nodiscard]] static size_t monthsDifference(const DateGP& date1, const DateGP& date2);

  // Non-throwing alternatives for callers that cannot use exceptions.
  [[nodiscard]] static std::optional<DateGP> tryParse(std::string_view date) noexcept;
  [[nodiscard]] static std::optional<DateGP> tryCreate(size_t y, size_t m, size_t d) noexcept;

private:
  std::optional<std::chrono::year_month_day> date_;
};



} // namespace kellerberrin


#endif // KEL_DATE_TIME_H
