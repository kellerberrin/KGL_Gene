//
// Created by kellerberrin on 23/9/21.
//

#include "kel_date_time.h"
#include "kel_exec_env.h"



#include <chrono>
#include <format>
#include <locale>
#include <sstream>



namespace kel = kellerberrin;


namespace {

// Tries to parse common delimited forms: 2020-01-01, 2020/1/1,
// 2001-Feb-28, 1-Jan-2020, etc.  Uses the classic (C) locale so month
// abbreviations are always English, matching the original output format.
std::optional<std::chrono::year_month_day> parseDate(std::string_view sv) noexcept {

  std::string s(sv);
  std::istringstream iss(s);
  iss.imbue(std::locale::classic());

  std::chrono::year_month_day ymd;

  auto tryPattern = [&](const char* pattern) -> bool {
    iss.clear();
    iss.seekg(0);
    iss.str(s);
    iss >> std::chrono::parse(pattern, ymd);
    return !iss.fail() && ymd.ok();
  };

  if (tryPattern("%Y-%m-%d")) return ymd; // 2020-01-01
  if (tryPattern("%Y/%m/%d")) return ymd; // 2020/1/1
  if (tryPattern("%Y-%b-%d")) return ymd; // 2001-Feb-28
  if (tryPattern("%Y/%b/%d")) return ymd; // 2001/Feb/28
  if (tryPattern("%d-%b-%Y")) return ymd; // 1-Jan-2020
  if (tryPattern("%d/%b/%Y")) return ymd; // 1/Jan/2020

  return std::nullopt;

}

} // anonymous namespace

kel::DateGP::DateGP(std::string_view date) {

  if (auto parsed = parseDate(date)) {

    date_ = *parsed;

  } else {

    ExecEnv::log().warn("DateGP: unable to parse date string '{}'", date);

  }

}

kel::DateGP::DateGP(size_t year, size_t month, size_t day) {

  std::chrono::year_month_day ymd{
    std::chrono::year{static_cast<int>(year)},
    std::chrono::month{static_cast<unsigned>(month)},
    std::chrono::day{static_cast<unsigned>(day)}
  };

  if (!ymd.ok()) {

    ExecEnv::log().warn("DateGP: invalid date {}-{}-{}", year, month, day);

  }

  date_ = ymd;

}

std::optional<kel::DateGP> kel::DateGP::tryParse(std::string_view date) noexcept {

  if (auto ymd = parseDate(date)) {

    DateGP result;
    result.date_ = *ymd;
    return result;

  }

  return std::nullopt;

}

std::optional<kel::DateGP> kel::DateGP::tryCreate(size_t y, size_t m, size_t d) noexcept {

  std::chrono::year_month_day ymd{
    std::chrono::year{static_cast<int>(y)},
    std::chrono::month{static_cast<unsigned>(m)},
    std::chrono::day{static_cast<unsigned>(d)}
  };

  if (!ymd.ok()) return std::nullopt;

  DateGP result;
  result.date_ = ymd;
  return result;

}

bool kel::DateGP::ok() const noexcept {

  return date_.has_value() && date_->ok();

}

size_t kel::DateGP::year() const {

  if (!ok()) {

    ExecEnv::log().warn("DateGP::year() on uninitialized date");
    return 0;

  }

  return static_cast<size_t>(static_cast<int>(date_->year()));

}

size_t kel::DateGP::month() const {

  if (!ok()) {

    ExecEnv::log().warn("DateGP::month() on uninitialized date");
    return 0;

  }

  return static_cast<unsigned>(date_->month());

}

size_t kel::DateGP::day() const {

  if (!ok()) {

    ExecEnv::log().warn("DateGP::day() on uninitialized date");
    return 0;

  }

  return static_cast<unsigned>(date_->day());

}

std::string kel::DateGP::text() const {

  if (!ok()) {

    ExecEnv::log().warn("DateGP::text() on uninitialized date");
    return {};

  }
  // Classic locale gives English month abbreviations ("Jan", "Feb", ...).
  return std::format(std::locale::classic(), "{:%Y-%b-%d}", *date_);

}

void kel::DateGP::setToday() {

  try {

    auto local = std::chrono::current_zone()->to_local(std::chrono::system_clock::now());
    date_ = std::chrono::year_month_day{floor<std::chrono::days>(local)};

  } catch (const std::exception& e) {

    ExecEnv::log().warn("DateGP::setToday; unable to determine local time zone ({}), using UTC", e.what());
    setUTCDate();

  }

}

void kel::DateGP::setUTCDate() {

  date_ = std::chrono::year_month_day{floor<std::chrono::days>(std::chrono::system_clock::now())};

}

size_t kel::DateGP::daysDifference(const DateGP& date1, const DateGP& date2) {

  if (!date1.ok() || !date2.ok()) {

    ExecEnv::log().warn("DateGP::daysDifference() on uninitialized date");
    return 0;

  }

  auto d1 = std::chrono::sys_days{*date1.date_};
  auto d2 = std::chrono::sys_days{*date2.date_};

  return static_cast<size_t>(std::chrono::abs(d1 - d2).count());

}

size_t kel::DateGP::monthsDifference(const DateGP& date1, const DateGP& date2) {

  if (!date1.ok() || !date2.ok()) {

    ExecEnv::log().warn("DateGP::monthsDifference() on uninitialized date");
    return 0;

  }

  auto months1 = (static_cast<int>(date1.date_->year()) * 12) + static_cast<unsigned>(date1.date_->month());
  auto months2 = (static_cast<int>(date2.date_->year()) * 12) + static_cast<unsigned>(date2.date_->month());

  auto month_diff = months1 > months2 ? months1 - months2 : months2 - months1;

  return static_cast<size_t>(month_diff);

}
