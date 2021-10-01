//
// Created by kellerberrin on 23/9/21.
//

#include "kel_date_time.h"

#include "kel_exec_env.h"

#include "boost/date_time/gregorian/gregorian.hpp"

// Define namespace alias
namespace kel = kellerberrin;
namespace bg = boost::gregorian;

// Must be a delimited date e.g. "1/1/2020" or "1-Jan-2020"
kel::DateGP::DateGP(const std::string& date_str) {

  try {

    bg::date boostdate(bg::from_string(date_str));
    year_ = boostdate.year();
    month_ = boostdate.month();
    day_ = boostdate.day();

  } catch(std::exception& e) {

    ExecEnv::log().error("DateGP::DateGP; unable to parse date string: '{}', exception: {}", date_str, e.what());

  }

}

kel::DateGP::DateGP(size_t year, size_t month, size_t day) {

  try {

    bg::date boostdate(year, month, day);
    year_ = boostdate.year();
    month_ = boostdate.month();
    day_ = boostdate.day();

  } catch(std::exception& e) {

    ExecEnv::log().error("DateGP::DateGP; invalid date: {}-{}-{}, exception: {}", day, month, year, e.what());

  }

}

// The date object is initialized to the current local date.
void kel::DateGP::setToday() {

  bg::date boostdate (bg::day_clock::local_day());

  year_ = boostdate.year();
  month_ = boostdate.month();
  day_ = boostdate.day();

}

// The date object is initialized to the current UTC (Greenwich) date.
void kel::DateGP::setUTCDate() {

  bg::date boostdate (bg::day_clock::universal_day());

  year_ = boostdate.year();
  month_ = boostdate.month();
  day_ = boostdate.day();

}

std::string kel::DateGP::text() const {

  bg::date boostdate(year_, month_, day_);
  return to_simple_string(boostdate);

}

bool kel::DateGP::operator<(const DateGP& cmp) const {

  bg::date thisdate(year_, month_, day_);
  bg::date cmpdate(cmp.year_, cmp.month_, cmp.day_);

  return thisdate < cmpdate;

}

bool kel::DateGP::operator==(const DateGP& cmp) const {

  bg::date thisdate(year_, month_, day_);
  bg::date cmpdate(cmp.year_, cmp.month_, cmp.day_);

  return thisdate == cmpdate;

}

size_t kel::DateGP::daysDifference(const DateGP& date1, const DateGP& date2) {

  bg::date d1(date1.year_, date1.month_, date1.day_);
  bg::date d2(date2.year_, date2.month_, date2.day_);

  if (d1 < d2) {

    bg::date_duration dd = d2 - d1;
    return dd.days();

  } else {

    bg::date_duration dd = d1 - d2;
    return dd.days();

  }

}


size_t kel::DateGP::monthsDifference(const DateGP& date1, const DateGP& date2) {

  bg::date d1(date1.year_, date1.month_, date1.day_);
  bg::date d2(date2.year_, date2.month_, date2.day_);

  if (d1 < d2) {

    size_t month_diff = (date2.year_ - date1.year_) * 12;
    month_diff += date2.month_ - date1.month_;
    return month_diff;

  } else {

    size_t month_diff = (date1.year_ - date2.year_) * 12;
    month_diff += date1.month_ - date2.month_;
    return month_diff;

  }

}


