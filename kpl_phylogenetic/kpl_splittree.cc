//
// Created by kellerberrin on 12/12/19.
//


#include "kpl_splittree.h"

namespace kpl = kellerberrin::phylogenetic;

// member function bodies go here


kpl::Split::Split() {

  _mask = 0L;
  _nleaves = 0;
  _bits_per_unit = (CHAR_BIT) * sizeof(Split::split_unit_t);
  clear();
  //std::cout << "Constructing a Split" << std::endl;

}

void kpl::Split::clear() {

  for (auto &u : _bits) {

    u = 0L;

  }

}



bool kpl::Split::operator==(const Split & other) const {

  return (_bits == other._bits);

}


bool kpl::Split::operator<(const Split & other) const {

  assert(_bits.size() == other._bits.size());

  return (_bits < other._bits);

}


void kpl::Split::resize(unsigned nleaves) {

  _nleaves = nleaves;
  unsigned nunits = 1 + ((nleaves - 1)/_bits_per_unit);
  _bits.resize(nunits);

  // create mask used to select only those bits used in final unit
  unsigned num_unused_bits = nunits*_bits_per_unit - _nleaves;
  unsigned num_used_bits = _bits_per_unit - num_unused_bits;
  _mask = 0L;
  split_unit_t unity = 1;

  for (unsigned i = 0; i < num_used_bits; i++) {

    _mask |= (unity << i);

  }

  clear();

}


void kpl::Split::setBitAt(unsigned leaf_index) {

  unsigned unit_index = leaf_index/_bits_per_unit;
  unsigned bit_index = leaf_index - unit_index * _bits_per_unit;

  split_unit_t unity = 1;
  split_unit_t bit_to_set = unity << bit_index;
  _bits[unit_index] |= bit_to_set;

}


kpl::Split::split_unit_t kpl::Split::getBits(unsigned unit_index) const {

  assert(unit_index < _bits.size());

  return _bits[unit_index];

}


bool kpl::Split::getBitAt(unsigned leaf_index) const {

  unsigned unit_index = leaf_index/_bits_per_unit;
  unsigned bit_index = leaf_index - unit_index*_bits_per_unit;
  split_unit_t unity = 1;
  split_unit_t bit_to_set = unity << bit_index;

  return (bool)(_bits[unit_index] & bit_to_set);

}


void kpl::Split::addSplit(const Split & other) {

  unsigned nunits = (unsigned)_bits.size();
  assert(nunits == other._bits.size());

  for (unsigned i = 0; i < nunits; ++i) {

    _bits[i] |= other._bits[i];

  }

}


std::string kpl::Split::createPatternRepresentation() const {

  std::string s;
  unsigned ntax_added = 0;

  for (unsigned i = 0; i < _bits.size(); ++i) {

    for (unsigned j = 0; j < _bits_per_unit; ++j) {

      split_unit_t bitmask = ((split_unit_t)1 << j);
      bool bit_is_set = ((_bits[i] & bitmask) > (split_unit_t)0);

      if (bit_is_set) {

        s += '*';

      }
      else {

        s += '-';

      }

      if (++ntax_added == _nleaves) {

        break;

      }

    }

  }

  return s;

}


bool kpl::Split::isEquivalent(const Split & other) const {

  unsigned nunits = (unsigned)_bits.size();
  assert(nunits > 0);

  // polarity 1 means root is on the same side of both splits
  // polarity 2 means they are inverted relative to one another
  unsigned polarity = 0;
  for (unsigned i = 0; i < nunits; ++i) {
    split_unit_t a = _bits[i];
    split_unit_t b = other._bits[i];
    bool a_equals_b = (a == b);
    bool a_equals_inverse_b = (a == ~b);
    if (i == nunits - 1) {

      a_equals_inverse_b = (a == (~b & _mask));

    }
    bool ok = (a_equals_b || a_equals_inverse_b);
    if (ok) {

      if (polarity == 0) {
        // First unit examined determines polarity
        if (a_equals_b) {

          polarity = 1;

        }
        else {

          polarity = 2;

        }

      }
      else {
        // Polarity determined by first unit used for all subsequent units
        if (polarity == 1 && !a_equals_b ) {

          return false;

        }
        else if (polarity == 2 && !a_equals_inverse_b ) {

          return false;

        }

      }

    }
    else {

      return false;

    }

  }

  // All of the units were equivalent, so that means the splits are equivalent
  return true;

}


bool kpl::Split::isCompatible(const Split & other) const {

  for (unsigned i = 0; i < _bits.size(); ++i) {

    split_unit_t a = _bits[i];
    split_unit_t b = other._bits[i];
    split_unit_t a_and_b = (a & b);
    bool equals_a   = (a_and_b == a);
    bool equals_b   = (a_and_b == b);

    if (a_and_b && !(equals_a || equals_b)) {
      // A failure of any unit to be compatible makes the entire split incompatible
      return false;
    }

  }

  // None of the units were incompatible, so that means the splits are compatible
  return true;

}


bool kpl::Split::conflictsWith(const Split & other) const {

  return !isCompatible(other);

}
