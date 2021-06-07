//
// Created by kellerberrin on 10/10/17.
//

#ifndef KGL_PATTERNS_H
#define KGL_PATTERNS_H

#include <map>

// Generic code patterns for the project.

namespace kellerberrin {   //  organization level namespace


// Pattern to delete duplicates from an iterable container.
// This is a surprisingly tricky function. So let us
// set it up as a pattern so there can be no mistake.

template<class T>
long deleteIterableDuplicates(T& iterableContainer) {

  long duplicates_removed = 0;
  auto iterator = iterableContainer.begin();

  while (iterator != iterableContainer.end()) {

    auto duplicate_iterator = iterator;
    ++duplicate_iterator;

    while(duplicate_iterator != iterableContainer.end()) {

      if (*iterator == *duplicate_iterator) {

        // C++ 11 erase() returns next iterator to next valid (or end())
        duplicate_iterator = iterableContainer.erase(duplicate_iterator);
        ++duplicates_removed;

      } else {

        ++duplicate_iterator;

      } // if erase

    } // while duplicate_iterator

    ++iterator;

  } // while iterator

  return duplicates_removed;

}

// Pattern to check for duplicates from an iterable container.
// This is a surprisingly tricky function. So let us
// set it up as a pattern so there can be no mistake.

template<class T>
long checkDuplicates(T& iterableContainer) {

  long duplicates = 0;

  for (auto iterator = iterableContainer.begin(); iterator != iterableContainer.end(); ++iterator) {

    auto duplicate_iterator = iterator;
    ++duplicate_iterator;

    while(duplicate_iterator != iterableContainer.end()) {

      if (*iterator == *duplicate_iterator) {

        ++duplicates;

      }

      ++duplicate_iterator;

    } // while duplicate_iterator

  } // while iterator

  return duplicates;

}


}   // end namespace


#endif // KGL_PATTERNS_H
