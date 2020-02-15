// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
// Created by kellerberrin on 10/10/17.
//

#ifndef KGL_PATTERNS_H
#define KGL_PATTERNS_H

#include <map>

// Generic code patterns for the project.

namespace kellerberrin::genome {   //  organization level namespace


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

// Pattern to delete elements from an iterable container
// A template boolean function, generally a lambda function,
// evaluates each element to determine if it is to be deleted (true).

template<class T, typename F>
long predicateIterableDelete(T& iterableContainer, F bool_pred) {

  size_t elements_removed = 0;
  auto iterator = iterableContainer.begin();
  while(iterator != iterableContainer.end()) {

    if (bool_pred(iterator)) {

      // C++ 11 erase() returns next iterator to next valid (or end())
      iterator = iterableContainer.erase(iterator);
      ++elements_removed;

    } else {

      ++iterator;

    } // if erase

  } // while duplicate_iterator

  return elements_removed;

}


// Pattern to delete items that match a predicate in a reference container.
// Note that this is quadratic in time and can be slow.
template<class T, typename P>
long deleteIterableReference(T& modified_container, const T& reference_container, P bool_pred) {

  long items_removed = 0;

  for (auto ref_it = reference_container.begin(); ref_it != reference_container.end(); ++ref_it) {

    auto mod_iterator = modified_container.begin();
    while(mod_iterator != modified_container.end()) {

      if (bool_pred(ref_it, mod_iterator)) {

        // C++ 11 erase() returns next iterator to next valid (or end())
        mod_iterator = modified_container.erase(mod_iterator);
        ++items_removed;

      } else {

        ++mod_iterator;

      } // if erase

    } // while mod_iterator

  } // for ref_it

  return items_removed;

}

// Pattern to delete items that do not match a predicate in a reference container.
template<class T, typename P>
long intersectIterable(T& modified_container, const T& reference_container, P bool_pred) {

  long items_removed = 0;

  auto mod_iterator = modified_container.begin();
  while(mod_iterator != modified_container.end()) {

    bool found = false;

    for (auto ref_it = reference_container.begin(); ref_it != reference_container.end(); ++ref_it) {

      if (bool_pred(ref_it, mod_iterator)) {

        found = true;
        break; // element found in the reference container so keep the element

      }

    } // for ref_it

    if (not found) {

      // C++ 11 erase() returns next iterator to next valid (or end())
      mod_iterator = modified_container.erase(mod_iterator);
      ++items_removed;

    } else {

      ++mod_iterator;

    }

  } // while mod_iterator

  return items_removed;

}


}   // end namespace


#endif // KGL_PATTERNS_H
