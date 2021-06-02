//
// Created by kellerberrin on 28/3/21.
//

#ifndef KKOL_ONTOLOGYTYPES_H
#define KOL_ONTOLOGYTYPES_H

#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>

namespace kellerberrin::ontology {


template<class T> using OntologySetType = std::unordered_set<T>;
template<class T, class V> using OntologyMapType = std::unordered_map<T, V>;


//template<class T> using OntologySetType = std::set<T>;
//template<class T, class V> using OntologyMapType = std::map<T, V>;



} // namespace

#endif //KOL_ONTOLOGYTYPES_H
