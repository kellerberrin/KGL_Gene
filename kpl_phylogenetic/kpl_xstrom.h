//
// Created by kellerberrin on 10/12/19.
//

#ifndef KPL_XSTROM_H
#define KPL_XSTROM_H


#include <boost/format.hpp>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace



class XStrom : public std::exception {

public:

  XStrom() throw() {}

  XStrom(const std::string s) throw() : _msg() { _msg = s; }

  XStrom(const boost::format &f) throw() : _msg() { _msg = boost::str(f); }

  virtual             ~XStrom() throw() {}

  const char *what() const throw() { return _msg.c_str(); }

private:

  std::string _msg;

};


}  // phylogenetic
} // kellerberrin


#endif // KPL_XSTROM_H
