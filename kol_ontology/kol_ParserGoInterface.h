#ifndef KOL_GO_PARSER_INTERFACE
#define KOL_GO_PARSER_INTERFACE

#include "kol_GoGraph.h"

namespace kellerberrin::ontology {


/*! \class ParserGoInterface
	\brief An interface class to define go graph parsers

	This class defines the interface of a go graph parser. Pure virtual methods
	  require that parsers implement these methods.

*/
class ParserGoInterface {

public:

  ParserGoInterface() = default;

  virtual ~ParserGoInterface() = default;

  //! A pure virtual method for parsing the file and returning a GoGraphImpl object.
  /*!
    This pure virtual method requires any parser to have a method that takes
      a filename string and returns a GoGraphImpl object pointer.
  */
  [[nodiscard]] virtual std::shared_ptr<GoGraph> parseGoFile(const std::string &fileName) const = 0;

  //! A pure virtual method for parsing the file and returning a GoGraphImpl object.
  /*!
  This pure virtual method requires any parser to have a method that takes
  a filename string and returns a GoGraphImpl object pointer.
  */
  [[nodiscard]] virtual bool isFileGood(const std::string &filename) const = 0;


};


} // namespace


#endif
