#ifndef GO_PARSER_INTERFACE
#define GO_PARSER_INTERFACE

#include <GoGraph.hpp>


/*! \class GoParserInterface
	\brief An interface class to define go graph parsers

	This class defines the interface of a go graph parser. Pure virtual methods
	  require that parsers implement these methods.

*/
class GoParserInterface {

public:

  GoParserInterface() = default;
  virtual ~GoParserInterface() = default;

	//! A pure virtual method for parsing the file and returning a GoGraph object.
	/*!
		This pure virtual method requires any parser to have a method that takes
		  a filename string and returns a GoGraph object pointer.
	*/
	[[nodiscard]] virtual std::unique_ptr<GoGraph> parseGoFile(const std::string& fileName) const = 0;

	//! A pure virtual method for parsing the file and returning a GoGraph object.
	/*!
	This pure virtual method requires any parser to have a method that takes
	a filename string and returns a GoGraph object pointer.
	*/
	[[nodiscard]] virtual bool isFileGood(const std::string& filename) const = 0;

	//! A pure virtual clone function for factory pattern
	/*!
		This pure virtual method returns an instance of this interface. This method
		  is used in a factory class to have the ability to decide the parser at runtime.
	*/
	[[nodiscard]] virtual std::unique_ptr<GoParserInterface> clone() const = 0;

};
#endif
