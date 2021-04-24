/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_ANNOTATION_PARSER_INTERFACE
#define KGL_ANNOTATION_PARSER_INTERFACE

#include "kol_AnnotationData.h"

namespace kellerberrin::ontology {


/*! \class AnnotationParserInterface
	\brief An interface class to define annotation parsers

	This class defines the interface of an annotation parser. Pure virtual methods
	  require that parsers implement these methods.

*/
class AnnotationParserInterface {
public:

  AnnotationParserInterface() = default;

  virtual ~AnnotationParserInterface() = default;

  //! A pure virtual method for parsing the file and returning an AnnotationData object.
  /*!
    This pure virtual method requires any parser to have a method that takes
      a filename string and returns an AnnotationData object pointer.
  */
  [[nodiscard]] virtual std::unique_ptr<AnnotationData> parseAnnotationFile(const std::string &fileName) const = 0;

  //! A pure virtual method for checking if a file exists and is formatted correctly.
  /*!
    This pure virtual function delegates format checking to the implementing class.
  */
  [[nodiscard]] virtual bool isFileGood(const std::string &fileName) const = 0;

  //! A pure virtual clone function for factory pattern
  /*!
    This pure virtual method returns an instance of this interface. This method
      is used in a factory class to have the ability to decide the parser at runtime.
  */
  [[nodiscard]] virtual std::unique_ptr<AnnotationParserInterface> clone() const = 0;

};

} // namespace


#endif
