//
// Created by kellerberrin on 22/9/21.
//

#ifndef KGL_LITERATURE_FILTER_VIRTUAL_H
#define KGL_LITERATURE_FILTER_VIRTUAL_H


#include <map>
#include <memory>
#include <vector>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract_ LiteratureFilter class uses the visitor pattern.
// Concrete Literature filters are defined in kgl_literature_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization level namespace


class PublicationSummary;   // Forward decl of the publication object.


class LiteratureFilter {

public:

  LiteratureFilter() = default;
  virtual ~LiteratureFilter() = default;

  [[nodiscard]] virtual bool applyFilter(const PublicationSummary& publication) const = 0;

  [[nodiscard]] virtual std::unique_ptr<LiteratureFilter> clone() const = 0;


private:


};



} // namespace


#endif //KGL_LITERATURE_FILTER_VIRTUAL_H
