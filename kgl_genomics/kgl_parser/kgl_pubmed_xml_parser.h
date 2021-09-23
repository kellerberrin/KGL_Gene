//
// Created by kellerberrin on 16/8/21.
//

#ifndef KGL_PUBMED_XML_PARSER_H
#define KGL_PUBMED_XML_PARSER_H

#include "kgl_pubmed_resource.h"
#include <rapidxml.h>

namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A local XML parser implementation to assist in parsing the complex Pubmed publication XML.
// These objects contain 3rd party rapidxml implementation details and should ONLY be included into the source
// file "kgl_pubmed_api.cpp". Keep it OUT of the general code base.
// The objects are functional only and cannot be created.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse citations/references.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ParseCitationXMLImpl {

public:

  ParseCitationXMLImpl() = delete;
  ~ParseCitationXMLImpl() = delete;

  // Parse citations/references. Returns false if any parsing problems encountered.
  [[nodiscard]] static std::pair<bool, LitCitationMap> parseCitationXML(const std::string& citation_text);

private:

  // Citation/reference XML nodes
  constexpr static const char* CITATION_ROOT_NODE_{"eLinkResult"};
  constexpr static const char* CITATION_RECORD_{"LinkSet"};
  constexpr static const char* CITATION_PMID_{"IdList"};
  constexpr static const char* CITATION_LINK_DB_{"LinkSetDb"};
  constexpr static const char* CITATION_LINK_SET_{"Link"};

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse publication summaries.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ParsePublicationXMLImpl {

public:

  ParsePublicationXMLImpl() = delete;
  ~ParsePublicationXMLImpl() = delete;

  // Parse publications. Returns false if any parsing problems encountered.
  [[nodiscard]] static std::pair<bool, APIPublicationMap> parsePublicationXML(const std::string& publication_xml_text);

private:

  // High level XML tags
  constexpr static const char* PUBLICATION_ROOT_NODE_{"PubmedArticleSet"};
  constexpr static const char* PUBLICATION_NODE_{"PubmedArticle"};

  // Article XML nodes
  constexpr static const char* MEDLINE_NODE_{"MedlineCitation"};
  constexpr static const char* PMID_NODE_{"PMID"};
  constexpr static const char* ARTICLE_NODE_{"Article"};
  constexpr static const char* PUBMED_NODE_{"PubmedData"};

  // Article field nodes
  constexpr static const char* ABSTRACT_NODE_{"Abstract"};
  constexpr static const char* ABSTRACT_TEXT_NODE_{"AbstractText"};
  constexpr static const char* ARTICLE_TITLE_NODE_{"ArticleTitle"};

  // Author field nodes
  constexpr static const char* AUTHOR_LIST_NODE_{"AuthorList"};
  constexpr static const char* AUTHOR_NODE_{"Author"};
  constexpr static const char* AUTHOR_SURNAME_NODE_{"LastName"};
  constexpr static const char* AUTHOR_COLLECTIVE_NODE_{"CollectiveName"};
  constexpr static const char* AUTHOR_INITIALS_NODE_{"Initials"};

  // Journal XML nodes.
  constexpr static const char* JOURNAL_NODE_{"Journal"};
  constexpr static const char* JOURNAL_ISSN_NODE_{"ISSN"};
  constexpr static const char* JOURNAL_ISSUE_NODE_{"JournalIssue"};
  constexpr static const char* VOLUME_NODE_{"Volume"};
  constexpr static const char* ISSUE_NODE_{"Issue"};
  constexpr static const char* JOURNAL_TITLE_NODE_{"Title"};

  // Chemical XML Nodes.
  constexpr static const char* CHEMICAL_LIST_NODE_{"ChemicalList"};
  constexpr static const char* CHEMICAL_NODE_{"Chemical"};
  constexpr static const char* NAME_SUBSTANCE_NODE_{"NameOfSubstance"};

  // MeSH XML Nodes.
  constexpr static const char* MESH_LIST_NODE_{"MeshHeadingList"};
  constexpr static const char* MESH_NODE_{"MeshHeading"};
  constexpr static const char* MESH_ID_ATTRIBUTE_{"UI"};
  constexpr static const char* MESH_DESCRIPTOR_NODE_{"DescriptorName"};

  // Doi XML Nodes.
  constexpr static const char* DOI_NODE_{"ELocationID"};
  constexpr static const char* DOI_ATTRIBUTE_{"EIdType"};
  constexpr static const char* DOI_ATTRIBUTE_VALUE_{"doi"};

  // Reference XML Nodes.
  constexpr static const char* REFERENCE_LIST_NODE_{"ReferenceList"};
  constexpr static const char* REFERENCE_NODE_{"Reference"};
  constexpr static const char* ARTICLE_ID_LIST_{"ArticleIdList"};
  constexpr static const char* REFERENCE_CITATION_NODE_{"Citation"};
  constexpr static const char* ARTICLE_ID_{"ArticleId"};
  constexpr static const char* ARTICLE_ID_ATTRIBUTE_{"IdType"};
  constexpr static const char* ARTICLE_ID_ATTRIBUTE_VALUE_{"pubmed"};

  // Date XML nodes
  constexpr static const char* PUB_DATE_NODE_{"PubDate"};
  constexpr static const char* YEAR_NODE_{"Year"};
  constexpr static const char* MONTH_NODE_{"Month"};
  constexpr static const char* DAY_NODE_{"Day"};
  constexpr static const char* HISTORY_NODE_{"History"};
  constexpr static const char* PUBMED_DATE_NODE_{"PubMedPubDate"};


  // Parse a publication.
  [[nodiscard]] static PublicationSummary parsePubmedArticleXML(rapidxml::xml_node<> * pubmed_article_node);
  // Parse publication sub-sections.
  static void parseJournalArticleXML(rapidxml::xml_node<> * journal_article_node, PublicationSummary& publication_details);
  static void parseArticleFieldsXML(rapidxml::xml_node<> * journal_article_node, PublicationSummary& publication);
  static void parseDoiXML(rapidxml::xml_node<> * journal_article_node, rapidxml::xml_node<> * pubmed_node, PublicationSummary& publication);
  static void parseAuthorsXML(rapidxml::xml_node<> * journal_article_node, PublicationSummary& publication);
  static void parseChemicalsXML(rapidxml::xml_node<> * medline_node, PublicationSummary& publication);
  static void parseMeSHXML(rapidxml::xml_node<> * medline_node, PublicationSummary& publication);
  static void parseReferencesXML(rapidxml::xml_node<> * pubmed_node, PublicationSummary& publication);
  static void parseXMLDate(rapidxml::xml_node<> * journal_article_node, rapidxml::xml_node<> * pubmed_node, PublicationSummary& publication);

  // Utility routines.
  [[nodiscard]] static rapidxml::xml_node<> * validSubNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name, const std::string& pmid);
  [[nodiscard]] static rapidxml::xml_attribute<> * validAttribute(rapidxml::xml_node<> * node_ptr, const char* attribute, const std::string& pmid);
  [[nodiscard]] static std::string validOptionalNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name);
  [[nodiscard]] static std::string parseTextWithEmbeddedNodes(rapidxml::xml_node<> * node_ptr, const char* sub_node_name, const std::string& pmid);

};





} // namespace


#endif // KGL_PUBMED_XML_PARSER_H
