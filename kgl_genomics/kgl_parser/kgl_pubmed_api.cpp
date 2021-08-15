//
// Created by kellerberrin on 11/8/21.
//

#include "kgl_pubmed_api.h"
#include "kel_exec_env.h"

#include <rapidxml.h>

#include <set>

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Citations/References for each publication.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::LitCitationMap kgl::PubmedRequester::getCitations(const std::vector<std::string>& pmid_vector) {

  return getCitationArgs(pmid_vector, PUBMED_ARTICLE_CITEDBY_ARGS_);

}

kgl::LitCitationMap kgl::PubmedRequester::getReferences(const std::vector<std::string>& pmid_vector) {

  return getCitationArgs(pmid_vector, PUBMED_ARTICLE_REFERENCES_ARGS_);

}


kgl::LitCitationMap kgl::PubmedRequester::getCitationArgs(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) {

  LitCitationMap citation_map;
  std::set<std::string> unique_pmids;

  // Allow trivial requests.
  if (pmid_vector.empty()) {

    return citation_map;

  }

  // For efficiency ensure that all pmids are unique.
  for (auto const& pmid : pmid_vector) {

    unique_pmids.insert(pmid);

  }

  std::vector<std::string> batch_pmids;
  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();
  size_t time_period_requests{0};
  for (auto const& pmid : unique_pmids) {

    batch_pmids.push_back(pmid);
    if (batch_pmids.size() >= PMID_PER_REQUEST_) {

      // Get api results.
      citation_map.merge(citationBatch(batch_pmids, cite_type_args));
      batch_pmids.clear();
      ++time_period_requests;
      if (time_period_requests >= REQUESTS_PER_SECOND_) {

        std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();
        auto batch_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch-begin_batch);

        if (batch_milliseconds < BATCH_TIME_INTERVAL_) {

          auto sleep_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(BATCH_TIME_INTERVAL_ - batch_milliseconds);
          std::this_thread::sleep_for(sleep_milliseconds);

        }

        time_period_requests = 0;
        begin_batch = std::chrono::steady_clock::now();

      }

    }

  }

  if (not batch_pmids.empty()) {

    citation_map.merge(citationBatch(batch_pmids, cite_type_args));
    batch_pmids.clear();

  }

  return citation_map;

}

kgl::LitCitationMap kgl::PubmedRequester::citationBatch(const std::vector<std::string>& pmid_vector, const std::string& cite_type_args) {

  LitCitationMap citation_map;
  std::string request_string;
  const std::string id_prefix{"&id="};
  for (auto const& pmid : pmid_vector) {

    request_string += id_prefix + pmid;

  }

  std::string citation_request = cite_type_args + request_string + API_KEY;

  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();

  auto [citation_result, citation_text] = pubmed_rest_api_.synchronousRequest(PUBMED_ELINK_URL_, citation_request);

  std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();

  if (not citation_result) {

    ExecEnv::log().error("PubmedRequester::getCitationArgs; problem with Pubmed API: {}", citation_text);

  } else {

    citation_map = parseCitationXML(citation_text);

  }

  auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch - begin_batch);

  ExecEnv::log().info("Batch Execution time: {} (ms)", batch_duration.count());

  return citation_map;

}


kgl::LitCitationMap kgl::PubmedRequester::parseCitationXML(const std::string& citation_text) {

  LitCitationMap citation_map;

  // All this raw pointer stuff is very nasty, is rapidxml the right library?
  rapidxml::xml_document<> doc;
  std::unique_ptr<char[]> text_buffer_ptr(std::make_unique<char[]>(citation_text.size() + 1));
  std::memcpy(text_buffer_ptr.get(), citation_text.c_str(), citation_text.size() + 1);

  // Parse the buffer using the xml file parsing library into doc
  doc.parse<0>(text_buffer_ptr.get());

  rapidxml::xml_node<>* root_node = doc.first_node(CITATION_ROOT_NODE_.c_str());

  if (root_node != nullptr) {

    rapidxml::xml_node<> * citation_node = root_node->first_node(CITATION_RECORD_.c_str());
    while (citation_node != nullptr)
    {

      auto citation_pmid = citation_node->first_node(CITATION_PMID_.c_str());
      if (citation_pmid != nullptr) {

        auto pmid_node = citation_pmid->first_node();
        if (pmid_node != nullptr) {

          std::vector<std::string> citation_pmid_vector;
          auto link_db = citation_node->first_node(CITATION_LINK_DB_.c_str());
          if (link_db != nullptr) {

            auto link_node = link_db->first_node(CITATION_LINK_SET_.c_str());
            while (link_node != nullptr) {

              auto cite_pmid_node = link_node->first_node();
              if (cite_pmid_node != nullptr) {

                citation_pmid_vector.emplace_back(cite_pmid_node->value());

              } else {

                ExecEnv::log().error("PubmedRequester::parseCitationXML; No citation node found");

              }

              link_node = link_node->next_sibling();

            }

          }

          auto [iter, result] = citation_map.try_emplace(pmid_node->value(), citation_pmid_vector);
          if (not result) {

            ExecEnv::log().error("PubmedRequester::parseCitationXML; expected duplicate for pmid_: {}", pmid_node->value());

          }

        } else {

          ExecEnv::log().error("PubmedRequester::parseCitationXML; Pmid Id Attribute: {} does not exist", CITATION_PMID_);

        }

      } else {

        ExecEnv::log().error("PubmedRequester::parseCitationXML; Pmid Id Attribute: {} does not exist", CITATION_PMID_);

      }

      citation_node = citation_node->next_sibling();

    }

  } else  {

    ExecEnv::log().error("PubmedRequester::parseCitationXML; Citation root node: {} is null", CITATION_ROOT_NODE_);

  }

  return citation_map;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Publication details.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::LitPublicationMap kgl::PubmedRequester::getPublicationDetails(const std::vector<std::string>& pmid_vector) {

  LitPublicationMap publication_map;
  std::set<std::string> unique_pmids;

  // Allow trivial requests.
  if (pmid_vector.empty()) {

    return publication_map;

  }

  // For efficiency ensure that all pmids are unique.
  for (auto const& pmid : pmid_vector) {

    unique_pmids.insert(pmid);

  }

  std::vector<std::string> batch_pmids;
  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();
  size_t time_period_requests{0};
  for (auto const& pmid : unique_pmids) {

    batch_pmids.push_back(pmid);
    if (batch_pmids.size() >= PMID_PER_REQUEST_) {

      // Get api results.
      publication_map.merge(publicationBatch(batch_pmids));
      batch_pmids.clear();
      ++time_period_requests;
      if (time_period_requests >= REQUESTS_PER_SECOND_) {

        std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();
        auto batch_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch-begin_batch);

        if (batch_milliseconds < BATCH_TIME_INTERVAL_) {

          auto sleep_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(BATCH_TIME_INTERVAL_ - batch_milliseconds);
          std::this_thread::sleep_for(sleep_milliseconds);

        }

        time_period_requests = 0;
        begin_batch = std::chrono::steady_clock::now();

      }

    }

  }

  if (not batch_pmids.empty()) {

    publication_map.merge(publicationBatch(batch_pmids));
    batch_pmids.clear();

  }

  return publication_map;

}


kgl::LitPublicationMap kgl::PubmedRequester::publicationBatch(const std::vector<std::string>& pmid_vector) {

  LitPublicationMap publication_map;
  std::string request_string;
  const std::string id_prefix{"&id="};
  for (auto const& pmid : pmid_vector) {

    request_string += id_prefix + pmid;

  }

  std::string citation_request = PUBMED_ARTICLE_DETAIL_ARGS_ + request_string + API_KEY;

  std::chrono::steady_clock::time_point begin_batch = std::chrono::steady_clock::now();

  auto [publication_result, publication_text] = pubmed_rest_api_.synchronousRequest(PUBMED_EFETCH_URL_, citation_request);

  std::chrono::steady_clock::time_point end_batch = std::chrono::steady_clock::now();

  if (not publication_result) {

    ExecEnv::log().error("PubmedRequester::publicationBatch; problem with Pubmed API: {}", publication_text);

  } else {

    ExecEnv::log().info("Publication text: {}", publication_text);
    publication_map = parsePublicationXML(publication_text);

  }

  auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_batch - begin_batch);

  ExecEnv::log().info("Batch Execution time: {} (ms)", batch_duration.count());

  return publication_map;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A local XML parser implementation to assist in parsing the complex Pubmed publication XML.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ParsePublicationXMLImpl {

public:

  ParsePublicationXMLImpl() = delete;
  ~ParsePublicationXMLImpl() = delete;

  static kgl::PubMedPublicationSummary parsePubmedArticleXML(rapidxml::xml_node<> * pubmed_article_node);


private:

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
  constexpr static const char* AUTHOR_INITIALS_NODE_{"Initials"};

  // Journal XML nodes.
  constexpr static const char* JOURNAL_NODE_{"Journal"};
  constexpr static const char* JOURNAL_ISSUE_NODE_{"JournalIssue"};
  constexpr static const char* VOLUME_NODE_{"Volume"};
  constexpr static const char* ISSUE_NODE_{"Issue"};
  constexpr static const char* PUB_DATE_NODE_{"PubDate"};
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

  // Reference XML Nodes.
  constexpr static const char* REFERENCE_LIST_NODE_{"ReferenceList"};
  constexpr static const char* REFERENCE_NODE_{"Reference"};
  constexpr static const char* ARTICLE_ID_LIST_{"ArticleIdList"};
  constexpr static const char* ARTICLE_ID_{"ArticleId"};
  constexpr static const char* ARTICLE_ID_ATTRIBUTE_{"IdType"};
  constexpr static const char* ARTICLE_ID_ATTRIBUTE_VALUE_{"pubmed"};


  // Date XML nodes
  constexpr static const char* YEAR_NODE_{"Year"};
  constexpr static const char* MONTH_NODE_{"Month"};
  constexpr static const char* DAY_NODE_{"Day"};

  // Parse publication sub-sections.
  static void parseJournalArticleXML(rapidxml::xml_node<> * journal_article_node, kgl::PubMedPublicationSummary& publication_details);
  static void parseArticleFieldsXML(rapidxml::xml_node<> * journal_article_node, kgl::PubMedPublicationSummary& publication);
  static void parseAuthorsXML(rapidxml::xml_node<> * journal_article_node, kgl::PubMedPublicationSummary& publication);
  static void parseChemicalsXML( rapidxml::xml_node<> * medline_node, kgl::PubMedPublicationSummary& publication);
  static void parseMeSHXML( rapidxml::xml_node<> * medline_node, kgl::PubMedPublicationSummary& publication);
  static void parseReferencesXML( rapidxml::xml_node<> * pubmed_node, kgl::PubMedPublicationSummary& publication);

  // Utility routines.
  [[nodiscard]] static std::string parseXMLDate(rapidxml::xml_node<> * date_node, const std::string& pmid);
  [[nodiscard]] static rapidxml::xml_node<> * validSubNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name, const std::string& pmid);
  [[nodiscard]] static rapidxml::xml_attribute<> * validAttribute(rapidxml::xml_node<> * node_ptr, const char* attribute, const std::string& pmid);
  [[nodiscard]] static std::string validOptionalNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name);

};



kgl::LitPublicationMap kgl::PubmedRequester::parsePublicationXML(const std::string& publication_xml_text) {

  LitPublicationMap publication_map;

  // Normally I avoid try/catch blocks but these XML documents are just too complex for if/then error handling.
  try {

    rapidxml::xml_document<> doc;
    std::unique_ptr<char[]> text_buffer_ptr(std::make_unique<char[]>(publication_xml_text.size() + 1));
    std::memcpy(text_buffer_ptr.get(), publication_xml_text.c_str(), publication_xml_text.size() + 1);

    // Parse the buffer using the xml file parsing library into doc
    doc.parse<0>(text_buffer_ptr.get());

    rapidxml::xml_node<>* root_node = doc.first_node(PUBLICATION_ROOT_NODE_.c_str());
    if (root_node == nullptr) {

      std::string error_message = "Publication Root Node: " + PUBLICATION_ROOT_NODE_ + " not found";
      throw std::runtime_error(error_message);

    }

    // Loop through the articles.
    rapidxml::xml_node<> * article_node = root_node->first_node(PUBLICATION_NODE_.c_str());
    while (article_node != nullptr) {

      auto publication = ParsePublicationXMLImpl::parsePubmedArticleXML(article_node);
      if (publication.pmid_.empty()) {

        std::string error_message = "publication Pubmed pmid_ not defined";
        throw std::runtime_error(error_message);

      }

      auto [iter, result] = publication_map.try_emplace(publication.pmid_, publication);
      if (not result) {

        std::string error_message = "cannot add duplicate publication: " + publication.pmid_;
        throw std::runtime_error(error_message);

      }

      // Next article.
      article_node = article_node->next_sibling();

    }

    ExecEnv::log().info("PubmedRequester::parsePublicationXML; article count {}", publication_map.size());

  } catch(std::exception& e) {

    ExecEnv::log().error("PubmedRequester::parsePublicationXML; error parsing XML Pubmed Article; {}", e.what());
    publication_map.clear();

  }

  return publication_map;

}


kgl::PubMedPublicationSummary ParsePublicationXMLImpl::parsePubmedArticleXML(rapidxml::xml_node<> * pubmed_article_node) {

  kgl::PubMedPublicationSummary publication;

  rapidxml::xml_node<>* medline_node = pubmed_article_node->first_node(MEDLINE_NODE_);
  if (medline_node == nullptr) {

    std::string error_message = "Publication Node: " + std::string(MEDLINE_NODE_) + " not found";
    throw std::runtime_error(error_message);

  }

  rapidxml::xml_node<>* pmid_node = medline_node->first_node(PMID_NODE_);
  if (pmid_node == nullptr) {

    std::string error_message = "Publication Node: " + std::string(PMID_NODE_) + " not found";
    throw std::runtime_error(error_message);

  }


  std::string pmid = pmid_node->value();
  publication.pmid_ = pmid;

  parseChemicalsXML( medline_node, publication);
  parseMeSHXML( medline_node, publication);

  auto journal_article_node = validSubNode(medline_node, ARTICLE_NODE_, publication.pmid_);
  parseArticleFieldsXML(journal_article_node, publication);
  parseJournalArticleXML(journal_article_node, publication);
  parseAuthorsXML(journal_article_node, publication);

  rapidxml::xml_node<>* pubmed_node = pubmed_article_node->first_node(PUBMED_NODE_);
  if (pubmed_node == nullptr) {

    std::string error_message = "Publication Node: " + std::string(PUBMED_NODE_) + " not found";
    throw std::runtime_error(error_message);

  }

  parseReferencesXML(pubmed_node, publication);

  return publication;

}

// Unpack the various XML tags for the article.
void ParsePublicationXMLImpl::parseAuthorsXML( rapidxml::xml_node<> * journal_article_node,
                                               kgl::PubMedPublicationSummary& publication) {

  auto author_list_node = validSubNode(journal_article_node, AUTHOR_LIST_NODE_, publication.pmid_);
  auto author_node = author_list_node->first_node(AUTHOR_NODE_);
  while(author_node != nullptr)  {

    auto surname = validSubNode(author_node, AUTHOR_SURNAME_NODE_, publication.pmid_);
    auto initials = validSubNode(author_node, AUTHOR_INITIALS_NODE_, publication.pmid_);

    std::string author_name = surname->value() + std::string(", ") + initials->value();
    publication.authors_.push_back(author_name);

    // Next author.
    author_node = author_node->next_sibling();

  }

}


// Unpack the various XML tags for the article.
void ParsePublicationXMLImpl::parseArticleFieldsXML(rapidxml::xml_node<> * journal_article_node,
                                                    kgl::PubMedPublicationSummary& publication) {

  auto abstract_node = validSubNode(journal_article_node, ABSTRACT_NODE_, publication.pmid_);

  auto abtract_text = validSubNode(abstract_node, ABSTRACT_TEXT_NODE_, publication.pmid_);

  publication.abstract_ = abtract_text->value();

  auto article_title = validSubNode(journal_article_node, ARTICLE_TITLE_NODE_, publication.pmid_);

  publication.title_ = article_title->value();

}


// Unpack the various XML tags to define the Publication Journal
void ParsePublicationXMLImpl::parseJournalArticleXML(rapidxml::xml_node<> * journal_article_node,
                                                     kgl::PubMedPublicationSummary& publication) {

  auto journal_node = validSubNode(journal_article_node, JOURNAL_NODE_, publication.pmid_);

  auto journal_title = validSubNode(journal_node, JOURNAL_TITLE_NODE_, publication.pmid_);

  publication.journal_ = journal_title->value();

  auto journal_issue = validSubNode(journal_node, JOURNAL_ISSUE_NODE_, publication.pmid_);

  publication.journal_issue_ = validOptionalNode(journal_issue, ISSUE_NODE_);

  publication.journal_volume_ = validOptionalNode(journal_issue, VOLUME_NODE_);

  auto pub_date = validSubNode(journal_issue, PUB_DATE_NODE_, publication.pmid_);

  publication.publication_date_ = parseXMLDate(pub_date, publication.pmid_);

}

// Unpack the various XML tags for chemicals_.
void ParsePublicationXMLImpl::parseChemicalsXML( rapidxml::xml_node<> * medline_node,
                                                 kgl::PubMedPublicationSummary& publication) {

  auto chem_list_node = medline_node->first_node(CHEMICAL_LIST_NODE_);
  if (chem_list_node == nullptr) {

    return;

  }
  auto chem_node = chem_list_node->first_node(CHEMICAL_NODE_);
  while(chem_node != nullptr)  {

    auto desc_node = validSubNode(chem_node, NAME_SUBSTANCE_NODE_, publication.pmid_);

    std::string chemical_description = desc_node->value();
    auto mesh_attr = validAttribute(desc_node, MESH_ID_ATTRIBUTE_, publication.pmid_);
    std::string mesh_id = mesh_attr->name();

    publication.chemicals_.emplace_back(mesh_id, chemical_description);

    // Next chemical.
    chem_node = chem_node->next_sibling();

  }

}


// Unpack the various XML tags for MeSH descriptors.
void ParsePublicationXMLImpl::parseMeSHXML( rapidxml::xml_node<> * medline_node,
                                            kgl::PubMedPublicationSummary& publication) {

  auto mesh_list_node = medline_node->first_node(MESH_LIST_NODE_);
  if (mesh_list_node == nullptr) {

    return;

  }
  auto mesh_node = mesh_list_node->first_node(MESH_NODE_);
  while(mesh_node != nullptr)  {

    auto desc_node = validSubNode(mesh_node, MESH_DESCRIPTOR_NODE_, publication.pmid_);

    std::string mesh_description = desc_node->value();
    auto mesh_attr = validAttribute(desc_node, MESH_ID_ATTRIBUTE_, publication.pmid_);
    std::string mesh_id = mesh_attr->name();

    publication.MeSHcodes_.emplace_back(mesh_id, mesh_description);

    // Next chemical.
    mesh_node = mesh_node->next_sibling();

  }

}


void ParsePublicationXMLImpl::parseReferencesXML( rapidxml::xml_node<> * pubmed_node,
                                                  kgl::PubMedPublicationSummary& publication) {

  auto reference_list_node = pubmed_node->first_node(REFERENCE_LIST_NODE_);
  if (reference_list_node == nullptr) {

    return;

  }
  auto reference_node = reference_list_node->first_node(REFERENCE_NODE_);
  while(reference_node != nullptr)  {

    auto article_list_node = validSubNode(reference_node, ARTICLE_ID_LIST_, publication.pmid_);
    auto id_node = validSubNode(article_list_node, ARTICLE_ID_, publication.pmid_);
    auto id_attrib = validAttribute(id_node, ARTICLE_ID_ATTRIBUTE_, publication.pmid_);
    std::string id_type = id_attrib->value();
    if (id_type != ARTICLE_ID_ATTRIBUTE_VALUE_) {

      kel::ExecEnv::log().warn("ParsePublicationXMLImpl::parseReferencesXML; publication pmid_: {} has reference id type: {}",
                               publication.pmid_, id_type);

    }

    std::string reference_pmid = id_node->value();

    auto [iter, result] = publication.references_.insert(reference_pmid);
    if (not result) {

      kel::ExecEnv::log().warn("ParsePublicationXMLImpl::parseReferencesXML; publication pmid_: {} has duplicate reference: {}",
                               publication.pmid_, reference_pmid);

    }

    // Next reference.
    reference_node = reference_node->next_sibling();

  }

}


std::string ParsePublicationXMLImpl::parseXMLDate(rapidxml::xml_node<> * date_node, const std::string& pmid) {

  auto year = validSubNode(date_node, YEAR_NODE_, pmid);
  auto month = validSubNode(date_node, MONTH_NODE_, pmid);
  auto day =  validOptionalNode(date_node, DAY_NODE_);

  std::string date_string;
  if (not day.empty()) {

    date_string = day + "-" + std::string(month->value()) + "-" + std::string(year->value());

  } else {

    date_string = std::string(month->value()) + "-" + std::string(year->value());

  }

  return date_string;

}

rapidxml::xml_node<> * ParsePublicationXMLImpl::validSubNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name, const std::string& pmid) {

  rapidxml::xml_node<>* sub_node_ptr = node_ptr->first_node(sub_node_name);
  if (sub_node_ptr == nullptr) {

    std::string error_message = "Pmid: " + pmid + " Node: '" + sub_node_name + "' not found";
    throw std::runtime_error(error_message);

  }

  return sub_node_ptr;

}

rapidxml::xml_attribute<> * ParsePublicationXMLImpl::validAttribute(rapidxml::xml_node<> * node_ptr, const char* attribute, const std::string& pmid) {

  rapidxml::xml_attribute<>* attrib_ptr = node_ptr->first_attribute(attribute);
  if (attrib_ptr == nullptr) {

    std::string error_message = "Pmid: " + pmid + " Attribute: '" + attribute + "' not found";
    throw std::runtime_error(error_message);

  }

  return attrib_ptr;

}


std::string ParsePublicationXMLImpl::validOptionalNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name) {

  rapidxml::xml_node<>* sub_node_ptr = node_ptr->first_node(sub_node_name);
  if (sub_node_ptr == nullptr) {

    return {""};

  }

  return sub_node_ptr->value();

}
