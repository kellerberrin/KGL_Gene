//
// Created by kellerberrin on 16/8/21.
//

#include "kgl_pubmed_xml_parser.h"
#include "kel_exec_env.h"

#include <cstring>
#include <thread>

namespace kgl = kellerberrin::genome;



std::pair<bool, kgl::LitCitationMap> kgl::ParseCitationXMLImpl::parseCitationXML(const std::string& citation_text, const DateGP& download_date) {

  LitCitationMap citation_map;
  bool parse_result{true};


  try {


    // All this raw pointer stuff is very nasty, is rapidxml the right library?
    rapidxml::xml_document<> doc;
    std::unique_ptr<char[]> text_buffer_ptr(std::make_unique<char[]>(citation_text.size() + 1));
    std::memcpy(text_buffer_ptr.get(), citation_text.c_str(), citation_text.size() + 1);

    // Parse the buffer using the xml file parsing library into doc
    doc.parse<0>(text_buffer_ptr.get());

    rapidxml::xml_node<>* root_node = doc.first_node(CITATION_ROOT_NODE_);
    if (root_node == nullptr) {

      ExecEnv::log().error("PubmedRequester::parseCitationXML; Citation root node: {} is null", CITATION_ROOT_NODE_);
      return {false, citation_map};

    }

    rapidxml::xml_node<> * citation_node = root_node->first_node(CITATION_RECORD_);
    while (citation_node != nullptr)
    {

      auto citation_pmid = citation_node->first_node(CITATION_PMID_);
      if (citation_pmid != nullptr) {

        auto pmid_node = citation_pmid->first_node();
        if (pmid_node != nullptr) {

          std::set<std::string> citation_pmid_set;
          auto link_db = citation_node->first_node(CITATION_LINK_DB_);
          if (link_db != nullptr) {

            auto link_node = link_db->first_node(CITATION_LINK_SET_);
            while (link_node != nullptr) {

              auto cite_pmid_node = link_node->first_node();
              if (cite_pmid_node != nullptr) {

                auto [iter, result] = citation_pmid_set.insert(cite_pmid_node->value());
                if (not result) {

                  ExecEnv::log().error("PubmedRequester::parseCitationXML; Cannot add duplicate citation: {}", cite_pmid_node->value());

                }

              } else {

                ExecEnv::log().error("PubmedRequester::parseCitationXML; No citation node found");
                parse_result = false;

              }

              link_node = link_node->next_sibling();

            }

          }

          auto [iter, result] = citation_map.try_emplace(pmid_node->value(), std::pair<DateGP, std::set<std::string>>{download_date, citation_pmid_set});
          if (not result) {

            ExecEnv::log().error("PubmedRequester::parseCitationXML; expected duplicate for pmid_: {}", pmid_node->value());
            parse_result = false;

          }

        } else {

          ExecEnv::log().error("PubmedRequester::parseCitationXML; Pmid Id Attribute: {} does not exist", CITATION_PMID_);
          parse_result = false;

        }

      } else {

        ExecEnv::log().error("PubmedRequester::parseCitationXML; Pmid Id Attribute: {} does not exist", CITATION_PMID_);
        parse_result = false;

      }

      citation_node = citation_node->next_sibling();

    }

  } catch(const std::exception& e) {

    ExecEnv::log().error("PubmedRequester::parseCitationXML; RapidXML Fails: {}", e.what());
    ExecEnv::log().info("PubmedRequester::parseCitationXML; Text: \n{}", citation_text);
    return {false, LitCitationMap() };

  }

  return {parse_result, citation_map};

}



std::pair<bool, kgl::APIPublicationMap> kgl::ParsePublicationXMLImpl::parsePublicationXML( const std::string& publication_xml_text,
                                                                                           const DateGP& download_date) {

  // Comment out after testing.
  // ExecEnv::log().info("PubmedRequester::parsePublicationXML; text:\n{}", publication_xml_text);

  APIPublicationMap publication_map;
  bool parse_result{true};

  try {

    // Place into a char buffer as rapidxml modifies the contents of the buffer.
    rapidxml::xml_document<> doc;
    std::unique_ptr<char[]> text_buffer_ptr(std::make_unique<char[]>(publication_xml_text.size() + 1));
    std::memcpy(text_buffer_ptr.get(), publication_xml_text.c_str(), publication_xml_text.size() + 1);

    // Parse the char buffer,
    doc.parse<0>(text_buffer_ptr.get());

    // Get the root node.
    rapidxml::xml_node<>* root_node = doc.first_node(PUBLICATION_ROOT_NODE_);
    if (root_node == nullptr) {

      ExecEnv::log().error("PubmedRequester::parsePublicationXML; error parsing publication Root Node: {}", PUBLICATION_ROOT_NODE_);
      return {false, publication_map}; // ReturnType the empty map.

    }

    // Get the first article and loop through the articles.
    rapidxml::xml_node<> * article_node = root_node->first_node(PUBLICATION_NODE_);
    while (article_node != nullptr) {

      // Parse the article in try/catch block to manage the complex error XML error handling.
      try {

        auto publication = ParsePublicationXMLImpl::parsePubmedArticleXML(article_node);
        if (publication.pmid().empty()) {

          std::string error_message = "publication Pubmed pmid not defined";
          throw std::runtime_error(error_message);

        }

        publication.downloadDate(download_date);

        auto [iter, result] = publication_map.try_emplace(publication.pmid(), std::make_shared<PublicationSummary>(publication));
        if (not result) {

          std::string error_message = "cannot add duplicate publication: " + publication.pmid();
          throw std::runtime_error(error_message);

        }

      } catch(std::exception& e) {

        // Issue an error message and skip to the next article.
        ExecEnv::log().error("PubmedRequester::parsePublicationXML; error parsing XML Pubmed Article; {}", e.what());
        parse_result = false;


      }

      // Next article.
      article_node = article_node->next_sibling();

    }

  } catch(const std::exception& e) {

    ExecEnv::log().error("PubmedRequester::parsePublicationXML; RapidXML Fails: {}", e.what());
    ExecEnv::log().info("PubmedRequester::parsePublicationXMLL; Text: \n{}", publication_xml_text);
    return {false, APIPublicationMap() };

  }



  return {parse_result, publication_map};

}


kgl::PublicationSummary kgl::ParsePublicationXMLImpl::parsePubmedArticleXML(rapidxml::xml_node<> * pubmed_article_node) {


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
  kgl::PublicationSummary publication(pmid);

  parseChemicalsXML( medline_node, publication);
  parseMeSHXML( medline_node, publication);

  auto journal_article_node = validSubNode(medline_node, ARTICLE_NODE_, publication.pmid());

  parseArticleFieldsXML(journal_article_node, publication);
  parseJournalArticleXML(journal_article_node, publication);
  parseAuthorsXML(journal_article_node, publication);

  rapidxml::xml_node<>* pubmed_node = pubmed_article_node->first_node(PUBMED_NODE_);
  if (pubmed_node == nullptr) {

    std::string error_message = "Publication Node: " + std::string(PUBMED_NODE_) + " not found";
    throw std::runtime_error(error_message);

  }

  parseReferencesXML(pubmed_node, publication);
  parseDoiXML(journal_article_node, pubmed_node, publication);
  parseXMLDate(journal_article_node, pubmed_node, publication);

  return publication;

}

// Unpack the various XML tags for the article.
void kgl::ParsePublicationXMLImpl::parseAuthorsXML(rapidxml::xml_node<> * journal_article_node,
                                                   PublicationSummary& publication) {

//  auto author_list_node = validSubNode(journal_article_node, AUTHOR_LIST_NODE_, publication.pmid());
  auto author_list_node = journal_article_node->first_node(AUTHOR_LIST_NODE_);
  if (author_list_node == nullptr) {

    return;

  }
  auto author_node = author_list_node->first_node(AUTHOR_NODE_);
  while(author_node != nullptr)  {

    auto surname = author_node->first_node(AUTHOR_SURNAME_NODE_);
    if (surname != nullptr) {

      std::string author_initials;
      auto initials = author_node->first_node(AUTHOR_INITIALS_NODE_);
      if (initials != nullptr) {

        author_initials = initials->value();

      }

      publication.addAuthor(surname->value(), author_initials);

    } else {

      auto collective = validSubNode(author_node, AUTHOR_COLLECTIVE_NODE_, publication.pmid());
      publication.addAuthor(collective->value(), "");

    }

    // Next author.
    author_node = author_node->next_sibling();

  }

}


// Unpack the various XML tags for the article.
void kgl::ParsePublicationXMLImpl::parseArticleFieldsXML(rapidxml::xml_node<> * journal_article_node,
                                                         PublicationSummary& publication) {

  auto abstract_node = journal_article_node->first_node(ABSTRACT_NODE_);
  if (abstract_node != nullptr) {

    auto abstract = parseTextWithEmbeddedNodes(abstract_node, ABSTRACT_TEXT_NODE_, publication.pmid());
    publication.abstract(abstract);

  } else {

    publication.abstract("");

  }

  auto title = parseTextWithEmbeddedNodes(journal_article_node, ARTICLE_TITLE_NODE_, publication.pmid());
  publication.title(title);

}

void kgl::ParsePublicationXMLImpl::parseDoiXML(rapidxml::xml_node<> * journal_article_node,
                                               rapidxml::xml_node<> * pubmed_node,
                                               PublicationSummary& publication) {


  rapidxml::xml_node<>* doi_node = journal_article_node->first_node(DOI_NODE_);
  if (doi_node != nullptr) {

    auto doi_attr = validAttribute(doi_node, DOI_ATTRIBUTE_, publication.pmid());
    if (std::string(doi_attr->value()) == std::string(DOI_ATTRIBUTE_VALUE_)) {

      publication.doi(doi_node->value());
      return;

    }

  }

  rapidxml::xml_node<>* article_id_node = pubmed_node->first_node(ARTICLE_ID_LIST_);
  if (article_id_node != nullptr) {

    auto article_id = article_id_node->first_node(ARTICLE_ID_);
    while(article_id != nullptr) {

      auto doi_attr = validAttribute(article_id, ARTICLE_ID_ATTRIBUTE_, publication.pmid());
      if (std::string(doi_attr->value()) == std::string(DOI_ATTRIBUTE_VALUE_)) {

        publication.doi(article_id->value());

      }

      article_id = article_id->next_sibling();

    }

  }

}


// Unpack the various XML tags to define the Publication Journal
void kgl::ParsePublicationXMLImpl::parseJournalArticleXML(rapidxml::xml_node<> * journal_article_node,
                                                          PublicationSummary& publication) {

  auto journal_node = validSubNode(journal_article_node, JOURNAL_NODE_, publication.pmid());

  auto journal_title = validSubNode(journal_node, JOURNAL_TITLE_NODE_, publication.pmid());

  publication.journal(journal_title->value());

  publication.journalISSN(validOptionalNode(journal_node, JOURNAL_ISSN_NODE_));

  auto journal_issue = validSubNode(journal_node, JOURNAL_ISSUE_NODE_, publication.pmid());

  publication.journalIssue(validOptionalNode(journal_issue, ISSUE_NODE_));

  publication.journalVolume(validOptionalNode(journal_issue, VOLUME_NODE_));


}

// Unpack the various XML tags for chemicals_.
void kgl::ParsePublicationXMLImpl::parseChemicalsXML(rapidxml::xml_node<> * medline_node,
                                                     PublicationSummary& publication) {

  auto chem_list_node = medline_node->first_node(CHEMICAL_LIST_NODE_);
  if (chem_list_node == nullptr) {

    return;

  }
  auto chem_node = chem_list_node->first_node(CHEMICAL_NODE_);
  while(chem_node != nullptr)  {

    auto desc_node = validSubNode(chem_node, NAME_SUBSTANCE_NODE_, publication.pmid());

    std::string chemical_description = desc_node->value();
    auto mesh_attr = validAttribute(desc_node, MESH_ID_ATTRIBUTE_, publication.pmid());
    std::string mesh_id = mesh_attr->value();

    publication.addChemical(mesh_id, chemical_description);

    // Next chemical.
    chem_node = chem_node->next_sibling();

  }

}


// Unpack the various XML tags for MeSH descriptors.
void kgl::ParsePublicationXMLImpl::parseMeSHXML( rapidxml::xml_node<> * medline_node,
                                                 kgl::PublicationSummary& publication) {

  auto mesh_list_node = medline_node->first_node(MESH_LIST_NODE_);
  if (mesh_list_node == nullptr) {

    return;

  }
  auto mesh_node = mesh_list_node->first_node(MESH_NODE_);
  while(mesh_node != nullptr)  {

    auto desc_node = validSubNode(mesh_node, MESH_DESCRIPTOR_NODE_, publication.pmid());

    std::string mesh_description = desc_node->value();
    auto mesh_attr = validAttribute(desc_node, MESH_ID_ATTRIBUTE_, publication.pmid());
    std::string mesh_id = mesh_attr->value();

    publication.addMeSHCode(mesh_id, mesh_description);

    // Next chemical.
    mesh_node = mesh_node->next_sibling();

  }

}


void kgl::ParsePublicationXMLImpl::parseReferencesXML(rapidxml::xml_node<> * pubmed_node,
                                                      PublicationSummary& publication) {

  auto reference_list_node = pubmed_node->first_node(REFERENCE_LIST_NODE_);
  if (reference_list_node == nullptr) {

    return;

  }
  auto reference_node = reference_list_node->first_node(REFERENCE_NODE_);
  while(reference_node != nullptr)  {

    std::string reference_pmid;
    rapidxml::xml_node<>* article_list_node = reference_node->first_node(ARTICLE_ID_LIST_);
    if (article_list_node != nullptr) {

      rapidxml::xml_node<>* id_node = article_list_node->first_node(ARTICLE_ID_);
      while (id_node != nullptr) {

        auto id_attrib = validAttribute(id_node, ARTICLE_ID_ATTRIBUTE_, publication.pmid());
        std::string id_type = id_attrib->value();
        if (id_type == ARTICLE_ID_ATTRIBUTE_VALUE_) {

          reference_pmid = id_node->value();
          break;

        }

        id_node = id_node->next_sibling(ARTICLE_ID_);

      }

    }

    std::string citation_text;
    rapidxml::xml_node<>* citation_node = reference_node->first_node(REFERENCE_CITATION_NODE_);
    if (citation_node != nullptr) {

      citation_text = citation_node->value();

    }

    if (not reference_pmid.empty() or not citation_text.empty()) {

      publication.addReference(reference_pmid, citation_text);

    } else {

      ExecEnv::log().warn("ParsePublicationXMLImpl::parseReferencesXML; reference publication pmid and citation text both empty",
                          publication.pmid(), reference_pmid);

    }

    // Next reference.
    reference_node = reference_node->next_sibling();

  }

}

void kgl::ParsePublicationXMLImpl::parseXMLDate(rapidxml::xml_node<> * journal_article_node,
                                                rapidxml::xml_node<> * pubmed_node,
                                                PublicationSummary& publication) {


  auto journal_node = validSubNode(journal_article_node, JOURNAL_NODE_, publication.pmid());
  auto journal_issue = validSubNode(journal_node, JOURNAL_ISSUE_NODE_, publication.pmid());
  auto pub_date = validSubNode(journal_issue, PUB_DATE_NODE_, publication.pmid());

  auto year_node = pub_date->first_node(YEAR_NODE_);
  auto month_node = pub_date->first_node(MONTH_NODE_);

  if (year_node != nullptr and month_node != nullptr) {

    auto day_node = pub_date->first_node(DAY_NODE_);
    if (day_node != nullptr) {

      std::string date_string = std::string(year_node->value()) + "-" + std::string(month_node->value()) + "-" + std::string(day_node->value());
      publication.publicationDate(DateGP(date_string));

    }else {

      std::string date_string = std::string(year_node->value()) + "-" + std::string(month_node->value()) + "-01";
      publication.publicationDate(DateGP(date_string));

    }

  } else {

    auto history_node = validSubNode(pubmed_node, HISTORY_NODE_, publication.pmid());
    auto pubmed_date_node = validSubNode(history_node, PUBMED_DATE_NODE_, publication.pmid());

    auto pm_year = pubmed_date_node->first_node(YEAR_NODE_);
    auto pm_month = pubmed_date_node->first_node(MONTH_NODE_);

    if (pm_year != nullptr and pm_month != nullptr) {

      auto pm_day = pubmed_date_node->first_node(DAY_NODE_);
      if (pm_day != nullptr) {

        std::string date_string = std::string(pm_year->value()) + "-" + std::string(pm_month->value()) + "-" + std::string(pm_day->value());
        publication.publicationDate(DateGP(date_string));

      }else {

        std::string date_string = std::string(pm_year->value()) + "-" + std::string(pm_month->value()) + "-01";
        publication.publicationDate(DateGP(date_string));

      }


    } else {

      if (year_node != nullptr) {

        std::string date_string = std::string(year_node->value()) + "-Jun-01";
        publication.publicationDate(DateGP(date_string));

      } else if (pm_year != nullptr) {

        std::string date_string = std::string(pm_year->value()) + "-Jun-01";
        publication.publicationDate(DateGP(date_string));

      } else {

        ExecEnv::log().warn("ParsePublicationXMLImpl::parseXMLDate; could not find a publication date for pmid: {}", publication.pmid());

      }

    }

  }

}

std::string kgl::ParsePublicationXMLImpl::parseTextWithEmbeddedNodes(rapidxml::xml_node<> * node_ptr,
                                                                     const char* sub_node_name,
                                                                     const std::string& pmid) {

  std::string sub_node_text;

  auto sub_node_ptr = validSubNode( node_ptr, sub_node_name, pmid);

  auto child_node_ptr = sub_node_ptr->first_node();

  while (child_node_ptr != nullptr) {

    sub_node_text += child_node_ptr->value();
    child_node_ptr = child_node_ptr->next_sibling();

  }

  return sub_node_text;

}


rapidxml::xml_node<> * kgl::ParsePublicationXMLImpl::validSubNode( rapidxml::xml_node<> * node_ptr,
                                                                   const char* sub_node_name,
                                                                   const std::string& pmid) {

  rapidxml::xml_node<>* sub_node_ptr = node_ptr->first_node(sub_node_name);
  if (sub_node_ptr == nullptr) {

    std::string error_message = "Pmid: " + pmid + " Node: '" + sub_node_name + "' not found";
    throw std::runtime_error(error_message);

  }

  return sub_node_ptr;

}

rapidxml::xml_attribute<> * kgl::ParsePublicationXMLImpl::validAttribute( rapidxml::xml_node<> * node_ptr,
                                                                          const char* attribute,
                                                                          const std::string& pmid) {

  rapidxml::xml_attribute<>* attrib_ptr = node_ptr->first_attribute(attribute);
  if (attrib_ptr == nullptr) {

    std::string error_message = "Pmid: " + pmid + " Attribute: '" + attribute + "' not found";
    throw std::runtime_error(error_message);

  }

  return attrib_ptr;

}


std::string kgl::ParsePublicationXMLImpl::validOptionalNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name) {

  rapidxml::xml_node<>* sub_node_ptr = node_ptr->first_node(sub_node_name);
  if (sub_node_ptr == nullptr) {

    return {""};

  }

  return sub_node_ptr->value();

}
