//
// Created by kellerberrin on 16/8/21.
//

#include "kgl_pubmed_xml_parser.h"
#include "kel_exec_env.h"


namespace kgl = kellerberrin::genome;



kgl::LitCitationMap kgl::ParseCitationXMLImpl::parseCitationXML(const std::string& citation_text) {

  LitCitationMap citation_map;

  // All this raw pointer stuff is very nasty, is rapidxml the right library?
  rapidxml::xml_document<> doc;
  std::unique_ptr<char[]> text_buffer_ptr(std::make_unique<char[]>(citation_text.size() + 1));
  std::memcpy(text_buffer_ptr.get(), citation_text.c_str(), citation_text.size() + 1);

  // Parse the buffer using the xml file parsing library into doc
  doc.parse<0>(text_buffer_ptr.get());

  rapidxml::xml_node<>* root_node = doc.first_node(CITATION_ROOT_NODE_);

  if (root_node != nullptr) {

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

              }

              link_node = link_node->next_sibling();

            }

          }

          auto [iter, result] = citation_map.try_emplace(pmid_node->value(), citation_pmid_set);
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



kgl::LitPublicationMap kgl::ParsePublicationXMLImpl::parsePublicationXML(const std::string& publication_xml_text) {

  LitPublicationMap publication_map;

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
    return publication_map; // Return the empty map.

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

      auto [iter, result] = publication_map.try_emplace(publication.pmid(), publication);
      if (not result) {

        std::string error_message = "cannot add duplicate publication: " + publication.pmid();
        throw std::runtime_error(error_message);

      }

    } catch(std::exception& e) {

      // Issue an error message and skip to the next article.
      ExecEnv::log().error("PubmedRequester::parsePublicationXML; error parsing XML Pubmed Article; {}", e.what());
      // Comment out after testing.
      // ExecEnv::log().error("PubmedRequester::parsePublicationXML; text:\n{}", publication_xml_text);


    }

    // Next article.
    article_node = article_node->next_sibling();

  }

  ExecEnv::log().info("PubmedRequester::parsePublicationXML; article count {}", publication_map.size());

  return publication_map;

}


kgl::PubMedPublicationSummary kgl::ParsePublicationXMLImpl::parsePubmedArticleXML(rapidxml::xml_node<> * pubmed_article_node) {


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
  kgl::PubMedPublicationSummary publication(pmid);

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
void kgl::ParsePublicationXMLImpl::parseAuthorsXML( rapidxml::xml_node<> * journal_article_node,
                                                    kgl::PubMedPublicationSummary& publication) {

  auto author_list_node = validSubNode(journal_article_node, AUTHOR_LIST_NODE_, publication.pmid());
  auto author_node = author_list_node->first_node(AUTHOR_NODE_);
  while(author_node != nullptr)  {

    auto surname = author_node->first_node(AUTHOR_SURNAME_NODE_);
    if (surname != nullptr) {

      auto initials = validSubNode(author_node, AUTHOR_INITIALS_NODE_, publication.pmid());
      std::string author_name = surname->value() + std::string("&") + initials->value();
      publication.addAuthor(author_name);

    } else {

      auto collective = validSubNode(author_node, AUTHOR_COLLECTIVE_NODE_, publication.pmid());
      std::string collective_name = collective->value() + std::string("&");
      publication.addAuthor(collective_name);

    }

    // Next author.
    author_node = author_node->next_sibling();

  }

}


// Unpack the various XML tags for the article.
void kgl::ParsePublicationXMLImpl::parseArticleFieldsXML( rapidxml::xml_node<> * journal_article_node,
                                                          kgl::PubMedPublicationSummary& publication) {

  auto abstract_node = journal_article_node->first_node(ABSTRACT_NODE_);
  if (abstract_node != nullptr) {

    auto abtract_text = validSubNode(abstract_node, ABSTRACT_TEXT_NODE_, publication.pmid());
    publication.abstract(abtract_text->value());

  } else {

    publication.abstract("");

  }

  auto article_title = validSubNode(journal_article_node, ARTICLE_TITLE_NODE_, publication.pmid());

  publication.title(article_title->value());

}

void kgl::ParsePublicationXMLImpl::parseDoiXML( rapidxml::xml_node<> * journal_article_node,
                                                rapidxml::xml_node<> * pubmed_node,
                                                PubMedPublicationSummary& publication) {


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
void kgl::ParsePublicationXMLImpl::parseJournalArticleXML( rapidxml::xml_node<> * journal_article_node,
                                                           kgl::PubMedPublicationSummary& publication) {

  auto journal_node = validSubNode(journal_article_node, JOURNAL_NODE_, publication.pmid());

  auto journal_title = validSubNode(journal_node, JOURNAL_TITLE_NODE_, publication.pmid());

  publication.journal(journal_title->value());

  auto journal_issue = validSubNode(journal_node, JOURNAL_ISSUE_NODE_, publication.pmid());

  publication.journalIssue(validOptionalNode(journal_issue, ISSUE_NODE_));

  publication.journalVolume(validOptionalNode(journal_issue, VOLUME_NODE_));


}

// Unpack the various XML tags for chemicals_.
void kgl::ParsePublicationXMLImpl::parseChemicalsXML( rapidxml::xml_node<> * medline_node,
                                                      kgl::PubMedPublicationSummary& publication) {

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
                                                 kgl::PubMedPublicationSummary& publication) {

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


void kgl::ParsePublicationXMLImpl::parseReferencesXML( rapidxml::xml_node<> * pubmed_node,
                                                        kgl::PubMedPublicationSummary& publication) {

  auto reference_list_node = pubmed_node->first_node(REFERENCE_LIST_NODE_);
  if (reference_list_node == nullptr) {

    return;

  }
  auto reference_node = reference_list_node->first_node(REFERENCE_NODE_);
  while(reference_node != nullptr)  {

    auto article_list_node = validSubNode(reference_node, ARTICLE_ID_LIST_, publication.pmid());
    auto id_node = validSubNode(article_list_node, ARTICLE_ID_, publication.pmid());
    auto id_attrib = validAttribute(id_node, ARTICLE_ID_ATTRIBUTE_, publication.pmid());
    std::string id_type = id_attrib->value();
    if (id_type != ARTICLE_ID_ATTRIBUTE_VALUE_) {

      ExecEnv::log().warn("ParsePublicationXMLImpl::parseReferencesXML; publication pmid_: {} has reference id type: {}",
                          publication.pmid(), id_type);

    }

    std::string reference_pmid = id_node->value();

    if (not publication.addReference(reference_pmid)) {

      ExecEnv::log().warn("ParsePublicationXMLImpl::parseReferencesXML; publication pmid_: {} has duplicate reference: {}",
                          publication.pmid(), reference_pmid);

    }

    // Next reference.
    reference_node = reference_node->next_sibling();

  }

}

void kgl::ParsePublicationXMLImpl::parseXMLDate(rapidxml::xml_node<> * journal_article_node,
                                                rapidxml::xml_node<> * pubmed_node,
                                                PubMedPublicationSummary& publication) {


  auto journal_node = validSubNode(journal_article_node, JOURNAL_NODE_, publication.pmid());
  auto journal_issue = validSubNode(journal_node, JOURNAL_ISSUE_NODE_, publication.pmid());
  auto pub_date = validSubNode(journal_issue, PUB_DATE_NODE_, publication.pmid());

  auto year = pub_date->first_node(YEAR_NODE_);
  auto month = pub_date->first_node(MONTH_NODE_);

  std::string date_string;
  if (year != nullptr and month != nullptr) {

    auto day =  validOptionalNode(pub_date, DAY_NODE_);

    if (not day.empty()) {

      date_string = day + "-" + std::string(month->value()) + "-" + std::string(year->value());

    } else {

      date_string = std::string(month->value()) + "-" + std::string(year->value());

    }

  } else {

    auto history_node = validSubNode(pubmed_node, HISTORY_NODE_, publication.pmid());
    auto pubmed_date_node = validSubNode(history_node, PUBMED_DATE_NODE_, publication.pmid());

    auto pm_year = pubmed_date_node->first_node(YEAR_NODE_);
    auto pm_month = pubmed_date_node->first_node(MONTH_NODE_);

    if (pm_year != nullptr and pm_month != nullptr) {

      auto day =  validOptionalNode(pubmed_date_node, DAY_NODE_);

      if (not day.empty()) {

        date_string = day + "-" + std::string(pm_month->value()) + "-" + std::string(pm_year->value());

      } else {

        date_string = std::string(pm_month->value()) + "-" + std::string(pm_year->value());

      }

    }

  }

  publication.publicationDate(date_string);

}



rapidxml::xml_node<> * kgl::ParsePublicationXMLImpl::validSubNode(rapidxml::xml_node<> * node_ptr, const char* sub_node_name, const std::string& pmid) {

  rapidxml::xml_node<>* sub_node_ptr = node_ptr->first_node(sub_node_name);
  if (sub_node_ptr == nullptr) {

    std::string error_message = "Pmid: " + pmid + " Node: '" + sub_node_name + "' not found";
    throw std::runtime_error(error_message);

  }

  return sub_node_ptr;

}

rapidxml::xml_attribute<> * kgl::ParsePublicationXMLImpl::validAttribute(rapidxml::xml_node<> * node_ptr, const char* attribute, const std::string& pmid) {

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
