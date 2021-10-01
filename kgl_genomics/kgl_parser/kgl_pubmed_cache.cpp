//
// Created by kellerberrin on 20/8/21.
//

#include "kgl_pubmed_cache.h"
#include "kgl_pubmed_xml_parser.h"  // Only used in implementation files.
#include <fstream>


namespace kgl = kellerberrin::genome;



kgl::LitPublicationMap kgl::PubmedAPICache::readCachedPublications() const {

  LitPublicationMap cached_pub_map;

  auto uncited_pub_map = readPublicationCache();

  // Merge the citations.
  auto citation_map = readCitationCache();
  for (auto& [pmid, publication_ptr] : uncited_pub_map) {

    auto result = citation_map.find(pmid);
    if (result != citation_map.end()) {

      auto const& [cite_pmid, cite_date_set] = *result;
      auto const& [cite_date, cite_set] = cite_date_set;

      // Check the download dates.
      if (cite_date != publication_ptr->downloadDate()) {

        ExecEnv::log().warn("PubmedAPICache::readCachedPublications; citation download date: {} does not match publication download date: {}",
                            cite_date.text(), publication_ptr->downloadDate().text());

      }

      publication_ptr->citations(cite_set);
      cached_pub_map.emplace(pmid, publication_ptr);

    }

  }

  ExecEnv::log().info("Cached Pubmed Cited publications: {}, Publications: {}, Citations: {}",
                      cached_pub_map.size(), uncited_pub_map.size(), citation_map.size());

  return cached_pub_map;

}


kgl::APIPublicationMap kgl::PubmedAPICache::readPublicationCache() const {

  APIPublicationMap publication_map;

  std::ifstream publication_xml_file(publication_cache_file_);

  if (not publication_xml_file.good()) {

    ExecEnv::log().warn("PubmedAPICache::readPublicationCache; unable to open cache file: {}", publication_cache_file_);
    return publication_map;

  }

  std::string cached_records;
  DateGP download_date;
  while(readCacheRecord(publication_xml_file, cached_records, download_date)) {

    auto [parse_result, detail_map] = ParsePublicationXMLImpl::parsePublicationXML(cached_records, download_date);
    if (parse_result) {

      publication_map.merge(detail_map);

    }

  }

  ExecEnv::log().info( "Found XML publication details: {} in cache file: {}",
                       publication_map.size(), publication_cache_file_);

  return publication_map;

}


kgl::LitCitationMap kgl::PubmedAPICache::readCitationCache() const {

  LitCitationMap citation_map;

  std::ifstream citation_xml_file(citation_cache_file_);

  if (not citation_xml_file.good()) {

    ExecEnv::log().warn("PubmedAPICache::readCitationCache; unable to open cache file: {}", citation_cache_file_);
    return citation_map;

  }

  std::string cached_records;
  DateGP download_date;
  while(readCacheRecord(citation_xml_file, cached_records, download_date)) {

    auto [parse_result, cite_map] = ParseCitationXMLImpl::parseCitationXML(cached_records, download_date);
    if (parse_result) {

      citation_map.merge(cite_map);

    }

  }

  ExecEnv::log().info( "Found XML citation details: {} in cache file: {}",
                       citation_map.size(), citation_cache_file_);

  return citation_map;

}


// Physical file operations.
bool kgl::PubmedAPICache::writeCitationCache(const std::string& xml_cache_record) const {

  std::ofstream citation_xml_file(citation_cache_file_, std::ios::app);

  if (not citation_xml_file.good()) {

    ExecEnv::log().error("PubmedAPICache::writeCitationCache; unable open cache file: {}", citation_cache_file_);
    return false;

  }
  DateGP todays_date;
  todays_date.setToday();

  // Wrap the XML text in pseudo begin and end tags
  citation_xml_file << std::string(START_CACHE_NODE_);
  citation_xml_file << "\"" << std::to_string(xml_cache_record.size()) << "\"";
  citation_xml_file << ATTRIBUTE_SEPARATOR_;
  citation_xml_file << std::string(START_CACHE_DOWNLOAD_DATE_);
  citation_xml_file << "\"" << todays_date.text() << "\"";
  citation_xml_file << START_CACHE_NODE_END_;
  citation_xml_file << xml_cache_record;
  citation_xml_file << std::string(END_CACHE_NODE_);

  return true;

}


bool kgl::PubmedAPICache::writePublicationCache(const std::string& xml_cache_record) const {

  std::ofstream detail_xml_file(publication_cache_file_, std::ios::app);

  if (not detail_xml_file.good()) {

    ExecEnv::log().error("PubmedAPICache::writePublicationCache; unable open cache file: {}", publication_cache_file_);
    return false;

  }
  DateGP todays_date;
  todays_date.setToday();

  // Wrap the XML text in pseudo begin and end tags
  detail_xml_file << std::string(START_CACHE_NODE_);
  detail_xml_file << "\"" << std::to_string(xml_cache_record.size()) << "\"";
  detail_xml_file << ATTRIBUTE_SEPARATOR_;
  detail_xml_file << std::string(START_CACHE_DOWNLOAD_DATE_);
  detail_xml_file << "\"" << todays_date.text() << "\"";
  detail_xml_file << START_CACHE_NODE_END_;
  detail_xml_file << xml_cache_record;
  detail_xml_file << std::string(END_CACHE_NODE_);

  return true;

}



bool kgl::PubmedAPICache::readCacheRecord(std::istream& input, std::string& record_string, DateGP& download_date) const {

  record_string.clear();

  if (input.eof()) {

    return false;

  }

  char node_buffer[128];
  // Read Cache Node:
  input.read(node_buffer, START_CACHE_NODE_.size());
  if (input.eof()) {

    return false;

  }

  node_buffer[START_CACHE_NODE_.size()] = '\0';

  if (std::string(node_buffer) != START_CACHE_NODE_) {

    ExecEnv::log().error("PubmedAPICache::readCacheRecord; unexpected start node: {}, XML cache is corrupt", std::string(node_buffer));
    return false;

  }

  std::string xml_size;
  bool found_comma{false};
  while(auto ch = input.get()) {

    if (ch == START_CACHE_NODE_END_) {

      break;

    }

    if (ch == ATTRIBUTE_SEPARATOR_) {

      found_comma = true;
      break;

    }

    if (input.eof()) {

      ExecEnv::log().error("PubmedAPICache::readCacheRecord; unexpected EOF encountered, XML cache is corrupt");
      return false;

    }

    if (ch != '\"') {

      xml_size += static_cast<char>(ch);

    }

  }

  if (found_comma) {

    input.read(node_buffer, START_CACHE_DOWNLOAD_DATE_.size());
    if (input.eof()) {

      ExecEnv::log().error("PubmedAPICache::readCacheRecord; unexpected EOF encountered, XML cache is corrupt");
      return false;

    }

    node_buffer[START_CACHE_DOWNLOAD_DATE_.size()] = '\0';

    if (std::string(node_buffer) != START_CACHE_DOWNLOAD_DATE_) {

      ExecEnv::log().error("PubmedAPICache::readCacheRecord; unexpected start node attribute: {}, XML cache is corrupt", std::string(node_buffer));
      return false;

    }

    // Read the Date format "YYYY-MMM-DD"
    char date_buffer[32];
    input.read(date_buffer, DateGP::TEXTSIZE+2);

    if (input.eof()) {

      ExecEnv::log().error("PubmedAPICache::readCacheRecord; unexpected EOF encountered, XML cache is corrupt");
      return false;

    }

    // Remove the quote characters.
    date_buffer[DateGP::TEXTSIZE+1] = '\0';
    DateGP attribute_date(&date_buffer[1]);
    if (attribute_date.notInitialized()) {

      ExecEnv::log().error("PubmedAPICache::readCacheRecord; invalid download date");
      return false;

    }

    download_date = attribute_date;

    // Get the node closing bracket '>'
    auto ch = input.get();
    if (ch != START_CACHE_NODE_END_) {

      ExecEnv::log().error("PubmedAPICache::readCacheRecord; expected '>' , XML cache is corrupt");
      return false;

    }

  }

  size_t record_size{0};
  try {

    record_size = std::stoll(xml_size);

  } catch(const std::exception& e) {

    ExecEnv::log().error("PubmedAPICache::readCacheRecord; invalid XML record size text: {}, exception: {}", xml_size, e.what());
    return false;

  }

  if (record_size <= MIN_CACHE_SIZE_ or record_size >= MAX_CACHE_SIZE_) {

    ExecEnv::log().error("PubmedAPICache::readCacheRecord; invalid XML cache size: {}, XML cache is corrupt", record_size);
    return false;

  }

  // Read the xml record into the string provided.
  record_string.resize(record_size);
  input.read(record_string.data(), record_string.size());

  if (input.eof()) {

    ExecEnv::log().error("PubmedAPICache::readCacheRecord; unexpected EOF encountered, XML cache is corrupt");
    return false;

  }

  input.read(node_buffer, END_CACHE_NODE_.size());
  node_buffer[END_CACHE_NODE_.size()] = '\0';

  if (std::string(node_buffer) != END_CACHE_NODE_) {

    ExecEnv::log().error("PubmedAPICache::readCacheRecord; unexpected end node: {}, XML cache is corrupt", std::string(node_buffer));
    return false;

  }

  return true;

}

