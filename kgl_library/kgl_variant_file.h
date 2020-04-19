//
// Created by kellerberrin on 15/4/20.
//

#ifndef KGL_VARIANT_FILE_H
#define KGL_VARIANT_FILE_H

#include "kgl_genome_types.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <string>
#include <vector>
#include <fstream>

namespace bio = boost::iostreams;
namespace kellerberrin::genome {   //  organization::project level namespace

// Basic VCF Record LIne.
////////////////////////////////////////////////////////////////////////////////////////////////////////

class VcfRecord
{
public:

  // Default constructor.
  VcfRecord() : offset(INVALID_POS), qual(MISSING_QUAL) {}

  static std::unique_ptr<VcfRecord> EOF_RECORD() { return std::unique_ptr<VcfRecord>(); }

  // Numeric id of the reference sequence.
  ContigId_t contig_id;
  // Position on the reference.
  ContigOffset_t offset;
  // Textual identifier of the variant.
  std::string id;
  // Bases in the reference.
  std::string ref;
  // Bases in the alternatives, comma-separated.
  std::string alt;
  // Quality
  double qual;
  // Value of FILTER field.
  std::string filter;
  // Value of INFO field.
  std::string info;
  // Value of FORMAT field.
  std::string format;
  // The genotype infos.
  std::vector<std::string> genotypeInfos;


private:
  // Constant for invalid position.
  static constexpr const ContigOffset_t INVALID_POS = std::numeric_limits<ContigOffset_t>::max();
  // Undefined quality number.
  static constexpr const double MISSING_QUAL = std::numeric_limits<double>::lowest();

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plug in raw IO objects to read text and gzipped files.

class TextStreamIO {

public:

  TextStreamIO() = default;
  TextStreamIO(const TextStreamIO &) = delete;
  ~TextStreamIO() = default;

  bool open(const std::string &file_name);
  bool readLine(std::string &text_line) { return not std::getline(file_, text_line).eof(); }

private:

  std::ifstream file_;

};


inline bool TextStreamIO::open(const std::string &file_name) {

  try {

    // Open input file.

    file_.open(file_name);
    if (not file_.good()) {

      ExecEnv::log().error("TextStreamIO; I/O error; could not open file: {}", file_name);
      return false;

    }
  }
  catch (std::exception const &e) {

    ExecEnv::log().error("TextStreamIO; Opening file: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GZStreamIO {

public:

  GZStreamIO() = default;
  GZStreamIO(const GZStreamIO &) = delete;
  ~GZStreamIO() = default;

  bool open(const std::string &file_name);
  bool readLine(std::string &text_line) { return not std::getline(gz_file_, text_line).eof(); }

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream gz_file_;

};

inline bool GZStreamIO::open(const std::string &file_name) {

  try {

    // Open input file.

    file_.open(file_name, std::ios_base::in | std::ios_base::binary);

    if (not file_.good()) {

      ExecEnv::log().error("GZStreamIO; I/O error; could not open file: {}", file_name);
      return false;

    }

    gz_file_.push(bio::gzip_decompressor());
    gz_file_.push(file_);

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("GZStreamIO; Opening file: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parses the VCF header for Contig and Genome infromation.

class VCFParseHeader {

public:

  VCFParseHeader() = default;

  ~VCFParseHeader() = default;

  template<class IOStream> bool parseHeader(const std::string& file_name);

  [[nodiscard]] const std::vector<std::string>& getGenomes() const { return vcf_genomes_; }

private:

  std::vector<std::string> vcf_genomes_;                // Field (genome) names for each VCF record

  // Parser constants.
  static constexpr const char* CONTIG_NAME_FRAGMENT_{"##contig"};
  static constexpr const size_t CONTIG_NAME_FRAGMENT_LENGTH_{8};
  static constexpr const char* CONTIG_INFO_START_{"=<"};
  static constexpr const char* CONTIG_INFO_END_{">"};
  static constexpr const char* FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.


};

template<class IOStream>
bool VCFParseHeader::parseHeader(const std::string& vcf_file_name) {

  IOStream vcf_file;

  // Open input file.
  if (not vcf_file.open(vcf_file_name)) {

    ExecEnv::log().critical("I/O error; could not open VCF file: {}", vcf_file_name);

  }

  try {

    long counter = 0;
    bool found_header = false;

    while (true) {

      std::string record_str;
      if (not vcf_file.readLine(record_str)) break;

      std::string contig_prefix = record_str.substr(0, CONTIG_NAME_FRAGMENT_LENGTH_);
      if (contig_prefix == CONTIG_NAME_FRAGMENT_) {

        std::string contig_string = Utility::trimAllWhiteSpace(record_str);

      }

      std::string line_prefix = record_str.substr(0, FIELD_NAME_FRAGMENT_LENGTH_);
      if (line_prefix == FIELD_NAME_FRAGMENT_) {

        found_header = true;
        std::vector<std::string> field_vector;
        if (not ParseVCFMiscImpl::tokenize(record_str, "\t",  field_vector)) {

          ExecEnv::log().error("Unable to parse VCF header line: {}", record_str);

        }
        size_t field_count = 0;
        for(auto const& field :field_vector) {

          if (field_count >= SKIP_FIELD_NAMES_) {

            vcf_genomes_.push_back(field);

          }

          ++field_count;

        }

        break; // #CHROM is the last field in the VCF header so stop processing.

      }

      ++counter;

    }

    if (not found_header) {

      ExecEnv::log().error("VCF Genome Names Not Found");

    } else {

      ExecEnv::log().info("{} Genomes in VCF Header {}, Header lines processed: {}", vcf_genomes_.size(), vcf_file_name, counter);

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCFParseHeader::parseHeader; VCF file: {}, unexpected I/O exception: {}", vcf_file_name, e.what());

  }

  return true;

}

} // end namespace


#endif //KGL_VARIANT_FILE_H
