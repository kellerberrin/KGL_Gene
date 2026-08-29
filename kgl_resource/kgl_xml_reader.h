//
// Created by kellerberrin on 26/8/26.
//

#ifndef KGL_XML_READER_H
#define KGL_XML_READER_H

#include "kel_exec_env.h"
#include "kel_property_tree.h"

#include <string>
#include <string_view>
#include <optional>
#include <vector>


namespace kellerberrin::genome {   //  organization::project level namespace


/// A fluent, chainable reader over a PropertyTree sub-tree.
/// Each read logs a contextual error on failure and records the failure so the
/// caller can short-circuit with a single `if (not reader.ok()) return std::nullopt;`.
class XmlReader {

public:

  explicit XmlReader(const PropertyTree& tree, std::string_view context) : tree_(tree), context_(context) {}

  /// Read a required string property. On failure logs and marks the reader failed.
  XmlReader& requiredString(std::string_view tag, std::string& out) {

    if (ok_ and not tree_.getProperty(std::string(tag), out)) {

      logError(tag, "string property");

    }
    return *this;

  }

  /// Read an optional string property; returns the value wrapped in std::optional.
  [[nodiscard]] std::optional<std::string> optionalString(std::string_view tag) {

    if (not ok_) return std::nullopt;
    std::string value;
    if (tree_.getProperty(std::string(tag), value)) {

      return value;

    }
    return std::nullopt;

  }

  /// Read an optional string property into an out-parameter; returns true if present.
  [[nodiscard]] bool optionalString(std::string_view tag, std::string& out) {

    if (not ok_) return false;
    return tree_.getProperty(std::string(tag), out);

  }

  /// Read a required file path (resolved against the work directory).
  XmlReader& requiredFile(std::string_view tag, std::string_view work_dir, std::string& out) {

    if (ok_ and not tree_.getFileProperty(std::string(tag), std::string(work_dir), out)) {

      logError(tag, "file");

    }
    return *this;

  }

  /// Read an optional file path; returns true if present.
  [[nodiscard]] bool optionalFile(std::string_view tag, std::string_view work_dir, std::string& out) {

    if (not ok_) return false;
    return tree_.getOptionalFileProperty(std::string(tag), std::string(work_dir), out);

  }

  /// Read a file path, creating it if absent.
  XmlReader& createFile(std::string_view tag, std::string_view work_dir, std::string& out) {

    if (ok_ and not tree_.getFileCreateProperty(std::string(tag), std::string(work_dir), out)) {

      logError(tag, "creatable file");

    }
    return *this;

  }

  /// Read a required child node vector.
  XmlReader& requiredNodes(std::string_view tag, std::vector<SubPropertyTree>& out) {

    if (ok_ and not tree_.getPropertyTreeVector(std::string(tag), out)) {

      logError(tag, "node list");

    }
    return *this;

  }

  /// True if all required reads so far succeeded.
  [[nodiscard]] bool ok() const { return ok_; }

private:

  void logError(std::string_view tag, std::string_view kind) {

    ExecEnv::log().error("{}; missing {} '{}'", context_, kind, tag);
    ok_ = false;

  }

  const PropertyTree& tree_;
  std::string_view context_;
  bool ok_{true};

};


} // namespace


#endif //KGL_XML_READER_H
