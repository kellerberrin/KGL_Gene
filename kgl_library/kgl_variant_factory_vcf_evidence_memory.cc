//
// Created by kellerberrin on 23/5/20.
//

#include "kgl_variant_factory_vcf_evidence_memory.h"
#include "kgl_variant_factory_vcf_evidence.h"



namespace kgl = kellerberrin::genome;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// data_item_count is the number of data items, for a float array of size 3 this is 1 for 1 array. For a string array of size 3 it is 3.
// data_item_size is the number of underlying data items. For float array of size 3 this is 3.
// For a string array this will be the total number of characters in all strings. For a string array of size 3 data_item_count=3, but data_item_size is size of
// all three strings added together. The total characters required for all 3 strings.
size_t kgl::InfoDataUsageCount::staticIncrementAndAllocate(InfoEvidenceIntern internal_type) {

  size_t item_offset{0};
  switch (internal_type) {

    case InfoEvidenceIntern::intern_char:  {

      item_offset = char_count_;
      ++char_count_;

    }// Size is fixed (boolean variables) and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_unity_integer_array:
    case InfoEvidenceIntern::intern_integer: {

      item_offset = integer_count_;
      ++integer_count_;

    }// Size is fixed and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_unity_float_array:
    case InfoEvidenceIntern::intern_float: {

      item_offset = float_count_;
      ++float_count_;  // should be 1

    }// Size is fixed and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_string: {

      item_offset = string_count_;
      ++string_count_;   // allocate a std::string_view

    } // Size varies between records.
      return item_offset;

    case InfoEvidenceIntern::intern_integer_array: {

      item_offset = array_count_;
      ++array_count_;

    } // Size is fixed and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_float_array: {

      item_offset = array_count_;
      ++array_count_;

    }
      return item_offset;


    case InfoEvidenceIntern::intern_string_array: {

      item_offset = array_count_;
      ++array_count_;

    } // Size varies between records.
      return item_offset;

    case InfoEvidenceIntern::intern_unity_string_array: {

      item_offset = array_count_;
      ++array_count_;

    }   // Size varies between records.
      return item_offset;

    case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
    default:
      return item_offset;

  }

}


bool kgl::InfoDataUsageCount::dynamicIncrementAndAllocate(InfoEvidenceIntern internal_type, const InfoParserToken& token) {

//   = subscribed_field.evidenceType().InternalInfoType();

  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      if (token.second != 0) {

        ExecEnv::log().warn("InfoDataUsageCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}", std::string(token.first), token.second);
        return false;
      }
      break;

    case InfoEvidenceIntern::intern_integer:
    case InfoEvidenceIntern::intern_float:
      if (token.second != 1) {

        ExecEnv::log().warn("InfoDataUsageCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}", std::string(token.first), token.second);
        return false;
      }
      break;

    case InfoEvidenceIntern::intern_unity_integer_array:
      if (token.second > 1) {

        ++unity_array_count_;
        integer_count_ += token.second;

      }
      break;

    case InfoEvidenceIntern::intern_unity_float_array:
      if (token.second > 1) {

        ++unity_array_count_;
        float_count_ += token.second;

      }
      break;

    case InfoEvidenceIntern::intern_string: {

      char_count_ += token.first.size(); // total char size of all strings.

    }
      break;

    case InfoEvidenceIntern::intern_integer_array: {

      integer_count_ += token.second; // size of the array

    } // Size is fixed and known at Info data subscription time.
      break;

    case InfoEvidenceIntern::intern_float_array: {

      float_count_ += token.second; // size of the array

    } // Size is fixed and known at Info data subscription time.
      break;

    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array: {

      char_count_ += token.first.size() - (token.second - 1); // total char size of all strings, less the delimiter chars.
      string_count_ += token.second;   // number strings, allocate a vector of std::string_views

    } // Size varies between records.
      break;

    case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
    default:
      break;

  }

  return true;

}


bool kgl::InfoDataUsageCount::operator==(const InfoDataUsageCount& cmp) const {

  bool result = unityArrayCount() == cmp.unityArrayCount() and
                arrayCount() == cmp.arrayCount() and
                floatCount() == cmp.floatCount() and
                integerCount() == cmp.integerCount() and
                stringCount() == cmp.stringCount() and
                charCount() == charCount();

  return result;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Resource allocation object.


std::optional<const kgl::InfoResourceHandle> kgl::InfoResourceAllocator::resourceRequest(DataResourceType resource_type, DataDynamicType dynamic_type, size_t data_size) {

  // Refuse allocation if a type mis-match.
  if (resource_type != resource_type_) {

    ExecEnv::log().error("InfoResourceAllocator::resourceRequest, mis-matched resource allocation requested, no data allocated");
    return std::nullopt;

  }

  // Field identifier.
  size_t unique_id = uniqueHandleId();

  // The resource object.
  InfoResourceHandle resource_handle(unique_id, fixed_resource_ptr_->allocateResource(data_size), data_size, dynamic_type, resource_type);

  // A dynamic array requested, queue until runtime size information available.
  if (dynamic_type == DataDynamicType::DynamicData) {

    // No memory offset calculated yet.
    dynamic_resource_ptr_->queueResource(InfoArrayIndex(unique_id, 0, data_size));

  }

  return resource_handle;

}


std::optional<const kgl::InfoResourceHandle> kgl::InfoStringAllocator::resourceRequest(DataResourceType resource_type, DataDynamicType dynamic_type, size_t data_size) {

  // Refuse allocation if a type mis-match.
  if (resource_type_ != DataResourceType::String) {

    ExecEnv::log().error("InfoStringAllocator::resourceRequest, mis-matched resource allocation requested, no data allocated");
    return std::nullopt;

  }

  // Field identifier.
  size_t unique_id = uniqueHandleId();

  // Note that the character resource is allocated at runtime when the size of the string is known.
  // The resource object. We allocate string view resources (generally 0 or 1).
  InfoResourceHandle resource_handle(unique_id, view_resource_ptr_->allocateResource(data_size), data_size, dynamic_type, resource_type);

  // A dynamic array requested, queue until runtime size information available.
  if (dynamic_type == DataDynamicType::DynamicData) {

    // No memory offset calculated yet.
    dynamic_resource_ptr_->queueResource(InfoArrayIndex(unique_id, 0, data_size));

  }

  return resource_handle;

}


void kgl::InfoStringAllocator::allocateStringChars(size_t size, const InfoParserToken& token)
{

  size_t size_chars;
  if (size > 1) {

    size_chars = token.first.size() - (size - 1); // Adjust for delimiter chars ','

  } else {

    size_chars = token.first.size();

  }

  char_resource_ptr_->allocateResource(size_chars);

}


// The initial resource request.
std::optional<const kgl::InfoResourceHandle> kgl::InfoMemoryResource::resourceRequest( DataResourceType resource_type,
                                                                                       DataDynamicType dynamic_type,
                                                                                       size_t data_size) {

  switch (resource_type) {

    case DataResourceType::Boolean:
      return bool_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

    case DataResourceType::Integer:
      return integer_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

    case DataResourceType::Float:
      return float_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

    case DataResourceType::String:
      return string_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

  }

  return std::nullopt;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the InfoMemoryResource object that controls memory resource allocation for all objects.

kgl::InfoMemoryResource::InfoMemoryResource() {

  unique_ident_ = std::make_shared<FixedResourceInstance>("Unique resource identifier", 0);

  char_memory_  = std::make_shared<FixedResourceInstance>("Allocated char size", 0);
  integer_memory_ = std::make_shared<FixedResourceInstance>("Allocated integer size", 0);
  float_memory_ = std::make_shared<FixedResourceInstance>("Allocated float size", 0);
  view_memory_ = std::make_shared<FixedResourceInstance>("Allocated string_view size", 0);
  array_memory_ = std::make_shared<DynamicResourceInstance>("Allocated dynamic_array size", std::vector<InfoArrayIndex>());

  bool_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Boolean, unique_ident_, char_memory_, array_memory_);
  integer_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Integer, unique_ident_, integer_memory_, array_memory_);
  float_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Float, unique_ident_, float_memory_, array_memory_);
  string_allocator_= std::make_unique<InfoStringAllocator>(unique_ident_, view_memory_, char_memory_, array_memory_);

}

kgl::InfoMemoryResource::InfoMemoryResource(const InfoMemoryResource& copy) {

  unique_ident_ = std::make_shared<FixedResourceInstance>(*copy.unique_ident_);

  char_memory_  = std::make_shared<FixedResourceInstance>(*copy.char_memory_);
  integer_memory_ = std::make_shared<FixedResourceInstance>(*copy.integer_memory_);
  float_memory_ = std::make_shared<FixedResourceInstance>(*copy.float_memory_);
  view_memory_ = std::make_shared<FixedResourceInstance>(*copy.view_memory_);
  array_memory_ = std::make_shared<DynamicResourceInstance>(*copy.array_memory_);

  bool_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Boolean, unique_ident_, char_memory_, array_memory_);
  integer_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Integer, unique_ident_, integer_memory_, array_memory_);
  float_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Float, unique_ident_, float_memory_, array_memory_);
  string_allocator_= std::make_unique<InfoStringAllocator>(unique_ident_, view_memory_, char_memory_, array_memory_);

}


kgl::InfoMemoryResource& kgl::InfoMemoryResource::operator=(const InfoMemoryResource& copy) {

  unique_ident_ = std::make_shared<FixedResourceInstance>(*copy.unique_ident_);

  char_memory_  = std::make_shared<FixedResourceInstance>(*copy.char_memory_);
  integer_memory_ = std::make_shared<FixedResourceInstance>(*copy.integer_memory_);
  float_memory_ = std::make_shared<FixedResourceInstance>(*copy.float_memory_);
  view_memory_ = std::make_shared<FixedResourceInstance>(*copy.view_memory_);
  array_memory_ = std::make_shared<DynamicResourceInstance>(*copy.array_memory_);

  return *this;

}


bool kgl::InfoMemoryResource::resolveAllocation(const InfoResourceHandle& item_resource_handle, const InfoParserToken& token) {

  if (item_resource_handle.dynamicType() == DataDynamicType::FixedData) {

    if (item_resource_handle.initialDataSize() != token.second) {

      // A boolean with a token size of 0 is OK.
      if (not (item_resource_handle.resourceType() == DataResourceType::Boolean and token.second == 0)) {

        // A fixed string with embedded ',' characters causing the parser to mis-count is OK.
        if (not (item_resource_handle.resourceType() == DataResourceType::String
            and std::string(token.first).find(VCFInfoParser::INFO_VECTOR_DELIMITER_) != std::string::npos)) {

          ExecEnv::log().warn("InfoMemoryResource::resolveAllocation; Mismatch between allocated size: {} and data size: {}, data value: {}"
          , item_resource_handle.initialDataSize(), token.second, std::string(token.first));
          return false;

        }

      }

    }

    // All strings must allocate character space at run time.
    if (item_resource_handle.resourceType() == DataResourceType::String) {

      string_allocator_->allocateStringChars(item_resource_handle.initialDataSize(), token);

    }

  } else {


    if (not resolveDynamic(item_resource_handle, token)) {

      return false;

    }

    // Find dynamic resource and allocate at runtime. No need to allocate if data storage pre-allocated.

  }

  return true;

}


bool kgl::InfoMemoryResource::resolveDynamic(const InfoResourceHandle& item_resource_handle, const InfoParserToken& token) {

  // Check if the space has been speculatively pre-allocated, this is done with 'A' alternate allele variables for efficiency.
  // As these variables can be arrays but are most often just a scalar.
  if (not (item_resource_handle.initialDataSize() != 0 and item_resource_handle.initialDataSize() == token.second)) {

    // Find the dynamic block. Just a simple linear search for now.
    bool found = false;
    for (auto& array_item : array_memory_->dynamicAllocation()) {

      // Found the resource.
      if (array_item.infoVariableIndex() == item_resource_handle.handleId()) {

        found = true;
        switch (item_resource_handle.resourceType()) {

          case DataResourceType::Boolean:
          // Boolean arrays are not defined. Complain and exit
            ExecEnv::log().error("InfoMemoryResource::resolveDynamic, Unexpected Boolean Array, resource ident: {}, token size: {}, token value: {}",
                                 item_resource_handle.handleId(), token.second, std::string(token.first));
            return false;

          case DataResourceType::Integer:
            array_item.infoOffset(integer_allocator_->allocateResource(token.second));
            array_item.infoSize(token.second);
            break;

          case DataResourceType::Float:
            array_item.infoOffset(float_allocator_->allocateResource(token.second));
            array_item.infoSize(token.second);
            break;

          case DataResourceType::String:
            array_item.infoOffset(string_allocator_->allocateViews(token.second));
            array_item.infoSize(token.second);
            break;

        } // switch

      } // if found

    } // for

    if (not found) {
      // If we fall through the loop then the dynamic item was not found.
      ExecEnv::log().error(
      "InfoMemoryResource::resolveDynamic, Unable to locate Dynamic resource ident: {}, token size: {}, token value: {}",
      item_resource_handle.handleId(), token.second, std::string(token.first));
      return false;

    } // if not found

  } // if not pre-allocated.

  // All strings must allocate character space at run time.
  if (item_resource_handle.resourceType() == DataResourceType::String) {

    string_allocator_->allocateStringChars(item_resource_handle.initialDataSize(), token);

  }

  return true;

}

