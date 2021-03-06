//
// Created by kellerberrin on 23/5/20.
//

#include "kgl_variant_factory_vcf_evidence_memory.h"
#include "kgl_variant_factory_vcf_evidence.h"



namespace kgl = kellerberrin::genome;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Resource allocation objects.


size_t kgl::StringResourceInstance::allocateViews(size_t views) {

  size_t offset = string_allocation_.size();

  for (size_t index = 0; index < views; ++index) {

    string_allocation_.emplace_back(std::string_view());

  }

  return offset;

}


bool kgl::DynamicResourceInstance::queueResource(InfoArrayIndex requested) {

  auto result = dynamic_allocation_.try_emplace(requested.infoVariableIndex(), requested);

  return result.second;

}




const kgl::InfoResourceHandle kgl::InfoResourceAllocator::resourceRequest(DataResourceType resource_type, DataDynamicType dynamic_type, size_t data_size) {

  // Refuse allocation if a type mis-match.
  if (resource_type != resource_type_) {

    ExecEnv::log().error("InfoResourceAllocator::resourceRequest, mis-matched resource allocation requested, no data allocated");

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


const kgl::InfoResourceHandle kgl::InfoStringAllocator::resourceRequest(DataResourceType resource_type, DataDynamicType dynamic_type, size_t data_size) {

  // Refuse allocation if a type mis-match.
  if (resource_type_ != DataResourceType::String) {

    ExecEnv::log().error("InfoStringAllocator::resourceRequest, mis-matched resource allocation requested, no data allocated");

  }

  // Field identifier.
  size_t unique_id = uniqueHandleId();

  // Note that the character resource is allocated at runtime when the size of the string is known.
  // The resource object. We allocate string view resources (generally 0 or 1).
  InfoResourceHandle resource_handle(unique_id, view_resource_ptr_->allocateViews(data_size), data_size, dynamic_type, resource_type);

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
const kgl::InfoResourceHandle kgl::InfoMemoryResource::resourceRequest( DataResourceType resource_type,
                                                                        DataDynamicType dynamic_type,
                                                                        size_t data_size) {

  switch (resource_type) {

    case DataResourceType::Boolean:
      return bool_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

    case DataResourceType::Integer:
      return integer_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

    case DataResourceType::Float:
      return float_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

    default:
    case DataResourceType::String:
      return string_allocator_->resourceRequest(resource_type, dynamic_type, data_size);

  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the InfoMemoryResource object that controls memory resource allocation for all objects.

kgl::InfoMemoryResource::InfoMemoryResource() {

  unique_ident_ = std::make_shared<FixedResourceInstance>("Unique resource identifier");

  bool_memory_  = std::make_shared<FixedResourceInstance>("Allocated bool char size");
  char_memory_  = std::make_shared<FixedResourceInstance>("Allocated string char size");
  integer_memory_ = std::make_shared<FixedResourceInstance>("Allocated integer size");
  float_memory_ = std::make_shared<FixedResourceInstance>("Allocated float size");
  view_memory_ = std::make_shared<StringResourceInstance>("Allocated string_view size");
  array_memory_ = std::make_shared<DynamicResourceInstance>("Allocated dynamic_array size");

  bool_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Boolean, unique_ident_, bool_memory_, array_memory_);
  integer_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Integer, unique_ident_, integer_memory_, array_memory_);
  float_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Float, unique_ident_, float_memory_, array_memory_);
  string_allocator_= std::make_unique<InfoStringAllocator>(unique_ident_, view_memory_, char_memory_, array_memory_);

}

kgl::InfoMemoryResource::InfoMemoryResource(const InfoMemoryResource& copy) {

  unique_ident_ = std::make_shared<FixedResourceInstance>(*copy.unique_ident_);

  bool_memory_  = std::make_shared<FixedResourceInstance>(*copy.bool_memory_);
  char_memory_  = std::make_shared<FixedResourceInstance>(*copy.char_memory_);
  integer_memory_ = std::make_shared<FixedResourceInstance>(*copy.integer_memory_);
  float_memory_ = std::make_shared<FixedResourceInstance>(*copy.float_memory_);
  view_memory_ = std::make_shared<StringResourceInstance>(*copy.view_memory_);
  array_memory_ = std::make_shared<DynamicResourceInstance>(*copy.array_memory_);

  bool_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Boolean, unique_ident_, bool_memory_, array_memory_);
  integer_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Integer, unique_ident_, integer_memory_, array_memory_);
  float_allocator_= std::make_unique<InfoResourceAllocator>(DataResourceType::Float, unique_ident_, float_memory_, array_memory_);
  string_allocator_= std::make_unique<InfoStringAllocator>(unique_ident_, view_memory_, char_memory_, array_memory_);

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

  } else {  // Dynamic or FixedDynamic.

    // Only allocate dynamic data if pre-allocated space not available (this is a major space saving).
    if (item_resource_handle.initialDataSize() != token.second) {
      // Dynamic data.
      if (not resolveDynamic(item_resource_handle, token)) {

        return false;

      }

    }

  }

  return true;

}


bool kgl::InfoMemoryResource::resolveDynamic(const InfoResourceHandle& item_resource_handle, const InfoParserToken& token) {

  // Check if the space has been speculatively pre-allocated, this is done with 'A' alternate allele variables for efficiency.
  // As these variables can be arrays but are most often just a scalar.
  if (not (item_resource_handle.initialDataSize() != 0 and item_resource_handle.initialDataSize() == token.second)) {

    // If the resource is FixedDynamic and the space pre-allocated for it mis-matches the data runtime size
    // then we treat it as a dynamic datatype and request a dynamic resource.
    if (item_resource_handle.dynamicType() == DataDynamicType::FixedDynamic) {

      requestDynamic(InfoArrayIndex(item_resource_handle.handleId(), 0, token.second));

    }

    // Finds the dynamic resource and allocates data space for it.
    if (not findUpdateDynamic(item_resource_handle, token)) {

      return false;

    }

  } // if not pre-allocated.

  // All strings must allocate character space at run time.
  if (item_resource_handle.resourceType() == DataResourceType::String) {

    string_allocator_->allocateStringChars(item_resource_handle.initialDataSize(), token);

  }

  return true;

}


bool kgl::InfoMemoryResource::findUpdateDynamic(const InfoResourceHandle& item_resource_handle, const InfoParserToken& token) {


  // Find the dynamic block.

  auto array_item_ptr = array_memory_->dynamicAllocation().find(item_resource_handle.handleId());

  if (array_item_ptr == array_memory_->dynamicAllocation().end()) {
    // If the dynamic item was not found.
    ExecEnv::log().error( "InfoMemoryResource::resolveDynamic, Unable to locate Dynamic resource identifier: {}, token size: {}, token value: {}",
                          item_resource_handle.handleId(), token.second, std::string(token.first));
    return false;

  } // if not found


  switch (item_resource_handle.resourceType()) {

    case DataResourceType::Boolean:
      // Allow for the possibility of boolean vectors.
      array_item_ptr->second.infoOffset(bool_allocator_->allocateResource(token.second));
      array_item_ptr->second.infoSize(token.second);
      break;

    case DataResourceType::Integer:
      array_item_ptr->second.infoOffset(integer_allocator_->allocateResource(token.second));
      array_item_ptr->second.infoSize(token.second);
      break;

    case DataResourceType::Float:
      array_item_ptr->second.infoOffset(float_allocator_->allocateResource(token.second));
      array_item_ptr->second.infoSize(token.second);
      break;

    case DataResourceType::String:
      array_item_ptr->second.infoOffset(string_allocator_->allocateViews(token.second));
      array_item_ptr->second.infoSize(token.second);
      break;

  } // switch


  return true;

}

