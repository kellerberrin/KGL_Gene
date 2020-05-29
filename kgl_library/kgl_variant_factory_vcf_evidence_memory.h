//
// Created by kellerberrin on 23/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEMORY_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEMORY_H

#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_evidence_data.h"



namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object manages the storage of individual Info fields from the raw data allocated below.

struct InfoDataUsageCount {

public:

  InfoDataUsageCount() = default;
  ~InfoDataUsageCount() = default;

  bool operator==(const InfoDataUsageCount& cmp) const;

  [[nodiscard]] size_t unityArrayCount() const { return unity_array_count_; }
  [[nodiscard]] size_t arrayCount() const { return array_count_; }
  [[nodiscard]] size_t floatCount() const { return float_count_; }
  [[nodiscard]] size_t integerCount() const { return integer_count_; }
  [[nodiscard]] size_t stringCount() const { return string_count_; }
  [[nodiscard]] size_t charCount() const { return char_count_; }

  void unityArrayCountAdd(size_t count) { unity_array_count_ += count; }
  void arrayCountAdd(size_t count) { array_count_ += count; }
  void floatCountAdd(size_t count) { float_count_ += count; }
  void integerCountAdd(size_t count) { integer_count_ += count; }
  void stringCountAdd(size_t count) { string_count_ += count; }
  void charCountAdd(size_t count) { char_count_ += count; }

  // Notionally allocates data on the 5 data arrays and returns a data block
  // The function is to be used sequentially on all subscribed variables.
  // The final values are used to actually allocate memory in the InfoDataBlock object.

  // Set up the indexes and pre-allocate fixed sized fields (run once for all fields)
  [[nodiscard]] size_t staticIncrementAndAllocate(InfoEvidenceIntern internal_type);

  // Allocate additional memory space at runtime run for every Info field parsed.
  // Sets up all the array indexes and verifies the size of the total allocated memory.
  [[nodiscard]] bool dynamicIncrementAndAllocate(size_t data_size, InfoEvidenceIntern internal_type, const InfoParserToken &token);

private:

  size_t unity_array_count_{0};
  size_t array_count_{0};
  size_t float_count_{0};
  size_t integer_count_{0};
  size_t string_count_{0};
  size_t char_count_{0};

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Internal array index structure (can use 64 bits for size efficiency).

// An array block, 1 for each data array. VCF data arrays are generally small, so this can be a memory overhead.
// Small is 64 bits total, large is 2 x 64 = 128 bits, generally large is fine (unless space is tight and a lot of arrays).
//#define INFO_DATA_SMALL_ARRAY_SIZE 1  // Uncomment for 64 bits per array, else 128 bits.


struct InfoArrayIndex {

public:

  InfoArrayIndex() = default;
  InfoArrayIndex(size_t variable_index, size_t element_offset, size_t element_count) {

    infoVariableIndex(variable_index);
    infoOffset(element_offset);
    infoSize(element_count);

  }
  InfoArrayIndex(const InfoArrayIndex&) = default;
  ~InfoArrayIndex() = default;

  [[nodiscard]] size_t infoVariableIndex() const { return static_cast<size_t>(info_variable_index_); }
  [[nodiscard]] size_t infoOffset() const { return static_cast<size_t>(info_element_offset_); }
  [[nodiscard]] size_t infoSize() const { return static_cast<size_t>(info_element_count_); }

  void infoVariableIndex(size_t info_variable_index) { info_variable_index_ = static_cast<VariableIndexImpl>(info_variable_index); }
  void infoOffset(size_t info_element_offset) { info_element_offset_ = static_cast<InfoOffsetImpl>(info_element_offset); }
  void infoSize(size_t info_element_count) { info_element_count_ = static_cast<InfoCountImpl>(info_element_count); }


private:

#ifdef INFO_DATA_SMALL_ARRAY_SIZE
  // Minimize storage size (64 bits).
  using VariableIndexImpl = uint16_t;
  using InfoOffsetImpl = uint32_t;
  using InfoCountImpl = uint16_t;
#else
  using VariableIndexImpl = uint32_t;
  using InfoOffsetImpl = uint64_t;
  using InfoCountImpl = uint32_t;
#endif


  VariableIndexImpl info_variable_index_{0};  // corresponds to the variable index in the header.
  InfoOffsetImpl info_element_offset_{0};
  InfoCountImpl info_element_count_{0};

};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Resource allocation objects.
// The ResourceHandle object is granted to a resource requester.
// The handle is granted by the resource allocators below.

// FixedDynamic is pre-allocated and can change data size at runtime
enum class DataDynamicType { FixedData, FixedDynamic, DynamicData};
// Data class
enum class DataResourceType { Boolean, Integer, Float, String };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The info handle is held is held permanently by the InfoSubscribedField after subscribing to the Info data.
// The handle is constant and uses indirect indexes to access each instance of an Info Data Block.
class InfoResourceHandle {

public:

  InfoResourceHandle(size_t unique_handle_id,
                     size_t initial_data_offset,
                     size_t initial_data_size,
                     DataDynamicType dynamic_type,
                     DataResourceType resource_type) : unique_handle_id_(unique_handle_id),
                                                      initial_data_offset_(initial_data_offset),
                                                      initial_data_size_(initial_data_size),
                                                      dynamic_type_(dynamic_type),
                                                      resource_type_(resource_type) {}

  ~InfoResourceHandle() = default;

  [[nodiscard]] size_t handleId() const { return unique_handle_id_; }
  [[nodiscard]] size_t initialDataOffset() const { return initial_data_offset_; }
  [[nodiscard]] size_t initialDataSize() const { return initial_data_size_; }
  [[nodiscard]] DataDynamicType dynamicType() const { return dynamic_type_; }
  [[nodiscard]] DataResourceType resourceType() const { return resource_type_; }


private:

  size_t unique_handle_id_;
  size_t initial_data_offset_;
  size_t initial_data_size_;
  DataDynamicType dynamic_type_;
  DataResourceType resource_type_;


};

// The named resource objects.
// The fixed (known at definition time) resource object.
class FixedResourceInstance {

public:

  FixedResourceInstance(std::string resource_name, size_t resource_value) : resource_name_(std::move(resource_name)), resource_value_(resource_value) {}
  ~FixedResourceInstance() = default;

  [[nodiscard]] const std::string& resourceName() const { return resource_name_; }
  [[nodiscard]] size_t resourceValue() const { return resource_value_; }
  // returns the current value and increments the resource.
  size_t allocateResource(size_t requested) { size_t current = resource_value_; resource_value_+= requested; return current; }

private:

  std::string resource_name_;
  size_t resource_value_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The dynamic resource object, this is only known at runtime.
class DynamicResourceInstance {

public:

  DynamicResourceInstance(std::string resource_name,  std::vector<InfoArrayIndex> dynamic_allocation)
  : resource_name_(std::move(resource_name)), dynamic_allocation_(std::move(dynamic_allocation)) {}
  ~DynamicResourceInstance() = default;

  [[nodiscard]] const std::string& resourceName() const { return resource_name_; }
  [[nodiscard]] std::vector<InfoArrayIndex>& dynamicAllocation() { return dynamic_allocation_; }
  void queueResource(const InfoArrayIndex& requested) {  dynamic_allocation_.push_back(requested); }

private:

  std::string resource_name_;
  std::vector<InfoArrayIndex> dynamic_allocation_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple resource allocator
// Allocates requested memory blocks, both fixed and dynamic.
class InfoResourceAllocator {

public:

  InfoResourceAllocator( DataResourceType resource_type,
                         std::shared_ptr<FixedResourceInstance> unique_handle_ptr,
                         std::shared_ptr<FixedResourceInstance> fixed_resource_ptr,
                         std::shared_ptr<DynamicResourceInstance> dynamic_resource_ptr)
                        : resource_type_(resource_type),
                          unique_handle_ptr_(std::move(unique_handle_ptr)),
                          fixed_resource_ptr_(std::move(fixed_resource_ptr)),
                          dynamic_resource_ptr_(std::move(dynamic_resource_ptr)) {}
  ~InfoResourceAllocator() = default;

  [[nodiscard]] std::optional<const InfoResourceHandle> resourceRequest(DataResourceType resource_type, DataDynamicType dynamic_type, size_t data_size);

  [[nodiscard]] size_t allocateResource(size_t size) { return  fixed_resource_ptr_->allocateResource(size); }



private:

  const DataResourceType resource_type_;
  std::shared_ptr<FixedResourceInstance> unique_handle_ptr_;
  std::shared_ptr<FixedResourceInstance> fixed_resource_ptr_;
  std::shared_ptr<DynamicResourceInstance> dynamic_resource_ptr_;

  size_t uniqueHandleId() { return unique_handle_ptr_->allocateResource(1);  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The String resource allocator.
class InfoStringAllocator {

public:

  InfoStringAllocator( std::shared_ptr<FixedResourceInstance> unique_handle_ptr,
                       std::shared_ptr<FixedResourceInstance> view_resource_ptr,
                       std::shared_ptr<FixedResourceInstance> char_resource_ptr,
                       std::shared_ptr<DynamicResourceInstance> dynamic_resource_ptr)
    : unique_handle_ptr_(std::move(unique_handle_ptr)),
      view_resource_ptr_(std::move(view_resource_ptr)),
      char_resource_ptr_(std::move(char_resource_ptr)),
      dynamic_resource_ptr_(std::move(dynamic_resource_ptr))  {}
  ~InfoStringAllocator() = default;

  [[nodiscard]] std::optional<const InfoResourceHandle> resourceRequest(DataResourceType resource_type, DataDynamicType dynamic_type, size_t data_size);

  void allocateStringChars(size_t size, const InfoParserToken& token);
  [[nodiscard]] size_t allocateViews(size_t size) { return view_resource_ptr_->allocateResource(size); }

private:

  const DataResourceType resource_type_{ DataResourceType::String };
  std::shared_ptr<FixedResourceInstance> unique_handle_ptr_;
  std::shared_ptr<FixedResourceInstance> view_resource_ptr_;
  std::shared_ptr<FixedResourceInstance> char_resource_ptr_;
  std::shared_ptr<DynamicResourceInstance> dynamic_resource_ptr_;

  size_t uniqueHandleId() { return unique_handle_ptr_->allocateResource(1); }

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object manges the memory for the Info Data Block.

class InfoMemoryResource {

public:

  InfoMemoryResource();
  InfoMemoryResource(const InfoMemoryResource& copy);
  InfoMemoryResource& operator=(const InfoMemoryResource& copy);
  ~InfoMemoryResource() = default;

// The initial resource request from a subscribed data item (this is only done once).
  [[nodiscard]] std::optional<const InfoResourceHandle> resourceRequest( DataResourceType resource_type, DataDynamicType dynamic_type, size_t data_size);
// FixedDynamic fields can request a dynamic data block if runtime data size mis-matches pre-allocated data size.
  void requestDynamic(size_t requestor_id) { array_memory_->queueResource(InfoArrayIndex(requestor_id, 0, 0)); }
// Resolves static data sizes with runtime data sizes.
  bool resolveAllocation(const InfoResourceHandle& item_resource_handle, const InfoParserToken& token);
// Raw memory audit functions.
  [[nodiscard]] size_t charSize() const { return char_memory_->resourceValue(); }
  [[nodiscard]] size_t integerSize() const { return integer_memory_->resourceValue(); }
  [[nodiscard]] size_t floatSize() const { return float_memory_->resourceValue(); }
  [[nodiscard]] size_t viewSize() const { return view_memory_->resourceValue(); }
  [[nodiscard]] size_t arraySize() const { return array_memory_->dynamicAllocation().size(); }

private:

  // Used to generate a unique identifier for resource requesters.
  std::shared_ptr<FixedResourceInstance> unique_ident_;    // The unique identifier resource.

  // These objects keep track of allocated memory.
  // The data allocation map is contained in these objects.
  std::shared_ptr<FixedResourceInstance> char_memory_;    // The number of characters allocated in the Info block
  std::shared_ptr<FixedResourceInstance> integer_memory_;    // The size of the integer array in the Info block
  std::shared_ptr<FixedResourceInstance> float_memory_;    // The size of the floating point array in the Info block
  std::shared_ptr<FixedResourceInstance> view_memory_;    // The size of the std::string_view array in the Info block
  std::shared_ptr<DynamicResourceInstance> array_memory_;    // The size of the array lookup in the Info block

  // These objects keep track of allocated data types and translate to a memory footprint.
  // These objects have a virtual view of allocation. The actual resource allocation is held in the ResourceInstance
  // objects held above.
  std::unique_ptr<InfoResourceAllocator> bool_allocator_;
  std::unique_ptr<InfoResourceAllocator> integer_allocator_;
  std::unique_ptr<InfoResourceAllocator> float_allocator_;
  std::unique_ptr<InfoStringAllocator> string_allocator_;

  bool resolveDynamic(const InfoResourceHandle& item_resource_handle, const InfoParserToken& token);
  bool findUpdateDynamic(const InfoResourceHandle& item_resource_handle, const InfoParserToken& token);


};



} //namespace.


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEMORY_H
