//
// Created by kellerberrin on 14/2/21.
//

#ifndef KEL_MEM_ALLOC_H
#define KEL_MEM_ALLOC_H


#include <memory_resource>


namespace kellerberrin {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Overload allocation to monitor memory allocation/de-allocation
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MonitorAllocate : public std::pmr_memory_resource {

  MonitorAllocate() = default;
  ~MonitorAllocate() override = default;

  [[nodiscard]] allocatedBytes() const { return allocated_bytes_; }
  [[nodiscard]] deallocatedBytes() const { return deallocated_bytes_; }

private:

  [[nodiscard]] void *do_allocate(size_t bytes, size_t align) override {

    allocated_bytes_ += bytes;

    return std::pmr::new_delete_resource()->allocate(bytes, align);

  }

  void do_deallocate(void* ptr, size_t bytes, size_t align) override {

    deallocated_bytes_ += bytes;

    std::pmr::new_delete_resource()->deallocate(ptr, bytes, align);

  }

  [[nodiscard]] bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {

    return std::pmr::new_delete_resource()->is_other(other);

  }


  std::atomic<size_t> allocated_bytes_{0};
  std::atomic<size_t> deallocated_bytes_{0};

};


} // namespace.


#endif //KEL_MEM_ALLOC_H
