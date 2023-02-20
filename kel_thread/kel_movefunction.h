// Copyright 2023 Kellerberrin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
// OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//

#ifndef KEL_MOVEFUNCTION_H
#define KEL_MOVEFUNCTION_H


#include <functional>


namespace kellerberrin {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A simple implementation of a std::function object for std::move_constructable callable objects.
// Callable objects are typically std::move_constructable if they have stateful arguments that are std::move_constructable.
// Thus the callable object returned from 'std::bind_front(std::forward<F>(f), std::forward<Args>(args)...)'
// where one or more of the args in 'std::forward<Args>(args)...' is a 'unique_ptr<>' is std::move_constructable.
// Typically, pushing and poping elements off the thread safe queues is most efficiently done if they are held by 'std::unique_ptr<>'.
//
// Note that this object will become obsolete when std::move_from_function becomes available in C++ 23.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> class MoveFunction; // Not defined.

template<typename ReturnType, typename... Args> class MoveFunction<ReturnType(Args...)> {

  struct InvokeBase {

    InvokeBase() = default;
    virtual ~InvokeBase() = default;

    virtual ReturnType invoke(Args &&... args) = 0;

  }; // InvokeBase

  template<typename Callable> struct InvokeImpl : InvokeBase {

    explicit InvokeImpl(Callable callable) : callable_(std::move(callable)) {}
    ~InvokeImpl() override = default;

    [[nodiscard]] ReturnType invoke(Args &&... args) override { return std::invoke(callable_, std::forward<Args>(args)...); }

    Callable callable_;

  }; // InvokeImpl

public:

  template<typename CallableObject> explicit MoveFunction(CallableObject callable_object)
      : callable_ptr_(std::make_unique<InvokeImpl<CallableObject>>(std::move(callable_object))) {}

  MoveFunction() = default;
  MoveFunction(const MoveFunction &) = delete;
  MoveFunction(MoveFunction &&function) noexcept: callable_ptr_(std::move(function.callable_ptr_)) {}
  ~MoveFunction() = default;

  [[nodiscard]] ReturnType operator()(Args &&... args) { return callable_ptr_->invoke(std::forward<Args>(args)...); }

private:

  std::unique_ptr<InvokeBase> callable_ptr_;

};


} // namespace


#endif //KEL_MOVEFUNCTION_H
