// File: InvPosition.H Define a 1/x^n operator for the BSpline package.

#ifndef BSPLINE_OPERATORS_INV_POSITION_H
#define BSPLINE_OPERATORS_INV_POSITION_H
 
#include <bspline/internal/misc.h>
#include <bspline/operators/GenericOperators.h>
 
namespace bspline::operators {
namespace internal = bspline::internal;
 
template <size_t n>
class InvPosition final : public Operator {
 private:
  template <typename T>
  static std::array<T, n + 1> expandPower(const T &xm) {
    std::array<T, n + 1> retVal;
    T power_of_xm = static_cast<T>(1);
    for (size_t i = 0; i < n + 1; i++) {
      retVal[n - i] = internal::binomialCoefficient<T>(n, i) * power_of_xm;
      power_of_xm /= xm;
    }
    return retVal;
  }
 
 public:
  static constexpr size_t outputOrder(size_t inputOrder) {
    return inputOrder + n;
  }
 
  template <typename T, size_t size>
  std::array<T, outputOrder(size - 1) + 1> transform(
      const std::array<T, size> &input, const support::Grid<T> &grid,
      size_t intervalIndex) const {
    constexpr size_t OUTPUT_SIZE = size + n;
 
    const T xm =
        (grid[intervalIndex] + grid[intervalIndex + 1]) / static_cast<T>(2);
 
    const std::array<T, n + 1> expanded = expandPower<T>(xm);
 
    std::array<T, OUTPUT_SIZE> retVal;
    retVal.fill(static_cast<T>(0));
 
    for (size_t i = 0; i < input.size(); i++) {
      for (size_t j = 0; j < expanded.size(); j++) {
        retVal[i + j] += expanded[j] * input[i];
      }
    }
    return retVal;
  }
};
 
template <size_t n>
using InvX = InvPosition<n>;
 
}  // namespace bspline::operators
#endif  // BSPLINE_OPERATORS_INV_POSITION_H