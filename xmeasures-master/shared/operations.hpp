//! \brief Basic common operations.
//! The Dao (Deterministic Agglomerative Overlapping) of Clustering library:
//! Robust & Fine-grained Deterministic Clustering for Large Networks.
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2015-08-06

#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <cstddef>  // ptrdiff_t
#include <cmath>  // fabs(), isinf()
#include <limits>  // Type limits
#include <type_traits>  // type check, enable_if_t
#include <cassert>
#include <utility>  // declval, move, forward
#include <iterator>  // iterator_traits
#include <stdexcept>  // Exception (for Arguments processing)
#include <string>  // to_string

#include "macrodef.h"  // TRACE, VALIDATE

#if TRACE >= 2
#include <cstdio>
#endif // TRACE

namespace daoc {

// Internal MACROSES:
//	- CONFIDENCE  - confidence of the operating values for the comparison
//	operations to work in the same confidence range. Float, E (0, 1]
//  NOTE: Affects only less() and equal() operations, not the memory address
//  related ordering.
//	- INSORTED_NONUNIQUE  - insorted() does not validate that the inserting
//	element is not present
//
// NOTE:
// - undefined maro definition is interpreted as having value 0
// - constexpr are noexcept since the noexcept operator always returns true for
//	a constant expression

using std::ptrdiff_t;
using std::numeric_limits;
using std::iterator_traits;
using std::is_same;
using std::is_floating_point;
using std::is_integral;
using std::is_scalar;
using std::is_reference;
using std::is_base_of;
using std::is_const;
using std::is_pointer;
using std::declval;
using std::max;
using std::min;
using std::enable_if_t;
using std::conditional_t;
using std::remove_pointer_t;
using std::distance;
using std::invalid_argument;
using std::string;


// Global parameters ----------------------------------------------------------
//! Confidence of the operating values for the comparison operations to work in
//! the same confidence range. Float, E (0, 1].
//! NOTE: Affects only less() and equal() operations, not the memory address related
//! ordering.
#ifndef CONFIDENCE
	#define CONFIDENCE  1
#endif // CONFIDENCE
static_assert(CONFIDENCE > 0 && CONFIDENCE <= 1, "Macrodef, values CONFIDENCE should E (0, 1]");
using ConfidanceInt = uint32_t;  //!< Integral type of Confidence
//! Confidence Denominator, Max possible value of the confidence
static constexpr ConfidanceInt  CONFIDENCE_DEN = 1E9;
//! Confidence Multiplier, E [1, CONFIDENCE_DEN]
static constexpr ConfidanceInt  CONFIDENCE_MUL = ceil(CONFIDENCE * CONFIDENCE_DEN);

// Tracing file
//!\brief Global trace log
//! Client may reassign ftrace, for example when the hierarchy is output to the
//! stdout so as this channel should not be used for the tracing.
//! \note Mainly stdout or stderr are assumed to used
#ifdef FTRACE_GLOBAL
extern FILE* ftrace;
#else
static FILE* ftrace =
// Note: TRACE <= 1 has vary few output
#if defined(DEBUG) || TRACE <= 1
	stderr
#else
	stdout  // Release
#endif // DEBUG, TRACE
	;
#endif // FTRACE_GLOBAL

// Compile time tracing functions ---------------------------------------------
//! \brief Trace specified type in compile time
//!
//! \return void
template <typename T>
constexpr void traceType()
{
	static_assert(T()._, "traceType() termination");
	// Or use RTTI typeid()
}

// Meta Definitions of Base Types ---------------------------------------------
//! \brief Type T by value or const reference depending of the type size
// Note: sizeof(T&) == sizeof(T)
template <typename T>
using ValCRef = conditional_t<sizeof(T) <= sizeof(void*) || is_scalar<T>::value, T, const T&>;

//! \brief Member type by value or const reference depending of the type size
template <typename ItemsT>
using MemberValCRef = ValCRef<typename ItemsT::value_type>;

//! \brief Fetches (depointered) type of the iterator value
//!
//! \tparam IterT  - iterator type of either a value or value pointer type
template <typename IterT>
using IterPureBaseT = remove_pointer_t<typename iterator_traits<IterT>::value_type>;

// Accessory type identifiers -------------------------------------------------
//! \brief Whether the type is an iterator of the specified category
//! \tparam T  - evaluating type
//! \tparam Category=std::input_iterator_tag  - target base category to be checked
//!
//! \return constexpr bool  - the type is iterator
template <typename T, typename Category=std::input_iterator_tag>
constexpr enable_if_t<is_base_of<Category, typename iterator_traits<T>::iterator_category>::value
	, bool> is_iterator()
{ return true; }

template <typename T, typename Category=std::input_iterator_tag>
constexpr enable_if_t<!is_base_of<Category, typename iterator_traits<T>::iterator_category>::value
	, bool> is_iterator()
{ return false; }

template <typename T, typename Category=std::input_iterator_tag>
constexpr bool is_iterator_v = is_iterator<T, Category>();


//! \brief Whether the type is an unordered associative container
//! \tparam T  - evaluating type
template <typename T, typename Accessory=void>
struct is_hashContainer: std::false_type  {};

template <typename T>
struct is_hashContainer<T, enable_if_t<std::is_object<typename T::hasher>::value>>
: std::true_type  {};

template <typename T>
constexpr bool is_hashContainer_v = is_hashContainer<T>::value;

// Accessory conditional operators --------------------------------------------
//! \brief Dereference iterator or its pointer to the underlying value or its reference
//!
//! \tparam ItemXT  - scalar value or pointer type
//!
//! \param el  - value or its pointer to be dereferenced
//! \return el Item  - resulting value or its reference
template <typename ItemXT>
constexpr enable_if_t<is_pointer<ItemXT>::value, remove_pointer_t<ItemXT>&> derefItemX(ItemXT el)
{
	static_assert(sizeof(ItemXT) <= sizeof(void*), "derefItemX(), ItemXT type size validation failed");
	return *el;
}

template <typename ItemXT>
constexpr enable_if_t<!is_pointer<ItemXT>::value, ItemXT> derefItemX(ItemXT el)
{
	static_assert(sizeof(ItemXT) <= sizeof(void*), "derefItemX(), ItemXT type size validation failed");
	return el;
}

// Accessory math functions ---------------------------------------------------
//! \brief Precision limit (error) for the specified floating type
// NOTE: at least GCC provides numeric_limits<ValT>::epsilon() / ValT(2) > 0
template <typename ValT>
constexpr ValT precision_limit()
{
	static_assert(is_floating_point<ValT>::value, "precision_limit() type should be fractional");
	// Note: much lower epsilon does not impact convergence,
	// but usage of sqrt(epsilon) allows to skip the size in comparison functions
	// ATTENTION: sqrt(eps) causes inclusion into mcands non-max neighbors in large dense networks
	// ATTENTION: values > eps cause lose of accuracy and order-dependence
	//return numeric_limits<ValT>::epsilon() * 5;  // 5 is half of the base 10
	//return numeric_limits<ValT>::epsilon() * 2;  // *2 to consider arithmetic error of +/- operations
	// ATTENTION: accurate results are provided by val ~< eps:
	//return pow(numeric_limits<ValT>::epsilon(), 1.25);  // *2 to consider arithmetic error of +/- operations
	return numeric_limits<ValT>::epsilon() / 2;  // /2 to have less the distinguishable value
}

//! \brief Strict less for floating point numbers
//! \note Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
//! 	CONFIM < CONFIDENCE_DEN provides less strict (looser) equal and stricter less comparisons
//!
//! \param a ValT  - val1
//! \param b ValT  - val2
//! \return bool  - val1 < val2
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_floating_point<ValT>::value, bool>
less(const ValT a, const ValT b=ValT(0))
{
	static_assert(is_floating_point<ValT>::value, "less(), value type should be fractional");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "less(), values CONFIDENCE should E (0, 1]");
	// a < b   <= For x > 0:  a + x < b =>  For y > 0:  a - b + y * fabs(a+b) < 0
	// ATTENTION:
	// - "precision_limit" must be the first multiplier to have correct evaluation for huge values
	// - fabs(a + b) should not be used in the normalization, it diminishes accuracy
	// Exact solution:
	//return a - b + precision_limit<ValT>() * max(fabs(a), fabs(b)) < 0;
	// NOTE: Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
	// NOTE: CONFID multiplier in the end to not affect accuracy if sizeof(ValT) > sizeof(CONFID)
	// NOTE: either <= with precision_limit<ValT>() or < with epsilon<ValT>()
	// Faster computation with sufficient accuracy:
	return (CONFIM < CONFIDENCE_DEN ? static_cast<ValT>(2 * CONFIM) / CONFIDENCE_DEN : 2)
		* (a - b) / (fabs(a) + fabs(b) + precision_limit<ValT>())
		+ precision_limit<ValT>() < 0;
}

//! \brief Strict less for integral numbers
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_integral<ValT>::value, bool>
less(const ValT a, const ValT b=ValT(0))
{
	static_assert(is_integral<ValT>::value, "less(), value type should be integral");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "less(), values CONFIDENCE should E (0, 1]");
	return CONFIM == CONFIDENCE_DEN ? a < b : less<double, CONFIM>(a, b);
}

//! \brief Strict less for floating point numbers
//! \pre size should be positive (!= 0)
//! \note Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
//! 	CONFIM < CONFIDENCE_DEN provides less strict (looser) equal and stricter less comparisons
//!
//! \param a ValT  - val1
//! \param b ValT  - val2
//! \param size=1 float  - average number of elements in a and b
//! 	to consider accumulation error
//! \return bool  - val1 < val2
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_floating_point<ValT>::value, bool>
lessx(const ValT a, const ValT b=ValT(0), const float size=1)
{
	static_assert(is_floating_point<ValT>::value, "less(), value type should be fractional");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "less(), values CONFIDENCE should E (0, 1]");
#if VALIDATE >= 2
	assert(size && "less(), positive size is expected for the imprecise values");
#endif // VALIDATE
	// a < b   <= For x > 0:  a + x < b =>  For y > 0:  a - b + y * fabs(a+b) < 0
	// ATTENTION:
	// - "1 +" is required to handle accumulation error correctly for float numbers < 1
	// - "precision_limit" must be the first multiplier to have correct evaluation for huge values
	// - size should not be used as a direct multiplier, because for the huge number of items
	// 	it will diminish precision and make the comparison invalid
	// - fabs(a + b) should not be used in the normalization, it diminishes accuracy
	// Exact solution:
	//return a - b + precision_limit<ValT>() * (1 + log2f(size)) * max(fabs(a), fabs(b)) < 0;
	// NOTE: Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
	// NOTE: CONFID multiplier in the end to not affect accuracy if sizeof(ValT) > sizeof(CONFID)
	// NOTE: either <= with precision_limit<ValT>() or < with epsilon<ValT>()
	// Faster computation with sufficient accuracy:
	return (CONFIM < CONFIDENCE_DEN ? static_cast<ValT>(2 * CONFIM) / CONFIDENCE_DEN : 2)
		* (a - b) / (fabs(a) + fabs(b) + precision_limit<ValT>())
		+ precision_limit<ValT>() * (1 + log2f(size)) < 0;
}

//! \brief Strict less for integral numbers
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_integral<ValT>::value, bool>
lessx(const ValT a, const ValT b=ValT(0), const float size=1)
{
	static_assert(is_integral<ValT>::value, "less(), value type should be integral");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "less(), values CONFIDENCE should E (0, 1]");
	return CONFIM == CONFIDENCE_DEN ? a < b : less<double, CONFIM>(a, b);
}

//! \brief Check equality of 2 floating point numbers
//! \note Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
//! 	CONFIM < CONFIDENCE_DEN provides less strict (looser) equal and stricter less comparisons
//!
//! \param a ValT  - val1
//! \param b ValT  - val2
//! \return bool  - whether val1 ~= val2
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_floating_point<ValT>::value, bool>
equal(const ValT a, const ValT b=ValT(0))
{
	static_assert(is_floating_point<ValT>::value, "equal(), value type should be fractional");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "equal(), values CONFIDENCE should E (0, 1]");
//#if VALIDATE >= 2
//#if TRACE >= 3
//	fprintf(ftrace, "equal()  a: %G, b: %G, precLim: %G (from %G)"
//		"\n\tres: %d, val: %G, lim: %G (lmul: %G)\n"
//		, a, b, precision_limit<ValT>(), numeric_limits<ValT>::epsilon()
//		, fabs(a - b) <= precision_limit<ValT>() * (1 + fabs(a) + fabs(b))
//		, fabs(a - b), precision_limit<ValT>() * (1 + fabs(a) + fabs(b))
//		, (1 + fabs(a) + fabs(b)));
//#endif // TRACE_EXTRA
//#endif // VALIDATE
	// ~ a == b with strictly defined equality margin
	// ATTENTION:
	// - "1 +" is required to handle accumulation error correctly for float numbers < 1
	// - "precision_limit" must be the first multiplier to have correct evaluation for huge values
	//// - "fabs(a + b)" is not appropriate here in case of a = -b && b > 0  => fabs(a) + fabs(b)
	//return fabs(a - b) <= precision_limit<ValT>() * max(fabs(a), fabs(b));
	// NOTE: Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
	// NOTE: CONFID multiplier in the end to not affect accuracy if sizeof(ValT) > sizeof(CONFID)
	// NOTE: either <= with precision_limit<ValT>() or < with epsilon<ValT>()
	return (CONFIM < CONFIDENCE_DEN ? static_cast<ValT>(2 * CONFIM) / CONFIDENCE_DEN : 2)
		* fabs(a - b) / (fabs(a) + fabs(b) + precision_limit<ValT>())
		<= precision_limit<ValT>();
}

//! \brief Strict equal for integral numbers
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_integral<ValT>::value, bool>
equal(const ValT a, const ValT b=ValT(0))
{
	static_assert(is_integral<ValT>::value, "equal(), value type should be integral");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "equal(), values CONFIDENCE should E (0, 1]");
	return CONFIM == CONFIDENCE_DEN ? a == b : equal<double, CONFIM>(a, b);
}

//! \brief Check equality of 2 floating point numbers
//! \pre size should be positive (!= 0)
//! \note Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
//! 	CONFIM < CONFIDENCE_DEN provides less strict (looser) equal and stricter less comparisons
//!
//! \param a ValT  - val1
//! \param b ValT  - val2
//! \param size=1 float  - average number of elements in a and b
//! 	to consider accumulation error
//! \return bool  - whether val1 ~= val2
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_floating_point<ValT>::value, bool>
equalx(const ValT a, const ValT b=ValT(0), const float size=1)
{
	static_assert(is_floating_point<ValT>::value, "equal(), value type should be fractional");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "equal(), values CONFIDENCE should E (0, 1]");
#if VALIDATE >= 2
	assert(size && "equal(), positive size is expected for the imprecise values");
#endif // VALIDATE
//#if VALIDATE >= 2
//	// Allow stand-alone nodes
//	//assert(size >= 1 && "equal(), size should be positive");
//#if TRACE >= 3
//	fprintf(ftrace, "equal()  a: %G, b: %G, cnt: %lu, precLim: %G (from %G)"
//		"\n\tres: %d, val: %G, lim: %G (lmul: %G)\n"
//		, a, b, size, precision_limit<ValT>(), numeric_limits<ValT>::epsilon()
//		, fabs(a - b) <= precision_limit<ValT>() * size * (1 + fabs(a) + fabs(b))
//		, fabs(a - b), precision_limit<ValT>() * size * (1 + fabs(a) + fabs(b))
//		, (1 + fabs(a) + fabs(b)));
//#endif // TRACE_EXTRA
//#endif // VALIDATE
	// ~ a == b with strictly defined equality margin
	// ATTENTION:
	// - "1 +" is required to handle accumulation error correctly for float numbers < 1
	// - "precision_limit" must be the first multiplier to have correct evaluation for huge values
	// - size should not be used as a direct multiplier, because for the huge number of items
	// 	it will diminish precision and make the comparison invalid
	//// - "fabs(a + b)" is not appropriate here in case of a = -b && b > 0  => fabs(a) + fabs(b)
	//return fabs(a - b) <= precision_limit<ValT>() * (1 + log2f(size)) * max(fabs(a), fabs(b));
	// NOTE: Exact Evaluations with Floating Point Numbers: https://goo.gl/A1DSwn
	// NOTE: CONFID multiplier in the end to not affect accuracy if sizeof(ValT) > sizeof(CONFID)
	// NOTE: either <= with precision_limit<ValT>() or < with epsilon<ValT>()
	return (CONFIM < CONFIDENCE_DEN ? static_cast<ValT>(2 * CONFIM) / CONFIDENCE_DEN : 2)
		* fabs(a - b) / (fabs(a) + fabs(b) + precision_limit<ValT>())
		<= precision_limit<ValT>() * (1 + log2f(size));
}

//! \brief Strict equal for integral numbers
template <typename ValT, ConfidanceInt CONFIM=CONFIDENCE_MUL>
constexpr enable_if_t<is_integral<ValT>::value, bool>
equalx(const ValT a, const ValT b=ValT(0), const float size=1)
{
	static_assert(is_integral<ValT>::value, "equal(), value type should be integral");
	static_assert(CONFIM > 0 && CONFIM <= CONFIDENCE_DEN, "equal(), values CONFIDENCE should E (0, 1]");
	return CONFIM == CONFIDENCE_DEN ? a == b : equal<double, CONFIM>(a, b);
}

// Comparison functions -------------------------------------------------------
//! Comparison function for the arguments
template <typename ValT>
constexpr bool cmpBase(const ValT a, const ValT b)
{
//	using ValT = common_type<Val1T, Val2T>::type;
//	// Note: scalar includes pointer
//	static_assert((is_scalar<ValT>::value || is_reference<ValT>::value)
	static_assert(!is_floating_point<ValT>::value
		&& (sizeof(ValT) <= sizeof(void*)  || is_reference<ValT>::value)
		, "cmpBase(), arguments should be integral and either small objects or references");
	return a < b;
}

//! Comparison function for objects sorting by dest attribute
//! \note Uses cmpBase()
template <typename T>
constexpr bool cmpDest(const T& a, const T& b)
{ return cmpBase(a.dest, b.dest); }

//! Comparison function for the objects weight
template <typename ObjsT, typename WeightT>
constexpr bool cmpObjsWeight(typename ObjsT::const_reference a, typename ObjsT::const_reference b)
// Note: strict < affects the results making them non-deterministic
{ return less<WeightT>(a.weight, b.weight); }

// Binary search --------------------------------------------------------------
// Min size of an array when bin search is faster than a linear scan, 11 or 9
constexpr int  BINSEARCH_MARGIN = 11;

//! \brief Direct value comparison (binary find callback)
//! 	Note: Items in the container are assumed to be ordered in acs order
//! \attention must be synced with cmpBase()
//!
//! \param cv T  - container value
//! \param v T  - control value
//! \return ptrdiff_t  - comparison result: cv - v
template <typename T>
constexpr enable_if_t<sizeof(T) <= sizeof(void*), ptrdiff_t>
bsVal(const T cv, const T v)
{
	//using ValT = common_type<Val1T, Val2T>::type;
	static_assert(!is_floating_point<T>::value
	//	//&& sizeof(T) <= sizeof(ptrdiff_t)
		, "bsVal(), argument type constraints failed");
	// ATTENTION: ptrdiff_t type cast gives invalid comparison when the second value
	// is large negative (high memory addresses) and the fist one is large positive
	//return reinterpret_cast<ptrdiff_t>(cv) - reinterpret_cast<ptrdiff_t>(v);
	return (cv != v) * ((cv < v) * -2 + 1);
}

template <typename T>
constexpr enable_if_t<(sizeof(T) > sizeof(void*)), ptrdiff_t>
bsVal(const T& cv, const T& v)
{
	static_assert(!is_floating_point<T>::value
		, "bsVal(), argument type constraints failed");
	return (cv != v) * ((cv < v) * -2 + 1);
}

//template <typename T>
//constexpr enable_if_t<is_floating_point<T>::value, ptrdiff_t>
//bsVal(const T cv, const T v, const size_t size=1)
//{
//	return !equal(cv, v, size) * (less(cv, v, size) * -2 + 1);
//}

////! \brief Dereferenced value address comparison (binary find callback)
////! 	Useful to process a container of iterators
////! 	Note: Items in the container are assumed to be ordered in acs order
////! \attention must be synced with cmpBase()
////!
////! \param cv T  - container value
////! \param v T  - control value
////! \return ptrdiff_t  - comparison result: cv - v
//template <typename T>  //, enable_if_t<sizeof(T) <= sizeof(ptrdiff_t)>* = nullptr>
//constexpr ptrdiff_t bsDerefAddr(const T& cv, const T& v)
//// ATTENTION: without casting difference will be incorrect for the types with size
//// different from ptrdiff_t and for the unsigned types
//{
//	static_assert(sizeof(T) > sizeof(ptrdiff_t)
//		, "bsDerefAddr(), heavy dereferencable type is expected");
//	return bsVal(&*cv, &*v);
//}

//! Binary sort comparison function for objects sorting by dest attribute
//! \note Uses bsVal()
//!
//! \param cv T  - container value
//! \param v T  - control value
//! \return ptrdiff_t  - comparison result
template <typename T>
constexpr ptrdiff_t bsDest(const T& a, const T& b)
{ return bsVal(a.dest, b.dest); }

//! \brief Elements dest comparison (binary find callback)
//! \note Items in the container are assumed to be ordered in acs order by bsVal().
//!
//! \tparam T  - container type
//! \param obj typename T::const_reference  - object
//! \param dest const typename T::value_type::DestT*  - dest attribute
//! \return ptrdiff_t  - comparison result
template <typename T>
constexpr ptrdiff_t bsObjsDest(typename T::const_reference obj
	, enable_if_t<(sizeof(typename T::value_type::DestT) > sizeof(void*))
		, const typename T::value_type::DestT*> dest)
{ return bsVal<const typename T::value_type::DestT*>(obj.dest, dest); }

//! \copydoc bsObjsDest<typename T>
template <typename T>
constexpr ptrdiff_t bsObjsDest(typename T::const_reference obj
	, enable_if_t<sizeof(typename T::value_type::DestT) <= sizeof(void*)
		, const typename T::value_type::DestT> dest)
{ return bsVal<const typename T::value_type::DestT>(obj.dest, dest); }

//! \brief Elements operator() comparison (binary find callback)
//! \note Items in the container are assumed to be ordered in acs order by bsVal().
//!
//! \tparam T  - container type
//! \param obj typename T::const_reference  - object
//! \param val const typename T::value_type::CallT  - comparing value
//! \return ptrdiff_t  - comparison result: cv - v
template <typename T>
constexpr ptrdiff_t bsObjsOp(typename T::const_reference obj
	, const typename T::value_type::CallT val)  // noexcept(typename T::value_type()())
// ATTENTION: without casting difference will be incorrect for the types with size
// different from ptrdiff_t and for the unsigned types
{
	static_assert(sizeof(typename T::value_type::CallT) <= sizeof(ptrdiff_t)
		, "bsObjsOp(), argument type constraints failed");
	return bsVal<const typename T::value_type::CallT>(obj(), val);
}

//! \brief Element operator() comparison (binary find callback)
//! \note Items in the container are assumed to be ordered in acs order by bsVal().
//!
//! \tparam T  - object type
//! \param obj const T&  - object
//! \param val const typename T::CallT  - comparing value
//! \return ptrdiff_t  - comparison result: cv - v
template <typename T>
constexpr ptrdiff_t bsObjOp(const T& obj, const typename T::CallT val)
// ATTENTION: without casting difference will be incorrect for the types with size
// different from ptrdiff_t and for the unsigned types
{
	static_assert(sizeof(typename T::CallT) <= sizeof(ptrdiff_t)
		, "bsObjOp(), argument type constraints failed");
	return bsVal<const typename T::CallT>(obj(), val);
}

////! \brief Element dest comparison (binary find callback)
////! \note Items in the container are assumed to be ordered in acs order by bsVal().
////!
////! \param obj typename T::const_reference  - object
////! \param dest const typename T::value_type::DestT*  - dest attribute
////! \return ptrdiff_t  - comparison result
//template <typename T>
//constexpr ptrdiff_t bsObjDest(T obj
//	, const typename T::DestT* dest)
//{ return bsVal<const typename T::DestT*>(obj.dest, dest); }

//! \brief Binary sort comparison function for objects sorting by get() member function
//! \note Uses bsVal()
//!
//! \param cv T  - container value
//! \param v T  - control value
//! \return ptrdiff_t  - comparison result
template <typename T>
constexpr ptrdiff_t bsObjGetF(const T& cv, const T& v)
{ return bsVal(cv.get(), v.get()); }

//! \brief Lambda function of struct comparison by the attr
//! \pre The attribute should be an object of non-floating point type
//!
//! \param T  - type of the struct to be compared
//! \param attr  - comparison field (attribute)
#define BS_STRUCT(T, attr)  \
	[](const T& cv, const T& v) noexcept -> ptrdiff_t \
	{  \
		static_assert(!is_floating_point<decltype(cv.attr)>::value, \
			"BS_STRUCT, non-floating point attribute type is expected"); \
		return bsVal<ValCRef<decltype(cv.attr)>>(cv.attr, v.attr);  \
	}
// 	~ bsVal<decltype(v.attr)>(cv.attr, v.attr);  // But bsVal complicates compile-time error analysis

//! \brief Lambda function of struct binary search the attr
//! \pre The attribute should be non-floating point object
//!
//! \param T  - type of the struct
//! \param attr  - comparison field (attribute)
#define BS_OBJ_ATTR(T, attr)  \
	[](const T& el, const ValCRef<decltype(declval<T>().attr)> ea) noexcept -> ptrdiff_t \
	{  \
		static_assert(!is_floating_point<decltype(el.attr)>::value, \
			"BS_OBJ_ATTR, non-floating point attribute type is expected"); \
		return bsVal<ValCRef<decltype(el.attr)>>(el.attr, ea);  \
	}

//! \brief Binary search in the ordered range with random iterators
//! 	for the closest following elemnt
//! \attention Te returned iterator is not necessary point to the end if the item is not present
//! \pre The elements are ordered by bsVal()
//! \post Returned iterator always has not less value than required
//!
//! \tparam RandIT  - random iterator
//! \tparam CompareF  - NONSTANDARD comparison function:
//! 		int cmp(const RandIT::reference iterVal, const T val);
//! 			returns negative if iterVal < val, 0 if iterVal == val
//! 			and positive otherwise
//! \return random iterator of the value position if exists, otherwise
//! 	the iterator on the following value / end
template <typename RandIT, typename Val, typename CompareF>
enable_if_t<is_iterator<RandIT, std::random_access_iterator_tag>(), RandIT>
binary_ifind(RandIT begin, const RandIT end, const Val val
	, CompareF cmp=bsVal<MemberValCRef<RandIT>>) noexcept
{
	static_assert(sizeof(val) <= sizeof(void*)
		, "binary_ifind(), T should not exceed the address type size");
	static_assert(is_same<decltype(cmp(*begin, val)), decltype(bsVal(nullptr, nullptr))>::value
		, "binary_ifind(), cmp() must return the same type as bsVal()");
	static_assert(is_iterator<RandIT, std::random_access_iterator_tag>()
		, "binary_ifind(), RandIT must be a random iterator type");
	RandIT iend = end;
	for(RandIT  pos = begin + (iend - begin) / 2; pos != iend
	; pos = begin + (iend - begin) / 2) {
		int  cres = cmp(*pos, val);

		if(cres < 0) {
			begin = pos + 1;
			continue;
		}
		if(cres > 0) {
			iend = pos;
			continue;
		}
		return pos;

		//if(!((cres < 0 && ((begin = pos + 1), true)) || (cres > 0 && ((iend = pos), true))))
		//	return pos;
	}
#if VALIDATE >= 2
//#if TRACE >= 3
//	if(!(iend == end || cmp(*iend, val) <= 0)) {
//		// Note: iend != end is always true here
//		fprintf(ftrace, "  >>>>> binary_ifind(), iend: %s, val: %s. Container: "
//			, iend->idstr().c_str(), val->idstr().c_str());
//		for(auto ic = begin; ic != end; ++ic)
//			fprintf(ftrace, " %s", ic->idstr().c_str());
//		fputc('\n', ftrace);
//	}
//#endif // TRACE
	assert((iend == end || cmp(*iend, val) >= 0)
	&& "binary_ifind(), iterator verification failed");
#endif // VALIDATE
	return iend;
}

//! Binary search in the ordered range with random iterators
template <typename RandIT, typename Val, typename CompareF>
inline enable_if_t<is_iterator<RandIT, std::random_access_iterator_tag>(), RandIT>
binary_find(RandIT begin, const RandIT end, const Val val
	, CompareF cmp=bsVal<MemberValCRef<RandIT>>) noexcept
{
	//static_assert(is_iterator<RandIT, std::random_access_iterator_tag>()
	//	, "binary_find(), RandIT must be a random iterator type");
	static_assert(sizeof(Val) <= sizeof(void*)
		, "binary_find(), T should not exceed the address type size");
	begin = binary_ifind(begin, end, val, cmp);
	return begin != end && !cmp(*begin, val) ? begin : end;
}

// Linear search --------------------------------------------------------------
//! \brief Linear search in the specified range for the closest following element
//! \pre The elements are ordered by bsVal().
//! 	The number of elements is small (typically up to BINSEARCH_MARGIN)
//!
//! \tparam IT  - iterator
//! \tparam Val  - value type
//! \tparam CompareF  - NONSTANDARD comparison function:
//! 		int cmp(const IT::reference iterVal, const T val);
//! 			returns negative if iterVal < val, 0 if iterVal == val
//! 			and positive otherwise
//! \return iterator of the value position or end
template <typename IT, typename Val, typename CompareF>
enable_if_t<sizeof(Val) <= sizeof(void*), IT>
linear_ifind(IT begin, const IT end, const Val val
	, CompareF cmp=bsVal<MemberValCRef<IT>>) noexcept
{
	static_assert(sizeof(Val) <= sizeof(void*)
		, "linear_ifind(), Val should not exceed the address type size");
	static_assert(is_same<decltype(cmp(*begin, val)), decltype(bsVal(nullptr, nullptr))>::value
		, "linear_ifind(), cmp() must return the same type as bsVal()");
	static_assert(is_iterator<IT>(), "linear_ifind(), IT must be a forward iterator type");
#if VALIDATE >= 2
	assert(distance(begin, end) <= BINSEARCH_MARGIN // * 2
		&& "linear_ifind(), a small number of items is expected");
#endif // VALIDATE
	while(begin != end && cmp(*begin, val) < 0)
		++begin;
	return begin;
}

template <typename IT, typename Val, typename CompareF>
enable_if_t<(sizeof(Val) > sizeof(void*)), IT>
linear_ifind(IT begin, const IT end, const Val& val
	, CompareF cmp=bsVal<MemberValCRef<IT>>) noexcept
{
	static_assert(sizeof(Val) > sizeof(void*)
		, "linear_ifind() obj, compound type Val is expected");
	static_assert(is_same<decltype(cmp(*begin, val)), decltype(bsVal(nullptr, nullptr))>::value
		, "linear_ifind() obj, cmp() must return the same type as bsVal()");
	static_assert(is_iterator<IT>(), "linear_ifind(), IT must be a forward iterator type");
#if VALIDATE >= 2
	assert(distance(begin, end) <= BINSEARCH_MARGIN // * 2
		&& "linear_ifind() obj, a small number of items is expected");
#endif // VALIDATE
	while(begin != end && cmp(*begin, val) < 0)
		++begin;
	return begin;
}

//! Linear search in the specified range
// Note: inline function anyway does not pass any arguments, so there is no
// sense to use const ref args
template <typename IT, typename Val, typename CompareF>
inline IT linear_find(IT begin, const IT end, const Val val
	, CompareF cmp=bsVal<MemberValCRef<IT>>) noexcept
{
	static_assert(sizeof(Val) <= sizeof(void*)
		, "linear_find(), Val should not exceed the address type size");
	begin = linear_ifind(begin, end, val, cmp);
	return begin != end && !cmp(*begin, val) ? begin : end;
}

//! Linear index find in the sorted dataset
//! \attention Te returned iterator is not necessary point to the end if the item is not present
// Note: Container is passed by non-const ref for better cache reuse to
// not employ const iterators not used by the nearby code
template <typename Container, typename Val, typename CompareF>
inline conditional_t<is_const<Container>::value, typename Container::const_iterator, typename Container::iterator>
linear_ifind(Container& cnt, const Val val
	, CompareF cmp=bsVal<MemberValCRef<Container>>) noexcept
{
	return linear_ifind(cnt.begin(), cnt.end(), val, cmp);
}

//! Fast index find in the sorted dataset using either binary or linear search
//! \attention Te returned iterator is not necessary point to the end if the item is not present
template <typename RandIT, typename Val, typename CompareF>
inline enable_if_t<is_iterator<RandIT, std::random_access_iterator_tag>(), RandIT>
fast_ifind(RandIT begin, const RandIT end, const Val val
	, CompareF cmp=bsVal<MemberValCRef<RandIT>>) noexcept
{
	static_assert(sizeof(Val) <= sizeof(void*)
		, "fast_ifind(), Val should not exceed the address type size");
	return (end - begin) < BINSEARCH_MARGIN
		? linear_ifind(begin, end, val, cmp)
		: binary_ifind(begin, end, val, cmp);
}

//! Fast index find in the sorted dataset using either binary or linear search
// Note: Container is passed by non-const ref for better cache reuse to
// not employ const iterators not used by the nearby code
template <typename Container, typename Val, typename CompareF>
inline enable_if_t<is_iterator<typename Container::iterator, std::random_access_iterator_tag>()
	, conditional_t<is_const<Container>::value, typename Container::const_iterator, typename Container::iterator>>
fast_ifind(Container& cnt, const Val val
	, CompareF cmp=bsVal<MemberValCRef<Container>>) noexcept
{
	return fast_ifind(cnt.begin(), cnt.end(), val, cmp);
}

//! Fast find in the sorted dataset using either binary or linear search
template <typename RandIT, typename Val, typename CompareF>
inline enable_if_t<is_iterator<RandIT, std::random_access_iterator_tag>(), RandIT>
fast_find(RandIT begin, const RandIT end, const Val val
	, CompareF cmp=bsVal<MemberValCRef<RandIT>>) noexcept
{
	static_assert(sizeof(Val) <= sizeof(void*)
		, "fast_find(), Val should not exceed the address type size");
	return (end - begin) < BINSEARCH_MARGIN
		? linear_find(begin, end, val, cmp)
		: binary_find(begin, end, val, cmp);
}

//! Fast find in the sorted dataset using either binary or linear search
// Note: Container is passed by non-const ref for better cache reuse to
// not employ const iterators not used by the nearby code
template <typename Container, typename Val, typename CompareF>
inline enable_if_t<is_iterator<typename Container::iterator, std::random_access_iterator_tag>()
	, conditional_t<is_const<Container>::value, typename Container::const_iterator, typename Container::iterator>>
fast_find(Container& cnt, const Val val
	, CompareF cmp=bsVal<MemberValCRef<Container>>) noexcept
{
	static_assert(is_iterator<typename Container::iterator, std::random_access_iterator_tag>()
		, "fast_find(), Container iterator must be a random iterator type");
	return fast_find(cnt.begin(), cnt.end(), val, cmp);
}

// Operations with containers -------------------------------------------------
////! \brief Extends ordered nonempty elements
////! \pre els should be not empty and el should not be already present in els
////! \post iel points to the element FOLLOWING the inserted one after the execution
////!
////! \param els ContainerT&  - elements to be extended considering the ordering
//////! \param iel typename ContainerT::iterator  - begin iterator in the els to
//////! 	speedup incremental inserts.
//////! 	Points to the elements FOLLOWING the inserted one after the execution.
////! \param el typename ContainerT::value_type  - element to be added to els
////! \param cmp=bsVal CompareF  - comparison function
////! 		int cmp(const ContainerT::value_type& iterVal
////!				, const ContainerT::value_type& val);
////!
////! 			returns negative if iterVal < val, 0 if iterVal == val
////! 			and positive otherwise
////! \return void
//template <typename ContainerT, typename CompareF>
//void addSorted(ContainerT& els  //, typename ContainerT::iterator iel
//	, const typename ContainerT::value_type el
//	, CompareF cmp=bsVal<MemberValCRef<ContainerT>>)
//{
//	static_assert(sizeof(typename ContainerT::value_type) <= sizeof(void*)
//		, "addSorted(), value_type should not exceed the address type size");
//	assert(!els.empty() && "addSorted(), els should be not empty");
//
//	auto iel = els.begin();
//	iel = cmp(*iel, el) < 0 ? fast_ifind(++iel, els.end(), el, cmp) : iel;
//	assert((iel == els.end() || cmp(*iel, el) != 0) && "addSorted(), elements must be unique");
//	els.insert(iel, el);
////	iel = els.insert(iel, el);
////	++iel;  // Note: required to speedup incremental inserts
//}

//! \brief Check whether the items are sorted (ordered by the cmp func)
//!
//! \tparam IterT  - Iterator type
//!
//! \param begin IterT  - begin iterator of the checking range
//! \param end const IterT  - end iterator of the checking range
//! \param cmp=bsVal CompareF  - comparison function
//! \param unique=true bool  - all items are unique or allow duplicates
//! \return bool  - the items are ordered by cmp func
template <typename IterT, typename CompareF>
bool sorted(IterT begin, const IterT end
	, CompareF cmp=bsVal<MemberValCRef<IterT>>, bool unique=true) noexcept // noexcept(CompareF())
{
	static_assert(is_same<decltype(cmp(*begin, *begin)), decltype(bsVal(nullptr, nullptr))>::value
		, "sorted(), cmp() must return the same type as bsVal()");
	static_assert(is_iterator<IterT>(), "sorted(), IterT must be a forward iterator type");
#if VALIDATE >= 2
	assert(begin != end && "sorted(), empty container is being checked");
#endif // VALIDATE
	if(begin == end)
		return true;
	// The container is sorted if it empty or has only one item
	while(++begin != end) {
		auto res = cmp(*begin, *(begin - 1));
		if(res < 0 || (res == 0 && unique))
			return false;
	}
	return true;
}

//! \brief Index of the specified element in the sorted elements
//! \pre Elements are ordered and unique
//! \note If not INSORTED_NONUNIQUE and VALIDATE >= 1 then the elements are validated
//! 	to not contain the searching el
//!
//! \tparam ContainerT  - container type
//! \tparam ItemT  - anchor type
//! \tparam CompareF  - type of the binary search comparison function
//!
//! \param els ContainerT&  - elements to search at
//! \param el const ItemT  - anchor whose position is being defined
//! \param cmp=bsVal<ItemT> CompareF  - comparison function
//! \return typename ContainerT::iterator  - iterator to els where dest should be inserted/emplaced
template <typename ContainerT, typename ItemT, typename CompareF>
typename ContainerT::iterator insorted(ContainerT& els, const ItemT el, CompareF cmp=bsVal<ItemT>)
#if defined(INSORTED_NONUNIQUE)
	noexcept(CompareF())
#endif // INSORTED_NONUNIQUE
{
	static_assert(sizeof(ItemT) <= sizeof(void*)
		, "insorted(), ItemT should not exceed the address type size");
#if VALIDATE >= 3
	assert(sorted(els.begin(), els.end(), cmp) && "insorted(), els should be sorted");
#endif // VALIDATE
	auto iel = els.begin();
	if(!els.empty() && cmp(*iel, el) < 0)
		iel = fast_ifind(++iel, els.end(), el, cmp);
#if TRACE >= 2
	if(!(iel == els.end() || cmp(*iel, el) != 0))
		fprintf(ftrace, "  >>>>>>> insorted(), #%u (%p) is already among %lu members\n"
			, el->id, el, els.size());
#endif // TRACE
#ifndef INSORTED_NONUNIQUE
	//assert((iel == els.end() || cmp(*iel, el) != 0) && "insorted(), elements must be unique");
	if(iel != els.end() && !cmp(*iel, el))
		throw invalid_argument(string("insorted(), elements must be unique: #")
			.append(std::to_string(derefItemX(el).id)) += '\n');
#endif // INSORTED_NONUNIQUE
	return iel;
}

//! \brief Index of the specified element in the SMALL number of sorted elements
//! \pre Elements are ordered, unique and does not contain the searching el.
//! 	The number of elements is small (typically up to BINSEARCH_MARGIN)
//!
//! \tparam ContainerT  - container type
//! \tparam ItemT  - anchor type
//! \tparam CompareF  - type of the binary search comparison function
//!
//! \param els ContainerT&  - elements to search at
//! \param el const ItemT  - anchor whose position is being defined
//! \param cmp=bsVal<ItemT> CompareF  - comparison function
//! \return typename ContainerT::iterator  - iterator to els where dest should be inserted/emplaced
template <typename ContainerT, typename ItemT, typename CompareF>
inline typename ContainerT::iterator insortedLight(ContainerT& els, const ItemT el
	, CompareF cmp=bsVal<ItemT>)
#if defined(INSORTED_NONUNIQUE)
	noexcept(CompareF())
#endif // INSORTED_NONUNIQUE
{
	static_assert(sizeof(ItemT) <= sizeof(void*)
		, "insortedLight(), ItemT should not exceed the address type size");
#if VALIDATE >= 2
	assert(els.size() <= BINSEARCH_MARGIN // * 2
		&& "insortedLight(), a small number of items is expected");
#endif // VALIDATE
//#if VALIDATE >= 3
//	assert(sorted(els.begin(), els.end(), cmp) && "insorted(), els should be not sorted");  // Note: issues occur in case of cmp is bsObjsDest
//#endif // VALIDATE
	auto iel = els.begin();
	if(!els.empty() && cmp(*iel, el) < 0)
		iel = linear_ifind(++iel, els.end(), el, cmp);
#if TRACE >= 2
	if(!(iel == els.end() || cmp(*iel, el) != 0))
		fprintf(ftrace, "  >>>>>>> insortedLight(), #%u (%p) is already among %lu members\n"
			, el->id, el, els.size());
#endif // TRACE
#ifndef INSORTED_NONUNIQUE
	//assert((iel == els.end() || cmp(*iel, el) != 0) && "insortedLight(), elements must be unique");
	if(iel != els.end() && !cmp(*iel, el))
		throw invalid_argument(string("insortedLight(), elements must be unique: #")
			.append(std::to_string(derefItemX(el).id)) += '\n');
#endif // INSORTED_NONUNIQUE
	return iel;
}

//! \brief Find min value among the specified arguments of the same type
//!
//! \param val ValT  - first argument
//! \param vs Ts...  - remained arguments
//! \return ValtT  - minimal argument
template <typename ValT, typename ...Ts>
ValT minval(ValT val, Ts... vs) noexcept
{
	static_assert(sizeof(ValT) <= sizeof(void*), "minval(), input types validation failed");
	return min(val, minval(vs...));  // Note: Tail call optimization is expected here instead of the recursion
}

//template <typename ValT, typename ...Ts>
//enable_if_t<!(sizeof(ValT) <= sizeof(void*)), ValT&> minval(ValT& val, Ts... vs) noexcept
//{
//	return min(val, minval(vs...));
//}

}  // daoc

#endif // OPERATIONS_H
