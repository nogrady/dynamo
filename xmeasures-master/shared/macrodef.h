//! \brief Global macro definitions.
//! The Dao (Deterministic Agglomerative Overlapping) of Clustering library:
//! Robust & Fine-grained Deterministic Clustering for Large Networks.
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2016-07-25

#ifndef MACRODEF_H
#define MACRODEF_H

// Global MACROSES:
//	- VALIDATE  - use alternative evaluations to validate results
//		- 0  - turn off heavy validation
//		- 1  - default value for the heavy validation
//		- 2  - extra heavy validation (might duplicate already performed heavy validation)
//		- 3  - cross validation of functions (executed on each call, but only once is enough)
//
//	- TRACE, TRACE_EXTRA  - detailed tracing under debug (trace nodes weights)
//		- 0  - turn off the tracing
//		- 1  - brief tracing that can be used in release to show warnings, etc.
//		- 2  - detailed tracing for DEBUG
//		- 3  - extra detailed tracing
//
//	- FTRACE_GLOBAL  - use global ftrace file for the whole project, or "shared/" headers
//		define it locally
//
//	- UTEST  - build [also] unit tests, requires installation and linking of the unit test library.
//
// NOTE: undefined maro definition is interpreted as having value 0

#ifndef TRACE
#ifdef DEBUG
	#define TRACE 2
#elif !defined(NDEBUG)  // RELEASE, !NDEBUG
	#define TRACE 1
//#else  // RELEASE, NDEBUG
//	#define TRACE 0
#endif // DEBUG
#endif // TRACE

#ifndef VALIDATE
#ifdef DEBUG
	#define VALIDATE 2
#elif !defined(NDEBUG)  // RELEASE, !NDEBUG
	#define VALIDATE 1
//#else  // ELEASE, NDEBUG
//	#define VALIDATE 0
#endif // DEBUG
#endif // VALIDATE

// SWIG related macro definitions
// Swig 3.0.12 does not understand some structures, workarounds are applied
// Note: defined only for SWIG interfaces
#ifdef SWIG
	// Just skip the static assert
    #define static_assert(a, b)
#endif // SWIG

// Note: SWIG_VERSION is not defined for SWIGJAVA and SWIGCSHARP
// Note: defined both for the SWIG interfaces and implementation
#if defined(SWIG_VERSION) || defined(SWIGJAVA) || defined(SWIGCSHARP)
	// Defined automatically when any SWIG processing is performed
	// (either the included as SWIG interface or implementation)
	#define DAOC_SWIGPROC
#endif // SWIG processing

// Define macros for the case when SWIG supports functions overloading
#if defined(SWIGCSHARP) || defined(SWIGD) || defined(SWIGJAVA)
	#define SWIG_OVERLOADS
#endif // OVERLOADS


#endif // MACRODEF_H
