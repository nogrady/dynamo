//! \brief Extrinsic measures evaluation interface implementation.
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2017-12-15

#include <cstdio>
//#include <bitset>
#include <errno.h>

#include "operations.hpp"
#include "interface.h"


using std::overflow_error;
using std::invalid_argument;
using namespace daoc;


// Omega Index related types and functions -------------------------------------
Id mutualnum(const RawClusterPtrs* a, const RawClusterPtrs* b, const Id nmax) noexcept
{
#if VALIDATE >= 2
	assert(a && b && "mutualnum(), valid containers are expected");
#endif // VALIDATE
	Id num = 0;
	if(b->size() < a->size()) {
		auto t = b;
		b = a;
		a = t;
	}
	const auto  eb = b->end();
	auto  ib = b->begin();
	for(auto acp: *a) {
		while(ib != eb && cmpBase(*ib, acp))
			++ib;
		if(ib == eb
		|| (*ib == acp && ++num >= nmax))
			break;
	}
	return num;
}

Id mutualnum(const RawClusterPtrs* a, const RawClusterPtrs* b) noexcept
{
#if VALIDATE >= 2
	assert(a && b && "mutualnum(), valid containers are expected");
#endif // VALIDATE
	Id num = 0;
	if(b->size() < a->size()) {
		auto t = b;
		b = a;
		a = t;
	}
	const auto  eb = b->end();
	auto  ib = b->begin();
	for(auto acp: *a) {
		while(ib != eb && cmpBase(*ib, acp))
			++ib;
		if(ib == eb)
			break;
		if(*ib == acp)
			++num;
	}
	return num;
}

// Other Measures related functions --------------------------------------------
//string to_string(Evaluation eval, bool bitstr)
//{
//	static_assert(sizeof(Evaluation) == sizeof(EvalBase)
//		, "to_string(), Evaluation type must be the same size as EvalBase");
//	// Convert to bit string
//	if(bitstr)
//		return bitset<sizeof(Evaluation) * 8>(static_cast<EvalBase>(eval))
//			.to_string().insert(0, "0b");
//
//	// Convert to semantic string
//	string  val;
//	switch(eval) {
//	case Evaluation::MULTIRES:
//		val = "MULTIRES";
//		break;
//	case Evaluation::OVERLAPPING:
//		val = "OVERLAPPING";
//		break;
//	case Evaluation::MULRES_OVP:
//		val = "MULRES_OVP";
//		break;
//	case Evaluation::NONE:
//	default:
//		val = "NONE";
//	}
//	return val;
//}

string to_string(F1 f1)
{
	// Convert to semantic string
	string  val;
	switch(f1) {
	case F1::PARTPROB:
		val = "PARTPROB";
		break;
	case F1::HARMONIC:
		val = "HARMONIC";
		break;
	case F1::AVERAGE:
		val = "AVERAGE";  // Suggested by Leskovec
		break;
	case F1::NONE:
	default:
		val = "NONE";
	}
	return val;
}

string to_string(Match mkind)
{
	// Convert to semantic string
	string  val;
	switch(mkind) {
	case Match::WEIGHTED:
		val = "WEIGHTED";
		break;
	case Match::UNWEIGHTED:
		val = "UNWEIGHTED";
		break;
	case Match::COMBINED:
		val = "COMBINED";
		break;
	case Match::NONE:
	default:
		val = "NONE";
	}
	return val;
}

bool xwmatch(Match m) noexcept
{
	return m == Match::WEIGHTED || m == Match::COMBINED;
}


bool xumatch(Match m) noexcept
{
	return m == Match::UNWEIGHTED || m == Match::COMBINED;
}

NodeBase NodeBase::load(const char* filename, float membership
	, ::AggHash* ahash, size_t cmin, size_t cmax, bool verbose)
{
	NodeBase  nb;  // Return using NRVO optimization
	NamedFileWrapper  finp(filename, "r");
	if(finp)
		static_cast<UniqIds&>(nb) = loadNodes<Id, AccId>(finp, membership
			, ahash, cmin, cmax, verbose);
	else perror((string("WARNING load(), can't open ") += filename).c_str());

	return nb;
}

// Accessory functions ---------------------------------------------------------
Id  parseId(char* str)
{
#if VALIDATE >= 2
	assert(!errno && "Initial errno should be zero");
#endif // VALIDATE
	auto nid = strtoul(str, nullptr, 10);
	static_assert(sizeof(nid) >= sizeof(Id), "Parsing value type is too small for Id");
	if(nid > numeric_limits<Id>::max() || (!nid && errno != 0)) {
		if(nid > numeric_limits<Id>::max())
			throw overflow_error("Loaded value of id is too large: " + std::to_string(nid) + "\n");
		else if(errno != 0)
			throw invalid_argument(string("Conversion to id can't be performed: ").append(str)
				+ ", errno: " + std::to_string(errno).append("\n"));
	}
	return nid;
}

AccProb hmean(AccProb a, AccProb b) noexcept
{
	static_assert(is_floating_point<AccProb>::value, "AccProb should be a floating point type");
	// Note: both a = b = 0 and a = -b are considered and yield 0
	return a + b != 0 ? 2 * a / (a + b) * b : 0;
}

AccProb gmean(AccProb a, AccProb b) noexcept
{
#ifdef DEBUG
	assert(a >= 0 && b >= 0 && "gmean(), the probabilities should E [0, 1]");
#endif // DEBUG
	return sqrt(a * b);
}

AccProb amean(AccProb a, AccProb b) noexcept
{
	return (a + b) / 2;
}
