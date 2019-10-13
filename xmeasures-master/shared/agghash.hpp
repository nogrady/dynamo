//! \brief AggHash simple (Aggregating Order Invariant Hashing) of the DAOC clustering library
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2017-02-21

#ifndef CODING_HPP
#define CODING_HPP

#include <cstdint>  // uintX_t
//#include <cstddef>  // size_t
#include <string>  // uintX_t
#include <functional>  // hash
//#include <cstring>  // memcmp
#include <type_traits>  // is_integral
#include <limits>  // numeric_limits
#include <stdexcept> // numeric_limits


namespace daoc {

using std::string;
using std::is_integral;
using std::numeric_limits;
using std::domain_error;

// Type Declarations ---------------------------------------------------
//! \brief Aggregation hash of ids
//! \pre Template types should be integral
//!
//! \tparam Id  - type of the member ids
//! \tparam AccId  - type of the accumulated Ids and accumulated squares of Ids
//! should have at least twice magnitude of the Id type (i.e. squared)
template <typename Id=uint32_t, typename AccId=uint64_t>
class AggHash {
	static_assert(is_integral<Id>::value && is_integral<AccId>::value
		&& sizeof(AccId) >= 2*sizeof(Id), "AggHash, types constraints are violated");

	// ATTENTION: type of the m_size should not be less than of m_idsum to
	// avoid gaps filled with trash on memory alignment
	// Note: size should be first as the most discriminative attribute, which
	// can be potentially used for the ordering
	// Note: the size is redundant and does not have any impact except for the structured ordering
	// if the sum does not increase AccId_MAX or if zero value of id is NOT allowed. The size is
	// necessary if id=0 may be present in the clusters.
	AccId  m_size;  //!< Size of the container
	AccId  m_idsum;  //!< Sum of the member ids
	AccId  m_id2sum;  //!< Sum of the squared member ids
protected:
	//! Id correction to prevent collisions
	constexpr static Id  idcor = sqrt(numeric_limits<Id>::max());
public:
	// Export the template parameter types
	using IdT = Id;  //!< Type of the member ids
	using AccIdT = AccId;  //!< Type of the accumulated Ids and accumulated squares of Ids

	//! \brief Default constructor
	AggHash() noexcept
	: m_size(0), m_idsum(0), m_id2sum(0) {}

	//! \brief Add id to the aggregation
	//! \note In case correction is used and id becomes out of range (initial id > IDMAX - IDCORR)
	//! 	then an exception is thrown, which crashes the whole application, which is OK
	//!
	//! \param id Id  - id to be included into the hash
	//! \return void
	void add(Id id) noexcept;

	//! \brief Clear/reset the aggregation
	//!
	//! \return void
	void clear() noexcept;

	//! \brief Number of the aggregated ids
	//!
	//! \return size_t  - number of the aggregated ids
	size_t size() const noexcept  { return m_size; }

	//! \brief Sum of the aggregated ids
	//!
	//! \return size_t  - sum of the aggregated ids
	size_t idsum() const noexcept  { return m_idsum; }

	//! \brief Sum of squares of the aggregated ids
	//!
	//! \return size_t  - sum of squares of the aggregated ids
	size_t id2sum() const noexcept  { return m_id2sum; }

//    //! \brief The hash is empty
//    //!
//    //! \return bool  - the hash is empty
//	bool empty() const noexcept  { return !m_size; }

	//! \brief Evaluate hash of the aggregation
	//!
	//! \return size_t  - resulting hash
	size_t hash() const;

	//! \brief Operator less
	//!
	//! \param ah const AggHash&  - comparing object
	//! \return bool operator  - result of the comparison
	inline bool operator <(const AggHash& ah) const noexcept;

	//! \brief Operator less or equal
	//!
	//! \param ah const AggHash&  - comparing object
	//! \return bool operator  - result of the comparison
	inline bool operator <=(const AggHash& ah) const noexcept;

	//! \brief Operator greater
	//!
	//! \param ah const AggHash&  - comparing object
	//! \return bool operator  - result of the comparison
	bool operator >(const AggHash& ah) const noexcept  { return !(*this <= ah); }

	//! \brief Operator greater or equal
	//!
	//! \param ah const AggHash&  - comparing object
	//! \return bool operator  - result of the comparison
	bool operator >=(const AggHash& ah) const noexcept  { return !(*this < ah); }

	//! \brief Operator equal
	//!
	//! \param ah const AggHash&  - comparing object
	//! \return bool operator  - result of the comparison
	inline bool operator ==(const AggHash& ah) const noexcept;

	//! \brief Operator unequal (not equal)
	//!
	//! \param ah const AggHash&  - comparing object
	//! \return bool operator  - result of the comparison
	bool operator !=(const AggHash& ah) const noexcept  { return !(*this == ah); }
};

// Type Definitions ----------------------------------------------------
template <typename Id, typename AccId>
void AggHash<Id, AccId>::add(Id id) noexcept
{
	id += idcor;  // Correct id to prevent collisions (see AgordiHash for details)
	// Check for the overflow after the correction
    // Note: the exception will crash the whole app since noexcept is used but it is fine
	if(id < idcor)
		throw domain_error(string("The corrected value of ").append(std::to_string(id))
			.append(" is too large and causes the overflow\n"));
	++m_size;
	m_idsum += id;
	m_id2sum += id * id;
}

template <typename Id, typename AccId>
void AggHash<Id, AccId>::clear() noexcept
{
	m_size = 0;
	m_idsum = 0;
	m_id2sum = 0;
}

template <typename Id, typename AccId>
size_t AggHash<Id, AccId>::hash() const
{
	// ATTENTION: requires filling with zero memory alignment trash or avoid the padding
	return std::hash<string>()(string(reinterpret_cast<const char*>(this), sizeof *this));
}

template <typename Id, typename AccId>
bool AggHash<Id, AccId>::operator <(const AggHash& ah) const noexcept
{
	return m_size < ah.m_size || (m_size == ah.m_size
		&& (m_idsum < ah.m_idsum || (m_idsum == ah.m_idsum && m_id2sum < ah.m_id2sum)));
}

template <typename Id, typename AccId>
bool AggHash<Id, AccId>::operator <=(const AggHash& ah) const noexcept
{
	return m_size < ah.m_size || (m_size == ah.m_size
		&& (m_idsum < ah.m_idsum || (m_idsum == ah.m_idsum && m_id2sum <= ah.m_id2sum)));
}

template <typename Id, typename AccId>
bool AggHash<Id, AccId>::operator ==(const AggHash& ah) const noexcept
{
	return m_size == ah.m_size && m_idsum == ah.m_idsum && m_id2sum == ah.m_id2sum;
	//return !memcmp(this, &ah, sizeof(AggHash));  // Note: memcmp returns 0 on full match
}

}  // daoc

#endif // CODING_HPP
