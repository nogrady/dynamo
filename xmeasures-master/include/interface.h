//! \brief Extrinsic measures evaluation interface.
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2017-02-13

#ifndef INTERFACE_H
#define INTERFACE_H

#include <unordered_map>
#include <memory>  // unique_ptr
#include <string>
#include <type_traits>
#include <limits>
#if VALIDATE >= 1
#include <stdexcept>
#endif // VALIDATE

#define INCLUDE_STL_FS
#include "fileio.hpp"
#if VALIDATE >= 2
#include "operations.hpp"
#endif // VALIDATE 2


using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::unique_ptr;
using std::string;
using std::pair;
using std::is_integral;
using std::is_pointer;
using std::is_floating_point;
using std::is_arithmetic;
using std::is_same;
//using std::enable_if;
using std::enable_if_t;
using std::conditional_t;
using std::numeric_limits;
#if VALIDATE >= 2
using std::invalid_argument;
#endif // VALIDATE

// Data Types ------------------------------------------------------------------
using Id = uint32_t;  //!< Node id type
// Note: Size should a magnitude larger than Id to hold Id*Id
using AccId = uint64_t;  //!< Accumulated Id type

using Prob = float;  //!< Probability
using AccProb = double;  //!< Accumulated Probability

//! Aggregated Hash of the loading cluster member ids
using AggHash = daoc::AggHash<Id, AccId>;

using RawIds = vector<Id>;  //!< Node ids, unordered

// Omega Index related types and functions -------------------------------------
using RawCluster = RawIds;  //!< Raw cluster of member node ids
using RawClusters = vector<RawCluster>;  //!< Raw clustering, container of the raw clusters
using RawClusterPtrs = vector<RawCluster*>;
using NodeRClusters = unordered_map<Id, pair<RawClusterPtrs, RawClusterPtrs>>;  //!< Raw node membership in the clusters

//! \brief Omega Index evaluation
//!
//! \tparam EXT bool  - extended Omega Index, which does not excessively penalize
//! 	distinct node shares
//!
//! \param ndrcs const NodeRClusters&  - node raw clusters relations
//! \param cls1 const RawClusters&  - clusters of the first collection
//! \param cls2 const RawClusters&  - clusters of the second collection
//! \return Prob  - omega index
template <bool EXT=false>
Prob omega(const NodeRClusters& ndrcs, const RawClusters& cls1, const RawClusters& cls2);

//! \brief Evaluate the number of mutual raw cluster pointers in the containers
//!
//! \pre Input raw clusters pointer containers are ordered by the cmpBase<RawCluster*>
//!
//! \param a const RawClusterPtrs*  - first raw cluster pointers
//! \param b const RawClusterPtrs*  - second raw cluster pointers
//! \param nmax const Id  - max number of matches for the early termination,
//! 	0 is allowed but senseless.
//! \return Id  - the number of mutual members
Id mutualnum(const RawClusterPtrs* a, const RawClusterPtrs* b, const Id nmax) noexcept;

Id mutualnum(const RawClusterPtrs* a, const RawClusterPtrs* b) noexcept;

// F1 & NMI related data types -------------------------------------------------
template <typename Count>
struct Cluster;

//! Cluster matching counter
//! \note Required only for F1 evaluation
//! \tparam Count  - arithmetic counting type
template <typename Count>
class Counter {
public:
	static_assert(is_arithmetic<Count>::value
		, "Counter(), Count should be an arithmetic type");
	using CountT = Count;  //!< Count type, arithmetic
	using ClusterT = Cluster<Count>;
private:
	// Note: it's OK to copy this pointer on assignment since it is never
	// allocated in this class
	ClusterT*  m_orig;  //!<  Originator cluster
	CountT  m_count;  //!<  Occurrences counter, <= members size
public:
    //! Default constructor
	Counter(): m_orig(nullptr), m_count(0)  {}

    //! \brief Update the counter from the specified origin
    //!
    //! \param orig ClusterT*  - counter origin
    //! \param cont Count  - contribution or share, actual only for the floating point counter
    //! \return void
	void operator()(ClusterT* orig, Count cont)
#if VALIDATE < 2
	noexcept
#endif // VALIDATE
	{
		if(m_orig != orig) {
			m_orig = orig;
			m_count = 0;
		}
		if(is_integral<CountT>::value)
			++m_count;
		else {
			static_assert(!is_floating_point<CountT>::value || sizeof(m_count) >= sizeof(double)
				, "operator(), types validation failed");
#if VALIDATE >= 2
			if(cont <= 0 || cont > 1)
				throw invalid_argument("operator(), cont should E (0, 1]\n");
#endif // VALIDATE
			m_count += cont;
		}
	}

    //! \brief Get counted value
    //!
    //! \return CountT  - counted value
	CountT operator()() const noexcept  { return m_count; }

    //! \brief Get counter origin
    //!
    //! \return ClusterT*  - counter origin
	ClusterT* origin() const noexcept  { return m_orig; }

    //! \brief Clear (reset) the counter
	void clear() noexcept
	{
		m_orig = nullptr;
		m_count = 0;
	}
};

//! Cluster
//! \tparam Count  - nodes contribution counter type
template <typename Count>
struct Cluster {
	static_assert(is_arithmetic<Count>::value
		, "Counter(), Count should be an arithmetic type");
	using CountT = Count;  //!< Count type, arithmetic

	RawIds  members;  //!< Node ids, unordered
	// Note: used by F1 only and always
	Counter<Count>  counter;  //!< Cluster matching counter
	////! Accumulated contribution
	//using AccCont = conditional_t<m_overlaps, Count, AccId>;
	//!< Contribution from members
	// Note: used only in case of a) overlaps (by all measures) or
	// b) multiple resolutions (by NMI only)
	Count  mbscont;
	static_assert(!is_floating_point<Count>::value || sizeof(mbscont) >= sizeof(double)
		, "operator(), types validation failed");

    //! Default constructor
	Cluster();

    //! \brief F1 measure
    //! \pre Clusters should be valid, i.e. non-empty
    //!
    //! \param matches Count  - the number of matched members
    //! \param capacity Count  - contributions capacity of the matching foreign cluster
    //! \return AccProb  - resulting F1
	AccProb f1(Count matches, Count capacity) const
#if VALIDATE < 2
	noexcept
#endif // VALIDATE
	{
		// F1 = 2 * pr * rc / (pr + rc)
		// pr = m / c1
		// rc = m / c2
		// F1 = 2 * m/c1 * m/c2 / (m/c1 + m/c2) = 2 * m / (c2 + c1)
		// ATTENTION: F1 compares clusters per-pair, so it is much simpler and has another
		// semantics of contribution for the multi-resolution case
		const Count  contrib = is_floating_point<Count>::value ? cont() : members.size();
#if VALIDATE >= 2
		if(matches < 0 || daoc::less<conditional_t<is_floating_point<Count>::value
		, Prob, Count>>(capacity, matches) || contrib <= 0)
			throw invalid_argument(string("f1(), both clusters should be non-empty, matches: ")
				.append(std::to_string(matches)).append(", capacity: ").append(std::to_string(capacity))
				.append(", contrib: ").append(std::to_string(contrib)) += '\n');
#endif // VALIDATE
		return 2 * matches / AccProb(capacity + contrib);  // E [0, 1]
		// Note that partial probability (non-normalized to the remained matches,
		// it says only how far this match from the full match) of the match is:
		// P = AccProb(matches * matches) / AccProb(size * members.size()),
		// where nodes contribution instead of the size should be use for overlaps.
		// The probability is more discriminative than F1 for high values
	}

    //! \brief Partial probability of the match (non-normalized to the other matches)
    //! \pre Clusters should be valid, i.e. non-empty
    //!
    //! \param matches Count  - the number of matched members
    //! \param capacity Count  - contributions capacity of the matching foreign cluster
    //! \return AccProb  - resulting probability
	AccProb pprob(Count matches, Count capacity) const
#if VALIDATE < 2
	noexcept
#endif // VALIDATE
	{
		// P = P1 * P2 = m/c1 * m/c2 = m*m / (c1*c2),
		// where nodes contribution instead of the size should be used for overlaps.
		// ATTENTION: F1 compares clusters per-pair, so it is much simpler and has another
		// semantics of contribution for the multi-resolution case comparing to NMI
		// that also uses cont()
		constexpr bool  floating = is_floating_point<Count>::value;
		const Count  contrib = floating ? cont() : members.size();
#if VALIDATE >= 2
		if(matches < 0 || daoc::less<conditional_t<floating, Prob, Count>>
		(capacity, matches) || contrib <= 0)
			throw invalid_argument(string("pprob(), both clusters should be non-empty, matches: ")
				.append(std::to_string(matches)).append(", capacity: ").append(std::to_string(capacity))
				.append(", contrib: ").append(std::to_string(contrib)) += '\n');
#endif // VALIDATE
		return floating ? static_cast<AccProb>(matches) * matches / (static_cast<AccProb>(capacity) * contrib)
			: static_cast<AccProb>(static_cast<AccId>(matches) * matches)
				/ (static_cast<AccId>(capacity) * contrib);  // E [0, 1]
	}

    //! \brief Cluster members contribution
    //!
    //! \return Count  - total contribution from the members
	Count cont() const noexcept
	{
//		return is_same<decltype(mbscont), EmptyStub>::value ? members.size() : mbscont;
		return mbscont;
	}
};

//! Automatic storage for the Cluster;
//! \tparam Count  - arithmetic counting type
template <typename Count>
using ClusterHolder = unique_ptr<Cluster<Count>>;

//! Cluster pointers, unordered
//! \tparam Count  - arithmetic counting type
template <typename Count>
using ClusterPtrs = vector<Cluster<Count>*>;

//! Node to clusters relations
//! \tparam Count  - arithmetic counting type
template <typename Count>
using NodeClusters = unordered_map<Id, ClusterPtrs<Count>>;

//! Resulting greatest matches for 2 input collections of clusters in a single direction
using Probs = vector<Prob>;

// Label-related types --------------------------------------------------------
//! Clusters Labels, Labels are ORDERED by cmpBase
template <typename Count>
using ClustersLabels = unordered_map<Cluster<Count>*, ClusterPtrs<Count>>;

// F1-related types -----------------------------------------------------------
using F1Base = uint8_t;

//! \brief F1 kind
enum struct F1: F1Base {
	//! Not initialized
	NONE = 0,
	//! Harmonic mean of the [weighted] average of the greatest (maximal) match
	//! by partial probabilities
	PARTPROB,
	//! Harmonic mean of the [weighted] average of the greatest (maximal) match by F1s
	HARMONIC,
	//! Arithmetic mean (average) of the [weighted] average of the greatest (maximal)
	//! match by F1s, i.e. F1-Score
	AVERAGE  // Suggested by Leskovec
};

//! \brief String representation of the F1
//! \relates F1
//!
//! \param f1 F1  - the value to be converted
//! \return string  - string value
string to_string(F1 f1);

// NMI-related types -----------------------------------------------------------
//! Internal element of the Sparse Matrix with Vector Rows
//! \tparam Index  - index (of the column) in the row
//! \tparam Value  - value type
template <typename Index, typename Value>
struct RowVecItem {
	static_assert(is_integral<Index>::value || is_pointer<Index>::value
		, "RowVecItem, Index should be an integral type");

	using CallT = Index;  //!< Type of the functor call

	Index  pos;  //!< Position (index) in the row
	Value  val;  //!< Target value (payload)

	//! Constructor in case of the simple value
    //!
    //! \param i=Index() Index  - index of value in the row
    //! \param v=Value() Value  - payload value
	template <typename T=Value, enable_if_t<sizeof(T) <= sizeof(void*)>* = nullptr>
	RowVecItem(Index i=Index(), Value v=Value()) noexcept(Value())
	: pos(i), val(v)  {}

	//! Constructor in case of the compound value
    //!
    //! \param i=Index() Index  - index of value in the row
    //! \param v=Value() Value  - payload value
	template <typename T=Value, enable_if_t<(sizeof(T) > sizeof(void*)), bool>* = nullptr>
	RowVecItem(Index i=Index(), Value&& v=Value()) noexcept(Value())
	: pos(i), val(move(v))  {}

    //! \brief Functor (call) operator
    //!
    //! \return CallT  - index of the value
	// Note: required to call obj()
	CallT operator()() const noexcept  { return pos; }

//	// Note: required for the comparison operations with index
//	operator CallT() const noexcept  { return this }
};

//! Row vector for the SparseMatrix
template <typename Index, typename Value>
using SparseMatrixRowVec = vector<RowVecItem<Index, Value>>;

//! Base type of the SparseMatrix (can be unordered_map, map, vector)
template <typename Index, typename Value>
using SparseMatrixBase = unordered_map<Index, SparseMatrixRowVec<Index, Value>>;

//! Sparse Matrix
//! \tparam Index  - index type
//! \tparam Value  - value type
template <typename Index, typename Value>
struct SparseMatrix: SparseMatrixBase<Index, Value> {
	static_assert((is_integral<Index>::value || is_pointer<Index>::value)
		&& is_arithmetic<Value>::value, "SparseMatrix(), invalid parameter types");

	using IndexT = Index;  //!< Indexes type, integral
	using ValueT = Value;  //!< Value type, arithmetic
	using BaseT = SparseMatrixBase<IndexT, ValueT>;  //!< SparseMatrixBase type
	using RowT = typename BaseT::mapped_type;  //!< Matrix row type
	//! Matrix row element type, which contains the value and might have
	//! additional attributes
	using RowItemT = typename RowT::value_type;

    //! \brief Default constructor
    //!
    //! \param rows=0 Id  - initial number of rows
	SparseMatrix(Id rows=0);

    //! \brief Access specified element inserting it if not exists
    //!
    //! \param i Index  - row index
    //! \param j Index  - column index
    //! \return Value& operator  - value of the element to be set
	Value& operator ()(Index i, Index j);

    //! \brief Access specified element without bounds checking
    //! \note fast, but unsafe
    //!
    //! \param i Index  - row index
    //! \param j Index  - column index
    //! \return Value& operator  - value of the element
	template <typename T=Value, enable_if_t<sizeof(T) <= sizeof(void*)>* = nullptr>
	Value operator ()(Index i, Index j) const noexcept; //  { return this->at(i) }

    //! \brief Access specified element without bounds checking
    //! \note fast, but unsafe
    //!
    //! \param i Index  - row index
    //! \param j Index  - column index
    //! \return Value& operator  - value of the element
	template <typename T=Value, enable_if_t<(sizeof(T) > sizeof(void*)), bool>* = nullptr>
	const Value& operator ()(Index i, Index j) const noexcept; //  { return this->at(i) }

    //! \brief Access specified element checking the bounds
    //!
    //! \param i Index  - row index
    //! \param j Index  - column index
    //! \return Value& operator  - value of the element
	template <typename T=Value, enable_if_t<sizeof(T) <= sizeof(void*)>* = nullptr>
	Value at(Index i, Index j); //  { return this->at(i) }

    //! \brief Access specified element checking the bounds
    //!
    //! \param i Index  - row index
    //! \param j Index  - column index
    //! \return Value& operator  - value of the element
	template <typename T=Value, enable_if_t<(sizeof(T) > sizeof(void*)), bool>* = nullptr>
	const Value& at(Index i, Index j); //  { return this->at(i) }

	using BaseT::at;  //!< Provide direct access to the matrix row
};

//using EvalBase = uint8_t;  //!< Base type for the Evaluation
//
////! \brief Evaluation type
//enum struct Evaluation: EvalBase {
//	NONE = 0,
////	HARD = 0
//	MULTIRES = 1,  //!< Multi-resolution non-overlapping clusters, compatible with hard partitioning
//	OVERLAPPING = 2,  //!< Overlapping clusters, compatible with hard partitioning
//	MULRES_OVP = 3  //!< Multi-resolution clusters with possible overlaps on each resolution level
//};
//
////! \brief Convert Evaluation to string
////! \relates Evaluation
////!
////! \param flag Evaluation  - the flag to be converted
////! \param bitstr=false bool  - convert to bits string or to Evaluation captions
////! \return string  - resulting flag as a string
//string to_string(Evaluation eval, bool bitstr=false);

struct RawNmi {
	Prob  mi;  //!< Mutual information of two collections
	Prob  h1;  //!< Information content of the 1-st collection
	Prob  h2;  //!< Information content of the 2-nd collection
	//Evaluation  eval;  //!< Evaluation type

	static_assert(is_floating_point<Prob>::value, "RawNmi, Prob should be a floating point type");
	RawNmi() noexcept: mi(0), h1(numeric_limits<Prob>::quiet_NaN())
		, h2(numeric_limits<Prob>::quiet_NaN())  {}

	void operator() (Prob mutinf, Prob cn1h, Prob cn2h) noexcept
	{
		mi = mutinf;
		h1 = cn1h;
		h2 = cn2h;
	};
};

// Collection ------------------------------------------------------------------
//! Node base interface
struct NodeBaseI {
    //! \brief Default virtual destructor
	virtual ~NodeBaseI()=default;

    //! \brief Whether the node base is actual (non-empty)
    //!
    //! \return bool  - the node base is non-empty
	operator bool() const noexcept  { return ndsnum(); };

    //! \brief The number of nodes
    //!
    //! \return Id  - the number of nodes in the collection
	virtual Id ndsnum() const noexcept = 0;

    //! \brief Whether exists the specified node
    //!
    //! \param nid  - node id
    //! \return bool  - specified node id exists
	virtual bool nodeExists(Id nid) const noexcept = 0;
};

//! Unique ids (node ids)
using UniqIds = unordered_set<Id>;

//! Node base
struct NodeBase: protected UniqIds, NodeBaseI {
	using UniqIds::clear;

	//! \copydoc NodeBaseI::nodeExists(Id nid) const noexcept
	Id ndsnum() const noexcept  { return size(); }

	//! \copydoc NodeBaseI::nodeExists(Id nid) const noexcept
	bool nodeExists(Id nid) const noexcept  { return count(nid); }

	//! \brief Load all unique nodes from the CNL file with optional filtering by the cluster size
	//!
	//! \param filename const char*  - name of the input file
    //! \param ahash=nullptr AggHash*  - resulting aggregated hash of the loaded
    //! node ids if not nullptr
	//! \param membership=1 float  - expected membership of the nodes, >0, typically >= 1.
	//! Used only for the node container preallocation to estimate the number of nodes
	//! if not specified in the file header
	//! \param cmin=0 size_t  - min allowed cluster size
	//! \param cmax=0 size_t  - max allowed cluster size, 0 means any size
    //! \param verbose=false bool  - print intermediate results to the stdout
    //! \return bool  - the collection is loaded successfully
	static NodeBase load(const char* filename, float membership=1
		, AggHash* ahash=nullptr, size_t cmin=0, size_t cmax=0, bool verbose=false);
};

//! Collection matching kind base
using MatchBase = uint8_t;

//! \brief Collection matching kind
enum struct Match: MatchBase {
	NONE = 0,  //!< Note initialized
	WEIGHTED,  //!< Weighted matching by the number of members in each cluster (macro weighting)
	UNWEIGHTED,  //!< Unweighted matching of each cluster (micro weighting)
	COMBINED  //!< Combined of macro and micro weightings using geometric mean
};

//! \brief String representation of the Match
//! \relates Match
//!
//! \param mkind Match  - the value to be converted
//! \return string  - string value
string to_string(Match mkind);

//! \brief The matching includes weighted match
//! \relates Match
//!
//! \param m Match  - matching kind
//! \return bool  - weighted matching included
bool xwmatch(Match m) noexcept;

//! \brief The matching includes unweighted match
//! \relates Match
//!
//! \param m Match  - matching kind
//! \return bool  - unweighted matching included
bool xumatch(Match m) noexcept;

//! Precision and recall
struct PrcRec {
	Prob prc;  //!< Precision
	Prob rec;  //!< Recall

	// Explicit members initialization by value to avoid uninitialized members
	PrcRec(Prob prc=0, Prob rec=0): prc(prc), rec(rec)  {}
};

//! Collection describing cluster-node relations
//! \tparam Count  - arithmetic counting type
template <typename Count>
class Collection: public NodeBaseI {
public:
	using CollectionT = Collection<Count>;
	//! Overlaps / multi-resolutions evaluation flag
	constexpr static bool  m_overlaps = is_floating_point<Count>::value;
	//! Accumulated contribution
	using AccCont = conditional_t<m_overlaps, Count, AccId>;
	//! Clusters matching matrix
	using ClustersMatching = SparseMatrix<Cluster<Count>*, AccCont>;  // Used only for NMI
	using ClsLabels = ClustersLabels<Count>;
private:
	// ATTENTNION: Collection manages the memory of the m_cls
	ClusterPtrs<Count>  m_cls;  //!< Clusters
	NodeClusters<Count>  m_ndcs;  //!< Node clusters relations
	size_t  m_ndshash;  //!< Nodes hash (of unique node ids only, not all members), 0 means was not evaluated
	//mutable bool  m_dirty;  //!< The cluster members contribution is not zero (should be reseted on reprocessing)
	//! Sum of contributions of all members in each cluster
	mutable AccCont  m_contsum;  // Used by NMI only, marked also by overlapping F1
protected:
    //! Default constructor
	Collection(): m_cls(), m_ndcs(), m_ndshash(0), m_contsum(0)  {}  //, m_dirty(false)  {}

	// Note: Actual for NMI and overlapping F1
	//! \brief Initialized cluster members contributions
	//!
	//! \param cn const CollectionT&  - target collection to initialize cluster
	//! members contributions
	//! \return void
	static void initconts(const CollectionT& cn) noexcept;
public:
	~Collection();

    //! \brief The number of clusters
    //!
    //! \return Id  - the number of clusters in the collection
	Id clsnum() const noexcept  { return m_cls.size(); }

    //! \brief The number of nodes
    //!
    //! \return Id  - the number of nodes in the collection
	Id ndsnum() const noexcept  { return m_ndcs.size(); }

	//! \copydoc NodeBaseI::nodeExists(Id nid) const noexcept
	bool nodeExists(Id nid) const noexcept  { return m_ndcs.count(nid); }

	//! \brief Load collection from the CNL file
	//! \pre All clusters in the file are expected to be unique and not validated for
	//! the mutual match until makeunique is set
	//!
	//! \param filename const char*  - name of the input file
	//! \param makeunique=false bool  - ensure that clusters contain unique members by
	//! 	removing the duplicates
	//! \param membership=1 float  - expected membership of the nodes, >0, typically >= 1.
	//! Used only for the node container preallocation to estimate the number of nodes
	//! if not specified in the file header
    //! \param ahash=nullptr AggHash*  - resulting hash of the loaded
    //! member ids base (unique ids only are hashed, not all ids) if not nullptr
	//! \param const nodebase=nullptr NodeBaseI*  - node base to filter-out nodes if required
	//! \param lostcls=nullptr RawIds*  - indices of the lost clusters during the node base
	//! synchronization
	//! \param verbose=false bool  - print the number of loaded nodes to the stdout
    //! \return CollectionT  - the collection is loaded successfully
	static CollectionT load(const char* filename, bool makeunique=false
		, float membership=1, AggHash* ahash=nullptr, const NodeBaseI* nodebase=nullptr
		, RawIds* lostcls=nullptr, bool verbose=false);

    //! \brief Transfer collection data
    //! \post This collection becomes empty
    //!
    //! \tparam FIRST bool  - fill first of second node clusters relations container
    //!
    //! \param cls RawClusters&  - raw clusters to be extended
    //! \param nds NodeRClusters&  - node raw clusters relations to be extended
    //! \return void
    template <bool FIRST>
	void transfer(RawClusters& cls, NodeRClusters& ndrcs);

    //! \brief Clear cluster counters
    //!
    //! \return void
	void clearcounts() const noexcept;

    //! \brief Label collection clusters according to the ground-truth cluster indices
    //!
    //! \param gt const CollectionT&  - ground-truth cluster collection
    //! \param cn const CollectionT&  - processing cluster collection
	//! \param lostcls const RawIds&  - indices of the lost clusters during the node base
	//! synchronization
    //! \param prob bool  - Partial Probabilities or F1 (harmonic) matching policy
    //! \param weighted=true bool  - weight labels by the number of instances or
    //! treat each label equally
    //! \param flname=nullptr const char*  - resulting label indices filename (.cll format)
    //! \param verbose=false bool  - print intermediate results to the stdout
    //! \return PrcRec  - resulting precision and recall for the labeled items
	static PrcRec label(const CollectionT& gt, const CollectionT& cn, const RawIds& lostcls
		, bool prob, bool weighted=true, const char* flname=nullptr, bool verbose=false);

	//! \brief Specified F1 evaluation of the Greatest (Max) Match for the
	//! multi-resolution clustering with possibly unequal node base
	//!
	//! Supported F1 measures are F1p <= F1h <= F1s, where:
	//! - F1p  - Harmonic mean of the [weighted] average of partial probabilities,
	//! 	the most discriminative and satisfies the largest number of the Formal
	//! 	Constraints (homogeneity, completeness, rag bag,  size/quantity, balance);
	//! - F1h  - Harmonic mean of the [weighted] average of F1s;
	//! - F1a  - Average F1-Score, i.e. arithmetic mean (average) of the [weighted]
	//! 	average of F1s, the least discriminative and satisfies the lowest number
	//! 	of the Formal Constraints.
	//!
	//! of the Greatest (Max) Match [Weighted] Average Harmonic Mean evaluation
	//! \note Undirected (symmetric) evaluation
	//!
	//! \param cn1 const CollectionT&  - first collection
	//! \param cn2 const CollectionT&  - second collection
    //! \param kind F1  - kind of F1 to be evaluated
    //! \param rec Prob&  - recall of cn2 relative to the ground-truth cn1 or
    //! 0 if the matching strategy does not have the precision/recall notations
    //! \param prc Prob&  - precision of cn2 relative to the ground-truth cn1 or
    //! 0 if the matching strategy does not have the precision/recall notations
    //! \param mkind=Match::WEIGHTED Match  - matching kind
    //! \param verbose=false bool  - print intermediate results to the stdout
	//! \return Prob  - resulting F1_gm
	static Prob f1(const CollectionT& cn1, const CollectionT& cn2, F1 kind
		, Prob& rec, Prob& prc, Match mkind=Match::WEIGHTED, bool verbose=false);

	//! \brief NMI evaluation
	//! \note Undirected (symmetric) evaluation
	//!
	//! \param cn1 const CollectionT&  - first collection
	//! \param cn2 const CollectionT&  - second collection
    //! \param expbase=false bool  - use ln (exp base) or log2 (Shannon entropy, bits)
    //! for the information measuring
    //! \param verbose=false bool  - perform additional verification and print details
	//! \return RawNmi  - resulting NMI
	static RawNmi nmi(const CollectionT& cn1, const CollectionT& cn2, bool expbase=false
		, bool verbose=false);
protected:
	// Label related functions -------------------------------------------------
    //! \brief Mark clusters of the argument collection with the labels
    //! \note For EACH label the best matching cluster is identified. Mutual match
    //! is not applied to guarantee coverage of the all ground-truth clusters to
    //! have meaningful F1
    //!
    //! \param cn const CollectionT&  - the collection to be labeled
    //! \param prob bool  - match labels by the Partial Probabilities or F1;
    //! prob maximizes gain otherwise loss is minimized and F1 is maximized
    //! \param weighted=true bool  - weight labels by the number of instances or
    //! treat each label equally
    //! \param csls=nullptr ClsLabels*  - resulting labels as clusters of the
    //! ground-truth collection if not nullptr
    //! \return PrcRec  - resulting average over all labels Precision and Recall
    //! for all nodes of the marked clusters, where each label can be assigned
    //! to multiple cn clusters and then all nodes of that clusters are matched
    //! to the ground truth cluster (label) nodes
	PrcRec mark(const CollectionT& cn, bool prob, bool weighted=true, ClsLabels* csls=nullptr) const;

	// F1-related functions ----------------------------------------------------
    //! \brief Average of the maximal matches (by F1 or partial probabilities)
    //! relative to the specified collection FROM this one
    //! \note External cn collection can have unequal node base and overlapping
    //! clusters on multiple resolutions. Small collection relative to the average
    //! or average relative to huge might yield best matching F1 equal to 1, but
    //! then the back match should be small.
    //! \attention Directed (non-symmetric) evaluation
    //!
    //! \param gmats const Probs&  - greatest (max) matching with another collection
    //! \param weighted bool  - weighted average by cluster size
    //! \return AccProb  - resulting max average match value from this collection
    //! to the specified one (DIRECTED)
	inline AccProb avggms(const Probs& gmats, bool weighted) const;  // const CollectionT& cn

    //! \brief Greatest (Max) matching value (F1 or partial probability) for each cluster
    //! to the corresponding clusters of the specified collection
    //! \note External cn collection can have unequal node base and overlapping
    //! clusters on multiple resolutions
    //! \attention Directed (non-symmetric) evaluation
    //! \post Modifies internal state of the collection
    //!
    //! \param cn const CollectionT&  - collection to compare with
    //! \param prob bool  - evaluate partial probability instead of F1
    //! \return Probs - resulting max F1 or partial probability for cluster
    //! (all member nodes are considered in the cluster)
	Probs gmatches(const CollectionT& cn, bool prob) const;

	// NMI-related functions ---------------------------------------------------
	//! \brief NMI evaluation considering overlaps, multi-resolution and possibly
	//! unequal node base
	//! \note Undirected (symmetric) evaluation
    //!
    //! \param cn const CollectionT&  - collection to compare with
    //! \param expbase bool  - use ln (exp base) or log2 (Shannon entropy, bits)
    //! for the information measuring
    //! \return RawNmi  - resulting NMI
	RawNmi nmi(const CollectionT& cn, bool expbase) const;

    //! \brief Clear contributions in each cluster and optionally
    //! evaluate the clusters matching
    //!
    //! \param cn const CollectionT&  - foreign collection to be processed with this one
    //! \param[out] clsmm=nullptr ClustersMatchingT*  - clusters matching matrix to be filled
    //! \return AccCont  - sum of all values of the clsmm matrix if specified
	AccCont evalconts(const CollectionT& cn, ClustersMatching* clsmm=nullptr) const;

    //! \brief Clear contributions in each cluster
    //!
    //! \return void
	void clearconts() const noexcept;
};

// Accessory functions ---------------------------------------------------------
//! \brief Compile time pair selector
//!
//! \tparam FIRST bool  - select .first or .second
//!
//! \param pr P&  - pair
//! \return selected field
template <bool FIRST, typename P>
enable_if_t<FIRST, typename P::first_type>& pairsel(P& pr) noexcept  { return pr.first; }

template <bool FIRST, typename P>
enable_if_t<!FIRST, typename P::second_type>& pairsel(P& pr) noexcept  { return pr.second; }

//! \brief Parse decimal c-string as id
//!
//! \param str char*  - id string
//! \return Id  - id value
Id  parseId(char* str);

//! \brief Harmonic mean
//! \note a + b = 0 are threated correctly resulting 0
//!
//! \param a AccProb  - first item
//! \param b AccProb  - second item
//! \return AccProb  - resulting mean
AccProb hmean(AccProb a, AccProb b) noexcept;

//! \brief Geometric mean
//!
//! \param a AccProb  - first item
//! \param b AccProb  - second item
//! \return AccProb  - resulting mean
AccProb gmean(AccProb a, AccProb b) noexcept;

//! \brief Arithmetic mean (average)
//!
//! \param a AccProb  - first item
//! \param b AccProb  - second item
//! \return AccProb  - resulting mean
AccProb amean(AccProb a, AccProb b) noexcept;

#endif // INTERFACE_H
