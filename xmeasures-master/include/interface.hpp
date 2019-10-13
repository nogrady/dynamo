//! \brief Extrinsic measures evaluation interface.
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2017-02-13

#include <cstring>  // strlen, strtok
//#include <cmath>  // sqrt
#include <algorithm>

#include "operations.hpp"
#include "interface.h"


using std::out_of_range;
using std::overflow_error;
using std::invalid_argument;
using std::to_string;
//using std::bitset;
using std::min;
using std::max;
using std::sort;
using std::unique;
//using std::adjacent_find;
#if VALIDATE >= 1
using std::domain_error;
#endif // VALIDATE 1
using namespace daoc;


template <bool EXT=false>
Prob omega(const NodeRClusters& ndrcs, const RawClusters& cls1, const RawClusters& cls2)
{
	//using AccContrib = conditional_t<EXT, AccProb, AccId>;  // Accumulated contribution of the pair to OI
	AccId  oh = 0;  // Observed contribution of the high (top right) part of the matrix
	RawIds  icount(cls1.size(), 0);
	RawIds  jcount(cls2.size(), 0);
	const auto encs = ndrcs.end();
	for(auto iincs = ndrcs.begin(); iincs != encs;) {
		auto& incs = *iincs;
		// Note: diagonal items for ovp always equal to 1 (or contrib_max) and will be added later
		for(auto ijncs = ++iincs; ijncs != encs; ++ijncs) {
			const Id  inum = mutualnum(&incs.second.first, &ijncs->second.first);
			const Id  jnum = mutualnum(&incs.second.second, &ijncs->second.second);
			oh += inum == jnum;
			++icount[inum];
			++jcount[jnum];
		}
	}
	AccId  eh = 0;  // Expected contribution of the high (top right) part of the matrix
	const Id  csize = min(jcount.size(), icount.size());
	for(Id i = 0; i < csize; ++i) {
		//printf("> omega() #%u: %lu, %u, %u\n"
		//	, i, AccId(icount[i]) * jcount[i], icount[i], jcount[i]);
		eh += AccId(icount[i]) * jcount[i];
	}
	// The number of pairs = nodes_num * (nodes_num - 1) / 2
	const AccId  npairs = ndrcs.size() * (ndrcs.size() - 1) >> 1;
	const AccProb  enorm = AccProb(eh) / npairs;
	//printf("> omega() observed: %lu, expected normalized: %G, %lu nodes\n"
	//	, oh, enorm, ndrcs.size());
	return (oh - enorm) / (npairs - enorm);
}

template <>
Prob omega<true>(const NodeRClusters& ndrcs, const RawClusters& cls1, const RawClusters& cls2)
{
	AccProb  oh = 0;  // Observed contribution of the high (top right) part of the matrix
	RawIds  icount(cls1.size(), 0);
	RawIds  jcount(cls2.size(), 0);
	const auto encs = ndrcs.end();
	for(auto iincs = ndrcs.begin(); iincs != encs;) {
		auto& incs = *iincs;
		// Note: diagonal items for ovp always equal to 1 (or contrib_max) and will be added later
		for(auto ijncs = ++iincs; ijncs != encs; ++ijncs) {
			const Id  inum = mutualnum(&incs.second.first, &ijncs->second.first);
			const Id  jnum = mutualnum(&incs.second.second, &ijncs->second.second);
			oh += inum == jnum ? 1 : AccProb(min(inum, jnum)) / max(inum, jnum);
			//eh += AccProb(max<Id>(inum, 1)) * max<Id>(jnum, 1);
			++icount[inum];
			++jcount[jnum];
		}
	}
	AccProb  eh = 0;  // Expected contribution of the high (top right) part of the matrix
	const Id  csize = min(jcount.size(), icount.size());
	for(Id i = 0; i < csize; ++i)
		eh += AccId(icount[i]) * jcount[i];
	// Consider remained accumulated counts
	RawIds& rcount = icount.size() > csize ? icount : jcount;
	const Id  rcsize = rcount.size();
	for(Id i = csize; i < rcsize; ++i)
		eh += rcount[i];
	// The number of pairs = nodes_num * (nodes_num - 1) / 2
	const AccId  npairs = ndrcs.size() * (ndrcs.size() - 1) >> 1;
	const AccProb  enorm = AccProb(eh) / npairs;
	//printf("> omega() observed: %G, expected normalized: %G, %lu nodes\n"
	//	, oh, enorm, ndrcs.size());
	return (oh - enorm) / (npairs - enorm);
}

// Cluster definition ----------------------------------------------------------
template <typename Count>
Cluster<Count>::Cluster(): members(), counter(), mbscont()
{}

// SparseMatrix definitions ----------------------------------------------------
template <typename Index, typename Value>
SparseMatrix<Index, Value>::SparseMatrix(Id rows)
{
	if(rows)
		BaseT::reserve(rows);  // Preallocate hash map
}

template <typename Index, typename Value>
Value& SparseMatrix<Index, Value>::operator ()(Index i, Index j)
{
	auto& rowi = (*this)[i];
	auto ir = fast_ifind(rowi, j, bsObjOp<RowVecItem<Index, Value>>);
	if(ir == rowi.end() || ir->pos != j)
		ir = rowi.emplace(ir, j);  // Transparently Insert a new element
	return ir->val;
}

template <typename Index, typename Value>
template <typename T, enable_if_t<sizeof(T) <= sizeof(void*)>*>
Value SparseMatrix<Index, Value>::operator()(Index i, Index j) const noexcept
{
#if VALIDATE >= 2
	auto irowi = BaseT::find(i);
	if(irowi == BaseT::end())
		fprintf(stderr, "ERROR operator(), row #%u does not exist\n", i);
	auto ie = fast_ifind(*irowi, j, bsObjOp<RowVecItem<Index, Value>>)->val;
	if(irowi == BaseT::end())
		fprintf(stderr, "ERROR operator(), element #u in the row #%u does not exist\n", j, i);
	return ie->val;
#else
	return fast_ifind(*find(i), j, bsObjOp<RowVecItem<Index, Value>>)->val;
#endif // VALIDATE
}

template <typename Index, typename Value>
template <typename T, enable_if_t<(sizeof(T) > sizeof(void*)), bool>*>
const Value& SparseMatrix<Index, Value>::operator()(Index i, Index j) const noexcept
{
#if VALIDATE >= 2
	auto irowi = BaseT::find(i);
	if(irowi == BaseT::end())
		fprintf(stderr, "ERROR operator() 2, row #%u does not exist\n", i);
	auto ie = fast_ifind(*irowi, j, bsObjOp<RowVecItem<Index, Value>>)->val;
	if(irowi == BaseT::end())
		fprintf(stderr, "ERROR operator(), element #u in the row #%u does not exist\n", j, i);
	return ie->val;
#else
	return fast_ifind(*find(i), j, bsObjOp<RowVecItem<Index, Value>>)->val;
#endif // VALIDATE
}

template <typename Index, typename Value>
template <typename T, enable_if_t<sizeof(T) <= sizeof(void*)>*>
Value SparseMatrix<Index, Value>::at(Index i, Index j)
{
	auto& rowi = BaseT::at(i);
	auto ie = fast_ifind(rowi, j, bsObjOp<RowVecItem<Index, Value>>);
	if(ie == rowi.end() || ie->pos != j)
		throw out_of_range("at(), attempt to access nonexistent element #"
			+ to_string(j) + " at the row #" + to_string(i) + "\n");
	return ie->val;
}

template <typename Index, typename Value>
template <typename T, enable_if_t<(sizeof(T) > sizeof(void*)), bool>*>
const Value& SparseMatrix<Index, Value>::at(Index i, Index j)
{
	auto& rowi = BaseT::at(i);
	auto ie = fast_ifind(rowi, j, bsObjOp<RowVecItem<Index, Value>>);
	if(ie == rowi.end() || ie->pos != j)
		throw out_of_range("at() 2, element #" + to_string(j) + " in the row #"
			+ to_string(i) + " does not exist\n");
	return ie->val;
}

// Collection definitions ------------------------------------------------------
template <typename Count>
Collection<Count>::~Collection()
{
	// Note: exception on double memory release occurs if m_cls contain non-unique pointers,
	// which should never happen
    for(auto cl: m_cls)
		delete cl;
	m_cls.clear();
}

template <typename Count>
Collection<Count> Collection<Count>::load(const char* filename, bool makeunique, float membership
	, ::AggHash* ahash, const NodeBaseI* nodebase, RawIds* lostcls, bool verbose)
{
	Collection  cn;  // Return using NRVO, named return value optimization

	// Open file
	NamedFileWrapper  file(filename, "r");
	if(!file) {
		perror(string("ERROR load(), failed on opening ").append(filename).c_str());
		return cn;
	}

	constexpr size_t  FILESIZE_INVALID = size_t(-1);
	const size_t fsize = file.size();
	if(!fsize) {
		fputs(("WARNING load(), the file '" + file.name()
			+ " is empty, skipped\n").c_str(), stderr);
		return cn;
	}

	// Load clusters
	// Note: CNL [CSN] format only is supported
	size_t  csnum = 0;  // The number of clusters
	size_t  nsnum = 0;  // The number of nodes
	// Note: strings defined out of the cycle to avoid reallocations
	StringBuffer  line;  // Reading line
	// Parse header and read the number of clusters if specified
	parseCnlHeader(file, line, csnum, nsnum);

	// Estimate the number of nodes in the file if not specified
	if(!nsnum) {
		if(fsize != FILESIZE_INVALID) {  // File length fetching failed
			nsnum = estimateCnlNodes(fsize, membership);
#if TRACE >= 2
			fprintf(stderr, "load(), %lu estimated nodes from %lu bytes\n", nsnum, fsize);
#endif // TRACE
		} else if(csnum)
			nsnum = 2 * csnum; // / membership;  // Note: use optimistic estimate instead of pessimistic (square / membership) to not overuse the memory
	}
	// Estimate the number of clusters in the file if not specified
	if(!csnum && nsnum)
		csnum = estimateClusters(nsnum, membership);
#if TRACE >= 2
	fprintf(stderr, "load(), expected %lu clusters, %lu nodes from %lu input bytes\n"
		, csnum, nsnum, fsize);
	if(nodebase)
		fprintf(stderr, "load(), nodebase provided with %u nodes\n", nodebase->ndsnum());
#endif // TRACE

	// Preallocate space for the clusters and nodes
	if(cn.m_cls.capacity() < csnum)  //  * cn.m_cls.max_load_factor()
		cn.m_cls.reserve(csnum);
	if(cn.m_ndcs.bucket_count() * cn.m_ndcs.max_load_factor() < nsnum)
		cn.m_ndcs.reserve(nsnum);

	// Parse clusters
	// ATTENTION: without '\n' delimiter the terminating '\n' is read as an item
	constexpr char  mbdelim[] = " \t\n";  // Delimiter for the members
	// Estimate the number of chars per node, floating number
	const float  ndchars = nsnum ? (fsize != FILESIZE_INVALID
		? fsize / float(nsnum) : log10(float(nsnum) + 1))  // Note: + 1 to consider the leading space
		: 1.f;
#if VALIDATE >= 2
	//fprintf(stderr, "load(), ndchars: %.4G\n", ndchars);
	assert(ndchars >= 1 && "load(), ndchars invalid");
#endif // VALIDATE
	::AggHash  mbhash;  // Nodes hash (only unique nodes, not all the members)
	ClusterHolder<Count>  chd(new Cluster<Count>());
	do {
		// Skip cluster id if specified and parse first node id
		char *tok = strtok(line, mbdelim);  // const_cast<char*>(line.data())

		// Skip comments
		if(!tok || tok[0] == '#')
			continue;
		// Skip the cluster id if present
		if(tok[strlen(tok) - 1] == '>') {
			const char* cidstr = tok;
			tok = strtok(nullptr, mbdelim);
			// Skip empty clusters, which actually should not exist
			if(!tok) {
				fprintf(stderr, "WARNING load(), empty cluster"
					" exists: '%s', skipped\n", cidstr);
				continue;
			}
		}

		// Parse remained node ids and load cluster members
		Cluster<Count>* const  pcl = chd.get();
		auto& members = pcl->members;
		members.reserve(line.length() / ndchars);  // Note: strtok() does not affect line.length()
		do {
			// Note: only the node id is parsed, share part is skipped if exists,
			// but potentially can be considered in NMI and F1 evaluation.
			// In the latter case abs diff of shares instead of co occurrence
			// counting should be performed.
			Id  nid = parseId(tok);
//			Id  nid = strtoul(tok, nullptr, 10);
//#if VALIDATE >= 2
//			if(!nid && tok[0] != '0') {
//				fprintf(stderr, "WARNING load(), conversion error of '%s' into 0: %s\n"
//					, tok, strerror(errno));
//				continue;
//			}
//#endif // VALIDATE
			// Filter out nodes if required
			if(nodebase && !nodebase->nodeExists(nid))
				continue;
			members.push_back(nid);
			auto& ncs = cn.m_ndcs[nid];
			// Update hash if required
			if(ncs.empty())
				mbhash.add(nid);
			ncs.push_back(pcl);
		} while((tok = strtok(nullptr, mbdelim)));
		if(!members.empty()) {
			members.shrink_to_fit();  // Free over reserved space
#if VALIDATE <= 1
			if(makeunique)
#endif // VALIDATE
			{
				// Ensure or validate that members are unique
				sort(members.begin(), members.end());
				const auto im = unique(members.begin(), members.end());
				//const auto im = adjacent_find(members.begin(), members.end());
				if(im != members.end()) {
					fprintf(stderr, "WARNING load(), #%lu cluster contained %lu duplicated members, corrected.\n"
						, cn.m_cls.size(), distance(im, members.end()));
					// Remove associated clusters
					for(auto jm = im; jm != members.end(); ++jm)
						cn.m_ndcs[*jm].pop_back();
					// Remove the tail of duplicated node ids
					members.erase(im, members.end());
					//fprintf(stderr, "WARNING load(), #%lu cluster contains duplicated member #%lu: %u\n"
					//	, cn.m_cls.size(), distance(members.begin(), im), *im);
					//throw invalid_argument("load(), the cluster contains duplicated members\n");
				}
			}
			//for(auto v: members)
			//	printf(" %u", v);
			//puts("");
			cn.m_cls.push_back(chd.release());
			// Start filling a new cluster
			chd.reset(new Cluster<Count>());
		} else if(lostcls)
			lostcls->push_back(lostcls->size() + cn.m_cls.size());
	} while(line.readline(file));
	// Save some space if it is essential
	if(cn.m_cls.size() < cn.m_cls.capacity() / 2)
		cn.m_cls.shrink_to_fit();
	// Rehash the clusters and nodes for faster traversing if required
	//if(cn.m_cls.size() < cn.m_cls.bucket_count() * cn.m_cls.max_load_factor() / 2)
	//	cn.m_cls.reserve(cn.m_cls.size());
	if(cn.m_ndcs.size() < cn.m_ndcs.bucket_count() * cn.m_ndcs.max_load_factor() / 2)
		cn.m_ndcs.reserve(cn.m_ndcs.size());
	// Assign hash to the results
	cn.m_ndshash = mbhash.hash();  // Note: required to identify the unequal node base in the processing collections
	if(ahash)
		*ahash = move(mbhash);
#if TRACE >= 2
	printf("load(), loaded %lu clusters (capacity: %lu, overhead: %0.2f %%) and"
		" %lu nodes (reserved %lu buckets, overhead: %0.2f %%) with hash %lu from %s\n"
		, cn.m_cls.size(), cn.m_cls.capacity()
		, cn.m_cls.size() ? float(cn.m_cls.capacity() - cn.m_cls.size()) / cn.m_cls.size() * 100
			: numeric_limits<float>::infinity()
		, cn.m_ndcs.size(), cn.m_ndcs.bucket_count()
		, cn.m_ndcs.size() ? float(cn.m_ndcs.bucket_count() - cn.m_ndcs.size()) / cn.m_ndcs.size() * 100
			: numeric_limits<float>::infinity()
		, cn.m_ndshash, file.name().c_str());
#elif TRACE >= 1
	if(verbose)
		printf("load(), loaded %lu clusters %lu nodes from %s\n", cn.m_cls.size()
			, cn.m_ndcs.size(), file.name().c_str());
#endif

	return cn;
}

template <typename Count>
template <bool FIRST>
void Collection<Count>::transfer(RawClusters& cls, NodeRClusters& ndrcs)
{
	// Clear node clusters relations of the collection
	m_ndcs.clear();
	m_ndcs.rehash(0);
	// Prepare target containers
	cls.reserve(m_cls.size());
	ndrcs.reserve(m_ndcs.size());

	for(auto cl: m_cls) {
		cls.push_back(move(cl->members));
		delete cl;
		auto& craw = cls.back();
		for(auto nd: craw)
			pairsel<FIRST>(ndrcs[nd]).push_back(&craw);
	}
	// Order node clusters
	for(auto& val: ndrcs)
		sort(pairsel<FIRST>(val.second).begin(), pairsel<FIRST>(val.second).end(), cmpBase<RawCluster*>);

	// Clear collection clusters
	m_cls.clear();
	m_cls.shrink_to_fit();
}

template <typename Count>
void Collection<Count>::clearcounts() const noexcept
{
	for(auto cl: m_cls)
		cl->counter.clear();
}

template <typename Count>
PrcRec Collection<Count>::label(const CollectionT& gt, const CollectionT& cn, const RawIds& lostcls
	, bool prob, bool weighted, const char* flname, bool verbose)
{
	// Initialized accessory data for evaluations if has not been done yet
	// (nmi also initializes mbscont)
	if(is_floating_point<Count>::value && !(gt.m_contsum && cn.m_contsum)) {  // Note: strict ! is fine here
		// Evaluate members contributions
		initconts(gt);
		initconts(cn);
	}

	// Label cn2 clusters from the ground-truth clusters cn1 and evaluate F1
	// including precision and recall
#if TRACE >= 3
	fputs("label(), Labeling the target collection\n", stderr);
#endif // TRACE
	// ATTENTION: mark() (and gmatches() changes internal state of the collection parameter, so
	// it should be called only once for each collection and with the same value of prob
	// Note: it's more convenient for the subsequent processing to assign labels to the clusters
	ClsLabels  csls;  // Clusters labels to be outputted
	auto pr = gt.mark(cn, prob, weighted, &csls);

	// Evaluate labels for each node by the node clusters and cluster labels

	// Output the resulting clusters labels to the specified file if required
	if(flname) {
		// Create the output file
		NamedFileWrapper  flbs(flname, "w");
		if(flbs) {
			fprintf(flbs, "# Clusters: %lu, Labels: %lu\n", csls.size(), gt.m_cls.size());
			// Map labels to their indices
			using LabelIndices = unordered_map<Cluster<Count>*, Id>;
			LabelIndices  lids;
			lids.reserve(gt.m_cls.size());
            for(size_t i = 0; i < gt.m_cls.size(); ++i)
				lids[gt.m_cls[i]] = i;
			// Output clusters marked with label indices
			for(auto cl: cn.m_cls) {
				auto iclbs = csls.find(cl);
				if(iclbs == csls.end()) {
					fprintf(flbs, "-\n");
					continue;
				}
				for(auto lb: iclbs->second)
					fprintf(flbs, "%u ", lids.at(lb));
				fputs("\n", flbs);
			}
		} else fprintf(stderr, "WARNING label(), labels output is omitted"
			": '%s' file can't be created\n", flname);
	}

	return pr;
}

template <typename Count>
PrcRec Collection<Count>::mark(const CollectionT& cn, bool prob, bool weighted, ClsLabels* csls) const
{
	// Reserve space for all clusters of the collection (but not more than the number of labels) to avoid reallocations
	if(csls)
		csls->reserve(min(m_cls.size(), cn.m_cls.size()));
	// Function evaluating value of the match
	auto fmatch = prob ? &Cluster<Count>::pprob : &Cluster<Count>::f1;
	// Traverse all gt clusters (labels)
	using MarkCands = unordered_set<Cluster<Count>*>;
	MarkCands  mcands;  // Marking candidates
	// Marked nodes of multiple clusters merged to the flat set
	using MarkNodes = unordered_set<Id>;
	MarkNodes  mnds;
	// Aggregated precision and recall
	AccProb  prc = 0;
	AccProb  rec = 0;
#if TRACE >= 2
	Id  nmlbs = 0;  // The number of labels with multiple clusters
#endif // TRACE
	// The number of missed (non-matched) labels, which is possible only when
	// the node base is not synchronized
	Id lbmissed = 0;
	AccProb  accw = 0;  // Accumulated weight of the labels or just their number (if !weighted)

	for(auto gtc: m_cls) {
		Prob  gmatch = 0; // Greatest value of the match (F1 or partial probability)
		// Traverse all members (node ids)
		for(auto nid: gtc->members) {
			// Find Matching clusters (containing the same member node id) in the foreign collection
			const auto imcls = cn.m_ndcs.find(nid);
			// Consider the case of unequal node base, i.e. missed node
			if(imcls == cn.m_ndcs.end())
				continue;
			for(auto mcl: imcls->second) {
				// Greatest matches (Max F1 or partial probability) for each ground-truth cluster
				// [of this collection, self] (label);
				if(m_overlaps)
					// In case of overlap contributes the smallest share (of the largest number of owners)
					mcl->counter(gtc, AccProb(1)
						/ max(m_ndcs.at(nid).size(), cn.m_ndcs.at(nid).size()));
				else mcl->counter(gtc, 1);
				// Note: only the max value for match is sufficient
				// ATTENTION: F1 compares clusters per-pair, so it is much simpler and
				// has another semantics of contribution for the multi-resolution case
				const Prob  match = (gtc->*fmatch)(mcl->counter(), m_overlaps
					? mcl->cont() : mcl->members.size());
				if(!less<Prob>(match, gmatch)) {
					if(!equal<Prob>(match, gmatch)) {
						gmatch = match;
						mcands.clear();
					}
					mcands.insert(mcl);
				}
			}
		}
#if TRACE >= 3
		// Note: sqrt() is used to provide semantic values, geometric mean of the Precision
		// and Recall, which is >= harmonic mean and <= arithmetic mean
		fprintf(stderr, "  %p (%lu) => %lu cands: %.3G", gtc, gtc->members.size(), mcands.size(), sqrt(gmatch));
#endif // TRACE
		// ATTENTION: update before the iteration skipping
		accw += weighted ? gtc->members.size() : 1;
		// Note: mcands can be empty only if the node base is not synchronized
		//assert(mcands.size() >= 1 && "mark(), each label should be matched to at least one cluster");
		if(mcands.empty()) {
			++lbmissed;
			continue;
		}
		// For the rare case of matching single label to multiple cn clusters, merge nodes
		// of that clusters to evaluate Precision and Recall of the aggregated match of the cluster nodes
		if(mcands.size() >= 2 && mnds.bucket_count() / mnds.max_load_factor()
		< max((*mcands.begin())->members.size(), (*++mcands.begin())->members.size()))
			mnds.reserve(max((*mcands.begin())->members.size(), (*++mcands.begin())->members.size()) * 1.5);
		// Mark candidate clusters with the labels, keep the labels ordered in each cluster
		// Note: on each iteration distinct label (ground-truth cluster) is provided
		if(csls) {
			for(auto cl: mcands) {
				auto& lbs = (*csls)[cl];
				// Note: typically the number of labels is small (1 or 2 even when the clusters are heavily overlapping)
				lbs.insert(linear_ifind(lbs, cl, bsVal<Cluster<Count>*>), gtc);
			}
		}
		// Aggregate clusters nodes
		if(mcands.size() >= 2) {
			for(auto cl: mcands)
				mnds.insert(cl->members.begin(), cl->members.end());
		}
		// Evaluate precision and recall
		if(mnds.empty()) {
			// The label marked a single cluster
			auto& mcl = **mcands.begin();
			AccProb  gm = mcl.counter();  // Matches
			if(weighted)
				gm *= gtc->members.size();
			prc += gm / static_cast<AccProb>(m_overlaps ? mcl.cont() : mcl.members.size());
			rec += gm / static_cast<AccProb>(m_overlaps ? gtc->cont() : gtc->members.size());
#if TRACE >= 3
			printf("  > mark(), gmatch: %G, gm: %G, accumulated prc: %G (mcl cont: %G), rec: %G (gtc cont: %G)\n"
				, gmatch, gm, prc, static_cast<AccProb>(m_overlaps ? mcl.cont() : mcl.members.size())
				, rec, static_cast<AccProb>(m_overlaps ? gtc->cont() : gtc->members.size()));
#endif // TRACE
		} else {
			// The label marked multiple clusters, compared it's nodes with the
			// nodes of the merged respective clusters
#if TRACE >= 2
			++nmlbs;
#endif // TRACE
			AccProb  accgm = 0;  // Matches
			if(m_overlaps) {
				AccProb  accont = 0;  // Accumulated contribution
				for(auto cl: mcands) {
					accgm += cl->counter();
					accont += cl->cont();
				}
#if TRACE >= 2
				assert(!less<Prob>(min<AccProb>(accont, gtc->cont()), accgm)
					&& "mark(), accgm ovp validation failed");
#endif // TRACE
				if(weighted)
					accgm *= gtc->members.size();
				prc += accgm / static_cast<AccProb>(accont);
				// Note: mnds.size() <= mcands.size()
				rec += accgm / static_cast<AccProb>(gtc->cont());
			} else {
				// Evaluate the number of matched nodes from the aggregated clusters (<= sum(cls_matches))
				for(auto nid: gtc->members)
					accgm += mnds.count(nid);
#if TRACE >= 2
				assert(accgm <= min(mnds.size(), gtc->members.size())
					&& "mark(), accgm multires validation failed");
#endif // TRACE
				if(weighted)
					accgm *= gtc->members.size();
				prc += accgm / static_cast<AccProb>(mnds.size());
				rec += accgm / static_cast<AccProb>(gtc->members.size());
			}
			mnds.clear();
		}
		mcands.clear();
	}
#if TRACE >= 3
	fputs("\n", stderr);
#endif // TRACE
	if(lbmissed)
		fprintf(stderr, "WARNING mark(), the number of non-matched labels: %u"
			" (possible only when the node base is not synchronized)\n", lbmissed);
#if TRACE >= 2
	fprintf(stderr, "  >> mark(), multi-cluster labels %u / %lu\n", nmlbs, m_cls.size());
#endif // TRACE
#if VALIDATE >= 2
	if(!weighted)
		assert(equal<Prob>(accw, m_cls.size()) && "mark(), total weight on unweighted eval"
			" should be equal to the number of labels");
#endif // VALIDATE
	return PrcRec(prc / accw, rec / accw);
}

template <typename Count>
void Collection<Count>::initconts(const CollectionT& cn) noexcept
{
	// Skip already initialized collection
	if(cn.m_contsum)
		return;

	for(auto& ndcs: cn.m_ndcs) {
		// ATTENTION: in case of fuzzy (unequal) overlaps the shares are unequal and
		// should be stored in the Collection (ncs clusters member or a map)
		const AccProb  ndshare = AccProb(1) / ndcs.second.size();
		for(auto cl: ndcs.second)
			cl->mbscont += ndshare;
	}
	// Mark that mbscont of clusters are used
	cn.m_contsum = -1;
}

template <typename Count>
Prob Collection<Count>::f1(const CollectionT& cn1, const CollectionT& cn2, F1 kind
	, Prob& rec, Prob& prc, Match mkind, bool verbose)
{
	if(kind == F1::NONE || mkind == Match::NONE) {
		fputs("WARNING f1(), f1 or match kind is not specified, the evaluation is skipped\n", stderr);
		return 0;
	}

	// Initialized accessory data for evaluations if has not been done yet
	// (nmi also initializes mbscont)
	if(is_floating_point<Count>::value && !(cn1.m_contsum && cn2.m_contsum)) {  // Note: strict ! is fine here
		// Evaluate members contributions
		initconts(cn1);
		initconts(cn2);
	}

	const bool  prob = kind == F1::PARTPROB;  // Evaluate by the partial probabilities
#if TRACE >= 3
	fputs("f1(), F1 Max Avg of the first collection\n", stderr);
#endif // TRACE
	// ATTENTION: gmatches() changes internal state of the collection parameter, so
	// it should be called only once for each collection and with the same value of prob
	const auto  gmats1 = cn1.gmatches(cn2, prob);
	const AccProb  f1ga1 = cn1.avggms(gmats1, mkind==Match::WEIGHTED);
	prc = f1ga1;  // cn1 (ground-truth) relative to cn2
#if TRACE >= 3
	fputs("f1(), F1 Max Avg of the second collection\n", stderr);
#endif // TRACE
	const auto  gmats2 = cn2.gmatches(cn1, prob);
	const AccProb  f1ga2 = cn2.avggms(gmats2, mkind==Match::WEIGHTED);
	if(kind != F1::AVERAGE)
		rec = f1ga2;  // cn2 relative to cn1 (ground-truth)
	else prc = rec = 0;  // Note: Precision and recall are not defined for the MF1a
	const AccProb  res = kind != F1::AVERAGE
		? hmean(f1ga1, f1ga2)
		: (f1ga1 + f1ga2) / 2;
#if TRACE <= 1
	if(verbose)
#endif // TRACE
	fprintf(verbose ? stdout : stderr, "f1(), f1ga: %G (f1ga1: %G, f1ga2: %G)\n", res, f1ga1, f1ga2);

	if(mkind == Match::COMBINED) {
		prc = rec = 0;  // There are no precision and recall notations for combined matching strategy
		// ATTENTION: gmats already evaluated and should be reused, their reevaluation would
		// affect internal states of the collections (counters) and requires reset of the state
		const AccProb  f1ga1w = cn1.avggms(gmats1, true);
		const AccProb  f1ga2w = cn2.avggms(gmats2, true);
#if TRACE <= 1
		if(verbose)
#endif // TRACE
		fprintf(verbose ? stdout : stderr, "f1(), f1ga1w: %G, f1ga2w: %G\n", f1ga1w, f1ga2w);
		const AccProb  resw = kind != F1::AVERAGE
			? hmean(f1ga1w, f1ga2w)
			: (f1ga1w + f1ga2w) / 2;
		fprintf(verbose ? stdout : stderr, "f1(), f1gaw: %G\n", resw);
		// Note: geometric mean >= harmonic mean for [0, 1] and yields more indicative values,
		// dropping the value not so much when (usually) the weighted match is larger
		return gmean(res, resw);
	}
	return res;
}

template <typename Count>
AccProb Collection<Count>::avggms(const Probs& gmats, bool weighted) const  // const CollectionT& cn,
{
	AccProb  accgm = 0;
	// ATTENTION: gmatches() changes internal state of the collection parameter, so it should be
	// called only once for each collection and with the same value of prob
	//const auto  gmats = gmatches(cn, prob);

	if(weighted) {
#if VALIDATE >= 1
		assert(gmats.size() == m_cls.size()
			&& "avggms(), matches are not synchronized with the clusters");
#endif // VALIDATE
		AccCont  csizesSum = 0;
		auto icl = m_cls.begin();
		for(auto gm: gmats) {
			// Evaluate members considering their shared contributions
			// ATTENTION: F1 compares clusters per-pair, so it is much simpler and
			// has another semantics of contribution for the multi-resolution case
			AccCont  ccont = m_overlaps ? (*icl)->cont() : (*icl)->members.size();
#if VALIDATE >= 2
			assert(ccont > 0 && "avggms(), the contribution should be positive");
#endif // VALIDATE
			++icl;
			accgm += gm * ccont;
			csizesSum += ccont;
#if VALIDATE >= 3
			fprintf(stderr, "	wgm %G, gm: %G, w: %G\n", gm * ccont, gm, AccProb(ccont));
#endif // VALIDATE
		}
		// Note: this assert is not the case for the collections containing multiple resolutions
		//assert(!less<AccCont>(csizesSum, csnum, csnum) && "avggms(), invalid sum of the cluster sizes");
#if VALIDATE >= 2
		fprintf(stderr, "avggms(), accgm: %G (%G / %G)\n", accgm / csizesSum, accgm, AccProb(csizesSum));
#endif // VALIDATE
		accgm /= csizesSum;
	} else {
		for(auto gm: gmats) {
#if VALIDATE >= 3
			if(less<decltype(gm)>(1, gm, gmats.size()))
				throw overflow_error("avggms(), gm E (0, 1] is out of range: " + to_string(gm) + "\n");
#endif // VALIDATE
			accgm += gm;
		}
		accgm /= gmats.size();
	}
#if VALIDATE >= 1
	if(lessx<AccProb>(1, accgm, gmats.size()))
		throw overflow_error("avggms(), accgm is invalid (> 1: " + std::to_string(accgm)
			+ "), probably invalid dataset is supplied\n");
#endif // VALIDATE
	return accgm;
}

template <typename Count>
Probs Collection<Count>::gmatches(const CollectionT& cn, bool prob) const
{
	// Greatest matches (Max F1 or partial probability) for each cluster [of this collection, self];
	Probs  gmats;  // Uses NRVO return value optimization
	gmats.reserve(m_cls.size());  // Preallocate the space

	// Function evaluating value of the match
	auto fmatch = prob ? &Cluster<Count>::pprob : &Cluster<Count>::f1;

	// Traverse all clusters in the collection
	for(auto cl: m_cls) {
		Prob  gmatch = 0; // Greatest value of the match (F1 or partial probability)
		//fprintf(stderr, "> gmatches() %#x (counter: %u (%#x), mbscont: %u) %lu mbs\n"
		//	, cl, cl->counter(), cl->counter.origin(), cl->mbscont, cl->members.size());
		//for(auto v: cl->members)
		//	printf(" %u", v);
		//puts("");
		// Traverse all members (node ids)
		for(auto nid: cl->members) {
			// Find Matching clusters (containing the same member node id) in the foreign collection
			const auto imcls = cn.m_ndcs.find(nid);
			// Consider the case of unequal node base, i.e. missed node
			if(imcls == cn.m_ndcs.end())
				continue;
			for(auto mcl: imcls->second) {
				//fprintf(stderr, ">> gmatches() #%u nid, mcl %#x (counter: %u (%#x), mbscont: %u) %lu mbs\n"
				//	, nid, mcl, mcl->counter(), mcl->counter.origin(), mcl->mbscont, mcl->members.size());
				if(m_overlaps)
					// In case of overlap contributes the smallest share (of the largest number of owners)
					mcl->counter(cl, AccProb(1)
						/ max(m_ndcs.at(nid).size(), cn.m_ndcs.at(nid).size()));
				else mcl->counter(cl, 1);
				// Note: only the max value for match is sufficient
				// ATTENTION: F1 compares clusters per-pair, so it is much simpler and
				// has another semantics of contribution for the multi-resolution case
				const Prob  match = (cl->*fmatch)(mcl->counter(), m_overlaps
					? mcl->cont() : mcl->members.size());
				if(gmatch < match)  // Note: <  usage is fine here
					gmatch = match;
			}
		}
		// Note: sqrt() is used to provide semantic values, geometric mean of the Precision
		// and Recall, which is >= harmonic mean and <= arithmetic mean
		if(prob)
			gmatch = sqrt(gmatch);
		gmats.push_back(gmatch);
#if TRACE >= 3
		fprintf(stderr, "  %p (%lu): %.3G", cl, cl->members.size(), gmatch);
#endif // TRACE
	}
#if TRACE >= 3
	fputs("\n", stderr);
#endif // TRACE
	return gmats;
}

template <typename Count>
RawNmi Collection<Count>::nmi(const CollectionT& cn1, const CollectionT& cn2, bool expbase, bool verbose)
{
	RawNmi  rnmi1;
	if(!cn1.clsnum() || !cn2.clsnum())
		return rnmi1;

	rnmi1 = cn1.nmi(cn2, expbase);
#if VALIDATE >= 1
#if VALIDATE < 2
	if(verbose)
#endif // VALIDATE 2
	{
		// Check NMI value for the inverse order of collections
		auto rnmi2 = cn2.nmi(cn1, expbase);
		fprintf(stderr, "nmi(),  mi1: %G, mi2: %G,  dmi: %G\n", rnmi1.mi, rnmi2.mi
			, rnmi1.mi - rnmi2.mi);
		assert((rnmi1.mi - rnmi2.mi) < precision_limit<Prob>() * 2
			&& "nmi(), rnmi is not symmetric, most likely overlaps are present and not considered but this implementation");
	}
#endif // VALIDATE 1
	return rnmi1;
}

template <typename Count>
RawNmi Collection<Count>::nmi(const CollectionT& cn, bool expbase) const
{
	ClustersMatching  clsmm;  //  Clusters matching matrix
	AccCont  cmmsum = evalconts(cn, &clsmm);  // Sum of all values of the clsmm

	RawNmi  rnmi;  // Resulting raw nmi, initially equal to 0
	if(clsmm.empty()) {
		fputs("WARNING nmi(), collection nodes have no any intersection"
			", the collections are totally different", stderr);
		return rnmi;
	}

	// Evaluate probabilities using the clusters matching matrix
	AccProb (*const clog)(AccProb) = expbase ? log : log2;  // log to be used

    //! \brief Information content
    //! \pre val > 0, capacity >= val
    //!
    //! \param val AccProb  - contributing value
    //! \param capacity AccProb  - total capacity
    //! \return AccProb  - resulting information content H(val, capacity)
	auto infocont = [clog](AccProb val, AccProb capacity) -> AccProb {
#if VALIDATE >= 2
		assert(val > 0 && capacity >= val  // && prob > 0 && prob <= 1
			&& "infocont(), invalid input values");
#endif // VALIDATE
		//if(!val)
		//	return 0;
		const AccProb  prob = val / capacity;
		return !equal<Prob>(prob, 1) ? prob * clog(prob) : -1;
	};

	AccProb  h12 = 0;  // Accumulated mutual probability over the matrix, I(cn1, cn2) in exp base
	AccProb  h1 = 0;  // H(cn1) - information size of the cn1 in exp base
	// Traverse matrix of the mutual information evaluate total mutual probability
	// and mutual probabilities for each collection
#if VALIDATE >= 2
	AccProb  psum = 0;  // Sum of probabilities over the whole matrix
#endif // VALIDATE
#if TRACE >= 2
	fprintf(stderr, "nmi(), nmi is evaluating with base %c, clsmatch matrix sum: %.3G\n"
		, expbase ? 'e' : '2', AccProb(cmmsum));

	#define TRACING_CLSMM_  // Note: local / private macroses are ended with '_'
	fprintf(stderr, "nmi(), clsmm:\n");
#endif // TRACE

	// TODO: To evaluate overlapping or multi-resolution NMI the clusters matching
	// matrix should be renormalized considering best matches and remaining parts
	// (the last part is missed in onmi, which caused it's inapplicability for
	// the multi-resolution evaluations).
	// Order collections by the number of clusters, find best matches to the
	// collection with the largest number of clusters and renormalize the matrix
	// to maximize the match resolving overlaps (or retaining if required).
	// The idea is to avoid NMI penalization of the overlaps by renormalizining
	// overlaps A,B : A',B' as A:A', B:B' moving to them contribution from A:B'
	// and B:A'.
	// The hard part is renormalization in case of the partial match with overlaps
	// in each part.

	// Evaluate NMI (standard NMI if cmmsum is not renormalized)
#if VALIDATE >= 2
	assert(m_contsum > 0 && cn.m_contsum > 0
		&& "nmi(), collection clusters contribution is invalid");
#endif // VALIDATE
	for(const auto& icm: clsmm) {
		// Evaluate information size (content) of the current cluster in the cn1
		// infocont(Accumulated value of the current cluster from cn1, the number of nodes)
		h1 -= infocont(icm.first->cont(), m_contsum);  // ndsnum(), cmmsum

		// Travers row
#ifdef TRACING_CLSMM_
		fprintf(stderr, "%.3G:  ", AccProb(icm.first->cont()));
#endif // TRACING_CLSMM_
		for(auto& icmr: icm.second) {
#if VALIDATE >= 2
			assert(icmr.val > 0 && "nmi(), matrix of clusters matching should contain only positive values");
#endif // VALIDATE
#ifdef TRACING_CLSMM_
			fprintf(stderr, " %G[%.3G]", AccProb(icmr.val), AccProb(icmr.pos->cont()));
#endif // TRACING_CLSMM_
#if VALIDATE >= 2
			// Evaluate mutual probability of the cluster (divide by multiplication of counts of both clusters)
			psum += AccProb(icmr.val) / cmmsum;
#endif // VALIDATE
			// Accumulate total normalized mutual information
			// Note: e base is used instead of 2 to have absolute entropy instead of bits length
			//const auto lval = icmr.val / (icmr.pos->cont() * cprob);  // cprob; icm.first->mbscont
			//mi += mcprob * clog(lval);  // mi = h1 + h2 - h12;  Note: log(a/b) = log(a) - log(b)

			// Note: in the original NMI: AccProb(icmr.val) / nodesNum [ = cmmsum]
			h12 -= infocont(icmr.val, cmmsum);
		}
#ifdef TRACING_CLSMM_
		fputs("\n", stderr);
#endif // TRACING_CLSMM_
	}
#if VALIDATE >= 2
	fprintf(stderr, "nmi(), psum: %G, h12: %G\n", psum, h12);
	assert(equalx(psum, 1., cmmsum) && "nmi(), total probability of the matrix should be 1");
#endif // VALIDATE

	// Evaluate information size cn2 clusters
	AccProb  h2 = 0;  // H(cn2) - information size of the cn1 in exp base
	for(const auto& c2: cn.m_cls)
		h2 -= infocont(c2->cont(), cn.m_contsum);  // cn.ndsnum(), cmmsum

	Prob  mi = h1 + h2 - h12;
#if VALIDATE >= 2
	assert(mi >= -precision_limit<Prob>() && "nmi(), mi is invalid");
#endif // VALIDATE
	if(fabs(mi) < precision_limit<Prob>())
		mi = 0;
	rnmi(mi, h1, h2);  // h1 + h2 - h12;  h12
	//rnmi(h12, h1, h2);  // h1 + h2 - h12;  h12  // ATTENTION: this approach is invalid, it shows higher values for worse matching, but [and] awards overlaps

	// VVV Tracks only shapes of clusters, but not the clusters [members] matching between the collections => need h12
	//rnmi(min(h1, h2), h1, h2);  // h1 + h2 - h12;  h12  // ATTENTION: this approach is invalid
#if TRACE >= 2
	Prob  nmix = rnmi.mi / max(h1, h2);
	fprintf(stderr, "nmi(),  mi: %G (h12: %G),  h1: %G, h2: %G,  NMI_max: %G\n"
		, rnmi.mi, h12, h1, h2, nmix);
#endif // TRACE

	return rnmi;
}

template <typename Count>
auto Collection<Count>::evalconts(const CollectionT& cn, ClustersMatching* pclsmm) const -> AccCont
{
	// Skip evaluations if they already performed (the case of evaluating
	// overlapping F1 after NMI)
	if(!pclsmm && m_contsum && cn.m_contsum) {  // Note: strict values usage is fine here
#if TRACE >= 2
		fputs("evalconts(), contributions evaluation skipped as already performed\n", stderr);
#endif // TRACE
		return 0;
	}

#if VALIDATE >= 2
	// Note: the processing is also fine (but redundant) for the empty collections
	assert(clsnum() && cn.clsnum() && "evalconts(), non-empty collections expected");
#endif // VALIDATE

	// Reset member contributions if not zero
	clearconts();
	cn.clearconts();

    //! \brief Contribution from the member of clusters
    //!
    //! \param owners  - the number of owner clusters of the member
    //! \return AccProb - resulting contribution to a single owner
    // Note: in contribute equally to each upper cluster on each resolution, where
    // only one relevant cluster on each resolution should exists. Multi-resolution
    // structure should have the same node base of each resolution, otherwise
    // overlapping evaluation is more fair. In fact overlaps is generalization of
    // the multi-resolution evaluation, but can't be directly applied to the
    // case of multi-resolution overlapping case, where 2-level (by levels and
	// by overlaps in each level) evaluations are required.
	auto mbcont = [](size_t owners) -> AccCont {
		static_assert(!m_overlaps || is_floating_point<AccCont>::value
			, "mbcont(), invalid types");
		return m_overlaps ? AccCont(1) / owners : AccCont(1);
	};

    //! \brief  Update contribution of to the member-related clusters
    //!
    //! \param cls const ClusterPtrs&  - clusters to be updated
    //! \return AccProb  - contributing value (share in case of overlaps in a single resolution)
    // ATTENTION: Such evaluation is applicable only for non-fuzzy overlaps (equal sharing)
	auto updateCont = [mbcont](const ClusterPtrs<Count>& cls) -> AccCont {
		const AccCont  share = mbcont(cls.size());
		for(auto cl: cls)
			cl->mbscont += share;
		return share;
	};

	ClustersMatching  clsmm = ClustersMatching(clsnum());  // Note: default max_load_factor is 1
	// Total sum of all values of the clsmm matrix, i.e. the number of
	// member nodes in both collections
	AccCont  cmmsum = 0;
	// Evaluate contributions to the clusters of each collection and to the
	// mutual matrix of clusters matching
	// Note: in most cases collections has the same node base, or it was synchronized
	// explicitly in advance, so handle unequal node base as a rare case
	//
	// ATTENTION: For 2+ levels cmmsum equals mbsnum * levs_num, in case of overlaps
	// it also > mbsnum.
	// =>  .mbscont is ALWAYS required except the case of non-overlapping clustering
	// on a single resolution

	// Consider the case of unequal node base, contribution from the missed nodes
	AccCont  econt1 = 0;  // Extra contribution from this collection
	for(auto& ncs: m_ndcs) {
		// Note: always evaluate contributions to the clusters of this collection
		const auto  cls1num = ncs.second.size();  // Note: equals to the number of resolutions for !m_overlaps
		const AccCont  share1 = mbcont(cls1num);
		// Evaluate contribution to the second collection if any
		auto incs2 = cn.m_ndcs.find(ncs.first);
		const auto  cls2num = incs2 != cn.m_ndcs.end() ? incs2->second.size() : 0;
		if(cls2num) {
			AccCont  cont = (m_overlaps ? updateCont(incs2->second)
				: mbcont(incs2->second.size())) * share1;  // Note: shares already divided by clsXnum
			// ATTENTION: share1 != cont * cls2num for !m_overlaps (cls1num - the number of resolutions)
			const AccCont  cont1sum = cont * cls2num;  // Total accumulative contribution from the cl
			cmmsum += cont1sum * cls1num;
			// Update clusters matching matrix
			for(auto cl: ncs.second) {
				cl->mbscont += cont1sum;
				for(auto cl2: incs2->second) {
					clsmm(cl, cl2) += cont;  // Note: contains only POSITIVE values
					if(!m_overlaps)
						cl2->mbscont += cont;
				}
			}
		} else {
			// Note: in this case cls2num and share2 are zero
			econt1 += share1 * cls1num;
			for(auto cl: ncs.second)
				cl->mbscont += share1;
		}

	}
	// Consider the case of unequal node base, contribution from the missed nodes
	AccCont  econt2 = 0;  // Extra contribution from the cn
	if((m_ndshash && cn.m_ndshash && m_ndshash != cn.m_ndshash) || ndsnum() != cn.ndsnum()) {
		for(auto& ncs: cn.m_ndcs) {
			// Skip processed nodes
			if(m_ndcs.count(ncs.first))
				continue;
			// Note: in this case cls1num and share1 are zero
			const auto  cls2num = ncs.second.size();  // Note: equals to the number of resolutions for !m_overlaps
			const AccCont  share2 = mbcont(cls2num);
			econt2 += share2 * cls2num;
			for(auto cl: ncs.second)
				cl->mbscont += share2;
		}
	} else if(!m_ndshash || !cn.m_ndshash)
		fputs("WARNING evalconts(), collection(s) hashes were not evaluated (%lu, %lu)"
			", so some unequal nodes might be skipped on evaluation, which cause"
			" approximate results\n", stderr);

	// Set contributions
	m_contsum = cmmsum + econt1;
	cn.m_contsum = cmmsum + econt2;

#if VALIDATE >= 1
	// Validate if there is anything
	if(!clsmm.empty()) {
#if VALIDATE >= 2
	// Validate sum of vector counts to evaluate probabilities
	// Note: to have symmetric values normalization should be done by the max values in the row / col
	//assert(cmmsum % 2 == 0 && "nmi(), cmmsum should always be even");
#if TRACE >= 2
	#define TRACING_CLSCOUNTS_
	fprintf(stderr, "evalconts(), cls1 counts (%lu): ", m_cls.size());
#endif // TRACE
	m_contsum = econt1;
	for(auto c1: m_cls) {
		m_contsum += c1->cont();
#ifdef TRACING_CLSCOUNTS_
		fprintf(stderr, " %.3G", AccProb(c1->cont()));
#endif // TRACING_CLSCOUNTS_
	}

#ifdef TRACING_CLSCOUNTS_
	fprintf(stderr, "\nevalconts(), cls2 counts (%lu): ", cn.m_cls.size());
#endif // TRACING_CLSCOUNTS_
	cn.m_contsum = econt2;
	for(const auto& c2: cn.m_cls) {
		cn.m_contsum += c2->cont();
#ifdef TRACING_CLSCOUNTS_
		fprintf(stderr, " %.3G", AccProb(c2->cont()));
#endif // TRACING_CLSCOUNTS_
	}
#ifdef TRACING_CLSCOUNTS_
	fputs("\n", stderr);
#endif // TRACING_CLSCOUNTS_
	if(m_overlaps
	&& !(equalx<AccCont>(m_contsum - econt1, ndsnum(), ndsnum())
	&& equalx<AccCont>(cn.m_contsum - econt2, cn.ndsnum(), cn.ndsnum()))) {  // consum equals to the number of nodes for the overlapping case
		fprintf(stderr, "evalconts(), c1csum: %.3G (- %.3G lacked), nds1num: %u"
			", c2csum: %.3G (- %.3G lacked), nds2num: %u,  cmmsum: %.3G\n"
			, AccProb(m_contsum), AccProb(econt1), ndsnum()
			, AccProb(cn.m_contsum), AccProb(econt2), cn.ndsnum(), AccProb(cmmsum));
		assert(0 && "evalconts(), consum validation failed");
	}
#endif // VALIDATE 2
	const bool match1 = equalx<AccProb>(m_contsum - econt1, cmmsum, m_cls.size());
	const bool match2 = equalx<AccProb>(cn.m_contsum - econt2, cmmsum, cn.m_cls.size());
	if((m_ndshash == cn.m_ndshash && m_ndshash && m_ndshash && (!match1 || !match2))  // The same node base
	|| (m_ndshash != cn.m_ndshash && (!match1 && !match2))   // Distinct node base
	) {  // Note: cmmsum should match to either of the sums
		fprintf(stderr, "evalconts(), c1csum: %.3G (- %.3G lacked), c2csum: %.3G (- %.3G lacked), cmmsum: %.3G\n"
			, AccProb(m_contsum), AccProb(econt1), AccProb(cn.m_contsum)
			, AccProb(econt2), AccProb(cmmsum));
		throw domain_error("nmi(), rows accumulation is invalid");
	}
#if TRACE >= 2
	fprintf(stderr, "evalconts(), c1csum: %.3G (- %.3G lacked), c2csum: %.3G (- %.3G lacked), cmmsum: %.3G\n"
		, AccProb(m_contsum), AccProb(econt1), AccProb(cn.m_contsum)
		, AccProb(econt2), AccProb(cmmsum));
#endif // TRACE
	}
#endif // VALIDATE 1
	// Set results
	if(pclsmm)
		*pclsmm = move(clsmm);

	return cmmsum;
}

template <typename Count>
void Collection<Count>::clearconts() const noexcept
{
	// Reset member contributions if not zero
	if(!m_contsum)  // Note: ! is fine here
		return;
	for(auto cl: m_cls)
		cl->mbscont = 0;
	m_contsum = 0;
}
