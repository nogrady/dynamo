//! \brief File IO utils
//!
//!	Interface macro definitions:
//! INCLUDE_STL_FS  - include STL filesystem library under fs namespace. This macros is
//! 	defined to avoid repetitive conditional inclusion of the STL FS.
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2017-02-13

#ifndef FILEIO_H
#define FILEIO_H

#include <cstdint>  // uintX_t
#include <cstdio>  // FILE
#include <utility>  // move
#include <string>
#include <vector>
#include <unordered_set>
// For the template definitions
#include <cstring>  // strtok
#include <cmath>  // sqrt

#ifdef INCLUDE_STL_FS
#if defined(__has_include) && __has_include(<filesystem>)
	#include <filesystem>
	namespace fs = std::filesystem;
#elif defined(__has_include) && __has_include(<experimental/filesystem>)
	#include <experimental/filesystem>
	namespace fs = std::experimental::filesystem;
#else
	#error "STL filesystem is not available. The native alternative is not implemented."
#endif // __has_include
#endif // INCLUDE_STL_FS

#include "agghash.hpp"

//#include "types.h"


namespace daoc {

using std::move;
using std::string;
using std::vector;
using std::unordered_set;

// File Wrapping Types ---------------------------------------------------------
//! \brief Wrapper around the FILE* to prevent hanging file descriptors
class FileWrapper {
    FILE*  m_dsc;
    bool  m_tidy;
public:
    //! \brief Constructor
    //!
    //! \param fd FILE*  - the file descriptor to be held
    //! \param cleanup=true bool  - close the file descriptor on destruction
    //! 	(typically false if stdin/out is supplied)
    FileWrapper(FILE* fd=nullptr, bool cleanup=true) noexcept
    : m_dsc(fd), m_tidy(cleanup)  {}

    //! \brief Copy constructor
    //! \note Any file descriptor should have a single owner
    FileWrapper(const FileWrapper&)=delete;

	//! \brief Move constructor
	// ATTENTION: fw.m_dsc is not set to nullptr by the default move operation
	// ATTENTION: std::vector will move their elements if the elements' move constructor
	// is noexcept, and copy otherwise (unless the copy constructor is not accessible)
    FileWrapper(FileWrapper&& fw) noexcept
    : FileWrapper(fw.m_dsc, fw.m_tidy)
    {
    	fw.m_dsc = nullptr;
    }

    //! \brief Copy assignment
    //! \note Any file descriptor should have the single owner
    FileWrapper& operator= (const FileWrapper&)=delete;

	//! \brief Move assignment
	// ATTENTION: fw.m_dsc is not set to nullptr by the default move operation
    FileWrapper& operator= (FileWrapper&& fw) noexcept
    {
    	reset(fw.m_dsc, fw.m_tidy);
    	fw.m_dsc = nullptr;
    	return *this;
    }

    //! \brief Destructor
    ~FileWrapper()  // noexcept by default
    {
        if(m_dsc && m_tidy) {
            fclose(m_dsc);
            m_dsc = nullptr;
        }
    }

    //! \brief Implicit conversion to the file descriptor
    //!
    //! \return FILE*  - self as a file descriptor
    operator FILE*() const noexcept  { return m_dsc; }

    //! \brief Reset the wrapper
    //!
    //! \param fd FILE*  - the file descriptor to be held
    //! \param cleanup=true bool  - close the file descriptor on destruction
    //! 	(typically false if stdin/out is supplied)
    //! \return void
	void reset(FILE* fd=nullptr, bool cleanup=true) noexcept
	{
        if(m_dsc && m_tidy && m_dsc != fd)
            fclose(m_dsc);
    	m_dsc = fd;
    	m_tidy = cleanup;
	}

    //! \brief Release ownership of the holding file
    //!
    //! \return FILE*  - file descriptor
    FILE* release() noexcept
    {
    	auto fd = m_dsc;
    	m_dsc = nullptr;
		return fd;
    }
};

//! \brief Wrapper around the FILE* that holds also the filename giving ability
//! to reopen it and perform meaningful
// Note: we can't inherit from the FileWrapper because semantic of reset differs
class NamedFileWrapper {
	FileWrapper  m_file;  //!< File descriptor
	string  m_name;  //!< File name
public:
    //! \brief Default Constructor
    // Note: Required tor return empty objects using NRVO optimization
	NamedFileWrapper() noexcept: m_file(), m_name()  {}

    //! \brief Constructor
    //! \pre Parent directory must exists
    //!
    //! \param filename const char*  - new file name to be opened
    //! \param mode const char*  - opening mode, the same as fopen() has
	NamedFileWrapper(const char* filename, const char* mode)
	: m_file(filename && mode ? fopen(filename, mode) : nullptr)
	, m_name(filename ? filename : "")  {}

    //! \brief Copy constructor
    //! \note Any file descriptor should have a single owner
    NamedFileWrapper(const NamedFileWrapper&)=delete;

	//! \brief Move constructor
	// ATTENTION: std::vector will move their elements if the elements' move constructor
	// is noexcept, and copy otherwise (unless the copy constructor is not accessible)
    NamedFileWrapper(NamedFileWrapper&& fw) noexcept
    : m_file(move(fw.m_file)), m_name(move(fw.m_name))  {}

    //! \brief Copy assignment
    //! \note Any file descriptor should have the single owner
    NamedFileWrapper& operator= (const NamedFileWrapper&)=delete;

	//! \brief Move assignment
    NamedFileWrapper& operator= (NamedFileWrapper&& fw) noexcept
    {
		m_file = move(fw.m_file);
		m_name = move(fw.m_name);
		return *this;
    }

    //! \brief File name
    //!
    //! \return const string&  - file name
    const string& name() const noexcept  { return m_name; }

    //! \brief File size
    //!
    //! \return size_t  - file size or -1 on error
    size_t size() const noexcept;

    //! \brief Implicit conversion to the file descriptor
    //!
    //! \return FILE*  - file descriptor
    operator FILE*() const noexcept  { return m_file; }

    //! \brief Reopen the file under another mode
    //!
    //! \param mode const char*  - the mode of operations, the same as in fopen()
    //! \return NamedFileWrapper&  - the reopened file or closed (if can't be opened)
    NamedFileWrapper& reopen(const char* mode)
    {
		m_file.reset(freopen(nullptr, mode, m_file));  // m_name.c_str()
		return *this;
    }

    //! \brief Reset the file, closes current file and opens another one
    //! \pre Parent directory must exists
    //!
    //! \param filename const char*  - new file name to be opened
    //! \param mode const char*  - opening mode, the same as fopen() has
    //! \return NamedFileWrapper&  - the newly opened file or just the old one closed
	NamedFileWrapper& reset(const char* filename, const char* mode);

    //! \brief Release ownership of the holding file
    //!
    //! \return FILE*  - file descriptor
    FILE* release() noexcept  { return m_file.release(); }
};

// File Reading Types ----------------------------------------------------------
//! \brief Base of the StringBuffer
using StringBufferBase = vector<char>;

//! \brief String buffer to real file by lines using c-strings
//! \note The last symbol in the string is always set to 0 automatically
class StringBuffer: protected StringBufferBase {
	constexpr static size_t  spagesize = 4096;  // Small page size on x64

	size_t  m_cur;  //! Current position for the writing
	size_t  m_length;  //! Current length of the holding c-string
//protected:
//	StringBufferBase::size();
public:
    //! \brief
    //! \post the allocated buffer will have size >= 2
    //!
    //! \param size=spagesize size_t  - size of the buffer
    // Note: can throw bad_alloc
	StringBuffer(size_t size=spagesize);

    //! \brief Reset the string and it's shrink the allocated buffer
    //!
    //! \param size=spagesize size_t  - new initial size of the string buffer
    //! \return void
	void reset(size_t size=spagesize);

    //! \brief Length of the string including the terminating '\n' if present,
    //! 	but without the terminating '0'
    //!
    //! \return size_t  - length of the holding c-string without the null terminator
	size_t length() const noexcept  { return m_length; }

    //! \brief Whether the string is empty or starts with the newline symbol
    //! \attention empty() is true for '\n' when length() == 1
    //!
    //! \return bool  - the string is empty or starts with the '\n'
	bool empty() const
#if VALIDATE < 2
		noexcept
#endif // VALIDATE
	;

    //! \brief C-string including '\n' if it was present in the file
	operator char*() noexcept  { return data(); }

    //! \brief Const C-string including '\n' if it was present in the file
	operator const char*() const noexcept  { return data(); }

    //! \brief Make public indexing operators
	using StringBufferBase::operator[];
	using StringBufferBase::at;

    //! \brief Read line from the file and store including the terminating '\n' symbol
    //! \attention The read string contains the trailing '\n' if exist in the file
    //! \note The buffer might contain [part of] the read line on reading error
    //!
    //! \param input FILE*  - processing file
    //! \return bool  - whether the current line is read without any errors or
    //! the all lines already read (and the current one is empty)
	bool readline(FILE* input);
};

// File I/O functions declaration ----------------------------------------------
//! \brief Ensure existence of the specified directory
//!
//! \param dir const string&  - directory to be created if has not existed
//! \return void
void ensureDir(const string& dir);

//! \brief  Parse the header of CNL file and validate the results
//! \post clsnum <= ndsnum if ndsnum > 0. 0 means not specified
//!
//! \param fcls NamedFileWrapper&  - the reading file
//! \param line StringBuffer&  - processing line (string, header) being read from the file
//! \param[out] clsnum size_t&  - resulting number of clusters if specified, 0 in case of parsing errors
//! \param[out] ndsnum size_t&  - resulting number of nodes if specified, 0 in case of parsing errors
//! \param verbose=false bool  - print information about the header parsing issue to the stdout
//! \return void
void parseCnlHeader(NamedFileWrapper& fcls, StringBuffer& line, size_t& clsnum
	, size_t& ndsnum, bool verbose=false);

//! \brief Load all unique nodes from the CNL file with optional filtering by the cluster size
//!
//! \tparam Id  - Node id type
//! \tparam AccId  - Accumulated node ids type
//!
//! \param file NamedFileWrapper&  - input collection of clusters in the CNL format
//! \param membership=1 float  - expected membership of the nodes, >0, typically >= 1.
//! Used only for the node container preallocation to estimate the number of nodes
//! if not specified in the file header
//! \param ahash=nullptr AggHash<Id, AccId>*  - resulting aggregated hash of the loaded
//! node ids if not nullptr
//! \param cmin=0 size_t  - min allowed cluster size
//! \param cmax=0 size_t  - max allowed cluster size, 0 means any size
//! \param verbose=true bool  - print the number of loaded nodes to the stdout
//! \return bool  - the collection is loaded successfully
template <typename Id, typename AccId>
unordered_set<Id> loadNodes(NamedFileWrapper& file, float membership=1
	, AggHash<Id, AccId>* ahash=nullptr, size_t cmin=0, size_t cmax=0, bool verbose=true);

//! \brief Estimate the number of nodes from the CNL file size
//!
//! \param filesize size_t  - the number of bytes in the CNL file
//! \param membership=1.f float  - average membership of the node,
//! 	> 0, typically ~= 1
//! \return size_t  - estimated number of nodes
size_t estimateCnlNodes(size_t filesize, float membership=1.f) noexcept;

//! \brief Estimate the number of clusters from the number of nodes
//!
//! \param ndsnum size_t - the number of nodes
//! \param membership=1.f float  - average membership of the node,
//! 	> 0, typically ~= 1
//! \return size_t  - estimated number of clusters
size_t estimateClusters(size_t ndsnum, float membership=1.f) noexcept;

//! \brief Convert value to yes/no c-string
//!
//! \param val bool  - value to be converted
//! \return constexpr const char*  - resulting c-string
constexpr const char* toYesNo(bool val) noexcept  { return val ? "yes" : "no"; }

// File I/O templates definition -----------------------------------------------
template <typename Id, typename AccId>
unordered_set<Id> loadNodes(NamedFileWrapper& file, float membership
	, AggHash<Id, AccId>* ahash, size_t cmin, size_t cmax, bool verbose)
{
	unordered_set<Id>  nodebase;  // Node base;  Note: returned using NRVO optimization

	if(!file)
		return nodebase;

	// Note: CNL [CSN] format only is supported
	size_t  clsnum = 0;  // The number of clusters
	size_t  ndsnum = 0;  // The number of nodes

	// Note: strings defined out of the cycle to avoid reallocations
	StringBuffer  line;  // Reading line
	// Parse header and read the number of clusters if specified
	// Note: line includes terminating '\n'
	parseCnlHeader(file, line, clsnum, ndsnum, verbose);

	// Estimate the number of nodes in the file if not specified
	if(!ndsnum) {
		size_t  cmsbytes = file.size();
		if(cmsbytes != size_t(-1))  // File length fetching failed
			ndsnum = estimateCnlNodes(cmsbytes, membership);
		else if(clsnum)
			ndsnum = 2 * clsnum; // / membership;  // Note: use optimistic estimate instead of pessimistic (square / membership) to not overuse the memory
#if TRACE >= 2
		fprintf(stderr, "loadNodes(), estimated %lu nodes\n", ndsnum);
#endif // TRACE
	}
#if TRACE >= 2
	else fprintf(stderr, "loadNodes(), specified %lu nodes\n", ndsnum);
#endif // TRACE

	// Preallocate space for nodes
	if(ndsnum)
		nodebase.reserve(ndsnum);

	// Load clusters
	// ATTENTION: without '\n' delimiter the terminating '\n' is read as an item
	constexpr char  mbdelim[] = " \t\n";  // Delimiter for the members
	vector<Id>  cnds;  // Cluster nodes. Note: a dedicated container is required to filter clusters by size
	cnds.reserve(sqrt(ndsnum));  // Note: typically cluster size does not increase the square root of the number of nodes
#if TRACE >= 2
	size_t  totmbs = 0;  // The number of read member nodes from the file including repetitions
	size_t  fclsnum = 0;  // The number of read clusters from the file
#endif // TRACE
	do {
#if TRACE >= 3
		// Note: line includes terminating '\n'
		fprintf(stderr, "%lu> %s", fclsnum, static_cast<const char*>(line));
#endif // TRACE
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
				fprintf(stderr, "WARNING loadNodes(), empty cluster"
					" exists: '%s', skipped\n", cidstr);
				continue;
			}
		}
		do {
			// Note: only node id is parsed, share part is skipped if exists,
			// but potentially can be considered in NMI and F1 evaluation.
			// In the latter case abs diff of shares instead of co occurrence
			// counting should be performed.
			Id  nid = strtoul(tok, nullptr, 10);
#if VALIDATE >= 2
			if(!nid && tok[0] != '0') {
				fprintf(stderr, "WARNING loadNodes(), conversion error of '%s' into 0: %s\n"
					, tok, strerror(errno));
				continue;
			}
#endif // VALIDATE
#if TRACE >= 2
			++totmbs;  // Update the total number of read members
#endif // TRACE
			cnds.push_back(nid);
		} while((tok = strtok(nullptr, mbdelim)));
#if TRACE >= 2
		++fclsnum;  // The number of valid read lines, i.e. clusters
#endif // TRACE

		// Filter read cluster by size
		if(cnds.size() >= cmin && (!cmax || cnds.size() <= cmax))
			nodebase.insert(cnds.begin(), cnds.end());
		// Prepare outer vars for the next iteration
		cnds.clear();
	} while(line.readline(file));
//	// Rehash the nodes decreasing the allocated space if required
//	if(nodebase.size() <= nodebase.bucket_count() * nodebase.max_load_factor() / 3)
//		nodebase.reserve(nodebase.size());
#if TRACE >= 2
	printf("loadNodes(), the loaded base has %lu nodes from the input %lu members of %lu clusters\n"
		, nodebase.size(), totmbs, fclsnum);
#else
	if(verbose)
		printf("loadNodes(), nodebase nodes loaded: %lu\n", nodebase.size());
#endif // TRACE 2

	// Evaluate nodes hash if required
	if(ahash && nodebase.size()) {
		AggHash<Id, AccId>  ndsh;
		for(auto nid: nodebase)
			ndsh.add(nid);
		*ahash = move(ndsh);
	}

	return nodebase;
}

}  // daoc

#endif // FILEIO_H
