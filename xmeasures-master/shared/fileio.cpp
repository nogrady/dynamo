//! \brief File IO utils
//!
//! \license Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0.html
//! > 	Simple explanation: https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)
//!
//! Copyright (c)
//! \authr Artem Lutov
//! \email luart@ya.ru
//! \date 2017-02-13

#include <cassert>
#include <system_error>  // error_code
//#include <stdexcept>

#ifdef __unix__
#include <sys/stat.h>
#endif // __unix__

#define INCLUDE_STL_FS
#include "fileio.hpp"


using std::error_code;
using std::to_string;
using fs::path;
using fs::create_directories;
using fs::is_directory;
using fs::exists;
using fs::status;
using std::logic_error;
using namespace daoc;

// File IO Types definitions ---------------------------------------------------
size_t NamedFileWrapper::size() const noexcept
{
	size_t  cmsbytes = -1;  // Return -1 on error
#ifdef __unix__  // sqrt(cmsbytes) lines => linebuf = max(4-8Kb, sqrt(cmsbytes) * 2) with dynamic realloc
	struct stat  filest;
	int fd = fileno(m_file);
	if(fd != -1 && !fstat(fd, &filest))
		return filest.st_size;
#endif // __unix
	error_code  err;
	cmsbytes = fs::file_size(m_name, err);
	if(cmsbytes == size_t(-1))
		fprintf(stderr, "WARNING size(), file size evaluation failed: %s\n"
			, err.message().c_str());

//	// Get length of the file
//	fseek(m_file, 0, SEEK_END);
//	cmsbytes = ftell(m_file);  // The number of bytes in the input communities
//	if(cmsbytes == size_t(-1))
//		perror("WARNING size(), file size evaluation failed");
//	//fprintf(stderr, "  %s: %lu bytes\n", fname, cmsbytes);
//	rewind(m_file);  // Set position to the begin of the file

	return cmsbytes;
}

NamedFileWrapper& NamedFileWrapper::reset(const char* filename, const char* mode)
{
	if(filename) {
		m_file.reset(fopen(filename, mode));
		m_name = filename;
	} else m_file.reset();
	return *this;
}

// File Reading Types ----------------------------------------------------------
StringBuffer::StringBuffer(size_t size)
: StringBufferBase(size), m_cur(0), m_length(0)
{
	if(size <= 2)
		size = 2;
	*data() = 0;  // Set first element to 0
	data()[size-2] = 0;  // Set prelast reserved element to 0
	// Note: data()[size-1] is set to 0 automatically on file read if
	// the reading data size >= size - 1 bytes
}

void StringBuffer::reset(size_t size)
{
	// Reset writing position
	m_cur = 0;
	m_length = 0;
	// Reset the buffer
	resize(size);  // Note: can throw bad_alloc
	shrink_to_fit();  // Free reserved memory
	*data() = 0;  // Set first element to 0
	data()[size-2] = 0;  // Set prelast reserved element to 0
	// Note: data()[size-1] is set to 0 automatically on file read if
	// the reading data size >= size - 1 bytes
}

//size_t StringBuffer::length() const
//#if VALIDATE < 2
//	noexcept
//#endif // VALIDATE
//{
//#if VALIDATE >= 2
//	const auto slen = strlen(data());
//	if(m_length != slen) {
//#if TRACE >= 2
//		fprintf(stderr, "length(), string: %s\n", data());
//#endif // TRACE
//		throw logic_error("ERROR length(), m_length (" + to_string(m_length)
//			+ ") != actual string length (" + to_string(slen) + ")\n");
//	}
//#endif // VALIDATE
//	return m_length;
//}

bool StringBuffer::empty() const
#if VALIDATE < 2
	noexcept
#endif // VALIDATE
{
#if VALIDATE >= 2
	if((!front() || front() == '\n') && m_length >= 2)
		throw logic_error("ERROR empty(), m_length (" + to_string(m_length)
			+ ") != actual string length (" + to_string(int(front() != 0)) + ")\n");
#endif // VALIDATE
	return !front() || front() == '\n';
}

bool StringBuffer::readline(FILE* input)
{
#if VALIDATE >= 2
	assert(input && !m_cur
		&& "readline(), valid file stream should be specified and have initial m_cur = 0");
#endif // VALIDATE
	*data() = 0;  // Set first element to 0 as an initialization to have the empty string on errors
	const auto ibeg = ftell(input);
	// Read data from file until the string is read or an error occurs
	while(fgets(data() + m_cur, size() - m_cur, input) && data()[size()-2]) {
#if TRACE >= 3  // Verified
		fprintf(stderr, "readline(), resizing buffer of %lu bytes, %lu pos: %s\n"
			, size(), m_cur, data());
#endif // TRACE
		m_cur = size() - 1;  // Start overwriting ending '0' of the string
		resize(size() + (size() / (spagesize * 2) + 1) * spagesize);
		data()[size() - 2] = 0;  // Set prelast element to 0
	}
	const auto iend = ftell(input);
#if VALIDATE >= 2
	if(iend == -1 || ibeg == -1)
		perror("ERROR, file position reading error");
	const size_t  slen = strlen(data());
	if(!((!m_cur || slen >= m_cur) && size_t(iend - ibeg) == slen)) {
		fprintf(stderr, "readline(), m_cur: %lu, slen: %lu, dpos: %li,  str: %s\n"
			, m_cur, slen, iend - ibeg, data());
		assert(0 && "readline(), string size validation failed");
	}
#endif // VALIDATE
	m_cur = 0;  // Reset the writing (appending) position
	// Note: prelast and last elements of the buffer will be always zero

	// Set string length
	m_length = iend != -1 && ibeg != -1 ? iend - ibeg : strlen(data());

	// Check for errors
	if((!m_length && feof(input)) || ferror(input)) {
		if(ferror(input))
			perror("ERROR readline(), file reading error");
		return false;  // No more lines can be read
	}

	return true;  // More lines can be read
}

// File I/O functions ----------------------------------------------------------
namespace daoc {

void ensureDir(const string& dir)
{
#if TRACE >= 3
	fprintf(stderr, "ensureDir(), ensuring existence of: %s\n", dir.c_str());
#endif // TRACE
	// Check whether the output directory exists and create it otherwise
	path  outdir = dir;
	if(!exists(outdir)) {
		error_code  err;
		if(!create_directories(outdir, err))
//			fputs(string("ERROR ensureDir(), target directory '").append(dir)
//				.append("' can't be created: ").append(err.message())
//				.append("\n").c_str(), stderr);
			throw std::ios_base::failure(string("ERROR ensureDir(), target directory '")
				.append(dir).append("' can't be created: ") += err.message());
	} else if(!is_directory(outdir))
//		fputs(string("ERROR ensureDir(), target entry '").append(dir)
//			.append("' already exists as a non-directory path\n").c_str(), stderr);
		throw std::ios_base::failure(string("ERROR ensureDir(), target entry '").append(dir)
			+= "' already exists as a non-directory path\n");
}

void parseCnlHeader(NamedFileWrapper& fcls, StringBuffer& line, size_t& clsnum
	, size_t& ndsnum, bool verbose)
{
    //! Parse count value
    //! \return  - id value of 0 in case of parsing errors
	auto parseCount = []() noexcept -> size_t {
		char* tok = strtok(nullptr, " \t,");  // Note: the value can't be ended with ':'
		//errno = 0;
		const auto val = strtoul(tok, nullptr, 10);
		if(errno)
			perror(string("WARNING parseCount(), id value parsing error for the tok '")
				.append(tok).append("'").c_str());
		return val;
	};

	errno = 0;  // Reset errno
	// Process the header, which is a special initial comment
	// The target header is:  # Clusters: <cls_num>[,] Nodes: <cls_num>
	constexpr char  clsmark[] = "clusters";
	constexpr char  ndsmark[] = "nodes";
	constexpr char  attrnameDelim[] = " \t:,";
#if TRACE >= 2
	size_t  lnum = 0;  // The number of lines read
#endif // TRACE
	while(line.readline(fcls)) {
#if TRACE >= 2
		++lnum;
#endif // TRACE
		// Skip empty lines
		if(line.empty())
			continue;
		// Consider only subsequent comments
		if(line[0] != '#')
			break;

		// Tokenize the line
		char *tok = strtok(line + 1, attrnameDelim);  // Note: +1 to skip the leading '#'
		// Skip comment without the string continuation and continuous comment
		if(!tok || tok[0] == '#')
			continue;
		uint8_t  attrs = 0;  // The number of read attributes
		do {
			// Lowercase the token
			for(char* pos = tok; *pos; ++pos)
				*pos = tolower(*pos);

			// Identify the attribute and read it's value
			if(!strcmp(tok, clsmark)) {
				clsnum = parseCount();
				++attrs;
#if TRACE >= 2
				fprintf(stderr, "parseCnlHeader(), clusters: %lu\n", clsnum);
#endif // TRACE
			} else if(!strcmp(tok, ndsmark)) {
				ndsnum = parseCount();
				++attrs;
#if TRACE >= 2
				fprintf(stderr, "parseCnlHeader(), nodes: %lu\n", ndsnum);
#endif // TRACE
			} else {
#if TRACE >= 1
#if TRACE < 2
			if(verbose)
#endif // TRACE 2
				fprintf(
#if TRACE >= 2
				stderr
#else
				stdout
#endif // TRACE 2
				, "WARNING parseCnlHeader(), the header parsing is omitted"
					" because of the unexpected attribute: %s\n", tok);
#endif // TRACE 1
				break;
			}
		} while((tok = strtok(nullptr, attrnameDelim)) && attrs < 2);

		// Validate and correct the number of clusters if required
		// Note: it's better to reallocate a container a few times than too much overconsume the memory
		if(ndsnum && clsnum > ndsnum) {
			fprintf(stderr, "WARNING parseCnlHeader(), clsnum (%lu) typically should be"
				" less than ndsnum (%lu)\n", clsnum, ndsnum);
			clsnum = ndsnum;
			//assert(0 && "parseCnlHeader(), clsnum typically should be less than ndsnum");
		}
		// Get following line for the unified subsequent processing
		line.readline(fcls);
		break;
	}
#if TRACE >= 2
	fprintf(stderr, "parseCnlHeader(), processed %lu lines of '%s'\n"
		, lnum, fcls.name().c_str());
#endif // TRACE
}

size_t estimateCnlNodes(size_t filesize, float membership) noexcept
{
	if(membership <= 0) {
		fprintf(stderr, "WARNING estimateCnlNodes(), invalid membership = %G specified"
			", reseted to 1\n", membership);
		membership = 1;
		//throw invalid_argument("estimateCnlNodes(), membership = "
		//	+ to_string(membership) + " should be positive\n");
	}

	size_t  ndsnum = 0;  // The estimated number of nodes
	if(filesize) {
		size_t  magn = 10;  // Decimal ids magnitude
		unsigned  img = 1;  // Index of the magnitude (10^1)
		size_t  reminder = filesize % magn;  // Reminder in bytes
		ndsnum = reminder / ++img;  //  img digits + 1 delimiter for each element
		while(filesize >= magn) {
			magn *= 10;
			ndsnum += (filesize - reminder) % magn / ++img;
			reminder = filesize % magn;
		}
	}
	return ndsnum / membership;
}

size_t estimateClusters(size_t ndsnum, float membership) noexcept
{
	if(membership <= 0) {
		fprintf(stderr, "WARNING estimateClusters(), invalid membership = %G specified"
			", reseted to 1\n", membership);
		membership = 1;
		//throw invalid_argument("estimateClusters(), membership = "
		//	+ to_string(membership) + " should be positive\n");
	}

	size_t  clsnum = 0;  // The estimated number of clusters
	// Usually the number of clusters does not increase square root of the number of nodes
	// Note: do not estimate in case the number of nodes is not specified
	if(ndsnum)
		clsnum = sqrt(ndsnum * membership) + 1;  // Note: +1 to consider rounding down
	return clsnum;
}

}  // daoc
