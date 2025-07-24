#ifndef UTIL_H
#define UTIL_H

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>

#include <htslib/sam.h>

#include "edlib.h"

namespace breaktracer
{

  inline uint32_t
  infixStart(EdlibAlignResult const& cigar) {
    int32_t tIdx = cigar.endLocations[0];
    for (int32_t i = 0; i < cigar.alignmentLength; i++) {
      if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
    }
    if (tIdx >= 0) return tIdx + 1;
    else return 0;
  }

  inline uint32_t
  infixEnd(EdlibAlignResult const& cigar) {
    return cigar.endLocations[0];
  }
  
  inline void
  printAlignmentPretty(std::string const& query, std::string const& target, EdlibAlignMode const modeCode, EdlibAlignResult const& align) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = align.endLocations[0];
        for (int32_t i = 0; i < align.alignmentLength; i++) {
            if (align.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
        }
    }
    std::cerr << std::endl;
    for (int start = 0; start < align.alignmentLength; start += 50) {
      std::cerr << "T: ";
      int32_t startTIdx = -1;
      for (int32_t j = start; ((j < start + 50) && (j < align.alignmentLength)); ++j) {
	if (align.alignment[j] == EDLIB_EDOP_INSERT) std::cerr << "-";
	else std::cerr << target[++tIdx];
	if (j == start) startTIdx = tIdx;
      }
      std::cerr << " (" << std::max(startTIdx, 0) << " - " << tIdx << ")" << std::endl;

      // match / mismatch
      std::cerr << ("   ");
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_MATCH) std::cerr <<  "|";
	else std::cerr << " ";
      }
      std::cerr << std::endl;

      // query
      std::cerr << "Q: ";
      int32_t startQIdx = qIdx;
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_DELETE) std::cerr << "-";
	else std::cerr << query[++qIdx];
	if (j == start) startQIdx = qIdx;
      }
      std::cerr << " ("<< std::max(startQIdx, 0) << " - " << qIdx << ")" << std::endl;
      std::cerr << std::endl;
    }
  }

  inline void
  printAlignment(std::string const& seqI, std::string const& seqJ, EdlibAlignMode const modeCode, EdlibAlignResult const& cigar) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) seqJ.size()) missingEnd = seqJ.size() - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) std::cerr << '-';
      }
    }
    // seqI
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) std::cerr << '-';
      else std::cerr << seqI[++qIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) std::cerr << '-';
      }
    }
    std::cerr << std::endl;
    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) std::cerr << seqJ[j];
      }
    }
    // seqJ
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) std::cerr << '-';
      else std::cerr << seqJ[++tIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) std::cerr << seqJ[++tIdx];
      }
    }
    std::cerr << std::endl;
  }

  
  template<typename TConfig>
  inline void
  checkSampleNames(TConfig& c) {
    uint32_t ucount = 0;
    std::set<std::string> snames;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      while (snames.find(c.sampleName[file_c]) != snames.end()) {
	std::cerr << "Warning: Duplicate sample names: " << c.sampleName[file_c] << std::endl;
	c.sampleName[file_c] += "_" + boost::lexical_cast<std::string>(ucount++);
	std::cerr << "Warning: Changing sample name to " << c.sampleName[file_c] << std::endl;
      }
      snames.insert(c.sampleName[file_c]);
    }
  }

  inline bool
  is_gz(boost::filesystem::path const& f) {
    std::ifstream bfile(f.string().c_str(), std::ios_base::binary | std::ios::ate);
    bfile.seekg(0, std::ios::beg);
    char byte1;
    bfile.read(&byte1, 1);
    char byte2;
    bfile.read(&byte2, 1);
    bfile.close();
    if ((byte1 == '\x1F') && (byte2 == '\x8B')) return true;
    else return false;
  }

  // Output directory/file checks
  inline bool
  _outfileValid(boost::filesystem::path const& outfile) {
    try {
      boost::filesystem::path outdir;
      if (outfile.has_parent_path()) outdir = outfile.parent_path();
      else outdir = boost::filesystem::current_path();
      if (!boost::filesystem::exists(outdir)) {
	std::cerr << "Output directory does not exist: " << outdir << std::endl;
	return false;
      } else {
	boost::filesystem::file_status s = boost::filesystem::status(outdir);
	boost::filesystem::ofstream file(outfile.string());
	file.close();
	if (!(boost::filesystem::exists(outfile) && boost::filesystem::is_regular_file(outfile))) {
	  std::cerr << "Fail to open output file " << outfile.string() << std::endl;
	  std::cerr << "Output directory permissions: " << s.permissions() << std::endl;
	  return false;
	} else {
	  boost::filesystem::remove(outfile.string());
	}
      }
    } catch (boost::filesystem::filesystem_error const& e) {
      std::cerr << e.what() << std::endl;
      return false;
    }
    return true;
  }


  inline uint32_t sequenceLength(bam1_t const* rec) {
    const uint32_t* cigar = bam_get_cigar(rec);
    uint32_t slen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CINS) || (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) slen += bam_cigar_oplen(cigar[i]);
    return slen;
  }

  inline int32_t
  readLength(bam1_t const* rec) {
    //int32_t slen = rec->core.l_qseq;  # Incorrect for seq. with hard-clips
    return sequenceLength(rec);
  }
    
  inline uint32_t alignmentLength(bam1_t const* rec) {
    const uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }
  
  inline std::size_t hash_lr(bam1_t* rec) {
    boost::hash<std::string> string_hash;
    std::string qname = bam_get_qname(rec);
    std::size_t seed = hash_string(qname.c_str());
    boost::hash_combine(seed, string_hash(qname));
    return seed;
  }

  inline std::size_t hash_lr(std::string const& qname) {
    boost::hash<std::string> string_hash;
    std::size_t seed = hash_string(qname.c_str());
    boost::hash_combine(seed, string_hash(qname));
    return seed;
  }
  
  inline void
  reverseComplement(std::string& sequence) {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }

  inline void
  getSMTag(std::string const& header, std::string const& fileName, std::string& sampleName) {
    std::set<std::string> smIdentifiers;
    std::string delimiters("\n");
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of(delimiters));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool rgPresent = false;
    for(;itH!=itHEnd; ++itH) {
      if (itH->find("@RG")==0) {
	std::string delim("\t");
	TStrParts keyval;
	boost::split(keyval, *itH, boost::is_any_of(delim));
	TStrParts::const_iterator itKV = keyval.begin();
	TStrParts::const_iterator itKVEnd = keyval.end();
	for(;itKV != itKVEnd; ++itKV) {
	  size_t sp = itKV->find(":");
	  if (sp != std::string::npos) {
	    std::string field = itKV->substr(0, sp);
	    if (field == "SM") {
	      rgPresent = true;
	      std::string rgSM = itKV->substr(sp+1);
	      smIdentifiers.insert(rgSM);
	    }
	  }
	}
      }
    }
    if (!rgPresent) {
      sampleName = fileName;
    } else if (smIdentifiers.size() == 1) {
      sampleName = *(smIdentifiers.begin());
    } else if (smIdentifiers.size() > 1) {
      sampleName = *(smIdentifiers.begin());
      std::cerr << "Warning: Multiple sample names (@RG:SM) present in the BAM file!" << std::endl;
    }
  }

}

#endif
