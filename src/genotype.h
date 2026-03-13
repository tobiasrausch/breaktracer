#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace breaktracer
{
  struct JunctionCount {
    std::vector<uint8_t> ref;
    std::vector<uint8_t> alt;
  };

  inline int32_t
  _editDistanceHW(std::string const& query, std::string const& target) {
    EdlibAlignResult align = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
    int32_t ed = align.editDistance;
    edlibFreeAlignResult(align);
    return ed;
  }
  
  inline int32_t
  _findSeqBp(bam1_t* rec, uint32_t const pos) {
    if (!rec->core.n_cigar) return -1;
    uint32_t rp = rec->core.pos; // reference pointer
    uint32_t sp = 0; // sequence pointer

    const uint32_t* cigar = bam_get_cigar(rec);
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
      uint32_t op = bam_cigar_op(cigar[i]);
      uint32_t oplen = bam_cigar_oplen(cigar[i]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	for (uint32_t k = 0; k < oplen; ++k, ++rp, ++sp) {
	  if (rp >= pos) return (int32_t)sp;
	}
      } else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	rp += oplen;
	if (rp >= pos) return (int32_t)sp;
      } else if (op == BAM_CINS) {
	sp += oplen;
      } else if (op == BAM_CSOFT_CLIP) {
	sp += oplen;
      } else if (op == BAM_CHARD_CLIP) {
      } else {
	std::cerr << "Unknown CIGAR operation" << std::endl;
      }
    }
    // Last aligned pos
    if (bam_cigar_op(cigar[rec->core.n_cigar - 1]) == BAM_CSOFT_CLIP) {
      return (int32_t)(sp - bam_cigar_oplen(cigar[rec->core.n_cigar - 1]));
    }
    return (int32_t)sp;
  }

  template<typename TConfig, typename TJunctionMap>
  inline void
  genotype(TConfig& c, std::vector<BrInTrace>& svs, TJunctionMap& jctMap) {
    if (svs.empty()) return;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> THeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    THeader hdr(c.files.size());
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
    }

    // Resize junction counts
    for (unsigned int file_c = 0; file_c < c.files.size(); ++file_c) jctMap[file_c].resize(svs.size());

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.hasDumpFile) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.dumpfile.string(), std::ios_base::out | std::ios_base::binary));
      dumpOut << "#svid\tbam\tqname\tchr\tpos\tmapq\ttype" << std::endl;
    }

    // Genotype SVs
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotyping" << std::endl;

    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr[0]->n_targets; ++refIndex) {
      // Fetch breakpoints
      typedef std::multimap<int32_t, int32_t> TBreakpointMap;
      TBreakpointMap bpMap;
      for(uint32_t i = 0; i < svs.size(); ++i) {
	if (svs[i].chr == refIndex) bpMap.insert(std::make_pair(svs[i].pos, i));
      }
      if (bpMap.empty()) continue;
      
      // Load sequence
      int32_t seqlen = -1;
      std::string tname(hdr[0]->target_name[refIndex]);
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr[0]->target_len[refIndex], &seqlen);

      // Iterate files
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[file_c]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

	  // Any SV breakpoint?
	  typedef std::set<int32_t> TSVSet;
	  TSVSet process;
	  int32_t rStart = rec->core.pos;
	  int32_t rEnd = bam_endpos(rec);
	  {
	    TBreakpointMap::const_iterator itBegin = bpMap.lower_bound(rStart);
	    TBreakpointMap::const_iterator itEnd = bpMap.upper_bound(rEnd);
	    for(; (itBegin != itEnd) && (itBegin != bpMap.end()); ++itBegin) process.insert(itBegin->second);
	  }

	  // Genotype SVs
	  std::string sequence;
	  for (typename TSVSet::const_iterator it = process.begin(); it != process.end(); ++it) {
	    int32_t svid = *it;
	    int32_t pos  = svs[svid].pos;
	    if ((svs[svid].chr != refIndex) || (pos < rStart) || (pos > rEnd)) continue;
	    if (svs[svid].consensus.empty()) continue;

	    // Load read sequence
	    if (sequence.empty()) {
	      sequence.resize(rec->core.l_qseq);
	      const uint8_t* seqptr = bam_get_seq(rec);
	      for (int ik = 0; ik < rec->core.l_qseq; ++ik) sequence[ik] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, ik)];
	    }
	    int32_t spBp = _findSeqBp(rec, pos);
	    if ((spBp < 0) || (spBp >= (int32_t)sequence.size())) continue;

	    // Adaptive flank
	    int32_t inslen = svs[svid].inslen;
	    int32_t maxFlank = std::min(500, std::max(inslen + 100, 200));
	    int32_t minFlank = std::max(30, std::min(inslen / 4 + 30, 100));
	    int32_t leftAvail = spBp;
	    int32_t rightAvail = (int32_t)sequence.size() - spBp;
	    if ((leftAvail < minFlank) || (rightAvail < minFlank)) continue;
	    int32_t leftFlank = std::min(leftAvail, maxFlank);
	    int32_t rightFlank = std::min(rightAvail, maxFlank);

	    // Probe
	    std::string probe = sequence.substr(spBp - leftFlank, leftFlank + rightFlank);

	    // Reference allele
	    int32_t tFlank   = 500;
	    int32_t refLeft  = std::max(0, pos - tFlank);
	    int32_t refRight = std::min(seqlen - 1, pos + tFlank);
	    std::string ref = boost::to_upper_copy(std::string(seq + refLeft, seq + refRight));

	    // Alt allele
	    std::string alt = boost::to_upper_copy(std::string(seq + refLeft, seq + pos)) + svs[svid].consensus + boost::to_upper_copy(std::string(seq + pos, seq + refRight));

	    // Align probe (query) to each allele (target)
	    int32_t edRef = _editDistanceHW(probe, ref);
	    int32_t edAlt = _editDistanceHW(probe, alt);
	    double normEdRef = (double)edRef / (double)probe.size();
	    double normEdAlt = (double)edAlt / (double)probe.size();	    
	    double maxErrRate = 0.15;
	    if ((normEdRef >= maxErrRate) && (normEdAlt >= maxErrRate)) continue;
	    if (normEdRef <= normEdAlt) jctMap[file_c][svid].ref.push_back(rec->core.qual);
	    else {
	      jctMap[file_c][svid].alt.push_back(rec->core.qual);
	      // Output insertion supporting reads
	      if (c.hasDumpFile) {
		std::string padNumber = std::to_string(svid);
		padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
		dumpOut << "BRINS" + padNumber << "\t" << c.files[file_c].string() << "\t" << bam_get_qname(rec) << "\t" << hdr[file_c]->target_name[rec->core.tid] << "\t" << rec->core.pos << "\t" << (int32_t)rec->core.qual << "\tSR" << std::endl;
	      }
	    }
	  }
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
      // Clean-up chromosome sequence
      if (seq != NULL) free(seq);
    }
    // Clean-up
    fai_destroy(fai);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      bam_hdr_destroy(hdr[file_c]);	  
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
    if (c.hasDumpFile) {
      dumpOut.pop();
      dumpOut.pop();
    }
  }

}

#endif
