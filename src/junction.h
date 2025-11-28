#ifndef JUNCTION_H
#define JUNCTION_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/dynamic_bitset.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace breaktracer
{

  struct TraceCandidate {
    uint32_t kLow;
    uint32_t kHigh;
    int32_t sStart;
    int32_t sEnd;
    int32_t len;
    int32_t score;
    double pid;

    TraceCandidate(uint32_t const kL, uint32_t const kH, int32_t const sSt, int32_t const sEn, int32_t const ln, int32_t const s, double p) : kLow(kL), kHigh(kH), sStart(sSt), sEnd(sEn), len(ln), score(s), pid(p) {}
    
    bool operator<(const TraceCandidate& other) const {
      return score > other.score;
    }
  };

  template<typename TReadBp>
  inline void
    _insertJunction(TReadBp& readBp, std::size_t const seed, bam1_t* rec, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = true;
    if (rec->core.flag & BAM_FREVERSE) fw = false;
    typedef typename TReadBp::mapped_type TJunctionVector;
    typename TReadBp::iterator it = readBp.find(seed);
    int32_t seqlen = readLength(rec);
    if (sp <= seqlen) {
      if (rec->core.flag & BAM_FREVERSE) {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, rp, seqlen - sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, rp, seqlen - sp, rec->core.qual))));
      } else {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, rp, sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, rp, sp, rec->core.qual))));
      }
    }
  }

  template<typename TConfig, typename TReadBp>
  inline void
  findJunctions(TConfig const& c, TReadBp& readBp) {
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Split-read scanning" << std::endl;

    // Load optional genome mask
    faidx_t* faiMap = NULL;
    if (c.hasGenomeMask) faiMap = fai_load(c.mask.string().c_str());
    
    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      // Load optional genome mask
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet masked;
      char* seq = NULL;
      if (c.hasGenomeMask) {
	masked.resize(hdr->target_len[refIndex], false);
	std::string tname(hdr->target_name[refIndex]);
	int32_t seqlen = faidx_seq_len(faiMap, tname.c_str());
	if (seqlen != - 1) {
	  int32_t seqout = -1;
	  seq = faidx_fetch_seq(faiMap, tname.c_str(), 0, seqlen, &seqout);
	  for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	    if (seq[i] == 'N') masked[i] = true;
	  }
	}
      }

      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Read alignments
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {

	  // Keep secondary alignments
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	  
	  std::size_t seed = hash_lr(rec);
	  //std::cerr << bam_get_qname(rec) << '\t' << seed << std::endl;
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	  
	  // Parse the CIGAR
	  const uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      sp += bam_cigar_oplen(cigar[i]);
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) {
		if ((!c.hasGenomeMask) || (!masked[rp])) _insertJunction(readBp, seed, rec, rp, sp, false);
	      }
	      rp += bam_cigar_oplen(cigar[i]);
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) { // Try look-ahead
		uint32_t spOrig = sp;
		uint32_t rpTmp = rp;
		uint32_t spTmp = sp;
		uint32_t dlen = bam_cigar_oplen(cigar[i]);
		for (std::size_t j = i + 1; j < rec->core.n_cigar; ++j) {
		  if ((bam_cigar_op(cigar[j]) == BAM_CMATCH) || (bam_cigar_op(cigar[j]) == BAM_CEQUAL) || (bam_cigar_op(cigar[j]) == BAM_CDIFF)) {
		    spTmp += bam_cigar_oplen(cigar[j]);
		    rpTmp += bam_cigar_oplen(cigar[j]);
		    if ((double) (spTmp - sp) / (double) (dlen + (rpTmp - rp)) > c.indelExtension) break;
		  } else if (bam_cigar_op(cigar[j]) == BAM_CDEL) {
		    rpTmp += bam_cigar_oplen(cigar[j]);
		    if (bam_cigar_oplen(cigar[j]) > c.minRefSep) {
		      // Extend deletion
		      dlen += (rpTmp - rp);
		      rp = rpTmp;
		      sp = spTmp;
		      i = j;
		    }
		  } else if (bam_cigar_op(cigar[j]) == BAM_CINS) {
		    if (bam_cigar_oplen(cigar[j]) > c.minRefSep) break; // No extension
		    spTmp += bam_cigar_oplen(cigar[j]);
		  } else break; // No extension
		}
		if ((!c.hasGenomeMask) || (!masked[rp])) _insertJunction(readBp, seed, rec, rp, spOrig, true);
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) {
		if ((!c.hasGenomeMask) || (!masked[rp])) _insertJunction(readBp, seed, rec, rp, sp, false);
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) { // Try look-ahead
		uint32_t rpOrig = rp;
		uint32_t rpTmp = rp;
		uint32_t spTmp = sp;
		uint32_t ilen = bam_cigar_oplen(cigar[i]);
		for (std::size_t j = i + 1; j < rec->core.n_cigar; ++j) {
		  if ((bam_cigar_op(cigar[j]) == BAM_CMATCH) || (bam_cigar_op(cigar[j]) == BAM_CEQUAL) || (bam_cigar_op(cigar[j]) == BAM_CDIFF)) {
		    spTmp += bam_cigar_oplen(cigar[j]);
		    rpTmp += bam_cigar_oplen(cigar[j]);
		    if ((double) (rpTmp - rp) / (double) (ilen + (spTmp - sp)) > c.indelExtension) break;
		  } else if (bam_cigar_op(cigar[j]) == BAM_CDEL) {
		    if (bam_cigar_oplen(cigar[j]) > c.minRefSep) break; // No extension
		    rpTmp += bam_cigar_oplen(cigar[j]);
		  } else if (bam_cigar_op(cigar[j]) == BAM_CINS) {
		    spTmp += bam_cigar_oplen(cigar[j]);
		    if (bam_cigar_oplen(cigar[j]) > c.minRefSep) {
		      // Extend insertion
		      ilen += (spTmp - sp);
		      rp = rpTmp;
		      sp = spTmp;
		      i = j;
		    }
		  } else {
		    break; // No extension
		  }
		}
		if ((!c.hasGenomeMask) || (!masked[rp])) _insertJunction(readBp, seed, rec, rpOrig, sp, true);
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	      int32_t finalsp = sp;
	      bool scleft = false;
	      if (sp == 0) {
		finalsp += bam_cigar_oplen(cigar[i]); // Leading soft-clip / hard-clip
		scleft = true;
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      //std::cerr << bam_get_qname(rec) << ',' << rp << ',' << finalsp << ',' << scleft << std::endl;
	      if (bam_cigar_oplen(cigar[i]) > c.minClip) {
		if ((!c.hasGenomeMask) || (!masked[rp])) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
	      }
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
      if (c.hasGenomeMask) {
	if (seq != NULL) free(seq);
      }
    }

    // Sort junctions
    for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) std::sort(it->second.begin(), it->second.end());

    // Clean-up
    if (c.hasGenomeMask) fai_destroy(faiMap);
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }


  template<typename TConfig, typename TReadBp>
  inline void
  process_read(TConfig const& c, std::string const& searchseq, TReadBp& readBp, bam1_t* rec, std::vector<TraceRecord>& traces) {
    int32_t maxFragSize = searchseq.size() + 0.15 * searchseq.size();
    std::size_t seed = hash_lr(rec);
    std::string sequence;
    sequence.resize(rec->core.l_qseq);
    const uint8_t* seqptr = bam_get_seq(rec);
    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);

    // Get candidates
    std::vector<TraceCandidate> candidates;
    for(uint32_t i = 0; i < readBp[seed].size(); ++i) {
      for(uint32_t k = i + 1; k < readBp[seed].size(); ++k) {
	// Breakpoint coords on read sequence
	int32_t sStart = readBp[seed][i].seqpos;
	int32_t sEnd = readBp[seed][k].seqpos;
	int32_t fragsize = sEnd - sStart;
	if (sEnd < (c.minSeedAlign + c.cropSize)) continue; 
	if (sStart + (c.minSeedAlign + c.cropSize) > rec->core.l_qseq) continue;
	if ((fragsize > (c.minSeedAlign + c.cropSize)) && (fragsize < maxFragSize)) {
	  double pid = 0;
	  bool candidate = false;
	  
	  // Full-length approach
	  std::string fullseq = sequence.substr(sStart, fragsize);
	  EdlibAlignResult cigarFull = edlibAlign(fullseq.c_str(), fullseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	  double pIdFull = 1.0 - ( (double) (cigarFull.editDistance) / (double) (fragsize) );
	  edlibFreeAlignResult(cigarFull);
	  
	  if (pIdFull > c.pctThres) {
	    candidate = true;
	    pid = pIdFull;
	  } else {
	    // Partial approach (scan seeds)
	    int32_t validSeeds = 0;
	    int32_t nCount = 0;
	    double percentIdentity = 0;
	    for(int32_t sCoord = sStart + c.cropSize; sCoord + c.minSeedAlign <= (sEnd - c.cropSize); sCoord += c.minSeedAlign) {
	      std::string subseq = sequence.substr(sCoord, c.minSeedAlign);
	      EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	      double pIdSeed = 1.0 - ( (double) (cigar.editDistance) / (double) (c.minSeedAlign) );
	      edlibFreeAlignResult(cigar);
	      percentIdentity += pIdSeed;
	      if (pIdSeed > c.pctThres) ++validSeeds;
	      ++nCount;
	    }
	    if (validSeeds > 0) {
	      percentIdentity /= (double) nCount;
	      if (percentIdentity > c.pctThres) {
		candidate = true;
		pid = percentIdentity;
	      }
	    }
	  }
	  if (candidate) candidates.push_back(TraceCandidate(i, k, sStart, sEnd, fragsize, (int32_t) (fragsize * pid), pid));
	}
      }
    }

    // Sort by score
    std::sort(candidates.begin(), candidates.end());

    // Flag forced mappings (if sequence is present in reference)
    typedef std::pair<int32_t, int32_t> TRefPos;
    std::set<TRefPos> filterRefPos;
    for(uint32_t i = 0; i < candidates.size(); ++i) {
      if (readBp[seed][candidates[i].kLow].refidx == readBp[seed][candidates[i].kHigh].refidx) {
	int32_t reflen = std::abs(readBp[seed][candidates[i].kHigh].refpos - readBp[seed][candidates[i].kLow].refpos);
	int32_t seqlen = candidates[i].len;
	double frac;
	if (reflen < seqlen) frac = (double) reflen / (double) seqlen;
	else frac = (double) seqlen / (double) reflen;
	if (frac > 0.8) {
	  filterRefPos.insert(std::make_pair(readBp[seed][candidates[i].kLow].refidx, readBp[seed][candidates[i].kLow].refpos));
	  filterRefPos.insert(std::make_pair(readBp[seed][candidates[i].kHigh].refidx, readBp[seed][candidates[i].kHigh].refpos));
	}
      }
    }

    // Process candidates
    std::vector<std::pair<int32_t, int32_t>> selectedIntervals;
    for(const auto& cand : candidates) {
      // check if one of the breakpoints got flagged
      bool flagged = false;
      if (filterRefPos.find(std::make_pair(readBp[seed][cand.kLow].refidx, readBp[seed][cand.kLow].refpos)) != filterRefPos.end()) flagged = true;
      if (filterRefPos.find(std::make_pair(readBp[seed][cand.kHigh].refidx, readBp[seed][cand.kHigh].refpos)) != filterRefPos.end()) flagged = true;
      if (flagged) continue;

      // Debug
      //std::cerr << bam_get_qname(rec) << '\t' << readBp[seed][cand.kLow].refidx << ':' << readBp[seed][cand.kLow].refpos << '\t' << readBp[seed][cand.kHigh].refidx << ':' << readBp[seed][cand.kHigh].refpos << '\t' << cand.len << ',' << cand.pid << ',' << cand.score << '\t' << cand.sStart << ',' << cand.sEnd << '\t' << cand.kLow << ',' << cand.kHigh << '\t' << flagged << std::endl;      

      // check if it overlaps a processed interval
      bool overlap = false;
      for(const auto& interval : selectedIntervals) {
	if (std::max(cand.sStart, interval.first) < std::min(cand.sEnd, interval.second)) {
	  overlap = true;
	  break;
	}
      }
      if (overlap) continue;

      // check for clean left flank
      bool validFlanks = true;
      int32_t leftAvail = cand.sStart;
      if (leftAvail > 0) {
	int32_t checkLen = std::min((int32_t) leftAvail, (int32_t) maxFragSize); 
	if (checkLen > 50) { 
	  std::string leftSeq = sequence.substr(cand.sStart - checkLen, checkLen);
	  EdlibAlignResult cigar = edlibAlign(leftSeq.c_str(), leftSeq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	  double pId = 1.0 - ( (double) (cigar.editDistance) / (double) (leftSeq.size()) );
	  edlibFreeAlignResult(cigar);
	  if (pId > c.pctThres) validFlanks = false;
	} else validFlanks = false;
      }

      // right flank
      int32_t rightAvail = sequence.size() - cand.sEnd;
      if ((validFlanks) && (rightAvail > 0)) {
	int32_t checkLen = std::min((int32_t) rightAvail, (int32_t) maxFragSize);
	if (checkLen > 50) {
	  std::string rightSeq = sequence.substr(cand.sEnd, checkLen);
	  EdlibAlignResult cigar = edlibAlign(rightSeq.c_str(), rightSeq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	  double pId = 1.0 - ( (double) (cigar.editDistance) / (double) (rightSeq.size()) );
	  edlibFreeAlignResult(cigar);
	  if (pId > c.pctThres) validFlanks = false;
	} else validFlanks = false;
      }

      if (validFlanks) {
	TraceRecord tr;
	// Canonical ordering
	if ((readBp[seed][cand.kLow].refidx < readBp[seed][cand.kHigh].refidx) || ( ( readBp[seed][cand.kLow].refidx == readBp[seed][cand.kHigh].refidx ) && (readBp[seed][cand.kLow].refpos <= readBp[seed][cand.kHigh].refpos) ) ) {
	  tr = TraceRecord(readBp[seed][cand.kLow].refidx, readBp[seed][cand.kLow].refpos, readBp[seed][cand.kLow].seqpos, readBp[seed][cand.kHigh].refidx, readBp[seed][cand.kHigh].refpos, readBp[seed][cand.kHigh].seqpos, (int32_t) (cand.pid * 100), cand.len, seed);
	} else {
	  tr = TraceRecord(readBp[seed][cand.kHigh].refidx, readBp[seed][cand.kHigh].refpos, readBp[seed][cand.kHigh].seqpos, readBp[seed][cand.kLow].refidx, readBp[seed][cand.kLow].refpos, readBp[seed][cand.kLow].seqpos, (int32_t) (cand.pid * 100), cand.len, seed);
	}
	traces.push_back(tr);
	selectedIntervals.push_back(std::make_pair(cand.sStart, cand.sEnd));
      }
    }
  }
  

  template<typename TConfig, typename TReadBp>
  inline void
  findBrIn(TConfig const& c, TReadBp& readBp, std::vector<TraceRecord>& tr) {
    std::string searchseq;
    if (c.insmode == 0) {
      std::string faname = "";
      if (!loadSingleFasta(c, faname, searchseq)) return;
    } else {
      if (c.insmode == 1) searchseq = MEI::alu;
      else if (c.insmode == 3) searchseq = MEI::sva;
      else if (c.insmode == 4) searchseq = MEI::numt; 
      else searchseq = MEI::line1;
      // Add poly-A for MEIs
      if (c.insmode != 4) searchseq += MEI::polyA;
    }

    // Augment with reverse complement
    std::string revseq = searchseq;
    reverseComplement(revseq);
    searchseq += revseq;

    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      hts_set_fai_filename(samfile[file_c], c.genome.string().c_str());
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    
    // Parse genome chr-by-chr
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Insertion search" << std::endl;

    // Insertion candidate read search
    std::set<std::size_t> clusteredReads;
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      typedef std::vector<uint16_t> TCounter;
      TCounter bp(hdr->target_len[refIndex], 0);
      uint8_t maxVal = std::numeric_limits<uint8_t>::max();

      // Cluster breakpoints +/-5bp
      for(uint32_t z = 0; z < 2; ++z) {
	for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	  if (it->second.size() > 1) {
	    for(uint32_t i = 0; i < it->second.size(); ++i) {
	      if (it->second[i].refidx == refIndex) {
		for(int32_t k = std::max(0, it->second[i].refpos - 5); k < std::min((int32_t) hdr->target_len[refIndex], it->second[i].refpos + 5); ++k) {
		  if (z) {
		    // 2nd pass
		    if (bp[k] >= c.minCliqueSize) clusteredReads.insert(it->first);
		  } else {
		    // 1st pass
		    if (bp[k] < maxVal) ++bp[k];
		  }
		}
	      }
	    }
	  }
	}
      }
      //std::cerr << refIndex << ',' << clusteredReads.size() << std::endl;
    }

    // Find insertion
    if (!clusteredReads.empty()) {

      // Threads
      ThreadPool pool(std::max<std::size_t>(1, c.maxThreads));
      std::vector<std::future<void>> futures;
      std::mutex tr_mutex;
      std::vector<bam1_t*> batch;
      
      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
	  // Read alignments
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    // Keep only primary alignments
	    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	    std::size_t seed = hash_lr(rec);
	    if (clusteredReads.find(seed) != clusteredReads.end()) {
	      if (readBp[seed].size()>1) {
		batch.push_back(bam_dup1(rec));
		if (batch.size() == c.batchSize) {
		  futures.push_back(pool.enqueue([&, batch] {
		    std::vector<TraceRecord> local_tr;
		    for(auto* srd : batch) {
		      std::vector<TraceRecord> singleTRs;
		      process_read(c, searchseq, readBp, srd, singleTRs);
		      local_tr.insert(local_tr.end(), singleTRs.begin(), singleTRs.end());
		    }
		    // Merge into shared vector
		    std::lock_guard<std::mutex> lock(tr_mutex);
		    tr.insert(tr.end(), local_tr.begin(), local_tr.end());
		    // Clean-up
		    for(auto* srd : batch) bam_destroy1(srd);
		  }));
		  batch.clear();
		}
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }
      if (!batch.empty()) {
	futures.push_back(pool.enqueue([&, batch] {
	  std::vector<TraceRecord> local_tr;
	  for(auto* srd : batch) {
	    std::vector<TraceRecord> singleTRs;
	    process_read(c, searchseq, readBp, srd, singleTRs);
	    local_tr.insert(local_tr.end(), singleTRs.begin(), singleTRs.end());
	  }
	  // Merge into shared vector
	  std::lock_guard<std::mutex> lock(tr_mutex);
	  tr.insert(tr.end(), local_tr.begin(), local_tr.end());
	  // Clean-up
	  for(auto* srd : batch) bam_destroy1(srd);
	}));
	batch.clear();
      }
      // Wait for threads
      pool.waitAll();
      for(auto& fut : futures) fut.get();
    }    

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }

    // Sort traces
    std::sort(tr.begin(), tr.end());
  }

  template<typename TConfig>
  inline void
  brInTraces(TConfig const& c, std::vector<TraceRecord>& tr) {
    // Breakpoints
    typedef std::vector<Junction> TJunctionVector;
    typedef std::map<std::size_t, TJunctionVector> TReadBp;
    TReadBp readBp;
    findJunctions(c, readBp);

    // Search breakpoint insertions
    findBrIn(c, readBp, tr);
  }

}

#endif
