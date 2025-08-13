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

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {

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
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, false);
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
		_insertJunction(readBp, seed, rec, rp, spOrig, true);
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      if (bam_cigar_oplen(cigar[i]) > c.minRefSep) _insertJunction(readBp, seed, rec, rp, sp, false);
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
		_insertJunction(readBp, seed, rec, rpOrig, sp, true);
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
	      if (bam_cigar_oplen(cigar[i]) > c.minClip) _insertJunction(readBp, seed, rec, rp, finalsp, scleft);
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }

    // Sort junctions
    for(typename TReadBp::iterator it = readBp.begin(); it != readBp.end(); ++it) std::sort(it->second.begin(), it->second.end());

    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
  }


  template<typename TConfig, typename TReadBp>
  inline bool
  process_read(TConfig const& c, std::string const& searchseq, TReadBp& readBp, bam1_t* rec, TraceRecord& singleTR) {
    int32_t maxFragSize = searchseq.size() + 0.15 * searchseq.size();
    std::size_t seed = hash_lr(rec);
    std::string sequence;
		
    // Get all sequences embedded in breakpoints
    uint32_t insAlignLength = 0;
    uint32_t kLow = 0;
    uint32_t kHigh = 0;
    double pid = 0;
    double pIdStart = 0;
    double pIdEnd = 0;
    for(uint32_t k = 1; k < readBp[seed].size(); ++k) {
      if (readBp[seed][k].seqpos < (c.minSeedAlign + c.cropSize) ) continue;
      if (readBp[seed][k].seqpos + (c.minSeedAlign + c.cropSize) > rec->core.l_qseq) continue;
      if (!insAlignLength) kLow = k - 1;  // No previous hit
      int32_t fragsize = readBp[seed][k].seqpos - readBp[seed][kLow].seqpos;
      if ((fragsize > (c.minSeedAlign + c.cropSize)) && (fragsize < maxFragSize)) {
	// Insertion fragment?
	if (sequence.empty()) {
	  sequence.resize(rec->core.l_qseq);
	  const uint8_t* seqptr = bam_get_seq(rec);
	  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  // Reverse complement if necessary
	  if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);
	}

	// Full-length approach (good for transductions)
	std::string fullseq = sequence.substr(readBp[seed][kLow].seqpos, fragsize);
	EdlibAlignResult cigarFull = edlibAlign(fullseq.c_str(), fullseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	double pIdFull = 1.0 - ( (double) (cigarFull.editDistance) / (double) (fragsize) );
	edlibFreeAlignResult(cigarFull);
	if ((pIdFull > c.pctThres)  &&  ( (double) (fragsize) > ((0.25 * (double) (insAlignLength)) + insAlignLength))) {
	  kHigh = k;
	  insAlignLength = fragsize;
	  pid = pIdFull;
	  continue;
	}
	// Partial approach (good for rearranged insertions)
	
	// Check infix start (if not checked yet)
	if (!insAlignLength) {
	  std::string subseq = sequence.substr(readBp[seed][kLow].seqpos + c.cropSize, c.minSeedAlign);
	  
	  // Align to Insertions
	  EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	  pIdStart = 1.0 - ( (double) (cigar.editDistance) / (double) (c.minSeedAlign) );
	  //printAlignment(subseq, searchseq, EDLIB_MODE_HW, cigar);
	  edlibFreeAlignResult(cigar);
	  if (pIdStart <= c.pctThres) continue; // No hit
	}
	
	// Always check infix end
	int32_t sCoord = readBp[seed][k].seqpos - (c.minSeedAlign + c.cropSize);
	std::string subseq = sequence.substr(sCoord, c.minSeedAlign);
	
	// Align to Insertions
	EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	pIdEnd = 1.0 - ( (double) (cigar.editDistance) / (double) (c.minSeedAlign) );
	//printAlignment(subseq, searchseq, EDLIB_MODE_HW, cigar);
	edlibFreeAlignResult(cigar);
	if (pIdEnd <= c.pctThres) continue; // No hit
	
	// At least 25% larger hit?
	if ( (double) (fragsize) > ((0.25 * (double) (insAlignLength)) + insAlignLength) ) {
	  
	  // Align remaining internal segments
	  double percentIdentity = pIdStart + pIdEnd;
	  int32_t nCount = 2;
	  for(int32_t sIter = readBp[seed][kLow].seqpos + c.minSeedAlign + c.cropSize; sIter + c.minSeedAlign < (readBp[seed][k].seqpos - (c.minSeedAlign + c.cropSize)); sIter += c.minSeedAlign) {
	    std::string subseq = sequence.substr(sIter, c.minSeedAlign);
	    
	    // Align to Insertions
	    EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	    double pIdInfix = 1.0 - ( (double) (cigar.editDistance) / (double) (c.minSeedAlign) );
	    //printAlignment(subseq, searchseq, EDLIB_MODE_HW, cigar);
	    edlibFreeAlignResult(cigar);
	    percentIdentity += pIdInfix;
	    ++nCount;
	    //std::cerr << bam_get_qname(rec) << '\t' << readBp[seed][kLow].seqpos << ',' << sIter << ',' << pIdInfix << ',' << (percentIdentity / (double) nCount) << ',' << readBp[seed][k].seqpos << std::endl;
	  }
	  percentIdentity /= (double) nCount;
	  
	  // Entire sequence highly-similar?
	  if (percentIdentity > c.pctThres) {
	    kHigh = k;
	    insAlignLength = fragsize;
	    pid = percentIdentity;
	  }
	}
      }
    }
    bool validRead = false;
    if (insAlignLength) {
      // Make sure the leading and trailing sequence is NOT insertion sequence
      validRead = true;
      for(int32_t bp = 0; ((bp < 2) && (validRead)); ++bp) {
	int32_t fragsize = readBp[seed][kLow].seqpos;
	int32_t sCoord = 0;
	if (bp) {
	  sCoord = readBp[seed][kHigh].seqpos;
	  fragsize = sequence.size() - sCoord;
	}
	if (fragsize < maxFragSize) {   // Otherwise clearly more than insertion content
	  std::string subseq = sequence.substr(sCoord, fragsize);
	  
	  // Align to insertion
	  EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	  double percentIdentity = 1.0 - ( (double) (cigar.editDistance) / (double) (fragsize) );
	  //printAlignment(subseq, searchseq, EDLIB_MODE_HW, cigar);
	  edlibFreeAlignResult(cigar);		      
	  if (percentIdentity > c.pctThres) validRead = false;
	}
      }
      
      // Store read
      if (validRead) {
	// Canonical ordering
	if ((readBp[seed][kLow].refidx < readBp[seed][kHigh].refidx) || ( ( readBp[seed][kLow].refidx == readBp[seed][kHigh].refidx ) && (readBp[seed][kLow].refpos <= readBp[seed][kHigh].refpos) ) ) singleTR = TraceRecord(readBp[seed][kLow].refidx, readBp[seed][kLow].refpos, readBp[seed][kLow].seqpos, readBp[seed][kHigh].refidx, readBp[seed][kHigh].refpos, readBp[seed][kHigh].seqpos, (int32_t) (pid * 100), insAlignLength, seed);
	else singleTR = TraceRecord(readBp[seed][kHigh].refidx, readBp[seed][kHigh].refpos, readBp[seed][kHigh].seqpos, readBp[seed][kLow].refidx, readBp[seed][kLow].refpos, readBp[seed][kLow].seqpos, (int32_t) (pid * 100), insAlignLength, seed);
      }
    }
    return validRead;
  }
  

  template<typename TConfig, typename TReadBp>
  inline void
  findBrIn(TConfig const& c, TReadBp& readBp, std::vector<TraceRecord>& tr) {
    std::string searchseq;
    if (c.insmode == 0) {
      std::string faname = "";
      if (!loadSingleFasta(c, faname, searchseq)) return;
    }
    else if (c.insmode == 1) searchseq = MEI::alu;
    else if (c.insmode ==3) searchseq = MEI::sva;
    else searchseq = MEI::line1;

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
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Insertion serach" << std::endl;

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
      std::vector<bam1_t*> batch;
      uint32_t batch_size = 100;
      
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
		if (batch.size() == batch_size) {
		  for(auto* srd : batch) {
		    TraceRecord singleTR;
		    if (process_read(c, searchseq, readBp, srd, singleTR)) tr.push_back(singleTR);
		  }
		  for(auto* srd : batch) bam_destroy1(srd);
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
	for(auto* srd : batch) {
	  TraceRecord singleTR;
	  if (process_read(c, searchseq, readBp, srd, singleTR)) tr.push_back(singleTR);
	}
	for(auto* srd : batch) bam_destroy1(srd);
      }
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



  template<typename TConfig, typename TReadBp>
  inline void
  findBrInAll(TConfig const& c, TReadBp& readBp, std::vector<TraceRecord>& tr) {
    int32_t maxFragSize = 10000;

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
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Insertion serach" << std::endl;

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
      for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
	// Collect reads from all samples
	for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	  // Read alignments
	  hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	    // Keep only primary alignments
	    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	    std::size_t seed = hash_lr(rec);

	    if (clusteredReads.find(seed) != clusteredReads.end()) {
	      if (readBp[seed].size()>1) {
		// Get all sequences embedded in breakpoints
		uint32_t insAlignLength = 0;
		uint32_t kLow = 0;
		uint32_t kHigh = 0;
		double pid = 0;
		for(uint32_t k = 1; k < readBp[seed].size(); ++k) {
		  if (readBp[seed][k].seqpos < (c.minSeedAlign + c.cropSize) ) continue;
		  if (readBp[seed][k].seqpos + (c.minSeedAlign + c.cropSize) > rec->core.l_qseq) continue;
		  if (!insAlignLength) kLow = k - 1;  // No previous hit
		  int32_t fragsize = readBp[seed][k].seqpos - readBp[seed][kLow].seqpos;
		  if ((fragsize > (c.minSeedAlign + c.cropSize)) && (fragsize < maxFragSize)) {
		    // At least 25% larger hit?
		    if ( (double) (fragsize) > ((0.25 * (double) (insAlignLength)) + insAlignLength) ) {
		      if ( (readBp[seed][kLow].refidx != readBp[seed][k].refidx) || ( (readBp[seed][kLow].refidx == readBp[seed][k].refidx) && (std::abs(readBp[seed][kLow].refpos - readBp[seed][k].refpos) > 100000) ) ) {
			kHigh = k;
			insAlignLength = fragsize;
			pid = 1.0;
		      }
		    }
		  }
		}
		if (insAlignLength) {
		  // Canonical ordering
		  if ((readBp[seed][kLow].refidx < readBp[seed][kHigh].refidx) || ( ( readBp[seed][kLow].refidx == readBp[seed][kHigh].refidx ) && (readBp[seed][kLow].refpos <= readBp[seed][kHigh].refpos) ) ) tr.push_back(TraceRecord(readBp[seed][kLow].refidx, readBp[seed][kLow].refpos, readBp[seed][kLow].seqpos, readBp[seed][kHigh].refidx, readBp[seed][kHigh].refpos, readBp[seed][kHigh].seqpos, (int32_t) (pid * 100), insAlignLength, seed));
		  else tr.push_back(TraceRecord(readBp[seed][kHigh].refidx, readBp[seed][kHigh].refpos, readBp[seed][kHigh].seqpos, readBp[seed][kLow].refidx, readBp[seed][kLow].refpos, readBp[seed][kLow].seqpos, (int32_t) (pid * 100), insAlignLength, seed));
		}
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }
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
