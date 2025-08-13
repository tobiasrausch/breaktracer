#ifndef CONSENSUS_H
#define CONSENSUS_H

#include <iostream>

#include "edlib.h"

namespace breaktracer
{

  template<typename TConfig, typename TAlign>
  inline void
  consensus(TConfig const& c, TAlign const& align, std::string& gapped, std::string& cs) {
    typedef typename TAlign::index TAIndex;

    // Calculate coverage
    typedef boost::multi_array<bool, 2> TFlag;
    TFlag fl(boost::extents[align.shape()[0]][align.shape()[1]]);
    typedef std::vector<int> TCoverage;
    TCoverage cov;
    cov.resize(align.shape()[1], 0);
    for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
      int start = 0;
      int end = -1;
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	fl[i][j] = false;
	if (align[i][j] != '-') end = j;
	else if (end == -1) start = j + 1;
      }
      for(TAIndex j = start; j<=end; ++j) {
	++cov[j];
	fl[i][j] = true;
      }
    }
    
    int covThreshold = c.minCliqueSize;
    TAIndex j = 0;
    std::vector<char> cons(align.shape()[1], '-');
    for(typename TCoverage::const_iterator itCov = cov.begin(); itCov != cov.end(); ++itCov, ++j) {
      int32_t maxIdx = 4;  // Leading/trailing gaps until min. coverage is reached
      if (*itCov >= covThreshold) {
	// Get consensus letter
	std::vector<int32_t> count(5, 0); // ACGT-
	for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	  if (fl[i][j]) {
	    if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	    else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	    else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	    else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	    else ++count[4];
	  }
	}
	maxIdx = 0;
	int32_t maxCount = count[0];
	for(uint32_t i = 1; i<5; ++i) {
	  if (count[i] > maxCount) {
	    maxCount = count[i];
	    maxIdx = i;
	  }
	}
      }
      switch (maxIdx) {
      case 0: cons[j] = 'A'; break;
      case 1: cons[j] = 'C'; break;
      case 2: cons[j] = 'G'; break;
      case 3: cons[j] = 'T'; break;
      default: break;
      }
    }
    gapped = std::string(cons.begin(), cons.end());
    for(uint32_t i = 0; i<cons.size(); ++i) {
      if (cons[i] != '-') cs.push_back(cons[i]);
    }
  }

  
  template<typename TAlign>
  inline void
  convertAlignment(std::string const& query, TAlign& align, EdlibAlignMode const modeCode, EdlibAlignResult const& cigar) {
    // Input alignment
    TAlign alignIn;
    alignIn.resize(boost::extents[align.shape()[0]][align.shape()[1]]);
    for(uint32_t i = 0; i < align.shape()[0]; ++i) {
      for(uint32_t j = 0; j < align.shape()[1]; ++j) {
	alignIn[i][j] = align[i][j];
      }
    }
	
    // Create new alignment
    uint32_t seqPos = alignIn.shape()[0];
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if (modeCode == EDLIB_MODE_HW) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) alignIn.shape()[1]) missingEnd = alignIn.shape()[1] - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
    }
    align.resize(boost::extents[alignIn.shape()[0]+1][missingStart + cigar.alignmentLength + missingEnd]);

    // infix alignment, fix start
    if (modeCode == EDLIB_MODE_HW) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) {
	  for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][j];
	  align[seqPos][j] = '-';
	}
      }
    }
    
    // target
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j + missingStart] = '-';
      } else {
	++tIdx;
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j + missingStart] = alignIn[seqIdx][tIdx];
      }
    }

    // query
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) align[seqPos][j + missingStart] = '-';
      else align[seqPos][j + missingStart] = query[++qIdx];
    }

    // infix alignment, fix end
    if (modeCode == EDLIB_MODE_HW) {
      if (missingEnd) {
	for (uint32_t j = cigar.alignmentLength + missingStart; j < cigar.alignmentLength + missingStart + missingEnd; ++j) {
	  ++tIdx;
	  for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][tIdx];
	  align[seqPos][j] = '-';
	}
      }
    }
  }

  template<typename TAlign>
  inline void
  consensusEdlib(TAlign const& align, std::string& cons) {
    typedef typename TAlign::index TAIndex;

    cons.resize(align.shape()[1]);
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      std::vector<int32_t> count(5, 0); // ACGT-
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	else ++count[4];
      }
      uint32_t maxIdx = 0;
      uint32_t sndIdx = 1;
      if (count[maxIdx] < count[sndIdx]) {
	maxIdx = 1;
	sndIdx = 0;
      }
      for(uint32_t i = 2; i<5; ++i) {
	if (count[i] > count[maxIdx]) {
	  sndIdx = maxIdx;
	  maxIdx = i;
	}
	else if (count[i] > count[sndIdx]) {
	  sndIdx = i;
	}
      }
      if (2 * count[sndIdx] < count[maxIdx]) {
	switch (maxIdx) {
	case 0: cons[j] = 'A'; break;
	case 1: cons[j] = 'C'; break;
	case 2: cons[j] = 'G'; break;
	case 3: cons[j] = 'T'; break;
	default: cons[j] = '-'; break;
	}
      } else {
	uint32_t k1 = maxIdx;
	uint32_t k2 = sndIdx;
	if (k1 > k2) {
	  k1 = sndIdx;
	  k2 = maxIdx;
	}
	// ACGT-
	if ((k1 == 0) && (k2 == 1)) cons[j] = 'M';
	else if ((k1 == 0) && (k2 == 2)) cons[j] = 'R';
	else if ((k1 == 0) && (k2 == 3)) cons[j] = 'W';
	else if ((k1 == 0) && (k2 == 4)) cons[j] = 'B';
	else if ((k1 == 1) && (k2 == 2)) cons[j] = 'S';
	else if ((k1 == 1) && (k2 == 3)) cons[j] = 'Y';
	else if ((k1 == 1) && (k2 == 4)) cons[j] = 'D';
	else if ((k1 == 2) && (k2 == 3)) cons[j] = 'K';
	else if ((k1 == 2) && (k2 == 4)) cons[j] = 'E';
	else if ((k1 == 3) && (k2 == 4)) cons[j] = 'F';
	else cons[j] = '-';
      }
    }
  }


  template<typename TConfig, typename TSplitReadSet>
  inline int
  msaEdlib(TConfig const& c, TSplitReadSet& sps, std::string& cs) {
    // Pairwise scores
    std::vector<int32_t> edit(sps.size() * sps.size(), 0);
    for(uint32_t i = 0; i < sps.size(); ++i) {
      for(uint32_t j = i + 1; j < sps.size(); ++j) {
	EdlibAlignResult align = edlibAlign(sps[i].c_str(), sps[i].size(), sps[j].c_str(), sps[j].size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	edit[i * sps.size() + j] = align.editDistance;
	edit[j * sps.size() + i] = align.editDistance;
	edlibFreeAlignResult(align);
      }
    }

    // Find best sequence to start alignment
    uint32_t bestIdx = 0;
    int32_t bestVal = sps[0].size();
    for(uint32_t i = 0; i < sps.size(); ++i) {
      std::vector<int32_t> dist(sps.size());
      for(uint32_t j = 0; j < sps.size(); ++j) dist[j] = edit[i * sps.size() + j];
      std::sort(dist.begin(), dist.end());
      if (dist[sps.size()/2] < bestVal) {
	bestVal = dist[sps.size()/2];
	bestIdx = i;
      }
    }
    
    // Align to best sequence
    std::vector<std::pair<int32_t, int32_t> > qscores;
    qscores.push_back(std::make_pair(0, bestIdx));
    std::string revc = sps[bestIdx];
    reverseComplement(revc);
    for(uint32_t j = 0; j < sps.size(); ++j) {
      if (j != bestIdx) {
	EdlibAlignResult align = edlibAlign(revc.c_str(), revc.size(), sps[j].c_str(), sps[j].size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	if (align.editDistance < edit[bestIdx * sps.size() + j]) {
	  reverseComplement(sps[j]);
	  qscores.push_back(std::make_pair(align.editDistance, j));
	} else qscores.push_back(std::make_pair(edit[bestIdx * sps.size() + j], j));
	edlibFreeAlignResult(align);
      }
    }
    std::sort(qscores.begin(), qscores.end());
    
    // Drop poorest 20% and order by centroid
    std::vector<uint32_t> selectedIdx;
    uint32_t lastIdx = (uint32_t) (0.8 * qscores.size());
    if (lastIdx < 3) lastIdx = 3;
    for(uint32_t i = 0; ((i < qscores.size()) && (i < lastIdx)); ++i) selectedIdx.push_back(qscores[i].second);
    
    // Extended IUPAC code
    EdlibEqualityPair additionalEqualities[20] = {{'M', 'A'}, {'M', 'C'}, {'R', 'A'}, {'R', 'G'}, {'W', 'A'}, {'W', 'T'}, {'B', 'A'}, {'B', '-'}, {'S', 'C'}, {'S', 'G'}, {'Y', 'C'}, {'Y', 'T'}, {'D', 'C'}, {'D', '-'}, {'K', 'G'}, {'K', 'T'}, {'E', 'G'}, {'E', '-'}, {'F', 'T'}, {'F', '-'}};

    // Incrementally align sequences    
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    align.resize(boost::extents[1][sps[selectedIdx[0]].size()]);
    uint32_t ind = 0;
    for(typename std::string::const_iterator str = sps[selectedIdx[0]].begin(); str != sps[selectedIdx[0]].end(); ++str) align[0][ind++] = *str;
    for(uint32_t i = 1; i < selectedIdx.size(); ++i) {
      // Convert to consensus
      std::string alignStr;
      consensusEdlib(align, alignStr);
      // Debug MSA
      //std::cerr << "Progressive MSA: " << i << '(' << align.shape()[0] << ':' << align.shape()[1] << ')' << std::endl;
      //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
      //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
      //std::cerr << std::endl;
      //}
      //std::cerr << "Consensus: " << std::endl;
      //std::cerr << alignStr << std::endl;
      //std::cerr << "ToBeAligned: " << sps[selectedIdx[i]] << std::endl;
      // Compute alignment
      EdlibAlignResult cigar = edlibAlign(sps[selectedIdx[i]].c_str(), sps[selectedIdx[i]].size(), alignStr.c_str(), alignStr.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, additionalEqualities, 20));
      convertAlignment(sps[selectedIdx[i]], align, EDLIB_MODE_NW, cigar);
      edlibFreeAlignResult(cigar);
    }
    
    // Debug MSA
    //std::cerr << "Output MSA " << '(' << align.shape()[0] << ':' << align.shape()[1] << ')' << std::endl;
    //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
    //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
    //std::cerr << std::endl;
    //}

    // Consensus
    std::string gapped;
    consensus(c, align, gapped, cs);
    //std::cerr << "Consensus:" << std::endl;
    //std::cerr << gapped << std::endl;

    // Return split-read support
    return align.shape()[0];
  }

  template<typename TConfig>
  inline void
  assemble(TConfig const& c, std::vector<TraceRecord>& tr, std::vector<BrInTrace>& sv) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "Local assembly" << std::endl;

    // Threads
    ThreadPool pool(std::max<std::size_t>(1, c.maxThreads));
    std::vector<std::future<void>> futures;
    
    // Link clustered reads back to SV 
    typedef std::set<std::size_t> TReadSet;
    typedef std::vector<uint32_t> TSVList; 
    typedef std::map<std::size_t, TSVList> TReadToSv;
    TReadToSv clusteredReads;
    for(uint32_t i = 0; i<sv.size(); ++i) {
      for(typename TReadSet::iterator itRead = sv[i].seeds.begin(); itRead != sv[i].seeds.end(); ++itRead) {
	typename TReadToSv::iterator itMap = clusteredReads.find((*itRead));
	if (itMap == clusteredReads.end()) clusteredReads.insert(std::make_pair((*itRead), TSVList()));
	clusteredReads[(*itRead)].push_back(i);
      }
    }
    
    // Link clustered reads back to trace record (simple 1:1 mapping, if a read has multiple insertions --> change!)
    typedef std::map<std::size_t, uint32_t> TReadToTrace;
    TReadToTrace trMap;
    for(uint32_t i = 0; i<tr.size(); ++i) trMap.insert(std::make_pair(tr[i].id, i));
    
    // Sequence store
    typedef std::vector<std::string> TSequences;
    typedef std::vector<TSequences> TSVSequences;
    TSVSequences seqStore(sv.size(), TSequences());

    // SV consensus done
    std::vector<bool> svcons(sv.size(), false);

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

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      // Collect reads from all samples
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  // Only primary alignments with the full sequence information
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	  
	  std::size_t seed = hash_lr(rec);
	  if (clusteredReads.find(seed) != clusteredReads.end()) {
	    // Get sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    const uint8_t* seqptr = bam_get_seq(rec);
	    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	    if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);

	    // Iterate all spanned SVs
	    for(uint32_t ri = 0; ri < clusteredReads[seed].size(); ++ri) {
	      int32_t svid = clusteredReads[seed][ri];
	      if ((!svcons[svid]) && (seqStore[svid].size() < c.maxReadPerSV)) {
		int32_t fragsize = std::abs(tr[trMap[seed]].seqpos - tr[trMap[seed]].seqpos2);
		int32_t sCoord = std::min(tr[trMap[seed]].seqpos, tr[trMap[seed]].seqpos2);
		// Debug
		//std::cerr << svid << ',' << sCoord << ',' << fragsize << std::endl;
		seqStore[svid].push_back(sequence.substr(sCoord, fragsize));

		// Enough split-reads?
		if ((seqStore[svid].size() == c.maxReadPerSV) || (seqStore[svid].size() == sv[svid].seeds.size())) {
		  futures.push_back(pool.enqueue([&, svid] {
		    if (seqStore[svid].size() > 1) {
		      std::string consensusStr;
		      TSequences seqMSA = seqStore[svid];
		      seqStore[svid].clear();
		      msaEdlib(c, seqMSA, consensusStr);
		      sv[svid].consensus = std::move(consensusStr);
		    }
		  }));
		  svcons[svid] = true;
		}
	      }
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }
    
    // Handle left-overs
    for(uint32_t svid = 0; svid < svcons.size(); ++svid) {	
      if (!svcons[svid]) {
	futures.push_back(pool.enqueue([&, svid] {
	  if (seqStore[svid].size() > 1) {
	    std::string consensusStr;
	    TSequences seqMSA = seqStore[svid];
	    seqStore[svid].clear();
	    msaEdlib(c, seqMSA, consensusStr);
	    sv[svid].consensus = std::move(consensusStr);
	  }
	}));
      }
    }

    // Threads
    pool.waitAll();
    for(auto& fut : futures) fut.get();
      
    // Clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
    
    // Clean-up unfinished SVs and failed local assemblies
    for(uint32_t svid = 0; svid < svcons.size(); ++svid) {
      if (!svcons[svid]) sv[svid].consensus = "";
      else {
	int32_t l1 = sv[svid].inslen;
	int32_t l2 = sv[svid].consensus.size();
	if (l1 > l2) {
	  int32_t tmpL = l1;
	  l1 = l2;
	  l2 = tmpL;
	}
	if ( ( (double) l1 / (double) l2 ) < 0.9 ) sv[svid].consensus = "";
      }
    }
  }


}

#endif
