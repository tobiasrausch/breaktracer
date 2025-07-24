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

  // Junction record
  struct Junction {
    bool forward;
    bool scleft;
    int32_t refidx;
    int32_t rstart;
    int32_t refpos;
    int32_t seqpos;
    uint16_t qual;

    Junction(bool const fw, bool const cl, int32_t const idx, int32_t const rst, int32_t const r, int32_t const s, uint16_t const qval) : forward(fw), scleft(cl), refidx(idx), rstart(rst), refpos(r), seqpos(s), qual(qval) {}

    bool operator<(const Junction& j2) const {
      return ((seqpos<j2.seqpos) || ((seqpos==j2.seqpos) && (refidx<j2.refidx)) || ((seqpos==j2.seqpos) && (refidx==j2.refidx) && (refpos<j2.refpos)) || ((seqpos==j2.seqpos) && (refidx==j2.refidx) && (refpos==j2.refpos) && (scleft < j2.scleft)));
    }
  };

  
  template<typename TReadBp>
  inline void
    _insertJunction(TReadBp& readBp, std::size_t const seed, bam1_t* rec, int32_t const rp, int32_t const sp, bool const scleft) {
    bool fw = true;
    if (rec->core.flag & BAM_FREVERSE) fw = false;
    int32_t readStart = rec->core.pos;
    if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) readStart = -1;
    typedef typename TReadBp::mapped_type TJunctionVector;
    typename TReadBp::iterator it = readBp.find(seed);
    int32_t seqlen = readLength(rec);
    if (sp <= seqlen) {
      if (rec->core.flag & BAM_FREVERSE) {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, readStart, rp, seqlen - sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, readStart, rp, seqlen - sp, rec->core.qual))));
      } else {
	if (it != readBp.end()) it->second.push_back(Junction(fw, scleft, rec->core.tid, readStart, rp, sp, rec->core.qual));
	else readBp.insert(std::make_pair(seed, TJunctionVector(1, Junction(fw, scleft, rec->core.tid, readStart, rp, sp, rec->core.qual))));
      }
    }
  }

  inline int32_t
  _selectReadStart(std::vector<Junction> const& jcvec) {
    for(uint32_t i = 0; i < jcvec.size(); ++i) {
      if (jcvec[i].rstart != -1) return jcvec[i].rstart;
    }
    return -1;
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
	for(int32_t refIndex = 0; refIndex < c.nchr; ++refIndex) {
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
  inline void
  findL1(TConfig const& c, TReadBp& readBp) {
    std::string line1("GGGGGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCGACGCAGAAGACGGTGATTTCTGCATTTCCATCTGAGGTACCGGGTTCATCTCACTAGGGAGTGCCAGACAGTGGGCGCAGGCCAGTGTGTGTGCGCACCGTGCGCGAGCCGAAGCAGGGCGAGGCATTGCCTCACCTGGGAAGCGCAAGGGGTCAGGGAGTTCCCTTTCTGAGTCAAAGAAAGGGGTGACGGTCGCACCTGGAAAATCGGGTCACTCCCACCCGAATATTGCGCTTTTCAGACCGGCTTAAGAAACGGCGCACCACGAGACTATATCCCACACCTGGCTCGGAGGGTCCTACGCCCACGGAATCTCGCTGATTGCTAGCACAGCAGTCTGAGATCAAACTGCAAGGCGGCAACGAGGCTGGGGGAGGGGCGCCCGCCATTGCCCAGGCTTGCTTAGGTAAACAAAGCAGCCGGGAAGCTCGAACTGGGTGGAGCCCACCACAGCTCAAGGAGGCCTGCCTGCCTCTGTAGGCTCCACCTCTGGGGGCAGGGCACAGACAAACAAAAAGACAGCAGTAACCTCTGCAGACTTAAGTGTCCCTGTCTGACAGCTTTGAAGAGAGCAGTGGTTCTCCCAGCACGCAGCTGGAGATCTGAGAACGGGCAGACAGACTGCCTCCTCAAGTGGGTCCCTGACTCCTGACCCCCGAGCAGCCTAACTGGGAGGCACCCCCCAGCAGGGGCACACTGACACCTCACACGGCAGGGTATTCCAACAGACCTGCAGCTGAGGGTCCTGTCTGTTAGAAGGAAAACTAACAACCAGAAAGGACATCTACACCGAAAACCCATCTGTACATCACCATCATCAAAGACCAAAAGTAGATAAAACCACAAAGATGGGGAAAAAACAGAACAGAAAAACTGGAAACTCTAAAACGCAGAGCGCCTCTCCTCCTCCAAAGGAACGCAGTTCCTCACCAGCAACGGAACAAAGCTGGATGGAGAATGATTTTGACGAGCTGAGAGAAGAAGGCTTCAGACGATCAAATTACTCTGAGCTACGGGAGGACATTCAAACCAAAGGCAAAGAAGTTGAAAACTTTGAAAAAAATTTAGAAGAATGTATAACTAGAATAACCAATACAGAGAAGTGCTTAAAGGAGCTGATGGAGCTGAAAACCAAGGCTCGAGAACTACGTGAAGAATGCAGAAGCCTCAGGAGCCGATGCGATCAACTGGAAGAAAGGGTATCAGCAATGGAAGATGAAATGAATGAAATGAAGCGAGAAGGGAAGTTTAGAGAAAAAAGAATAAAAAGAAATGAGCAAAGCCTCCAAGAAATATGGGACTATGTGAAAAGACCAAATCTACGTCTGATTGGTGTACCTGAAAGTGATGTGGAGAATGGAACCAAGTTGGAAAACACTCTGCAGGATATTATCCAGGAGAACTTCCCCAATCTAGCAAGGCAGGCCAACGTTCAGATTCAGGAAATACAGAGAACGCCACAAAGATACTCCTCGAGAAGAGCAACTCCAAGACACATAATTGTCAGATTCACCAAAGTTGAAATGAAGGAAAAAATGTTAAGGGCAGCCAGAGAGAAAGGTCGGGTTACCCTCAAAGGAAAGCCCATCAGACTAACAGTGGATCTCTCGGCAGAAACCCTACAAGCCAGAAGAGAGTGGGGGCCAATATTCAACATTCTTAAAGAAAAGAATTTTCAACCCAGAATTTCATATCCAGCCAAACTAAGCTTCATAAGTGAAGGAGAAATAAAATACTTTATAGACAAGCAAATGTTGAGAGATTTTGTCACCACCAGGCCTGCCCTAAAAGAGCTCCTGAAGGAAGCGCTAAACATGGAAAGGAACAACCGGTACCAGCCGCTGCAAAATCATGCCAAAATGTAAAGACCATCGAGACTAGGAAGAAACTGCATCAACTAATGAGCAAAATCACCAGCTAACATCATAATGACAGGATCAAATTCACACATAACAATATTAACTTTAAATATAAATGGACTAAATTCTGCAATTAAAAGACACAGACTGGCAAGTTGGATAAAGAGTCAAGACCCATCAGTGTGCTGTATTCAGGAAACCCATCTCACGTGCAGAGACACACATAGGCTCAAAATAAAAGGATGGAGGAAGATCTACCAAGCCAATGGAAAACAAAAAAAGGCAGGGGTTGCAATCCTAGTCTCTGATAAAACAGACTTTAAACCAACAAAGATCAAAAGAGACAAAGAAGGCCATTACATAATGGTAAAGGGATCAATTCAACAAGAGGAGCTAACTATCCTAAATATTTATGCACCCAATACAGGAGCACCCAGATTCATAAAGCAAGTCCTCAGTGACCTACAAAGAGACTTAGACTCCCACACATTAATAATGGGAGACTTTAACACCCCACTGTCAACATTAGACAGATCAACGAGACAGAAAGTCAACAAGGATACCCAGGAATTGAACTCAGCTCTGCACCAAGCAGACCTAATAGACATCTACAGAACTCTCCACCCCAAATCAACAGAATATACCTTTTTTTCAGCACCACACCACACCTATTCCAAAATTGACCACATAGTTGGAAGTAAAGCTCTCCTCAGCAAATGTAAAAGAACAGAAATTATAACAAACTATCTCTCAGACCACAGTGCAATCAAACTAGAACTCAGGATTAAGAATCTCACTCAAAGCCGCTCAACTACATGGAAACTGAACAACCTGCTCCTGAATGACTACTGGGTACATAACGAAATGAAGGCAGAAATAAAGATGTTCTTTGAAACCAACGAGAACAAAGACACCACATACCAGAATCTCTGGGACGCATTCAAAGCAGTGTGTAGAGGGAAATTTATAGCACTAAATGCCTACAAGAGAAAGCAGGAAAGATCCAAAATTGACACCCTAACATCACAATTAAAAGAACTAGAAAAGCAAGAGCAAACACATTCAAAAGCTAGCAGAAGGCAAGAAATAACTAAAATCAGAGCAGAACTGAAGGAAATAGAGACACAAAAAACCCTTCAAAAAATCAATGAATCCAGGAGCTGGTTTTTTGAAAGGATCAACAAAATTGATAGACCGCTAGCAAGACTAATAAAGAAAAAAAGAGAGAAGAATCAAATAGACACAATAAAAAATGATAAAGGGGATATCACCACCGATCCCACAGAAATACAAACTACCATCAGAGAATACTACAAACACCTCTACGCAAATAAACTAGAAAATCTAGAAGAAATGGATACATTCCTCGACACATACACTCTCCCAAGACTAAACCAGGAAGAAGTTGAATCTCTGAATAGACCAATAACAGGCTCTGAAATTGTGGCAATAATCAATAGTTTACCAACCAAAAAGAGTCCAGGACCAGATGGATTCACAGCCGAATTCTACCAGAGGTACATGGAGGAACTGGTACCATTCCTTCTGAAACTATTCCAATCAATAGAAAAAGAGGGAATCCTCCCTAACTCATTTTATGAGGCCAGCATCATTCTGATACCAAAGCCGGGCAGAGACACAACCAAAAAAGAGAATTTTAGACCAATATCCTTGATGAACATTGATGCAAAAATCCTCAATAAAATACTGGCAAACCGAATCCAGCAGCACATCAAAAAGCTTATCCACCATGATCAAGTGGGCTTCATCCCTGGGATGCAAGGCTGGTTCAATATACGCAAATCAATAAATGTAATCCAGCATATAAACAGAGCCAAAGACAAAAACCACATGATTATCTCAATAGATGCAGAAAAAGCCTTTGACAAAATTCAACAACCCTTCATGCTAAAAACTCTCAATAAATTAGGTATTGATGGGACGTATTTCAAAATAATAAGAGCTATCTATGACAAACCCACAGCCAATATCATACTGAATGGGCAAAAACTGGAAGCATTCCCTTTGAAAACCGGCACAAGACAGGGATGCCCTCTCTCACCGCTCCTATTCAACATAGTGTTGGAAGTTCTGGCCAGGGCAATCAGGCAGGAGAAGGAAATAAAGGGTATTCAATTAGGAAAAGAGGAAGTCAAATTGTCCCTGTTTGCAGACGACATGATTGTATATCTAGAAAACCCCATCGTCTCAGCCCAAAATCTCCTTAAGCTGATAAGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTACAAAAATCACAAGCATTCTTATACACCAACAACAGACAAACAGAGAGCCAAATCATGGGTGAACTCCCATTCGTAATTGCTTCAAAGAGAATAAAATACCTAGGAATCCAACTTACAAGGGATGTGAAGGACCTCTTCAAGGAGAACTACAAACCACTGCTCAAGGAAATAAAAGAGGACACAAACAAATGGAAGAACATTCCATGCTCATGGGTAGGAAGAATCAATATCGTGAAAATGGCCATACTGCCCAAGGTAATTTACAGATTCAATGCCATCCCCATCAAGCTACCAATGACTTTCTTCACAGAATTGGAAAAAACTACTTTAAAGTTCATATGGAACCAAAAAAGAGCCCGCATTGCCAAGTCAATCCTAAGCCAAAAGAACAAAGCTGGAGGCATCACACTACCTGACTTCAAACTATACTACAAGGCTACAGTAACCAAAACAGCATGGTACTGGTACCAAAACAGAGATATAGATCAATGGAACAGAACAGAGCCCTCAGAAATAATGCCGCATATCTACAACTATCTGATCTTTGACAAACCTGAGAAAAACAAGCAATGGGGAAAGGATTCCCTATTTAATAAATGGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATCAATTCAAGATGGATTAAAGATTTAAACGTTAAACCTAAAACCATAAAAACCCTAGAAGAAAACCTAGGCATTACCATTCAGGACATAGGCGTGGGCAAGGACTTCATGTCCAAAACACCAAAAGCAATGGCAACAAAAGACAAAATTGACAAATGGGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTACCATCAGAGTGAACAGGCAACCTACAACATGGGAGAAAATTTTCGCAACCTACTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAACTTAAACAAATTTACAAGAAAAAAACAAACAACCCCATCAAAAAGTGGGCGAAGGACATGAACAGACACTTCTCAAAAGAAGACATTTATGCAGCCAAAAAACACATGAAGAAATGCTCATCATCACTGGCCATCAGAGAAATGCAAATCAAAACCACTATGAGATATCATCTCACACCAGTTAGAATGGCAATCATTAAAAAGTCAGGAAACAACAGGTGCTGGAGAGGATGCGGAGAAATAGGAACACTTTTACACTGTTGGTGGGACTGTAAACTAGTTCAACCATTGTGGAAGTCAGTGTGGCGATTCCTCAGGGATCTAGAACTAGAAATACCATTTGACCCAGCCATCCCATTACTGGGTATATACCCAAATGAGTATAAATCATGCTGCTATAAAGACACATGCACACGTATGTTTATTGCGGCACTATTCACAATAGCAAAGACTTGGAACCAACCCAAATGTCCAACAATGATAGACTGGATTAAGAAAATGTGGCACATATACACCATGGAATACTATGCAGCCATAAAAAATGATGAGTTCATATCCTTTGTAGGGACATGGATGAAATTGGAAACCATCATTCTCAGTAAACTATCGCAAGAACAAAAAACCAAACACCGCATATTCTCACTCATAGGTGGGAATTGAACAATGAGATCACATGGACACAGGAAGGGGAATATCACACTCTGGGGACTGTGGTGGGGTCGGGGGAGGGGGGAGGGATAGCATTGGGAGATATACCTAATGCTAGATGACACATTAGTGGGTGCAGCGCACCAGCATGGCACATGTATACATATGTAACTAACCTGCACAATGTGCACATGTACCCTAAAACTTAGAGTATATACTCTAAGTTTTAGGGTACATGTGCACATTGTGCAGGTTAGTTACATATGTATACATGTGCCATGCTGGTGCGCTGCACCCACTAATGTGTCATCTAGCATTAGGTATATCTCCCAATGCTATCCCTCCCCCCTCCCCCGACCCCACCACAGTCCCCAGAGTGTGATATTCCCCTTCCTGTGTCCATGTGATCTCATTGTTCAATTCCCACCTATGAGTGAGAATATGCGGTGTTTGGTTTTTTGTTCTTGCGATAGTTTACTGAGAATGATGGTTTCCAATTTCATCCATGTCCCTACAAAGGATATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGTATATGTGCCACATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAAACATACGTGTGCATGTGTCTTTATAGCAGCATGATTTATACTCATTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACAGTCCCACCAACAGTGTAAAAGTGTTCCTATTTCTCCGCATCCTCTCCAGCACCTGTTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGATATCTCATAGTGGTTTTGATTTGCATTTCTCTGATGGCCAGTGATGATGAGCATTTCTTCATGTGTTTTTTGGCTGCATAAATGTCTTCTTTTGAGAAGTGTCTGTTCATGTCCTTCGCCCACTTTTTGATGGGGTTGTTTGTTTTTTTCTTGTAAATTTGTTTAAGTTCATTGTAGATTCTGGATATTAGCCCTTTGTCAGATGAGTAGGTTGCGAAAATTTTCTCCCATGTTGTAGGTTGCCTGTTCACTCTGATGGTAGTTTCTTTTGCTGTGCAGAAGCTCTTTAGTTTAATTAGATCCCATTTGTCAATTTTGTCTTTTGTTGCCATTGCTTTTGGTGTTTTGGACATGAAGTCCTTGCCCACGCCTATGTCCTGAATGGTAATGCCTAGGTTTTCTTCTAGGGTTTTTATGGTTTTAGGTTTAACGTTTAAATCTTTAATCCATCTTGAATTGATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACCATTTATTAAATAGGGAATCCTTTCCCCATTGCTTGTTTTTCTCAGGTTTGTCAAAGATCAGATAGTTGTAGATATGCGGCATTATTTCTGAGGGCTCTGTTCTGTTCCATTGATCTATATCTCTGTTTTGGTACCAGTACCATGCTGTTTTGGTTACTGTAGCCTTGTAGTATAGTTTGAAGTCAGGTAGTGTGATGCCTCCAGCTTTGTTCTTTTGGCTTAGGATTGACTTGGCAATGCGGGCTCTTTTTTGGTTCCATATGAACTTTAAAGTAGTTTTTTCCAATTCTGTGAAGAAAGTCATTGGTAGCTTGATGGGGATGGCATTGAATCTGTAAATTACCTTGGGCAGTATGGCCATTTTCACGATATTGATTCTTCCTACCCATGAGCATGGAATGTTCTTCCATTTGTTTGTGTCCTCTTTTATTTCCTTGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCACATCCCTTGTAAGTTGGATTCCTAGGTATTTTATTCTCTTTGAAGCAATTACGAATGGGAGTTCACCCATGATTTGGCTCTCTGTTTGTCTGTTGTTGGTGTATAAGAATGCTTGTGATTTTTGTACATTGATTTTGTATCCTGAGACTTTGCTGAAGTTGCTTATCAGCTTAAGGAGATTTTGGGCTGAGACGATGGGGTTTTCTAGATATACAATCATGTCGTCTGCAAACAGGGACAATTTGACTTCCTCTTTTCCTAATTGAATACCCTTTATTTCCTTCTCCTGCCTGATTGCCCTGGCCAGAACTTCCAACACTATGTTGAATAGGAGCGGTGAGAGAGGGCATCCCTGTCTTGTGCCGGTTTTCAAAGGGAATGCTTCCAGTTTTTGCCCATTCAGTATGATATTGGCTGTGGGTTTGTCATAGATAGCTCTTATTATTTTGAAATACGTCCCATCAATACCTAATTTATTGAGAGTTTTTAGCATGAAGGGTTGTTGAATTTTGTCAAAGGCTTTTTCTGCATCTATTGAGATAATCATGTGGTTTTTGTCTTTGGCTCTGTTTATATGCTGGATTACATTTATTGATTTGCGTATATTGAACCAGCCTTGCATCCCAGGGATGAAGCCCACTTGATCATGGTGGATAAGCTTTTTGATGTGCTGCTGGATTCGGTTTGCCAGTATTTTATTGAGGATTTTTGCATCAATGTTCATCAAGGATATTGGTCTAAAATTCTCTTTTTTGGTTGTGTCTCTGCCCGGCTTTGGTATCAGAATGATGCTGGCCTCATAAAATGAGTTAGGGAGGATTCCCTCTTTTTCTATTGATTGGAATAGTTTCAGAAGGAATGGTACCAGTTCCTCCATGTACCTCTGGTAGAATTCGGCTGTGAATCCATCTGGTCCTGGACTCTTTTTGGTTGGTAAACTATTGATTATTGCCACAATTTCAGAGCCTGTTATTGGTCTATTCAGAGATTCAACTTCTTCCTGGTTTAGTCTTGGGAGAGTGTATGTGTCGAGGAATGTATCCATTTCTTCTAGATTTTCTAGTTTATTTGCGTAGAGGTGTTTGTAGTATTCTCTGATGGTAGTTTGTATTTCTGTGGGATCGGTGGTGATATCCCCTTTATCATTTTTTATTGTGTCTATTTGATTCTTCTCTCTTTTTTTCTTTATTAGTCTTGCTAGCGGTCTATCAATTTTGTTGATCCTTTCAAAAAACCAGCTCCTGGATTCATTGATTTTTTGAAGGGTTTTTTGTGTCTCTATTTCCTTCAGTTCTGCTCTGATTTTAGTTATTTCTTGCCTTCTGCTAGCTTTTGAATGTGTTTGCTCTTGCTTTTCTAGTTCTTTTAATTGTGATGTTAGGGTGTCAATTTTGGATCTTTCCTGCTTTCTCTTGTAGGCATTTAGTGCTATAAATTTCCCTCTACACACTGCTTTGAATGCGTCCCAGAGATTCTGGTATGTGGTGTCTTTGTTCTCGTTGGTTTCAAAGAACATCTTTATTTCTGCCTTCATTTCGTTATGTACCCAGTAGTCATTCAGGAGCAGGTTGTTCAGTTTCCATGTAGTTGAGCGGCTTTGAGTGAGATTCTTAATCCTGAGTTCTAGTTTGATTGCACTGTGGTCTGAGAGATAGTTTGTTATAATTTCTGTTCTTTTACATTTGCTGAGGAGAGCTTTACTTCCAACTATGTGGTCAATTTTGGAATAGGTGTGGTGTGGTGCTGAAAAAAAGGTATATTCTGTTGATTTGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCTGCTTGGTGCAGAGCTGAGTTCAATTCCTGGGTATCCTTGTTGACTTTCTGTCTCGTTGATCTGTCTAATGTTGACAGTGGGGTGTTAAAGTCTCCCATTATTAATGTGTGGGAGTCTAAGTCTCTTTGTAGGTCACTGAGGACTTGCTTTATGAATCTGGGTGCTCCTGTATTGGGTGCATAAATATTTAGGATAGTTAGCTCCTCTTGTTGAATTGATCCCTTTACCATTATGTAATGGCCTTCTTTGTCTCTTTTGATCTTTGTTGGTTTAAAGTCTGTTTTATCAGAGACTAGGATTGCAACCCCTGCCTTTTTTTGTTTTCCATTGGCTTGGTAGATCTTCCTCCATCCTTTTATTTTGAGCCTATGTGTGTCTCTGCACGTGAGATGGGTTTCCTGAATACAGCACACTGATGGGTCTTGACTCTTTATCCAACTTGCCAGTCTGTGTCTTTTAATTGCAGAATTTAGTCCATTTATATTTAAAGTTAATATTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTGATTTTGCTCATTAGTTGATGCAGTTTCTTCCTAGTCTCGATGGTCTTTACATTTTGGCATGATTTTGCAGCGGCTGGTACCGGTTGTTCCTTTCCATGTTTAGCGCTTCCTTCAGGAGCTCTTTTAGGGCAGGCCTGGTGGTGACAAAATCTCTCAACATTTGCTTGTCTATAAAGTATTTTATTTCTCCTTCACTTATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTTTCTTTAAGAATGTTGAATATTGGCCCCCACTCTCTTCTGGCTTGTAGGGTTTCTGCCGAGAGATCCACTGTTAGTCTGATGGGCTTTCCTTTGAGGGTAACCCGACCTTTCTCTCTGGCTGCCCTTAACATTTTTTCCTTCATTTCAACTTTGGTGAATCTGACAATTATGTGTCTTGGAGTTGCTCTTCTCGAGGAGTATCTTTGTGGCGTTCTCTGTATTTCCTGAATCTGAACGTTGGCCTGCCTTGCTAGATTGGGGAAGTTCTCCTGGATAATATCCTGCAGAGTGTTTTCCAACTTGGTTCCATTCTCCACATCACTTTCAGGTACACCAATCAGACGTAGATTTGGTCTTTTCACATAGTCCCATATTTCTTGGAGGCTTTGCTCATTTCTTTTTATTCTTTTTTCTCTAAACTTCCCTTCTCGCTTCATTTCATTCATTTCATCTTCCATTGCTGATACCCTTTCTTCCAGTTGATCGCATCGGCTCCTGAGGCTTCTGCATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCTCCATCAGCTCCTTTAAGCACTTCTCTGTATTGGTTATTCTAGTTATACATTCTTCTAAATTTTTTTCAAAGTTTTCAACTTCTTTGCCTTTGGTTTGAATGTCCTCCCGTAGCTCAGAGTAATTTGATCGTCTGAAGCCTTCTTCTCTCAGCTCGTCAAAATCATTCTCCATCCAGCTTTGTTCCGTTGCTGGTGAGGAACTGCGTTCCTTTGGAGGAGGAGAGGCGCTCTGCGTTTTAGAGTTTCCAGTTTTTCTGTTCTGTTTTTTCCCCATCTTTGTGGTTTTATCTACTTTTGGTCTTTGATGATGGTGATGTACAGATGGGTTTTCGGTGTAGATGTCCTTTCTGGTTGTTAGTTTTCCTTCTAACAGACAGGACCCTCAGCTGCAGGTCTGTTGGAATACCCTGCCGTGTGAGGTGTCAGTGTGCCCCTGCTGGGGGGTGCCTCCCAGTTAGGCTGCTCGGGGGTCAGGAGTCAGGGACCCACTTGAGGAGGCAGTCTGTCTGCCCGTTCTCAGATCTCCAGCTGCGTGCTGGGAGAACCACTGCTCTCTTCAAAGCTGTCAGACAGGGACACTTAAGTCTGCAGAGGTTACTGCTGTCTTTTTGTTTGTCTGTGCCCTGCCCCCAGAGGTGGAGCCTACAGAGGCAGGCAGGCCTCCTTGAGCTGTGGTGGGCTCCACCCAGTTCGAGCTTCCCGGCTGCTTTGTTTACCTAAGCAAGCCTGGGCAATGGCGGGCGCCCCTCCCCCAGCCTCGTTGCCGCCTTGCAGTTTGATCTCAGACTGCTGTGCTAGCAATCAGCGAGATTCCGTGGGCGTAGGACCCTCCGAGCCAGGTGTGGGATATAGTCTCGTGGTGCGCCGTTTCTTAAGCCGGTCTGAAAAGCGCAATATTCGGGTGGGAGTGACCCGATTTTCCAGGTGCGACCGTCACCCCTTTCTTTGACTCAGAAAGGGAACTCCCTGACCCCTTGCGCTTCCCAGGTGAGGCAATGCCTCGCCCTGCTTCGGCTCGCGCACGGTGCGCACACACACTGGCCTGCGCCCACTGTCTGGCACTCCCTAGTGAGATGAACCCGGTACCTCAGATGGAAATGCAGAAATCACCGTCTTCTGCGTCGCTCACGCTGGGAGCTGTAGACCGGAGCTGTTCCTATTCGGCCATCTTGGCTCCTCCCCC");
    int32_t minSeedAlign = 130;
    int32_t cropSize = 20;  // To cover poly-A tail or micro-insertions at the breakpoint
    int32_t maxFragSize = 7000;
    float pctThres = 0.9;

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
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] " << "L1 serach" << std::endl;

    // L1 candidate read search
    std::set<std::size_t> clusteredReads;
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet bp(hdr->target_len[refIndex], false);
      TBitSet firstHit(hdr->target_len[refIndex], false);
      TBitSet secondHit(hdr->target_len[refIndex], false);

      // Cluster breakpoints +/-5bp
      for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	if (it->second.size() > 1) {
	  for(uint32_t i = 0; i < it->second.size(); ++i) {
	    if (it->second[i].refidx == refIndex) {
	      int32_t refStart = 0;
	      if (it->second[i].refpos > 5) refStart = it->second[i].refpos - 5;
	      int32_t refEnd = hdr->target_len[refIndex];
	      if (it->second[i].refpos + 5 < refEnd) refEnd = it->second[i].refpos + 5;
	      for(int32_t k = refStart; k < refEnd; ++k) {
		if (secondHit[k]) bp[k] = true;
		else {
		  if (firstHit[k]) secondHit[k] = true;
		  else firstHit[k] = true;
		}
	      }
	    }
	  }
	}
      }
      
      // 2nd pass to fetch reads
      for(typename TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	if (it->second.size() > 1) {
	  for(uint32_t i = 0; i < it->second.size(); ++i) {
	    if (it->second[i].refidx == refIndex) {
	      int32_t refStart = 0;
	      if (it->second[i].refpos > 5) refStart = it->second[i].refpos - 5;
	      int32_t refEnd = hdr->target_len[refIndex];
	      if (it->second[i].refpos + 5 < refEnd) refEnd = it->second[i].refpos + 5;
	      for(int32_t k = refStart; k < refEnd; ++k) {
		if (bp[k]) clusteredReads.insert(it->first);
	      }
	    }
	  }
	}
      }
      //std::cerr << refIndex << ',' << clusteredReads.size() << std::endl;
    }

    // Find L1
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
		std::string sequence;
		
		// Get all sequences embedded in breakpoints
		uint32_t l1AlignLength = 0;
		uint32_t kLow = 0;
		uint32_t kHigh = 0;
		double pid = 0;
		double pIdStart = 0;
		double pIdEnd = 0;
		for(uint32_t k = 1; k < readBp[seed].size(); ++k) {
		  if (readBp[seed][k].seqpos < (minSeedAlign + cropSize) ) continue;
		  if (readBp[seed][k].seqpos + (minSeedAlign + cropSize) > rec->core.l_qseq) continue;
		  if (!l1AlignLength) kLow = k - 1;  // No previous hit
		  int32_t fragsize = readBp[seed][k].seqpos - readBp[seed][kLow].seqpos;
		  if ((fragsize > (minSeedAlign + cropSize)) && (fragsize < maxFragSize)) {
		    // L1 fragment?
		    if (sequence.empty()) {
		      sequence.resize(rec->core.l_qseq);
		      const uint8_t* seqptr = bam_get_seq(rec);
		      for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		      // Reverse complement if necessary
		      if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);
		    }

		    // Full-length approach (good for transductions)
		    std::string fullseq = sequence.substr(readBp[seed][kLow].seqpos, fragsize);
                    EdlibAlignResult cigarFull = edlibAlign(fullseq.c_str(), fullseq.size(), line1.c_str(), line1.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		    double pIdFull = 1.0 - ( (double) (cigarFull.editDistance) / (double) (fragsize) );
		    edlibFreeAlignResult(cigarFull);
		    if ((pIdFull > pctThres)  &&  ( (double) (fragsize) > ((0.25 * (double) (l1AlignLength)) + l1AlignLength))) {
		      kHigh = k;
		      l1AlignLength = fragsize;
		      pid = pIdFull;
		      continue;
		    }
		    // Partial approach (good for rearranged L1s)

		    // Check infix start (if not checked yet)
		    if (!l1AlignLength) {
		      std::string subseq = sequence.substr(readBp[seed][kLow].seqpos + cropSize, minSeedAlign);
		      
		      // Align to L1
		      EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), line1.c_str(), line1.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		      pIdStart = 1.0 - ( (double) (cigar.editDistance) / (double) (minSeedAlign) );
		      //printAlignment(subseq, line1, EDLIB_MODE_HW, cigar);
		      edlibFreeAlignResult(cigar);
		      if (pIdStart <= pctThres) continue; // No hit
		    }
		    
		    // Always check infix end
		    int32_t sCoord = readBp[seed][k].seqpos - (minSeedAlign + cropSize);
		    std::string subseq = sequence.substr(sCoord, minSeedAlign);
		      
		      // Align to L1
		    EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), line1.c_str(), line1.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		    pIdEnd = 1.0 - ( (double) (cigar.editDistance) / (double) (minSeedAlign) );
		    //printAlignment(subseq, line1, EDLIB_MODE_HW, cigar);
		    edlibFreeAlignResult(cigar);
		    if (pIdEnd <= pctThres) continue; // No hit

		    // At least 25% larger hit?
		    if ( (double) (fragsize) > ((0.25 * (double) (l1AlignLength)) + l1AlignLength) ) {

		      // Align remaining internal segments
		      double percentIdentity = pIdStart + pIdEnd;
		      int32_t nCount = 2;
		      for(int32_t sIter = readBp[seed][kLow].seqpos + minSeedAlign + cropSize; sIter + minSeedAlign < (readBp[seed][k].seqpos - (minSeedAlign + cropSize)); sIter += minSeedAlign) {
			std::string subseq = sequence.substr(sIter, minSeedAlign);
		    
			// Align to L1
			EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), line1.c_str(), line1.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
			double pIdInfix = 1.0 - ( (double) (cigar.editDistance) / (double) (minSeedAlign) );
			//printAlignment(subseq, line1, EDLIB_MODE_HW, cigar);
			edlibFreeAlignResult(cigar);
			percentIdentity += pIdInfix;
			++nCount;
			//std::cerr << bam_get_qname(rec) << '\t' << readBp[seed][kLow].seqpos << ',' << sIter << ',' << pIdInfix << ',' << (percentIdentity / (double) nCount) << ',' << readBp[seed][k].seqpos << std::endl;
		      }
		      percentIdentity /= (double) nCount;

		      // Entire sequence highly-similar?
		      if (percentIdentity > pctThres) {
			kHigh = k;
			l1AlignLength = fragsize;
			pid = percentIdentity;
		      }
		    }
		  }
		}
		if (l1AlignLength) {
		  // Make sure the leading and trailing sequence is NOT L1 sequence
		  bool validRead = true;
		  for(int32_t bp = 0; ((bp < 2) && (validRead)); ++bp) {
		    int32_t fragsize = readBp[seed][kLow].seqpos;
		    int32_t sCoord = 0;
		    if (bp) {
		      sCoord = readBp[seed][kHigh].seqpos;
		      fragsize = sequence.size() - sCoord;
		    }
		    if (fragsize < maxFragSize) {   // Otherwise clearly more than L1 content
		      std::string subseq = sequence.substr(sCoord, fragsize);
		      
		      // Align to L1
		      EdlibAlignResult cigar = edlibAlign(subseq.c_str(), subseq.size(), line1.c_str(), line1.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		      double percentIdentity = 1.0 - ( (double) (cigar.editDistance) / (double) (fragsize) );
		      //printAlignment(subseq, line1, EDLIB_MODE_HW, cigar);
		      edlibFreeAlignResult(cigar);		      
		      if (percentIdentity > pctThres) validRead = false;
		    }
		  }

		  // Output read
		  if (validRead) {
		    int32_t offset = std::abs(readBp[seed][kLow].refpos - readBp[seed][kHigh].refpos);
		    if (readBp[seed][kLow].refidx != readBp[seed][kHigh].refidx) std::cout << "Inter-chromosomal SV with L1 fragment" << '\t';
		    else if (offset > 1000) std::cout << "Intra-chromosomal SV with L1 fragment" << '\t';
		    else std::cout << "L1 insertion" << '\t';
		    std::cout << bam_get_qname(rec) << '\t';
		    std::cout << hdr->target_name[readBp[seed][kLow].refidx] << ':' << readBp[seed][kLow].refpos << " (ReadPos: " << readBp[seed][kLow].seqpos << ')' << '\t';
		    std::cout << hdr->target_name[readBp[seed][kHigh].refidx] << ':' << readBp[seed][kHigh].refpos << " (ReadPos: " << readBp[seed][kHigh].seqpos << ')' << '\t';
		    std::cout << l1AlignLength << '\t';
		    std::cout << pid << std::endl;
		  }
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
  }   
  
  template<typename TConfig>
  inline void
  brInTraces(TConfig const& c) {
    // Breakpoints
    typedef std::vector<Junction> TJunctionVector;
    typedef std::map<std::size_t, TJunctionVector> TReadBp;
    TReadBp readBp;
    findJunctions(c, readBp);

    // Line1
    std::cout << "L1InsType\tReadName\tRefCoordBefore\tRefCoordAfter\tL1FragmentSize\tL1Similarity" << std::endl;
    findL1(c, readBp);    
  }


}

#endif
