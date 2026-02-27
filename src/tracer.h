#ifndef TRACER_H
#define TRACER_H

#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h>
#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include "edlib.h"
#include "util.h"
#include "junction.h"
#include "cluster.h"
#include "consensus.h"

namespace breaktracer {


  struct TracerConfig {
    bool hasGenomeMask;
    uint16_t insmode;
    uint16_t minMapQual;
    uint32_t minRefSep;
    uint32_t minClip;
    uint32_t graphPruning;
    uint32_t batchSize;
    uint32_t minCliqueSize;
    uint32_t maxReadPerSV;
    uint32_t maxThreads;
    int32_t nchr;
    int32_t minSeedAlign;
    int32_t cropSize;
    float pctThres;
    float indelExtension;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    boost::filesystem::path mask;
    boost::filesystem::path insseq;
    std::vector<std::string> sampleName;
  };

  struct InsAnno {
    bool isRC;
    bool td3Inv;
    bool td5Inv;
    int32_t posLeft;  // Breakpoint homology
    int32_t posRight;
    int32_t pos2Left;
    int32_t pos2Right;
    int32_t polyT;  // poly-T
    int32_t polyA;  // poly-A
    int32_t insStart;  // coordinates with respect to search seq
    int32_t insEnd;
    int32_t td5Len; // transduction
    int32_t td3Len;
    float consId;
    std::string cigar;
    std::string td5Seq;
    std::string td3Seq;
    
    InsAnno() : isRC(false), td3Inv(false), td5Inv(false), posLeft(0), posRight(0), pos2Left(0), pos2Right(0), polyT(0), polyA(0), insStart(0), insEnd(0), td5Len(0), td3Len(0), consId(0), cigar("*"), td5Seq("*"), td3Seq("*") {}
  };

  inline int32_t
  estimatePolyTail(std::string const& seq, bool const checkPolyA, float const minIdentity = 0.9) {
    int32_t len = 0;
    int32_t matches = 0;
    int32_t total = 0;
    if (checkPolyA) {
      for (size_t i = 0; i < seq.size(); ++i) {
        if (seq[i] == 'A') ++matches;
        ++total;
	if ((float) matches / (float) total >= minIdentity) len = total;
	else if ((total > 10) && ((float) matches / (float) total < 0.8)) break;
      }
    } else { 
      for (int i = seq.size() - 1; i >= 0; --i) {
        if (seq[i] == 'T') ++matches;
        ++total;
        if ((float) matches / (float) total >= minIdentity) len = total;
	else if ((total > 10) && ((float) matches / (float) total < 0.8)) break;
      }
    }
    return len;
  }

  template<typename TConfig>
  inline bool
  checkForInversion(TConfig const& c, std::string const& seq5, std::string const& meiseq) {
    EdlibAlignResult invRes = edlibAlign(seq5.c_str(), seq5.size(), meiseq.c_str(), meiseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
    double invId = 1.0 - ((double)invRes.editDistance / (double)seq5.size());
    edlibFreeAlignResult(invRes);
    if (invId > c.pctThres) return true;
    else return false;
  }
  
  template<typename TConfig>
  inline void
  annotateHomology(TConfig const& c, bam_hdr_t* hdr, faidx_t* fai, char const* seq, std::string const& searchseq, BrInTrace const& sv, InsAnno& ianno) {
    // Left homology length (pos)
    int minHomLen = 50;
    while (minHomLen) {
      if ((sv.pos >= minHomLen) && (sv.pos + minHomLen <= (int) hdr->target_len[sv.chr])) {
	std::string refseq = boost::to_upper_copy(std::string(seq + sv.pos - minHomLen, seq + sv.pos));
	EdlibAlignResult cigar = edlibAlign(refseq.c_str(), refseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	double pId = 1.0 - ( (double) (cigar.editDistance) / (double) (refseq.size()) );
	edlibFreeAlignResult(cigar);
	if (pId > c.pctThres) {
	  ianno.posLeft = minHomLen;
	  minHomLen += 10;
	} else break;
      } else break;
    }

    // Right homology length (pos)
    minHomLen = 50;
    while (minHomLen) {
      if ((sv.pos >= minHomLen) && (sv.pos + minHomLen <= (int) hdr->target_len[sv.chr])) {
	std::string refseq = boost::to_upper_copy(std::string(seq + sv.pos, seq + sv.pos + minHomLen));
	EdlibAlignResult cigar = edlibAlign(refseq.c_str(), refseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	double pId = 1.0 - ( (double) (cigar.editDistance) / (double) (refseq.size()) );
	edlibFreeAlignResult(cigar);
	if (pId > c.pctThres) {
	  ianno.posRight = minHomLen;
	  minHomLen += 10;
	} else break;
      } else break;
    }

    // Left homology length (pos2)
    minHomLen = 50;
    while (minHomLen) {
      std::string refseq;
      if ((sv.pos2 >= minHomLen) && (sv.pos2 + minHomLen <= (int) hdr->target_len[sv.chr2])) {
	if (sv.chr == sv.chr2) refseq = boost::to_upper_copy(std::string(seq + sv.pos2 - minHomLen, seq + sv.pos2));
	else {
	  // Different chromosome
	  int32_t seqlen = -1;
	  std::string chrom(hdr->target_name[sv.chr2]);
	  char* locseq = NULL;
	  locseq = faidx_fetch_seq(fai, chrom.c_str(), sv.pos2 - minHomLen, sv.pos2, &seqlen);
	  refseq = std::string(locseq);
	  if (locseq != NULL) free(locseq);
	}
      }
      if (!refseq.empty()) {
	EdlibAlignResult cigar = edlibAlign(refseq.c_str(), refseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	double pId = 1.0 - ( (double) (cigar.editDistance) / (double) (refseq.size()) );
	edlibFreeAlignResult(cigar);
	if (pId > c.pctThres) {
	  ianno.pos2Left = minHomLen;
	  minHomLen += 10;
	} else break;
      } else break;
    }

    // Right homology length (pos2)
    minHomLen = 50;
    while (minHomLen) {
      std::string refseq;
      if ((sv.pos2 >= minHomLen) && (sv.pos2 + minHomLen <= (int) hdr->target_len[sv.chr2])) {
	if (sv.chr == sv.chr2) refseq = boost::to_upper_copy(std::string(seq + sv.pos2, seq + sv.pos2 + minHomLen));
	else {
	  // Different chromosome
	  int32_t seqlen = -1;
	  std::string chrom(hdr->target_name[sv.chr2]);
	  char* locseq = NULL;
	  locseq = faidx_fetch_seq(fai, chrom.c_str(), sv.pos2, sv.pos2 + minHomLen, &seqlen);
	  refseq = std::string(locseq);
	  if (locseq != NULL) free(locseq);
	}
      }
      if (!refseq.empty()) {
	EdlibAlignResult cigar = edlibAlign(refseq.c_str(), refseq.size(), searchseq.c_str(), searchseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	double pId = 1.0 - ( (double) (cigar.editDistance) / (double) (refseq.size()) );
	edlibFreeAlignResult(cigar);
	if (pId > c.pctThres) {
	  ianno.pos2Right = minHomLen;
	  minHomLen += 10;
	} else break;
      } else break;
    }

    // Characterize insertion sequence
    if (!sv.consensus.empty()) {
      // Forward or reverse?
      int32_t halfLen = searchseq.size() / 2;
      std::string fwdseq = searchseq.substr(0, halfLen);
      double fwdId = 0.0;
      std::string revseq = searchseq.substr(halfLen, halfLen);
      double revId = 0.0;
      EdlibAlignResult resFwd = edlibAlign(sv.consensus.c_str(), sv.consensus.size(), fwdseq.c_str(), fwdseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
      if (resFwd.status == EDLIB_STATUS_OK) fwdId = 1.0 - ( (double) (resFwd.editDistance) / (double) (sv.consensus.size()) );
      EdlibAlignResult resRev = edlibAlign(sv.consensus.c_str(), sv.consensus.size(), revseq.c_str(), revseq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
      if (resRev.status == EDLIB_STATUS_OK) revId = 1.0 - ( (double) (resRev.editDistance) / (double) (sv.consensus.size()) );
      if (fwdId > revId) ianno.isRC = false;
      else ianno.isRC = true;
      edlibFreeAlignResult(resFwd);
      edlibFreeAlignResult(resRev);

      // Any transductions?
      std::string meiSeq;
      if (ianno.isRC) {
	meiSeq = revseq;
	if ((c.insmode > 0) && (c.insmode < 4)) meiSeq = meiSeq.substr(MEI::polyA.size()); // Trim poly-T
      } else {
	meiSeq = fwdseq;
	if ((c.insmode > 0) && (c.insmode < 4)) meiSeq = meiSeq.substr(0, meiSeq.size() - MEI::polyA.size()); // Trim poly-A
      }
      int32_t matchStart = 0;
      int32_t matchEnd = sv.consensus.size();
      EdlibAlignResult trRes = edlibAlign(meiSeq.c_str(), meiSeq.size(), sv.consensus.c_str(), sv.consensus.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
      //printAlignment(meiSeq, sv.consensus, EDLIB_MODE_HW, trRes);
      if ((trRes.status == EDLIB_STATUS_OK) && (trRes.numLocations > 0)) {
	matchStart = trRes.startLocations[0];
	matchEnd = trRes.endLocations[0];
	if (matchStart > 10) {
	  std::string seq5 = sv.consensus.substr(0, matchStart);
	  if (ianno.isRC) {
	    ianno.polyT = estimatePolyTail(seq5, false); // poly-T
	    seq5 = seq5.substr(0, seq5.size() - ianno.polyT);
	    if (seq5.size() > 10) {
	      ianno.td3Inv = checkForInversion(c, seq5, fwdseq);
	      ianno.td3Len = seq5.size();
	      ianno.td3Seq = seq5;
	    }
	  } else {
	    ianno.polyA = estimatePolyTail(seq5, true); // poly-A
	    seq5 = seq5.substr(ianno.polyA);
	    if (seq5.size() > 10) {
	      ianno.td5Inv = checkForInversion(c, seq5, revseq);
	      ianno.td5Len = seq5.size();
	      ianno.td5Seq = seq5;
	    }
	  }
	}
	if (matchEnd < (int) sv.consensus.size() - 11) {
	  std::string seq5 = sv.consensus.substr(matchEnd + 1);
	  if (ianno.isRC) {
	    ianno.polyT = estimatePolyTail(seq5, false); // poly-T
	    ianno.td5Inv = checkForInversion(c, seq5, fwdseq);
	    ianno.td5Len = (int) sv.consensus.size() - 1 - matchEnd;
	    ianno.td5Seq = seq5;
	  } else {
	    ianno.polyA = estimatePolyTail(seq5, true); // poly-A
	    ianno.td3Inv = checkForInversion(c, seq5, revseq);
	    ianno.td3Len = (int) sv.consensus.size() - 1 - matchEnd;
	    ianno.td3Seq = seq5;
	  }
	}
      }
      edlibFreeAlignResult(trRes);

      // Do proper alignment calculation without transductions
      std::string consBody = sv.consensus.substr(matchStart, (matchEnd - matchStart));
      if (!consBody.empty()) {
	EdlibAlignResult res = edlibAlign(consBody.c_str(), consBody.size(), meiSeq.c_str(), meiSeq.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	if (res.status == EDLIB_STATUS_OK) {
	  ianno.consId = 1.0 - ( (double) (res.editDistance) / (double) (consBody.size()) );
	  if (res.numLocations > 0) {
	    if (ianno.isRC) {
	      ianno.insStart = meiSeq.size() - 1 - res.endLocations[0];
	      ianno.insEnd = meiSeq.size() - 1 - res.startLocations[0];
	    } else {
	      ianno.insStart = res.startLocations[0];
	      ianno.insEnd = res.endLocations[0];
	    }
	  }
	}
	char* cig = edlibAlignmentToCigar(res.alignment, res.alignmentLength, EDLIB_CIGAR_STANDARD);
	ianno.cigar = std::string(cig);
	free(cig);
	edlibFreeAlignResult(res);
      }
    }
  }
  
  template<typename TConfig>
  inline void
  output(TConfig const& c, std::vector<BrInTrace>& sv) {
    // Sort by chromosome
    std::sort(sv.begin(), sv.end());

    // Load insertion reference sequence
    std::string searchseq;
    if (c.insmode == 0) {
      std::string faname = "";
      if (!loadSingleFasta(c, faname, searchseq)) return;
    } else {
      if (c.insmode == 1) searchseq = MEI::alu;
      else if (c.insmode == 3) searchseq = MEI::sva;
      else if (c.insmode == 4) searchseq = MEI::numt;
      else searchseq = MEI::line1;
      if (c.insmode != 4) searchseq += MEI::polyA;
    }
    int32_t minComplexOffset = searchseq.size() * 2;
    std::string revseq = searchseq;
    reverseComplement(revseq);
    searchseq += revseq;

    // ALT allele string based on insertion mode
    std::string altStr;
    if (c.insmode == 1) altStr = "<INS:ME:ALU>";
    else if (c.insmode == 2) altStr = "<INS:ME:L1>";
    else if (c.insmode == 3) altStr = "<INS:ME:SVA>";
    else if (c.insmode == 4) altStr = "<INS:MT>";
    else altStr = "<INS>";

    // Open BAM file for header and annotation
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.files[0].string().c_str());
    bam_hdr_t* bamhd = sam_hdr_read(samfile);
    faidx_t* fai = fai_load(c.genome.string().c_str());

    // Open output VCF/BCF (binary BCF to file, text VCF to stdout)
    std::string fmtout = "wb";
    if (c.outfile.string() == "-") fmtout = "w";
    htsFile* fp = hts_open(c.outfile.string().c_str(), fmtout.c_str());
    bcf_hdr_t* vcfhdr = bcf_hdr_init("w");

    // VCF header
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::gregorian::date today = now.date();
    std::string datestr("##fileDate=");
    datestr += boost::gregorian::to_iso_string(today);
    bcf_hdr_append(vcfhdr, datestr.c_str());
    bcf_hdr_append(vcfhdr, "##ALT=<ID=INS,Description=\"Insertion\">");
    bcf_hdr_append(vcfhdr, "##ALT=<ID=INS:ME:ALU,Description=\"ALU mobile element insertion\">");
    bcf_hdr_append(vcfhdr, "##ALT=<ID=INS:ME:L1,Description=\"LINE1 mobile element insertion\">");
    bcf_hdr_append(vcfhdr, "##ALT=<ID=INS:ME:SVA,Description=\"SVA mobile element insertion\">");
    bcf_hdr_append(vcfhdr, "##ALT=<ID=INS:MT,Description=\"Nuclear mitochondrial insertion\">");
    bcf_hdr_append(vcfhdr, "##FILTER=<ID=LowQual,Description=\"Low quality insertion.\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for POS2 in case of an inter-chromosomal event\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Genomic position for CHR2 in case of an inter-chromosomal event\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Estimated insertion fragment size\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=BRINTYPE,Number=1,Type=String,Description=\"Breakpoint insertion type\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=INSSEQID,Number=1,Type=Float,Description=\"Percent identity of consensus to insertion reference sequence (0-1)\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=INSSEQSTART,Number=1,Type=Integer,Description=\"Start coordinate in insertion reference sequence\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=INSSEQEND,Number=1,Type=Integer,Description=\"End coordinate in insertion reference sequence\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=INSSEQSTRAND,Number=1,Type=String,Description=\"Strand of insertion sequence\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=INSSEQCIGAR,Number=1,Type=String,Description=\"CIGAR of consensus alignment to insertion reference sequence\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=POLYT,Number=1,Type=Integer,Description=\"poly-T length\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=POLYA,Number=1,Type=Integer,Description=\"poly-A length\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=TD5LEN,Number=1,Type=Integer,Description=\"5-prime transduction length\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=TD5SEQ,Number=1,Type=String,Description=\"5-prime transduction sequence\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=TD5INV,Number=0,Type=Flag,Description=\"5-prime transduction inverted\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=TD3LEN,Number=1,Type=Integer,Description=\"3-prime transduction length\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=TD3SEQ,Number=1,Type=String,Description=\"3-prime transduction sequence\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=TD3INV,Number=0,Type=Flag,Description=\"3-prime transduction inverted\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=HOMLEN1L,Number=1,Type=Integer,Description=\"Bp1 left homology length\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=HOMLEN1R,Number=1,Type=Integer,Description=\"Bp1 right homology length\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=HOMLEN2L,Number=1,Type=Integer,Description=\"Bp2 left homology length\">");
    bcf_hdr_append(vcfhdr, "##INFO=<ID=HOMLEN2R,Number=1,Type=Integer,Description=\"Bp2 right homology length\">");
    std::string refloc("##reference=");
    refloc += c.genome.string();
    bcf_hdr_append(vcfhdr, refloc.c_str());
    for (int i = 0; i < bamhd->n_targets; ++i) {
      std::string refname("##contig=<ID=");
      refname += std::string(bamhd->target_name[i]) + ",length=" + std::to_string(bamhd->target_len[i]) + ">";
      bcf_hdr_append(vcfhdr, refname.c_str());
    }
    bcf_hdr_add_sample(vcfhdr, NULL);
    if (bcf_hdr_write(fp, vcfhdr) != 0) std::cerr << "Error: Failed to write VCF header!" << std::endl;

    // Iterate SVs and write VCF records
    char* seq = NULL;
    int lastIdx = -1;
    bcf1_t* rec = bcf_init();
    for (uint32_t i = 0; i < sv.size(); ++i) {
      // Lazy loading of reference genome
      if (sv[i].chr != lastIdx) {
	if (seq != NULL) { free(seq); seq = NULL; }
	int32_t seqlen = -1;
	std::string chrom(bamhd->target_name[sv[i].chr]);
	seq = faidx_fetch_seq(fai, chrom.c_str(), 0, bamhd->target_len[sv[i].chr], &seqlen);
	lastIdx = sv[i].chr;
      }

      // Insertion annotation
      InsAnno ianno;
      annotateHomology(c, bamhd, fai, seq, searchseq, sv[i], ianno);
      char strand = ianno.isRC ? '-' : '+';

      // Breakpoint insertion type
      int32_t offset = std::abs(sv[i].pos - sv[i].pos2);
      std::string brintype;
      if (sv[i].chr != sv[i].chr2) brintype = "InterChromosomalSVwithInsertion";
      else if (offset > minComplexOffset) brintype = "IntraChromosomalSVwithInsertion";
      else brintype = "PlainInsertion";

      // ID
      std::string padNumber = std::to_string(i);
      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
      std::string id = "BRINS" + padNumber;

      // Build record
      bcf_clear(rec);
      rec->rid = bcf_hdr_name2id(vcfhdr, bamhd->target_name[sv[i].chr]);
      rec->pos = sv[i].pos;  // 0-based (htslib internal)
      bcf_update_id(vcfhdr, rec, id.c_str());
      std::string alleles(1, (char) std::toupper((unsigned char) seq[sv[i].pos]));
      alleles += sv[i].consensus.empty() ? "," + altStr : "," + sv[i].consensus;
      bcf_update_alleles_str(vcfhdr, rec, alleles.c_str());
      int32_t qvalout = sv[i].mapq * sv[i].seeds.size();
      if (qvalout < 0) qvalout = 0;
      if (qvalout > 10000) qvalout = 10000;
      rec->qual = qvalout;
      int32_t fltPass = bcf_hdr_id2int(vcfhdr, BCF_DT_ID, "PASS");
      if (qvalout < 60) { fltPass = bcf_hdr_id2int(vcfhdr, BCF_DT_ID, "LowQual"); }
      bcf_update_filter(vcfhdr, rec, &fltPass, 1);

      // INFO fields
      bool precise = !sv[i].consensus.empty();
      if (precise) bcf_update_info_flag(vcfhdr, rec, "PRECISE", NULL, 1);
      else bcf_update_info_flag(vcfhdr, rec, "IMPRECISE", NULL, 1);
      bcf_update_info_string(vcfhdr, rec, "SVTYPE", "INS");
      std::string svmethod("BreakTracer-v");
      svmethod += breaktracerVersionNumber;
      bcf_update_info_string(vcfhdr, rec, "SVMETHOD", svmethod.c_str());

      int32_t tmpi;
      if (sv[i].chr == sv[i].chr2) {
	tmpi = sv[i].pos2 + 1;  // 1-based VCF INFO/END
	bcf_update_info_int32(vcfhdr, rec, "END", &tmpi, 1);
      } else {
	tmpi = sv[i].pos + 1;   // 1-based
	bcf_update_info_int32(vcfhdr, rec, "END", &tmpi, 1);
	bcf_update_info_string(vcfhdr, rec, "CHR2", bamhd->target_name[sv[i].chr2]);
	tmpi = sv[i].pos2 + 1;  // 1-based VCF INFO/POS2
	bcf_update_info_int32(vcfhdr, rec, "POS2", &tmpi, 1);
      }
      tmpi = sv[i].inslen;
      bcf_update_info_int32(vcfhdr, rec, "SVLEN", &tmpi, 1);
      bcf_update_info_string(vcfhdr, rec, "BRINTYPE", brintype.c_str());
      tmpi = (int32_t) sv[i].seeds.size();
      bcf_update_info_int32(vcfhdr, rec, "SR", &tmpi, 1);
      tmpi = sv[i].mapq;
      bcf_update_info_int32(vcfhdr, rec, "MAPQ", &tmpi, 1);
      float tmpf = (float) sv[i].qual / (float) 100.0;
      bcf_update_info_float(vcfhdr, rec, "INSSEQID", &tmpf, 1);

      if (precise) {
	tmpi = ianno.insStart;
	bcf_update_info_int32(vcfhdr, rec, "INSSEQSTART", &tmpi, 1);
	tmpi = ianno.insEnd;
	bcf_update_info_int32(vcfhdr, rec, "INSSEQEND", &tmpi, 1);
	std::string strandStr(1, strand);
	bcf_update_info_string(vcfhdr, rec, "INSSEQSTRAND", strandStr.c_str());
	bcf_update_info_string(vcfhdr, rec, "INSSEQCIGAR", ianno.cigar.c_str());
	tmpi = ianno.polyT;
	bcf_update_info_int32(vcfhdr, rec, "POLYT", &tmpi, 1);
	tmpi = ianno.polyA;
	bcf_update_info_int32(vcfhdr, rec, "POLYA", &tmpi, 1);
	if (ianno.td5Len > 0) {
	  tmpi = ianno.td5Len;
	  bcf_update_info_int32(vcfhdr, rec, "TD5LEN", &tmpi, 1);
	  bcf_update_info_string(vcfhdr, rec, "TD5SEQ", ianno.td5Seq.c_str());
	  if (ianno.td5Inv) bcf_update_info_flag(vcfhdr, rec, "TD5INV", NULL, 1);
	}
	if (ianno.td3Len > 0) {
	  tmpi = ianno.td3Len;
	  bcf_update_info_int32(vcfhdr, rec, "TD3LEN", &tmpi, 1);
	  bcf_update_info_string(vcfhdr, rec, "TD3SEQ", ianno.td3Seq.c_str());
	  if (ianno.td3Inv) bcf_update_info_flag(vcfhdr, rec, "TD3INV", NULL, 1);
	}
	tmpi = ianno.posLeft;
	bcf_update_info_int32(vcfhdr, rec, "HOMLEN1L", &tmpi, 1);
	tmpi = ianno.posRight;
	bcf_update_info_int32(vcfhdr, rec, "HOMLEN1R", &tmpi, 1);
	tmpi = ianno.pos2Left;
	bcf_update_info_int32(vcfhdr, rec, "HOMLEN2L", &tmpi, 1);
	tmpi = ianno.pos2Right;
	bcf_update_info_int32(vcfhdr, rec, "HOMLEN2R", &tmpi, 1);
      }
      bcf_write1(fp, vcfhdr, rec);
    }
    bcf_destroy(rec);
    if (seq != NULL) free(seq);
    fai_destroy(fai);

    // Clean-up
    bam_hdr_destroy(bamhd);
    hts_idx_destroy(idx);
    sam_close(samfile);
    bcf_hdr_destroy(vcfhdr);
    hts_close(fp);

    // Build index for binary BCF output
    if (c.outfile.string() != "-") bcf_index_build(c.outfile.string().c_str(), 14);
  }

  
  template<typename TConfig>
  inline int32_t
  runTracer(TConfig& c) {

#ifdef PROFILE
    ProfilerStart("tracer.prof");
#endif

    // Search breakpoint insertion traces
    typedef std::vector<TraceRecord> TTraceVector;
    TTraceVector tr;
    brInTraces(c, tr);
    //for(uint32_t i = 0; i < tr.size(); ++i) std::cerr << tr[i].chr << ':' << tr[i].pos << '\t' << tr[i].chr2 << ':' << tr[i].pos2 << std::endl;
    
    // Cluster reads
    typedef std::vector<BrInTrace> TBrInVector;
    TBrInVector sv;
    cluster(c, tr, sv);
    
    // Assemble
    assemble(c, tr, sv);
    
    // Output
    output(c, sv);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    
    return 0;
  }

  int tracer(int argc, char **argv) {
    TracerConfig c;
   
    // Parameter
    std::string mode;
    std::string instag;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("technology,y", boost::program_options::value<std::string>(&mode)->default_value("ont"), "seq. technology [pb, ont]")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "BCF output file")
      ("threads,t", boost::program_options::value<uint32_t>(&c.maxThreads)->default_value(8), "number of threads")
      ;
    
    boost::program_options::options_description disc("Split-read options");
    disc.add_options()
      ("minrefsep,m", boost::program_options::value<uint32_t>(&c.minRefSep)->default_value(30), "min. reference separation")
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("min-clique-size,z", boost::program_options::value<uint32_t>(&c.minCliqueSize)->default_value(3), "min. clique size")
      ("max-reads,p", boost::program_options::value<uint32_t>(&c.maxReadPerSV)->default_value(15), "max. reads for local assembly")
      ;
    
    boost::program_options::options_description brin("Breakpoint insertion options");
    brin.add_options()
      ("mask,k", boost::program_options::value<boost::filesystem::path>(&c.mask), "genome mask")
      ("cropsize,r", boost::program_options::value<int32_t>(&c.cropSize)->default_value(20), "leading/trailing crop size")
      ("seedlen,s", boost::program_options::value<int32_t>(&c.minSeedAlign)->default_value(130), "min. seed length")
      ("pctid,i", boost::program_options::value<float>(&c.pctThres)->default_value(0.9), "min. percent identity")
      ("instag,n", boost::program_options::value<std::string>(&instag)->default_value("L1"), "Type of insertion [ALU|L1|SVA|NUMT]")
      ("insseq,e", boost::program_options::value<boost::filesystem::path>(&c.insseq), "FASTA with insertion sequence [overrides -n]")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
      ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "graph pruning cutoff")
      ("batchsize,b", boost::program_options::value<uint32_t>(&c.batchSize)->default_value(100), "multi-threading batch size")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(disc).add(brin).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(disc).add(brin);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: breaktracer " << argv[0] << " [OPTIONS] -g <ref.fa> <sample1.sort.bam> <sample2.sort.bam> ..." << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }
    
    // Clique size
    if (c.minCliqueSize < 2) c.minCliqueSize = 2;
    
    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.genome.string().c_str()) == -1) {
	  std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.genome.string().c_str());
      }
      fai_destroy(fai);
    }
   
    // Check input files
    c.sampleName.resize(c.files.size());
    c.nchr = 0;
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
	std::cerr << "Alignment file is missing: " << c.files[file_c].string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.files[file_c].string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.files[file_c].string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      if (!c.nchr) c.nchr = hdr->n_targets;
      else {
	if (c.nchr != hdr->n_targets) {
	  std::cerr << "BAM files have different number of chromosomes!" << std::endl;
	  return 1;
	}
      }
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (!faidx_has_seq(fai, tname.c_str())) {
	  std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	  return 1;
	}
      }
      fai_destroy(fai);
      std::string sampleName = "unknown";
      getSMTag(std::string(hdr->text), c.files[file_c].stem().string(), sampleName);
      c.sampleName[file_c] = sampleName;
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
    checkSampleNames(c);
    
    // Genome mask
    if (vm.count("mask")) c.hasGenomeMask = true;
    else c.hasGenomeMask = false;
   
    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }
    
    // Long-read options
    if (mode == "pb") c.indelExtension = 0.7;
    else if (mode == "ont") c.indelExtension = 0.5;
    
   // Insertion mode, 0=FASTA, 1=ALU, 2=L1, 3=SVA
    if (vm.count("insseq")) {
      c.insmode = 0;
      if (!(boost::filesystem::exists(c.insseq) && boost::filesystem::is_regular_file(c.insseq) && boost::filesystem::file_size(c.insseq))) {
	std::cerr << "Insertion sequence file is missing: " << c.insseq.string() << std::endl;
	return 1;
      }
    }
    else if (instag == "ALU") c.insmode = 1;
    else if (instag == "SVA") c.insmode = 3;
    else if (instag == "NUMT") c.insmode = 4;
    else c.insmode = 2;
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "breaktracer ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    return runTracer(c);
  }



}

#endif
