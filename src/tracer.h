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
    uint16_t insmode;
    uint16_t minMapQual;
    uint32_t minRefSep;
    uint32_t minClip;
    uint32_t graphPruning;
    uint32_t minCliqueSize;
    uint32_t maxReadPerSV;
    int32_t nchr;
    int32_t minSeedAlign;
    int32_t cropSize;
    float pctThres;
    float indelExtension;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path genome;
    boost::filesystem::path insseq;
    std::vector<std::string> sampleName;
  };
  
  template<typename TConfig>
  inline void
  output(TConfig const& c, std::vector<BrInTrace>& sv) {
    // Open output file
    std::streambuf * buf;
    std::ofstream of;
    if(c.outfile.string() != "-") {
      of.open(c.outfile.string().c_str());
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream out(buf);
    out << "BrInId\tRefCoord1\tRefCoord2\tBrInEstFragmentSize\tPercentIdentity\tReadSupport\tBrInType\tBrInSeqSize\tBrInSequence" << std::endl;

    // Open file handles
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.files[0].string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    for(uint32_t i = 0; i < sv.size(); ++i) {
      std::string padNumber = boost::lexical_cast<std::string>(i);
      padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
      out << "BRINS" << padNumber << '\t';
      out << hdr->target_name[sv[i].chr] << ':' << sv[i].pos << '\t';
      out << hdr->target_name[sv[i].chr2] << ':' << sv[i].pos2 << '\t';
      out << sv[i].inslen << '\t';
      out << sv[i].qual << '\t';
      out << sv[i].seeds.size() << '\t';
      int32_t offset = std::abs(sv[i].pos - sv[i].pos2);
      //out << c.sampleName[file_c] << '\t';
      if (sv[i].chr != sv[i].chr2) out << "InterChromosomalSVwithInsertion";
      else if (offset > 1000) out << "IntraChromosomalSVwithInsertion";
      else out << "PlainInsertion";
      out << '\t';
      out << sv[i].consensus.size() << '\t';
      out << sv[i].consensus;
      out << std::endl;
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    if(c.outfile.string() != "-") of.close();
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
     ("cropsize,r", boost::program_options::value<int32_t>(&c.cropSize)->default_value(20), "leading/trailing crop size")
     ("seedlen,s", boost::program_options::value<int32_t>(&c.minSeedAlign)->default_value(130), "min. seed length")
     ("pctid,i", boost::program_options::value<float>(&c.pctThres)->default_value(0.9), "min. percent identity")
     ("instag,t", boost::program_options::value<std::string>(&instag)->default_value("L1"), "Type of insertion [ALU|L1|SVA]")
     ("insseq,e", boost::program_options::value<boost::filesystem::path>(&c.insseq), "FASTA with insertion sequence [overrides -t]")
     ;

   boost::program_options::options_description hidden("Hidden options");
   hidden.add_options()
     ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
     ("pruning,j", boost::program_options::value<uint32_t>(&c.graphPruning)->default_value(1000), "graph pruning cutoff")
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
