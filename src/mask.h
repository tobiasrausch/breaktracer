#ifndef MASK_H
#define MASK_H

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


  struct MaskConfig {
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
    boost::filesystem::path insseq;
    std::vector<std::string> sampleName;
  };
  
  template<typename TConfig>
  inline int32_t
  runMask(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("tracer.prof");
#endif
    
   
#ifdef PROFILE
    ProfilerStop();
#endif

    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    
    return 0;
  }

  int mask(int argc, char **argv) {
   MaskConfig c;
   
   // Parameter
   std::string mode;
   std::string instag;
   boost::program_options::options_description generic("Generic options");
   generic.add_options()
     ("help,?", "show help message")
     ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output genome mask in FASTA format")
     ;
   
   boost::program_options::options_description brin("Sequence options");
   brin.add_options()
     ("seedlen,s", boost::program_options::value<int32_t>(&c.minSeedAlign)->default_value(130), "min. seed length")
     ("pctid,i", boost::program_options::value<float>(&c.pctThres)->default_value(0.9), "min. percent identity")
     ("instag,n", boost::program_options::value<std::string>(&instag)->default_value("L1"), "Type of sequence [ALU|L1|SVA|NUMT]")
     ("insseq,e", boost::program_options::value<boost::filesystem::path>(&c.insseq), "FASTA with insertion sequence [overrides -n]")
     ;

   boost::program_options::options_description hidden("Hidden options");
   hidden.add_options()
     ("input-file", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
     ;
   
   boost::program_options::positional_options_description pos_args;
   pos_args.add("input-file", -1);
   
   boost::program_options::options_description cmdline_options;
   cmdline_options.add(generic).add(brin).add(hidden);
   boost::program_options::options_description visible_options;
   visible_options.add(generic).add(brin);
   boost::program_options::variables_map vm;
   boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
   boost::program_options::notify(vm);
   
   // Check command line arguments
   if ((vm.count("help")) || (!vm.count("input-file"))) {
     std::cerr << std::endl;
     std::cerr << "Usage: breaktracer " << argv[0] << " [OPTIONS] <ref.fa>" << std::endl;
     std::cerr << visible_options << "\n";
     return 0;
   }

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
   
   // Check outfile
   if (!vm.count("outfile")) c.outfile = "-";
   else {
     if (c.outfile.string() != "-") {
       if (!_outfileValid(c.outfile)) return 1;
     }
   }

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
   
   return runMask(c);
 }

}

#endif
