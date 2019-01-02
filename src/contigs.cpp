#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "minimap.h"
#include "/uufs/chpc.utah.edu/common/home/u0401321/TwoFus/bin/externals/minimap2/src/minimap2_project/kseq.h"
#include "zlib.h"


#include "contigs.hpp"
#include "input.hpp"
#include "util.hpp"


KSEQ_INIT(gzFile, gzread)


void contigs::findAllContigs(){
  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(i_.contigBamPath_)){
    std::cout << "Could not open bam file in contigs::findAllContigs() for " << i_.contigBamPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << i_.contigBamPath_ << " could not be opened in contigs::findAllContigs()" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  while(reader.GetNextAlignment(al)){
    contigVec_.push_back(al);
  }

}


void contigs::findInsertionContigs(){
}

void contigs::alignContigsToMEList(){
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;
  int n_threads = 3;

  mm_verbose = 3; // print to std out
  mm_set_opt(0, &iopt, &mopt); //initialize alignment parameters to default
  mopt.flag |= MM_F_CIGAR; // perform alignment                                                                                                                                                                                               

  gzFile f = gzopen(i_.contigFastqPath_.c_str(), "r");
  assert(f);
  kseq_t *ks = kseq_init(f);

  // open index reader
  if(! util::fileExists(i_.contigFastqPath_)){
    std::cout << "Could not open alu fasta file " << i_.contigFastqPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    exit (EXIT_FAILURE);
  }

  mm_idx_reader_t *r = mm_idx_reader_open(i_.contigFastqPath_.c_str(), &iopt, 0);
  mm_idx_t *mi;
  while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
    mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence
    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread                                                                                                                              

    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence                                                                                                                                                           
      //std::cout << "looping through kseq_reads" << std::endl; 

      mm_reg1_t *reg;
      int j, i, n_reg;

      reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query

      for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	mm_reg1_t *r = &reg[j];
	contigsAlignedToMEList_.insert({mi->seq[r->rid].name, int(r->mapq)});
	free(r->p);
	//std::cout << "found alu hit for contig: " << ks->name.s << std::endl;
      }
      free(reg);
    }
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
  }
  mm_idx_reader_close(r); 
  kseq_destroy(ks); // close the query file                                                                                                                                                                                                   
  gzclose(f);
}

void contigs::findMobileElementContigs(){
  contigs::alignContigsToMEList();
  std::cout << "aligned: " << contigsAlignedToMEList_.size() << " contigs to ME list" << std::endl;
  
  for(const auto & g : groupedContigsVec_){
    for(const auto & c: g){
      /*if(contigs::contigAligns(c).second != -1){
	groupedMobileElementContigVec_.push_back(g);
	break;
	}*/
    }
  }
}

void contigs::findTranslocationContigs(){
}

void contigs::groupNearbyContigs(){

  BamTools::BamAlignment previousContig;
  int32_t maxDist = 100;

  for(const auto & c : contigVec_){
    groupedContigs g;
    if(previousContig.RefID == c.RefID && std::abs(c.Position-previousContig.Position) < maxDist){
      g.push_back(c);
      g.push_back(previousContig);
    }
    else{
      g.push_back(c);
    }
    previousContig = c;
    groupedContigsVec_.push_back(g);
  }

}


contigs::contigs(const input & i) : i_(i){
  contigs::findAllContigs();
  contigs::groupNearbyContigs();
  contigs::findMobileElementContigs();
  
}

contigs::~contigs(){
}
