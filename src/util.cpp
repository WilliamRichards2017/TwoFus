#include "util.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <unistd.h>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

const bool util::breakpointOverlapsPeak(const std::pair<int32_t, int32_t> & indices, const int32_t & breakpoint){
  return(breakpoint >= indices.first-5 and breakpoint <= indices.second+5);
}

void replaceAll( string &s, const string &search, const string &replace ) {
  for( size_t pos = 0; ; pos += replace.length() ) {
    // Locate the substring to replace
    pos = s.find( search, pos );
    if( pos == string::npos ) break;
    // Replace by erasing and inserting
    s.erase( pos, search.length() );
    s.insert( pos, replace );
  }
}

const std::string util::replaceDoubleGenerator(str){

  std::string replace = replaceAll(str, "generator.generator", "generator");

  std::cout << "replace: " << replace << std::endl;
}

const maxPeak util::getMaxPeak(const std::vector<std::pair<int32_t, int32_t> > & peakVec, const BamTools::BamAlignment & al){
  maxPeak mp;
  mp.indices = std::make_pair(-1,-1);
  mp.value = -34;
  

  for(const auto & p : peakVec){
    if(int(al.Qualities[p.first])-33 > mp.value){
      mp.indices = p;
      mp.value = int(al.Qualities[p.first])-33;
    }
  }
  return mp;
}

const std::vector<int32_t> util::getPeakVector(const BamTools::BamAlignment & al){
  
  std::vector<int32_t> peakVec;
  //std::cout << "peakVector : ";
  for(auto c : al.Qualities){
    peakVec.push_back(int(c)-33); // ascii conversion
    //std::cout << peakVec.back() << ", ";
  }
  //std::cout << std::endl;
  return peakVec;  
}

const std::vector<std::pair<int32_t, int32_t> > util::getPeaks(const BamTools::BamAlignment & al) {

  
  auto amp = util::getPeakVector(al);
  int wideStart = -1;                 // The start of any current wide peak
  
  int grad = -1;                      // Sign of gradient (almost)
  //    =  1 for increasing
  //    =  0 for level AND PREVIOUSLY INCREASING (so potential wide peak)
  //    = -1 for decreasing OR level, but previously decreasing
  // A sharp peak is identified by grad=1 -> grad=-1
  // A wide  peak is identified by grad=0 -> grad=-1

  std::vector<std::pair<int32_t, int32_t> > peakCoords;

  for (int i = 0; i < amp.size() - 1; i++) {
    if(amp[i+1] < amp[i]){    
      if(grad == 1){
	//std::cout << "Sharp peak of " << amp[i] << " at i = " << i << '\n';
	peakCoords.push_back(std::make_pair(i,i));
      }
      else if(grad == 0){
	peakCoords.push_back(std::make_pair(wideStart,i));
	//std::cout << "Wide peak of " << amp[i] << " from i = " << wideStart << " to " << i << '\n';
      }
      grad = -1;
    }
    else if(amp[i+1] == amp[i]){   // Check for start of a wide peak
      if(grad == 1){
	wideStart = i;
	grad = 0;
      }
    }
    else{
      grad = 1;
    }
  }
  return peakCoords;
}


const bool util::breakpointHasSupport(const BamTools::BamAlignment &){

  std::string variant;
  std::vector<std::string> kmers;

}


const std::vector<BamTools::BamAlignment> util::filterOutPrimaryAlignment(const BamTools::BamAlignment & pAl, const std::vector<BamTools::BamAlignment> & allAl){
  std::vector<BamTools::BamAlignment> sa;

  for(const auto & al : allAl){
    if(al.Position != pAl.Position){
      sa.push_back(al);
    }
  }
  return sa;
}

const std::vector<BamTools::BamAlignment>  util::pullAllReadsWithName(const std::string & readName, 
								      const std::map<std::string, std::vector<BamTools::BamAlignment> > & SAMap){
  std::vector<BamTools::BamAlignment> nullVec;
  auto it = SAMap.find(readName);

  if(it != SAMap.end()){
    return it->second;
  }
  return nullVec;
}

const std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > util::checkIfSecondariesAreNearby(const std::vector<std::pair<BamTools::BamAlignment, std::vector<BamTools::BamAlignment> > > & contigs){
  
  std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment> > primaryContigs;


  for(unsigned i = 0; i < contigs.size()-1; ++i){
    auto pc1 = contigs[i].first;
    auto pc2 = contigs[i+1].second;

    auto sa1 = contigs[i].second;
    auto sa2 = contigs[i+1].second;


    if(sa1.size() > 0 and sa2.size() < 4){
      for(const auto s1 : sa1){
	for(const auto s2 : sa2){
	  if(util::isNearby(s1, s2)){
	    std::cout << "Found nearby supplemental alignments" << std::endl;	  

	    auto p1 = std::make_pair(pc1, s1);
	    auto p2 = std::make_pair(pc1, s2);
	    primaryContigs.push_back(p1);
	    primaryContigs.push_back(p2);
	    return primaryContigs;
	  }
	}
      }
    }
  }
  return primaryContigs;
}

const bool util::checkClipsConverge(const BamTools::BamAlignment & al1, const BamTools::BamAlignment & al2){
  bool rightBound1 = false;
  bool rightBound2 = false;

  BamTools::BamAlignment lAl;
  BamTools::BamAlignment rAl;

  clipCoords cc1 = clipCoords{al1};
  clipCoords cc2 = clipCoords{al2};

  if(cc1.leftPos_ + cc1.globalOffset_ < cc2.leftPos_ + cc2.globalOffset_){
    rAl = al1;
    lAl = al2;
    std::cout << "rAl.Name " << rAl.Name << std::endl;
    std::cout << "rAl.leftPos " << cc1.leftPos_ + cc1.globalOffset_ << std::endl;

    std::cout << "lAl.Name " << lAl.Name << std::endl;
    std::cout << "lAl.leftPos " << cc2.leftPos_ + cc2.globalOffset_ << std::endl;
  }
  else{
    rAl = al2;
    lAl = al1;

    std::cout << "rAl.Name " << rAl.Name << std::endl;
    std::cout << "rAl.leftPos " << cc2.leftPos_ + cc2.globalOffset_ << std::endl;
    std::cout << "lAl.Name " << lAl.Name << std::endl;
    std::cout << "lAl.leftPos " << cc1.leftPos_ + cc1.globalOffset_ << std::endl;
  }
  
  std::vector<BamTools::CigarOp> cig1 = lAl.CigarData;
  std::vector<BamTools::CigarOp> cig2 = rAl.CigarData;
  
  if(cig1[0].Type == 'S'){
    rightBound1 = true;
  }
  if(cig2[0].Type == 'S'){
    rightBound2 = true;
  }


  std::cout << "rightBound1 " << rightBound1 << std::endl;
  std::cout << "rightBound2 " << rightBound2 << std::endl;
    return !rightBound1 and rightBound2;
}

const bool util::isNearby(const BamTools::BamAlignment & al1, const BamTools::BamAlignment & al2){
  int32_t maxDist = 1000;
  if(al1.RefID == al2.RefID and std::abs(al1.Position-al2.Position) < maxDist){
    return true;
  }
  return false;
}

const std::pair<BamTools::BamAlignment, BamTools::BamAlignment> util::findLeftAndRightContigs(const std::vector<BamTools::BamAlignment> & groupedContigs){
  BamTools::BamAlignment leftContig;
  BamTools::BamAlignment rightContig;
  
  if(groupedContigs.size() < 2){
    std::cerr << "Inside util::populateLeftAndRightContigs() without atleast two contigs..." << std::endl;
    std::cerr << "Exiting run with non-zero status..." << std::endl;
    exit (EXIT_FAILURE);
  }

  else if(groupedContigs.size() > 2){
    std::cerr << "Warning: group has more than two contigs in util::populateLeftAndRightContigs()" << std::endl;
    std::cerr << "Proceding using left most and right most contig" << std::endl;
  }

  int32_t leftMostPos = std::numeric_limits<int32_t>::max();
  int32_t rightMostPos = -1;
  int32_t leftIndex = -1;
  int32_t rightIndex = -1;
  for(unsigned i = 0; i < groupedContigs.size(); ++i){

    //std::cout << "populating left and right contigs for contig: " << groupedContigs[i].Name << " at position " << groupedContigs[i].Position <<  std::endl;
    //std::cout << "comparing if " << groupedContigs[i].Position << " < " << leftMostPos << std::endl;
    //std::cout << "comparing if " << groupedContigs[i].Position << " > " << rightMostPos << std::endl;
                                                                            
    if(groupedContigs[i].Position < leftMostPos){
      leftMostPos = groupedContigs[i].Position;
      leftIndex = i;
    }
    if(groupedContigs[i].Position > rightMostPos){
      rightMostPos = groupedContigs[i].Position;
      rightIndex = i;
    }
  }

  if(leftIndex != -1){
    leftContig = groupedContigs[leftIndex];
    //std::cout << "leftContig_.Name is: " << leftContig_.Name << std::endl;
  }
  if(rightIndex != -1){
    rightContig = groupedContigs[rightIndex];
    //std::cout << "rightContig_.Name is: " << rightContig_.Name << std::endl;
  }
  return std::make_pair(leftContig, rightContig);
}


const int32_t util::countMinKmerDepth(const std::vector<std::pair<std::string, int32_t> > & kmers){
  std::vector<int32_t> kmerCounts;

  for(const auto & k : kmers){
    if(k.second > 0){
      kmerCounts.push_back(k.second);
    }
  }

  if(kmerCounts.size() == 0){
    return 0;
  }

  return *std::min_element(kmerCounts.begin(), kmerCounts.end());
}


const std::unordered_map<std::string, int32_t> util::countKmersFromJhash(const std::string & jhashPath, const std::vector<std::string> & kmers){

  std::string jellyfishPath = "../bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish";
  
  std::unordered_map<std::string, int32_t> ret;
  for (const auto & kmer : kmers){

    std::string queryCommand = jellyfishPath + " query " + jhashPath + " " + kmer;
    //std::cout << "executing command: " << queryCommand << std::endl;

    std::string queryOutput = util::exec(queryCommand.c_str());
    //std::cout << "queryOutput is " << queryOutput << std::endl;
    


    std::istringstream iss(queryOutput);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

    if(kmerCount.size() == 2){
      ret.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
    }
  }
  return ret;
}

const std::vector<std::pair<std::string, int32_t> > util::countKmersFromText(const std::string & textPath, const std::vector<std::string> & kmers){
  std::ifstream file(textPath);
  std::string line;

  std::vector<std::pair<std::string, int32_t> > kmerCounts;
  std::map<std::string, int32_t> kmerMap;

  while(std::getline(file, line)){
    std::istringstream iss(line);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
    if(kmerCount.size() == 2){
      kmerMap.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
    }
  }
  for(auto k : kmers){
    //std::cout << "looking up kmer: " << k << std::endl;
    auto it = kmerMap.find(k);
    auto revIt = kmerMap.find(util::revComp(k));
    if(it != kmerMap.end()){
      kmerCounts.push_back(std::make_pair(it->first, it->second));
      //std::cout << "found kmer in map" << std::endl;
    }
    else if(revIt != kmerMap.end()){
      kmerCounts.push_back(std::make_pair(k, revIt->second));
      //std::cout << "found revComp(kmer) in map" << std::endl;
    }
    else{
      //std::cout << "else statement" << std::endl;
      kmerCounts.push_back(std::make_pair(k, 0));
    }
  }
  //std::cout << "Returning kmerCounts with size of" << kmerCounts.size() << std::endl;
  return kmerCounts;
}

const std::vector<std::string> util::kmerize(const std::string & sequence, const int32_t & kmerSize){
  int32_t kmercount = 0;
  std::vector<std::string> kmers;

  while(kmercount + kmerSize <= sequence.length()){
    std::string kmer = sequence.substr(kmercount, kmerSize);
    if(kmer.size() == kmerSize){
      kmers.push_back(kmer);
    }
    ++kmercount;
  }
  return kmers;
}

const std::string util::getChromosomeFromRefID(const int32_t & id, const std::vector<BamTools::RefData> & refData){
  std::string ret = "";
  if(id == -1) {
    ret = "unmapped";
    return ret;
  }
  ret = refData[id].RefName;
  return ret;
}

const std::string util::pullRefSequenceFromRegion(const breakpoint & leftBP, const int32_t & rightPos,
						  const std::string & refPath, const std::vector<BamTools::RefData> & refData){

  std::string fastahackPath = "../bin/externals/fastahack/src/fastahack_project-build/tools/fastahack";

  std::string cmd = fastahackPath + " -r " + util::getChromosomeFromRefID(leftBP.refID, refData) + ':' + std::to_string(leftBP.position) + ".." + std::to_string(rightPos) + ' ' + refPath;

  std::cout << "Executing command: " << cmd << std::endl;
  
  return util::exec(cmd.c_str());
}

std::vector<BamTools::RefData> util::populateRefData(const std::string & bamPath){
  BamTools::BamReader reader = util::openBamFile(bamPath);
  return reader.GetReferenceData();
}

const std::vector<std::string> util::split(const std::string & line, const char delim)
{
  std::vector<std::string> tokens;
  std::stringstream lineStream(line);
  std::string token;
  while(getline(lineStream, token, delim)){
    tokens.push_back(token);
  }
  return tokens;
}

const float util::calculateStrandBiasFromContigName(const std::string & contigName){

    std::vector<std::string> tokenizedName = util::split(contigName, ':');

  std::pair<std::pair<int32_t, int32_t>, float> f;

  if(tokenizedName.size() > 2){
    f.first = std::make_pair(std::atoi(tokenizedName[1].c_str()), std::atoi(tokenizedName[2].c_str()));
    f.second = float(f.first.first)/(float(f.first.first) + float(f.first.second));

  }
  return f.second;
}

const float util::calculateStrandBiasFromContigNames(const std::vector<std::string> & contigNames){
  float forwardCount = 0;
  float reverseCount = 0;
  for(const auto & cn : contigNames){
    std::vector<std::string> tokenizedName = util::split(cn, ':');
    forwardCount += std::atoi(tokenizedName[1].c_str());
    reverseCount += std::atoi(tokenizedName[2].c_str());
    
    std::cout << "forwardCount: " << forwardCount << std::endl;
    std::cout << "reverseCount: " << reverseCount << std::endl;
  }

  std::cout << "SB from contigNames is: " << forwardCount/(forwardCount+reverseCount) << std::endl;
  return forwardCount/(forwardCount+reverseCount);
}

BamTools::BamReader util::openBamFile(const std::string & bamPath){

  BamTools::BamReader reader;  

  if(!reader.Open(bamPath)){
    std::cout << "could not open bamPath " << bamPath << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);
  }
  
  reader.LocateIndex();

  if(!reader.HasIndex()){
    std::cout << "Index for " << bamPath << " could not be located" << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit(EXIT_FAILURE);
  }
  return reader;
}

const bool util::fileExists(const std::string & Filename){
  return access(Filename.c_str(), 0) == 0;
}

const bool util::isReadLeftBound(const std::vector<BamTools::CigarOp> & cigOps){
  if(cigOps[0].Type == 'S'){
    return true;
  }
  return false;
}

const std::vector<int32_t> util::getInsertionVec(const BamTools::BamAlignment & al){
  const std::vector<BamTools::CigarOp> cig = al.CigarData;
  std::vector<int32_t> insertionVec;
  int32_t indel = 0;
  for(auto c : cig){
    if(c.Type =='S'){
      insertionVec.push_back(indel);
      indel = 0;
    }
    else if(c.Type == 'D'){
      indel -= c.Length;
    }
  }
  return insertionVec;
}

const std::vector<std::string> util::getClipSeqs(const BamTools::BamAlignment & al){
  std::vector<std::string> clipSeqs;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  al.GetSoftClips(clipSizes, readPositions, genomePositions);

  const std::vector<int32_t> insertionVec = util::getInsertionVec(al);
  for(int i = 0; i < readPositions.size(); ++i){
    if(util::isReadLeftBound(al.CigarData)){
      clipSeqs.push_back(al.QueryBases.substr(0, clipSizes[i]));
    }
    
    else{  
      //std::cout << "Clipped seq for read is: " << al.QueryBases.substr(readPositions[i], clipSizes[i]) << std::endl;
      clipSeqs.push_back(al.QueryBases.substr(readPositions[i]+insertionVec[i], clipSizes[i]));
    }
  }
  return clipSeqs;
}

std::string util::exec(char const* cmd) {
  char buffer[128];
  std::string result = "";
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
      if (fgets(buffer, 128, pipe) != NULL)
	result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  return result;
}

const std::string util::revComp (std::string sequence){
  auto lambda = [](const char c) {
    switch (c) {
    case 'A':
    case 'a':
    return 'T';
    case 'G':
    case 'g':
    return 'C';
    case 'C':
    case 'c':
    return 'G';
    case 'T':
    case 't':
    return 'A';
    case 'N':
    case 'n': 
    return 'N';
    default:
    std::cout << "Invalid nucleotide in revComp: " << c << std::endl;
    }
  };

  std::transform(sequence.cbegin(), sequence.cend(), sequence.begin(), lambda);
  return sequence;
}
