#ifndef __SRC_INPUT_H__
#define __SRC_INPUT_H__

#include <string>
#include <vector>

class input{

 public:
  input();
  input(const input &);
  input(int, char **);
  ~input();

  std::string probandBamPath_;
  std::string probandMutBamPath_;
  std::string contigFastqPath_;
  std::string contigBamPath_;
  std::string mobileElementFastaPath_;
  std::string mobileElementIndexPath_;
  std::string referencePath_;
  std::string referenceIndexPath_;
  std::string vcfOutPath_;
  std::string kmerPath_;
  std::string hashListPath_;
  std::vector<std::string> parentBamPaths_;

  std::string probandRefPath_;
  std::string probandAltPath_;
  std::vector<std::string> parentRefPaths_;
  std::vector<std::string> parentAltPaths_;

  void printArgs();

 private:
  
  int argc_;
  char ** argv_;

  void parseArgs();
  void populateKmerPaths();

};

#endif // __SRC_INPUT_H__
