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
  std::string contigFastqPath_;
  std::string contigBamPath_;
  std::string mobileElementFastaPath_;
  std::string mobileElementIndexPath_;
  std::string referencePath_;
  std::string referenceIndexPath_;
  std::string vcfOutPath_;

  std::vector<std::string> parentBamPaths_;

  void printArgs();

 private:
  
  int argc_;
  char ** argv_;

  void parseArgs();

};


#endif // __SRC_INPUT_H__
