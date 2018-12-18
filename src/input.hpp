#ifndef __SRC_INPUT_H__
#define __SRC_INPUT_H__

#include <string>
#include <vector>

class input{

 public:
  input(const int, const char **);
  ~input();

  std::string probandBamPath_;
  std::string contigFastaPath_;
  std::string contigBamPath_;
  std::string mobileElementFastqPath_;
  std::string mobileElementIndexPath_;
  std::string referencePath_;
  std::string referenceIndexPath_;
  std::string vcfOutPath_;

  const std::vector<std::string> parentBamPaths_;

 private:
  
  const int argc_;
  const char ** argv_;

  void parseArgs();

};


#endif // __SRC_INPUT_H__
