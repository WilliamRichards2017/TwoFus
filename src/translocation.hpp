#ifndef __SRC_TRANSLOCATION_HPP__
#define __SRC_TRANSLOCATION_HPP__

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "input.hpp"

class translocation{
public:
  translocation();
  translocation(const std::vector<BamTools::BamAlignment> &, const std::map<std::string, std::vector<BamTools::BamAlignment> > &, const input &);
  translocation(const translocation &);
  ~translocation();


  const std::map<std::string, std::vector<BamTools::BamAlignment> > & getSAMap();

  const BamTools::BamAlignment & getLeftContig();
  const BamTools::BamAlignment & getRightContig();

private:

  input i_;
  std::vector<BamTools::BamAlignment> groupedContigs_;
  std::map<std::string, std::vector<BamTools::BamAlignment> > SAMap_;

  BamTools::BamAlignment leftContig_;
  BamTools::BamAlignment rightContig_;


  void populateLeftAndRightContigs();
  const bool isTrans();
  std::vector<BamTools::BamAlignment> pullAllReadsWithName(const std::string &);

};

#endif //__SRC_TRANSLOCATION_HPP__
