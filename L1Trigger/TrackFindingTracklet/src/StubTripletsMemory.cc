#include "L1Trigger/TrackFindingTracklet/interface/StubTripletsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/L1TStub.h"
#include "L1Trigger/TrackFindingTracklet/interface/Stub.h"
#include <iomanip>
#include <filesystem>

using namespace std;
using namespace trklet;

StubTripletsMemory::StubTripletsMemory(string name, Settings const& settings) : MemoryBase(name, settings) {}

void StubTripletsMemory::writeST(bool first, unsigned int iSector) {
  iSector_ = iSector;

  const string dirSP = settings_.memPath() + "StubPairs/";
  openFile(first, dirSP, "StubTriplets_");

  for (unsigned int j = 0; j < stubs1_.size(); j++) {
    string stub1index = stubs1_[j]->stubindex().str();
    string stub2index = stubs2_[j]->stubindex().str();
    string stub3index = stubs3_[j]->stubindex().str();
    out_ << hexstr(j) << " " << stub1index << "|" << stub2index << "|" << stub3index << endl;
  }
  out_.close();
}
