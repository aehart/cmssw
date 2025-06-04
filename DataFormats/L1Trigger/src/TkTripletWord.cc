// Class to store the 128-bit TkTriplet word for L1 Track Trigger.
// Author: George Krathanasis, CU Boulder (December 2023)

#include "DataFormats/L1Trigger/interface/TkTripletWord.h"

namespace l1t {
  TkTripletWord::TkTripletWord(valid_t valid,
                               WMass_t mass,
                               unassigned_t unassigned) {
    setTkTripletWord(valid, mass, unassigned);
  }

  template <class packVarType>
  inline void TkTripletWord::packIntoWord(unsigned int& currentOffset,
                                          unsigned int wordChunkSize,
                                          packVarType& packVar) {
    for (unsigned int b = currentOffset; b < (currentOffset + wordChunkSize); ++b) {
      tkTripletWord_.set(b, packVar[b - currentOffset]);
    }
    currentOffset += wordChunkSize;
  }

  void TkTripletWord::setTkTripletWord(valid_t valid,
                                       WMass_t mass,
                                       unassigned_t unassigned) {
    // pack the TkTriplet word
    unsigned int offset = 0;
    packIntoWord(offset, TkTripletBitWidths::kValidSize, valid);
    packIntoWord(offset, TkTripletBitWidths::kMassSize, mass);
  }

}  //namespace l1t
