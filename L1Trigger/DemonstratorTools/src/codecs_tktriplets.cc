#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/DemonstratorTools/interface/codecs/tktriplets.h"

namespace l1t::demo::codecs {

  ap_uint<64> encodeTriplet(const l1t::TkTripletWord& t) {
    ap_uint<1> valid = (t.validWord());
    ap_uint<64 - (l1t::TkTripletWord::kMassSize + 1)> unassigned = 0;
    ap_uint<64> tripletWord = (unassigned, t.massWord(), valid);
    return tripletWord;
  }

  // Encodes etsum collection onto 1 output link (or 2+ depending on size of word)
  std::array<std::vector<ap_uint<64>>, 1> encodeTriplets(const edm::View<l1t::TkTripletWord>& triplets) {
    std::vector<ap_uint<64>> tripletWords;

    for (const auto& triplet : triplets)
      tripletWords.push_back(encodeTriplet(triplet));

    std::array<std::vector<ap_uint<64>>, 1> linkData;

    for (size_t i = 0; i < linkData.size(); i++) {
      // Pad etsum vectors -> full packet length (48 frames, but only 1 etsum max)
      tripletWords.resize(1, 0);
      linkData.at(i) = tripletWords;
    }

    return linkData;
  }

} // namespace l1t::demo::codecs

  