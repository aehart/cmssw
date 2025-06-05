#ifndef L1Trigger_L1TTrackmatch_TkTripletEmuAlgo_HH
#define L1Trigger_L1TTrackmatch_TkTripletEmuAlgo_HH

#include <ap_int.h>
#include <ap_fixed.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1Trigger/interface/TkTripletWord.h"

// Namespace that defines constants and types used by the TkTriplet collection

namespace l1ttripletemu {

    const unsigned int kValidSize{1};
    const unsigned int kMassSize{23};
    const unsigned int kMassMagSize{17};
    const unsigned int kUnassignedSize{64 - (kMassSize + kValidSize)};

    typedef ap_ufixed<kMassSize, kMassMagSize, AP_RND_CONV, AP_SAT> WMass_t;

    // Output definition as per interface document, only used when creating output format
    const double kMaxMass = 1 << kMassMagSize;
    const double kStepMass = kMaxMass / (1 << kMassSize);

    struct TkTriplet {
        WMass_t mW;
    };

}  // namespace l1ttripletemu
#endif