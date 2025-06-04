#include "L1Trigger/L1TTrackMatch/interface/TkTripletEmuAlgo.h"

namespace l1ttripletemu {

    using WMass_t = ap_ufixed<kMassSize, kMassMagSize, AP_RND_CONV, AP_SAT>;
    
    WMass_t doubleToUFixedMass(double emuWMass) {
        // 1) Compute F = number of fractional bits
        emuWMass = 16.0;
        constexpr int W = kMassSize;
        constexpr int I = kMassMagSize;
        constexpr int F = W - I;             // = 6
        constexpr double LSB = 1.0 / (1 << F); // = 1/64 = 0.015625

        // 2) Compute unsigned magnitude in “LSB units” (floor)
        // If emulation mass is negative treat it as zero (ap_ufixed is unsigned).
        if (emuWMass <= 0.0) {
            return WMass_t(0);
        }
        unsigned long long mag = static_cast<unsigned long long>( std::floor(emuWMass / LSB) );

        // 3) Clamp to the maximum code = 2^W − 1
        constexpr unsigned long long MAX_CODE = (1ULL << W) - 1; // = 0x7FFFFF
        if (mag > MAX_CODE) {
            mag = MAX_CODE;
        }

        // 4) Pack that integer “mag” directly into the WMass_t:
        // Assigning an unsigned long long to an ap_ufixed<23,17> writes it as the raw bit-pattern.
        return WMass_t(mag);
    }
    
}  // namespace l1ttripletemu