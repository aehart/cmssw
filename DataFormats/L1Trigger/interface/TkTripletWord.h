#ifndef FIRMWARE_TkTripletWord_h
#define FIRMWARE_TkTripletWord_h

// Class to store the 128-bit TkTriplet word for the L1TrackTriggerMatch.
// Original author: George Karathanasis (Dec 2023)

#include <vector>
#include <ap_int.h>
#include <cassert>
#include <cmath>
#include <bitset>
#include <string>
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"

namespace l1t {

  class TkTripletWord {
  public:
    // ----------constants, enums and typedefs ---------
    static constexpr double MAX_MASS = 1000.;
    static constexpr double MAX_ETA = 8.;
    static constexpr double MAX_CHARGE = 3.;
    static constexpr double MAX_Z0 = 25.;

    enum TkTripletBitWidths {
      // The sizes of the triplet word components and total word size
      kValidSize = 1,           // Width of the valid bit
      kMassSize = 23,           // Width of the entire triplet mass
      kMassMagSize = 17,        // Width of the triplet mass magnitude
      kUnassignedSize = 40,
      kTkTripletWordSize = kValidSize + kMassSize + kUnassignedSize,
    };

    enum TkTripletBitLocations {
      // The location of the least significant bit (LSB) and most significant bit (MSB) in the triplet word for different fields
      kValidLSB = 0,
      kValidMSB = kValidLSB + TkTripletBitWidths::kValidSize - 1,
      kMassLSB = kValidMSB + 1,
      kMassMSB = kMassLSB + TkTripletBitWidths::kMassSize - 1,
      kUnassignedLSB = kMassMSB + 1,
      kUnassignedMSB = kUnassignedLSB + TkTripletBitWidths::kUnassignedSize - 1,
    };

    // ap parameter types
    typedef ap_uint<kValidSize> valid_t;                                             //valid
    typedef ap_ufixed<kMassSize, kMassMagSize, AP_RND_CONV, AP_SAT> WMass_t;          //triplet mass
    typedef ap_uint<TkTripletBitWidths::kUnassignedSize> unassigned_t;
    typedef std::bitset<TkTripletBitWidths::kTkTripletWordSize> tktripletword_bs_t;
    typedef ap_uint<TkTripletBitWidths::kTkTripletWordSize> tktripletword_t;

  public:
    // ----------Constructors --------------------------
    TkTripletWord() {}
    TkTripletWord(valid_t valid,
                  WMass_t mass,
                  unassigned_t unassigned);

    ~TkTripletWord() {}

    // ----------copy constructor ----------------------
    TkTripletWord(const TkTripletWord& word) { tkTripletWord_ = word.tkTripletWord_; }

    // ----------operators -----------------------------
    TkTripletWord& operator=(const TkTripletWord& word) {
      tkTripletWord_ = word.tkTripletWord_;
      return *this;
    }

    // ----------member functions (getters) ------------
    // These functions return arbitarary precision words (lists of bits) for each quantity
    valid_t validWord() const {
      return tkTripletWord()(TkTripletBitLocations::kValidMSB, TkTripletBitLocations::kValidLSB);
    }
    WMass_t massWord() const {
      WMass_t ret;
      ret.V = tkTripletWord()(TkTripletBitLocations::kMassMSB, TkTripletBitLocations::kMassLSB);
      return ret;
    }
    unassigned_t unassignedWord() const {
      return tkTripletWord()(TkTripletBitLocations::kUnassignedMSB, TkTripletBitLocations::kUnassignedLSB);
    }
    tktripletword_t tkTripletWord() const { return tktripletword_t(tkTripletWord_.to_string().c_str(), 2); }

    // These functions return the packed bits in integer format for each quantity
    // Signed quantities have the sign enconded in the left-most bit.
    unsigned int validBits() const { return validWord().to_uint(); }
    unsigned int massBits() const { return massWord().to_uint(); }
    unsigned int unassignedBits() const { return unassignedWord().to_uint(); }

    // These functions return the unpacked and converted values
    // These functions return real numbers converted from the digitized quantities by unpacking the 64-bit vertex word
    bool valid() const { return validWord().to_bool(); }
    float mass() const {
      return unpackSignedValue(
          massWord(), TkTripletBitWidths::kMassSize, MAX_MASS / (1 << TkTripletBitWidths::kMassSize));
    }
    unsigned int unassigned() const { return unassignedWord().to_uint(); }

    // ----------member functions (setters) ------------
    void setTkTripletWord(valid_t valid,
                          WMass_t wMass,
                          unassigned_t unassigned);

    template <class packVarType>
    inline void packIntoWord(unsigned int& currentOffset, unsigned int wordChunkSize, packVarType& packVar);

  private:
    // ----------private member functions --------------
    double unpackSignedValue(unsigned int bits, unsigned int nBits, double lsb) const {
      int isign = 1;
      unsigned int digitized_maximum = (1 << nBits) - 1;
      if (bits & (1 << (nBits - 1))) {  // check the sign
        isign = -1;
        bits = (1 << (nBits + 1)) - bits;  // if negative, flip everything for two's complement encoding
      }
      return (double(bits & digitized_maximum) + 0.5) * lsb * isign;
    }

    // ----------member data ---------------------------
    tktripletword_bs_t tkTripletWord_;
  };

  typedef std::vector<l1t::TkTripletWord> TkTripletWordCollection;

}  // namespace l1t

#endif
