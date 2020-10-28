# Poor-man enum class with string conversion
class _Enum:
    def __init__(self, **values):
        import six
        self._reverse = {}
        for key, value in six.iteritems(values):
            setattr(self, key, value)
            if value in self._reverse:
                raise Exception("Value %s is already used for a key %s, tried to re-add it for key %s" % (value, self._reverse[value], key))
            self._reverse[value] = key

    def toString(self, val):
        return self._reverse[val]

SubDet = _Enum(
    BPix = 1,
    FPix = 2,
    TIB = 3,
    TID = 4,
    TOB = 5,
    TEC = 6
)

# Needs to be kept consistent with
# DataFormats/TrackReco/interface/TrackBase.h
Algo = _Enum(
    undefAlgorithm = 0, ctf = 1,
    duplicateMerge = 2, cosmics = 3,
    initialStep = 4,
    lowPtTripletStep = 5,
    pixelPairStep = 6,
    detachedTripletStep = 7,
    mixedTripletStep = 8,
    pixelLessStep = 9,
    tobTecStep = 10,
    jetCoreRegionalStep = 11,
    conversionStep = 12,
    muonSeededStepInOut = 13,
    muonSeededStepOutIn = 14,
    exoRegionalStep = 15,
    outInEcalSeededConv = 16, inOutEcalSeededConv = 17,
    nuclInter = 18,
    standAloneMuon = 19, globalMuon = 20, cosmicStandAloneMuon = 21, cosmicGlobalMuon = 22,
    # Phase1
    highPtTripletStep = 23, lowPtQuadStep = 24, detachedQuadStep = 25,
    reservedForUpgrades1 = 26, reservedForUpgrades2 = 27,
    bTagGhostTracks = 28,
    beamhalo = 29,
    gsf = 30,
    # HLT algo name
    hltPixel = 31,
    # steps used by PF
    hltIter0 = 32,
    hltIter1 = 33,
    hltIter2 = 34,
    hltIter3 = 35,
    hltIter4 = 36,
    # steps used by all other objects @HLT
    hltIterX = 37,
    # steps used by HI muon regional iterative tracking
    hiRegitMuInitialStep = 38,
    hiRegitMuLowPtTripletStep = 39,
    hiRegitMuPixelPairStep = 40,
    hiRegitMuDetachedTripletStep = 41,
    hiRegitMuMixedTripletStep = 42,
    hiRegitMuPixelLessStep = 43,
    hiRegitMuTobTecStep = 44,
    hiRegitMuMuonSeededStepInOut = 45,
    hiRegitMuMuonSeededStepOutIn = 46,
    algoSize = 47
)

# Needs to kept consistent with
# DataFormats/TrackReco/interface/TrajectoryStopReasons.h
StopReason = _Enum(
  UNINITIALIZED = 0,
  MAX_HITS = 1,
  MAX_LOST_HITS = 2,
  MAX_CONSECUTIVE_LOST_HITS = 3,
  LOST_HIT_FRACTION = 4,
  MIN_PT = 5,
  CHARGE_SIGNIFICANCE = 6,
  LOOPER = 7,
  MAX_CCC_LOST_HITS = 8,
  NO_SEGMENTS_FOR_VALID_LAYERS = 9,
  SEED_EXTENSION = 10,
  SIZE = 12,
  NOT_STOPPED = 255
)

# Need to be kept consistent with
# DataFormats/TrackReco/interface/SeedStopReason.h
SeedStopReason = _Enum(
  UNINITIALIZED = 0,
  NOT_STOPPED = 1,
  SEED_CLEANING = 2,
  NO_TRAJECTORY = 3,
  SEED_REGION_REBUILD = 4,
  FINAL_CLEAN = 5,
  SMOOTHING_FAILED = 6,
  SIZE = 7
)

# to be kept is synch with enum HitSimType in TrackingNtuple.py
HitSimType = _Enum(
    Signal = 0,
    ITPileup = 1,
    OOTPileup = 2,
    Noise = 3,
    Unknown = 99
)

