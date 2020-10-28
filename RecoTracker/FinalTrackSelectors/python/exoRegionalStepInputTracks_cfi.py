import FWCore.ParameterSet.Config as cms
from RecoTracker.FinalTrackSelectors.TrackCollectionMerger_cfi import *
from RecoTracker.FinalTrackSelectors.trackAlgoPriorityOrder_cfi import trackAlgoPriorityOrder

import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
exoRegionalStepInputTracks =  TrackCollectionMerger.clone()
exoRegionalStepInputTracks.trackProducers = ['earlyGeneralTracks',
                                             'muonSeededTracksInOut',
                                             'muonSeededTracksOutIn'
                                     ]
exoRegionalStepInputTracks.inputClassifiers =["earlyGeneralTracks",
                                              "muonSeededTracksInOutClassifier",
                                              "muonSeededTracksOutInClassifier"
                                      ]
from Configuration.Eras.Modifier_trackingLowPU_cff import trackingLowPU
trackingLowPU.toModify(exoRegionalStepInputTracks,
    trackProducers = [
        'earlyGeneralTracks',
        'muonSeededTracksInOut',
        'muonSeededTracksOutIn'
    ],
    inputClassifiers = [
        "earlyGeneralTracks",
        "muonSeededTracksInOutClassifier",
        "muonSeededTracksOutInClassifier"
    ]
)
from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
_forPhase1 = dict(
    trackProducers = [
        'earlyGeneralTracks',
        'muonSeededTracksInOut',
        'muonSeededTracksOutIn'
    ],
    inputClassifiers = [
        "earlyGeneralTracks",
        "muonSeededTracksInOutClassifier",
        "muonSeededTracksOutInClassifier"
    ],
)
trackingPhase1.toModify(exoRegionalStepInputTracks, **_forPhase1)

# For Phase2PU140
from Configuration.Eras.Modifier_trackingPhase2PU140_cff import trackingPhase2PU140
from RecoTracker.FinalTrackSelectors.trackListMerger_cfi import trackListMerger as _trackListMerger
trackingPhase2PU140.toReplaceWith(exoRegionalStepInputTracks, _trackListMerger.clone(
    TrackProducers =['initialStepTracks',
                     'highPtTripletStepTracks',
                     'lowPtQuadStepTracks',
                     'lowPtTripletStepTracks',
                     'detachedQuadStepTracks',
                     'pixelPairStepTracks',
                    ],
    hasSelector = [1,1,1,1,1,1],
    indivShareFrac = [1.0,0.16,0.095,0.09,0.09,0.09],
    selectedTrackQuals = cms.VInputTag(cms.InputTag("initialStepSelector","initialStep"),
                                       cms.InputTag("highPtTripletStepSelector","highPtTripletStep"),
                                       cms.InputTag("lowPtQuadStepSelector","lowPtQuadStep"),
                                       cms.InputTag("lowPtTripletStepSelector","lowPtTripletStep"),
                                       cms.InputTag("detachedQuadStep"),
                                       cms.InputTag("pixelPairStepSelector","pixelPairStep"),
                                       ),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4,5), pQual=cms.bool(True) )
                             ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
    )
)
