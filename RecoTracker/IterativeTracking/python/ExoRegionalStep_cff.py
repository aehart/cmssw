import FWCore.ParameterSet.Config as cms
import RecoTracker.IterativeTracking.iterativeTkConfig as _cfg
from Configuration.Eras.Modifier_fastSim_cff import fastSim

#for dnn classifier
from Configuration.ProcessModifiers.trackdnn_cff import trackdnn

#######################################################################
# Very large impact parameter tracking using TOB + TEC ring 5 seeding #
#######################################################################

exoRegionalStepClusters = _cfg.clusterRemoverForIter("ExoRegionalStep")
for _eraName, _postfix, _era in _cfg.nonDefaultEras():
    _era.toReplaceWith(exoRegionalStepClusters, _cfg.clusterRemoverForIter("ExoRegionalStep", _eraName, _postfix))

# TRIPLET SEEDING LAYERS
from RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi import *
exoRegionalStepSeedLayersTripl = cms.EDProducer("SeedingLayersEDProducer",
    layerList = cms.vstring(
    #TOB
    'TOB1+TOB2+MTOB3','TOB1+TOB2+MTOB4',
    #TOB+MTEC
    'TOB1+TOB2+MTEC1_pos','TOB1+TOB2+MTEC1_neg',
    ),
    TOB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         skipClusters   = cms.InputTag('exoRegionalStepClusters')
    ),
    MTOB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
         skipClusters   = cms.InputTag('exoRegionalStepClusters'),
         rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
    ),
    MTEC = cms.PSet(
        rphiRecHits    = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag('exoRegionalStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
        minRing = cms.int32(6),
        maxRing = cms.int32(7)
    )
)

# Triplet TrackingRegion
from RecoTracker.FinalTrackSelectors.exoRegionalStepInputTracks_cfi import exoRegionalStepInputTracks
from HiggsLongLived.TreeMaker.goodTrackProducer_cfi import goodTrackProducer as _goodTrackProducer
from HiggsLongLived.TreeMaker.generalV0Candidates_cfi import generalV0Candidates as _generalV0Candidates
from HiggsLongLived.TreeMaker.regionOfInterestProducer_cfi import regionOfInterestProducer as _regionOfInterestProducer
from HiggsLongLived.TreeMaker.regionOfInterestSelector_cfi import regionOfInterestSelector as _regionOfInterestSelector
from HiggsLongLived.TreeMaker.vertexOfInterestProducer_cfi import vertexOfInterestProducer as _vertexOfInterestProducer
from HiggsLongLived.TreeMaker.trueRegionOfInterestProducer_cfi import trueRegionOfInterestProducer as _trueRegionOfInterestProducer
exoRegionalStepTrackingRegionsStepZero = _goodTrackProducer.clone(tracks = cms.InputTag("exoRegionalStepInputTracks"))
exoRegionalStepTrackingRegionsStepOne = _generalV0Candidates.clone(
    trackRecoAlgorithm = cms.InputTag("exoRegionalStepTrackingRegionsStepZero", ""),
    doFit = cms.bool(False),
    useRefTracks = cms.bool(False),
    vtxDecayXYCut = cms.double(1.),
    ssVtxDecayXYCut = cms.double(5.),
    innerTkDCACut = cms.double(0.2),
    outerTkDCACut = cms.double(1.0),
    innerHitPosCut = cms.double(4.0),
    maxV0sCut = cms.int32(-1),
)
exoRegionalStepTrackingRegionsStepTwo = _regionOfInterestProducer.clone(
    minRadius = cms.double (2.0),
    trackClusters = cms.InputTag("exoRegionalStepTrackingRegionsStepOne", "Kshort"),
)
exoRegionalStepTrackingRegionsStepThree = _regionOfInterestSelector.clone(
    discriminatorCut = cms.double(0.5),
    genMatch = cms.string("None"),
    #genMatch = cms.string("HToSSTobbbb"),
    #genMatch = cms.string("stopToBottom"),
    graph_path = cms.FileInPath('HiggsLongLived/DeepSets/data/FullData_Phi-64-128-256_16-32-64_F-128-64-32_Model.pb'),
    input_names = cms.vstring(
        'phi_0',
        'phi_1'
    ),
    maxROIs = cms.int32(-1),
    minNConstituents = cms.uint32(0),
    minRadius = cms.double(2.0),
    nThreads = cms.uint32(1),
    output_names = cms.vstring('model_5/activation_10/Softmax'),
    regionsOfInterest = cms.InputTag("exoRegionalStepTrackingRegionsStepTwo", ""),
    singleThreadPool = cms.string('no_threads'),
)
exoRegionalStepTrackingRegionsStepFour = _vertexOfInterestProducer.clone(regionsOfInterest = cms.InputTag("exoRegionalStepTrackingRegionsStepThree", ""))
exoRegionalStepTrackingRegionsStepFive = _trueRegionOfInterestProducer.clone(physics = cms.string("HToSSTobbbb"))
#exoRegionalStepTrackingRegionsStepFive = _trueRegionOfInterestProducer.clone(physics = cms.string("stopToBottom"))

from RecoTracker.TkTrackingRegions.globalTrackingRegionWithVertices_cfi import globalTrackingRegionWithVertices as _globalTrackingRegionWithVertices
exoRegionalStepTrackingRegionsTripl = _globalTrackingRegionWithVertices.clone(RegionPSet = dict(
    originRadius = 1.0,
    fixedError = 1.0,
    VertexCollection = "exoRegionalStepTrackingRegionsStepFour",
    useFakeVertices = True,
    ptMin = 0.55,
    allowEmpty = True
))

from Configuration.Eras.Modifier_pp_on_XeXe_2017_cff import pp_on_XeXe_2017
from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
from RecoTracker.IterativeTracking.MixedTripletStep_cff import _mixedTripletStepTrackingRegionsCommon_pp_on_HI
(pp_on_XeXe_2017 | pp_on_AA_2018).toReplaceWith(exoRegionalStepTrackingRegionsTripl, 
                _mixedTripletStepTrackingRegionsCommon_pp_on_HI.clone(RegionPSet=dict(
                    ptMinScaling4BigEvts= False,
                    fixedError = 5.0,
                    ptMin = 2.0,
                    originRadius = 3.5
                )                                                                      )
)

# Triplet seeding
from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import ClusterShapeHitFilterESProducer as _ClusterShapeHitFilterESProducer
exoRegionalStepClusterShapeHitFilter = _ClusterShapeHitFilterESProducer.clone(
    ComponentName = 'exoRegionalStepClusterShapeHitFilter',
    doStripShapeCut = cms.bool(False),
    clusterChargeCut = dict(refToPSet_ = 'SiStripClusterChargeCutTight')
)

from RecoTracker.TkHitPairs.hitPairEDProducer_cfi import hitPairEDProducer as _hitPairEDProducer
exoRegionalStepHitDoubletsTripl = _hitPairEDProducer.clone(
    seedingLayers = "exoRegionalStepSeedLayersTripl",
    trackingRegions = "exoRegionalStepTrackingRegionsTripl",
    maxElement = 50000000,
    produceIntermediateHitDoublets = True,
)
from RecoTracker.TkSeedGenerator.multiHitFromChi2EDProducer_cfi import multiHitFromChi2EDProducer as _multiHitFromChi2EDProducer
exoRegionalStepHitTripletsTripl = _multiHitFromChi2EDProducer.clone(
    doublets = "exoRegionalStepHitDoubletsTripl",
    extraPhiKDBox = 0.01,
)
from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsEDProducer_cff import seedCreatorFromRegionConsecutiveHitsEDProducer as _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer
from RecoPixelVertexing.PixelLowPtUtilities.StripSubClusterShapeSeedFilter_cfi import StripSubClusterShapeSeedFilter as _StripSubClusterShapeSeedFilter
_exoRegionalStepSeedComparitorPSet = dict(
    ComponentName = 'CombinedSeedComparitor',
    mode = cms.string("and"),
    comparitors = cms.VPSet(
        cms.PSet(# FIXME: is this defined in any cfi that could be imported instead of copy-paste?
            ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
            FilterAtHelixStage = cms.bool(True),
            FilterPixelHits = cms.bool(False),
            FilterStripHits = cms.bool(True),
            ClusterShapeHitFilterName = cms.string('exoRegionalStepClusterShapeHitFilter'),
            ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache") # not really needed here since FilterPixelHits=False
        ),
        _StripSubClusterShapeSeedFilter.clone()
    )
)
exoRegionalStepSeedsTripl = _seedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer.clone(#empirically better than 'SeedFromConsecutiveHitsTripletOnlyCreator'
    seedingHitSets = "exoRegionalStepHitTripletsTripl",
    SeedComparitorPSet = _exoRegionalStepSeedComparitorPSet,
)
#fastsim
import FastSimulation.Tracking.TrajectorySeedProducer_cfi
_fastSim_exoRegionalStepSeedsTripl = FastSimulation.Tracking.TrajectorySeedProducer_cfi.trajectorySeedProducer.clone(
    trackingRegions = "exoRegionalStepTrackingRegionsTripl",
    hitMasks = cms.InputTag("exoRegionalStepMasks"),
)
from FastSimulation.Tracking.SeedingMigration import _hitSetProducerToFactoryPSet
_fastSim_exoRegionalStepSeedsTripl.seedFinderSelector.MultiHitGeneratorFactory = _hitSetProducerToFactoryPSet(exoRegionalStepHitTripletsTripl)
_fastSim_exoRegionalStepSeedsTripl.seedFinderSelector.MultiHitGeneratorFactory.SeedComparitorPSet=cms.PSet(  ComponentName = cms.string( "none" ) )
_fastSim_exoRegionalStepSeedsTripl.seedFinderSelector.MultiHitGeneratorFactory.refitHits = False
_fastSim_exoRegionalStepSeedsTripl.seedFinderSelector.layerList = exoRegionalStepSeedLayersTripl.layerList.value()
fastSim.toReplaceWith(exoRegionalStepSeedsTripl,_fastSim_exoRegionalStepSeedsTripl)

# PAIR SEEDING LAYERS
exoRegionalStepSeedLayersPair = cms.EDProducer("SeedingLayersEDProducer",
    layerList = cms.vstring('TOB1+TEC1_pos','TOB1+TEC1_neg', 
                            'TEC1_pos+TEC2_pos','TEC1_neg+TEC2_neg', 
                            'TEC2_pos+TEC3_pos','TEC2_neg+TEC3_neg', 
                            'TEC3_pos+TEC4_pos','TEC3_neg+TEC4_neg', 
                            'TEC4_pos+TEC5_pos','TEC4_neg+TEC5_neg', 
                            'TEC5_pos+TEC6_pos','TEC5_neg+TEC6_neg', 
                            'TEC6_pos+TEC7_pos','TEC6_neg+TEC7_neg'),
    TOB = cms.PSet(
         TTRHBuilder    = cms.string('WithTrackAngle'), clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         skipClusters   = cms.InputTag('exoRegionalStepClusters')
    ),
    TEC = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag('exoRegionalStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight')),
        minRing = cms.int32(5),
        maxRing = cms.int32(5)
    )
)
# Pair TrackingRegion
exoRegionalStepTrackingRegionsPair = _globalTrackingRegionWithVertices.clone(RegionPSet = dict(
    originRadius = 1.0,
    fixedError = 1.0,
    VertexCollection = "exoRegionalStepTrackingRegionsStepFour",
    useFakeVertices = True,
    ptMin = 0.6,
    allowEmpty = True
))

(pp_on_XeXe_2017 | pp_on_AA_2018).toReplaceWith(exoRegionalStepTrackingRegionsPair, 
                _mixedTripletStepTrackingRegionsCommon_pp_on_HI.clone(RegionPSet=dict(
                    ptMinScaling4BigEvts= False,
                    fixedError = 7.5,
                    ptMin = 2.0,
                    originRadius = 6.0
                )                                                                      )
)


# Pair seeds
exoRegionalStepHitDoubletsPair = _hitPairEDProducer.clone(
    seedingLayers = "exoRegionalStepSeedLayersPair",
    trackingRegions = "exoRegionalStepTrackingRegionsPair",
    produceSeedingHitSets = True,
    maxElementTotal = 12000000,
)
from RecoTracker.TkSeedGenerator.seedCreatorFromRegionConsecutiveHitsEDProducer_cff import seedCreatorFromRegionConsecutiveHitsEDProducer as _seedCreatorFromRegionConsecutiveHitsEDProducer
exoRegionalStepSeedsPair = _seedCreatorFromRegionConsecutiveHitsEDProducer.clone(
    seedingHitSets = "exoRegionalStepHitDoubletsPair",
    SeedComparitorPSet = _exoRegionalStepSeedComparitorPSet,
)
#fastsim
import FastSimulation.Tracking.TrajectorySeedProducer_cfi
fastSim.toReplaceWith(exoRegionalStepSeedsPair,
                      FastSimulation.Tracking.TrajectorySeedProducer_cfi.trajectorySeedProducer.clone(
        trackingRegions = "exoRegionalStepTrackingRegionsPair",
        hitMasks = cms.InputTag("exoRegionalStepMasks"),
        seedFinderSelector = dict(layerList = exoRegionalStepSeedLayersPair.layerList.value())
        )
)


# Combined seeds
import RecoTracker.TkSeedGenerator.GlobalCombinedSeeds_cfi
exoRegionalStepSeeds = RecoTracker.TkSeedGenerator.GlobalCombinedSeeds_cfi.globalCombinedSeeds.clone()
exoRegionalStepSeeds.seedCollections = cms.VInputTag(cms.InputTag('exoRegionalStepSeedsTripl'),cms.InputTag('exoRegionalStepSeedsPair'))

# LowPU
from Configuration.Eras.Modifier_trackingLowPU_cff import trackingLowPU
trackingLowPU.toModify(exoRegionalStepHitDoubletsPair, seedingLayers = 'exoRegionalStepSeedLayers')
trackingLowPU.toReplaceWith(exoRegionalStepSeeds, _seedCreatorFromRegionConsecutiveHitsEDProducer.clone(
    seedingHitSets = "exoRegionalStepHitDoubletsPair",
))


# QUALITY CUTS DURING TRACK BUILDING (for inwardss and outwards track building steps)
import TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff
_exoRegionalStepTrajectoryFilterBase = TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff.CkfBaseTrajectoryFilter_block.clone(
    maxLostHits = 0,
    minimumNumberOfHits = 5,
    minPt = 0.1,
    minHitsMinPt = 3
    )
exoRegionalStepTrajectoryFilter = _exoRegionalStepTrajectoryFilterBase.clone(
    seedPairPenalty = 1,
)
trackingLowPU.toReplaceWith(exoRegionalStepTrajectoryFilter, _exoRegionalStepTrajectoryFilterBase.clone(
    minimumNumberOfHits = 6,
))
for e in [pp_on_XeXe_2017, pp_on_AA_2018]:
    e.toModify(exoRegionalStepTrajectoryFilter, minPt=2.0)

exoRegionalStepInOutTrajectoryFilter = exoRegionalStepTrajectoryFilter.clone(
    minimumNumberOfHits = 4,
)


import RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi
exoRegionalStepChi2Est = RecoTracker.MeasurementDet.Chi2ChargeMeasurementEstimator_cfi.Chi2ChargeMeasurementEstimator.clone(
    ComponentName = cms.string('exoRegionalStepChi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(16.0),
    clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTight'))
)
trackingLowPU.toModify(exoRegionalStepChi2Est,
    clusterChargeCut = dict(refToPSet_ = 'SiStripClusterChargeCutTiny')
)

# TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi
exoRegionalStepTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.clone(
    MeasurementTrackerName = '',
    trajectoryFilter = cms.PSet(refToPSet_ = cms.string('exoRegionalStepTrajectoryFilter')),
    inOutTrajectoryFilter = cms.PSet(refToPSet_ = cms.string('exoRegionalStepInOutTrajectoryFilter')),
    useSameTrajFilter = False,
    minNrOfHitsForRebuild = 4,
    alwaysUseInvalidHits = False,
    maxCand = 2,
    estimator = cms.string('exoRegionalStepChi2Est'),
    #startSeedHitsInRebuild = True
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7)
    )
# Important note for LowPU: in RunI_ExoRegionalStep the
# inOutTrajectoryFilter parameter is spelled as
# inOutTrajectoryFilterName, and I suspect it has no effect there. I
# chose to "fix" the behaviour here, so the era is not fully
# equivalent to the customize. To restore the customize behaviour,
# uncomment the following lines
#trackingLowPU.toModify(exoRegionalStepTrajectoryBuilder,
#    inOutTrajectoryFilter = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.inOutTrajectoryFilter.clone(),
#    inOutTrajectoryFilterName = cms.PSet(refToPSet_ = cms.string('exoRegionalStepInOutTrajectoryFilter'))
#)

# MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
exoRegionalStepTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('exoRegionalStepSeeds'),
    clustersToSkip = cms.InputTag('exoRegionalStepClusters'),
    ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),

    TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('exoRegionalStepTrajectoryBuilder')),
    doSeedingRegionRebuilding = True,
    useHitsSplitting = True,
    cleanTrajectoryAfterInOut = True,
    TrajectoryCleaner = 'exoRegionalStepTrajectoryCleanerBySharedHits'
)
import FastSimulation.Tracking.TrackCandidateProducer_cfi
fastSim.toReplaceWith(exoRegionalStepTrackCandidates,
                      FastSimulation.Tracking.TrackCandidateProducer_cfi.trackCandidateProducer.clone(
        MinNumberOfCrossedLayers = 3,
        src = cms.InputTag("exoRegionalStepSeeds"),
        hitMasks = cms.InputTag("exoRegionalStepMasks")
        )
                      )


from TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi import trajectoryCleanerBySharedHits
exoRegionalStepTrajectoryCleanerBySharedHits = trajectoryCleanerBySharedHits.clone(
    ComponentName = cms.string('exoRegionalStepTrajectoryCleanerBySharedHits'),
    fractionShared = cms.double(0.09),
    allowSharedFirstHit = cms.bool(True)
    )
trackingLowPU.toModify(exoRegionalStepTrajectoryCleanerBySharedHits, fractionShared = 0.19)

# TRACK FITTING AND SMOOTHING OPTIONS
import TrackingTools.TrackFitters.RungeKuttaFitters_cff
exoRegionalStepFitterSmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.KFFittingSmootherWithOutliersRejectionAndRK.clone(
    ComponentName = 'exoRegionalStepFitterSmoother',
    EstimateCut = 30,
    MinNumberOfHits = 7,
    Fitter = cms.string('exoRegionalStepRKFitter'),
    Smoother = cms.string('exoRegionalStepRKSmoother')
    )
trackingLowPU.toModify(exoRegionalStepFitterSmoother, MinNumberOfHits = 8)

exoRegionalStepFitterSmootherForLoopers = exoRegionalStepFitterSmoother.clone(
    ComponentName = 'exoRegionalStepFitterSmootherForLoopers',
    Fitter = cms.string('exoRegionalStepRKFitterForLoopers'),
    Smoother = cms.string('exoRegionalStepRKSmootherForLoopers')
)

# Also necessary to specify minimum number of hits after final track fit
exoRegionalStepRKTrajectoryFitter = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectoryFitter.clone(
    ComponentName = cms.string('exoRegionalStepRKFitter'),
    minHits = 7
)
trackingLowPU.toModify(exoRegionalStepRKTrajectoryFitter, minHits = 8)

exoRegionalStepRKTrajectoryFitterForLoopers = exoRegionalStepRKTrajectoryFitter.clone(
    ComponentName = cms.string('exoRegionalStepRKFitterForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

exoRegionalStepRKTrajectorySmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectorySmoother.clone(
    ComponentName = cms.string('exoRegionalStepRKSmoother'),
    errorRescaling = 10.0,
    minHits = 7
)
trackingLowPU.toModify(exoRegionalStepRKTrajectorySmoother, minHits = 8)

exoRegionalStepRKTrajectorySmootherForLoopers = exoRegionalStepRKTrajectorySmoother.clone(
    ComponentName = cms.string('exoRegionalStepRKSmootherForLoopers'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
)

import TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi
exoRegionalFlexibleKFFittingSmoother = TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi.FlexibleKFFittingSmoother.clone(
    ComponentName = cms.string('exoRegionalFlexibleKFFittingSmoother'),
    standardFitter = cms.string('exoRegionalStepFitterSmoother'),
    looperFitter = cms.string('exoRegionalStepFitterSmootherForLoopers'),
)


# TRACK FITTING
import RecoTracker.TrackProducer.TrackProducer_cfi
exoRegionalStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'exoRegionalStepTrackCandidates',
    AlgorithmName = cms.string('exoRegionalStep'),
    #Fitter = 'exoRegionalStepFitterSmoother',
    Fitter = 'exoRegionalFlexibleKFFittingSmoother',
    )
fastSim.toModify(exoRegionalStepTracks, TTRHBuilder = 'WithoutRefit')


# TRACK SELECTION AND QUALITY FLAG SETTING.
from RecoTracker.FinalTrackSelectors.TrackMVAClassifierPrompt_cfi import *
from RecoTracker.FinalTrackSelectors.TrackMVAClassifierDetached_cfi import *
exoRegionalStepClassifier1 = TrackMVAClassifierDetached.clone()
exoRegionalStepClassifier1.src = 'exoRegionalStepTracks'
exoRegionalStepClassifier1.mva.GBRForestLabel = 'MVASelectorIter6_13TeV'
exoRegionalStepClassifier1.qualityCuts = [-0.6,-0.45,-0.3]
fastSim.toModify(exoRegionalStepClassifier1, vertices = "firstStepPrimaryVerticesBeforeMixing")

exoRegionalStepClassifier2 = TrackMVAClassifierPrompt.clone()
exoRegionalStepClassifier2.src = 'exoRegionalStepTracks'
exoRegionalStepClassifier2.mva.GBRForestLabel = 'MVASelectorIter0_13TeV'
exoRegionalStepClassifier2.qualityCuts = [0.0,0.0,0.0]
fastSim.toModify(exoRegionalStepClassifier2,vertices = "firstStepPrimaryVerticesBeforeMixing")

from RecoTracker.FinalTrackSelectors.ClassifierMerger_cfi import *
exoRegionalStep = ClassifierMerger.clone()
exoRegionalStep.inputClassifiers=['exoRegionalStepClassifier1','exoRegionalStepClassifier2']

from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
trackingPhase1.toReplaceWith(exoRegionalStep, exoRegionalStepClassifier1.clone(
     mva = dict(GBRForestLabel = 'MVASelectorTobTecStep_Phase1'),
     qualityCuts = [-0.6,-0.45,-0.3]
))

from RecoTracker.FinalTrackSelectors.TrackLwtnnClassifier_cfi import *
from RecoTracker.FinalTrackSelectors.trackSelectionLwtnn_cfi import *
trackdnn.toReplaceWith(exoRegionalStep, TrackLwtnnClassifier.clone(
     src = 'exoRegionalStepTracks',
     qualityCuts = [-0.4, -0.25, -0.1]
))
(trackdnn & fastSim).toModify(exoRegionalStep,vertices = "firstStepPrimaryVerticesBeforeMixing")

pp_on_AA_2018.toModify(exoRegionalStep, qualityCuts = [-0.6,-0.3,0.7])

import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
trackingLowPU.toReplaceWith(exoRegionalStep, RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.multiTrackSelector.clone(
    src = 'exoRegionalStepTracks',
    useAnyMVA = cms.bool(False),
    GBRForestLabel = cms.string('MVASelectorIter6'),
    trackSelectors = [
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
            name = 'exoRegionalStepLoose',
            chi2n_par = 0.4,
            res_par = ( 0.003, 0.001 ),
            minNumberLayers = 5,
            maxNumberLostLayers = 1,
            minNumber3DLayers = 2,
            d0_par1 = ( 2.0, 4.0 ),
            dz_par1 = ( 1.8, 4.0 ),
            d0_par2 = ( 2.0, 4.0 ),
            dz_par2 = ( 1.8, 4.0 )
        ),
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.tightMTS.clone(
            name = 'exoRegionalStepTight',
            preFilterName = 'exoRegionalStepLoose',
            chi2n_par = 0.3,
            res_par = ( 0.003, 0.001 ),
            minNumberLayers = 5,
            maxNumberLostLayers = 0,
            minNumber3DLayers = 2,
            d0_par1 = ( 1.5, 4.0 ),
            dz_par1 = ( 1.4, 4.0 ),
            d0_par2 = ( 1.5, 4.0 ),
            dz_par2 = ( 1.4, 4.0 )
        ),
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.highpurityMTS.clone(
            name = 'QualityMasks',
            preFilterName = 'exoRegionalStepTight',
            chi2n_par = 0.2,
            res_par = ( 0.003, 0.001 ),
            minNumberLayers = 5,
            maxNumberLostLayers = 0,
            minNumber3DLayers = 2,
            d0_par1 = ( 1.4, 4.0 ),
            dz_par1 = ( 1.3, 4.0 ),
            d0_par2 = ( 1.4, 4.0 ),
            dz_par2 = ( 1.3, 4.0 )
        ),
    ] #end of vpset
)) #end of clone



ExoRegionalStepTask = cms.Task(exoRegionalStepClusters,
                          exoRegionalStepSeedLayersTripl,
                          exoRegionalStepInputTracks,
                          exoRegionalStepTrackingRegionsStepZero,
                          exoRegionalStepTrackingRegionsStepOne,
                          exoRegionalStepTrackingRegionsStepTwo,
                          exoRegionalStepTrackingRegionsStepThree,
                          exoRegionalStepTrackingRegionsStepFour,
                          exoRegionalStepTrackingRegionsStepFive,
                          exoRegionalStepTrackingRegionsTripl,
                          exoRegionalStepHitDoubletsTripl,
                          exoRegionalStepHitTripletsTripl,
                          exoRegionalStepSeedsTripl,
                          exoRegionalStepSeedLayersPair,
                          exoRegionalStepTrackingRegionsPair,
                          exoRegionalStepHitDoubletsPair,
                          exoRegionalStepSeedsPair,
                          exoRegionalStepSeeds,
                          exoRegionalStepTrackCandidates,
                          exoRegionalStepTracks,
                          exoRegionalStepClassifier1,exoRegionalStepClassifier2,
                          exoRegionalStep)
ExoRegionalStep = cms.Sequence(ExoRegionalStepTask)


### Following are specific for LowPU, they're collected here to
### not to interfere too much with the default configuration
# SEEDING LAYERS
exoRegionalStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    layerList = cms.vstring('TOB1+TOB2', 
        'TOB1+TEC1_pos', 'TOB1+TEC1_neg', 
        'TEC1_pos+TEC2_pos', 'TEC2_pos+TEC3_pos', 
        'TEC3_pos+TEC4_pos', 'TEC4_pos+TEC5_pos', 
        'TEC5_pos+TEC6_pos', 'TEC6_pos+TEC7_pos', 
        'TEC1_neg+TEC2_neg', 'TEC2_neg+TEC3_neg', 
        'TEC3_neg+TEC4_neg', 'TEC4_neg+TEC5_neg', 
        'TEC5_neg+TEC6_neg', 'TEC6_neg+TEC7_neg'),
    TOB = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag('exoRegionalStepClusters'),
        TTRHBuilder = cms.string('WithTrackAngle'), clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTiny'))
    ),
    TEC = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag('exoRegionalStepClusters'),
        #    untracked bool useSimpleRphiHitsCleaner = false
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'), clusterChargeCut = cms.PSet(refToPSet_ = cms.string('SiStripClusterChargeCutTiny')),
        minRing = cms.int32(5),
        maxRing = cms.int32(5)
    )
)

trackingLowPU.toReplaceWith(ExoRegionalStepTask, 
    cms.Task(
    exoRegionalStepClusters,
    exoRegionalStepSeedLayers,
    exoRegionalStepTrackingRegionsPair,
    exoRegionalStepHitDoubletsPair,
    exoRegionalStepSeeds,
    exoRegionalStepTrackCandidates,
    exoRegionalStepTracks,
    exoRegionalStep
    )
)

#fastsim
import FastSimulation.Tracking.FastTrackerRecHitMaskProducer_cfi
exoRegionalStepMasks = FastSimulation.Tracking.FastTrackerRecHitMaskProducer_cfi.maskProducerFromClusterRemover(exoRegionalStepClusters)
fastSim.toReplaceWith(ExoRegionalStepTask,
                      cms.Task(exoRegionalStepMasks
                                   ,exoRegionalStepTrackingRegionsTripl
                                   ,exoRegionalStepSeedsTripl
                                   ,exoRegionalStepTrackingRegionsPair
                                   ,exoRegionalStepSeedsPair
                                   ,exoRegionalStepSeeds
                                   ,exoRegionalStepTrackCandidates
                                   ,exoRegionalStepTracks
                                   ,exoRegionalStepClassifier1,exoRegionalStepClassifier2
                                   ,exoRegionalStep
                                   )
)
