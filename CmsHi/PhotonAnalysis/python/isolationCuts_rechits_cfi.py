import FWCore.ParameterSet.Config as cms

import RecoEgamma.PhotonIdentification.isolationCalculator_cfi

isolationSumsCalculator100MeV  = RecoEgamma.PhotonIdentification.isolationCalculator_cfi.isolationSumsCalculator.clone(
    EcalRecHitThreshEA_Barrel       = cms.double(0.10),
    EcalRecHitThreshEB_Barrel       = cms.double(0.10)
    )

isolationSumsCalculator150MeV  = RecoEgamma.PhotonIdentification.isolationCalculator_cfi.isolationSumsCalculator.clone(
    EcalRecHitThreshEA_Barrel       = cms.double(0.15),
    EcalRecHitThreshEB_Barrel       = cms.double(0.15)
    )
isolationSumsCalculator200MeV  = RecoEgamma.PhotonIdentification.isolationCalculator_cfi.isolationSumsCalculator.clone(
    EcalRecHitThreshEA_Barrel       = cms.double(0.20),
    EcalRecHitThreshEB_Barrel       = cms.double(0.20)
    )







