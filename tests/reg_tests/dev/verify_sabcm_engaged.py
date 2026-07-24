import os, sys, copy
sys.path.insert(0, '/home/mdo/MDOLab_3_v2/adflow_sabcm/tests/reg_tests')
from adflow import ADFLOW
from reg_default_options import adflowDefOpts
from reg_bcm import ap_bcm_tut_wing, bcmBaseOptionsSmooth, bcmPlainSAOptions, bcmAeroDVs, getStateBlocks, maskStateVector
import numpy

restartfile = "/home/mdo/MDOLab_3_v2/adflow_sabcm/tests/reg_tests/output_files/mdo_tutorial_bcm_smooth_000_vol.cgns"

def build(options):
    o = copy.copy(adflowDefOpts)
    o["outputdirectory"] = "/home/mdo/MDOLab_3_v2/adflow_sabcm/tests/reg_tests/output_files"
    o.update(copy.deepcopy(options))
    o["restartfile"] = restartfile
    o["ncycles"] = 1
    o["usenksolver"] = False
    o["useanksolver"] = False
    s = ADFLOW(options=o, debug=True)
    ap = copy.deepcopy(ap_bcm_tut_wing)
    for dv in bcmAeroDVs:
        ap.addDV(dv)
    res = s.getResidual(ap)
    return s, ap, res

sSABCM, apSABCM, resSABCM = build(bcmBaseOptionsSmooth)
sSA, apSA, resSA = build(bcmPlainSAOptions)

nw, blocks = getStateBlocks(sSABCM)
diffNorm = numpy.sqrt(numpy.sum((resSABCM - resSA) ** 2))
resSABCMNorm = numpy.sqrt(numpy.sum(resSABCM ** 2))
resSANorm = numpy.sqrt(numpy.sum(resSA ** 2))
print("residual(use_SABCM=True) norm:", resSABCMNorm)
print("residual(use_SABCM=False, plain SA) norm:", resSANorm)
print("||residual(SABCM) - residual(plain SA)||:", diffNorm)

# nuTilde-row-only comparison -- SABCM's production multiplier acts through this row
rSABCM_nu = maskStateVector(resSABCM, nw, blocks["nuTilde"])
rSA_nu = maskStateVector(resSA, nw, blocks["nuTilde"])
print("||resid[nuTilde](SABCM)||:", numpy.sqrt(numpy.sum(rSABCM_nu**2)))
print("||resid[nuTilde](plain SA)||:", numpy.sqrt(numpy.sum(rSA_nu**2)))
print("||resid[nuTilde](SABCM) - resid[nuTilde](plain SA)||:", numpy.sqrt(numpy.sum((rSABCM_nu-rSA_nu)**2)))

wDot = sSABCM.getStatePerturbation(314)
dRdw_SABCM = sSABCM.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
dRdw_SA = sSA.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
print("||dR/dw*wDot (SABCM)||:", numpy.sqrt(numpy.sum(dRdw_SABCM**2)))
print("||dR/dw*wDot (plain SA)||:", numpy.sqrt(numpy.sum(dRdw_SA**2)))
print("||dR/dw*wDot (SABCM) - (plain SA)||:", numpy.sqrt(numpy.sum((dRdw_SABCM-dRdw_SA)**2)))
