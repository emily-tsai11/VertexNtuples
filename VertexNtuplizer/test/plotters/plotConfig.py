import sys
sys.path.insert(0, "../cmsstyle/src/cmsstyle")
import cmsstyle as CMS

from genVertices import *
from secVertices import *
from vtxMatch import *
from efficiencies import *

# fname_TTnoPU = "/eos/user/e/etsai/workspace/VertexNtuples_CMSSW_14_1_0_pre3/src/VertexNtuples/VertexNtuplizer/test/histo_TTToHadronicNoPU.root"
# fname_TTPU200 = "/eos/user/e/etsai/workspace/VertexNtuples_CMSSW_14_1_0_pre3/src/VertexNtuples/VertexNtuplizer/test/histo_TTToHadronicPU200.root"
fname_TTnoPU = "/eos/cms/store/group/phys_btag/etsai/VertexNtuples/TTToHadronic_noPU.root"
fname_TTPU200 = "/eos/cms/store/group/phys_btag/etsai/VertexNtuples/TTToHadronic_PU200.root"
# savePath = "/eos/home-b/btvweb/www/Offline/Phase2/etsai/Phase2/"
save_path = "/eos/home-b/btvweb/www/Offline/Phase2/etsai/Run2/"

COLORS = [CMS.p10.kBlue, CMS.p10.kYellow, CMS.p10.kRed, CMS.p10.kGray, CMS.p10.kViolet, CMS.p10.kBrown, CMS.p10.kOrange, CMS.p10.kGreen, CMS.p10.kAsh, CMS.p10.kCyan]

ttnoPU = "t\\bar{t} no PU"
ttPU200 = "t\\bar{t} PU 200"

# match_vars = [
#   "xres",
#   "yres",
#   "zres",
#   "xpull",
#   "ypull",
#   "zpull",
#   "trk_deltaR",
#   "trk_ptresnorm",
#   "matchdxy",
#   "matchd3d",
#   "matchdxyerr",
#   "matchd3derr",
#   "matchdxysig",
#   "matchd3dsig"
# ]
# match_labels = [
#   "x resolution [cm]",
#   "y resolution [cm]",
#   "z resolution [cm]",
#   "x pull",
#   "y pull",
#   "z pull",
#   "SV track to GV daughter #Delta R",
#   "|GV daughter p_{T} - SV track p_{T}| / (GV daughter p_{T} + SV track p_{T})",
#   "SV to GV d_{xy} [cm]",
#   "SV to GV d_{3D} [cm]",
#   "SV to GV d_{xy} error [cm]",
#   "SV to GV d_{3D} error [cm]",
#   "SV to GV d_{xy} significance",
#   "SV to GV d_{3D} significance"
# ]

# gv_sv_vars = gv_vars.copy()
# gv_sv_vars.extend(match_vars)
# gv_sv_labels = gv_labels.copy()
# gv_sv_labels.extend(match_labels)
# sv_gv_vars = sv_vars.copy()
# sv_gv_vars.extend(match_vars)
# sv_gv_labels = sv_labels.copy()
# sv_gv_labels.extend(match_labels)
