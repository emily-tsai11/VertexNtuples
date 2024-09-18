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
