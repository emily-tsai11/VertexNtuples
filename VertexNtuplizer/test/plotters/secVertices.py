sv_names = ["svInclusive", "svMerged", "svSlimmed"]

sv_vars = [
  "trk_tval",
  "trk_terr",
  "trk_tsig",
  # "trk_tqual",
  "trk_pt",
  "trk_pt2",
  "trk_eta",
  "trk_phi",
  "trk_dxy",
  "trk_dz",
  "trk_d3d",
  "trk_dxyerr",
  "trk_dzerr",
  "trk_d3derr",
  "trk_dxysig",
  "trk_dzsig",
  "trk_d3dsig",
  "trk_charge",
  "trk_chi2",
  "trk_ndof",
  "trk_chi2dof",
  "x",
  "y",
  "z",
  "xerr",
  "yerr",
  "zerr",
  "dxy",
  "dz",
  "d3d",
  "dxyerr",
  "dzerr",
  "d3derr",
  "dxysig",
  "dzsig",
  "d3dsig",
  "pt",
  "pt2",
  "eta",
  "phi",
  "tavg",
  "trange",
  "chi2",
  "ndof",
  "chi2dof",
  "ntrk"
]
sv_var_labels = [
  "Track time [ps]",
  "Track time error [ps]",
  "Track time significance",
  # "Track time quality",
  "Track p_{T} [GeV]",
  "Track p_{T}^{2} [GeV^{2}]",
  "Track #eta",
  "Track #phi",
  "Track d_{xy} [cm]",
  "Track d_{z} [cm]",
  "Track d_{3D} [cm]",
  "Track d_{xy} error [cm]",
  "Track d_{z} error [cm]",
  "Track d_{3D} error [cm]",
  "Track d_{xy} significance",
  "Track d_{z} significance",
  "Track d_{3D} significance",
  "Track charge",
  "Track #chi^{2}",
  "Track N.D.O.F.",
  "Track #chi^{2} / N.D.O.F.",
  "x [cm]",
  "y [cm]",
  "z [cm]",
  "x error [cm]",
  "y error [cm]",
  "z error [cm]",
  "d_{xy} [cm]",
  "d_{z} [cm]",
  "d_{3D} [cm]",
  "d_{xy} error [cm]",
  "d_{z} error [cm]",
  "d_{3D} error [cm]",
  "d_{xy} significance",
  "d_{z} significance",
  "d_{3D} significance",
  "p_{T} [GeV]",
  "p_{T}^{2} [GeV^{2}]",
  "#eta",
  "#phi",
  "<t> [ps]",
  "#Delta t [ps]",
  "#chi^{2}",
  "N.D.O.F.",
  "#chi^{2} / N.D.O.F.",
  "# of tracks"
]

sv_2d_names = [
  "svInclusiveMTDPV",
  "svMergedMTDPV",
]

sv_var2d = [
  "trk_pt_tval",
  "trk_pt_terr",
  "trk_pt_tsig",
  # "trk_pt_tqual",
  "trk_eta_tval",
  "trk_eta_terr",
  "trk_eta_tsig",
  # "trk_eta_tqual",
]
sv_var2d_labels = [
  ("p_{T} [GeV]", "Track time [ps]"),
  ("p_{T} [GeV]", "Track time error [ps]"),
  ("p_{T} [GeV]", "Track time significance"),
  # ("p_{T} [GeV]", "Track time quality"),
  ("#eta", "Track time [ps]"),
  ("#eta", "Track time error [ps]"),
  ("#eta", "Track time significance"),
  # ("#eta", "Track time quality"),
]
