import ROOT
import numpy as np

# Make global style changes
ROOT.gStyle.SetOptStat(0) # Disable the statistics box
ROOT.gStyle.SetTextFont(42)

# Create a canvas
c = ROOT.TCanvas('c', 'my canvas', 800, 600)

# Create a histogram with some dummy data and draw it
TFile myFile("MuonInfo.root")
mtree->Draw("muon_e")

# Save as png file and show interactively
c.SaveAs('dummy_data.png')
c.Draw()
