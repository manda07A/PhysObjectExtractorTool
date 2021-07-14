# Implementation of the plotting step of the analysis
#
# The plotting combines the histograms to plots which allow us to study the
# inital dataset based on observables motivated through physics.

import ROOT
ROOT.gROOT.SetBatch(True)


# Declare a human-readable label for each variable
labels = {
        "npv": "Number of primary vertices",
        }

# Specify the color for each process
colors = {
        "ggH": ROOT.TColor.GetColor("#BF2229"),
        "qqH": ROOT.TColor.GetColor("#00A88F"),
        "TT": ROOT.TColor.GetColor(155, 152, 204),
        "W": ROOT.TColor.GetColor(222, 90, 106),
        "QCD":  ROOT.TColor.GetColor(250, 202, 255),
        "ZLL": ROOT.TColor.GetColor(100, 192, 232),
        "ZTT": ROOT.TColor.GetColor(248, 206, 104),
        }

# Retrieve a histogram from the input file based on the process and the variable
# name
def getHistogram(tfile, name, variable, tag=""):
    print (name)
    print (variable)
    name = "%s_%s%s" % (name, variable, tag)
    h = tfile.Get(name)
    if not h:
        raise Exception("Failed to load histogram %s." % (name))
    return h

# Main function of the plotting step
#
# The major part of the code below is dedicated to define a nice-looking layout.
# The full version of this analysis continas the estimation of QCD background
def main(variable):
    tfile = ROOT.TFile("histograms.root", "READ")

    # Styles
    ROOT.gStyle.SetOptStat(0)

    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetCanvasDefW(600)
    ROOT.gStyle.SetCanvasDefX(0)
    ROOT.gStyle.SetCanvasDefY(0)

    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    ROOT.gStyle.SetEndErrorSize(2)
    ROOT.gStyle.SetMarkerStyle(20)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleXOffset(1.00)
    ROOT.gStyle.SetTitleYOffset(1.60)

    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")

    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetPaperSize(20., 20.)
    ROOT.gStyle.SetHatchesLineWidth(5)
    ROOT.gStyle.SetHatchesSpacing(0.05)

    #ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "Y")


    # Simulation
    ZLL = getHistogram(tfile, "ZLL", variable)

    # Data
    data = getHistogram(tfile, "dataRunB", variable)
    dataRunC = getHistogram(tfile, "dataRunC", variable)
    data.Add(dataRunC)

    # Draw histograms
    data.SetMarkerStyle(20)
    data.SetLineColor(ROOT.kBlack)
    
    for x, l in [(ZLL, "ZLL")]:
        x.SetLineWidth(0)
        x.SetFillColor(colors[l])
        
    stack = ROOT.THStack("", "")
    for x in [ZLL]:
        stack.Add(x)

    c = ROOT.TCanvas("", "", 600, 600)
    stack.Draw("hist")
    name = data.GetTitle()
    if name in labels:
        title = labels[name]
    else:
        title = name
    stack.GetXaxis().SetTitle(title)
    stack.GetYaxis().SetTitle("N_{Events}")
    stack.SetMaximum(max(stack.GetMaximum(), data.GetMaximum()) * 1.4)
    stack.SetMinimum(1.0)

    data.Draw("E1P SAME")
    
    # Add legend
    legend = ROOT.TLegend(0.4, 0.73, 0.90, 0.88)
    legend.SetNColumns(2)
    legend.AddEntry(ZLL, "Z#rightarrowll", "f")
    legend.AddEntry(data, "Data", "lep")
    legend.SetBorderSize(0)
    legend.Draw()

    # Add title
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatex(0.6, 0.935, "11.5 fb^{-1} (2012, 8 TeV)")
    latex.DrawLatex(0.16, 0.935, "#bf{CMS Open Data}")

    # Save
    c.SaveAs("%s.pdf" % (variable))
    c.SaveAs("%s.png" % (variable))

# Loop over all variable names and make a plot for each
if __name__ == "__main__":
    for variable in labels.keys():
        main(variable)
