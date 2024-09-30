#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include <Math/Functor.h>
#include <Fit/Fitter.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TF1.h>
#include <TH1F.h>
#include <TArc.h>

#define PDD pair<double, double>

using namespace std;

vector<PDD> propagate(double momentum, double magneticField, double detectorGap, int numDetectors);                                               // calculate the position when the particle hits the detectors
vector<PDD> gaussianSmearing(vector<PDD> &hits, double resolution);                                                                               // apply gaussian smearing to the hits to simulate the detector resolution
pair<PDD, double> calculateCircle(vector<PDD> &smearedHits);                                                                                      // calculate the circle that fits the smeared hits for the reconstruction
double fitCircle(vector<PDD> &smearedHits, double resolution);                                                                                    // fit a circle to the smeared hits
double reconstruction(vector<PDD> &smearedHits, double resolution, double magneticField, double detectorGap);                                     // reconstruct the momentum of the particle
double simulate(vector<PDD> &hits, double resolution, double magneticField, double detectorGap);                                                  // simulate the whole process
vector<PDD> positionMomentumResolution(vector<double> &momentums, double resolution, double magneticField, double detectorGap, int numDetectors); // simulate the whole process for different momentums

bool debug = false;    // Debug mode
bool detailed = true; // Detailed mode

// Constants
const double c = 299792458; // speed of light in m/s

// Particle properties (proton)
double charge = 1;    // Charge of the particle in e
double mass0 = 938e6; // Rest mass of the particle in eV/c^2

// Simulation properties
int N = debug ? 1 : 10000; // Number of simulations

int main()
{
    // Particle properties
    vector<double> momentums;
    for (int i = -5; i <= 25; i++)
        momentums.push_back(TMath::Power(10, i / 10.));
    if (debug)
        momentums.resize(1);

    // Magnetic field properties
    double magneticField = 1; // Magnetic field strength in Tesla

    // Detector properties
    int numDetectors = 5;                          // The number of detectors
    double detectorGap = 15e-2;                    // The gap between each pair of detectors
    vector<double> resolutions = {5, 10, 30, 100}; // The resolutions of the detectors in micrometers
    if (debug)
        resolutions.resize(1), resolutions[0] = 100;

    // Draw the relative resolution as a function of the momentum
    TCanvas *relativeResolutionCanvas = new TCanvas("relativeResolutionCanvas", "Relative Resolution", 1414, 1000);
    relativeResolutionCanvas->cd();
    relativeResolutionCanvas->SetGrid();

    TMultiGraph *relativeResolutionMultiGraph = new TMultiGraph();
    TLegend *relativeResolutionLegend = new TLegend(0.1, 0.7, 0.5, 0.9);

    int color[] = {kRed, kBlue, kGreen, kBlack};
    int colorIndex = 0;

    if (debug || detailed)
        relativeResolutionCanvas->Print("positionMomentum.pdf[");
    for (double resolution : resolutions)
    {
        vector<PDD> results = positionMomentumResolution(momentums, resolution, magneticField, detectorGap, numDetectors);

        TGraph *relativeResolutionGraph = new TGraph(results.size());
        for (int i = 0; i < results.size(); i++)
            relativeResolutionGraph->SetPoint(i, results[i].first, results[i].second);
        relativeResolutionGraph->SetLineColor(color[colorIndex %= 4]);
        relativeResolutionGraph->SetMarkerStyle(20);
        relativeResolutionGraph->SetMarkerSize(1);
        relativeResolutionGraph->SetMarkerColor(color[colorIndex++]);

        relativeResolutionMultiGraph->Add(relativeResolutionGraph);
        relativeResolutionLegend->AddEntry(relativeResolutionGraph, Form("Resolution = %g #mum", resolution));
    }
    if (debug || detailed)
        relativeResolutionCanvas->Print("positionMomentum.pdf]");

    gPad->SetLogx();
    gPad->SetLogy();
    relativeResolutionMultiGraph->SetTitle("Relative Resolution;Momentum (GeV/c);Relative Resolution");
    relativeResolutionMultiGraph->Draw("APL");
    relativeResolutionMultiGraph->GetXaxis()->SetNdivisions(505);

    relativeResolutionLegend->Draw();

    relativeResolutionCanvas->Print("relativeResolution.png");

    delete relativeResolutionCanvas;
    delete relativeResolutionMultiGraph;
    delete relativeResolutionLegend;

    return 0;
}

vector<PDD> propagate(double momentum, double magneticField, double detectorGap, int numDetectors)
{
    // Calculate the radius of the particle's motion
    double R = 1e9 * momentum / (abs(charge) * magneticField) / c; // Radius of the particle's motion in meters

    if (debug)
        cout << "Calculated radius: " << R << endl;

    // Calculate the position when the particle hits the detectors
    vector<PDD> hits;
    for (int i = 0; i < numDetectors; i++)
    {
        double xHit = detectorGap * i;
        double yHit = R - TMath::Sqrt(R * R - xHit * xHit);
        hits.push_back({xHit, yHit});
    }

    return hits;
}

vector<PDD> gaussianSmearing(vector<PDD> &hits, double resolution)
{
    TRandom3 random(0);
    vector<PDD> smearedHits;
    for (auto hit : hits)
    {
        double x = hit.first;
        double y = hit.second;
        double ySmeared = random.Gaus(y, resolution);
        smearedHits.push_back({x, ySmeared});
    }

    return smearedHits;
}

pair<PDD, double> calculateCircle(vector<PDD> &smearedHits)
{
    int n = smearedHits.size();
    double x1 = smearedHits[0].first;
    double y1 = smearedHits[0].second;
    double x2 = smearedHits[n / 2].first;
    double y2 = smearedHits[n / 2].second;
    double x3 = smearedHits[n - 1].first;
    double y3 = smearedHits[n - 1].second;

    double A = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
    double B = (x1 * x1 + y1 * y1) * (y3 - y2) + (x2 * x2 + y2 * y2) * (y1 - y3) + (x3 * x3 + y3 * y3) * (y2 - y1);
    double C = (x1 * x1 + y1 * y1) * (x2 - x3) + (x2 * x2 + y2 * y2) * (x3 - x1) + (x3 * x3 + y3 * y3) * (x1 - x2);
    double D = (x1 * x1 + y1 * y1) * (x3 * y2 - x2 * y3) + (x2 * x2 + y2 * y2) * (x1 * y3 - x3 * y1) + (x3 * x3 + y3 * y3) * (x2 * y1 - x1 * y2);

    double x = -B / (2 * A);
    double y = -C / (2 * A);
    double r = sqrt((B * B + C * C - 4 * A * D) / (4 * A * A));

    return {{x, y}, r};
}

double fitCircle(vector<PDD> &smearedHits, double resolution)
{
    auto chi2Function = [&](const double *par)
    {
        double chi2 = 0;
        for (int i = 0; i < smearedHits.size(); i++)
        {
            double dx = smearedHits[i].first - par[0];
            double dy = smearedHits[i].second - par[1];
            double r = TMath::Sqrt(dx * dx + dy * dy);
            double k = 1 / r;
            double dk = par[2] - k;
            double err = (dy * dy) / TMath::Power(r, 3) * resolution;
            chi2 += dk * dk / (err * err);
        }

        return chi2;
    };

    ROOT::Math::Functor chi2Functor(chi2Function, 3); // warp the chi2 function into a functor
    ROOT::Fit::Fitter fitter;
    // fitter.Config().MinimizerOptions().SetTolerance(1e2);
    if (debug) // set the debug level
        fitter.Config().MinimizerOptions().SetPrintLevel(3);

    pair<PDD, double> circle = calculateCircle(smearedHits);
    double pStart[3] = {circle.first.first, circle.first.second, 1 / circle.second};
    fitter.SetFCN(chi2Functor, pStart);

    bool ok = fitter.FitFCN();
    if (!ok)
    {
        if (debug)
            cout << "Fit failed" << endl;
        return 1e100;
    }

    const double *pFit = fitter.Result().GetParams();
    double radius = 1 / pFit[2] * pFit[1] / abs(pFit[1]);

    if (debug)
    {
        cout << "Fitted radius: " << radius << endl;

        // Plot the smeared hits and the fitted circle
        TCanvas *fitCircleCanvas = new TCanvas("fitCircleCanvas", "Fit Circle", 600, 600);
        fitCircleCanvas->cd();
        fitCircleCanvas->SetGrid();
        fitCircleCanvas->DrawFrame(-1.1, -0.1, 1.1, 2.1);
        TGraph *fitCircleGraph = new TGraph(smearedHits.size());
        for (int i = 0; i < smearedHits.size(); i++)
            fitCircleGraph->SetPoint(i, smearedHits[i].first, smearedHits[i].second);
        fitCircleGraph->SetMarkerStyle(20);
        fitCircleGraph->SetMarkerSize(1);
        fitCircleGraph->Draw("P");

        auto result = fitter.Result();
        TArc *fitCircleArc = new TArc(result.Parameter(0), result.Parameter(1), 1 / result.Parameter(2));
        fitCircleArc->SetLineColor(kRed);
        fitCircleArc->SetLineWidth(4);
        fitCircleArc->SetFillStyle(0);
        fitCircleArc->Draw();

        fitCircleCanvas->Print("fitCircle.png");

        // Print the fit result
        result.Print(std::cout);

        delete fitCircleCanvas;
        delete fitCircleGraph;
        delete fitCircleArc;
    }

    return radius;
}

double reconstruction(vector<PDD> &smearedHits, double resolution, double magneticField, double detectorGap)
{
    // Calculate the radius of the particle's motion
    double R = fitCircle(smearedHits, resolution);

    if (R == 1e100)
        return R;

    // Calculate the momentum of the particle
    double momentum = R * abs(charge) * magneticField * c / 1e9;

    return momentum;
}

double simulate(vector<PDD> &hits, double resolution, double magneticField, double detectorGap)
{
    resolution *= 1e-6; // Convert the resolution to meters

    vector<PDD> smearedHits = gaussianSmearing(hits, resolution);

    if (debug)
    {
        cout << "Smeared hits:" << endl;
        for (auto hit : smearedHits)
            cout << hit.first << " " << hit.second << endl;
    }

    return reconstruction(smearedHits, resolution, magneticField, detectorGap);
}

vector<PDD> positionMomentumResolution(vector<double> &momentums, double resolution, double magneticField, double detectorGap, int numDetectors)
{
    vector<PDD> results;

    TCanvas *positionMomentumCanvas = new TCanvas("positionMomentumCanvas", "Position-Momentum", 1414, 1000);

    for (double momentum : momentums)
    {
        vector<PDD> hits = propagate(momentum, magneticField, detectorGap, numDetectors);
        if (debug)
        {
            cout << "Hits:" << endl;
            for (auto hit : hits)
                cout << hit.first << " " << hit.second << endl;
        }

        TH1F *positionMomentumHistogram = new TH1F("positionMomentum", Form("Momentum Distribution  (Resolution = %g #mum);1/Momentum (c/GeV);Counts", resolution),
                                                   100, 0, 0);

        for (int i = 0; i < N; ++i)
        {
            double reconstructedMomentum = simulate(hits, resolution, magneticField, detectorGap);
            if (reconstructedMomentum == 1e100)
                continue;
            positionMomentumHistogram->Fill(1 / reconstructedMomentum);
        }
        if (debug)
            positionMomentumHistogram->Fit("gaus");
        else
            positionMomentumHistogram->Fit("gaus", "Q");

        positionMomentumHistogram->GetXaxis()->SetNdivisions(505);

        positionMomentumCanvas->cd();
        positionMomentumHistogram->Draw();
        if (debug || detailed)
            positionMomentumCanvas->Print("positionMomentum.pdf");

        double mean = positionMomentumHistogram->GetFunction("gaus")->GetParameter(1);
        double sigma = positionMomentumHistogram->GetFunction("gaus")->GetParameter(2);

        results.push_back({momentum, sigma / mean});

        delete positionMomentumHistogram;
    }

    delete positionMomentumCanvas;

    return results;
}
