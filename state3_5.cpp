#include <iostream>
#include <vector>
#include <unordered_map>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
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

using namespace std;
using PDD = pair<double, double>;
using PVV = pair<TVector2, TVector2>;

TRandom3 *Random = new TRandom3(0); // Random number generator

vector<TVector3> propagate(double momentum);                                          // calculate the position when the particle hits the detectors
vector<TVector3> gaussianSmearing(vector<TVector3> &hits, double resolution);         // apply gaussian smearing to the hits to simulate the detector resolution
pair<PDD, double> calculateCircle(vector<TVector3> &smearedHits);                     // calculate the circle that fits the smeared hits for the reconstruction
double fitCircle(vector<TVector3> &smearedHits, double resolution);                   // fit a circle to the smeared hits
double fitLine(vector<TVector3> &smearedHits);                                        // fit a line to the smeared hits
double reconstruction(vector<TVector3> &smearedHits, double resolution);              // reconstruct the momentum of the particle
double simulate(double momentum, double resolution);                                  // simulate the whole process
vector<PDD> positionMomentumResolution(vector<double> &momentums, double resolution); // simulate the whole process for different momentums

bool debug = false;   // Debug mode
bool detailed = true; // Detailed mode

// Constants
const double c = 299792458; // speed of light in m/s

// Particle properties (proton)
double charge = 1;     // Charge of the particle in e
double mass0 = 938e-3; // Rest mass of the particle in GeV/c^2

// Magnetic field properties
double magneticField = 1; // Magnetic field strength in Tesla

// Simulation properties
int N = debug ? 1 : 10000; // Number of simulations

// Detector properties
int numDetectors = 5;       // The number of detectors
double detectorGap = 15e-2; // The gap between each pair of detectors
enum materials
{
    silicon,
    carbon,
    vacuum
};

int main()
{
    // Particle properties
    vector<double> momentums;
    for (int i = -3; i <= 25; i++)
        momentums.push_back(TMath::Power(10, i / 10.));
    if (debug)
        momentums.resize(1);

    // Detector properties
    vector<double> resolutions = {5, 10, 30, 100}; // The resolutions of the detectors in micrometers
    if (debug)
        resolutions.resize(1);

    // Draw the relative resolution as a function of the momentum
    TCanvas *relativeResolutionCanvas = new TCanvas("relativeResolutionCanvas", "Relative Resolution", 1410, 1000);
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
        vector<PDD> results = positionMomentumResolution(momentums, resolution);

        TGraph *relativeResolutionGraph = new TGraph(results.size());
        for (int i = 0; i < results.size(); i++)
            relativeResolutionGraph->SetPoint(i, results[i].first, results[i].second);
        relativeResolutionGraph->SetLineColor(color[colorIndex %= 4]);
        relativeResolutionGraph->SetMarkerStyle(20);
        relativeResolutionGraph->SetMarkerSize(1);
        relativeResolutionGraph->SetMarkerColor(color[colorIndex++]);

        relativeResolutionMultiGraph->Add(relativeResolutionGraph);
        relativeResolutionLegend->AddEntry(relativeResolutionGraph, Form("Resolution = %g #mum", resolution));

        static int count = 0;
        cout << "Resolution caculation progress: " << ++count * 100 / resolutions.size() << "%" << endl;
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

PVV propagateInMagneticField(TVector3 position, TVector3 &momentum, double distance)
{
    TVector2 momentum2D(momentum.Y(), momentum.Z()); // momentum of the particle in the y-z plane
    TVector2 position2D(position.Y(), position.Z()); // position of the particle in the y-z plane

    double mass = TMath::Sqrt(momentum.Mag2() + TMath::Power(mass0, 2));          // Mass of the particle in GeV/c^2
    double r = momentum2D.Mod() * 1e9 / (charge * magneticField * c);             // Radius of the particle's motion in meters
    double T = TMath::TwoPi() * r / (momentum2D.Mod() * c / mass);                // Period of the particle's motion in seconds
    TVector2 center = position2D + r * momentum2D.Rotate(TMath::Pi() / 2).Unit(); // Center of the particle's motion in the y-z plane

    double newZ = position.Z() + distance;
    double newY = center.X() + TMath::Sqrt(r * r - TMath::Power(newZ - center.Y(), 2));
    if (isnan(newY))
    {
        cerr << "[Error] The particle cannot be fully detected!" << endl;
        cerr << "[Error] The particle cannot be fully detected!" << endl;
        cerr << "[Error] The particle cannot be fully detected!" << endl;
        exit(1);
    }
    double dY = newY - position.Y();

    double d2 = dY * dY + distance * distance;
    double angle = TMath::ACos((r * r * 2 - d2) / (2 * r * r));

    double t = T * angle / TMath::TwoPi();
    double dX = t * momentum.X() * c / mass;

    return {TVector2(dX, dY), TVector2(angle, 0)};
}

double getStoppingPower(materials material, TVector3 &momentum, double distance)
{

    double Z, A, phi, I;
    switch (material)
    {
    case silicon:
        Z = 14, A = 28.0855, phi = 2.329, I = 563.446635e-9;
        break;
    case carbon:
        Z = 6, A = 12.0107, phi = 2.210, I = 171.684645e-9;
        break;
    case vacuum:
    default:
        return 1;
    }

    double mass = TMath::Sqrt(momentum.Mag2() + TMath::Power(mass0, 2)); // Mass of the particle in GeV/c^2
    double beta = momentum.Mag() / mass;                                 // Velocity of the particle in c
    double gamma = 1. / TMath::Sqrt(1 - beta * beta);                    // Lorentz factor of the particle
    double mass0e = 0.511e-6;                                            // Rest mass of the electron in GeV/c^2

    double K = 0.307075e-3; // The constant in GeV mol^-1 cm^2
    double z = charge;      // The atomic number of the incident particle

    double Wmax = 2 * mass0e * beta * beta * gamma * gamma / (1 + 2 * gamma * mass0e / mass0 + TMath::Power(mass0e / mass0, 2)); // Maximum energy transfer in GeV
    double hbarOmegap = TMath::Sqrt(phi * (Z / A)) * 28.816e-9;                                                                  // Plasmon energy in GeV
    double halfDeltaBetaGamma = TMath::Log(hbarOmegap / I) + TMath::Log(beta * gamma) - 0.5;                                     // The density effect correction

    double S = K * z * z * Z / (A * beta * beta) * (0.5 * TMath::Log(2 * mass0e * beta * beta * gamma * gamma * Wmax / (I * I)) - beta * beta - halfDeltaBetaGamma) * phi * distance * 100; // Stopping power in GeV

    double x = distance * 1e2 * phi; // Thickness of the material in g/cm^-2
    double sigma = 0.5 * K * z * z * x * Z / (A * beta * beta);
    double deltaP = sigma * (TMath::Log(2 * mass * beta * beta * gamma * gamma / I) + TMath::Log(sigma / I) + 0.2 - beta * beta - halfDeltaBetaGamma); // Fluctuation in energy loss in GeV

    // cout << S << endl; exit(0);

    double energyLoss = Random->Gaus(S, deltaP);

    return 1 - energyLoss / momentum.Mag();
}

PVV propagateByMultipleScattering(TVector3 &momentum, double distance, materials material)
{
    double X0;
    switch (material)
    {
    case silicon:
        X0 = 21.82 / 2.329 * 1e-2; // Radiation length of silicon in meters (X0 [g/cm^-2] / density [g/cm^-3])
        break;
    case carbon:
        X0 = 42.70 / 2.210 * 1e-2; // Radiation length of carbon in meters (X0 [g/cm^-2] / density [g/cm^-3])
        break;
    case vacuum:
    default:
        return {TVector2(0, 0), TVector2(0, 0)};
    }

    double mass = TMath::Sqrt(momentum.Mag2() + TMath::Power(mass0, 2)); // Mass of the particle in GeV/c^2
    TVector3 velocity = momentum * (1. / mass);                          // Velocity of the particle in c

    double thickness = distance / X0; // Thickness of the material in radiation lengths
    double theta0 = 13.6e-3 / (velocity.Mag() * momentum.Mag()) * charge * TMath::Sqrt(thickness) *
                    (1 + 0.038 * TMath::Log(thickness * charge * charge / velocity.Mag2())); // Multiple scattering angle in rad

    double z1x = Random->Gaus(0, 1), z1y = Random->Gaus(0, 1);
    double z2x = Random->Gaus(0, 1), z2y = Random->Gaus(0, 1);

    auto getD = [&](double z1, double z2)
    {
        return z1 * distance * theta0 / TMath::Sqrt(12) + z2 * distance * theta0 / 2.;
    };

    return {TVector2(getD(z1x, z2x), getD(z1y, z2y)), TVector2(z2x * theta0, z2y * theta0)};
}

vector<TVector3> propagate(double momentumz)
{
    // Simulate the particle's motion
    vector<TVector3> hits;

    TVector3 position(0, 0, 0);         // Initial position of the particle in meters
    TVector3 momentum(0, 0, momentumz); // Initial momentum of the particle in GeV/c

    for (int i = 0; i < numDetectors; i++)
    {
        double distance;

        auto prop_ = [&](double distance, materials material)
        {
            auto [d1, a1] = propagateInMagneticField(position, momentum, distance);
            auto [d2, a2] = propagateByMultipleScattering(momentum, distance, material);
            momentum *= getStoppingPower(material, momentum, distance);
            auto d = d1 + d2, a = a1 + a2;
            position += TVector3(d.X(), d.Y(), distance);
            momentum.RotateX(a.X());
            momentum.RotateY(a.Y());
        };

        // propagate in detectors
        prop_(3e-4, silicon); // The first silicon detector (detector the particle's x position)
        prop_(1e-2, carbon);  // The carbon support
        prop_(3e-4, silicon); // The second silicon detector (detector the particle's y position)
        hits.push_back(position);

        // propagate in detector gaps
        if (i < numDetectors - 1)
            prop_(detectorGap - 106e-4, vacuum); // The vacuum gap (0 means X0 is infinity)
    }

    return hits;
}

vector<TVector3> gaussianSmearing(vector<TVector3> &hits, double resolution)
{
    vector<TVector3> smearedHits;
    for (auto hit : hits)
    {
        double x = hit.X();
        double y = hit.Y();
        double xSmeared = Random->Gaus(x, resolution);
        double ySmeared = Random->Gaus(y, resolution);
        smearedHits.push_back({xSmeared, ySmeared, hit.Z()});
    }

    return smearedHits;
}

pair<PDD, double> calculateCircle(vector<TVector3> &smearedHits)
{
    int n = smearedHits.size();
    double x1 = smearedHits[0].Z();
    double y1 = smearedHits[0].Y();
    double x2 = smearedHits[n / 2].Z();
    double y2 = smearedHits[n / 2].Y();
    double x3 = smearedHits[n - 1].Z();
    double y3 = smearedHits[n - 1].Y();

    double A = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
    double B = (x1 * x1 + y1 * y1) * (y3 - y2) + (x2 * x2 + y2 * y2) * (y1 - y3) + (x3 * x3 + y3 * y3) * (y2 - y1);
    double C = (x1 * x1 + y1 * y1) * (x2 - x3) + (x2 * x2 + y2 * y2) * (x3 - x1) + (x3 * x3 + y3 * y3) * (x1 - x2);
    double D = (x1 * x1 + y1 * y1) * (x3 * y2 - x2 * y3) + (x2 * x2 + y2 * y2) * (x1 * y3 - x3 * y1) + (x3 * x3 + y3 * y3) * (x2 * y1 - x1 * y2);

    double x = -B / (2 * A);
    double y = -C / (2 * A);
    double r = sqrt((B * B + C * C - 4 * A * D) / (4 * A * A));

    return {{x, y}, r};
}

double fitCircle(vector<TVector3> &smearedHits, double resolution)
{
    auto chi2Function = [&](const double *par)
    {
        double chi2 = 0;
        for (auto hit : smearedHits)
        {
            double dx = hit.Z() - par[0];
            double dy = hit.Y() - par[1];
            double r = TMath::Sqrt(dx * dx + dy * dy);
            double k = 1 / r;
            double dk = par[2] - k;
            double err = dy / TMath::Power(r, 3) * resolution * 1e3;
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
    double radius = -1 / pFit[2] * pFit[1] / abs(pFit[1]);

    if (debug)
    {
        cout << "Fitted radius: " << radius << endl;

        // Plot the smeared hits and the fitted circle
        TCanvas *fitCircleCanvas = new TCanvas("fitCircleCanvas", "Fit Circle", 600, 600);
        fitCircleCanvas->cd();
        fitCircleCanvas->SetGrid();
        fitCircleCanvas->DrawFrame(-0.5, -0.8, 1.1, 0.8);
        TGraph *fitCircleGraph = new TGraph(smearedHits.size());
        for (int i = 0; i < smearedHits.size(); i++)
            fitCircleGraph->SetPoint(i, smearedHits[i].Z(), smearedHits[i].Y());
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

double fitLine(vector<TVector3> &smearedHits)
{
    TGraph *fitLineGraph = new TGraph(smearedHits.size());
    for (int i = 0; i < smearedHits.size(); i++)
        fitLineGraph->SetPoint(i, smearedHits[i].Z(), smearedHits[i].X());

    TF1 *fitLineFunction = new TF1("fitLineFunction", "[0] * x", 1);
    fitLineGraph->Fit(fitLineFunction, "Q");

    double slope = fitLineFunction->GetParameter(0);

    if (debug)
    {
        TCanvas *fitLineCanvas = new TCanvas("fitLineCanvas", "Fit Line", 600, 600);
        fitLineCanvas->cd();
        fitLineCanvas->SetGrid();
        fitLineCanvas->DrawFrame(-0.1, -0.04, 0.7, 0.04);
        fitLineGraph->SetMarkerStyle(20);
        fitLineGraph->SetMarkerSize(1);
        fitLineGraph->Draw("P");

        fitLineCanvas->Print("fitLine.png");

        delete fitLineCanvas;
    }

    delete fitLineGraph;
    delete fitLineFunction;

    return slope;
}

double reconstruction(vector<TVector3> &smearedHits, double resolution)
{
    // Calculate the radius of the particle's motion
    double r = fitCircle(smearedHits, resolution);

    if (r == 1e100)
        return r;

    // Calculate the momentum of the particle
    double momentumyz = r * charge * magneticField * c / 1e9;

    double k = fitLine(smearedHits);

    double momentum = TMath::Sqrt(1 + k * k) * momentumyz;

    // cout << momentum << " " << k * momentumyz << endl;

    return momentum;
}

double simulate(double momentum, double resolution)
{
    resolution *= 1e-6; // Convert the resolution to meters

    vector<TVector3> hits = propagate(momentum);

    vector<TVector3> smearedHits = gaussianSmearing(hits, resolution);

    if (debug)
    {
        cout << "Smeared hits:" << endl;
        for (auto hit : smearedHits)
            hit.Print();
    }

    return reconstruction(smearedHits, resolution);
}

vector<PDD> positionMomentumResolution(vector<double> &momentums, double resolution)
{
    vector<PDD> results;

    TCanvas *positionMomentumCanvas = new TCanvas("positionMomentumCanvas", "Position-Momentum", 1414, 1000);

    for (double momentum : momentums)
    {
        TH1F *positionMomentumHistogram = new TH1F("positionMomentum", Form("Momentum Distribution  (Resolution = %g #mum);Momentum (GeV/c);Counts", resolution),
                                                   100, 0, 0);

        for (int i = 0; i < N; ++i)
        {
            double reconstructedMomentum = simulate(momentum, resolution);
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
