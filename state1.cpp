#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLine.h>

int state1();

using namespace std;

int main()
{
    return state1();
}

int state1()
{
    // Constants
    const double c = 299792458; // Speed of light in m/s

    // Particle properties
    // Particle is a proton, moving in the positive X-axis
    double charge = 1;      // Charge of the particle in e
    double mass0 = 938e6;   // Rest mass of the particle in eV/c^2
    double momentum0 = 1e9; // Initial momentum of the particle in eV/c
    double x = 0;           // Initial x position in meters
    double y = 0;           // Initial y position in meters

    // Magnetic field properties
    // Magnetic field is in the positive Z-axis
    double magneticField = 1; // Magnetic field strength in Tesla

    // Detector properties
    double detectorGap = 15e-2; // The gap between each pair of detectors
    int numDetectors = 5;       // The number of detectors

    // Time properties
    double timeStep = 1e-6; // Time step for the simulation in seconds

    // Calculate
    double energy = TMath::Sqrt(TMath::Power(momentum0, 2) + TMath::Power(mass0, 2)); // Total energy of the particle in eV
    double mass = energy;                                                             // Mass of the particle in eV/c^2
    double gamma = energy / mass0;                                                    // Lorentz factor
    double momentumX = momentum0;                                                     // Momentum of the particle in eV/c
    double momentumY = 0;                                                             // Momentum of the particle in eV/c
    double velocityX = momentumX / mass;                                              // X component of the velocity in c
    double velocityY = momentumY / mass;                                              // Y component of the velocity in c
    double lastVelocity = TMath::Sqrt(TMath::Power(velocityX, 2) + TMath::Power(velocityY, 2));

    // Print the initial conditions
    printf("Initial velocityX: %g m/s, velocityY: %g m/s\n", velocityX * c, velocityY * c);

    // Simulate the particle's motion
    vector<double> xx, yy; // Vectors to store the position of the particle
    xx.push_back(x), yy.push_back(y);

    // Calculate the radius of the particle's motion
    double R = (energy - mass0) / (abs(charge) * magneticField) / c;
    printf("R: %g m\n", R);

    double maxX = 0, maxY = 0;
    int detectorIndex = 0;
    int WatchDog = 0;
    while (x <= detectorGap * (numDetectors - 1))
    // while (true)
    {
        // WatchDog
        WatchDog++;
        if (WatchDog > 1e6)
        {
            cout << "WatchDog triggered" << endl;
            break;
        }

        // Calculate the new position using recurrence relations
        double forceX = -c * charge * velocityY * magneticField;                                          // Lorentz force in the X direction in eV/m
        double forceY = c * charge * velocityX * magneticField;                                           // Lorentz force in the Y direction in eV/m
        double velocity = TMath::Sqrt(TMath::Power(velocityX, 2) + TMath::Power(velocityY, 2));           // Velocity of the last time step in c
        double deltaVelocity = velocity - lastVelocity;                                                   // Change in velocity in c
        double accelerationX = forceX / mass - velocity * deltaVelocity * velocityX / (c * c * timeStep); // Acceleration in the X direction in c^2/m
        double accelerationY = forceY / mass - velocity * deltaVelocity * velocityY / (c * c * timeStep); // Acceleration in the Y direction in c^2/m
        double lastX = x, lastY = y;
        x = x + velocityX * timeStep + 0.5 * accelerationX * timeStep * timeStep; // Update the X position
        y = y + velocityY * timeStep + 0.5 * accelerationY * timeStep * timeStep; // Update the Y position
        lastVelocity = velocity;                                                  // Update the last velocity
        velocityX = velocityX + accelerationX * timeStep;                         // Update the X component of the velocity
        velocityY = velocityY + accelerationY * timeStep;                         // Update the Y component of the velocity

        if (lastX <= detectorGap * detectorIndex && x >= detectorGap * detectorIndex)
        {
            printf("Detected particle at x = %g m, y = %g m, t = %g s\n", x, y, timeStep * WatchDog);
            detectorIndex++;
        }

        if (detectorIndex >= numDetectors)
        {
            break;
        }

        if (x > maxX)
            maxX = x;
        if (y > maxY)
            maxY = y;
        xx.push_back(x), yy.push_back(y); // Store the new position
    }

    // Create a TCanvas object to display the plot
    TCanvas *canvas = new TCanvas("canvas", "Particle Trajectory", 800, 600);

    // Create a TGraph object to plot the trajectory
    TGraph *graph = new TGraph(xx.size(), &xx[0], &yy[0]);

    graph->Draw("AL");
    graph->SetTitle(Form("Particle Trajectory (Momentum: %g GeV/c)", momentum0 / 1e9));
    graph->GetXaxis()->SetTitle("X position (m)");
    graph->GetYaxis()->SetTitle("Y position (m)");

    // Draw the detectors
    for (int i = 0; i < numDetectors; i++)
    {
        double detectorX = i * detectorGap;

        int t = upper_bound(xx.begin(), xx.end(), detectorX) - xx.begin();
        printf("Detected particle at x = %g m, y = %g m, t = %g s\n", xx[t], yy[t], t * timeStep);
        TLine *detectorLine = new TLine(detectorX, 0, detectorX, maxY * 1.1);
        detectorLine->SetLineColor(kRed);
        detectorLine->Draw();
    }

    TLine *RLine = new TLine(R, 0, R, maxY * 1.1);
    RLine->SetLineColor(kGreen);
    RLine->Draw();

    printf("R: %g m (Calculate), %g m (Simulate)\n", R, maxX);

    // Draw the canvas to display the plot on screen
    canvas->Draw();

    // Save the plot as a PDF file
    canvas->Print("trajectory.png");

    // // Clean up
    delete graph;
    delete canvas;

    return 0;
}
