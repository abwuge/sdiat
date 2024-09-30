#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TRandom.h>

int main()
{
    TH1D *hk = new TH1D("hk", ";k;Events", 100, 0, 0);
    TH1D *habsk = new TH1D("habsk", ";1/k;Events", 100, 0, 0);

    int n = 100000;
    double k;
    double mean = 1,
           sigma = 1;
    for (int i = 0; i < n; i++)
    {
        k = gRandom->Gaus(mean, sigma);
        hk->Fill(k);
        habsk->Fill(fabs(k));
    }
    TCanvas *c1 = new TCanvas("c1", "c1", 1414, 1000);

    hk->Draw();
    c1->SaveAs("test1.png");

    habsk->Draw();
    c1->SaveAs("test2.png");

    return 0;
}
