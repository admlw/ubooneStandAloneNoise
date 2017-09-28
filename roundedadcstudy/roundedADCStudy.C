#include "TH1D.h"
#include "TFile.h"
#include <iostream>

TH1D* roundAdcs(TH1D* th1d, int roundingVal, int waveform_size, TH1D* roundedAdcs){
    
    th1d->Scale(1./(double)roundingVal);

    for (int i = 1; i < waveform_size+1; i++){
        int roundedVal = std::floor(th1d->GetBinContent(i) + 0.5);
        roundedAdcs->SetBinContent(i, roundedVal);
    }

    roundedAdcs->Scale(roundingVal);
    return roundedAdcs;

}

void roundedADCStudy(){

    TFile fout("outputHistos.root", "update");
    TFile fin("noiseHistograms.root", "read");

    int roundingVal = 8;

    TString outputHistoName = Form("roundedto%iADCs", roundingVal);
    TH1D* summedFreqSpace = new TH1D("summedFreqSpace", "", 9594, 0, 9594);
    summedFreqSpace->SetName(outputHistoName);

    for (int i = 1; i < 501; i ++){

        TString histoName = Form("noiseTimeDom_1000cm;%i", i);
        fin.cd();
        TH1D* noiseHisto = (TH1D*)fin.Get(histoName);
        TH1D* roundedHisto = new TH1D("roundedHisto", "", 9594, 0, 9594); 
        roundedHisto = roundAdcs(noiseHisto, roundingVal, 9594, roundedHisto);

        TH1* freqDom = 0;
        freqDom = (TH1*)roundedHisto->FFT(freqDom, "MAG RTC M");

        summedFreqSpace->Add(freqDom);

        noiseHisto->Delete();
        roundedHisto->Delete();
        freqDom->Delete();
    }

    fout.cd();
    summedFreqSpace->Write();
    fout.Close();
    

}

int main(){

    roundedADCStudy();
    return 0;

}
