#include "standalonenoise.h"

//
// stand alone code to produce a noise spectrum which mimics the noise in the 
// MicroBooNE detector for a given wirelength
//
// MicroBooNE's effective number of bits (ENOB) is ~10.9, with a waveform size 
// of 9596 ticks. This corresponds to a sideband magnitude of ~60 in the 
// feequency domain
//
// Author: Adam Lister
// Email : a.lister1@lancaster.ac.uk
//

void GenNoisePostFilter(double wirelength, double ENOB, TFile* f_output)
{

    // 
    // setup
    //

    const size_t waveform_size = 9596; 

    // initialise histograms

    TH1D* noiseFreqDom = new TH1D("noiseFreqDom", "", waveform_size, 0, waveform_size); 
    TH1D* noiseTimeDom = new TH1D("noiseTimeDom", "", waveform_size, 0, waveform_size);
    TH1D* adcs = new TH1D("adcs", "", 20, -10, 10);
    
    TString wl = Form("_%gcm", wirelength);
    noiseTimeDom->SetName("noiseTimeDom"+wl);
    noiseFreqDom->SetName("noiseFreqDom"+wl);
    adcs->SetName("adcs"+wl);

    // vars

    std::vector<double> _pfn_rho_v;
    std::vector<double> _pfn_value_re;
    std::vector<double> _pfn_value_im;
    _pfn_rho_v.resize(waveform_size);
    _pfn_value_re.resize(waveform_size);
    _pfn_value_im.resize(waveform_size);
    Double_t fitpar[9] = {0.};
    Double_t wldparams[2] = {0.};

    // calculate FFT magnitude of noise from ENOB
    
    double baseline_noise = std::sqrt(waveform_size/12)*std::pow(2, 12 - ENOB);
 
    // wire length dependence function 
    
    TF1* _wld_f = new TF1("_wld_f", "[0] + [1]*x", 0.0, 470);

    // gain function
    
    TF1* _pfn_f1 = new TF1("_pfn_f1", "((([0]*1/(x*[8]/2 + 10) + [1]*exp(-0.5*(((x*[8]/2)-[2])/[3])**2)*exp(-0.5*pow(x*[8]/(2*[4]),[5])))*[6])+[7])", 0.0, (double)waveform_size);

    // custom poisson
    
    TF1* _poisson = new TF1("_poisson", "[0]**(x) * exp([-0]) / tgamma(x+1.)", 0.0, 30.0);
  
    //
    // set data-driven parameters
    //

    // posson mean
    const double params = 3.30762;

    _poisson->SetParameter(0, params);
    
    wldparams[0] = 0.395;
    wldparams[1] = 0.001304;

    _wld_f->SetParameters(wldparams);
    double wldValue = _wld_f->Eval(wirelength);

    fitpar[0] = 9.27790e+02;
    fitpar[1] = 1.20284e+07;
    fitpar[2] = 4.93692e+03;
    fitpar[3] = 1.03438e+03;
    fitpar[4] = 2.33306e+02;
    fitpar[5] = 1.36605e+00;
    fitpar[6] = wldValue;
    fitpar[7] = baseline_noise;
    fitpar[8] = waveform_size;

    _pfn_f1->SetParameters(fitpar);


    // 2MHz Digitization

    for(size_t i=0; i < waveform_size; i++){

        Double_t freq;
        if (i < waveform_size/2.){
   
            freq = (i)*2./waveform_size; 
        
        }
        else{
        
            freq = (waveform_size-i)*2./waveform_size;
        
        }

        double pfnf1val = _pfn_f1->Eval(freq);
     

        // define FFT parameters
        double randomizer = _poisson->GetRandom()/params;


        _pfn_rho_v[i] = pfnf1val * randomizer;

        Double_t rho = _pfn_rho_v[i];
        Double_t phi = gRandom->Uniform(0,1) * 2. * TMath::Pi();


        _pfn_value_re[i] = rho*cos(phi)/((double)waveform_size);
        _pfn_value_im[i] = rho*sin(phi)/((double)waveform_size);

        noiseFreqDom->SetBinContent(i, rho);

    }

    //
    // Inverse FFT
    //
    
    TVirtualFFT::SetTransform(0);
    int n = waveform_size;
    TVirtualFFT* _pfn_ifft = TVirtualFFT::FFT(1,&n,"MAG C2R M K");
    _pfn_ifft->SetPointsComplex(&_pfn_value_re[0],&_pfn_value_im[0]);
    _pfn_ifft->Transform();

    TH1 *fb = 0;
    fb = TH1::TransformHisto(_pfn_ifft,fb,"Re");

    //
    // Define noise spectrum
    //

    for(size_t i=0; i<waveform_size; ++i) {
        noiseTimeDom->SetBinContent(i, fb->GetBinContent(i));
        adcs->Fill(fb->GetBinContent(i));
    }

    delete fb;

    f_output->cd();
    adcs->Write();
    noiseTimeDom->Write();
    noiseFreqDom->Write();
}

int main(){

    TFile *f_output = new TFile("noiseHistograms.root", "recreate");

    double ENOB = 10.91;
   
    for (int i; i < 500; i++){
        GenNoisePostFilter(i, ENOB, f_output);
    }

    f_output->Close();
    
    return 0;

}
