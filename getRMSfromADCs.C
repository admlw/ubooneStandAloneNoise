void getRMSfromADCs() {

    TH2D* lengthrms = new TH2D("lengthrms", "", 470, 0, 470, 100, 1, 3);

    for (int i = 0; i <470; i++){

        TString plotname = Form("adcs_%icm", i);
        TH1D *h = (TH1D*)_file0->Get(plotname);

        //Manual calculation of RMS on h_adc
        double pars[3];
        double rms;
        if (h->GetSum()>0){
            double xc = 0.5-0.34;
            h->GetQuantiles(1, &pars[0], &xc);

            xc = 0.5;
            h->GetQuantiles(1, &pars[1], &xc);

            xc = 0.5+0.34;
            h->GetQuantiles(1, &pars[2], &xc);

            rms = sqrt((pow(pars[1]-pars[0],2) + pow(pars[2]-pars[1],2))/2.);

        }


        lengthrms->Fill(i, rms);

        h->Delete();
    }

    lengthrms->Draw();

}
