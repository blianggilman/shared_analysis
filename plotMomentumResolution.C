
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "insert.h"


//Global Variables
TTree *tree = (TTree*) file->Get("tracks");

float gpx, gpy, gpz, px, py, pz;

double etavals[15] = {-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
char etachars[15][25] = {"n35", "n30", "n25", "n20", "n15", "n10", "n05", "00", "05", "10", "15", "20", "25", "30", "35"};
char etabins[14][25] = {"-3.5 ;< #eta < -3.0", "-3.0 < #eta < -2.5", "-2.5 < #eta < -2.0", 
    "-2.0 < #eta < -1.5", "-1.5 < #eta < -1.0", "-1.0 < #eta < -0.5", "-0.5 < #eta < 0.0", 
    "0.0 < #eta < 0.5", "0.5 < #eta < 1.0", "1.0 < #eta < 1.5", "1.5 < #eta < 2.0", 
    "2.0 < #eta < 2.5", "2.5 < #eta < 3.0", "3.0 < #eta < 3.5"};

char pbins[10][20] = {"0 < p < 2", "2 < p < 4", "4 < p < 6", "6 < p < 8", "8 < p < 10", "10 < p < 12",
    "12 < p < 14", "14 < p < 16", "16 < p < 18", "18 < p < 20"};

//first need to initialize histogram
vector<vector<TH1F*>> mom_res_preliminary_hists(14);
vector<vector<TH1F*>> mom_res_hists(14);
vector<vector<TH1F*>> mom_res_hists_direct(14);

int canvas_ctr = 0;
// int num_eta_bins = 14;
// int num_pt_bins = 10;

void initializeHists(){
    for (int i=0; i<14; i++){
        mom_res_preliminary_hists[i] = vector<TH1F*>(10);
        mom_res_hists[i] = vector<TH1F*>(10);
        mom_res_hists_direct[i] = vector<TH1F*>(10);

        for(int j=0; j<10; j++){
            // cout << "iteration: " << i << ", " << j << endl;

            char title[1024];
            sprintf(title, "mom_res_eta%s-%s_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            char label[1024];
            sprintf(label, ";dp/p, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            mom_res_preliminary_hists[i][j] = new TH1F(title,label,100,-1,1);

        }
    }
}




//calculate the three vector momentum
float calculateP(float px, float py, float pz){
    float sum = pow(px,2) + pow(py,2) + pow(pz,2);
    return sqrt(sum);
}


//calculate eta from theta
float calculateEta(float px, float py, float pz){
    //to find theta: theta is the COM scattering angle = the angle between the z-axis and 2D momentum
    //cos theta = pz/p
    float p = calculateP(px, py, pz);
    float theta = acos(pz/p);

    //given theta, calculate eta = -ln(tan(theta/2))
    return -log(tan(theta/2.));
}

//calculate the momentum resolution
float calculateMomRes(float gpx, float gpy, float gpz, float px, float py, float pz){
    //take the magnitude of each momentum vector, then subtract and take the absolute value
    float p = calculateP(px, py, pz);
    float gp = calculateP(gpx, gpy, gpz);
    float num = p-gp;
    return (num/abs(gp));
}


//make gaussian fits
Double_t** makeFits(vector<vector<TH1F*>> hists, Double_t** st_dev=NULL){

    double lowerbound = -0.1;
    double upperbound = 0.1;

    // TF1* g1 = new TF1("g1", "gaus",lowerbound,upperbound);
    Double_t** st_dev_new = 0;
    st_dev_new = new double*[14]; //[eta range][p range] in order from min to max (neg->pos)
    for (int i=0; i<14; i++){
        st_dev_new[i] = new double[10];
        for (int j=0; j<10; j++){
            if (st_dev != NULL){
                lowerbound = st_dev[i][j]*-3.;
                upperbound = st_dev[i][j]*3.;
            }
            TF1* g1 = new TF1("g1", "gaus",lowerbound,upperbound);

            hists[i][j]->Fit(g1, "RQ");
            st_dev_new[i][j] = g1->GetParameter(2);
        }
    }


    // // print st_dev_new values
    // for (int i=0; i<10; i++){
    //     cout << st_dev_new[0][i] << endl;
    // }

    return st_dev_new;


}



//calculate PWG requirements
double* calculatePWGreqs(int i, double p[]){
    double A[] = {0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1};
    double B[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1., 1., 1., 2., 2.};
    double* dp_p = 0;
    dp_p = new double[10];
    
    for (int j=0; j<10; j++){
        double sum = pow(A[i]*p[j],2) + pow(B[i],2);
        dp_p[j] = sqrt(sum);
    }

    return dp_p;

}


//plot SD results
void plotSD(Double_t** st_dev, Double_t** st_dev_direct){
    
    //initialize array of TGraphs
    std::vector<TGraph*> pwg_req_eqs(14);
    std::vector<TGraph*> momPlots_by_p(14);
    std::vector<TGraph*> momPlots_by_p_direct(14);
    std::vector<TGraph*> momPlots_by_eta(14);
    const Int_t n = 10;
    double p[n] = {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.};
    double etaforplot[14] = {-3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25};
    

    //make eta labels
    std::vector<TLatex*> etalabels(14);
    // cout << "length of eta " << sizeof(etabins) << " " << sizeof(etalabels) << endl;
    for (int i=0; i<14; i++) {
        double* dp_p = calculatePWGreqs(i, p);
        double* max;

        cout << "terms in dp/p: ";
        for (int j=0; j<10; j++){
            // max = max_element(dp_p[0], dp_p[j]+10);
            cout << dp_p[j] << " ";
        }

        double width = 16.;
        double height = (dp_p[0]+dp_p[9])/2;
        etalabels[i] = new TLatex(width,height,etabins[i]);
        cout << "MAXIMUM!!!! " << *max << ", " << height << endl;

        // cout << "eta label " << i << ": " << etalabels[i] << endl;
    }

    //make p labels
    std::vector<TLatex*> plabels(10);
    for (int i=0; i<10; i++) {
        double width = -1.;
        double height = 2.;
        char name[1024];
        sprintf(name, "%.0f < p < %.0f", p[i], p[i+1]);
        // plabels[i] = new TLatex(width,height,name); //"? ;< p < ?"
        plabels[i] = new TLatex(width,height,pbins[i]);
    }



    //fill graphs and plot, also get PWG requirements
    TCanvas *c15 = new TCanvas("c15","c15",1200,900);

    //if the first p bin in each eta range has less than 100 events, exclude that eta range
    if (canvas_ctr <= 12) c15->Divide(3,4);
    else c15->Divide(3,5);
    
    for(int i=0; i<14; i++){ //CHANGE UPPER BOUND TO CANVAS_CTR
        // cout << "iteration: " << i << endl;
        // cout << st_dev[1][1] << ", " << st_dev[5][5] << endl;

        if (i==0 || i==13) continue;

        //convert st dev's into percentages
        for (int j=0; j<10; j++){
            st_dev[i][j] *= 100.0;
            st_dev_direct[i][j] *= 100.0;
        }
        
        momPlots_by_p[i] = new TGraph(n,p,st_dev[i]);
        momPlots_by_p_direct[i] = new TGraph(n,p,st_dev_direct[i]);

        double* dp_p = calculatePWGreqs(i, p);
        pwg_req_eqs[i] = new TGraph(n,p,dp_p);

        c15->cd(i); //CHANGED FROM ORIGINAL i+1
        TMultiGraph *mg = new TMultiGraph();
        mg->Add(momPlots_by_p[i]);
        mg->Add(pwg_req_eqs[i]);
        mg->Add(momPlots_by_p_direct[i]);
        mg->Draw("ALP"); //"ACP"

        momPlots_by_p[i]->SetLineColor(4);
        momPlots_by_p[i]->SetMarkerColor(4);
        momPlots_by_p[i]->SetMarkerStyle(2);
        momPlots_by_p_direct[i]->SetLineColor(2);
        momPlots_by_p_direct[i]->SetMarkerColor(2);
        momPlots_by_p_direct[i]->SetMarkerStyle(5);

        mg->GetXaxis()->SetTitle("Track p [GeV/c]");
        mg->GetYaxis()->SetTitle("#Deltap/p (%)");
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();

        etalabels[i]->SetTextFont(43); etalabels[i]->SetTextSize(12);
        etalabels[i]->Draw("same");

    }

    c15->Print("plots/mom_res_SD_by_p.pdf");



    //fill graphs and plot, also get PWG requirements
    TCanvas *c16 = new TCanvas("c16","c16",1200,900);
    c16->Divide(5,2);
    
    for(int i=0; i<10; i++){

        //save each column of st_dev 2d array into an array
        Double_t st_dev_col[14];
        for (int a=0; a<14; a++){
            st_dev_col[a] = st_dev[a][i];
            cout << st_dev[a][i] << ", ";
        }
        cout << endl;

        momPlots_by_eta[i] = new TGraph(14,etaforplot,st_dev_col);
        // for (int b=0; b<14; b++){
        //     cout << b << ": (" << etaforplot[b] << ", " << st_dev_col[b] << ")" << endl;
        // }
        if (i==6) {
            for (int b=0; b<14; b++){
                cout << b << ": (" << etaforplot[b] << ", " << st_dev_col[b] << ")" << endl;
            }
        }


        c16->cd(i+1); 
        momPlots_by_eta[i]->Draw("ALP"); //"ACP"

        momPlots_by_eta[i]->SetLineColor(4);
        momPlots_by_eta[i]->SetMarkerColor(2);
        momPlots_by_eta[i]->SetMarkerStyle(2);

        momPlots_by_eta[i]->GetXaxis()->SetTitle("Track #eta [rad]");
        momPlots_by_eta[i]->GetYaxis()->SetTitle("#Deltap/p (%)");
        momPlots_by_eta[i]->GetXaxis()->CenterTitle();
        momPlots_by_eta[i]->GetYaxis()->CenterTitle();

        plabels[i]->SetTextFont(43); plabels[i]->SetTextSize(12);
        plabels[i]->Draw("same");

    }

    c16->Print("plots/mom_res_SD_by_eta.pdf");







}



//plot histogram of dp/p in multiple bins of eta and p (see graph for specifics)
//eta: [-3.5,3.5] in bins of 0.5
//p: [0,10] in bins of 1 or 2? //try 1 first
void plotMomRes(int nEntries){    

    int ctr[14];
    int max[14];
    int savei[14];


    //fill histograms based on the bins
    for (int i=0; i<nEntries; i++) {
        if (i%1000000==0) cout << "jet: " << i << " out of: " << nEntries << endl;
        
        tree->GetEntry(i);
        float eta = calculateEta(px,py,pz);
        float p = calculateP(px,py,pz);
        float momRes = calculateMomRes(gpx,gpy,gpz,px,py,pz);
        

        if (i<10) cout << "p for " <<i << " is " << px << ", " << py << ", " << pz << " --> " << p << endl;
        if (i<10) cout << "eta for " <<i << " is " << eta << endl;
        if (i<10) cout << "theta for " <<i << " is " << acos(pz/p) << endl;

        //fill hists
        for (int i=0; i<14; i++){
            for (int j=0; j<10; j++){
                if (eta > etavals[i] && eta <= etavals[i+1]) {
                    ctr[i]++;
                    if (momRes > max[i]) {
                        max[i] = momRes;
                        savei[i] = i;
                    }
                    if (p > j*2 && p <= j*2+2) mom_res_preliminary_hists[i][j]->Fill(momRes); //does comparison have to be a double?
                }
            }
        }
        


    }

    //if the first p bin in each eta range has less than 100 events, exclude that eta range
    for (int i=0; i<14; i++){
        if (mom_res_preliminary_hists[i][0]->GetEntries() > 100) canvas_ctr++;
    }

    //print statements
    cout << "COUNTERS" << endl;
    for(int i=0; i<14; i++){
        cout << i << ": " << ctr[0] << endl;
    }
    cout << "MAX VALS PER ETA (event#: max resolution)" << endl;
    for(int i=0; i<14; i++){
        cout << savei[i] << ": " << max[i] << endl;
    }



    //calculate standard deviations
    Double_t** st_dev_direct = 0;
    st_dev_direct = new double*[14]; //[eta range][p range] in order from min to max (neg->pos)
    for (int i=0; i<14; i++){
        st_dev_direct[i] = new double[10];
        for (int j=0; j<10; j++){
            st_dev_direct[i][j] = mom_res_preliminary_hists[i][j]->GetStdDev();
        }
    }


    //make fits
    Double_t** st_dev_fromfit = makeFits(mom_res_preliminary_hists);

    //print statements
    cout << "compare st devs!!!!" << endl;
    for (int i=0; i<14; i++){
        for (int j=0; j<10; j++){
            cout << st_dev_direct[i][j] << " ?= " << st_dev_fromfit[i][j] << endl;
        }
    }
    

    //redefine histograms with new binning from std dev w fits and calculated directly
    for (int i=0; i<14; i++){
        for(int j=0; j<10; j++){
            char title[1024];
            sprintf(title, "mom_res_eta%s-%s_p%d-%d_rebinned", etachars[i], etachars[i+1], j*2, j*2+2);
            char label[1024];
            sprintf(label, ";dp/p, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            mom_res_hists[i][j] = new TH1F(title, label, 100, -3*st_dev_fromfit[i][j], 3*st_dev_fromfit[i][j]);

            //calculated directly
            char title_direct[1024];
            sprintf(title_direct, "mom_res_eta%s-%s_p%d-%d_direct_rebinned", etachars[i], etachars[i+1], j*2, j*2+2);
            mom_res_hists_direct[i][j] = new TH1F(title_direct, label, 100, -3*st_dev_direct[i][j], 3*st_dev_direct[i][j]);
        }
    }


    //fill histograms with new binning from std dev w fits
    for (int i=0; i<nEntries; i++) {
        if (i%1000000==0) cout << "jet: " << i << " out of: " << nEntries << endl;
        
        tree->GetEntry(i);
        float eta = calculateEta(px,py,pz);
        float p = calculateP(px,py,pz);
        float momRes = calculateMomRes(gpx,gpy,gpz,px,py,pz);

        for (int i=0; i<14; i++){
            for (int j=0; j<10; j++){
                if (eta > etavals[i] && eta <= etavals[i+1]) {
                    if (p > j*2 && p <= j*2+2) {
                        mom_res_hists[i][j]->Fill(momRes); //does comparison have to be a double?
                        mom_res_hists_direct[i][j]->Fill(momRes);
                    }
                }
            }
        }
    }

    Double_t** st_dev_final_fromfit = makeFits(mom_res_hists, st_dev_fromfit);
    Double_t** st_dev_final_direct = makeFits(mom_res_hists, st_dev_direct);



    //make canvas to draw
    vector<TCanvas*> canvas(14);
    // vector<TCanvas*> canvas_direct(14);
    for (int i=0; i<14; i++){
        char canvas_title[1024];
        sprintf(canvas_title, "c%d", i);
        canvas[i] = new TCanvas(canvas_title,canvas_title,1200,900);
        canvas[i]->Divide(5,2);
        // canvas_direct[i] = new TCanvas(canvas_title,canvas_title,1200,900);
        // canvas_direct[i]->Divide(5,2);

    }


    for (int i=0; i<14; i++){
        for (int j=0; j<10; j++){
            canvas[i]->cd(j+1);
            // mom_res_preliminary_hists[i][j]->Draw();
            mom_res_hists[i][j]->Draw();

            // canvas_direct[i]->cd(j+q);
            // mom_res_hists_direct[i]->Draw();

        }
        char outname[1024];
        sprintf(outname, "plots/mom_res_etan%s-%s_bins100_rebinned.pdf", etachars[i], etachars[i+1]);
        canvas[i]->Print(outname);
    }


    //plot SD
    plotSD(st_dev_final_fromfit, st_dev_final_direct);

}




//combine - main method
void plotMomentumResolution(){
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    initializeHists();

    tree->SetBranchAddress("gpx", &gpx);
    tree->SetBranchAddress("gpy", &gpy);
    tree->SetBranchAddress("gpz", &gpz);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);

    
    const int nEntries = tree->GetEntries();
    plotMomRes(nEntries); 

}