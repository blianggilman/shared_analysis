#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

#include "insert.h"


//Global Variables
TTree *tree = (TTree*) file->Get("tracks");

float gpx, gpy, gpz, px, py, pz;
float dca2d, pcaz, gvz;

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

vector<vector<TH1F*>> phi_res_preliminary_hists(14);
vector<vector<TH1F*>> theta_res_preliminary_hists(14);
vector<vector<TH1F*>> longDCA_res_preliminary_hists(14);
vector<vector<TH1F*>> transDCA_res_preliminary_hists(14);

vector<vector<TH1F*>> phi_res_hists(14);
vector<vector<TH1F*>> theta_res_hists(14);
vector<vector<TH1F*>> longDCA_res_hists(14);
vector<vector<TH1F*>> transDCA_res_hists(14);


// int num_eta_bins = 14;
// int num_pt_bins = 10;



void pretty_TGraph( TGraph * g, int color = 1, int marker = 20, TString xtitle = "", TString ytitle = "", TString title="");

void initializeHists(){
    for (int i=0; i<14; i++){
        mom_res_preliminary_hists[i] = vector<TH1F*>(10);
        mom_res_hists[i] = vector<TH1F*>(10);
        mom_res_hists_direct[i] = vector<TH1F*>(10);

        
        phi_res_preliminary_hists[i] = vector<TH1F*>(10);
        theta_res_preliminary_hists[i] = vector<TH1F*>(10);
        longDCA_res_preliminary_hists[i] = vector<TH1F*>(10);
        transDCA_res_preliminary_hists[i] = vector<TH1F*>(10);

        phi_res_hists[i] = vector<TH1F*>(10);
        theta_res_hists[i] = vector<TH1F*>(10);
        longDCA_res_hists[i] = vector<TH1F*>(10);
        transDCA_res_hists[i] = vector<TH1F*>(10);


        for(int j=0; j<10; j++){
            // cout << "iteration: " << i << ", " << j << endl;

            char title[1024];
            sprintf(title, "mom_res_eta%s-%s_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            char label[1024];
            sprintf(label, ";dp/p, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            mom_res_preliminary_hists[i][j] = new TH1F(title,label,100,-1,1);

            char phititle[1024];
            sprintf(phititle, "phi_res_eta%s-%s_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            char philabel[1024];
            sprintf(philabel, ";d#phi, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            phi_res_preliminary_hists[i][j] = new TH1F(phititle,philabel,100,-.1,.1);

            char thetatitle[1024];
            sprintf(thetatitle, "theta_res_eta%s-%s_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            char thetalabel[1024];
            sprintf(thetalabel, ";d#theta, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            theta_res_preliminary_hists[i][j] = new TH1F(thetatitle,thetalabel,100,-1,4);

            char longDCAtitle[1024];
            sprintf(longDCAtitle, "longDCA_res_eta%s-%s_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            char thetaDCAlabel[1024];
            sprintf(thetaDCAlabel, ";Longitudinal DCA, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            longDCA_res_preliminary_hists[i][j] = new TH1F(longDCAtitle,thetaDCAlabel,100,-175,175);

            char transDCAtitle[1024];
            sprintf(transDCAtitle, "transDCA_res_eta%s-%s_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            char transDCAlabel[1024];
            sprintf(transDCAlabel, ";Transverse DCA, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            transDCA_res_preliminary_hists[i][j] = new TH1F(transDCAtitle,transDCAlabel,100,-5,5);

        }
    }
}




//calculate the three vector momentum
float calculateP(float px, float py, float pz){
    float sum = pow(px,2) + pow(py,2) + pow(pz,2);
    return sqrt(sum);
}


//calculate phi
float calculatePhi(float px, float py, float pz){
    //phi = tan(py/px)???
    return atan(py/px);
}

//calculate theta
float calculateTheta(float px, float py, float pz){
    //to find theta: theta is the COM scattering angle = the angle between the z-axis and 2D momentum
    //cos theta = pz/p
    float p = calculateP(px, py, pz);
    float theta = acos(pz/p);

    return theta;
}


//calculate eta from theta
float calculateEta(float px, float py, float pz){
    //to find theta: theta is the COM scattering angle = the angle between the z-axis and 2D momentum
    //cos theta = pz/p
    float theta = calculateTheta(px, py, pz);

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

//calculate the angular theta resolution
float calculatePhiRes(float gpx, float gpy, float gpz, float px, float py, float pz){
    float phi = calculatePhi(px, py, pz);
    float gphi = calculatePhi(gpx, gpy, gpz);
    return phi-gphi;
}

//calculate the angular theta resolution
float calculateThetaRes(float gpx, float gpy, float gpz, float px, float py, float pz){
    float theta = calculateTheta(px, py, pz);
    float gtheta = calculateTheta(gpx, gpy, gpz);
    return theta-gtheta;
}

//calculate the longitudinal (z) DCA
float calculateLongDCA(float pcaz, float gvz){
    // float dvt = dca2d;
    float dvl = pcaz - gvz;
    return dvl;
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
        //convert pwg reqs into decimals (originally in percentages)
        dp_p[j] /= 100.0;
    }

    return dp_p;

}

void grrReadFile(){

    
}

//plot SD results
//flag 0 = mom res, 1 = phi res, 2 = theta res, 3 = longDCA res, 4 = transDCA
// void plotSD(Double_t** st_dev, Double_t** st_dev_direct, int canvas_ctr, int flag=0){
void plotSD(Double_t** st_dev, int canvas_ctr, int flag=0){
    
    /*//initialize array of TGraphs
    std::vector<TGraph*> pwg_req_eqs(14);
    std::vector<TGraph*> momPlots_by_p(14);
    std::vector<TGraph*> momPlots_by_p_direct(14);
    std::vector<TGraph*> momPlots_by_eta(14);*/
    std::vector<TGraph*> EICPlots(14);
    const Int_t n = 10;
    double p[n] = {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.};
    // double etaforplot[14] = {-3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25};
    static int canvas_ctr_static = canvas_ctr;
    double etaforplot[12] = {-2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75};


    if (flag==0){
        for (int i=0; i<14; i++) {

            if (i==0 || i==13) continue;

            double* dp_p = calculatePWGreqs(i, p);

            //save PWG requirements
            char newfilename[1024];
            sprintf(newfilename, "datafiles/%s/data_mom_res_PWG_req_eta%s-%s.dat", configuration, etachars[i], etachars[i+1]);

            ofstream outdata; // outdata is like cin
            outdata.open(newfilename); // opens the file
            if( !outdata ) { // file couldn't be opened
                cerr << "Error: file could not be opened" << endl;
                exit(1);
            }

            for (int j=0; j<n; j++)
                outdata << p[j] << " " << dp_p[j] << endl;
            outdata.close();
        }

    

        for(int i=0; i<14; i++){ //CHANGE UPPER BOUND TO CANVAS_CTR
            // cout << "iteration: " << i << endl;
            // cout << st_dev[1][1] << ", " << st_dev[5][5] << endl;

            if (i==0 || i==13) continue;

            Double_t* EIC_x = 0;
            EIC_x = new double[9];

            Double_t* EIC_y = 0;
            EIC_y = new double[9];

            char thingy[1024];
            sprintf(thingy, "dat%i.csv", i);
            cout << "file: " << thingy << endl;

            ifstream in(thingy);
            vector<vector<double>> fields;

            if (in) {
                string line;
                while (getline(in, line)) {
                    stringstream sep(line);
                    string field;
                    // cout << line << endl;

                    fields.push_back(vector<double>());

                    while (getline(sep, field, ',')) {
                        fields.back().push_back(stod(field));
                    }
                }
            }

            // for (auto row : fields) {
            //     for (auto field : row) {
            //         cout << field << ' '; //PRINTS
            //     }
            //     cout << '\n';
            // }

            for (int a=0; a<9; a++){
                EIC_x[a] = fields[a][0];
                EIC_y[a] = fields[a][1];
            }
            // EICPlots[i] = new TGraph(9, EIC_x, EIC_y);
            in.close();  
            in.clear();


            //save PWG requirements
            char ecceintnotefilename[1024];
            sprintf(ecceintnotefilename, "datafiles/data_mom_res_ECCE_intnote_eta%s-%s.dat", etachars[i], etachars[i+1]);

            ofstream outdata; // outdata is like cin
            outdata.open(ecceintnotefilename); // opens the file
            if( !outdata ) { // file couldn't be opened
                cerr << "Error: file could not be opened" << endl;
                exit(1);
            }

            for (int a=0; a<9; a++)
                outdata << fields[a][0] << " " << fields[a][1] << endl;
            outdata.close();


            /*
            momPlots_by_p[i] = new TGraph(n,p,st_dev[i]);
            momPlots_by_p_direct[i] = new TGraph(n,p,st_dev_direct[i]);

            double* dp_p = calculatePWGreqs(i, p);
            pwg_req_eqs[i] = new TGraph(n,p,dp_p);

            c15->cd(i); //CHANGED FROM ORIGINAL i+1
            TMultiGraph *mg = new TMultiGraph();
            mg->Add(momPlots_by_p[i]);
            mg->Add(EICPlots[i]);
            mg->Add(pwg_req_eqs[i]);
            mg->Add(momPlots_by_p_direct[i]);
            mg->Draw("ALP"); //"ACP"

            momPlots_by_p[i]->SetLineColor(4);
            momPlots_by_p[i]->SetMarkerColor(4);
            momPlots_by_p[i]->SetMarkerStyle(2);
            momPlots_by_p_direct[i]->SetLineColor(2);
            momPlots_by_p_direct[i]->SetMarkerColor(2);
            momPlots_by_p_direct[i]->SetMarkerStyle(5);
            EICPlots[i]->SetLineColor(8);
            EICPlots[i]->SetMarkerColor(8);
            EICPlots[i]->SetMarkerStyle(4);

            mg->GetXaxis()->SetTitle("Track p [GeV/c]");
            mg->GetYaxis()->SetTitle("#Deltap/p");
            mg->GetXaxis()->CenterTitle();
            mg->GetYaxis()->CenterTitle();
            mg -> GetXaxis() -> SetNdivisions(107);
            mg -> GetXaxis() -> SetTitleSize(0.05);
            mg -> GetXaxis() -> SetLabelSize(0.05);
            mg -> GetYaxis() -> SetNdivisions(107);
            mg -> GetYaxis() -> SetTitleSize(0.05);
            mg -> GetYaxis() -> SetLabelSize(0.05);

            etalabels[i]->SetTextFont(43); etalabels[i]->SetTextSize(18);
            etalabels[i]->Draw("same");
        
            auto legend = new TLegend(0.1,.7,0.48,0.9);
            legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
            legend->AddEntry(momPlots_by_p_direct[i],"Momentum Plots By P","l");
            legend->AddEntry(EICPlots[i],"ECCE Paper Plots","l");
            legend->AddEntry(pwg_req_eqs[i],"PWG Requirement Equations","l");
            legend->Draw();
            */
        }
    }



    for (int i=0; i<14; i++) {

        // if (i==0 || i==13) continue;


        //save data
        char configfilename[1024];
        if (flag == 0) sprintf(configfilename, "datafiles/%s/data_mom_res_sim_eta%s-%s.dat", configuration, etachars[i], etachars[i+1]);
        if (flag == 1) sprintf(configfilename, "datafiles/%s/data_phi_res_sim_eta%s-%s.dat", configuration, etachars[i], etachars[i+1]);
        if (flag == 2) sprintf(configfilename, "datafiles/%s/data_theta_res_sim_eta%s-%s.dat", configuration, etachars[i], etachars[i+1]);
        if (flag == 3) sprintf(configfilename, "datafiles/%s/data_longDCA_res_sim_eta%s-%s.dat", configuration, etachars[i], etachars[i+1]);
        if (flag == 4) sprintf(configfilename, "datafiles/%s/data_transDCA_res_sim_eta%s-%s.dat", configuration, etachars[i], etachars[i+1]);
        ofstream outdata;
        outdata.open(configfilename); // opens the file
        if( !outdata ) { // file couldn't be opened
            cerr << "Error: file could not be opened" << endl;
            exit(1);
        }

        for (int j=0; j<n; j++)
            outdata << p[j] << " " << st_dev[i][j] << endl;
        outdata.close();

    }
    


    /*
    c15->Print("plots/mom_res_SD_by_p.pdf");



    //fill graphs and plot, also get PWG requirements
    TCanvas *c16 = new TCanvas("c16","c16",1200,900);
    c16->Divide(5,2);
    
    for(int i=0; i<10; i++){

        //save each column of st_dev 2d array into an array
        Double_t st_dev_col[canvas_ctr];
        for (int a=0; a<14; a++){
            // cout << "Otirignal " << st_dev[a][i] <<", " << endl;
            if (a==0 || a ==13) continue;
            st_dev_col[a-1] = st_dev[a][i];
            // cout << st_dev[a][i] << ", " << endl;;
        }
        // cout << endl;

        momPlots_by_eta[i] = new TGraph(canvas_ctr,etaforplot,st_dev_col);

        pretty_TGraph(momPlots_by_eta[i],62,20,"Track #eta [rad]","#Deltap/p (%)","");
        momPlots_by_eta[i]->GetXaxis() -> SetTitleOffset(0.0);
        
        // for (int b=0; b<14; b++){
        //     cout << b << ": (" << etaforplot[b] << ", " << st_dev_col[b] << ")" << endl;
        // }
        if (i==6) {
            for (int b=0; b<canvas_ctr; b++){
                cout << b << ": (" << etaforplot[b] << ", " << st_dev_col[b] << ")" << endl;
            }
        }


        c16->cd(i+1);

        gPad -> SetLeftMargin(0.2);

        momPlots_by_eta[i]->Draw("ALP"); //"ACP"

        //momPlots_by_eta[i]->SetLineColor(4);
        //momPlots_by_eta[i]->SetMarkerColor(2);
        //momPlots_by_eta[i]->SetMarkerStyle(2);

        momPlots_by_eta[i]->GetXaxis()->SetTitle("Track #eta [rad]");
        momPlots_by_eta[i]->GetYaxis()->SetTitle("#Deltap/p");
        momPlots_by_eta[i]->GetXaxis()->CenterTitle();
        momPlots_by_eta[i]->GetYaxis()->CenterTitle();
	
   	 
        plabels[i]->SetTextFont(43); 
	    plabels[i]->SetTextSize(18);
        plabels[i]->Draw("same");
  
    }

    c16->Print("plots/mom_res_SD_by_eta.pdf");
    */







}

void pretty_TGraph( TGraph * g, int color = 1, int marker = 20, TString xtitle = "", TString ytitle = "", TString title=""){
    g -> SetLineColor(color);
    g -> SetMarkerColor(color);
    g -> SetLineWidth(2);
    g -> SetMarkerStyle(marker);

    g -> GetXaxis() -> SetTitle(xtitle);
    g -> GetXaxis() -> CenterTitle();
    g -> GetXaxis() -> SetNdivisions(107);
    g -> GetXaxis() -> SetTitleSize(0.05);
    g -> GetXaxis() -> SetLabelSize(0.05);

    g -> GetYaxis() -> SetTitle(ytitle);
    g -> GetYaxis() -> CenterTitle();
    g -> GetYaxis() -> SetNdivisions(107);
    g -> GetYaxis() -> SetTitleSize(0.05);
    g -> GetYaxis() -> SetLabelSize(0.05);

    g -> SetTitle(title);
}


//plot histogram of dp/p in multiple bins of eta and p (see graph for specifics)
//eta: [-3.5,3.5] in bins of 0.5
//p: [0,10] in bins of 1 or 2? //try 1 first
void plotMomRes(int nEntries){    

    int ctr[14];
    int max[14];
    int savei[14];
    int maxphi[14];
    int maxtheta[14];
    int maxlongDCA[14];
    int maxtransDCA[14];


    //fill histograms based on the bins
    for (int i=0; i<nEntries; i++) {
        if (i%1000000==0) cout << "jet: " << i << " out of: " << nEntries << endl;
        if (i==1000000) break;
        
        tree->GetEntry(i);
        float eta = calculateEta(px,py,pz);
        float p = calculateP(px,py,pz);
        float momRes = calculateMomRes(gpx,gpy,gpz,px,py,pz);
        float phiRes = calculatePhiRes(gpx,gpy,gpz,px,py,pz);
        float thetaRes = calculateThetaRes(gpx,gpy,gpz,px,py,pz);
        float longDCARes = calculateLongDCA(pcaz,gvz); //a DCA
        float transDCARes = dca2d; // r-phi DCA
        

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
                    if (phiRes > maxphi[i]) maxphi[i] = phiRes;
                    if (thetaRes > maxtheta[i]) maxtheta[i] = thetaRes;
                    if (longDCARes > maxlongDCA[i]) maxlongDCA[i] = longDCARes;
                    if (transDCARes > maxtransDCA[i]) maxtransDCA[i] = transDCARes;


                    if (p > j*2 && p <= j*2+2) {
                        mom_res_preliminary_hists[i][j]->Fill(momRes); //does comparison have to be a double?
                        phi_res_preliminary_hists[i][j]->Fill(phiRes);
                        theta_res_preliminary_hists[i][j]->Fill(thetaRes);
                        longDCA_res_preliminary_hists[i][j]->Fill(longDCARes);
                        transDCA_res_preliminary_hists[i][j]->Fill(transDCARes);
                    }
                }
            }
        }
        


    }

    //if the first p bin in each eta range has less than 100 events, exclude that eta range
    int canvas_ctr = 0;
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
        cout << "max phi " << i << ": " << maxphi[i] << endl;
        cout << "max theta " << i << ": " << maxtheta[i] << endl;
        cout << "max longDCA " << i << ": " << maxlongDCA[i] << endl;
        cout << "max transDCA " << i << ": " << maxtransDCA[i] << endl;
    }



    //calculate standard deviations
    Double_t** st_dev_direct = 0;
    Double_t** phi_st_dev_direct = 0;
    Double_t** theta_st_dev_direct = 0;
    Double_t** longDCA_st_dev_direct = 0;
    Double_t** transDCA_st_dev_direct = 0;
    st_dev_direct = new double*[14]; //[eta range][p range] in order from min to max (neg->pos)
    phi_st_dev_direct = new double*[14];
    theta_st_dev_direct = new double*[14];
    longDCA_st_dev_direct = new double*[14];
    transDCA_st_dev_direct = new double*[14];
    for (int i=0; i<14; i++){
        st_dev_direct[i] = new double[10];
        phi_st_dev_direct[i] = new double[10];
        theta_st_dev_direct[i] = new double[10];
        longDCA_st_dev_direct[i] = new double[10];
        transDCA_st_dev_direct[i] = new double[10];
        for (int j=0; j<10; j++){
            st_dev_direct[i][j] = mom_res_preliminary_hists[i][j]->GetStdDev();
            phi_st_dev_direct[i][j] = phi_res_preliminary_hists[i][j]->GetStdDev();
            theta_st_dev_direct[i][j] = theta_res_preliminary_hists[i][j]->GetStdDev();
            longDCA_st_dev_direct[i][j] = longDCA_res_preliminary_hists[i][j]->GetStdDev();
            transDCA_st_dev_direct[i][j] = transDCA_res_preliminary_hists[i][j]->GetStdDev();
        }
    }

    cout << "Preliminary histograms made" << endl;


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

            char phititle[1024];
            sprintf(phititle, "phi_res_eta%s-%s_p%d-%d_direct_rebinned", etachars[i], etachars[i+1], j*2, j*2+2);
            char philabel[1024];
            sprintf(philabel, ";d#phi, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            phi_res_hists[i][j] = new TH1F(phititle, philabel, 100, -3*phi_st_dev_direct[i][j], 3*phi_st_dev_direct[i][j]);
            char thetatitle[1024];
            sprintf(thetatitle, "theta_res_eta%s-%s_p%d-%d_direct_rebinned", etachars[i], etachars[i+1], j*2, j*2+2);
            char thetalabel[1024];
            sprintf(thetalabel, ";d#theta, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            theta_res_hists[i][j] = new TH1F(thetatitle, thetalabel, 100, -3*theta_st_dev_direct[i][j], 3*phi_st_dev_direct[i][j]);
            char longDCAtitle[1024];
            sprintf(longDCAtitle, "longDCA_res_eta%s-%s_p%d-%d_direct_rebinned", etachars[i], etachars[i+1], j*2, j*2+2);
            char longDCAlabel[1024];
            sprintf(longDCAlabel, ";z DCA, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            longDCA_res_hists[i][j] = new TH1F(longDCAtitle, longDCAlabel, 100, -3*longDCA_st_dev_direct[i][j], 3*phi_st_dev_direct[i][j]);
            char transDCAtitle[1024];
            sprintf(transDCAtitle, "transDCA_res_eta%s-%s_p%d-%d_direct_rebinned", etachars[i], etachars[i+1], j*2, j*2+2);
            char transDCAlabel[1024];
            sprintf(transDCAlabel, ";r-#phi DCA, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            transDCA_res_hists[i][j] = new TH1F(transDCAtitle, transDCAlabel, 100, -3*transDCA_st_dev_direct[i][j], 3*phi_st_dev_direct[i][j]);

        }
    }


    //fill histograms with new binning from std dev w fits
    for (int i=0; i<nEntries; i++) {
        if (i%1000000==0) cout << "jet: " << i << " out of: " << nEntries << endl;
        if (i==1000000) break;
        
        tree->GetEntry(i);
        float eta = calculateEta(px,py,pz);
        float p = calculateP(px,py,pz);
        float momRes = calculateMomRes(gpx,gpy,gpz,px,py,pz);
        float phiRes = calculatePhiRes(gpx,gpy,gpz,px,py,pz);
        float thetaRes = calculateThetaRes(gpx,gpy,gpz,px,py,pz);
        float longDCARes = calculateLongDCA(pcaz,gvz); //a DCA
        float transDCARes = dca2d; // r-phi DCA

        for (int i=0; i<14; i++){
            for (int j=0; j<10; j++){
                if (eta > etavals[i] && eta <= etavals[i+1]) {
                    if (p > j*2 && p <= j*2+2) {
                        mom_res_hists[i][j]->Fill(momRes); //does comparison have to be a double?
                        mom_res_hists_direct[i][j]->Fill(momRes);
                        phi_res_hists[i][j]->Fill(phiRes);
                        theta_res_hists[i][j]->Fill(thetaRes);
                        longDCA_res_hists[i][j]->Fill(longDCARes);
                        transDCA_res_hists[i][j]->Fill(transDCARes);
                    }
                }
            }
        }
    }

    Double_t** st_dev_final_fromfit = makeFits(mom_res_hists, st_dev_fromfit);
    Double_t** st_dev_final_direct = makeFits(mom_res_hists, st_dev_direct);

    Double_t** phi_st_dev_final = makeFits(phi_res_hists, phi_st_dev_direct);
    Double_t** theta_st_dev_final = makeFits(theta_res_hists, theta_st_dev_direct);
    Double_t** longDCA_st_dev_final = makeFits(longDCA_res_hists, longDCA_st_dev_direct);
    Double_t** transDCA_st_dev_final = makeFits(transDCA_res_hists, transDCA_st_dev_direct);


    cout << "About to draw histograms to canvas" << endl;

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
        sprintf(outname, "plots/%s/mom_res_eta%s-%s_bins100_rebinned.pdf", configuration, etachars[i], etachars[i+1]);
        canvas[i]->Print(outname);
    }


    //make other canvases to draw
    vector<TCanvas*> phicanvas(14);
    vector<TCanvas*> thetacanvas(14);
    vector<TCanvas*> longDCAcanvas(14);
    vector<TCanvas*> transDCAcanvas(14);
    for (int i=0; i<14; i++){
        char phi_canvas_title[1024]; sprintf(phi_canvas_title, "c%d phi", i);
        phicanvas[i] = new TCanvas(phi_canvas_title,phi_canvas_title,1200,900);
        phicanvas[i]->Divide(5,2);
        char theta_canvas_title[1024]; sprintf(theta_canvas_title, "c%d theta", i);
        thetacanvas[i] = new TCanvas(theta_canvas_title,theta_canvas_title,1200,900);
        thetacanvas[i]->Divide(5,2);
        char longDCA_canvas_title[1024]; sprintf(longDCA_canvas_title, "c%d longDCA", i);
        longDCAcanvas[i] = new TCanvas(longDCA_canvas_title,longDCA_canvas_title,1200,900);
        longDCAcanvas[i]->Divide(5,2);
        char transDCA_canvas_title[1024]; sprintf(transDCA_canvas_title, "c%d transDCA", i);
        transDCAcanvas[i] = new TCanvas(transDCA_canvas_title,transDCA_canvas_title,1200,900);
        transDCAcanvas[i]->Divide(5,2);

        for (int j=0; j<10; j++){
            phicanvas[i]->cd(j+1);
            phi_res_hists[i][j]->Draw();
            thetacanvas[i]->cd(j+1);
            theta_res_hists[i][j]->Draw();
            longDCAcanvas[i]->cd(j+1);
            longDCA_res_hists[i][j]->Draw();
            transDCAcanvas[i]->cd(j+1);
            transDCA_res_hists[i][j]->Draw();

        }
        char phioutname[1024];
        sprintf(phioutname, "plots/%s/phi_res_eta%s-%s_bins100_rebinned.pdf", configuration, etachars[i], etachars[i+1]);
        phicanvas[i]->Print(phioutname);
        char thetaoutname[1024];
        sprintf(thetaoutname, "plots/%s/theta_res_eta%s-%s_bins100_rebinned.pdf", configuration, etachars[i], etachars[i+1]);
        thetacanvas[i]->Print(thetaoutname);
        char longDCAoutname[1024];
        sprintf(longDCAoutname, "plots/%s/longDCA_res_eta%s-%s_bins100_rebinned.pdf", configuration, etachars[i], etachars[i+1]);
        longDCAcanvas[i]->Print(longDCAoutname);
        char transDCAoutname[1024];
        sprintf(transDCAoutname, "plots/%s/transDCA_res_eta%s-%s_bins100_rebinned.pdf", configuration, etachars[i], etachars[i+1]);
        transDCAcanvas[i]->Print(transDCAoutname);

    }


    // for (int i=0; i<14; i++){
    //     for (int j=0; j<10; j++){
    //         canvas[i]->cd(j+1);
    //         // mom_res_preliminary_hists[i][j]->Draw();
    //         mom_res_hists[i][j]->Draw();

    //         // canvas_direct[i]->cd(j+q);
    //         // mom_res_hists_direct[i]->Draw();

    //     }
    //     char outname[1024];
    //     sprintf(outname, "plots/%s/mom_res_etan%s-%s_bins100_rebinned.pdf", configuration, etachars[i], etachars[i+1]);
    //     canvas[i]->Print(outname);
    // }
    


    //plot SD
    // plotSD(st_dev_final_fromfit, st_dev_final_direct, canvas_ctr);
    plotSD(st_dev_final_direct, canvas_ctr);
    plotSD(phi_st_dev_final, canvas_ctr, 1);
    plotSD(theta_st_dev_final, canvas_ctr, 2);
    plotSD(longDCA_st_dev_final, canvas_ctr, 3);
    plotSD(transDCA_st_dev_final, canvas_ctr, 4);

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
    tree->SetBranchAddress("dca2d", &dca2d);
    tree->SetBranchAddress("pcaz", &pcaz);
    tree->SetBranchAddress("gvz", &gvz);

    
    const int nEntries = tree->GetEntries();
    plotMomRes(nEntries); 

}
