#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

#include "insert.h"





void test2(){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    std:vector<TGraph*> EICPlots(14);

    string dat[14];
    for (int i=0; i<14; i++){
       char thingy[1024];
       sprintf(thingy, "/project/projectdirs/alice/eyeats/out_ECCE/60610987/tracking_output/dat%i.csv", i);
       dat[i] = thingy;
    }

    cout << "CHECKPOINT 1!!!!!" << endl;
    for (int i=0; i<14; i++){
        ifstream src(dat[i]);
        string line;
        int counter = 0;

        Double_t* EIC_x = 0;
        EIC_x = new double[9];

        Double_t* EIC_y = 0;
        EIC_y = new double[9];

        // while(!src.eof()){
        while (std::getline(src, linestr)){
            istringstream iss(linestr);
            // parse the stringstream into your vectors - I wouldn't use getline here, but you could
        }
            src >> line;
            cout << line << " " << endl;
            if (counter%2 ==0){
                //cout << line.substr(0,line.find(",")) << " " << i << " " << counter << endl;
                EIC_x[i] = stod(line.substr(0,line.find(",")));
            } else {
                EIC_y[i] = stod(line);
            }
            counter++;
        for (int j=0; j<9; j++){
            cout << EIC_x[j] << " " << EIC_y[j] << endl;
        }
        EICPlots[i] = new TGraph(9, EIC_x, EIC_y);
    }

}


