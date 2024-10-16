#define run2analysis_cxx
#include "run2analysis.h"
#include <TMatrix.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>


enum TrackQuality
{
	undefQuality = -1,
	loose = 0,
	tight = 1,
	highPurity = 2,
	confirmed = 3,      // means found by more than one iteration
	goodIterative = 4,  // meaningless
	looseSetWithPV = 5,
	highPuritySetWithPV = 6,
	discarded = 7,  // because a better track found. kept in the collection for reference....
	qualitySize = 8
};



// Modification of the code in September 2021 to 
// align it with the use of xtalk inversion (only for cluster cleaning) and saturation 
//
void run2analysis::Loop()
{

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;


   double pi=acos(-1);
   double dEdxSF [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF


   TH2D* dEdX0stripVsP = new TH2D("dEdX0stripVsP", "dEdX:Momentum [GeV]", 250,0,50, 400, 0.,20.);
   TH2D* dEdX0stripVsP_lowp = new TH2D("dEdX0stripVsP_lowp", "dEdX:Momentum [GeV]", 200,0,5, 400, 0.,20.);
   TH2D* dEdX0stripVsP_charge = new TH2D("dEdX0stripVsP_charge", "dEdX:Charge*Momentum [GeV]", 400,-5,5, 400, 0.,20.);

   TH1D* eta_distribution = new TH1D("eta_distribution", "eta_distribution", 200, -5, 5);
   TH1D* eta_distribution_NoSelection = new TH1D("eta_distribution_NoSelection", "eta_distribution_NoSelection", 200, -5, 5);
   TH2D* mass_reconstruction = new TH2D("mass_reconstruction", "mass_reconstruction", 200, 0, 5, 100, 0, 2); 

   dEdX0stripVsP->Sumw2();
   dEdX0stripVsP_lowp->Sumw2();
   dEdX0stripVsP_charge->Sumw2();


   TString outputfilename="ROOT_histograms/EtaDistribution.root";
   TFile* OutputHisto = new TFile(outputfilename,"RECREATE");
   TTree *outputTree = new TTree("outputTree", "Output Tree");

   cout << "save file: " << outputfilename << endl;
   cout << " --- Run:  " << nentries << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%1000000 == 0 && jentry!=0) cout << jentry <<  " (" << (100.*jentry)/(1.*nentries) << " %)" <<endl;

		if (npv<1) continue;

		int ntracks1=0;
		for (int itr=0; itr<ntracks; itr++)
		{

			int index_of_the_track=itr;
			int presk= track_prescale[index_of_the_track];
			bool selection=true;
			if (track_pt[itr] < 0.5) selection = false;
			//if (track_pt[index_of_the_track]>10) selection=false;
			//if (abs(track_eta[index_of_the_track])>2.1) selection=false;

			if (track_nvalidhits[index_of_the_track]<8) selection=false;
			if (track_npixhits[index_of_the_track]<2) selection=false;
			if (track_validfraction[index_of_the_track]<0.8) selection=false;
				
			bool is_high_qual =  (track_qual[index_of_the_track] & (1 << TrackQuality::highPurity)) >> TrackQuality::highPurity ;
			if (!is_high_qual) selection=false;
			if (track_chi2[index_of_the_track]>5) selection=false;

			// which selection for dz and dyx w/r to the primary vertex ? here 1st PV to pass the selection
			if (abs(track_dxy[index_of_the_track])>0.5) selection=false;
			if (abs(track_dz[index_of_the_track])>0.5) selection=false;

			// no cut on the isolation
			// CUT ISO ICI 

			// cut on sigma pT/pT for signal  :
			if (track_pterr[index_of_the_track]/track_pt[index_of_the_track]>0.25) selection=false;

                        eta_distribution_NoSelection->Fill(track_eta[itr]);

			if (!selection) continue;

			std::vector <float> charge_corr;
			std::vector <float> pathlength;
			std::vector <int> subdetId;
			std::vector <UInt_t> detId;
			std::vector <int> moduleGeometry;
			std::vector <bool> bool_cleaning;
			std::vector <bool> mustBeInside;

			int nsatclu=0;
			for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++)
			{
				float ch1=dedx_charge[iclu];
				bool clean1=true;
				if (dedx_subdetid[iclu]>=3)
				{
					// strip
					// without any correction :
                    //ch1=sclus_charge[iclu];  // charge without any correction
                    // clean1=sclus_clusclean[iclu];  
					// with Saturation only (but no Xtalk inversion) :
					float check_charge=0;
					vector<int> Quncor;
					for (int istrip=sclus_index_strip[iclu]; istrip<sclus_index_strip[iclu]+sclus_nstrip[iclu]; istrip++)
					{
						check_charge+=strip_ampl[istrip];
						Quncor.push_back(strip_ampl[istrip]);
					}
					vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25);
					float newcharge =0;
					for (unsigned int inwc=0; inwc<Qcor.size(); inwc++) { newcharge+=Qcor[inwc]; }
					ch1=newcharge;
					clean1=sclus_clusclean2[iclu]; // clusterCleaning with Xtalk inversion and Saturation (September 22, 2021)
					// for the record : there is a possibility to recompute the cluster cleaning with the Saturation only : clusterCleaning(Qcor, 1, &exitCode)
					// but it is not what we want here
					if (dedx_isstrip[iclu] && clean1 && dedx_insideTkMod[iclu] && (sclus_sat254[iclu] || sclus_sat255[iclu])) nsatclu++;

					if (clean1 && dedx_insideTkMod[iclu])
					{
						// fill info
						charge_corr.push_back(ch1);
						pathlength.push_back(dedx_pathlength[iclu]);
						subdetId.push_back(dedx_subdetid[iclu]);
						detId.push_back(dedx_detid[iclu]);
						moduleGeometry.push_back(dedx_modulgeom[iclu]);
						mustBeInside.push_back(dedx_insideTkMod[iclu]);
						bool_cleaning.push_back(clean1);
					}  // end if cleaning & inside			     
				}  // end if Strip 
			}

						
			int nv=0;
			int ns=0;

        	double ih0_strip = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF,  NULL,2, 0.,  nv, ns);

			if(charge_corr.size() > 9 )
			{
				if (track_p[itr]<5) {
					dEdX0stripVsP_lowp->Fill(track_p[itr],ih0_strip,presk);
					dEdX0stripVsP_charge->Fill(track_p[itr]*track_charge[itr],ih0_strip,presk);

                                        eta_distribution->Fill(track_eta[itr]);
                                        if (ih0_strip < 3.05983) continue;
                                        else mass_reconstruction->Fill(track_p[itr], track_p[itr] * sqrt((ih0_strip-3.05983)/2.66455), presk);
				}
				dEdX0stripVsP->Fill(track_p[itr],ih0_strip,presk);
			}
      } //END LOOP ON ALL TRACKS

   }

   OutputHisto->cd();

   dEdX0stripVsP->Write();
   dEdX0stripVsP_lowp->Write();
   dEdX0stripVsP_charge->Write();
   eta_distribution->Write();
   eta_distribution_NoSelection->Write();
   mass_reconstruction->Write();
   OutputHisto->Close();

}


double run2analysis::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, int & nv, int & ns) {
  double result= getdEdX(charge, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, scaleFactors, templateHisto, n_estim, dropLowerDeDxValue, 0., nv, ns);
  return result;
}


double run2analysis::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, double dropHigherDeDxValue, int & nv, int & ns) {
  double result=-1;
     size_t MaxStripNOM=99;
     bool usePixel=true;
     bool useStrip=true;

     std::vector<double> vect;

     bool debugprint=false;
     unsigned int SiStripNOM = 0;
     ns=0;

     for(unsigned int h=0;h<charge.size();h++){
        if (debugprint) std::cout << "look on dedxHits in computedEdx " << h << std::endl;
        if(!usePixel && subdetId[h]<3)continue; // skip pixels
        if(!useStrip && subdetId[h]>=3)continue; // skip strips        
        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;

        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;
        if(useStrip && subdetId[h]>=3 && !mustBeInside[h])continue;
        if(useStrip && subdetId[h]>=3 && ++SiStripNOM > MaxStripNOM) continue; // skip remaining strips, but not pixel

        int ClusterCharge = charge[h];
        if (subdetId[h]>=3 && charge[h]>=254) ns++;

        double scaleFactor = scaleFactors[0];
        if (subdetId[h]<3) scaleFactor *= scaleFactors[1]; // add pixel scaling
        if (debugprint) std::cout << " after SF " << std::endl;

        if(templateHisto){  //save discriminator probability
           double ChargeOverPathlength = scaleFactor*ClusterCharge/(pathlength[h]*10.0*(subdetId[h]<3?265:1));
           int    BinX   = templateHisto->GetXaxis()->FindBin(moduleGeometry[h]);
           int    BinY   = templateHisto->GetYaxis()->FindBin(pathlength[h]*10.0); //*10 because of cm-->mm
           int    BinZ   = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
           double Prob   = templateHisto->GetBinContent(BinX,BinY,BinZ);
           vect.push_back(Prob); //save probability
           if (debugprint) std::cout << " after Prob vect.push_back " << std::endl;
        }else{
           double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
           double ChargeOverPathlength = scaleFactor*Norm*ClusterCharge/pathlength[h];
           vect.push_back(ChargeOverPathlength); //save charge
           if (debugprint) std::cout << " after ChargeOverPathlength vect.push_back " << std::endl;
        }
     }

     if(dropLowerDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
         int nTrunc = tmp.size()*dropLowerDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropLowerDeDxValue " << std::endl;

     if(dropHigherDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::less<double>() );
         int nTrunc = tmp.size()*dropHigherDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropHigherDeDxValue " << std::endl;






     int size = vect.size();
     nv = size;

     if(size>0){
        if(templateHisto){  //dEdx discriminator
          //Ias discriminator
          result = 1.0/(12*size);
           std::sort(vect.begin(), vect.end(), std::less<double>() );
           for(int i=1;i<=size;i++){
              result += vect[i-1] * pow(vect[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
           }
           result *= (3.0/size);
           if (debugprint) std::cout << " Ias discriminator " << result << std::endl;
        }else{  //dEdx estimator
           //harmonic2 estimator
           result=0;
//           double expo = -2;
           double expo = -1* n_estim;
           for(int i = 0; i< size; i ++){
              result+=pow(vect[i],expo);
           }
           result = pow(result/size,1./expo);
           if (debugprint) std::cout << " harmonic discriminator " << result << " with expo " << expo << std::endl;
        }
     }else{
        result = -1;
     }
     if (debugprint) std::cout << " ok finished computeDeDx " << std::endl;


  return result;
}


int run2analysis::GetLayerLabel(int subdetid_, UInt_t detid_, int year)
{
// from https://github.com/dapparu/HSCP/blob/c69cf1c71dd99289f72ab6d03077915776c85690/src/Cluster.cc
// and https://cmssdt.cern.ch/lxr/source/DataFormats/SiPixelDetId/interface/PXBDetId.h
// https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder
        if(subdetid_==1)
        {
             if (year==2016) {
                if(((detid_>>16)&0xF)==1) return 23;
                else if(((detid_>>16)&0xF)==2) return 24;
                else if(((detid_>>16)&0xF)==3) return 25;
                else if(((detid_>>16)&0xF)==4) return 26;  // do not exist in 2016
             }
             else {
                if(((detid_>>20)&0xF)==1) return 23;
                else if(((detid_>>20)&0xF)==2) return 24;
                else if(((detid_>>20)&0xF)==3) return 25;
                else if(((detid_>>20)&0xF)==4) return 26;
             }

        }
        else if(subdetid_==2)
        {
             if (year==2016) {
                if(((detid_>>16)&0xF)==1) return 27;
                else if(((detid_>>16)&0xF)==2) return 28;
                else if(((detid_>>16)&0xF)==3) return 29; // do not exist in 2016
             }
             else {
                if(((detid_>>18)&0xF)==1) return 27;
                else if(((detid_>>18)&0xF)==2) return 28;
                else if(((detid_>>18)&0xF)==3) return 29;
             }
/*  https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder for 2016
                cout << "  side " << int((detid_>>23)&0x3) <<
                        "  disk " << int((detid_>>16)&0xF) << 
                        "  blade "  << int((detid_>>10)&0x3F) <<
                        "  panel "  << int((detid_>>8)&0x3) <<
                        "  mod "  << int((detid_>>2)&0x3F) << endl;
*/
//                if (((detid_>>16)&0xF)==0) cout << " disk 0 ? " << std::endl;
        }
        else if(subdetid_==3)  // TIB
        {
                if(((detid_>>14)&0x7)==1) return 1;
                else if(((detid_>>14)&0x7)==2) return 2;
                else if(((detid_>>14)&0x7)==3) return 3;
                else if(((detid_>>14)&0x7)==4) return 4;
        }
        else if(subdetid_==5) // TOB
        {
                if(((detid_>>14)&0x7)==1) return 5;
                else if(((detid_>>14)&0x7)==2) return 6;
                else if(((detid_>>14)&0x7)==3) return 7;
                else if(((detid_>>14)&0x7)==4) return 8;
                else if(((detid_>>14)&0x7)==5) return 9;
                else if(((detid_>>14)&0x7)==6) return 10;
        }
        else if(subdetid_==4)  //TID
        {
                if(((detid_>>11)&0x3)==1) return 11;
                else if(((detid_>>11)&0x3)==2) return 12;
                else if(((detid_>>11)&0x3)==3) return 13;
        }
        else if(subdetid_==6) // TEC
        {
                if(((detid_>>14)&0xF)==1) return 14;
                else if(((detid_>>14)&0xF)==2) return 15;
                else if(((detid_>>14)&0xF)==3) return 16;
                else if(((detid_>>14)&0xF)==4) return 17;
                else if(((detid_>>14)&0xF)==5) return 18;
                else if(((detid_>>14)&0xF)==6) return 19;
                else if(((detid_>>14)&0xF)==7) return 20;
                else if(((detid_>>14)&0xF)==8) return 21;
                else if(((detid_>>14)&0xF)==9) return 22;
        }
        return -1;
}


std::vector<int> run2analysis::SaturationCorrection(const std::vector<int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat) {
  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);
  Double_t a=1-2*x1-2*x2;
//  TMatrix A(N,N);

//---  que pour 1 max bien net
 if(Q.size()<2 || Q.size()>8){
        for (unsigned int i=0;i<Q.size();i++){
                QII.push_back((int) Q[i]);
        }
        return QII;
  }
 if(way){
          vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())      ;
          if(*mQ>253){
                 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
                 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
                     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
          }
      else{
          return Q; // no saturation --> no x-talk inversion
      }
  }
//---
 // do nothing else
 return Q; 
}


bool run2analysis::clusterCleaning(const std::vector<int>&  Q, int crosstalkInv, uint8_t * exitCode)
{
     vector<int>  ampls = Q;
//   if(crosstalkInv==1)ampls = CrossTalkInv(ampls,0.10,0.04, true,20,25);


  // ----------------  COMPTAGE DU NOMBRE DE MAXIMA   --------------------------
  //----------------------------------------------------------------------------
  //
         Int_t NofMax=0; Int_t recur255=1; Int_t recur254=1;
         bool MaxOnStart=false;bool MaxInMiddle=false, MaxOnEnd =false;
         Int_t MaxPos=0;
        // D?but avec max
        if(ampls.size()!=1 && ((ampls[0]>ampls[1])
            || (ampls.size()>2 && ampls[0]==ampls[1] && ampls[1]>ampls[2] && ampls[0]!=254 && ampls[0]!=255)
            || (ampls.size()==2 && ampls[0]==ampls[1] && ampls[0]!=254 && ampls[0]!=255)) ){
          NofMax=NofMax+1;  MaxOnStart=true;  }

        // Maximum entour?
        if(ampls.size()>2){
          for (unsigned int i =1; i < ampls.size()-1; i++) {
                if( (ampls[i]>ampls[i-1] && ampls[i]>ampls[i+1])
                    || (ampls.size()>3 && i>0 && i<ampls.size()-2 && ampls[i]==ampls[i+1] && ampls[i]>ampls[i-1] && ampls[i]>ampls[i+2] && ampls[i]!=254 && ampls[i]!=255) ){
                 NofMax=NofMax+1; MaxInMiddle=true;  MaxPos=i;
                }
                if(ampls[i]==255 && ampls[i]==ampls[i-1]) {
                        recur255=recur255+1;
                        MaxPos=i-(recur255/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
                if(ampls[i]==254 && ampls[i]==ampls[i-1]) {
                        recur254=recur254+1;
                        MaxPos=i-(recur254/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
            }
        }
        // Fin avec un max
       if(ampls.size()>1){
          if(ampls[ampls.size()-1]>ampls[ampls.size()-2]
             || (ampls.size()>2 && ampls[ampls.size()-1]==ampls[ampls.size()-2] && ampls[ampls.size()-2]>ampls[ampls.size()-3] )
             ||  ampls[ampls.size()-1]==255){
           NofMax=NofMax+1;  MaxOnEnd=true;   }
         }
        // Si une seule strip touch?e
        if(ampls.size()==1){    NofMax=1;}

  // ---  SELECTION EN FONCTION DE LA FORME POUR LES MAXIMA UNIQUES ---------
  //------------------------------------------------------------------------
  //
  //               ____
  //              |    |____
  //          ____|    |    |
  //         |    |    |    |____
  //     ____|    |    |    |    |
  //    |    |    |    |    |    |____
  //  __|____|____|____|____|____|____|__
  //    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
  //

   bool shapecdtn=false;
   if (exitCode) *exitCode = 255;

      if(crosstalkInv==1){
        if(NofMax==1){shapecdtn=true; if (exitCode) *exitCode=0;}
        return shapecdtn;
      }

        Float_t C_M=0.0;        Float_t C_D=0.0;        Float_t C_Mn=10000;     Float_t C_Dn=10000;     Float_t C_Mnn=10000;    Float_t C_Dnn=10000;
        Int_t CDPos;
        Float_t coeff1=1.7;     Float_t coeff2=2.0;
        Float_t coeffn=0.10;    Float_t coeffnn=0.02; Float_t noise=4.0;

        if(NofMax==1){

                if(MaxOnStart==true){
                        C_M=(Float_t)ampls[0]; C_D=(Float_t)ampls[1];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[2] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=2;}
                                else if(ampls.size()>3){ C_Dn=(Float_t)ampls[2];  C_Dnn=(Float_t)ampls[3] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=3;
                                }
                }

                if(MaxOnEnd==true){
                        C_M=(Float_t)ampls[ampls.size()-1]; C_D=(Float_t)ampls[ampls.size()-2];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[0] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=4;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[ampls.size()-3] ; C_Dnn=(Float_t)ampls[ampls.size()-4] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=5;
                                }
                }
                if(MaxInMiddle==true){
                        C_M=(Float_t)ampls[MaxPos];
                        int LeftOfMaxPos=MaxPos-1;if(LeftOfMaxPos<=0)LeftOfMaxPos=0;
                        int RightOfMaxPos=MaxPos+1;if(RightOfMaxPos>=(int)ampls.size())RightOfMaxPos=ampls.size()-1;
                        if(ampls[LeftOfMaxPos]<ampls[RightOfMaxPos]){ C_D=(Float_t)ampls[RightOfMaxPos]; C_Mn=(Float_t)ampls[LeftOfMaxPos];CDPos=RightOfMaxPos;} else{ C_D=(Float_t)ampls[LeftOfMaxPos]; C_Mn=(Float_t)ampls[RightOfMaxPos];CDPos=LeftOfMaxPos;}
                        if(C_Mn<coeff1*coeffn*C_M+coeff2*coeffnn*C_D+2*noise || C_M==255){
                                if(ampls.size()==3) shapecdtn=true ;
                                else if(ampls.size()>3){
                                        if(CDPos>MaxPos){
                                                if(ampls.size()-CDPos-1==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1==1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1>1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=(Float_t)ampls[CDPos+2];
                                                }
                                                if(MaxPos>=2){
                                                        C_Mnn=(Float_t)ampls[MaxPos-2];
                                                }
                                                else if(MaxPos<2) C_Mnn=0;
                                        }
                                        if(CDPos<MaxPos){
                                                if(CDPos==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(CDPos==1){
                                                        C_Dn=(Float_t)ampls[0];
                                                        C_Dnn=0;
                                                }
                                                if(CDPos>1){
                                                        C_Dn=(Float_t)ampls[CDPos-1];
                                                        C_Dnn=(Float_t)ampls[CDPos-2];
                                                }
                                                if(ampls.size()-LeftOfMaxPos>1 && MaxPos+2<(int)(ampls.size())-1){
                                                        C_Mnn=(Float_t)ampls[MaxPos+2];
                                                }else C_Mnn=0;
                                        }
                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                           && C_Mnn<=coeff1*coeffn*C_Mn+coeff2*coeffnn*C_M+2*noise
                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise) {
                                                shapecdtn=true;
                                        }

                                }
                        } else if (exitCode) *exitCode=6;
                }
        }
        else if (NofMax>1 && exitCode) *exitCode = 1; // more than one maximum
        if(ampls.size()==1){shapecdtn=true;}
        if(shapecdtn && exitCode) *exitCode=0;

   return shapecdtn;
}




