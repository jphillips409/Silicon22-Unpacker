void Calibrate()
{


  TCanvas *mycan1 = (TCanvas*)gROOT->FindObject("mycan1");
  if(!mycan1)
    {
      mycan1 = new TCanvas("mycan1","",1200,1000);
      mycan1->Divide(4,4,0.0000001,0.0000001);
    }

  float energy[2][32][5] = {0.};
  float peak[2][32][5] = {0.};



  ifstream file("SiCalibPoints.txt");
  if(!file.is_open())
    {
      cout << "No Si Calib Data" << endl;
      return;
    }

  ofstream out("testSiRecalib.cal");

  int  itele = -1;
  int istrip = -1;

  for(int i = 0;i<64;i++)
    {

      file >> itele >> istrip;
      itele = itele-6;
      file >>energy[itele][istrip][0];
      file >> energy[itele][istrip][1] >> energy[itele][istrip][2];
      file >> energy[itele][istrip][3] >> energy[itele][istrip][4];
      file >> peak[itele][istrip][0] >> peak[itele][istrip][1];
      file >> peak[itele][istrip][2] >> peak[itele][istrip][3];
      file >> peak[itele][istrip][4];

    }

  float slope = 0.;
  float inter = 0.;
  for(int i = 0;i<16;i++)
    {
      mycan1->cd(i+1);
      TGraph *mygraph = new TGraph(5,peak[0][i],energy[0][i]);
      mygraph->SetMarkerStyle(20);
      mygraph->Draw("AP");

      TF1 *fit = new TF1("fit","pol1",0,500);
      mygraph->Fit("fit","RQ");

      slope = fit->GetParameter(1);
      inter = fit->GetParameter(0);
      cout << istrip <<" Slope = " << slope;
      cout << " inter = " << inter << endl;
      out << 0 << " " << i << " " << slope << " " << inter << endl;

    }

  TCanvas *mycan2 = (TCanvas*)gROOT->FindObject("mycan2");
  if(!mycan2)
    {
      mycan2 = new TCanvas("mycan2","",1200,1000);
      mycan2->Divide(4,4,0.0000001,0.0000001);
    }
  for(int i = 0;i<16;i++)
    {
      int istrip = i+16;
      mycan2->cd(i+1);
      TGraph *mygraph = new TGraph(5,peak[0][istrip],energy[0][istrip]);
      mygraph->SetTitle(Form("Strip %i",istrip));
      mygraph->SetMarkerStyle(20);
      mygraph->Draw("AP");

      TF1 *fit = new TF1("fit","pol1",0,500);
      mygraph->Fit("fit","RQ");

      cout << istrip << " Slope = " << fit->GetParameter(1);
      cout << " inter = " << fit->GetParameter(0) << endl;
      slope = fit->GetParameter(1);
      inter = fit->GetParameter(0);
      out << 0 << " " << istrip << " " << slope << " " << inter << endl;

    }
  TCanvas *mycan3 = (TCanvas*)gROOT->FindObject("mycan3");
  if(!mycan3)
    {
      mycan3 = new TCanvas("mycan3","",1200,1000);
      mycan3->Divide(4,4,0.0000001,0.0000001);
    }
  for(int i = 0;i<16;i++)
    {
      int istrip = i;
      mycan3->cd(i+1);
      TGraph *mygraph = new TGraph(5,peak[1][istrip],energy[1][istrip]);
      mygraph->SetTitle(Form("Strip %i",istrip));
      mygraph->SetMarkerStyle(20);
      mygraph->Draw("AP");

      TF1 *fit = new TF1("fit","pol1",0,500);
      mygraph->Fit("fit","RQ");

      cout << istrip << " Slope = " << fit->GetParameter(1);
      cout << " inter = " << fit->GetParameter(0) << endl;

      slope = fit->GetParameter(1);
      inter = fit->GetParameter(0);
      out << 1 << " " << istrip << " " << slope << " " << inter << endl;


    }
  TCanvas *mycan4 = (TCanvas*)gROOT->FindObject("mycan4");
  if(!mycan4)
    {
      mycan4 = new TCanvas("mycan4","",1200,1000);
      mycan4->Divide(4,4,0.0000001,0.0000001);
    }
  for(int i = 0;i<16;i++)
    {
      int istrip = i+16;
      mycan4->cd(i+1);
      TGraph *mygraph = new TGraph(5,peak[1][istrip],energy[1][istrip]);
      mygraph->SetTitle(Form("Strip %i",istrip));
      mygraph->SetMarkerStyle(20);
      mygraph->Draw("AP");

      TF1 *fit = new TF1("fit","pol1",0,500);
      mygraph->Fit("fit","RQ");

      cout << istrip << " Slope = " << fit->GetParameter(1);
      cout << " inter = " << fit->GetParameter(0) << endl;

      slope = fit->GetParameter(1);
      inter = fit->GetParameter(0);
      out << 1 << " " << istrip << " " << slope << " " << inter << endl;


    }


  return;
}
