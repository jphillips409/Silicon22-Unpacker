void fbcal()
{
  string name;
  name = "fbline.dat";
  ifstream ifile(name.c_str());
  ofstream ofile("fb.cal");
  if (!ifile.is_open()) 
    {
      cout << "could not open file " << name << endl;
    }
  for(int itele =0;itele<14;itele++)
    {
      int n =0;
      ifile >> n;  // number of points
      cout << n << endl;
      float *x = new float [n]; //Front
      float *y = new float [n]; //Back
      for (int i=0;i<n;i++) 
	ifile >> x[i] >> y[i];
      TGraph *mygraph = new TGraph(n,x,y);
      mygraph->Draw("AP");
      mygraph->SetMarkerStyle(20);
      TF1 *myfit = new TF1("myfit","pol2",0,75);
      mygraph->Fit(myfit);
      TF1 *myfit2 = new TF1("myfit2","pol1",0,75);
      myfit2->SetLineColor(kRed);
      mygraph->Fit(myfit2,"+");
      ofile << 0 << " " << itele << " " << myfit->GetParameter(1);
      ofile << " " << myfit->GetParameter(0) << " ";
      ofile << myfit->GetParameter(2) << endl;
      c1->Update();
      c1->Modified();

      TMarker * mark;
      mark=(TMarker*)c1->WaitPrimitive("TMarker");
      delete mark;


    }



  ifile.close();
  ofile.close();
  return;
}

