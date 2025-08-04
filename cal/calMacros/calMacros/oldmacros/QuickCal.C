{
  ifstream ifile("ceasar_Na22.dat");
  ifstream ifile2("ceasar_bkg.dat");
  ofstream ofile("ceasar.peaks");
  int i1 = 0;
  double ipos1[196] ={0.};
  double ipos2[196] ={0.};
  double ipos3[196] ={0.};
  double idum = 0.;
  double idum2 = 0.;


  for(;;)
    {
      if(ifile.eof()) break;
      ifile >> i1 >> idum;
      if(idum == 0) idum = -1;
      ipos1[i1] = idum;     

    }

  for(;;)
    {
      if(ifile2.eof()) break;
      ifile2 >> i1 >> idum >> idum2;
      ipos2[i1] = idum;     
      ipos3[i1] = idum2;     
    }

  for(int itele = 0; itele < 1;itele++)
    {
      for(int istrip = 0; istrip<196;istrip++)
	{
	  ofile << istrip << " " << ipos1[istrip] << " " << ipos2[istrip];
	    ofile << " " << ipos3[istrip] << endl;

	}
    }

  ifile.close();
  ofile.close();

}
