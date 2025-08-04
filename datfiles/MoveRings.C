void MoveRings()
{

    ifstream infile("S4Rings.txt");
    ofstream ofile("S4Rings_new.txt");
    
    int Ring[128][2]={0};
    int RingNew[128][2]={0};

    char dum[200];
    infile.getline(dum,200);
    int dat1,dat2,dat3;
    for(int i=0;i<128;i++)
      {
	infile >> dat1 >> dat2 >> dat3;
	Ring[dat3][0] = dat1;
	Ring[dat3][1] = dat2;
      }

    ofile << dum << endl;

    
    for(int i=0;i<128;i++)
      {
	int newring =0;
	if(i<32)
	  newring = 31-i;
	else if(i<64)
	  newring = 63-i+32;
	else if(i<96)
	  newring = 95-i+64;
	else
	  newring = 127-i+96;
	RingNew[newring][0] = Ring[i][0];
	RingNew[newring][1] = Ring[i][1];

	cout << i << " " << newring << endl;
      }
    for(int i=0;i<128;i++)
      {
	ofile << RingNew[i][0] << "\t"<< RingNew[i][1] <<"\t" << i << endl;
      }


    return;
}
