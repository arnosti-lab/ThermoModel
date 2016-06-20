
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
using namespace std;

int countnew;

void runMAST(int, int, double);
double** getMASTresults(int);
int** findOverlaps(double const * const * const, int);
double ** removeOverlaps(double **, int, int**);
void orderMAST(int, double **, int);
int** findFootprints(double const* const* const,int);
void createInput(double const* const* const,int,int);
void create_thermo_inputs(int, int, int, int, int, int, int, int);
void help_flags();

void help_flags()
{
	cout << "The flags used in ./create_inputs are the following:" << endl;
	cout << "-thresh: When the SAME threshold is used for MAST analysis on Dl, Twi, and Sna, this is the threshold used" << endl;
	cout << "-Dthresh: When DIFFERENT thresholds are used for MAST analysis on Dl, Twi, and Sna, this is the Dl threshold used"<< endl;
	cout << "-Tthresh: When DIFFERENT thresholds are used for MAST analysis on Dl, Twi, and Sna, this is the Twi threshold used"<< endl;
	cout << "-Sthresh: When DIFFERENT thresholds are used for MAST analysis on Dl, Twi, and Sna, this is the Sn threshold used"<< endl;
	cout << "-dlPWM: Dl PWM used for the MAST analysis (1 = FlyReg, 2 = Wolfe, 3 = Average)" << endl;
	cout << "-twPWM: Twi PWM used for the MAST analysis (1 = Zinzen, 2 = BDTNP, 3 = Average)" << endl;
	cout << "-snPWM: Sna PWM used for the MAST analysis (1 = BDTNP, 2 = Zinzen selex, 3 = FlyReg)" << endl;
	cout << "-inter: Type of Interactions (none=0,NN=1,allN=2)" << endl;
	cout << "-coop: Type of Cooperativity (fixed=0,binned=1,ProteinBinned=2,Log=3,Linear=4,Gauss=5)" << endl;
	cout << "NOTE: If -coop is set to 1 or 2, need to include 1 more flag, -bincoop" << endl;
	cout << "-bincoop: When -coop is set to 1 or 2 --> 2 values: number of bins and size of bins" << endl;
	cout << "-quench: Type of Quenching (fixed=0,binned=1,MSB=2,Log=3,Linear=4,Gauss=5)" << endl;
	cout << "NOTE: If -quench is set to 1, need to include 1 more flag, -binquench" << endl;
	cout << "-binquench: When -quench is set to 1 --> 2 values: number of bins and size of bins" << endl;
	
}

void runMAST(int construct, int protein, double threshold,int* PWMs)
{
	//if protein = 1 : get MAST results for Dorsal sites on rho
	//if protein = 2 : get MAST results for Twist sites on rho
	//if protein = 3 : get MAST results for Snail sites on rho
	
	int i;
	char command[100], enhancer[100];
	char buffer[20];
	char c[10];
	int n;
	n = sprintf(c, "%i", construct);
	
	
	//rho
	
	if(construct<10)
		strcpy(enhancer,"con0");
	else 
		strcpy(enhancer,"con");
	strcat(enhancer,c);
	strcat(enhancer,".txt");
	

		
	if(protein==1)
	{
		if(PWMs[0]==1)
			strcpy(command,"./mast dl_flyreg.txt ");
		else if(PWMs[0]==2)
			strcpy(command,"./mast dl_B1H_wolfe.txt ");
		else if(PWMs[0]==3)
			strcpy(command,"./mast dl_average.txt ");
		else
		{
			cout << "CANNOT RUN MAST OUTPUT!!!! Incorrect Dl PWM passed to runMAST()" << endl;
			exit(1);
		}
		strcat(command,enhancer);
		strcat(command," -hit_list -mt ");
		gcvt(threshold,6,buffer);
		strcat(command,buffer);
		strcat(command," -bfile bg.txt > MASToutputDorsal.txt");
	}
	else if(protein==2)
	{
		if(PWMs[1]==1)
			strcpy(command,"./mast twi_selex_zinzen.txt ");
		else if(PWMs[1]==2)
			strcpy(command,"./mast twi_selex_BDTNP.txt ");
		else if(PWMs[1]==3)
			strcpy(command,"./mast twi_average.txt ");
		else
		{
			cout << "CANNOT RUN MAST OUTPUT!!!! Incorrect Twi PWM passed to runMAST()" << endl;
			exit(1);
		}
		strcat(command,enhancer);
		strcat(command," -hit_list -mt ");
		gcvt(threshold,6,buffer);
		strcat(command,buffer);
		strcat(command," -bfile bg.txt > MASToutputTwist.txt");
	}
	else if(protein==3)
	{
		if(PWMs[2]==1)
			strcpy(command,"./mast sna_BDTNP.txt ");
		else if(PWMs[2]==2)
			strcpy(command,"./mast sna_selex.txt ");
		else if(PWMs[2]==3)
			strcpy(command,"./mast sna_flyreg_rup.txt ");
		else
		{
			cout << "CANNOT RUN MAST OUTPUT!!!! Incorrect Sna PWM passed to runMAST()" << endl;
			exit(1);
		}
		strcat(command,enhancer);
		strcat(command," -hit_list -mt ");
		gcvt(threshold,6,buffer);
		strcat(command,buffer);
		strcat(command," -bfile bg.txt > MASToutputSnail.txt");
	}
	else
	{
		cout << "CANNOT READ MAST OUTPUT!!!! Incorrect protein value passed to runMAST()" << endl;
		exit(1);
	}
	
	
	
	
	
	if (!system(NULL)) {
		cerr << "System(NULL) failed, unable to run external command" << endl;
        exit (1);
	}
	
	i=system(command);
}

double** getMASTresults(int protein)
{
	//if protein = 1 : get MAST results for Dorsal sites on rho
	//if protein = 2 : get MAST results for Twist sites on rho
	//if protein = 3 : get MAST results for Snail sites on rho
	
	string dummy;
	int i;
	char c;
	char file[80];
	double** MASTresults = new double*[4];
	countnew++;
	
	if(protein==1)
		strcpy(file,"MASToutputDorsal.txt");
	else if(protein==2)
		strcpy(file,"MASToutputTwist.txt");
	else if(protein==3)
		strcpy(file,"MASToutputSnail.txt");
	else
	{
		cout << "CANNOT READ MAST OUTPUT!!!! Incorrect value passed to getMASTresults()" << endl;
		exit(1);
	}
	
	ifstream mast(file);
	int countBS = 0;
	while ( mast >> c ){
		//mast >> c;
		if(c!='#')
		{
			countBS ++;
		}
		getline(mast,dummy);
	}
	//cout << "found " << countBS << " binding sites " << endl;
	

	MASTresults[0] = new double[1];
	countnew++;
	MASTresults[0][0] = double(countBS);
	MASTresults[1] = new double[countBS];
	countnew++;
	MASTresults[2] = new double[countBS];
	countnew++;
	MASTresults[3] = new double[countBS];
	countnew++;
	
	countBS = 0;
	mast.close();
	mast.open(file);
	while ( mast >> c ){
		if(c!='#')
		{
			mast >> dummy;
			mast >> i;
			mast >> MASTresults[1][countBS];
			mast >> MASTresults[2][countBS];
			mast >> MASTresults[3][countBS];
			//cout << "start = " << start[countBS];
			//cout << ", end = " << finish[countBS];
			//cout << ", score = " << score[countBS] << endl;
			countBS++;
		}
		getline(mast,dummy);
	}
	return MASTresults;
}

int** findOverlaps(double const * const * const sites, int numsites)
{
	int i,j;
	int ** overlaps = new int*[numsites];
	countnew++;
	
	for(i=0;i<numsites;i++)
	{
		overlaps[i] = new int[numsites];
		countnew++;
	}
	for(i=0;i<numsites;i++)
	{
		for(j=0;j<numsites;j++)
		{
			if(i==j)
				overlaps[i][j]=0;
			else if(sites[1][i]>=sites[1][j] && sites[1][i]<=sites[2][j])
			{	overlaps[i][j]=1;
				overlaps[j][i]=1;
			}
			else 
				overlaps[i][j]=0;

		}
	}
	return overlaps;
}

double ** removeOverlaps(double ** sites, int numsites, int** overlaps)
{
	int i,j,k, removed = 0;
	for(i=0;i<numsites;i++)
	{
		for(j=i+1;j<numsites;j++)
		{
			//if you find an overlap and neither site has been removed yet,
			//see if it is a multiple overlap or not
			if(overlaps[i][j]==1 && sites[4][j]>0 && sites[4][i]>0)
			{
				for(k=j+1;k<numsites;k++)
				{
					//if multi-overlap and none of the sites has been removed yet,
					//remove weakest site!!!!
					if(overlaps[i][k]==1 && overlaps[j][k]==1 && sites[4][j]>0 && sites[4][i]>0 && sites[4][k]>0)
					{
						removed++;
						if(sites[3][k]<=sites[3][i] && sites[3][k]<=sites[3][j])
							sites[4][k] = -1;
						else if(sites[3][j]<=sites[3][i] && sites[3][j]<=sites[3][k])
							sites[4][j] = -1;
						else if(sites[3][i]<=sites[3][j] && sites[3][i]<=sites[3][k])
							sites[4][i] = -1;
						else {
							cout << "PROBLEM: found multi-overlap but cannot remove it!!!" << endl;
						}

					}
				}
			}
		}
	}
	double ** BSs = new double*[6];
	countnew++;
	for(i=0;i<5;i++)
	{
		BSs[i] = new double[numsites-removed];
		countnew++;
	}
	BSs[5] = new double[1];
	countnew++;
	BSs[5][0] = numsites-removed;
	int countsites = 0;
	//copy over sites, removing those that should be removed
	for(i=0;i<numsites;i++)
	{
		if(sites[4][i]>0)
		{
			BSs[0][countsites] = sites[0][i];
			BSs[1][countsites] = sites[1][i];
			BSs[2][countsites] = sites[2][i];
			BSs[3][countsites] = sites[3][i];
			BSs[4][countsites] = sites[4][i];
			countsites++;
		}
		delete []overlaps[i];
		countnew--;
	}
	delete []overlaps;
	countnew--;
	for(i=0;i<5;i++)
	{
		delete []sites[i];
		countnew--;
	}
	delete []sites;
	countnew--;
	return BSs;
}
					 

void orderMAST(int construct,double ** sites, int numsites) 
{
	int i, j;
	
	int ** overlaps;
	overlaps = findOverlaps(sites,numsites);
	//cout << "numsites = " << numsites << endl;
	double **BSs = removeOverlaps(sites,numsites,overlaps);
	numsites = BSs[5][0];
	//cout << "new numsites = " << numsites << endl;
	double * middleInOrder = new double[numsites];
	countnew++;
	int * counted = new int[numsites];
	countnew++;
	double ** orderedBSs = new double*[5];
	countnew++;
	for(i=0;i<5;i++)
	{
		orderedBSs[i] = new double[numsites];
		countnew++;
	}
	for(i=0;i<numsites; i++)
	{
		middleInOrder[i] = BSs[0][i];
		counted[i] = 0;
	}
	sort(middleInOrder, middleInOrder+numsites);
	for(i=0;i<numsites;i++)
	{
		for(j=0;j<numsites;j++)
		{
			if(middleInOrder[i]==BSs[0][j] && counted[j]==0)
			{
				orderedBSs[0][i] = BSs[0][j];
				orderedBSs[1][i] = BSs[1][j];
				orderedBSs[2][i] = BSs[2][j];
				orderedBSs[3][i] = BSs[3][j];
				orderedBSs[4][i] = BSs[4][j];
				counted[j]=1;
				break;
			}
		}
	}
	for(i=0;i<6;i++)
	{
		delete []BSs[i];
		countnew--;
	}
	delete []BSs;
	countnew--;
	delete []counted;
	countnew--;
	delete []middleInOrder;
	countnew--;
	//createWTinput(orderedBSs,numsites);
	
	if(construct<35)
	{
		int** ff = findFootprints(orderedBSs,numsites);
		for(i=0;i<numsites;i++)
		{
			if(ff[0][i]>0 && construct == 1)
				cout << "found footprint " << ff[1][i] << ff[2][i] << endl;
		}
		
		createInput(orderedBSs,numsites,construct);
	}
	else {
		createInput(orderedBSs,numsites,construct);
	}

	
	for(i=0;i<5;i++)
	{
		delete []orderedBSs[i];
		countnew--;
	}
	delete []orderedBSs;
	countnew--;
}

int ** findFootprints(double const* const* const orderedBSs,int numsites)
{
	int i,j;
	int footprints[4][12];
	// start, end, type, #
	/*
	 D1 - 43-53
	 D2 - 179-189
	 D3 - 238-248
	 D4 - 285-295
	 
	 T1 - 213-222
	 T2 - 224-233
	 
	 S1 - 74-83
	 S2 - 144-153
	 S3 - 165-174 (Rev)
	 S4 - 225-234
	 
	 bHLH1 - 64-70
	 bHLH2 - 87-93
	 */
	//D1
	footprints[0][0] = 43; footprints[1][0] = 53; footprints[2][0] = 1; footprints[3][0] = 1;
	//D2
	footprints[0][1] = 179; footprints[1][1] = 189; footprints[2][1] = 1; footprints[3][1] = 2;
	//D3
	footprints[0][2] = 238; footprints[1][2] = 248; footprints[2][2] = 1; footprints[3][2] = 3;
	//D4
	footprints[0][3] = 285; footprints[1][3] = 295; footprints[2][3] = 1; footprints[3][3] = 4;
	//T1
	footprints[0][4] = 213; footprints[1][4] = 222; footprints[2][4] = 2; footprints[3][4] = 1;
	//T2
	footprints[0][5] = 224; footprints[1][5] = 233; footprints[2][5] = 2; footprints[3][5] = 2;
	//S1
	footprints[0][6] = 74; footprints[1][6] = 83; footprints[2][6] = 3; footprints[3][6] = 1;
	//S2
	footprints[0][7] = 144; footprints[1][7] = 153; footprints[2][7] = 3; footprints[3][7] = 2;
	//S3
	footprints[0][8] = 165; footprints[1][8] = 174; footprints[2][8] = 3; footprints[3][8] = 3;
	//S4
	footprints[0][9] = 225; footprints[1][9] = 234; footprints[2][9] = 3; footprints[3][9] = 4;
	//BHLH
	footprints[0][10] = 64; footprints[1][10] = 70; footprints[2][10] = 2; footprints[3][10] = 3;
	footprints[0][11] = 87; footprints[1][11] = 93; footprints[2][11] = 2; footprints[3][11] = 4;
	
	
	
	int ** foundfootprints = new int*[3];
	foundfootprints[0] = new int[numsites];
	foundfootprints[1] = new int[numsites];
	foundfootprints[2] = new int[numsites];
	// 1 if footprint, type, #
	
	for(i=0;i<numsites;i++)
	{
		foundfootprints[0][i]=0;
		for(j=0;j<12;j++)
		{
			//first make sure they are the same type of site
			if(orderedBSs[4][i]==footprints[2][j] && foundfootprints[0][i] == 0)
			{
				if(orderedBSs[1][i]<=footprints[0][j] && footprints[0][j] <=orderedBSs[2][i])
				{
					//found a footprint
					foundfootprints[0][i] = 1;
					foundfootprints[1][i] = orderedBSs[4][i];
					foundfootprints[2][i] = footprints[3][j];
				}
				else if(footprints[0][j]<=orderedBSs[1][i] && orderedBSs[1][i]<=footprints[1][j])
				{
					//found a footprint
					foundfootprints[0][i] = 1;
					foundfootprints[1][i] = orderedBSs[4][i];
					foundfootprints[2][i] = footprints[3][j];
				}
				else {
					//not a footprint
					foundfootprints[0][i] = 0;
					foundfootprints[1][i] = 0;
					foundfootprints[2][i] = 0;
				}
			}
			
		}
	}
	return foundfootprints;
}

void createInput(double const* const* const orderedBSs,int numsites,int construct)
{
	int i,j,numover=0;
	char c[10];
	int ** overlaps = findOverlaps(orderedBSs,numsites);
	char file[80];
	string n;
	
	n = sprintf(c, "%i", construct);
	strcpy(file,"Inputs_");
	strcat(file,c);
	strcat(file,".txt");
	
	
	ofstream inputs(file);
	inputs << "numBS: " << numsites << endl;
	inputs << "types(D=1,T=2,S=3): ";
	for(i=0;i<numsites;i++)
		inputs << orderedBSs[4][i] << " ";
	inputs << endl;
	inputs << "distances: ";
	for(i=0;i<numsites-1;i++)
	{
		if(overlaps[i][i+1]>0)
		{
			numover++;
			inputs << "0 ";
		}
		else
			inputs << orderedBSs[0][i+1]-orderedBSs[0][i] << " ";
	}
	inputs << endl;
	inputs << "binding_affinities: ";
	for(i=0;i<numsites;i++)
	{	inputs << orderedBSs[3][i] << " ";
//cout << orderedBSs[3][i] << endl;
		if(orderedBSs[3][i]<0)
		{	cout << "PROBLEM: negative Binding Affinity!!!!" << endl;
			if(orderedBSs[4][i]==1)
				cout << "problem is Dorsal" << endl;
			else if(orderedBSs[4][i]==2)
				cout << "problem is Twist" << endl;
			else if(orderedBSs[4][i]==3)
				cout << "problem is Snail" << endl;
			//exit(1);
		}
	}
	inputs << endl;
	inputs << "numOverlapping: " << numover << endl;
	inputs << "Overlapping_distances(from_center_to_center): ";
	for(i=0;i<numsites-1;i++)
	{
		if(overlaps[i][i+1]>0)
		{
			inputs << orderedBSs[0][i+1]-orderedBSs[0][i] << " ";
		}
	}
	inputs << endl;
	
	for(i=0;i<numsites;i++)
	{
		delete []overlaps[i];
		countnew--;
	}
	delete []overlaps;
	countnew--;
	
	inputs.close();
}

void create_thermo_inputs(int NumCons, int typeInter, int typeC, int numCbins,
						  int sizeCbins, int typeQ, int numQbins, int sizeQbins)
{
	int i;
	string lin1,lin2,lin3,lin4,lin5,lin6;
	char name[128];
	char c[10];
	int n;
	char file[80];
	
	
	
	for(i=0;i<NumCons;i++)
	{
		n = sprintf(c, "%i", i+1);
		strcpy(file,"Inputs_");
		strcat(file,c);
		strcat(file,".txt");
		
		ifstream paramsI(file);
		getline(paramsI,lin1);
		getline(paramsI,lin2);
		getline(paramsI,lin3);
		getline(paramsI,lin4);
		getline(paramsI,lin5);
		getline(paramsI,lin6);
		
		ofstream paramsO(file);
		paramsO << lin1 << endl;
		paramsO << lin2 << endl;
		paramsO << lin3 << endl;
		paramsO << lin4 << endl;
		paramsO << lin5 << endl;
		paramsO << lin6 << endl;
		paramsO << "type_Interactions(none=0,NN=1,allN=2): " << typeInter << endl;
		paramsO << "type_Cooperativity(fixed=0,binned=1,ProteinBinned=2,Log=3,Linear=4,Gauss=5): " << typeC << endl;
		if(typeC==1 || typeC==2)
		{
			paramsO << "If_binned_coop_numBins_and_sizeBins: " << numCbins << " " << sizeCbins << endl;
		}
		else
		{
			paramsO << "If_binned_coop_numBins_and_sizeBins: " << endl;
		}
		paramsO << "type_Quenching(fixed=0,binned=1,MSB=2,Log=3,Linear=4,Gauss=5): " << typeQ << endl;
		if(typeQ==1)
		{
			paramsO << "If_binned_quench_numBins_and_sizeBins: " << numQbins << " " << sizeQbins << endl;
		}
		else
		{
			paramsO << "If_binned_quench_numBins_and_sizeBins: " << endl;
		}
		paramsO << "scaling_factors:";
		paramsO << endl << "cooperativities:";
		
		paramsO << endl << "quenching:";
		
		paramsO.close();
	}
}
	
int main(int argc, char * argv[]) 
{
	//Parameters from the command line
	
	
	//-------------------------------------------------
	
	
	int i, construct, con, j, bins, totalsites;
	double threshold = 0,Dthreshold = 0,Tthreshold = 0,Sthreshold = 0;
	char name[128];
	string dummy;
	char c;
	int ncon[47] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,47,48,49,50,52,53,54,55,56,57,58,59,60};
	char file[80];
	int PWMs[3] = {0,0,0};
	double **x1, **x2, **x3, **sites;
	int typeC = -1, typeQ = -1, typeInter = 0, numCbins = 0, numQbins = 0, sizeCbins = 0, sizeQbins = 0;
	
	
	
	for (i = 1; i < argc; i++) {
		string arg = argv[i];
		//cout << "Param " << i << ": " << arg << endl;
		int n = 1;
		while (i+n < argc && argv[i+n][0] != '-') 
			n++;
	
		if(arg == "-h"){
			help_flags();
			exit(1);
		}
		else if (arg == "-thresh") {
			if (n == 2)
			{
				threshold = atof(argv[++i]);
				Dthreshold = threshold;
				Tthreshold = threshold;
				Sthreshold = threshold;
			}
			else 
				cerr << "-thresh requires a threshold number (double)" << endl;
		}
		else if (arg == "-Dthresh") {
			if (n == 2)
				Dthreshold = atof(argv[++i]);
			else 
				cerr << "-Dthresh requires a threshold number (double)" << endl;
		}
		else if (arg == "-Tthresh") {
			if (n == 2)
				Tthreshold = atof(argv[++i]);
			else 
				cerr << "-Tthresh requires a threshold number (double)" << endl;
		}
		else if (arg == "-Sthresh") {
			if (n == 2)
				Sthreshold = atof(argv[++i]);
			else 
				cerr << "-Sthresh requires a threshold number (double)" << endl;
		}
		else if (arg == "-dlPWM") {
			if (n==2) 
			{   
				PWMs[0] = atoi(argv[++i]);
				if(PWMs[0]<1 || PWMs[0]>3)
				{	cerr << "-dlPWM entered is invalid!" << endl;
					exit(1);
				}
			}
			else 
				cerr << "-dlPWM requires an integer" << endl;
		}
		else if (arg == "-twPWM") {
			if (n==2)
			{
				PWMs[1] = atoi(argv[++i]);
				if(PWMs[1]<1 || PWMs[1]>3)
				{	cerr << "-twPWM entered is invalid!" << endl;
					exit(1);
				}
			}
			else
				cerr << "-twPWM requires an integer" << endl;
		}
		else if (arg == "-snPWM") {
			if (n==2)
			{	PWMs[2] = atoi(argv[++i]);
				if(PWMs[2]<1 || PWMs[2]>3)
				{	cerr << "-snPWM entered is invalid!" << endl;
					exit(1);
				}
			}
			else 
				cerr << "-snPWM requires an integer" << endl;
		}
		else if (arg == "-inter") {
			if (n==2)
			{	typeInter = atoi(argv[++i]);
				if(typeInter<0 || typeInter>2)
				{	cerr << "-inter entered is invalid!" << endl;
					exit(1);
				}
			}
			else 
				cerr << "-inter requires an integer" << endl;
		}
		else if (arg == "-coop") {
			if (n==2)
			{	typeC = atoi(argv[++i]);
				if(typeC<0 || typeC>5)
				{	cerr << "-coop entered is invalid!" << endl;
					exit(1);
				}
			}
			else 
				cerr << "-coop requires an integer" << endl;
		}
		else if (arg == "-bincoop") {
			if (n==3)
			{	numCbins = atoi(argv[++i]);
				sizeCbins = atoi(argv[++i]); 
				if(numCbins<1 || sizeCbins<1)
				{	cerr << "-bincoop entered is invalid!" << endl;
					exit(1);
				}
			}
			else 
				cerr << "-bincoop requires 2 integers" << endl;
		}
		else if (arg == "-quench") {
			if (n==2)
			{	typeQ = atoi(argv[++i]);
				if(typeQ<0 || typeQ>5)
				{	cerr << "-quench entered is invalid!" << endl;
					exit(1);
				}
			}
			else 
				cerr << "-quench requires an integer" << endl;
		}
		else if (arg == "-binquench") {
			if (n==3)
			{	numQbins = atoi(argv[++i]);
				sizeQbins = atoi(argv[++i]); 
				if(numQbins<1 || sizeQbins<1)
				{	cerr << "-binquench entered is invalid!" << endl;
					exit(1);
				}
			}
			else 
				cerr << "-binquench requires 2 integers" << endl;
		}
		else {
			cerr << "Parameter \"" << arg << "\" is not recognized." << endl;
		}
	}
	if(Dthreshold<=0 || Tthreshold<=0 || Sthreshold<=0)
	{	
		cerr << "atleast 1 threshold not set" << endl;
		if (Dthreshold <= 0) cerr << "It's Dthreshold" << endl;
		if (Tthreshold <= 0) cerr << "It's Tthreshold" << endl;
		if (Sthreshold <= 0) cerr << "It's Sthreshold" << endl;
		exit(1);
	}
	if(PWMs[0]==0)
	{	
		cerr << "dlPWM not set" << endl;
		exit(1);
	}
	if(PWMs[1]==0)
	{	
		cerr << "twPWM not set" << endl;
		exit(1);
	}
	if(PWMs[2]<=0)
	{	
		cerr << "snPWM not set" << endl;
		exit(1);
	}
	if(typeInter==0)
	{	
		cerr << "type of Interaction not set" << endl;
		exit(1);
	}
	if(typeC<0)
	{	
		cerr << "type of Cooperativity not set" << endl;
		exit(1);
	}
	if(typeQ<0)
	{	
		cerr << "type of Quenching not set" << endl;
		exit(1);
	}
	if(typeC==1 || typeC==2)
	{	
		if(numCbins==0)
		{
			cerr << "number of bins for Cooperativity not set" << endl;
			exit(1);
		}
		if(sizeCbins==0)
		{
			cerr << "size of bins for Cooperativity not set" << endl;
			exit(1);
		}
	}
	if(typeQ==1)
	{	
		if(numQbins==0)
		{
			cerr << "number of bins for Quenching not set" << endl;
			exit(1);
		}
		if(sizeQbins==0)
		{
			cerr << "size of bins for Quenching not set" << endl;
			exit(1);
		}
	}
	
	
	for (con=0;con<59;con++){
		//construct = ncon[con];
		//cout << "construct = " << construct << endl;
		construct = con+1;
	
		runMAST(construct,1,Dthreshold,PWMs);
		runMAST(construct,2,Tthreshold,PWMs);
		runMAST(construct,3,Sthreshold,PWMs);
		x1 = getMASTresults(1);
		x2 = getMASTresults(2);
		x3 = getMASTresults(3);
		if(con==0){
		cout << "# of dorsal = " << x1[0][0] << endl;
		cout << "# of twist = " << x2[0][0] << endl;
		cout << "# of snail = " << x3[0][0] << endl;
		}
		totalsites = x1[0][0]+x2[0][0]+x3[0][0];
		sites = new double*[5];
		countnew++;
		for(i=0;i<5;i++)
		{
			sites[i] = new double[totalsites];
			countnew++;
		}
		
		for(i=0;i<x1[0][0]; i++)
		{
			j = (x1[1][i]+x1[2][i])/2;
			sites[0][i] = j;
			sites[1][i] = x1[1][i];
			sites[2][i] = x1[2][i];
			sites[3][i] = x1[3][i];
			sites[4][i] = 1;
		}
		for(i=0;i<x2[0][0]; i++)
		{
			j = (x2[1][i]+x2[2][i])/2;
			sites[0][(int)x1[0][0]+i] = j;
			sites[1][(int)x1[0][0]+i] = x2[1][i];
			sites[2][(int)x1[0][0]+i] = x2[2][i];
			sites[3][(int)x1[0][0]+i] = x2[3][i];
			sites[4][(int)x1[0][0]+i] = 2;
		}
		
		for(i=0;i<x3[0][0]; i++)
		{
			j=(x3[1][i]+x3[2][i])/2;
			sites[0][(int)x1[0][0]+(int)x2[0][0]+i] = j;
			sites[1][(int)x1[0][0]+(int)x2[0][0]+i] = x3[1][i];
			sites[2][(int)x1[0][0]+(int)x2[0][0]+i] = x3[2][i];
			sites[3][(int)x1[0][0]+(int)x2[0][0]+i] = x3[3][i];
			sites[4][(int)x1[0][0]+(int)x2[0][0]+i] = 3;
		}
		
		for(i=0;i<4;i++)
		{
			delete[]x1[i];
			countnew--;
			delete[]x2[i];
			countnew--;
			delete[]x3[i];
			countnew--;
		}
		delete[]x1;
		countnew--;
		delete[]x2;
		countnew--;
		delete[]x3;
		countnew--;
		
		orderMAST(construct,sites,totalsites);
	}
	create_thermo_inputs(60,typeInter,typeC,numCbins,sizeCbins,typeQ,numQbins,sizeQbins);
	return 0;
	
}