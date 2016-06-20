
#include "ThermoModel.h"
#include "All_Inputs.h"
#include <iostream>
#include <fstream>
#include <math.h>

double findBestError(void)
{
	string words;
	int i, j,count = 0;
	double num[13], best = 540;
	char name[1000];
	ifstream cmaes1;
	//read file "outcmaesfit.dat" Do not create the file
	//if it doesn't already exist
	
	cmaes1.open("outcmaesfit.dat");
	if(cmaes1.bad())
	{
		cout << "outcmaesfit.dat\n";
	}
	getline(cmaes1,words);
	//cout << words << endl;
	int atend = 0;
	while ( !cmaes1.eof() ){ // keep reading until end-of-file
	//for(j=0;j<1;j++)
	//{
		count++;
		//cout << "line = " << count << endl;
		for(i=0; i<13; i++)
		{
			cmaes1 >> num[i];
			//cout << "num = " << num << endl;
			/*if(i==0 && num[0] !=count){
				atend = 1;
				break;
			}
			/*	
			if(i==4 && num<best)
				best = num;
			//cout << "Best error so far = " << best << endl;*/
		}		
		
	}
	best = num[4];
	cout << "the Best error = " << best << endl;
	return best;
}

double* findBestParameters(double error, int numP)
{
	string words;
	int i, j,count = 0;
	double* x = new double[numP];
	double num;
	ifstream cmaes2;
	//read file "outcmaesxrecentbest.dat" Do not create the file
	//if it doesn't already exist
	
	cmaes2.open("outcmaesxrecentbest.dat");
	if(cmaes2.bad())
	{
		cout << "outcmaesxrecentbest.dat\n";
	}
	getline(cmaes2,words);
	//cout << words << endl;
	int foundit = 0;
	while ( foundit == 0 ) { // keep reading until end-of-file
		//for(j=0;j<1;j++)
		//{
		count++;
		//cout << "line = " << count << endl;
		for(i=0; i<(5+numP); i++)
		{
			cmaes2 >> num;
			//cout << "num = " << num << endl;
			if(i==4 && num==error){
				foundit = 1;
				for(j=0; j<numP; j++)
					cmaes2 >> x[j];
				break;
			}
			//cout << "Best error so far = " << best << endl;
		}		
	}
	cout << "the Best parameter set is = ";
	for(j=0; j<numP; j++){
		cout <<  x[j] << ", ";
	}
	
	cout << endl;
	
	return x;
	
}

int matrixExpression(void)
{
	char name[128];
	double express[17];
	char c[10];
	int i, j, n, construct;
	char file[80];
	
	int N, typeC,typeQ;
	double* x;
	
	
	ifstream inputs;
	ofstream outputs;
	//read file "Expresssion_#.txt" Do not create the file
	//if it doesn't already exist
	for(i=0; i<59; i++)
	{
		construct = i+1;
		n = sprintf(c, "%i", construct);
		strcpy(file,"Expression_");
		strcat(file,c);
		strcat(file,".txt");
		
		inputs.open(file);
		if(inputs.bad())
		{
			cout << file << " does not exist" << endl;
		}
		
		strcpy(file,"MExpression_");
		strcat(file,c);
		strcat(file,".txt");
		outputs.open(file);
		inputs >> name;
		for (j=0; j<17;j++) {
			inputs >>  express[j];
			outputs << express[j] << " ";
		}
		
		inputs.close();
		outputs.close();
		
	}
	
	return 0;
}



int main(int nNumberofArgs, char* pszArgs[])
{
	int i, construct, numC, numQ, n, con, j, bins;
	double totalSSE = 0;
	string lin1,lin2,lin3,lin4,lin5,lin6,lin7,lin8,lin9,lin10,lin11,lin12;
	char name[128];
	Params p[17];
	double dl[17], tw[17], sn[17], express[17];
	char c[10];
	char file[80];
	
	initialize_thermo();
	
	int N, typeC,typeQ;
	double* x;
	
	strcpy(file,"Inputs_1");
	strcat(file,".txt");
	 
	ifstream paramsI(file);
	getline(paramsI,lin1);
	getline(paramsI,lin2);
	getline(paramsI,lin3);
	getline(paramsI,lin4);
	getline(paramsI,lin5);
	getline(paramsI,lin6);
	getline(paramsI,lin7);
	paramsI >> name;
	paramsI >> typeC;
	getline(paramsI,lin8);
	paramsI >> name;
	if(typeC == 1 || typeC ==2)
		paramsI >> numC;
	getline(paramsI,lin9);
	paramsI >> name;
	paramsI >> typeQ;
	paramsI >> name;
	if(typeQ == 1)
		paramsI >> numQ;
	if(typeC==0 || typeC==3 || typeC==4 ||typeC==5)
		numC = 4;
	else if(typeC==1)
		numC = numC*2;
	else if(typeC==2)
		numC = numC*4;
	if(typeQ==0 || typeQ==3 || typeQ==4 || typeQ==5)
		numQ = 2;
	else if(typeQ==1) 
		numQ = numQ*2;
	else if(typeQ==2)
		numQ=0;
	N = 3+numC+numQ;
	cout << "typeC = " <<typeC <<", numC = " << numC << " ,typeQ = " << typeQ  << " and numQ = " << numQ<< endl;
	
	
	x = findBestParameters(findBestError(),N);
	
	/*try same as ThermoExample.c
	double *pos = new double[N];
	for(i=0; i<N; i++)
		pos[i] = x[i]*x[i];
	for (i=N-20; i<N; i++)
		pos[i] = exp(-pos[i]);
	cout << "Check ThermoExample.c: " << endl;
	for(i=0; i<N; i++)
		cout << pos[i] << endl;
	 */
	

	ofstream Bestp("parameters_found.txt");
	for(j=0; j<N; j++)
		x[j] =  x[j]*x[j];
	
	
	double lowerSF = 0, lowerC = 0;
	
	
	for(j=0; j<3;j++)
		x[j] = x[j]+lowerSF;
	for(j=3; j<3+numC;j++)
		x[j] = x[j]+lowerC;
	
	if(numQ>2){
		for(j=N-numQ; j<N; j++)
		{
			if(exp(-x[j])<1.0e-20)
				x[j] = 0;
			else
				x[j] =  exp(-x[j]);
		}
	}
	for(j=0; j<N; j++)
		Bestp << x[j] << endl;
	double check = runobjective(x,N,typeC,typeQ);
	cout << "check error = " << check << endl;
	
	ifstream inputs;
	//read file "Concentrations.txt" Do not create the file
	//if it doesn't already exist
	
	inputs.open("Concentrations.txt");
	if(inputs.bad())
	{
		cout << "Concentrations.txt\n";
	}
	//cout << "Dorsal: ";
	inputs >>  name;
	for (i=0; i<17;i++) {
		inputs >>  dl[i];
		//cout << dl[i] << "\t";
	}
	//cout << endl << "Twist: ";;
	inputs >>  name;
	for (i=0; i<17;i++) {
		inputs >>  tw[i];
		//cout << tw[i] << "\t";
	}
	//cout << endl << "Snail: ";
	inputs >>  name;
	for (i=0; i<17;i++) {
		inputs >>  sn[i];
		//cout << sn[i] << "\t";
	}
	totalSSE = 0;
	for(con=0;con<59;con++)
	{
		construct = con+1;
		cout << "construct = " << construct << endl;
		
		n = sprintf(c, "%i", construct);
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
		if(construct==1)
		{	getline(paramsI,lin7);
			getline(paramsI,lin8);
			getline(paramsI,lin9);
			getline(paramsI,lin10);
			getline(paramsI,lin11);
		}
		
		ofstream paramsO(file);
		paramsO << lin1 << endl;
		paramsO << lin2 << endl;
		paramsO << lin3 << endl;
		paramsO << lin4 << endl;
		paramsO << lin5 << endl;
		paramsO << lin6 << endl;
		paramsO << lin7 << endl;
		paramsO << lin8 << endl;
		paramsO << lin9 << endl;
		paramsO << lin10 << endl;
		paramsO << lin11 << endl;
		paramsO << "scaling_factors:";
		for(i=0;i<3;i++){
			paramsO << " " << x[i];
			thermo_parameters[con].Scalings[i] = x[i];
		}
		paramsO << endl << "cooperativities:";
		for(i=0;i<numC;i++){
			paramsO << " " << x[3+i];
			if(typeC ==1 || typeC ==2)
				thermo_parameters[con].Cs[i+3] = x[3+i];
			else 
				thermo_parameters[con].Cs[i+1] = x[3+i];
		}
		paramsO << endl << "quenching:";
		if(typeQ==0 || typeQ==3 || typeQ==4 || typeQ==5)
		{
			thermo_parameters[con].Qs[1] = x[3+numC];
			thermo_parameters[con].Qs[2] = x[3+numC+1];
		}
		else if(typeQ==1)
		{
			for(i=0;i<numQ;i++)
			{
				thermo_parameters[con].Qs[3+i] = x[3+numC+i];
			}
		}
		
		for(i=0;i<numQ;i++)
			paramsO << " " << x[3+numC+i];
		paramsO.close();
		
		//cout << endl << "Expression: " << endl;;
		
		for (i=0; i<17;i++) {
			p[i].SetParameters(construct);
			p[i].SetConcentrations(dl[i],tw[i],sn[i]);
			p[i].saveStates();
			express[i] = p[i].expression();
			delete []p[i].typeBS;
			//countnew--;
			delete []p[i].distances;
			//countnew--;
			delete []p[i].BAs;
			//countnew--;
			delete []p[i].overlapDistances;
			//countnew--;
			delete []p[i].Scalings;
			//countnew--;
			delete []p[i].Cs;
			//countnew--;
			delete []p[i].Qs;
			//countnew--;
			int maxstates = pow(2,p[i].n);
			for(j=0;j<maxstates;j++)
			{
				delete []p[i].allStates[j];
				//countnew--;
			}
			delete []p[i].allStates;
			//countnew--;
		}
		strcpy(file,"Prediction_");
		strcat(file,c);
		strcat(file,".txt");
		ofstream prediction(file);
		for (i=0; i<17; i++) 
			prediction << express[i] << " ";
		totalSSE = totalSSE + SSE(express,construct);
	}
	cout << "totalSSE = " << totalSSE << endl;
	//cout << "countnew = " << countnew << endl;*/
	matrixExpression();
	return 0;
	
}
