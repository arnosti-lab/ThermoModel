
#include "ThermoModel.h"
#include "All_Inputs.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>
using namespace std;


bool equals(double num1, double num2, double tol = numeric_limits<float>::epsilon() )
{
	if(fabs(num1-num2)<tol)
	{
		return true;
	}
	else {
		return false;
	}
}


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
	//cout << "the Best error = " << best << endl;
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
	while (cmaes2 && foundit == 0 ) { // keep reading until end-of-file
		//for(j=0;j<1;j++)
		//{
		count++;
		//cout << "line = " << count << endl;
		for(i=0; i<(5+numP); i++)
		{
			cmaes2 >> num;
			//cout << "num = " << num << endl;
			if(i==4 && equals(num,error)){
				foundit = 1;
				for(j=0; j<numP; j++)
					cmaes2 >> x[j];
				break;
			}
			//cout << "Best error so far = " << best << endl;
		}		
	}
	if(foundit == 0)
	{
		cout << "ERROR: could not find the best error in the file \"outcmaesxrecentbest.dat\"" << endl;
		exit(1);
	}
	//cout << "the Best parameter set is = ";
	//for(j=0; j<numP; j++){
//		cout <<  x[j] << ", ";
//	}
	
//	cout << endl;
	
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
	for(i=0; i<38; i++)
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
	
	double bestError;
	bestError = findBestError();
	//cout << "bestError = " << bestError << endl;
	x = findBestParameters(bestError,N);
	
	
	
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
	double check = runobjective(x,N,typeC,typeQ);
	//cout << "check error = " << check << endl;
	if(equals(check,bestError,0.000001))
		cout << "Everything checks out!" << endl;
	else
		cout << "PROBLEM!!!!!  Errors are not the same!!!!" << endl;
	
	return 0;
	
}