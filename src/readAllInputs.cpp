/*
 *  readAllInputs.cpp
 *  
 *
 *  Created by Jacqueline Dresch on 10/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

using namespace std;


int main(int nNumberofArgs, char* pszArgs[])
{
	int i, construct, numC, n, m, nn, con, j;
	double d;
	char name[128];
	char c[10];
	char file[80];
	double dl[17],tw[17],sn[17];
	
	
	ofstream outputs("All_Inputs.h");
	ifstream inputs;
	ifstream express;
	
	outputs << "#ifndef ALL_INPUTS_H" << endl;
	outputs << "#define ALL_INPUTS_H" << endl;
	outputs << "#include <string.h>" << endl;
	outputs << "#include <stdio.h>" << endl;
	outputs << "#include <stdlib.h>" << endl;
	outputs << "#include <iostream>" << endl;
	outputs << "#include <iomanip>" << endl;
	outputs << "#include <fstream>" << endl;
	outputs << "#include <math.h>" << endl;
	//outputs << "#include \"tbb/parallel_for.h\"" << endl;
	//outputs << "#include \"tbb/blocked_range.h\"" << endl;
	outputs << "#include \"ThermoModel.h\"" << endl;
	outputs << "Params* thermo_parameters = new Params[62];" << endl;
	//outputs << "extern double** Alloutcomes;" << endl;
	outputs << "double Concentrations[3][17];" << endl;
	outputs << "double Expressions[62][17];" << endl << endl;
	
	outputs << "using namespace std;" << endl << endl;
	outputs << "void initialize_thermo(void) {" << endl;
	
	
	for(con=0;con<60;con++)
	{
		construct = con+1;
		n = sprintf(c, "%i", construct);
		strcpy(file,"Inputs_");
		strcat(file,c);
		strcat(file,".txt");
		inputs.open(file);
		if(!inputs.is_open())
		{
			cout << "Couldn't find " << file << "\n";
			exit(1);
		}
		
		strcpy(file,"Expression_");
		strcat(file,c);
		strcat(file,".txt");
		express.open(file);
		if(!express.is_open())
		{
			cout << "CAN NOT find Expression.txt\n";
		}
		
	
		char name[128];
		outputs << "thermo_parameters[" << con << "].n = ";
		//input # of binding sites
		inputs >>  name >> n;
		outputs << n <<  ";" << endl;
		/*input the type of binding sites
		 1 = Dl
		 2 = Twi
		 3 = Sna
		 4 = Twi/Sna*/
		inputs >>  name;
		outputs << "thermo_parameters[" << con << "].typeBS = new int[" << n << "];" << endl;
		//cout << "type of binding sites: ";
		for (i=0; i<n; i++) {
			inputs >> m;
			outputs << "thermo_parameters[" << con << "].typeBS[" << i << "] = "<< m << ";" << endl;
		}
		//input the distances between binding sites
		outputs << "thermo_parameters[" << con << "].distances = new int[" << n-1 << "];" << endl;
		inputs >>  name;
		for (i=0; i<n-1; i++) {
			inputs >> m;
			outputs << "thermo_parameters[" << con << "].distances[" << i << "] = "<< m << ";" << endl;
		}
		//input Binding Affinity for each site
		outputs << "thermo_parameters[" << con << "].BAs = new double[" << n << "];" << endl;
		inputs >>  name;
		for (i=0; i<n; i++) {
			inputs >> d;
			outputs << "thermo_parameters[" << con << "].BAs[" << i << "] = " << d << ";" << endl;
		}
		//input the number of overlapping sites
		inputs >>  name >> m;
		outputs << "thermo_parameters[" << con << "].nOverlap = " << m << ";" << endl;
		//input distances between overlapping binding sites (from middle to middle)
		inputs >>  name;
		outputs << "thermo_parameters[" << con << "].overlapDistances = new int[" << m << "];" << endl;
		for (i=0; i<m; i++) {
			inputs >> nn;
			outputs << "thermo_parameters[" << con << "].overlapDistances[" << i << "] = " << nn << ";" << endl;
			
		}
		inputs >>  name;
		inputs >> m;
		outputs << "thermo_parameters[" << con << "].typeInteractions = " << m << ";" << endl;
		inputs >>  name;
		int typeC, typeQ,numCbins,numQbins,sizeCbins,sizeQbins;
		inputs >> typeC;
		inputs >>  name;
		if(typeC==1 || typeC==2)
		{
			inputs >> numCbins >> sizeCbins;
		}
		inputs >> name;
		inputs >> typeQ;
		inputs >>  name;
		if(typeQ==1)
		{
			inputs >> numQbins >> sizeQbins;
		}	
		
		if(typeC==0 || typeC==3 ||typeC==4 || typeC==5)
		{
			outputs << "thermo_parameters[" << con << "].Cs = new double[5];" << endl;
			outputs << "thermo_parameters[" << con << "].Cs[0] = " << typeC << ";" << endl;
		}
		//binned coop
		else if(typeC==1)
		{
			outputs << "thermo_parameters[" << con << "].Cs = new double[" << numCbins*2+3 << "];" << endl;
			outputs << "thermo_parameters[" << con << "].Cs[0] = " << typeC << ";" << endl;
			outputs << "thermo_parameters[" << con << "].Cs[1] = " << numCbins << ";" << endl;
			outputs << "thermo_parameters[" << con << "].Cs[2] = " << sizeCbins << ";" << endl;
		}
		//binned coop for EACH protein-protein interaction (Dl-Dl,Tw-Tw,Sn-Sn,Dl-Tw)
		else if(typeC==2)
		{
			outputs << "thermo_parameters[" << con << "].Cs = new double[" << numCbins*4+3 << "];" << endl;
			outputs << "thermo_parameters[" << con << "].Cs[0] = " << typeC << ";" << endl;
			outputs << "thermo_parameters[" << con << "].Cs[1] = " << numCbins << ";" << endl;
			outputs << "thermo_parameters[" << con << "].Cs[2] = " << sizeCbins << ";" << endl;
			
		}
		if(typeQ==0 || typeQ==3 || typeQ==4 || typeQ==5)
		{
			outputs << "thermo_parameters[" << con << "].Qs = new double[3];" << endl;
			outputs << "thermo_parameters[" << con << "].Qs[0] = " << typeQ << ";" << endl;
		}
		else if(typeQ==1)
		{
			outputs << "thermo_parameters[" << con << "].Qs = new double[" << numQbins*2+3 << "];" << endl;
			outputs << "thermo_parameters[" << con << "].Qs[0] = " << typeQ << ";" << endl;
			outputs << "thermo_parameters[" << con << "].Qs[1] = " << numQbins << ";" << endl;
			outputs << "thermo_parameters[" << con << "].Qs[2] = " << sizeQbins << ";" << endl;
		}
		else if(typeQ==2)
		{
			outputs << "thermo_parameters[" << con << "].Qs = new double[1];" << endl;
			outputs << "thermo_parameters[" << con << "].Qs[0] = " << typeQ << ";" << endl;
		}
		
		outputs << "thermo_parameters[" << con << "].Scalings = new double[3];" << endl;
		
		inputs.close();
	
		
		express >> name;
		for(i=0;i<17;i++)
		{
			express >> d;
			outputs << "Expressions[" << con << "][" << i << "] = " << d << ";" << endl;
		}
		express.close();
	}
	
	ifstream concent;
	//read file "Concentrations.txt" Do not create the file
	//if it doesn't already exist
	
	concent.open("Concentrations.txt");
	if(!concent.is_open())
	{
		cout << "Concentrations.txt NOT found\n";
	}
	concent >>  name;
	for (i=0; i<17;i++) {
		concent >>  dl[i];
		outputs << "Concentrations[0][" << i<< "] = " << dl[i] << ";" << endl;
	}
	concent >>  name;
	for (i=0; i<17;i++) {
		concent >>  tw[i];
		outputs << "Concentrations[1][" << i<< "] = " << tw[i] << ";" << endl;
	}
	//cout << endl << "Snail: ";
	concent >>  name;
	for (i=0; i<17;i++) {
		concent >>  sn[i];
		outputs << "Concentrations[2][" << i<< "] = " << sn[i] << ";" << endl;
	}
	concent.close();
	outputs << "}" << endl;
		
	
	outputs << endl << endl << "#endif" << endl;
	
	//cout << "// = " << // << endl;
	return 0;
}

