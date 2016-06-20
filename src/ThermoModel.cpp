#include "ThermoModel.h"

void Params::SetParameters(int construct) {
    n = thermo_parameters[construct-1].n;
    /*input the type of binding sites
    1 = Dl
    2 = Twi
    3 = Sna
    4 = Twi/Sna*/
    typeBS = new int[n];
    for (int i=0; i<n; i++) typeBS[i] = thermo_parameters[construct-1].typeBS[i];
    //input the distances between binding sites
    distances = new int[n-1];
    for (int i=0; i<n-1; i++) distances[i] = thermo_parameters[construct-1].distances[i];
    //input Binding Affinity for each site
    BAs = new double[n];
    for (int i=0; i<n; i++) BAs[i] = thermo_parameters[construct-1].BAs[i];
    nOverlap = thermo_parameters[construct-1].nOverlap;;
    //input distances between overlapping binding sites (from middle to middle)
    overlapDistances = new int[nOverlap];
    for (int i=0; i<nOverlap; i++) overlapDistances[i] = thermo_parameters[construct-1].overlapDistances[i];
    /*input type of interactions
    NO interactions (turn off all cooperativity and quenching
    typeInteractions = 0
    ONLY allow adjacent bound proteins to interact (coop or quench)
    typeInteractions = 1
    allow ANY proteins bound together to interact (coop or quench)
    typeInteractions = 2*/
    typeInteractions = thermo_parameters[construct-1].typeInteractions;
    int typeC, typeQ,numCbins,numQbins,sizeCbins,sizeQbins;
    typeC = thermo_parameters[construct-1].Cs[0];
    typeQ = thermo_parameters[construct-1].Qs[0];;
    Scalings = new double[3];
    for (int i=0; i<3; i++) Scalings[i] = thermo_parameters[construct-1].Scalings[i];
    /*input cooperativity parameters
    Cs[0] = Dl-Dl coop
    Cs[1] = Twi-Twi coop
    Cs[2] = Sna-Sna coop
    Cs[3] = Dl-Twi coop */
    //constant coop, logistic coop, linear coop, or gaussian coop
    switch(typeC) {
        case 0:
        case 3:
        case 4:
        case 5:
            Cs = new double[5];
            Cs[0] = typeC;
            for (int i=1; i<5; i++) Cs[i] = thermo_parameters[construct-1].Cs[i];
            break;
        //binned coop
        case 1:
            numCbins = thermo_parameters[construct-1].Cs[1];
            Cs = new double[numCbins*2+3];
            Cs[0] = typeC;
            Cs[1] = numCbins;
            Cs[2] = thermo_parameters[construct-1].Cs[2];
            for (int i=3; i<numCbins*2+3; i++) Cs[i] = thermo_parameters[construct-1].Cs[i];
            break;
        //binned coop for EACH protein-protein interaction (Dl-Dl,Tw-Tw,Sn-Sn,Dl-Tw)
        case 2:
            numCbins = thermo_parameters[construct-1].Cs[1];
            Cs = new double[numCbins*4+3];
            Cs[0] = typeC;
            Cs[1] = numCbins;
            Cs[2] = thermo_parameters[construct-1].Cs[2];
            for (int i=3; i<numCbins*4+3; i++) Cs[i] = thermo_parameters[construct-1].Cs[i];
            break;
    }
    switch(typeQ) {
        case 0:
        case 3:
        case 4:
        case 5:
            Qs = new double[3];
            Qs[0] = typeQ;
            Qs[1] = thermo_parameters[construct-1].Qs[1];
            Qs[2] = thermo_parameters[construct-1].Qs[2];
            break;
        case 1:
            numQbins = thermo_parameters[construct-1].Qs[1];
            Qs = new double[numQbins*2+3];
            Qs[0] = typeQ;
            Qs[1] = numQbins;
            Qs[2] = thermo_parameters[construct-1].Qs[2];
            for (int i=3; i<numQbins*2+3; i++) Qs[i] = thermo_parameters[construct-1].Qs[i];
            break;
        case 2:
            Qs = new double[1];
            Qs[0] = typeQ;
            break;
    }
    int maxstates = pow(2,n);
    allStates = new bool*[maxstates];
    for(int i=0;i<maxstates;i++) allStates[i] = new bool[n];
}

void Params::SetConcentrations(double dln, double twn, double snn) {	
    dl = dln;
    tw = twn;
    sn = snn;
}

double quench(int p1,int p2, double d, double* qs) {
    // inputs:
    //		p1 = type of protein 1
    //		p2 = type of protein 2
    //		d = distance between proteins
    int typeQ = qs[0];
    double Q;
    switch(typeQ) {
        case 0:
            //NO distance dependence for now.....but could be added later
            // Sn and Dl
            if(p1-p2>1 || p2-p1>1) Q = qs[1];
            // Sn and Tw
            else Q = qs[2];
            break;
        case 1:
            //Snail can quench Twist and Dorsal with different strengths, 
            //quenching is limited by numQbins and sizeQbins
            // Sn and Dl
            if(p1-p2>1 || p2-p1>1) Q = binnedQuench(1,d,qs);
            // Sn and Tw
            else Q = binnedQuench(2,d,qs);
            break;
        case 2:
            //Quenching has NO parameters.....is just the MSB (Fakhouri et al)
            //quenching function
            Q = -4.2494e-09*pow(d,5)+1.0199e-06*pow(d,4)-1.0199e-04*pow(d,3)+5.0099e-03*pow(d,2)-1.0731e-01*d+1.0;
            if(Q<0 || d>100) Q = 0;
            break;
        case 3:
            //logistic quenching
            // Sn and Dl
            if(p1-p2>1 || p2-p1>1) Q = 2/(1+exp(-d/qs[1]))-1;
            // Sn and Tw
            else Q = 2/(1+exp(-d/qs[2]))-1;
            Q = 1-Q;
            break;
        case 4:
            //linear quenching
            // Sn and Dl
            if(p1-p2>1 || p2-p1>1) Q = qs[1]*d;
            // Sn and Tw
            else Q = qs[2]*d;
            if(Q>1) Q = 1;
            Q = 1-Q;
            break;
        case 5:
            //Gaussian quenching
            // Sn and Dl
            if(p1-p2>1 || p2-p1>1) Q = exp(-d*d/qs[1]);
            // Sn and Tw
            else Q = exp(-d*d/qs[2]);
            break;
    }
    return Q;
}

double binnedQuench(int A, int d, double* quench) {
    // Sna-Dl: A=1
    // Sna-Twi: A=2
    double Q = 0;
    int numBins = quench[1];
    int sizeBins = quench[2];
    double Qd[numBins],Qt[numBins];
    for(int i=0;i<numBins;i++) {
        Qd[i] = quench[i+3];
        Qt[i] = quench[i+3+numBins];
    }
    for(int i=0;i<numBins;i++) {
        if(d>sizeBins*i && d<=sizeBins*(i+1)) {
            if(A==1) Q = Qd[i];
            else Q = Qt[i];
        }
    }
    return Q;
}

double coop(int p1, int p2, double d, double* coops) {
    // inputs:
    //		p1 = type of protein 1
    //		p2 = type of protein 2
    //		d = distance between proteins
    //		coops = array of cooperativities
    double C;
    int Ctype = coops[0];
    switch(Ctype) {
        case 0:
            // NO distance dependence for now.....but could be added later
            if(d<100) {
                // both Dl
                if(p1==1 && p2 ==1) C = coops[1];
                // both Tw
                else if(p1==2 && p2 ==2) C = coops[2];
                // both Sn
                else if(p1==3 && p2 ==3) C = coops[3];
                // one Dl, one Tw
                else C = coops[4];
            }
            else  C = 1.0;
            break;
        case 1: 
            //binned coop, last bin is INFINITE
            //homotypic
            if((p1-p2)==0) C = binnedCoop(0, d, coops);
            //heterotypic
            else C = binnedCoop(1, d, coops);
            break;
        case 2:
            //binned coop for EACH type of interaction
            //homotypic, Dl-Dl
            if(p1==1 && p2==1) C = ProteinbinnedCoop(1,d,coops);
            //homotypic, Twi-Twi
            else if(p1==2 && p2==2) C = ProteinbinnedCoop(2,d,coops);
            //homotypic, Sna-Sna
            else if(p1==3 && p2==3) C = ProteinbinnedCoop(3,d,coops);
            //heterotypic
            else C = ProteinbinnedCoop(4,d,coops);
            break;
        case 3:
            //logistic coop - NOTE: we don't know scale, so we have a vertical stretch
            //and a horizontal stretch!
            //homotypic
            if((p1-p2)==0) C = coops[2]*(1-(2/(1+exp(-d/coops[1]))-1));
            //heterotypic
            else C = coops[4]*(1-(2/(1+exp(-d/coops[3]))-1));
            break;
        case 4:
            //linear coop - NOTE: we don't know scale, so we have an intercept
            //and a negative slope!
            //homotypic
            if((p1-p2)==0) C = (-coops[1])*d+(-coops[2]);
            //heterotypic
            else C = (-coops[3])*d+(-coops[4]);
            if(C<0) C = 0;
            break;
        case 5:
            //Gaussian coop- NOTE: we don't know scale, so we have a vertical stretch
            //and a horizontal stretch!
            //homotypic
            if((p1-p2)==0) C = coops[2]*exp(-d*d/coops[1]);
            //heterotypic
            else C = coops[4]*exp(-d*d/coops[3]);
            break;
    }
    return C;
}
        
double binnedCoop(int HH, int d, double* coops) {
    // Homotypic: HH=0
    // Heterotypic: HH=1
    double C;
    int numBins = coops[1];
    int sizeBins = coops[2];
    for(int i=0;i<numBins-1;i++) if(d>sizeBins*i && d<=sizeBins*(i+1)) C = coops[i+3+numBins*HH];
    if(d>sizeBins*(numBins-1)) C = coops[numBins*(HH+1)-1+3];
    return C;
}

double ProteinbinnedCoop(int proteins, int d, double* coops) {
    //Dl-Dl: proteins = 1
    //Twi-Twi: proteins = 2
    //Sna-Sna: proteins = 3
    //Dl-Twi: proteins = 4
    double C;
    int numBins = coops[1];
    int sizeBins = coops[2];
    for(int i=0;i<numBins-1;i++) if(d>sizeBins*i&& d<=sizeBins*(i+1)) C = coops[i+3+numBins*(proteins-1)];
    if(d>sizeBins*(numBins-1)) C = coops[numBins*proteins-1+3];
    return C;
}

double Params::expression() {
    double conts[2] = {0.0, 0.0};
    contributionsAllStates(*this, conts);
    double result = conts[0]/(conts[1]+1);
    return result;
}

void Params::saveStates() {
    bool *v = new bool[n];
    bool *vcopy = new bool[n];
    for (int i=0;i<n;i++) {
        v[i] = 0;
        vcopy[i] = 0;
    }
    numStates = 0;
    States(n,v,vcopy);
    delete []v;
    delete []vcopy;
}

void Params::States(int nn, bool *v, bool  *originalv) {
    //records all possible states into p.allStates
    //2 possibilities: bound and unbound - but must beware of overlapping sites!!!
    //check if it is an overlapping site - if yes and overlapping site is bound, cannot bind!
    //overlapping site is bound
    if(nn<n && nn>0 && distances[nn-1]==0 && v[nn]==1) {
        //not bound
        for(int j=0;j<nn-1;j++) v[j]=originalv[j];
        //you are not at the very first binding site (meaning you have not gone through all)
        if(nn>1) States(nn-1,v,originalv);
        //you are at the very first binding site (meaning you have gone through all)
        else addState(v);
    }
    //either not overlapping or overlapping site not bound 
    else {
        //not bound
        for(int j=0;j<nn-1;j++) v[j]=originalv[j];
        //you are not at the very first binding site (meaning you have not gone through all)
        if(nn>1) States(nn-1,v,originalv);
        //you are at the very first binding site (meaning you have gone through all)
        else addState(v);
        //bound
        v[nn-1] = 1;
        for(int j=0;j<nn-1;j++) v[j]=originalv[j];
        //you are not at the very first binding site (meaning you have not gone through all)
        if(nn>1) States(nn-1,v,originalv);
        //you are at the very first binding site (meaning you have gone through all)
        else addState(v);
    }
}

void Params::addState(bool *v) {
    for(int i=0;i<n;i++) allStates[numStates][i] = v[i];
    numStates++;
}

void contributionsAllStates(Params const & p, double *final) {
    double *outcomes = new double[2];
    int states = p.numStates/64;
    for(int j=0; j<64; j++) {
        for(int i=j*states;i<(j+1)*states;i++) {
            getContributions(p.allStates[i],p, outcomes);
            final[0] += outcomes[0];
            final[1] += outcomes[1];
        }
    }
    delete[] outcomes;
}

void getContributions(bool const * const v,Params const & p, double *express) {
    if(p.typeInteractions==0) NoNcontributions(v,p, express);
    else if(p.typeInteractions==1) NearNcontributions(v,p, express);
    else if(p.typeInteractions==2) allNcontributions(v,p, express);
}


void NoNcontributions(bool const * const v,Params const & p, double *express) {
    int i, j;
    int numA = 0, numR = 0, num0 = 0;
    int nbound,protbound,pos,prot;
    double con_scaling[3], BA;
    con_scaling[0] = p.dl*p.Scalings[0];
    con_scaling[1] = p.tw*p.Scalings[1];
    con_scaling[2] = p.sn*p.Scalings[2];
    //check what is bound 
    for(i=0;i<p.n;i++) {
        //not bound
        if(v[i]==0) num0++;
        //we have an activator bound
        else if(v[i]>0 && p.typeBS[i]<3) numA++;
        //we have a repressor bound
        else if(v[i]>0 && p.typeBS[i]==3) numR++;
    }
    nbound = p.n-num0;
    // configuration.proteins is all 0's.....nothing is bound
    if(num0==p.n) {
        express[0] = 0;
        express[1] = 0;
    }
    //atleast 1 activator is bound
    else if(numA>0) {
        express[0] = 1.0;
        express[1] = 1.0;
        for(pos=0;pos<p.n;pos++) {
            //find BOUND proteins
            if(v[pos]==0) continue;
            //get type of protein
            prot = p.typeBS[pos];
            //need contribution of scaling and con_scalingration for protein
            express[0] *= con_scaling[prot-1]*p.BAs[pos];
            express[1] *= con_scaling[prot-1]*p.BAs[pos];
        }
    }
    //no activators bound BUT atleast 1 repressor
    else if(numA==0 && numR>0) {
        express[0] = 0.0;
        express[1] = 1.0;
        for(pos=0;pos<p.n;pos++) {
            //find BOUND proteins
            if(v[pos]==0) continue;
            //get type of protein
            prot = p.typeBS[pos];
            //need contribution of scaling and con_scalingration for protein
            express[1] *= con_scaling[prot-1]*p.BAs[pos];
        }
    }
}

void NearNcontributions(bool const * const v,Params const & p, double *express) {
    int i, j, countoverlap;
    int numA = 0, numR = 0, num0 = 0;
    int nbound,protbound,pos,pos1,pos2,prot1,prot2,dist12;
    double con_scaling[3];
    int *distances = new int[p.n-1];
    con_scaling[0] = p.dl*p.Scalings[0];
    con_scaling[1] = p.tw*p.Scalings[1];
    con_scaling[2] = p.sn*p.Scalings[2];
    countoverlap = 0;
    for(i=0; i<p.n-1; i++) if(p.distances[i]==0) {
        distances[i] = p.overlapDistances[countoverlap];
        countoverlap++;
    }
    else distances[i]=p.distances[i];
    //check what is bound 
    for(i=0;i<p.n;i++) {
        //not bound
        if(v[i]==0) num0++;
        //we have an activator bound
        else if(v[i]>0 && p.typeBS[i]<3) numA++;
        //we have a repressor bound
        else if(v[i]>0 && p.typeBS[i]==3 ) numR++;
    }
    nbound = p.n-num0;
    // configuration.proteins is all 0's.....nothing is bound
    if(num0==p.n) {
        express[0] = 0;
        express[1] = 0;
    }
    //if only 1 protein is bound
    else if(num0==p.n-1) {
        pos = 0;
        while(v[pos]==0) pos++;
        //get type of protein bound
        protbound = p.typeBS[pos];
        if(numA == 1) express[0] = con_scaling[protbound-1]*p.BAs[pos];
        else if(numR == 1) express[0] = 0.0;
        express[1] = con_scaling[protbound-1]*p.BAs[pos];
    }
    //more than 1 protein is bound and 1 is an activator
    else if(numA>0 && num0<p.n-1) {
        express[0] = 1.0;
        express[1] = 1.0;
        pos = 0;
        pos1 = 0;
        while(v[pos]==0 && pos < p.n) pos++;
        pos1 = pos;
        //need contribution of scaling and con_scalingration for FIRST
        //protein bound
        prot1 = p.typeBS[pos1];
        express[0] *= con_scaling[prot1-1]*p.BAs[pos1];
        express[1] *= con_scaling[prot1-1]*p.BAs[pos1];
        for(pos2=pos+1;pos2<p.n;pos2++) {
            //find adjacent BOUND proteins
            if(v[pos2]==0) continue;
            //get type of protein2
            prot2 = p.typeBS[pos2];
            //need contribution of scaling and con_scalingration for protein2
            express[0] *= con_scaling[prot2-1]*p.BAs[pos2];
            express[1] *= con_scaling[prot2-1]*p.BAs[pos2];
            dist12 = 0;
            for(j=pos1; j<pos2; j++) dist12 += distances[j];
            //if they are BOTH activators or BOTH repressors, cooperate!
            if((prot1<3 && prot2<3)||(prot1==3 && prot2==3)) {
                express[0] *= coop(prot1,prot2,dist12,p.Cs);
                express[1] *= coop(prot1,prot2,dist12,p.Cs);
            } 
            //if 1 is an activator and the other is a repressor, quench!
            else express[0] *= (1-quench(prot1,prot2,dist12,p.Qs));
            prot1 = prot2;
            pos1 = pos2;
        }
    }
    //no activators bound BUT atleast 2 repressors
    else if(numA==0 && numR>1) {
        express[0] = 0.0;
        express[1] = 1.0;
        pos = 0;
        pos1 = 0;
        while(v[pos]==0) pos++;
        pos1 = pos;
        //need contribution of scaling and con_scalingration for FIRST
        //protein bound
        prot1 = p.typeBS[pos1];
        express[1] *= con_scaling[prot1-1]*p.BAs[pos1];
        for(pos2=pos+1;pos2<p.n;pos2++) {
            //find adjacent BOUND proteins
            if(v[pos2]==0) continue;
            //get type of protein2
            prot2 = p.typeBS[pos2];
            //need contribution of scaling and con_scalingration for protein2
            express[1] *= con_scaling[prot2-1]*p.BAs[pos2];
            dist12 = 0;
            for(j=pos1; j<pos2; j++) dist12 += distances[j];
            //if they are BOTH repressors, cooperate!
            if(prot1==3 && prot2==3) express[1] *= coop(prot1,prot2,dist12,p.Cs);
            prot1 = prot2;
            pos1 = pos2;
        }
    }
    delete []distances;
}

void allNcontributions(bool const * const v,Params const & p, double *express) {
    int i, j, countoverlap;
    int numA = 0, numR = 0, num0 = 0;
    int nbound,protbound,pos,pos1,pos2,prot1,prot2,dist12;
    double con_scaling[3];
    int *distances = new int[p.n-1];
    con_scaling[0] = p.dl*p.Scalings[0];
    con_scaling[1] = p.tw*p.Scalings[1];
    con_scaling[2] = p.sn*p.Scalings[2];
    countoverlap = 0;
    for(i=0; i<p.n-1; i++) {
        if(p.distances[i]==0) {
            distances[i] = p.overlapDistances[countoverlap];
            countoverlap++;
        }
        else distances[i]=p.distances[i];
    }
    //check what is bound 
    for(i=0;i<p.n;i++) {
        //not bound
        if(v[i]==0) num0++;
        //we have an activator bound
        else if(v[i]>0 && p.typeBS[i]<3) numA++;
        //we have a repressor bound
        else if(v[i]>0 && p.typeBS[i]==3) numR++;
    }
    nbound = p.n-num0;
    // v is all 0's.....nothing is bound
    if(num0==p.n) {
        express[0] = 0;
        express[1] = 0;
    }
    //if only 1 protein is bound
    else if(num0==p.n-1) {
        pos = 0;
        while(v[pos]==0) pos++;
        //get type of protein bound
        protbound = p.typeBS[pos];
        if(numA == 1) express[0] = con_scaling[protbound-1]*p.BAs[pos];
        else if(numR == 1) express[0] = 0.0;	
        express[1] = con_scaling[protbound-1]*p.BAs[pos];
    }
    //atleast 1 protein is bound and 1 is an activator
    else if(numA>0 && num0<p.n-1) {
        express[0] = 1.0;
        express[1] = 1.0;
        for(pos1=0;pos1<p.n-1;pos1++) {
            if(v[pos1]==0) continue;
            //need contribution of scaling and con_scalingration for FIRST
            //protein bound
            prot1 = p.typeBS[pos1];
            express[0] *= con_scaling[prot1-1]*p.BAs[pos1];
            express[1] *= con_scaling[prot1-1]*p.BAs[pos1];
            for(pos2=pos1+1;pos2<p.n;pos2++) {
                //find adjacent BOUND proteins
                if(v[pos2]==0) continue;
                //get type of protein2
                prot2 = p.typeBS[pos2];
                //need contribution of scaling and con_scalingration for protein2
                express[0] = express[0]*con_scaling[prot2-1]*p.BAs[pos2];
                express[1] = express[1]*con_scaling[prot2-1]*p.BAs[pos2];
                dist12 = 0;
                for(j=pos1; j<pos2; j++) dist12 = dist12+distances[j];
                //if they are BOTH activators or BOTH repressors, cooperate!
                if((prot1<3 && prot2<3)||(prot1==3 && prot2==3)){
                    express[0] *= coop(prot1,prot2,dist12,p.Cs);
                    express[1] *= coop(prot1,prot2,dist12,p.Cs);
                }
                //if 1 is an activator and the other is a repressor, quench!
                else express[0] = express[0]*(1-quench(prot1,prot2,dist12,p.Qs));
            }
        }
    }
    //no activators bound BUT atleast 1 repressors
    else if(numA==0 && numR>1){
        express[0] = 0.0;
        express[1] = 1.0;
        for(pos1=0;pos1<p.n-1;pos1++) {
            if(v[pos1]==0) continue;
            //need contribution of scaling and con_scalingration for FIRST
            //protein bound
            prot1 = p.typeBS[pos1];
            express[0] *= con_scaling[prot1-1]*p.BAs[pos1];
            express[1] *= con_scaling[prot1-1]*p.BAs[pos1];
            for(pos2=pos1+1;pos2<p.n;pos2++) {
                //find adjacent BOUND proteins
                if(v[pos2]==0) continue;
                //get type of protein2
                prot2 = p.typeBS[pos2];
                //need contribution of scaling and con_scalingration for protein2
                express[1] *= con_scaling[prot2-1]*p.BAs[pos2];
                dist12 = 0;
                for(j=pos1; j<pos2; j++) dist12 = dist12+distances[j];
                //if they are BOTH activators or BOTH repressors, cooperate!
                if((prot1<3 && prot2<3)||(prot1==3 && prot2==3)) express[1] = express[1]*coop(prot1,prot2,dist12,p.Cs);
            }
        }
    }
    delete[] distances;
}

double SSE(double *model,int construct) {
    int i;
    char name[128];
    double E = 0, express[17];
    for (i=0; i<17;i++) {
        express[i] = Expressions[construct-1][i];
        E =  E + (express[i]-model[i])*(express[i]-model[i]);
    }
    return E;
}


double runobjective(double const *x, int dim, int typeC, int typeQ) {
    int i, construct, numC, n, con, j;
    double totalSSE = 0;
    char name[128];
    Params p[17];
    double dl[17], tw[17], sn[17], express[17];
    char c[10];
    char file[80]; 
int constructs[38] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38};;;
	int numcons = 38;
	for(int conindex=0; conindex<numcons; conindex++) {
 con=constructs[conindex];

    for(i=0;i<3;i++)
        thermo_parameters[con-1].Scalings[i] = x[i];
        switch(typeC) {
            case 0:
            case 3:
            case 4:
            case 5:
                numC = 4;
                break;
            case 1:
                numC = 2*thermo_parameters[con-1].Cs[1];
                break;
            case 2:
                numC = 4*thermo_parameters[con-1].Cs[1];
                break;
        }
        for(i=0;i<numC;i++){
            if(typeC ==1 || typeC ==2) thermo_parameters[con-1].Cs[i+3] = x[3+i];
            else  thermo_parameters[con-1].Cs[i+1] = x[3+i];
        }
        if(typeQ==0 || typeQ==3 || typeQ==4 || typeQ==5) {
            thermo_parameters[con-1].Qs[1] = x[3+numC];
            thermo_parameters[con-1].Qs[2] = x[3+numC+1];
        }
        else if(typeQ==1) {
            int numQ = 2*thermo_parameters[con-1].Qs[1];
            for(i=0;i<numQ;i++) thermo_parameters[con-1].Qs[3+i] = x[3+numC+i];
        }
        for (i=0; i<17;i++) {
            p[i].SetParameters(con);
            p[i].SetConcentrations(Concentrations[0][i],Concentrations[1][i],Concentrations[2][i]);
            p[i].saveStates();
            express[i] = p[i].expression();
            delete[] p[i].typeBS;
            delete[] p[i].distances;
            delete[] p[i].BAs;
            delete[] p[i].overlapDistances;
            delete[] p[i].Scalings;
            delete[] p[i].Cs;
            delete[] p[i].Qs;
            int maxstates = pow(2,p[i].n);
            for(j=0;j<maxstates;j++) delete[] p[i].allStates[j];
            delete[] p[i].allStates;
        }
        totalSSE += SSE(express,con);
    }
    return totalSSE;
}