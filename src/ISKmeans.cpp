#include "RCbridge.h"

#include <stdlib.h> 
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

double l2nV(vector<double>& vec){
    double accum = 0.0;
       for (int i = 0; i < vec.size(); ++i) {
           accum += vec[i] * vec[i];
       }
       return sqrt(accum);
}

double l2nDivision(vector<double>& numerator, vector<double>& denominator){
  // calculate 	l2n((numerator+lambda)/(denominator+lambda));	
	double sum2=0;
	double ad;
	for(int i=0;i<numerator.size();i++){
		ad = numerator[i]/denominator[i];
		sum2 += ad*ad;
	}
	return sqrt(sum2);	
}

double l2nDivisionLambda(vector<double>& numerator, vector<double>& denominator, double lambda){
  // calculate 	l2n((numerator+lambda)/(denominator+lambda));	
	double sum2=0;
	double ad;
	for(int i=0;i<numerator.size();i++){
		ad = numerator[i]/(denominator[i]+lambda);
		sum2 += ad*ad;
	}
	return sqrt(sum2);	
}

double BinarySearch_ADMM(vector<double>& numerator, vector<double>& denominator){
  // find a number such that l2n((numerator+lambda)/(denominator+lambda)) = 1;
  if(l2nDivision(numerator,denominator)<=1){
  	return 0.0;
  } 
  double stableE = 1.0/10000; 
  double maxNum = numerator[0];
  double minDen = denominator[0];
  double su;
  
  int n = numerator.size();
  
  for(int i=1;i<n;i++){
	  if(maxNum<numerator[i]){
	  	maxNum = numerator[i];
	  }
	  if(minDen>denominator[i]){
	  	minDen = denominator[i];
	  }	  
  }
  
  double lam1 = 0.0;
  double lam2 = sqrt(n) * maxNum - minDen + stableE;
  int iter = 0;
  while(iter<=20 && (lam2-lam1)>stableE){	  
    su = l2nDivisionLambda(numerator, denominator, (lam1+lam2)/2);
    if(su<1){
      lam2 = (lam1+lam2)/2;
    } else {
      lam1 = (lam1+lam2)/2;
    }
    iter = iter+1;
  }
  return (lam1+lam2)/2;
}

void updateX(double *x, double *y, double *z, int *groupLevel,int *genePos,double *coef, int J, int G, int L, double rho){
	int curStart = 0;
	int curPos;
	double l2na;
	
	// agroupLen: group size of one group
	int agroupLen;
	//double *a;
	
	// a: temp array to store results	
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}
		if(agroupLen==0)		
			cout<<"agroupLen"<<agroupLen<<endl;
  	    std::vector<double> a(agroupLen, 0);
		
		//a = (double*)malloc((agroupLen)*sizeof(double));		
		for(int ap=0; ap<agroupLen; ap++){
			curPos = curStart + ap;
			if(ap>=agroupLen){
				cout<<"ap error 1, "<<ap<<endl;				
			}
			if(curPos>=L){
				cout<<"L error 2, "<<L<<endl;				
				cout<<"curPos error 2, "<<curPos<<endl;				
			}
			if(genePos[curPos]>=J){
				cout<<"genePos[curPos] error 3, "<<genePos[curPos]<<endl;				
			}			
			a[ap] = coef[curPos] * z[genePos[curPos]] - y[curPos]/rho;
		}
		
		l2na = l2nV(a);
		//  cout<<"g:"<<g<<". agroupLen: "<<agroupLen<<endl;				
		// cout<<"g:"<<g<<". l2na: "<<l2na<<endl;				
	    if(l2na*rho<=1){
			for(int bp=curStart; bp<curStart+agroupLen; bp++){
				if(bp>=L){
					cout<<"L error 4, "<<L<<endl;				
					cout<<"bp error 4, "<<bp<<endl;				
				}			
				x[bp] = 0;
			}
		} else {
			for(int ap=0; ap<agroupLen; ap++){
				curPos = curStart + ap;			
				if(curPos>=L){
					cout<<"L error 5, "<<L<<endl;				
					cout<<"curPos error 5, "<<curPos<<endl;				
				}			
				x[curPos] = (1-1/l2na/rho)*a[ap];
			}			
		}
		curStart += agroupLen;
		//free(a);
	}
	
}

void copyZ(vector<double>& z_old,double *z){
	for(int j=0;j<z_old.size();j++){
		z_old[j] = z[j];
	}
}

void updateZ(double *x, double *y, double *z, double *r, int *groupLevel,int *genePos,double *coef, int J, int G, int L, double rho){
	int curStart = 0;
	double am;
	int geneCurPos;
	// agroupLen: group size of one group
	int agroupLen;
	
	// new double [totalLength]()
	// a: temp array to store results
    std::vector<double> b(J, 0);
    std::vector<double> c(J, 0);
    std::vector<double> rplusc(J, 0);
    std::vector<int> nonzero(J, 0);
	
    //double *b = (double*)malloc(J*sizeof(double));
    //double *c = (double*)malloc(J*sizeof(double));
	//double *rplusc = (double*)malloc(J*sizeof(double));
    //int *nonzero = (int*)malloc(J*sizeof(int));

	for(int g=1;g<=G;g++){
		agroupLen = 0;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}		
		if(agroupLen==0)		
			cout<<"agroupLen"<<agroupLen<<endl;
		
		for(int bp=curStart; bp<curStart + agroupLen; bp++){
			am = coef[bp];
			geneCurPos = genePos[bp];
			b[geneCurPos] += am*am*rho;
			c[geneCurPos] += (rho * x[bp] + y[bp])*am;
		}
		curStart += agroupLen;
	}
	
	int countSumTure=0;
	int nonZeroIndex=0;
	for(int j=0;j<J;j++){		
		rplusc[j] = r[j] + c[j];
		//cout<<"j: "<<j<<". rplusc[j]:"<<rplusc[j]<<". r[j]: "<<r[j]<<endl;
		nonzero[j] = (rplusc[j] > 0) ? 1: 0;
		//nonzero[j] = (rplusc[j] > 0 && fabs(rplusc[j])/r[j] > 1.0/1000) ? 1: 0;
		if(nonzero[j]) countSumTure++;		
	}
	/*
	if(countSumTure==0){
		cout<<"not enough penalty for lambda!"<<endl;
	}
	*/

    std::vector<double> nonZeroNum(countSumTure, 0);
    std::vector<double> nonZeroDen(countSumTure, 0);

    //double *nonZeroNum = (double*)malloc(countSumTure*sizeof(double));
    //double *nonZeroDen = (double*)malloc(countSumTure*sizeof(double));
	for(int j=0;j<J;j++){		
		if(nonzero[j]){
			nonZeroNum[nonZeroIndex] = rplusc[j];
			nonZeroDen[nonZeroIndex] = b[j];
			nonZeroIndex++;
		} else {
			z[j] = 0;
		}
	}
	
    double u = BinarySearch_ADMM(nonZeroNum, nonZeroDen);

	nonZeroIndex = 0;
	for(int j=0;j<J;j++){		
		if(nonzero[j]){
			z[j] = nonZeroNum[nonZeroIndex]/(nonZeroDen[nonZeroIndex] + u);
			nonZeroIndex++;			
		}
	}	
	/*
	free(b);
	free(c);
	free(rplusc);
	free(nonzero);
	free(nonZeroNum);
	free(nonZeroDen);	
	*/
}

void updateZbyX(double *x, double *z, int *groupLevel,int *genePos,double *coef, int J, int G){
	int curStart = 0;
	double am;
	int geneCurPos;
	// agroupLen: group size of one group
	int agroupLen;
			
    std::vector<int> markZtoZero(J, 0);
	
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}		
		if(agroupLen==0)		
			cout<<"agroupLen"<<agroupLen<<endl;
		
		for(int bp=curStart; bp<curStart + agroupLen; bp++){
			am = coef[bp];
			geneCurPos = genePos[bp];
			if(am==0 || x[bp]!=0){
				markZtoZero[geneCurPos] = 1;
			} 
		}				
		curStart += agroupLen;
	}
	for(int j=0;j<J;j++){
		if(markZtoZero[j]==0)
			z[j] = 0;		
	}	
}

void updateY(double *x, double *y, double *z, int *groupLevel,int *genePos,double *coef, int G, double rho){
	int curStart = 0;
	// agroupLen: group size of one group
	int agroupLen;

	for(int g=1;g<=G;g++){
		agroupLen = 0;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}
		if(agroupLen==0)		
			cout<<"agroupLen"<<agroupLen<<endl;
				
		for(int bp=curStart; bp<curStart + agroupLen; bp++){
			y[bp] += rho*(x[bp] - coef[bp] * z[genePos[bp]]);
		}	
		curStart += agroupLen;
	}
}

double updateR(double *x, double *z, int *groupLevel,int *genePos,double *coef, int G){
	double sumRP = 0;
	int curStart = 0;
	int curPos;
	// agroupLen: group size of one group
	int agroupLen;
	//double *a;
	
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		//cout<<"groupLevel[curStart+agroupLen]: "<<groupLevel[curStart+agroupLen]<<endl;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}
		if(agroupLen==0)		
			cout<<"agroupLen"<<agroupLen<<endl;
		
  	    std::vector<double> a(agroupLen, 0);
		
		//a = (double*)malloc((agroupLen)*sizeof(double));				
		for(int ap=0; ap<agroupLen; ap++){
			curPos = curStart + ap;
			a[ap] = x[curPos] - coef[curPos] * z[genePos[curPos]];
		}
		sumRP += l2nV(a);
		curStart += agroupLen;	
		//free(a);
	}
	return sumRP;
}

double getObj(double *r,double *z,int J,int *groupLevel,int *genePos,double *coef, int G){
	double obj = 0;
	for(int j=0;j<J;j++){
		obj -= r[j]*z[j];
	}
	
	int curStart = 0;
	int curPos;
	//double *a;
	
	// agroupLen: group size of one group
	int agroupLen;
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		//cout<<"groupLevel[curStart+agroupLen]: "<<groupLevel[curStart+agroupLen]<<endl;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}
		if(agroupLen==0)		
			cout<<"agroupLen"<<agroupLen<<endl;
  	    
		std::vector<double> a(agroupLen, 0);
		
		//a = (double*)malloc((agroupLen)*sizeof(double));		
		for(int ap=0; ap<agroupLen; ap++){
			curPos = curStart + ap;
			a[ap] = coef[curPos] * z[genePos[curPos]];
		}
		
		obj += l2nV(a);
		
		curStart += agroupLen;	
		//free(a);
	}
	return obj;
}

double getRd2(vector<double>& z_old,double *z,int *groupLevel,int *genePos,double *coef, int G, double rho){
	double Rd2=0;
    double a;
	
	int curStart = 0;
	int curPos;
	
	// agroupLen: group size of one group
	int agroupLen;
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		//cout<<"groupLevel[curStart+agroupLen]: "<<groupLevel[curStart+agroupLen]<<endl;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}
		if(agroupLen==0)		
			cout<<"agroupLen"<<agroupLen<<endl;
		
		for(int ap=0; ap<agroupLen; ap++){
			curPos = curStart + ap;
			a = coef[curPos] * (z[genePos[curPos]] - z_old[genePos[curPos]]);
			Rd2 += a*a;
		}
		curStart += agroupLen;		
	}
	
	Rd2 = rho*sqrt(Rd2);
	return Rd2;
}

double getError(vector<double>& z_old,double *z, int J){
	double numerator = 0;
	double denominator = 0;
	for(int j=0;j<J;j++){
		numerator += fabs(z_old[j] - z[j]);
		denominator += fabs(z[j]);
	}
	return numerator/denominator;
}

void ADMM_updatew(double *x, double *y, double *z, double *r, double *objective, int *groupLevel,int *genePos,double *coef, int J, int G, int L){
    bool stopCrit = false;
    int iter = 0;
	double error1 = 1.0/1000;
	int maxiterADMM = 1000;
	
	std::vector<double> z_old(z, z+J);
	
	//double *z_old = (double*)malloc(J*sizeof(double));
	double thisError;
	
	// primal residual sum of square;
	double sumRP;
	// dual residual sum of square;
	double sumDP;

	// parameters for ADMM;
	double rho = 1;
	double incr = 2;
  	double decr = 2;
  	double mu = 10;
	
	while(!stopCrit){
	    iter++;
	    copyZ(z_old,z);
		updateX(x, y, z, groupLevel,genePos,coef, J, G, L, rho);
		updateZ(x, y, z, r, groupLevel, genePos,coef, J, G, L, rho);
		updateY(x,y,z,groupLevel,genePos,coef,G,rho);
		
		sumRP = updateR(x, z, groupLevel,genePos,coef, G);
		*objective = getObj(r,z,J,groupLevel,genePos,coef,G);
		//cout<<"objective: "<<*objective<<endl;
		sumDP = getRd2(z_old,z,groupLevel,genePos,coef,G,rho);
	    thisError = getError(z_old,z,J);
    	
		//cout<<"l2 of primal residual is: "<<sumRP<<endl;
		//cout<<"l2 of dual residual is: "<<sumDP<<endl;
		//cout<<"objective is: "<<*objective<<endl;
		//cout<<"change in z is: "<<thisError<<endl;
		//cout<<"iter: "<<iter<<endl;

		if(sumRP > mu * sumDP) rho = rho * incr;
		if(sumDP > mu * sumRP) rho = rho / decr;								
	    if( (sumRP<error1 && sumDP<error1) || iter>maxiterADMM) stopCrit = true;		
	}
	
	// update z by x
	updateZbyX(x, z, groupLevel, genePos, coef, J, G);
    //free(z_old);
	
}

extern "C" {
	void ADMM_updatew_R(double *x, double *y, double *z, double *r, double *objective, int *groupLevel,int *genePos,double *coef, int *J, int *G, int *L){
		ADMM_updatew(x, y, z, r, objective, groupLevel, genePos, coef, *J, *G, *L);
	}
}

// how to compile
// g++ -c ADMM.cpp -o ADMM.o
// g++ -c -fPIC ADMM.cpp -o ADMM.o
// R CMD SHLIB ADMM.cpp
