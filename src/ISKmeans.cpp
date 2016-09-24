#include "RCbridge.h"

#include <stdlib.h> 
#include <iostream>
#include <math.h>
using namespace std;

double l2n(double *vec, int n){
    double accum = 0.;
       for (int i = 0; i < n; ++i) {
           accum += vec[i] * vec[i];
       }
       return sqrt(accum);
}

double l2nDivision(double *numerator, double *denominator,int n){
  // calculate 	l2n((numerator+lambda)/(denominator+lambda));	
	double sum2=0;
	double ad;
	for(int i=0;i<n;i++){
		ad = numerator[i]/denominator[i];
		sum2 += ad*ad;
	}
	return sqrt(sum2);	
}

double l2nDivisionLambda(double *numerator, double *denominator, double lambda, int n){
  // calculate 	l2n((numerator+lambda)/(denominator+lambda));	
	double sum2=0;
	double ad;
	for(int i=0;i<n;i++){
		ad = numerator[i]/(denominator[i]+lambda);
		sum2 += ad*ad;
	}
	return sqrt(sum2);	
}

double BinarySearch_ADMM(double *numerator, double *denominator,int n){
  // find a number such that l2n((numerator+lambda)/(denominator+lambda)) = 1;
  if(l2nDivision(numerator,denominator,n)<=1){
  	return 0;
  } 
  double stableE = 1.0/10000; 
  double maxNum = numerator[0];
  double minDen = denominator[0];
  double su;
  
  for(int i=1;i<n;i++){
	  if(maxNum<numerator[i]){
	  	maxNum = numerator[i];
	  }
	  if(minDen>denominator[i]){
	  	minDen = denominator[i];
	  }	  
  }
  
  double lam1 = 0;
  double lam2 = sqrt(n) * maxNum - minDen + stableE;
  int iter = 0;
  while(iter<=20 && (lam2-lam1)>stableE){	  
    su = l2nDivisionLambda(numerator, denominator, (lam1+lam2)/2, n);
    if(su<1){
      lam2 = (lam1+lam2)/2;
    } else {
      lam1 = (lam1+lam2)/2;
    }
    iter = iter+1;
  }
  return (lam1+lam2)/2;
}

void updateX(double *x, double *y, double *z, int *groupLevel,int *genePos,double *coef, int G, double rho){
	int curStart = 0;
	int curPos;
	double l2na;
	
	// agroupLen: group size of one group
	int agroupLen;
	
	// a: temp array to store results
    double *a;
	
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}		
		a = (double*)malloc((agroupLen)*sizeof(double));		
		for(int ap=0; ap<agroupLen; ap++){
			curPos = curStart + ap;
			a[ap] = coef[curPos] * z[genePos[curPos]] - y[curPos]/rho;
		}
		
		l2na = l2n(a,agroupLen);
		//  cout<<"g:"<<g<<". agroupLen: "<<agroupLen<<endl;				
		// cout<<"g:"<<g<<". l2na: "<<l2na<<endl;				
	    if(l2na*rho<=1){
			for(int bp=curStart; bp<curStart+agroupLen; bp++){
				x[bp] = 0;
			}
		} else {
			for(int ap=0; ap<agroupLen; ap++){
				curPos = curStart + ap;			
				x[curPos] = (1-1/l2na/rho)*a[ap];
			}			
		}
		curStart += agroupLen;
		free(a);
	}
	
}

void copyZ(double *z_old,double *z, int J){
	for(int j=0;j<J;j++){
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
    double *b = (double*)malloc(J*sizeof(double));
    double *c = (double*)malloc(J*sizeof(double));
	double *rplusc = (double*)malloc(J*sizeof(double));
    int *nonzero = (int*)malloc(J*sizeof(int));
	for(int j=0;j<J;j++){
		nonzero[j] = j;
	}

	for(int j=0;j<J;j++){
		b[j] = 0.0;
		c[j] = 0.0;
	}
	
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}		
		
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
		nonzero[j] = (rplusc[j] > 0 && fabs(rplusc[j])/r[j] > 1.0/1000) ? 1: 0;
		if(nonzero[j]) countSumTure++;		
	}
	if(countSumTure==0){
		cout<<"not enough penalty for lambda!"<<endl;
	}

    double *nonZeroNum = (double*)malloc(countSumTure*sizeof(double));
    double *nonZeroDen = (double*)malloc(countSumTure*sizeof(double));
	for(int j=0;j<J;j++){		
		if(nonzero[j]){
			nonZeroNum[nonZeroIndex] = rplusc[j];
			nonZeroDen[nonZeroIndex] = b[j];
			nonZeroIndex++;
		} else {
			z[j] = 0;
		}
	}
	
    double u = BinarySearch_ADMM(nonZeroNum, nonZeroDen,countSumTure);

	nonZeroIndex = 0;
	for(int j=0;j<J;j++){		
		if(nonzero[j]){
			z[j] = nonZeroNum[nonZeroIndex]/(nonZeroDen[nonZeroIndex] + u);
			nonZeroIndex++;			
		}
	}	
	
	free(b);
	free(c);
	free(rplusc);
	free(nonzero);
	free(nonZeroNum);
	free(nonZeroDen);	
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
				
		for(int bp=curStart; bp<curStart + agroupLen; bp++){
			y[bp] += rho*(x[bp] - coef[bp] * z[genePos[bp]]);
		}	
		curStart += agroupLen;
	}
}

double updateR(double *x, double *z, int *groupLevel,int *genePos,double *coef, int G){
	double sumRP = 0;
	double *a;
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
		
		a = (double*)malloc((agroupLen)*sizeof(double));		
		for(int ap=0; ap<agroupLen; ap++){
			curPos = curStart + ap;
			a[ap] = x[curPos] - coef[curPos] * z[genePos[curPos]];
		}
		sumRP += l2n(a,agroupLen);
		curStart += agroupLen;	
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
    double *a;
	
	// agroupLen: group size of one group
	int agroupLen;
	for(int g=1;g<=G;g++){
		agroupLen = 0;
		//cout<<"groupLevel[curStart+agroupLen]: "<<groupLevel[curStart+agroupLen]<<endl;
		while(groupLevel[curStart+agroupLen] == g){
			agroupLen++;		
		}
		
		a = (double*)malloc((agroupLen)*sizeof(double));		
		for(int ap=0; ap<agroupLen; ap++){
			curPos = curStart + ap;
			a[ap] = coef[curPos] * z[genePos[curPos]];
		}
		
		obj += l2n(a,agroupLen);
		
		curStart += agroupLen;	
		free(a);
	}
	return obj;
}

double getRd2(double *z_old,double *z,int *groupLevel,int *genePos,double *coef, int G, double rho){
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

double getError(double *z_old,double *z, int J){
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
	double *z_old = (double*)malloc(J*sizeof(double));
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
	    copyZ(z_old,z,J);
		updateX(x, y, z, groupLevel,genePos,coef, G, rho);
		updateZ(x, y, z, r, groupLevel, genePos,coef, J, G, L, rho);
		updateY(x,y,z,groupLevel,genePos,coef,G,rho);
		
		sumRP = updateR(x, z, groupLevel,genePos,coef, G);
		*objective = getObj(r,z,J,groupLevel,genePos,coef,G);
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
    free(z_old);
	
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
