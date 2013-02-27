
/* ************************************************************************ *
 * ************************************************************************ *

   File: hmm.C
   The class CHMM defines operations for HMM

  * ************************************************************************ *

   Authors: Daniel DeMenthon & Marc Vuilleumier
   building up on original C code by Tapas Kanungo
   Date:  2-18-99 

 * ************************************************************************ *

   Modification Log:
	4-14-99: Compute log(A) and log(pi) in respective classes
	4-16-99: Compute ViterbiLog with state duration probabilities

 * ************************************************************************ *
   Log for new ideas:
 * ************************************************************************ *
               Language and Media Processing
               Center for Automation Research
               University of Maryland
               College Park, MD  20742
 * ************************************************************************ *
 * ************************************************************************ */
 
//===============================================================================

//===============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>

using namespace std;

#include "utils.h"
#include "read.h"
#include "snp.h"
#include "obs.h"
#include "obsSeq.h"
#include "obsProb.h"
#include "discreteObsProb.h"
#include "gaussianObsProb.h"
#include "vectorObsProb.h"
#include "stateTrans.h"
#include "plainStateTrans.h"
#include "gammaProb.h"
#include "explicitDurationTrans.h"
#include "initStateProb.h"
#include "hmm.h"

extern vector<SNP*> snp_list;
extern vector<READ*> reads_list;
//int DELTA = 0.01;
extern int EM;
extern int flag_iter;
extern int errate;

//===============================================================================

//#define DEBUG
//#define FULLDEBUG

//===============================================================================

//===============================================================================

inline int Factorial(int x) {
  return (x == 1 || x == 0 ? 1 : x * Factorial(x - 1));
}

CHMM::CHMM(CStateTrans *a, CObsProb *b, CInitStateProb *pi)
{
	mA = a;
	mB = b;
	mPi = pi;

	mN = mA->GetN();
}

//===============================================================================

CHMM::~CHMM(void)
{
// Nothing for now
}

//===============================================================================

double CHMM::Forward(double **alpha, double *scale, CObs **obs, long T, boolean doLog)
     // Scaling is used to prevent roundoff errors
     // Same scaling is used for backward and forward procedures
     // so that the scales cancel out in the Baum Welch formula
     // Quantity returned is actually - log(P(O | model),
     // i.e. exponential of this quantity is 1 / P(O | model)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double logProb;
	double bi1, aij, bjt1;

// 1. Initialization
// alpha = f(i)
// bjt1 = emission
	for (i = 1; i <= mN; i++) {
	  bi1 = mB->at(i, 1);
	  alpha[1][i] = mPi->at(i) *  bi1;
	}
	scale[1] = Normalize(alpha[1], mN);

// 2. Induction
	for (t = 1; t <= T - 1; t++) {
	  for (j = 1; j <= mN; j++) {
	    sum = 0.0;
	    for (i = 1; i <= mN; i++){
	      aij = mA->at(i, j);
	      sum += alpha[t][i] * aij;
	    }
	    bjt1 =  mB->at(j, t+1);
	    alpha[t+1][j] = sum * bjt1;
	  }
// scale is the normalization factor s(i)
// alpha is the normalized ~f(i)
	  scale[t+1] = Normalize(alpha[t+1], mN);
	}
	logProb = 0.0;

// 3. Termination
	if(doLog){
	  for (t = 1; t <= T; t++){
	    logProb += log(scale[t]);
	  }
	}// endif

	return logProb;// zero returned if doLog is false
}
//===============================================================================
/*
double CHMM::ForwardAlgo(double **alpha, double *scale, CObs **obs, long T, boolean doLog)
     // Scaling is used to prevent roundoff errors
     // Same scaling is used for backward and forward procedures
     // so that the scales cancel out in the Baum Welch formula
     // Quantity returned is actually - log(P(O | model),
     // i.e. exponential of this quantity is 1 / P(O | model)
{
	int	i, j; 	// state indices 
	int	t;	// time index 

	double sum;	// partial sum 
	double logProb;
	double bi1, aij, bjt1;

// 1. Initialization
// alpha = f(i)
// bjt1 = emission
	for (i = 1; i <= mN; i++) {
// Experimental: Need to confirm the emission for the first read
	  bi1 = mB->at(i, 1);
	  alpha[1][i] = mPi->at(i) *  bi1;
	}
	scale[1] = Normalize(alpha[1], mN);

// 2. Induction
	int nread = 2;
        SNP **reads_snp_list = new SNP*[1000];
        int *index = new int[2000];

	for (t = 1; t <= T - 1; t++) {
		int common_snp_count = 0;
		int known_snp_count = 0;
		bool overlap;
		overlap = TRUE;
                double obslik[1000][2][3], genlik[1000][2][3];
                READ *pd = (((CFlexibleObs<READ*>*)(obs[t]))->Get(1));
                READ *nd = (((CFlexibleObs<READ*>*)(obs[t+1]))->Get(1));
                int rd_start = (*pd).GetPos() < (*nd).GetPos() ? (*nd).GetPos() : (*pd).GetPos();
                int rd_end = (*pd).GetPos()+(*pd).GetLen() > (*nd).GetPos()+(*nd).GetLen() ? (*nd).GetPos()+(*nd).GetLen()-1 : (*pd).GetPos()+(*pd).GetLen()-1;

                GetCommonSnpList(obs, reads_snp_list, &common_snp_count, &known_snp_count, index, t+1);
                if(common_snp_count==0) {
                        cout << "Read " << t << " and " << t+1 << " have no overlapping snps. Skipping to the next pair.." << endl << endl;
                        //exit(1);
                        // (nd)->assignHaplotype(1,0.5);
			overlap = FALSE;
                        // continue;
                }

		for (j = 1; j <= mN; j++) {
			sum = 0.0;
			for (i = 1; i <= mN; i++){
				aij = mA->at(i, j);
				// sum += alpha[nread-1][i] * aij;
				// This is the scaled alpha
				sum += alpha[t][i] * aij;
	    		}

#ifdef DEBUG
cout << "Read " << t << " and " << t+1 << ", " << rd_start << ", " << rd_end << " with " << common_snp_count << " overlapping snps and " << known_snp_count << " common snps." << endl << endl;
cout << "K SnpPos  R A L L\tAA_gen\tAA_obs\tAA_genob\tAB_gen\tAB_obs\tAB_genob\tBB_gen\tBB_obs\tBB_genob\tprob\tlogprob\n";
#endif
			double logBiOt = 0;
			for(int count=0; count<common_snp_count; count++) {
                                SNP *sp = reads_snp_list[count];
#ifdef DEBUG
cout << sp->GetKnown() << " " << sp->GetPos() << " " << sp->GetRef() << " " << sp->GetAlt() << " " << (*pd).GetAllele(index[2*count]) << " " << (*nd).GetAllele(index[2*count+1]);
#endif
                                double emission = compute_new_emission(reads_snp_list, count, obs, t+1, index, j, obslik[count][j-1], genlik[count][j-1]);
                                logBiOt += log(emission)/common_snp_count;
			}
			if(common_snp_count==0) {
				logBiOt = log(0.5);
			}
			bjt1 = pow(2.7182, logBiOt);
			mB->addEmission(j,t+1,bjt1);
	    		bjt1 =  mB->at(j, t+1);
#ifdef DEBUG
cout << "Total emission = " << bjt1 << endl << endl;
#endif
			// Computing emission here
	    		// alpha[nread][j] = sum * bjt1;
	    		// Unscaled yet
	    		alpha[t+1][j] = sum * bjt1;
	  	}
		// scale is the normalization factor s(i)
		// alpha is the normalized ~f(i)
	  	// scale[nread] = Normalize(alpha[nread], mN);
	  	scale[t+1] = Normalize(alpha[t+1], mN);
		nread++;
	}
	nread--;
	//logProb = 0.0;
	// Experimental: logProb contains P(x)
	logProb = 1.0;

// 3. Termination
	if(doLog){
	  for (t = 1; t <= T; t++){
	    // logProb += log(scale[nread-1]);
	    logProb += log(scale[t]);
	  }
	}// endif

	delete [] index;
	delete [] reads_snp_list;
	return logProb;// zero returned if doLog is false
}
*/
double CHMM::ForwardEmission(double **alpha, double *scale, CObs **obs, long T, boolean doLog)
     // Scaling is used to prevent roundoff errors
     // Same scaling is used for backward and forward procedures
     // so that the scales cancel out in the Baum Welch formula
     // Quantity returned is actually - log(P(O | model),
     // i.e. exponential of this quantity is 1 / P(O | model)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double logProb;
	double bi1, aij, bjt1;

// 1. Initialization
// alpha = f(i)
// bjt1 = emission
	for (i = 1; i <= mN; i++) {
	  bi1 = mB->at(i, 1);
	  alpha[1][i] = mPi->at(i) *  bi1;
	}
	scale[1] = Normalize(alpha[1], mN);

// 2. Induction
	int nread = 2;
        SNP **reads_snp_list = new SNP*[1000];
        int *index = new int[2000];

	for (t = 1; t <= T - 1; t++) {
		int common_snp_count = 0;
		bool overlap;
		overlap = TRUE;
                double obslik[1000][2][3], genlik[1000][2][3];
                READ *pd = (((CFlexibleObs<READ*>*)(obs[t]))->Get(1));
                READ *nd = (((CFlexibleObs<READ*>*)(obs[t+1]))->Get(1));
                int rd_start = (*pd).GetPos() < (*nd).GetPos() ? (*nd).GetPos() : (*pd).GetPos();
                int rd_end = (*pd).GetPos()+(*pd).GetLen() > (*nd).GetPos()+(*nd).GetLen() ? (*nd).GetPos()+(*nd).GetLen()-1 : (*pd).GetPos()+(*pd).GetLen()-1;

//                GetCommonSnpList(obs, reads_snp_list, &common_snp_count, index, t+1);
                if(common_snp_count==0) {
                        cout << "Read " << t << " and " << t+1 << " have no overlapping snps. Skipping to the next pair.." << endl << endl;
                        //exit(1);
                        // (nd)->assignHaplotype(1,0.5);
			// overlap = FALSE;
                        // continue;
                }

		for (j = 1; j <= mN; j++) {
			sum = 0.0;
			for (i = 1; i <= mN; i++){
				aij = mA->at(i, j);
				// sum += alpha[nread-1][i] * aij;
				// This is the scaled alpha
				sum += alpha[t][i] * aij;
	    		}
	if(sum!=0.5) {
		//printf("Sum = %d\t%f\n",t,sum);
	}

#ifdef FULLDEBUG
cout << "Read " << t << " and " << t+1 << ", " << rd_start << ", " << rd_end << " with " << common_snp_count << " overlapping snps" << endl << endl;
cout << "K SnpPos  R A L L\tAA_gen\tAA_obs\tAA_genob\tAB_gen\tAB_obs\tAB_genob\tBB_gen\tBB_obs\tBB_genob\tprob\tlogprob\n";
#endif
			double logBiOt = 0;
			for(int count=0; count<common_snp_count; count++) {
                                SNP *sp = reads_snp_list[count];
#ifdef FULLDEBUG
cout << sp->GetKnown() << " " << sp->GetPos() << " " << sp->GetRef() << " " << sp->GetAlt() << " " << (*pd).GetAllele(index[2*count]) << " " << (*nd).GetAllele(index[2*count+1]);
#endif
                                double emission = compute_new_emission(reads_snp_list, count, obs, t+1, index, j, obslik[count][j-1], genlik[count][j-1]);
                                //logBiOt += emission;
                                logBiOt += log(emission)/common_snp_count;
			}
			if(common_snp_count==0) {
				//logBiOt = 0.5;
				//sum = 0.5;
				logBiOt = log(0.5);
			} else {
				//bjt1 = logBiOt/common_snp_count;
			}
			bjt1 = pow(2.7182, logBiOt);
			// RECHECK: Do I update the emissions here or not?
			//mB->addEmission(j,t+1,bjt1);
	    		//bjt1 =  mB->at(j, t+1);
#ifdef DEBUG
cout << "Total emission = " << bjt1 << endl << endl;
#endif
			// Computing emission here
	    		// alpha[nread][j] = sum * bjt1;
	    		// Unscaled yet
	    		alpha[t+1][j] = sum * bjt1;
	  	}
		// scale is the normalization factor s(i)
		// alpha is the normalized ~f(i)
	  	// scale[nread] = Normalize(alpha[nread], mN);
	  	scale[t+1] = Normalize(alpha[t+1], mN);
		nread++;
	}
	nread--;

	logProb = 0.0;

// 3. Termination
	if(doLog){
	  for (t = 1; t <= T; t++){
	    // logProb += log(scale[nread-1]);
	    logProb += log(scale[t]);
	  }
	}// endif

	delete [] index;
	delete [] reads_snp_list;
	return logProb;// zero returned if doLog is false
}

double CHMM::ForwardLogEmission(double **alpha, double *scale, CObs **obs, long T, boolean doLog)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double logProb;
	double bi1, aij, bjt1;

// 1. Initialization
// alpha = f(i)
// bjt1 = emission
	for (i = 1; i <= mN; i++) {
	  bi1 = mB->at(i, 1);
	  if(bi1==0.0) {
		cout << "bi1 = 0... exiting.\n";
		exit(1);
	  }
	  //if(bi1==0.0) bi1 = pow(2.7183, -10.0);
	  alpha[1][i] = log(bi1);
	  //alpha[1][i] = log(mPi->at(i)) + log(bi1);
	}

// 2. Induction
	int nread = 2;
        SNP **reads_snp_list = new SNP*[2000];
        SNP **known_snp_list = new SNP*[100];
        int *common_index = new int[4000];
        int *known_index = new int[200];

	for (t = 1; t <= T - 1; t++) {
		int common_snp_count = 0;
		int known_snp_count = 0;
		bool overlap;
                double obslik[2000][2][3], genlik[2000][2][3];
                READ *pd = (((CFlexibleObs<READ*>*)(obs[t]))->Get(1));
                READ *nd = (((CFlexibleObs<READ*>*)(obs[t+1]))->Get(1));
                int rd_start = (*pd).GetPos() < (*nd).GetPos() ? (*nd).GetPos() : (*pd).GetPos();
                int rd_end = (*pd).GetPos()+(*pd).GetLen() > (*nd).GetPos()+(*nd).GetLen() ? (*nd).GetPos()+(*nd).GetLen()-1 : (*pd).GetPos()+(*pd).GetLen()-1;

                GetCommonSnpList(obs, reads_snp_list, known_snp_list, &common_snp_count, &known_snp_count, common_index, known_index, t+1);
		(*nd).AddKnownCount(known_snp_count);
#ifdef DEBUG
if(known_snp_count>0)
cout << "For read " << nd->GetPos() << " added known overlapping count of " << known_snp_count << endl;
#endif

		for (j = 1; j <= mN; j++) {
			sum = 0.0;
			for (i = 1; i <= mN; i++) {
				aij = mA->at(i, j);
				double aij1;
				if(i==1) {
					sum += (alpha[t][i] + log(aij));
					aij1 = aij;
				} else if(i==2) {
					double sumi = pow(2.7183, alpha[t][i] - alpha[t][i-1] + log(aij) - log(aij1));
					sum += log(1 + sumi);
				} else {
					cout << "Unexpected number of states (>2).\n";
					exit(1);
				}
	    		}
#ifdef FULLDEBUG
cout << "Read " << t << " and " << t+1 << ", " << rd_start << ", " << rd_end << " with " << common_snp_count << " overlapping snps and " << known_snp_count << " common snps." << endl << endl;
cout << "K SnpPos  R A L L\tAA_gen AA_obs AA_genob\tAB_gen AB_obs AB_genob\tBB_gen BB_obs BB_genob\tprob\tlogprob\n";
#endif
			double logBiOt = 0.0;
			if(flag_iter) {
				for(int count=0; count<common_snp_count; count++) {
                        	        SNP *sp = reads_snp_list[count];
#ifdef FULLDEBUG
cout << sp->GetKnown() << " " << sp->GetPos() << " " << sp->GetRef() << " " << sp->GetAlt() << " " << (*pd).GetAllele(common_index[2*count]) << " " << (*nd).GetAllele(common_index[2*count+1]);
#endif
					// j here is (hap1==hap2), not the value of hap itself.
					// absolute value of state is decided after backward algorithm is complete,
					// in the calculation of gamma values in the calling function.
					double emission = compute_new_emission(reads_snp_list, count, obs, t+1, common_index, j, obslik[count][j-1], genlik[count][j-1]);
					if(sp->GetKnown()==1)
                        	        	//logBiOt += log(emission);
                        	        	logBiOt += (emission*1000.0);
					else
                        	        	//logBiOt += log(emission)/1000.0;
                        	        	logBiOt += (emission);
				}
				if(common_snp_count==0) {
					//logBiOt = log(1.0/T);
					logBiOt = (1.0/T);
					bjt1 = logBiOt;
                        		cout << "Read " << t << " and " << t+1 << " have no overlapping snps. Skipping to the next pair.." << endl << endl;
				} else {
					// Ignoring nth root.
					bjt1 = logBiOt;
					bjt1 = logBiOt/(double)(common_snp_count);
					//bjt1 = logBiOt/(double)(common_snp_count+(999*known_snp_count));
				}
			} else {
				for(int count=0; count<known_snp_count; count++) {
                                	SNP *sp = known_snp_list[count];
#ifdef FULLDEBUG
cout << sp->GetKnown() << " " << sp->GetPos() << " " << sp->GetRef() << " " << sp->GetAlt() << " " << (*pd).GetAllele(known_index[2*count]) << " " << (*nd).GetAllele(known_index[2*count+1]);
#endif
					// j here is (hap1==hap2), not the value of hap itself.
					// absolute value of state is decided after backward algorithm is complete,
					// in the calculation of gamma values in the calling function.
					double emission = compute_new_emission(known_snp_list, count, obs, t+1, known_index, j, obslik[count][j-1], genlik[count][j-1]);
                                	logBiOt += log(emission);
				}
				if(known_snp_count==0) {
					logBiOt = log(1.0/T);
					bjt1 = logBiOt;
                        		cout << "Read " << t << " and " << t+1 << " have no overlapping snps. Skipping to the next pair.." << endl << endl;
				} else {
					// Ignoring nth root.
					bjt1 = logBiOt;
					bjt1 = logBiOt/known_snp_count;
				}
			}
#ifdef DEBUG
cout << "Total emission = " << bjt1 << endl << endl;
#endif
			// Computing emission here
			// REVISIT: Still need to ensure whether emissions need to be updated here for 
			// use in Backward algorithm
			//mB->addEmission(j,t+1,pow(2.7183,bjt1));
			mB->addEmission(j,t+1,bjt1);
	    		//alpha[t+1][j] = sum + (bjt1);
	    		alpha[t+1][j] = sum + log(bjt1);
	  	}
		nread++;
	}
	nread--;
	// Experimental: logProb contains P(x)
	logProb = 0.0;
	flag_iter=1;

// 3. Termination

	double full_prob = alpha[T][1] + log(0.5) + log(1.0+pow(2.7183,alpha[T][2]-alpha[T][1]));
	if(doLog){
	  for (t = 1; t <= T; t++){
	    // logProb += log(scale[nread-1]);
	    logProb += log(scale[t]);
	  }
	}// endif

	delete [] common_index;
	delete [] known_index;
	delete [] reads_snp_list;
	delete [] known_snp_list;
	//return logProb;// zero returned if doLog is false
	return full_prob;
}

//===============================================================================

void CHMM::Backward(double **beta, double *scale, CObs **obs, long T)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
	double sum;
 
 
// 1. Initialization
		for (i = 1; i <= mN; i++){
			beta[T][i] = 1.0/scale[T];
		}
 
// 2. Induction
// beta = b(i)
     for (t = T - 1; t >= 1; t--){
	  for (i = 1; i <= mN; i++){
	    sum = 0.0;
	    for (j = 1; j <= mN; j++){
	      sum += mA->at(i,j) * mB->at(j, t+1) * beta[t+1][j];
	    }
// beta is the normalized ~b(i)
	    beta[t][i] = sum/scale[t];
	    }
     }
}
//===============================================================================

void CHMM::BackwardAlgo(double **beta, double *scale, CObs **obs, long T)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
	double sum;
 
 
// 1. Initialization
		for (i = 1; i <= mN; i++){
			beta[T][i] = 1.0/scale[T];
		}
 
// 2. Induction
// beta = b(i)
     for (t = T - 1; t >= 1; t--){
	  for (i = 1; i <= mN; i++){
	    sum = 0.0;
	    for (j = 1; j <= mN; j++){
	      sum += mA->at(i,j) * mB->at(j, t+1) * beta[t+1][j];
	      // sum += mA->at(i,j) * mB->at(j, obs[t+1]) * beta[t+1][j];
	    }
// beta is the normalized ~b(i)
	    beta[t][i] = sum/scale[t];
	    }
     }
}

void CHMM::BackwardEmission(double **beta, double *scale, CObs **obs, long T)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
	double sum;
 
 
// 1. Initialization
		for (i = 1; i <= mN; i++){
			beta[T][i] = 1.0/scale[T];
		}
 
// 2. Induction
// beta = b(i)
// RECHECK: Verify beta assignments here, or before getting here
     for (t = T - 1; t >= 1; t--){
	  for (i = 1; i <= mN; i++){
	    sum = 0.0;
	    for (j = 1; j <= mN; j++){
		double trans = mA->at(i,j);
		double em = mB->at(j,t+1);
		double bet = beta[t+1][j];
	      sum += trans * em * bet;
	      // sum += mA->at(i,j) * mB->at(j, obs[t+1]) * beta[t+1][j];
	    }
// beta is the normalized ~b(i)
	    beta[t][i] = sum/scale[t];
	    }
     }
}

double CHMM::BackwardLogEmission(double **beta, double *scale, CObs **obs, long T)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
	double sum;
 
 
// 1. Initialization
		for (i = 1; i <= mN; i++){
			beta[T][i] = log(mA->at(i,1));
		}
 
// 2. Induction
// beta = b(i)
     for (t = T - 1; t >= 1; t--) {
	  for (i = 1; i <= mN; i++) {
	    sum = 0.0;
	    for (j = 1; j <= mN; j++) {
		double trans = mA->at(i,j);
		double em = mB->at(j,t+1);
		double bet = beta[t+1][j];
		if(j==1)
			sum += (log(trans) + log(em) + bet);
		else if(j==2) {
			double sumi = pow(2.7183,log(mB->at(j,t+1))-log(mB->at(j-1,t+1))+beta[t+1][j]-beta[t+1][j-1]);
	      		sum += log(1.0 + sumi);
		} else {
			cout << "Unexpected number of states (>2).\n";
			exit(1);
		}
	    }
// beta is the normalized ~b(i)
	    beta[t][i] = sum;
	    }
     }
    double bi1 = mB->at(1, 1);
    double bi2 = mB->at(2, 1);

    //if(bi1==0) bi1 = pow(2.7183, -10.0);
    //if(bi2==0) bi2 = pow(2.7183, -10.0);
    if(bi1==0||bi2==0) {
	cout << "bi2 = 0. Exiting..\n";
	exit(1);
    }

    return log(mA->at(1,2))+log(bi2)+beta[1][2] + log(1 + pow(2.7183,log(bi2)-log(bi1)+beta[1][2]-beta[1][1]));
}

//===============================================================================

double CHMM::Viterbi(CObs **obs, long T, int *q)
// returns sequence q of most probable states
// and probability of seeing that sequence
// if A, B, pi are already given
{
	int 	i, j;	// state indices
	int  	t;	// time index

	int	argmaxval;
	double	maxval, val;
	double pprob;
	double *prevDelta, *delta, *tmp;
	int **psi;
	
	delta = SetVector(mN);
	prevDelta = SetVector(mN);
	psi = SetIntMatrix(T, mN);

//1. Initialization
	for (i = 1; i <= mN; i++) {
		prevDelta[i] = mPi->at(i) * mB->at(i, obs[1]);
		psi[1][i] = 0;
	}	

// 2. Recursion
	for (t = 2; t <= T; t++) {
		for (j = 1; j <= mN; j++) {
			maxval = 0.0;
			argmaxval = 1;	
			for (i = 1; i <= mN; i++) {
				val = prevDelta[i] * mA->at(i, j);
				if (val > maxval) {
					maxval = val;	
					argmaxval = i;	
				}
			}
			delta[j] = maxval * mB->at(j, obs[t]);
			psi[t][j] = argmaxval; // reverse order of indices?
		}
		tmp = delta; delta = prevDelta; prevDelta = tmp;
	}

// 3. Termination
	pprob = 0.0;
	q[T] = 1;
	for (i = 1; i <= mN; i++) {
	  if (prevDelta[i] > pprob) {
			pprob = prevDelta[i];	
			q[T] = i;
	  }
	}

// 4. Path (state sequence) backtracking

	for (t = T - 1; t >= 1; t--){
		q[t] = psi[t+1][q[t+1]];
	}
	delete [] psi[1];
	delete [] psi;
	delete [] prevDelta;
	delete delta;
	
	return pprob;
}

//===============================================================================
double CHMM::ViterbiLog(CObs **obs, long T, int *q)
// returns sequence q of most probable states
// and probability of seeing that sequence
// if A, B, pi are already given
// This implementation uses logarithms to avoid underflows.
{

	int 	i, j;	// state indices
	int  	t;	// time index

	int	argmaxval;
	double	maxval, val, bVal, logVal;
	double logProb;
	double *prevDelta, *delta, *tmp;
	int **psi;
	double **logBiOt;
	int zeroProbCount;

/*
	delta = SetVector(mN);
	prevDelta = SetVector(mN);
	psi = SetIntMatrix(T, mN);

// We do not preprocess the logs for B
// because the data have variable length T
	logBiOt =  SetMatrix(mN, T);
	for (t = 1; t <= 1; t++){
            zeroProbCount = 0;
            for (i = 1; i <= mN; i++){ 
              bVal = mB->at(i, t); // obs[t] is an entire seq of snps in a read
                if(bVal<=0.0){
                    logVal = -1000.0;
                    zeroProbCount++;
                }
                else{
                    logVal = log(bVal);
                }
                logBiOt[i][t] = logVal;
            }// for i
            if(zeroProbCount == mN){// unseen obs, Viterbi decides which seen obs to use
                cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                for (i = 1; i <= mN; i++){		  
                    logBiOt[i][t] = log(1.0/mN);
                }
            }
 	}// for t
	// May have to move this closing brace to the end or rather
	// perform this log conversion in the recursion step.
// 1. Initialization
	for (i = 1; i <= mN; i++){
		prevDelta[i] = mPi->logAt(i) + logBiOt[i][1];
		psi[1][i] = 0; // What's this for?
		mA->InitViterbiDurations(i); // NO-OP
	}
 
// 2. Recursion
	for (t = 2; t <= T; t++) { // successive reads
		zeroProbCount = 0;
		int common_snp_count = 0;
		for (j = 1; j <= mN; j++) {
			maxval = prevDelta[1] + mA->logAt(1, j);
			argmaxval = 1;
			for (i = 1; i <= mN; i++) {// previous state
				val = prevDelta[i] + mA->logAt(i, j);
				if (val > maxval) {
					maxval = val;
					argmaxval = i;
				}
			}

		SNP **reads_snp_list = new SNP*[250];
		int *index = new int[500];
			common_snp_count = 0;
			int hap = (logBiOt[j][t-1] > logBiOt[(j%mN)+1][t-1]) ? 1 : 2; // 1 -> H1==H2, 2 -> H1!=H2

			GetCommonSnpList(obs, reads_snp_list, &common_snp_count, index, t);

			// I end up not assigning the actual emission matrix, but only the log matrix for the computations

			logBiOt[j][t] = 1.0;
			for(int count=0; count<common_snp_count; count++) {
				double emission = compute_new_emission(reads_snp_list, count, obs, t, index, hap);
				if(emission <= 0.0) {
					logVal = -1000.0;
                    			zeroProbCount++;
                		} else{
                    			logVal = log(emission);
                		}
                		logBiOt[j][t] *= logVal;
			}

			delta[j] = maxval + logBiOt[j][t]; 
			psi[t][j] = argmaxval; // What's this for?
			mA->UpdateViterbiDurations(argmaxval, j);// ***dfd 4-16-99
		}
		if(zeroProbCount > common_snp_count) { // need to check per snp
                	cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                	for (i = 1; i <= mN; i++){		  
                    		logBiOt[i][t] = log(1.0/mN);
                	}
		}
		tmp = delta; delta = prevDelta; prevDelta = tmp;
	}
 
// 3. Termination
	logProb = prevDelta[1];
	q[T] = 1;
	for (i = 1; i <= mN; i++) {
		if (prevDelta[i] > logProb) {
		      logProb = prevDelta[i];
		      q[T] = i;
		}
	} 
 
// 4. Path (state sequence) backtracking
	for (t = T - 1; t >= 1; t--){
		q[t] = psi[t+1][q[t+1]];
	}
	
	delete [] psi[1];
	delete [] psi;
	delete [] prevDelta;
	delete [] delta;
	delete [] logBiOt[1];
	delete [] logBiOt;
#if 0
	delete [] logA[1];
	delete [] logA;
	delete [] logPi;
#endif	
*/
	return logProb;
}
double CHMM::ViterbiLog(CObs **obs, long T, int *q, double *probarray)
// returns sequence q of most probable states
// and probability of seeing that sequence
// if A, B, pi are already given
// This implementation uses logarithms to avoid underflows.
{

	int 	i, j;	// state indices
	int  	t;	// time index

	int	argmaxval;
	double  prob;
	double	maxval, val, bVal, logVal;
	double logProb;
	double *prevDelta, *delta, *tmp;
	int **psi;
	double **logBiOt;
	int zeroProbCount;

	delta = SetVector(mN);
	prevDelta = SetVector(mN);
	psi = SetIntMatrix(T, mN);

// We do not preprocess the logs for B
// because the data have variable length T
	logBiOt =  SetMatrix(mN, T);
	for (t = 1; t <= 1; t++){
            zeroProbCount = 0;
            for (i = 1; i <= mN; i++){ 
              bVal = mB->at(i, t); // obs[t] is an entire seq of snps in a read
                if(bVal<=0.0){
                    logVal = -10.0;
                    zeroProbCount++;
                }
                else{
                    logVal = log(bVal);
                }
                logBiOt[i][t] = logVal;
            }// for i
            if(zeroProbCount == mN){// unseen obs, Viterbi decides which seen obs to use
                cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                for (i = 1; i <= mN; i++){		  
                    logBiOt[i][t] = log(1.0/mN);
                }
            }
 	}// for t

// 1. Initialization
// Initialization is performed for prevDelta and psi only.
// PrevDelta stores probability to stage k-1 for all states

	for (i = 1; i <= mN; i++){
		prevDelta[i] = mPi->logAt(i) + logBiOt[i][1];
		psi[1][i] = 0; // What's this for?
		mA->InitViterbiDurations(i); // NO-OP
	}
 
// 2. Recursion
	int nread = 2;
	SNP **reads_snp_list = new SNP*[2000];
	int *index = new int[2000];
	for (t = 2; t <= T; t++) { // successive reads
		zeroProbCount = 0;
		int common_snp_count = 0;
		double obslik[2000][2][3], genlik[2000][2][3];

		READ *pd = (((CFlexibleObs<READ*>*)(obs[t-1]))->Get(1));
		READ *nd = (((CFlexibleObs<READ*>*)(obs[t]))->Get(1));
		int rd_start = (*pd).GetPos() < (*nd).GetPos() ? (*nd).GetPos() : (*pd).GetPos();
		int rd_end = (*pd).GetPos()+(*pd).GetLen() > (*nd).GetPos()+(*nd).GetLen() ? (*nd).GetPos()+(*nd).GetLen()-1 : (*pd).GetPos()+(*pd).GetLen()-1;

//		GetCommonSnpList(obs, reads_snp_list, &common_snp_count, index, t);

		//Ignore this for now: RECHECK: There is a bug here. In case two consecutive reads with 0 overlapping snps are encountered, it may not
		//continue to work as expected. We might have to reset haplotype assumptions, or compare the last read with overlapping
		//snps to the next such one (again to which there might be very little chance)
		if(common_snp_count==0) {
			cout << "Read " << t-1 << " and " << t << " have no overlapping snps. Skipping to the next pair.." << endl << endl;
			//exit(1);
			(nd)->assignHaplotype(1,0.5);
			continue;
		}
		for (j = 1; j <= mN; j++) {
			// prevDelta[i] = v(k)(i)
			maxval = prevDelta[1] + mA->logAt(1, j);
			argmaxval = 1;
			prob = maxval;
			for (i = 1; i <= mN; i++) {// previous state
				val = prevDelta[i] + mA->logAt(i, j);
				if (val > maxval) {
					maxval = val;
					argmaxval = i;
					prob = maxval;
				}
			}
			// maxval = max(k) (v(k)(i)a(k)(l))
			// need to add emission here and find max state for next read

			// I end up not assigning the actual emission matrix, but only the log matrix for the computations

			logBiOt[j][nread] = 0.0;
			int hap = prevDelta[j] > prevDelta[(j%mN)+1] ? 1 : (prevDelta[j] < prevDelta[(j%mN)+1] ? 2 : j) ;
			int abs_hap = prevDelta[j] >= prevDelta[(j%mN)+1] ? j : (j%mN)+1;

#ifdef FULLDEBUG
cout << "Read " << t-1 << " and " << t << ", Hap " << hap << ", Abs Hap = " << abs_hap << ", " << rd_start << ", " << rd_end << " with " << common_snp_count << " overlapping snps" << endl << endl;
cout << "K SnpPos  R A L L\tAA_gen\tAA_obs\tAA_genob\tAB_gen\tAB_obs\tAB_genob\tBB_gen\tBB_obs\tBB_genob\tprob\tlogprob\n";
#endif
			for(int count=0; count<common_snp_count; count++) {
				SNP *sp = reads_snp_list[count];
//				if(sp->GetKnown()==1) {
#ifdef FULLDEBUG
cout << sp->GetKnown() << " " << sp->GetPos() << " " << sp->GetRef() << " " << sp->GetAlt() << "\t" << (*pd).GetAllele(index[2*count]) << " " << (*nd).GetAllele(index[2*count+1]);
#endif
				double emission = compute_new_emission(reads_snp_list, count, obs, t, index, hap, obslik[count][j-1], genlik[count][j-1]);
				if(emission <= 0.0) {
					logVal = -100.0;
                    			zeroProbCount++;
                		} else{
                    			logVal = log(emission);
                		}
#ifdef FULLDEBUG
cout << "\t" << emission << "\t" << logVal << endl;
#endif
				// Experimental: Do I divide by common_snp_count or not?
                			logBiOt[j][nread] += logVal/common_snp_count;
                			//logBiOt[j][nread] += logVal;
//				}
			}
			// logBiOt[j][nread] now contains the mN different emissions to be added to maxval to determine the max state at this stage
#ifdef FULLDEBUG
cout << "Total emission = " << logBiOt[j][nread] << "," << maxval << endl << endl;
#endif

			delta[j] = maxval + logBiOt[j][nread];
			psi[nread][j] = argmaxval; // What's this for?
			mA->UpdateViterbiDurations(argmaxval, j);// ***dfd 4-16-99
		}
		// By now all v(l)(i+1)s are computed from which the max has to be found
		if(zeroProbCount > common_snp_count) { // need to check per snp
                	cerr << "*** Unseen data, renormalizing logBiOt to equiprobs ***"<<endl;
                	for (i = 1; i <= mN; i++){		  
                    		logBiOt[i][nread] = log(1.0/mN);
                	}
		}

		// Updating posteriors here
		int ind_snp = 0;
		double hap_prob = 0.0;
		ind_snp = delta[1] > delta[2] ? 1 : 2;
		// RECHECK: hap_prob assignment seems incorrect
		//hap_prob = pow(2.7182, -(logBiOt[ind_snp][nread] - logBiOt[ind_snp%mN + 1][nread]) );
		double prob1 = pow(2.7182, logBiOt[ind_snp][nread]);
		double prob2 = pow(2.7182, logBiOt[ind_snp%mN+1][nread]);
		hap_prob = prob1/(prob1+prob2);
#ifdef FULLDEBUG
cout << "leading emission " << prob1 << ", trailing emission " << prob2 << endl;
cout << "Assigning read " << (*nd).GetPos() << " haplotype " << ind_snp << " and happrob " << hap_prob << endl << endl << endl;
#endif
		(nd)->assignHaplotype(ind_snp, hap_prob);
		nread++;
		// prevDelta updated to find the max from in the next read
		tmp = delta; delta = prevDelta; prevDelta = tmp;
	}
	nread--;
 
// 3. Termination
	logProb = prevDelta[1];
	//q[T] = 1;
	q[nread] = 1;
	for (i = 1; i <= mN; i++) {
		if (prevDelta[i] > logProb) {
		      logProb = prevDelta[i];
		      q[nread] = i;
		      //q[T] = i;
		}
	} 
 
// 4. Path (state sequence) backtracking
	//for (t = T - 1; t >= 1; t--){
	for (t = nread - 1; t >= 1; t--){
		q[t] = psi[t+1][q[t+1]];
	}
	
	delete [] psi[1];
	delete [] psi;
	delete [] prevDelta;
	delete [] delta;
	delete [] logBiOt[1];
	delete [] logBiOt;
#if 0
	delete [] logA[1];
	delete [] logA;
	delete [] logPi;
#endif	
	return logProb;
}

void CHMM::UpdateBaumWelchGenotypePosteriors(long snp_start, long snp_end)
{
	for(vector<SNP*>::iterator snp_it = (snp_list).begin(); snp_it != (snp_list).end(); snp_it++) {
		long pos = (*snp_it)->GetPos();
		if(pos<=snp_start)
			continue;
		else if(pos>snp_end)
			break;

		int g=0;
		int genotype = 0;
		double post[3];
		double prior[3];
		double *happrob = new double[3];
		double *genprob = new double[2];
		double genp = 1.0;
		double norm = 0.0;

		int ref_ct = (*snp_it)->GetRefCount();
		int alt_ct = (*snp_it)->GetAltCount();
		int err_ct = (*snp_it)->GetErrCount();

		prior[0] = (*snp_it)->GetPosteriors()[0];
		prior[1] = (*snp_it)->GetPosteriors()[1];
		prior[2] = (*snp_it)->GetPosteriors()[2];

		happrob[0] = 0.33;
		happrob[1] = 0.33;
		happrob[2] = 0.33;
		haplotypeProbability(snp_it, happrob);

		for(g=0; g<3; g++) {
			norm += prior[g]*happrob[g];
		}

		for(g=0; g<3; g++) {
			post[g] = prior[g]*happrob[g]/norm;
		}
#ifdef DEBUG
cout << "Posterior het prob for snp " << (*snp_it)->GetPos() << " = " << post[1] << endl;
#endif
    		(*snp_it)->add_posteriors(post);

		delete [] happrob;
		delete [] genprob;
	}
}

void CHMM::UpdateGenotypes(long snp_start, long snp_end)
{
	for(vector<SNP*>::iterator snp_it = (snp_list).begin(); snp_it != (snp_list).end(); snp_it++) {
		long pos = (*snp_it)->GetPos();
		if(pos<=snp_start)
			continue;
		else if(pos>snp_end)
			break;

		int g=0;
		int genotype = 0;
		double post[3];
		double prior[3];
		double *happrob = new double[3];
		double *genprob = new double[2];
		double genp = 1.0;
		double norm = 0.0;

		int ref_ct = (*snp_it)->GetRefCount();
		int alt_ct = (*snp_it)->GetAltCount();
		int err_ct = (*snp_it)->GetErrCount();

		// REVISIT: Was using genlik, replacing with posterior
		// should be using the latest as in updatebaumwelchgenotype..
		prior[0] = (*snp_it)->GetGenLik()[0];
		prior[1] = (*snp_it)->GetGenLik()[1];
		prior[2] = (*snp_it)->GetGenLik()[2];
		// prior[0] = (*snp_it)->GetPosteriors()[0];
		// prior[1] = (*snp_it)->GetPosteriors()[1];
		// prior[2] = (*snp_it)->GetPosteriors()[2];
		haplotypeProbability(snp_it, happrob);

#ifdef DEBUG
cout << "Returned happrob =" << "\t" << happrob[0] << "\t" << happrob[1] << "\t" << happrob[2] << endl;
#endif

		for(g=0; g<3; g++) {
			norm += prior[g]*happrob[g];
		}

		for(g=0; g<3; g++) {
			post[g] = prior[g]*happrob[g]/norm;
		}
		// post[0] = (*snp_it)->GetPosteriors()[0];
		// post[1] = (*snp_it)->GetPosteriors()[1];
		// post[2] = (*snp_it)->GetPosteriors()[2];
    		(*snp_it)->add_posteriors(post);

		if(post[1]>post[0]&&post[1]>post[2]) {
			genotypeProbability(snp_it, genprob);
			if(genprob[0]>genprob[1]) {
				genotype = 2;
				genp = genprob[0];
			} else if(genprob[1]>genprob[0]) {
				genotype = 3;
				genp = genprob[1];
			} else {
				cout << "Khao Britannia 50-50. Is it very tasty tasty?" << endl;
			}
		} else if(post[0]>post[1]&&post[0]>post[2]) {
			genotype = 1;
		} else if(post[2]>post[0]&&post[2]>post[1]) {
			genotype = 4;
		} else {
			cout << "Go learn your math from a high school kid and come back here" << endl;
		}
		(*snp_it)->assign_genotype(genotype, genp);

#ifdef DEBUG
		cout << "SNP: " << (*snp_it)->GetPos() << " " << (*snp_it)->GetKnown() << endl;
		cout << ref_ct << "\t" << alt_ct << "\t" << err_ct << "\t" << (*snp_it)->GetPos() << endl;

		cout << "Priors:\t";
		cout << ref_ct << "\t" << alt_ct << "\t" << err_ct << "\t";
		cout << "SNP:" << (*snp_it)->GetPos() << "\t" << (*snp_it)->GetKnown() << "\t";
		cout << prior[0] << "\t" << prior[1] << "\t" << prior[2] << endl;

		cout << "Happrobs:\t";
		cout << ref_ct << "\t" << alt_ct << "\t" << err_ct << "\t";
		cout << "SNP:" << (*snp_it)->GetPos() << "\t" << (*snp_it)->GetKnown() << "\t";
		cout << happrob[0] << "\t" << happrob[1] << "\t" << happrob[2] << endl;

		cout << "Posteriors:\t";
		cout << ref_ct << "\t" << alt_ct << "\t" << err_ct << "\t";
		cout << "SNP:" << (*snp_it)->GetPos() << "\t" << (*snp_it)->GetKnown() << "\t";
		cout << post[0] << "\t" << post[1] << "\t" << post[2] << endl;

		cout << "Genotype:\t" << genotype << "\t" << genp << endl << endl;
#endif
		delete [] happrob;
		delete [] genprob;
	}
}

//double CHMM::haplotypeProbability(int refct, int altct, int errct, int g)
void CHMM::haplotypeProbability(vector<SNP*>::iterator snp_it, double probs[3])
{
	char ref = (*snp_it)->GetRef();
	char alt = (*snp_it)->GetAlt();
	//double qual = (*snp_it)->GetQualScore();
	//double qualscore = pow(10.0,-(qual/10.0)); // P(alt==right)
	int real_count = 0, ref_ct = 0, alt_ct = 0;
	double errate = 0.02;

#ifdef FULLDEBUG
cout << "SNP: " << (*snp_it)->GetPos() << endl;
#endif

	probs[0] = 1.0;
	probs[1] = 1.0;
	probs[2] = 1.0;
	// Here I assume that if an observed allele is not ref (or alt in case of BB),
	// then it belongs to the other haplotype.
	for(int count=0; count<(*snp_it)->GetReadCount(); count++) {
		READ *rd = (*snp_it)->GetRead(count);
		char all;
		double qual, qualscore;
		for(int s_pos=0; s_pos<rd->GetSnpCount(); s_pos++) {
			if(rd->GetSnp(s_pos)->GetPos() == (*snp_it)->GetPos()) {
				all = rd->GetAllele(s_pos);
				qual = rd->GetQualScore(s_pos)-33;
				qualscore = pow(10.0,-(qual/10.0));
				break;
			}
		}
		//char all = rd->GetAllele(count);
		double haprob = rd->GetHapProb();
		int hap = rd->GetHap();
		int state = (*snp_it)->GetEmissionState();
		if(state!=hap)
			haprob = 1.0 - haprob;

		if(rd->GetKnownCount()==0) {
		//if(rd->GetSnpCount()==0) {
#ifdef DEBUG
cout << "For snp " << (*snp_it)->GetPos() << " skipping read " << rd->GetPos() << " for happrob calculation" << endl;
#endif
			// samtools support
			 continue;
		}

#ifdef DEBUG
cout << "Estimating snp " << (*snp_it)->GetPos() << ", qualscore: " << qualscore << " with read " << rd->GetPos() << " bearing haplotype " << hap << " with prob " << haprob << " and state " << state << endl;
#endif

		if(all == ref) {
			probs[0] *= (haprob); //Do nothing
			probs[1] *= (0.5*haprob);
			probs[2] *= ((haprob) * (qualscore));
			ref_ct++;
		} else if(all == alt) {
			probs[0] *= ((haprob) * (qualscore));
			probs[1] *= (0.5*haprob);
			probs[2] *= (haprob); //Do nothing
			alt_ct++;
		} else {
///*
			probs[0] *= (errate);
			probs[1] *= (errate);
			probs[2] *= (errate);
//*/
/*
			probs[0] *= haprob*(1.0-errate);
			probs[1] *= haprob*(errate);
			probs[2] *= haprob*(errate);
*/
		}
		real_count++;
#ifdef DEBUG
cout << "For snp " << (*snp_it)->GetPos() << " bearing allele " << all << ref << alt << " on read " << rd->GetPos() << " with haprob " << haprob << " for happrob calculation " << probs[0] << ":" << probs[1] << ":" << probs[2] << endl;
#endif
	}
	if(real_count > 0) {
		probs[1] *= (Factorial(ref_ct+alt_ct)/(Factorial(ref_ct)*Factorial(alt_ct)*pow(2,ref_ct+alt_ct)));
/*
		probs[0] /= real_count;
		probs[1] /= real_count;
		probs[2] /= real_count;
		probs[0] = pow(2.7183, probs[0]);
		probs[1] = pow(2.7183, probs[1]);
		probs[2] = pow(2.7183, probs[2]);
*/
	} else {
		probs[0] = probs[1] = probs[2] = 1.0;
	}
}

void CHMM::genotypeProbability(vector<SNP*>::iterator snp_it, double probs[2])
{
	probs[0] = probs[1] = 1.0;
	char ref = (*snp_it)->GetRef();
	char alt = (*snp_it)->GetAlt();
	int read_count = (*snp_it)->GetReadCount();
//	double errate = 0.01;

	for(int count=0; count<read_count; count++) {
		READ *rd = (*snp_it)->GetRead(count);
		char all = rd->GetAllele(count);
		int hap = rd->GetHap();
		double happrob = rd->GetHapProb();

		if(all==ref) {
			probs[0] = probs[0];
			probs[1] *= (1-happrob);
		} else if(all==alt) {
			probs[0] *= (1-happrob);
			probs[1] = probs[1];
		} else {
			probs[0] *= errate;
			probs[1] *= errate;
		}
	}
	// Experimental: Logs for revisit above
	probs[0] = pow(2.7182, log(probs[0])/read_count);
	probs[1] = pow(2.7182, log(probs[1])/read_count);
	probs[2] = pow(2.7182, log(probs[2])/read_count);
}

void CHMM::GetCommonSnpList(CObs**obs, SNP**reads_snp_list, SNP** known_snp_list, int *common_snp_count, int *known_snp_count, int *common_index, int *known_index, int t)
{
	READ *prev_read = (((CFlexibleObs<READ*>*)(obs[t-1]))->Get(1));
	int prev_read_snp_count = prev_read->GetSnpCount();
	SNP **prev_snp_list = prev_read->GetSnpList();
	READ *curr_read = (((CFlexibleObs<READ*>*)(obs[t]))->Get(1));
	int curr_read_snp_count = curr_read->GetSnpCount();
	SNP **curr_snp_list = curr_read->GetSnpList();
	int it1 = 0, it2 = 0, it3 = 0, it4 = 0;
#ifdef FULLDEBUG
cout << "Start positions of two reads: " << prev_read->GetPos() << "\t" << curr_read->GetPos() << endl;
#endif
	while(it1<prev_read_snp_count && it2<curr_read_snp_count) {
		SNP* snp1 = prev_snp_list[it1];
		SNP* snp2 = curr_snp_list[it2];

		if(snp1->GetPos() == snp2->GetPos()) {
			if(snp1->GetKnown()>=1) {
				known_snp_list[*known_snp_count] = snp1;
				known_index[2*(*known_snp_count)] = it1;
				known_index[(2*(*known_snp_count))+1] = it2;
				known_snp_list[*known_snp_count]->IncrKnownOverlapCount();
				(*known_snp_count)++;
			}
			reads_snp_list[*common_snp_count] = snp1;
			common_index[2*(*common_snp_count)] = it1;
			common_index[(2*(*common_snp_count))+1] = it2;
			reads_snp_list[*common_snp_count]->IncrOverlapCount();
			(*common_snp_count)++;it1++;it2++;
		} else if(snp1->GetPos() < snp2->GetPos()) {
			it1++;
		} else if(snp1->GetPos() > snp2->GetPos()) {
			it2++;
		}
	}

}

//===============================================================================


double CHMM::compute_new_emission(SNP **reads_snp_list, int count, CObs **obs, int t, int *index, int hap, double obslik[3], double genlik[3])
{
	//double err_rate = errate;
	double err_rate = 0.02;
	double known_snp_rate = 0.1;
	double novel_snp_rate = 0.001;
	double gen_prior[3];
	double snp_rate = novel_snp_rate;
	double gen_lik[3] = {0.9, 0.01, 0.09};
	double obs_lik[3] = {0.0, 0.0, 1.0};
	double prob = 0.0;
	double lik = 0.0;

	char ref = reads_snp_list[count]->GetRef();
	char alt = reads_snp_list[count]->GetAlt();
	char all1 = ((CFlexibleObs<READ*>*)(obs[t-1]))->Get(1)->GetAllele(index[2*count]);
	char all2 = ((CFlexibleObs<READ*>*)(obs[t]))->Get(1)->GetAllele(index[(2*count)+1]);
	int obt = (ref==all1) ? ((ref==all2) ? 1 : 2) : ((ref==all2) ? 3 : 4);
	if(reads_snp_list[count]->GetKnown() == 1)
		snp_rate = known_snp_rate;
	gen_prior[0] = 1 - snp_rate - snp_rate*snp_rate;
	gen_prior[1] = snp_rate;
	gen_prior[2] = snp_rate*snp_rate;
	double good = (1-err_rate)*(1-err_rate);
	double bad = (1-err_rate)*err_rate;
	double ugly = err_rate*err_rate;
	for(int j=0; j<3; j++) {
		//lik += (reads_snp_list[count]->GetGenLik())[j] * gen_prior[j];
		lik += (reads_snp_list[count]->GetPosteriors())[j] * gen_prior[j];
	}

	for(int i=0; i<3; i++) {
		//gen_lik[i] = ((reads_snp_list[count]->GetGenLik())[i] * gen_prior[i])/lik;
		gen_lik[i] = ((reads_snp_list[count]->GetPosteriors())[i] * gen_prior[i])/lik;

		switch(i) {
		case 0: // genotype = A/A
			if(obt==1&&hap==1)
				obs_lik[i] = good;
			else if(obt==1&&hap==2)
				obs_lik[i] = good;
			else if(obt==2&&hap==1)
				obs_lik[i] = bad;
			else if(obt==2&&hap==2)
				obs_lik[i] = bad;
			else if(obt==3&&hap==1)
				obs_lik[i] = bad;
			else if(obt==3&&hap==2)
				obs_lik[i] = bad;
			else if(obt==4&&hap==1)
				obs_lik[i] = ugly;
			else if(obt==4&&hap==2)
				obs_lik[i] = ugly;
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
//obs_lik[i] = 0;
		break;
		case 1: // genotype = A/B
// #define DIVIDE
			if(obt==1&&hap==1)
				obs_lik[i] = good + ugly;
			else if(obt==1&&hap==2)
				obs_lik[i] = 2*bad;
			else if(obt==2&&hap==1)
				obs_lik[i] = 2*bad;
			else if(obt==2&&hap==2)
				obs_lik[i] = ugly + good;
			else if(obt==3&&hap==1)
				obs_lik[i] = 2*bad;
			else if(obt==3&&hap==2)
				obs_lik[i] = good + ugly;
			else if(obt==4&&hap==1)
				obs_lik[i] = ugly + good;
			else if(obt==4&&hap==2)
				obs_lik[i] = 2*bad;
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
#ifdef DIVIDE
			obs_lik[i] /= 2;
#endif
		break;
		case 2: // genotype = B/B
			if(obt==1&&hap==1)
				obs_lik[i] = ugly;
			else if(obt==1&&hap==2)
				obs_lik[i] = ugly;
			else if(obt==2&&hap==1)
				obs_lik[i] = bad;
			else if(obt==2&&hap==2)
				obs_lik[i] = bad;
			else if(obt==3&&hap==1)
				obs_lik[i] = bad;
			else if(obt==3&&hap==2)
				obs_lik[i] = bad;
			else if(obt==4&&hap==1)
				obs_lik[i] = good;
			else if(obt==4&&hap==2)
				obs_lik[i] = good;
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		}

#ifdef FULLDEBUG
cout << "\t" << gen_lik[i] << " " << obs_lik[i] << " " << gen_lik[i] * obs_lik[i];
#endif
		// REVISIT: Normalize here alone?
		prob += gen_lik[i] * obs_lik[i];
		obslik[i] = obs_lik[i];
		genlik[i] = gen_lik[i];
	}
	reads_snp_list[count]->addEmission(hap,prob);
#ifdef FULLDEBUG
cout << "\t" << prob << "\t" << log(prob) << endl;
//cout << "Posterior:";
	for(int p=0;p<3;p++) {
//		cout << "\t" << (gen_lik[p]*obs_lik[p])/prob;
	}
//cout << endl;
#endif
	//return gen_lik[1]*obs_lik[1];
	return prob;
}

double CHMM::BaumWelchCore(CObs **obs, long T, double *gamma, double **xi,
				 boolean doLog)
// Operations on a single observation sequence in the Baum-Welch loop are grouped here
{
  int	i, j, t;
  double logProb;
  double bjBeta;
  double *scale;
  double **alpha, **beta;

  alpha = SetMatrix(T, mN);// different size every time
  beta = SetMatrix(T, mN);

  scale = SetVector(T);
  
  logProb = Forward(alpha, scale, obs, T, doLog);
  Backward(beta, scale, obs, T);

  for (t = 1; t <= T - 1; t++){
    for (j = 1; j <= mN; j++) {
      gamma[j] = alpha[t][j] * beta[t][j];
      bjBeta =  beta[t+1][j] * mB->at(j, obs[t+1]);
      for (i = 1; i <= mN; i++){ 
	 xi[i][j] = alpha[t][i] *  mA->at(i, j) * bjBeta;
      }
    }//end for j
    Normalize(xi, mN, mN);
    Normalize(gamma, mN);
    
    mA->BWSum(xi);
    mB->BWSum(gamma, obs[t]);

    if(t==1){
      mPi->BWSum(gamma);
    }
  }// end for t

  // Step for t = T
  for (j = 1; j <= mN; j++) {
    gamma[j] = alpha[T][j] * beta[T][j];
  }
  Normalize(gamma, mN);
  mB->BWSum(gamma, obs[T]);

  delete [] alpha[1];
  delete [] alpha;
  delete [] beta[1];
  delete [] beta;
  delete [] scale;
  
  return logProb;
}

double CHMM::BaumWelchCoreEmission(CObs **obs, long T, ostream &distanceOutput, double *gamma, double **xi,
				 boolean doLog, long snp_start, long snp_end)
// Operations on a single observation sequence in the Baum-Welch loop are grouped here
{
  int	i, j, t, ind;
  double logProb;
  double backProb;
  double bjBeta;
  double *scale;
  double **alpha, **beta;

  alpha = SetMatrix(T, mN);// different size every time
  beta = SetMatrix(T, mN);

  scale = SetVector(T);
 
  //logProb = ForwardEmission(alpha, scale, obs, T, doLog);
  //BackwardEmission(beta, scale, obs, T);
  logProb = ForwardLogEmission(alpha, scale, obs, T, doLog);
  backProb = BackwardLogEmission(beta, scale, obs, T);

#ifdef DEBUG
cout << "logProb = " << logProb << " and backProb = " << backProb << endl << endl;
#endif

  //distanceOutput.seekp((streamoff)0, ios_base::beg);
  // Calculation of gammas (and thus respective posterior probabilities) is tricky here.
  // The treatment of the j index in the alpha matrix does not correspond to j->haplotype
  // Rather, it is j-> h(l)==h(l-1) ? 1 : 2;
  // This is so because in the compute_emission function, we need to specify j as hap1==hap2, not hap val itself
  // However, we cannot do the above math over there because the posteriors (and thus the h(l-1)) can only be
  // calculated after the Backward algorithm has been run as well. Hence we do that in the loop that follows.
  // The first read remains an exception though. We know for sure for the first read that j->haplotype.
  // This is used to start the ball rolling for later reads.
  // At the same time, there is no such requirement for beta probability calculations and j->haplotype.
  // The distinction can be noted in the loop below.

// RESUME debugging on gamma values here.
    for (j = 1; j <= mN; j++) {
    	gamma[j] = alpha[1][j]+beta[1][j];
#ifdef DEBUG
cout << "gamma[" << j << "] = " << gamma[j] << endl;
#endif
    }

    READ *nd = (((CFlexibleObs<READ*>*)(obs[1]))->Get(1));
    int int_hap = gamma[1] > gamma[2] ? 1: 2;
    //assert(int_hap==1);
    // samtools support here
    if(nd->GetKnownCount()>0) {
    //if(nd->GetSnpCount()>0) {
	double normel = pow(2.7183, gamma[1] - logProb) + pow(2.7183, gamma[2] - logProb);
#ifdef DEBUG
cout << "gamma[int_hap] = " << gamma[int_hap] << endl;
cout << "Assigned read " << nd->GetPos() << " haplotype " << int_hap << " and probablity of " << pow(2.7183, gamma[int_hap] - logProb)/normel << endl;
#endif
    	nd->assignHaplotype(int_hap, pow(2.7183, gamma[int_hap] - logProb)/normel);
    } else {
#ifdef FULLDEBUG
cout << "Skipped assigning haplotype to read " << nd->GetPos() << endl;
#endif
    }
    	LogNormalize(gamma, mN);
    	mB->BWSum(gamma, 1);

  for (t = 2; t <= T; t++){
    for (j = 1; j <= mN; j++) {
      // int_hap always corresponds to the higher probable state for the previous read.
      if(j==int_hap)
        gamma[j] = alpha[t][1] + beta[t][j];
      else
        gamma[j] = alpha[t][2] + beta[t][j];
    }//end for j

  // UpdateGenotype in posterior field
    // Not going to assign the haplotype here, only haplotype probability for haplotype with higher posterior.
    READ *nd = (((CFlexibleObs<READ*>*)(obs[t]))->Get(1));
    int_hap = gamma[1] > gamma[2] ? 1: 2;
    // samtools support here
    if(nd->GetKnownCount()>0) {
    //if(nd->GetSnpCount()>0) {
	double normel = pow(2.7183, gamma[1] - logProb) + pow(2.7183, gamma[2] - logProb);
#ifdef DEBUG
cout << "gamma[int_hap] = " << gamma[int_hap] << endl;
cout << "Assigned read " << nd->GetPos() << " haplotype " << int_hap << " and probablity of " << pow(2.7183, gamma[int_hap] - logProb)/normel << endl;
#endif
    	nd->assignHaplotype(int_hap, pow(2.7183, gamma[int_hap] - logProb)/normel);
    } else {
#ifdef FULLDEBUG
cout << "Skipped assigning haplotype to read " << nd->GetPos() << endl;
#endif
    }
    	LogNormalize(gamma, mN);
    	mB->BWSum(gamma, t);
  }// end for t

  UpdateBaumWelchGenotypePosteriors(snp_start, snp_end);

  // Step for t = T
  /*
  for (j = 1; j <= mN; j++) {
    gamma[j] = alpha[T][j] + beta[T][j];
  }
  Normalize(gamma, mN);
  mB->BWSum(gamma, T);
  */

  delete [] alpha[1];
  delete [] alpha;
  delete [] beta[1];
  delete [] beta;
  delete [] scale;
  
  return logProb;
}

//===============================================================================

double CHMM::IterBaumWelch(CObsSeq *obsSeq, double *gamma, double **xi)
{
	  double deltaAB, deltaA, deltaB, deltaPi, delta;
	  double logProb;
	  int i;
	  const boolean NOLOG = FALSE;

	  mA->StartIter();// Zero sums used to cumulate results from each sequence
	  mB->StartIter();
	  mPi->StartIter();
  
	  for(i=1;i<=obsSeq->mNbSequences;i++){ // Loop over observation files:

	    logProb = BaumWelchCore(obsSeq->mObs[i], obsSeq->mNbObs[i], gamma, xi, NOLOG);

	  }// end loop over observation files

	  deltaA = mA->EndIter();
	  deltaB = mB->EndIter();
	  deltaPi = mPi->EndIter();
	  
	  deltaAB = deltaA > deltaB ? deltaA : deltaB;
	  delta = deltaAB > deltaPi ? deltaAB : deltaPi;

	  return delta;
}

double CHMM::IterEmissionBaumWelch(CObsSeq *obsSeq, ostream &distanceOutput, double *gamma, double **xi, long snp_start, long snp_end)
{
	  double deltaAB, deltaA, deltaB, deltaPi, delta;
	  double logProb;
	  int i;
	  const boolean NOLOG = FALSE;

	  mB->StartIter();
  
	  for(i=1;i<=obsSeq->mNbSequences;i++){ // Loop over observation files:

	    logProb = BaumWelchCoreEmission(obsSeq->mObs[i], obsSeq->mNbObs[i], distanceOutput, gamma, xi, NOLOG, snp_start, snp_end);

	  }// end loop over observation files

	  //mB->UnLog();
	  deltaB = mB->EndLogIter();
	  
	  return deltaB;
}

//===============================================================================

// Supposed to run forward and backward algorithms
// Obtain P(x) and thus posteriors
// Perform haplotype calling
// Perform genotype calling
// Perform EM

void CHMM::FindFBDistance(CObsSeq *obsSeq, ostream &outFile, long snp_start, long snp_end)
{
	double *scale;
	double **alpha, **beta;
	double logProb;
	int i, j, em, t, k;
	double **posterior;
	double post;

	for(i=1;i<=obsSeq->mNbSequences;i++){ // Loop over observation files:
		long T = obsSeq->mNbObs[i];
		CObs **obs = obsSeq->mObs[i];
		alpha = SetMatrix(T, mN);// different size every time
		beta = SetMatrix(T, mN);
		scale = SetVector(T);
		int *q = SetIntVector(T);
		posterior = SetMatrix(T,mN);

		// logProb contains log P(x)
		// Haplotype Calling
		//logProb = ForwardAlgo(alpha, scale, obs, T, TRUE);
		//BackwardAlgo(beta, scale, obs, T);
		logProb = Forward(alpha, scale, obs, T, FALSE);
		Backward(beta, scale, obs, T);

		for(t=1;t<=T;t++) {
			READ *nd = (((CFlexibleObs<READ*>*)(obs[t]))->Get(1));
			for(k=1;k<=mN;k++) {
				// Experimental: P(x) cancels scaling for alpha and beta
				// so posterior is just their product
				double st = 1.0;
				for(j=T; j>t; j--) {
					st *= scale[j];
				}
				post = (alpha[t][k]*beta[t][k]*scale[t]);
#ifdef FULLDEBUG
cout << endl;
cout << "Scale = " << scale[t] << endl;
cout << "Scaled alpha = " << alpha[t][k] << endl;
cout << "beta = " << beta[t][k] << endl;
cout << "Scaled beta = " << beta[t][k]*scale[t] << endl;
#endif
				posterior[t][k] = post;
			}
			// REVSIT: hardcoding #states.
			q[t] = posterior[t][1] > posterior[t][2] ? 1 : 2;
#ifdef FULLDEBUG
cout << "q[t] = " << q[t] << endl;
cout << "posterior = " << posterior[t][q[t]] << endl;
#endif
			nd->assignHaplotype(q[t],posterior[t][q[t]]);
		}

		//UpdateGenotypes(snp_start, snp_end);

		for(j=2;j<=T;j++) {
			outFile << q[j] << "\t" << ((reads_list)[reads_list.size()-T+j-1])->GetPos() << endl;
		}
		delete [] q;
	}
}

void CHMM::LearnBaumWelch(CObsSeq *obsSeq)
// Follows (6.110) p. 369 of Rabiner-Huang
// The Baum-Welch loop is done over all the sequences
{
	int  iCount = 0;
	double delta;

	double *gamma = SetVector(mN);
	double **xi = SetMatrix(mN, mN);
	
	mA->Start();// Allocate sums used to cumulate results from each sequence
	mB->Start();
	mPi->Start();

// Baum-Welch loop starts here
	do  {	
	  delta = IterBaumWelch(obsSeq, gamma, xi);

	  iCount++;
	  cout<< endl << "BW iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;

	}
	while(delta > DELTA);

	mPi->End();
	mB->End();
	mA->End();

	cout << endl << "num iterations " << iCount << endl;
//	cout << "logTotalProb: " << sumProbf << endl;
	
	delete [] xi[1];
	delete [] xi;
	delete [] gamma;
}

void CHMM::LearnEmissionBaumWelch(CObsSeq *obsSeq, ostream &distanceOutput, long snp_start, long snp_end)
// Follows (6.110) p. 369 of Rabiner-Huang
// The Baum-Welch loop is done over all the sequences
{
	int  iCount = 0;
	double delta;
int iter = 0;
	double *gamma = SetVector(mN);
	double **xi = SetMatrix(mN, mN);
	
	mB->Start();

// Baum-Welch loop starts here
	do  {
	  delta = IterEmissionBaumWelch(obsSeq, distanceOutput, gamma, xi, snp_start, snp_end);

	  iCount++;
	  cout<< endl << "BW iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;

	}
	//while(delta > DELTA);
	while(++iter<EM);

	mB->End();

	//UpdateGenotypes(snp_start, snp_end);
	// REVISIT: Save haplotype values in file

	cout << endl << "num iterations " << iCount << endl;
//	cout << "logTotalProb: " << sumProbf << endl;
	
	delete [] xi[1];
	delete [] xi;
	delete [] gamma;
}

//===============================================================================

double CHMM::SegmentalKMeansCore(CObs **obs, long T)
// Operations on a single observation sequence in the segmental K-means loop
{
  int	t;
  int thisQ, nextQ;
  double logProb;
  int *q;

  q = SetIntVector(T);// best state sequence
  
  logProb = ViterbiLog(obs, T, q);
  
  mPi->SKMSum(q[1]);

  nextQ = q[1];
  for (t = 1; t <= T - 1; t++){
    thisQ = nextQ;
    nextQ = q[t+1];
    mA->SKMSum(thisQ, nextQ);
    mB->SKMSum(thisQ, obs[t]);
  }// end for t

  // Step for t = T
  mB->SKMSum(nextQ, obs[T]);

  delete [] q;  
  return logProb;
}

//===============================================================================

double CHMM::IterSegmentalKMeans(CObsSeq *obsSeq)
{
	  double deltaAB, deltaA, deltaB, deltaPi, delta;
	  double logProb;
	  int i;

	  mA->StartIter();// Zero sums used to cumulate results from each sequence
	  mB->StartIter();
	  mPi->StartIter();
  
	  for(i=1;i<=obsSeq->mNbSequences;i++){ // Loop over observation files:

	    logProb = SegmentalKMeansCore(obsSeq->mObs[i], obsSeq->mNbObs[i]);

	  }// end loop over observation files

	  deltaA = mA->EndIter();
	  deltaB = mB->EndIter();
	  deltaPi = mPi->EndIter();
	  
	  deltaAB = deltaA > deltaB ? deltaA : deltaB;
	  delta = deltaAB > deltaPi ? deltaAB : deltaPi;

	  return delta;
}

//===============================================================================

void CHMM::LearnSegmentalKMeans(CObsSeq *obsSeq)
// Follows (6.15.2) p. 383 of Rabiner-Huang
// The Segmental K-Means loop is done over all the sequences
{
	int  iCount = 0;
	double delta;

	mA->Start();// Allocate sums used to cumulate results from each sequence
	mB->Start();
	mPi->Start();

// Baum-Welch loop starts here
	do  {	
	  delta = IterSegmentalKMeans(obsSeq);

	  iCount++;

	  cout<< endl << "SKM iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;

	}
	while(delta > DELTA);

	mPi->End();
	mB->End();
	mA->End();

	cout << endl << "num iterations " << iCount << endl;
//	cout << "logTotalProb: " << sumProbf << endl;
}

//===============================================================================

void CHMM::LearnHybridSKM_BW(CObsSeq *obsSeq)
// Follows (6.15.2) p. 383 of Rabiner-Huang
// Combine Segmental K-Means and Baum-Welch
{
	int  iCount = 0;
	double delta;
	double *gamma = SetVector(mN);
	double **xi = SetMatrix(mN, mN);

	mA->Start();// Allocate sums used to cumulate results from each sequence
	mB->Start();
	mPi->Start();

	for(int i=0;i<10;i++){
	  delta = IterBaumWelch(obsSeq, gamma, xi);
	  iCount++;
	  cout<< endl << "BW iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;
	}
// Baum-Welch loop starts here
	do  {	
	  delta = IterSegmentalKMeans(obsSeq);

	  iCount++;

	  cout<< endl << "SKM iteration no. "<< iCount << endl;
	  cout << "delta = " << delta <<endl<<endl;

	}
	while(delta > DELTA);

	mPi->End();
	mB->End();
	mA->End();

	cout << endl << "num iterations " << iCount << endl;
//	cout << "logTotalProb: " << sumProbf << endl;
}

//===============================================================================


CObsSeq* CHMM::GenerateSequences(long nbSequences, long nbObs, int seed)
{
  int i;
  int anyState = 1;
  CObs* obsType = mB->PickObservation(anyState);// to pass observation type
  CObsSeq* obsSeq = new CObsSeq(obsType, nbSequences, nbObs);

  MyInitRand(seed);

  for(i=1;i<=nbSequences;i++){
	   obsSeq->mObs[i] = GenerateObservations(nbObs);
  }
  return obsSeq;
}

//===============================================================================


CObs** CHMM::GenerateObservations(long T)
// Generate an observation sequence of length T using A, B and Pi
{
   int t = 1;
	int currentState;
	CObs** obs;

	obs = new CObs*[T+1];
	assert(obs != NULL);

        currentState = mPi->PickInitialState();
        mA->InitDuration(currentState);
        obs[1] = mB->PickObservation(currentState);
 
        for (t = 2; t <= T; t++) {
                currentState =  mA->PickNextState(currentState);
                obs[t] =  mB->PickObservation(currentState);
        }
	return obs;
}

//===============================================================================

void CHMM::PrintStatesAndExpectedObs(CObsSeq *obsSeq,
				     ostream& stateFile, ostream& bestObsFile)
// Print sequences of most probable states 
// and sequences of expected observations 
// corresponding to those states
{
	long i, j, t, T, nbSequences;
	double logProb;
 	int *q;// most probable state sequence
	CObs **stateToObsMap;// list of expected observations for each state
	CObs *expectedObs;

	obsSeq->PrintHeader(stateFile);
	obsSeq->PrintHeader(bestObsFile);
	
	stateToObsMap = mB->MapStateToObs();// expected observation for each state

	for(i=1;i<=obsSeq->GetNbSequences();i++){ // Loop over observation files:
		T = obsSeq->mNbObs[i];
		stateFile <<"T= "<< T << endl;
		bestObsFile <<"T= "<< T << endl;
		
  		q = SetIntVector(T);// best state sequence
		logProb = ViterbiLog(obsSeq->mObs[i], T, q);

		for (t=1; t <= T; t++){
			stateFile << q[t] << " ";
			
			expectedObs = stateToObsMap[q[t]];
			expectedObs->Print(bestObsFile);
		}
		stateFile << endl;
		bestObsFile << endl;
		delete [] q;
	}
#if 1
	for (j=1; j<= mN; j++){
		delete stateToObsMap[j];
	}
#endif
	delete [] stateToObsMap;
}

//===============================================================================

#if 0

CObsSeq* CHMM::ReadSequences(ifstream &inputSeqFile)
// Read observation file into a data structure to make training faster
{
  long i, j, T, nbSequences;
  char magicID[32];
  CObs* obs;
  CObsSeq* obsSeq = new CObsSeq;

  mB->ReadFileHeader(inputSeqFile);// P5 or P6

  inputSeqFile >> magicID;
  assert(strcmp(magicID, "nbSequences=")==0);
  inputSeqFile >> nbSequences;
  obsSeq->mNbSequences = nbSequences;
  obsSeq->mObsCount = 0;
  obsSeq->mObs = new CObs**[nbSequences+1];
  obsSeq->mNbObs = new long[nbSequences+1];

  for(i=1;i<=nbSequences;i++){
    cout <<"Sequence "<< i <<endl;
    inputSeqFile >> magicID;
    assert(strcmp(magicID, "T=")==0);

    inputSeqFile >> T;// nb of observations for each sequence
    obsSeq->mNbObs[i] =  T;
  	 obsSeq->mObsCount += T;
    obsSeq->mObs[i] = new CObs*[T+1];// This array is from 1 to T
    assert(obsSeq->mObs[i] != NULL);

    for(j=1; j <=T; j++){
      obs = mB->ReadObsFrom(inputSeqFile);// obs type depends on type of mB
      obsSeq->mObs[i][j] = obs;
    }
  }

  return obsSeq;
}

#endif

//===============================================================================

double CHMM::FindDistance(CObsSeq *obsSeq, ostream &outFile)
// Find product of probabilities for each sequence
{
  long i, T, nbSequences;

  double logProb;
  double sumProb = 0.0;
  double *scale;
  double **alpha, **beta;
  const boolean YESLOG = TRUE;
  long stepCount = 0;
  double distance;


  nbSequences = obsSeq->mNbSequences;
  //outFile << "nbSequences= " << nbSequences << endl;

  for(i=1;i<=nbSequences;i++){ // Loop over observation files:
    T = obsSeq->mNbObs[i];
    stepCount += T;

    alpha = SetMatrix(T, mN);// different size every time
    beta = SetMatrix(T, mN);

    scale = SetVector(T);
  
    logProb = Forward(alpha, scale, obsSeq->mObs[i], T, YESLOG);

    sumProb += logProb;

//    cout<<i<<". T = "<<T<< ", logProb = " <<logProb<<", -prob/step = "<<-logProb/T << endl;
   // outFile << -logProb/T << endl;
    cout << -logProb/T << endl;

    delete [] alpha[1];
    delete [] alpha;
    delete [] beta[1];
    delete [] beta;
    delete [] scale;
  }

  cout << endl << "-Log Prob of all " << nbSequences << " sequences: "<< -sumProb << endl;
  cout << "Nb of steps: " << stepCount << endl;

  distance = - sumProb/ stepCount; // average distance per step
  cout << "Average Forward distance: " << distance << endl;
  
#if 0
  double averageProb = exp(-distance);
  cout << endl << "Average Forward prob = "<< averageProb << endl;
#endif
 return distance;

}

//===============================================================================
//GLOBAL
//double CHMM::FindViterbiDistance(CObsSeq *obsSeq, ostream &outFile, vector<READ*> *reads_list, vector<SNP*> *snp_list)
double CHMM::FindViterbiDistance(CObsSeq *obsSeq, ostream &outFile, ostream &gtFile)
// Find for each sequence log prob corresponding to best state segmentation
{
  long i, T, nbSequences;

  double logProb;
  double sumProb = 0.0;
  int *q;// most probable state sequence
  double *prob;
  long stepCount = 0;
  double distance;


  nbSequences = obsSeq->mNbSequences;

  for(i=1;i<=nbSequences;i++){ // Loop over observation files:
    // T = #symbols
    T = obsSeq->mNbObs[i];
    stepCount += T;
  
    q = SetIntVector(T);// best state sequence
    prob = SetVector(T);

    for(int em=0; em<EM; em++) {
    	// Haplotype calling
    	logProb = ViterbiLog(obsSeq->mObs[i], T, q, prob);
    	sumProb += logProb;

    	// Genotype calling
    	UpdateGenotypes(0,0);
    }

    cout << -logProb/T << endl;
    int j;
    for(j=1;j<=T;j++) {
	//outFile << q[j] << ":" << prob[j] << ": ";
	outFile << q[j] << "\t" << ((reads_list)[j-1])->GetPos() << endl;
	int k_max = ((reads_list)[j-1])->GetSnpCount();
    }
    delete [] q;
  }

  cout << "-Log Prob of all " << nbSequences << " sequences: "<< -sumProb << endl;
  cout << "Nb of steps: " << stepCount << endl;

  distance = - sumProb/ stepCount; // average distance per step
  
#if 0
  double averageProb = exp(-distance);
  cout << endl << "Average Viterbi prob = "<< averageProb << endl;
#endif
 return distance;

}

//===============================================================================

double CHMM::FindCrossEntropyDistance(CObsSeq *obsSeq, ostream &outFile)
// Find for each sequence log prob corresponding to best state segmentation
{
  long i, T, nbSequences;

  const double kWeight = 1.0;
  double logProb, weightedLogProb;
  double sumProb = 0.0;
  int *q;// most probable state sequence
  long stepCount = 0;
  double distance;
  double logQToPiProb;


  nbSequences = obsSeq->mNbSequences;
  outFile << "nbSequences= " << nbSequences << endl;

  for(i=1;i<=nbSequences;i++){ // Loop over observation files:
    T = obsSeq->mNbObs[i];
    stepCount += T;
  
    q = SetIntVector(T);// best state sequence
    logProb = ViterbiLog(obsSeq->mObs[i], T, q);
    logQToPiProb = FindQToPiProb(T, q);
    weightedLogProb = (kWeight*logProb + logQToPiProb)/(1.0+kWeight);
    sumProb += weightedLogProb;

    outFile <<-weightedLogProb/T << endl;
    cout <<-logProb/T << ", "<< -logQToPiProb/T << endl;
    delete [] q;
  }

  cout << endl << "-Log Prob of all " << nbSequences << " sequences: "<< -sumProb << endl;
  cout << "Nb of steps: " << stepCount << endl;

  distance = - sumProb/ stepCount; // average distance per step
  cout << "Average cross-entropy distance: " << distance << endl;
  
#if 0
  double averageProb = exp(-distance);
  cout << endl << "Average cross-entropy prob = "<< averageProb << endl;
#endif
 return distance;

}

//===============================================================================

double CHMM::FindQToPiProb(long T, int *q)
// Find probability of seeing a Pi distribution given that we have a Qi distribution
{
  int i, t;
  double *probQ;
  double piProb, qProb;
  double qToPiProb = 0.0;

  probQ = SetVector(mN);// prob distribution of observed Qs
  SetToZero(probQ, mN);

  for (t = 1; t <= T; t++){
    probQ[q[t]]++;
  }
  NonZeroNormalizeRow(probQ, mN);

  for(i=1;i<=mN;i++){
    piProb = mPi->at(i);
    qProb = probQ[i];
    qToPiProb += piProb * log(qProb/piProb);
  }
  qToPiProb *= T;// each number of observations is T * mPi->at(i)

  delete [] probQ;

  return qToPiProb;
}

//===============================================================================

void CHMM::Print(ostream &outFile)
{
	mA->Print(outFile);
	mB->Print(outFile);
	mPi->Print(outFile);
}

//===============================================================================

//===============================================================================
//===============================================================================

