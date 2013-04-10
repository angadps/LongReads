
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/beta.hpp>

using namespace std;

#include "read.h"
#include "snp.h"
#include "hmm.h"

extern vector<SNP*> snp_list;
extern vector<READ*> reads_list;
extern int errate;
extern int coverage;

//===============================================================================

//#define DEBUG
//#define FULLDEBUG

//===============================================================================

//===============================================================================

inline long long int Factorial(int x) {
  return (x == 1 || x == 0 ? 1 : x * Factorial(x - 1));
}

long int nCr (int n, int r)
{
	int d = n-r;
	long int fact = 1;
	if(r==0)
		return fact;

	if(r>n-r) {
		r = n-r;
		d = n-r;
	}
	for(int i=n; i>d; i--)
		fact *= (long int)i;
	fact /= Factorial(r);
	return fact;
}

char GetPosAllele(int t, long pos)
{
	for(int i=0; i<reads_list[t]->GetSnpCount(); i++) {
		if(reads_list[t]->GetSnpList()[i]->GetPos()==pos)
			return reads_list[t]->GetAllele(i);
	}
	return 'N';
}

void NaiveBayes(int snp_start, int snp_end, int T)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double logProb;
	double bi1, aij, bjt1;

        SNP **known_snp_list = new SNP*[100];
        int *known_index = new int[200];
	int reads_list_size = reads_list.size();

	for (t = 1+reads_list_size-T; t <= reads_list_size-1; t++) {
		int known_snp_count = 0;
		bool overlap;
                double emission_list[3], haprob[3];
		double norm = 0.0;
		READ *pd = reads_list[t-1];
		READ *nd = reads_list[t];
		int hap = 0;
                int rd_start = (*pd).GetPos() < (*nd).GetPos() ? (*nd).GetPos() : (*pd).GetPos();
                int rd_end = (*pd).GetPos()+(*pd).GetLen() > (*nd).GetPos()+(*nd).GetLen() ? (*nd).GetPos()+(*nd).GetLen()-1 : (*pd).GetPos()+(*pd).GetLen()-1;

                GetSnpList(known_snp_list, &known_snp_count, known_index, t);
		(*nd).AddKnownCount(known_snp_count);
#ifdef DEBUG
if(known_snp_count>0)
cout << "For read " << nd->GetPos() << " added known overlapping count of " << known_snp_count << endl;
cout << "Read " << t << " and " << t+1 << ", " << rd_start << ", " << rd_end << " with " << known_snp_count << " known snps." << endl << endl;
cout << "SnpPos  R A LL..\tAA_gen AA_obs AA_genob\tAB_gen AB_obs AB_genob\tBB_gen BB_obs BB_genob\tprob\n";
#endif

		cout << "Throwing emissions from read: " << nd->GetPos() << endl;
		for(j=1;j<=2;j++) {
			cout << "Haplotype " << j << endl;
			emission_list[j] = 0;
			for(int count=0; count<known_snp_count; count++) {
				double emission;
                       	        SNP *sp = known_snp_list[count];
#ifdef DEBUG
cout << sp->GetPos() << " " << sp->GetRef() << " " << sp->GetAlt() << " " << (*nd).GetAllele(known_index[count]);
#endif
				emission = compute_new_emission(known_snp_list, count, t, known_index, j);
					emission_list[j] += emission;
			}
#ifdef DEBUG
cout << "Total Emission = " << emission_list[j] << endl;
#endif
	  	}
		norm = emission_list[1] + emission_list[2];
		if(norm>0.0) {
			haprob[1] = emission_list[1] / norm;
			haprob[2] = emission_list[2] / norm;
		} else {
			haprob[1] = 0.5;
			haprob[2] = 0.5;
		}
#ifdef DEBUG
cout << "Happrob[1] = " << haprob[1] << endl;
cout << "Happrob[2] = " << haprob[2] << endl;
#endif
		hap = haprob[1] > haprob[2] ? 1 : haprob == haprob ? t%2+1 : 2;
		double happrob = haprob[hap];
		hap = (*pd).GetHap() == hap ? 1 : 2;
		nd->assignHaplotype(hap, happrob);
	}

	UpdateBaumWelchGenotypePosteriors(snp_start, snp_end);
	FindSomaticMutations(snp_start, snp_end);

	delete [] known_index;
	delete [] known_snp_list;
}

void UpdateBaumWelchGenotypePosteriors(long snp_start, long snp_end)
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
		double prior[3], sam_prior[3];
		double *happrob = new double[3];
		double genp = 1.0;
		double norm = 0.0;
		double snp_rate = 0.001;

		if((*snp_it)->GetKnown()==1)
			snp_rate = 0.1;

		sam_prior[0] = (*snp_it)->GetPosteriors()[0];
		sam_prior[1] = (*snp_it)->GetPosteriors()[1];
		sam_prior[2] = (*snp_it)->GetPosteriors()[2];

		prior[0] = 1 - snp_rate - snp_rate/2.0;
		prior[1] = snp_rate;
		prior[2] = snp_rate/2.0;

#ifdef DEBUG
cout << "Prior het prob for snp " << (*snp_it)->GetPos() << " = " << prior[1] << endl;
#endif
		haplotypeProbability(snp_it, happrob);

		if(happrob[0]==happrob[1]==happrob[2]==1.0) {
			prior[0] = sam_prior[0];
			prior[1] = sam_prior[1];
			prior[2] = sam_prior[2];
		}

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
	}
}

void FindSomaticMutations(long snp_start, long snp_end)
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
		double happrob[3];
		double genp = 1.0;
		double norm = 0.0;
		double snp_rate = 0.001;

		if((*snp_it)->GetKnown()==1)
			snp_rate = 0.1;

		prior[0] = (*snp_it)->GetPosteriors()[0];
		prior[1] = (*snp_it)->GetPosteriors()[1];
		prior[2] = (*snp_it)->GetPosteriors()[2];

		prior[0] = 1 - snp_rate - snp_rate*snp_rate;
		prior[1] = snp_rate;
		prior[2] = snp_rate*snp_rate;

#ifdef DEBUG
cout << "Prior het prob for som " << (*snp_it)->GetPos() << " = " << prior[1] << endl;
#endif

		somaticHaplotypeProbability(snp_it, happrob);

		for(g=0; g<3; g++) {
			norm += prior[g]*happrob[g];
		}

		for(g=0; g<3; g++) {
			post[g] = (prior[g]*happrob[g])/norm;
		}

#ifdef DEBUG
		if(post[1] > 1.0)
			cout << "Posterior het for somatic " << (*snp_it)->GetPos() << " > 1.0" << endl;
#endif

#ifdef DEBUG
cout << "Posterior het prob for som " << (*snp_it)->GetPos() << " = " << post[1] << endl;
#endif
    		(*snp_it)->add_somatic_posteriors(post);
	}
}

void haplotypeProbability(vector<SNP*>::iterator snp_it, double probs[3])
{
	char ref = (*snp_it)->GetRef();
	char alt = (*snp_it)->GetAlt();
	//double qual = (*snp_it)->GetQualScore();
	//double qualscore = pow(10.0,-(qual/10.0)); // P(alt==right)
	int real_count = 0, ref_ct = 0, alt_ct = 0, known = 0;
	double qual, qualscore;
	double errate = 0.02;

	probs[0] = 1.0;
	probs[1] = 1.0;
	probs[2] = 1.0;
	// Here I assume that if an observed allele is not ref (or alt in case of BB),
	// then it belongs to the other haplotype.
	for(int count=0; count<(*snp_it)->GetReadCount(); count++) {
		READ *rd = (*snp_it)->GetRead(count);
		READ *pd = (*snp_it)->GetRead(count-1);
		char all, pall = 'N';
		int state = 0;
		for(int s_pos=0; s_pos<rd->GetSnpCount(); s_pos++) {
			if(rd->GetSnp(s_pos)->GetPos() == (*snp_it)->GetPos()) {
				all = rd->GetAllele(s_pos);
				qual = rd->GetQualScore(s_pos)-33;
				qualscore = pow(10.0,-(qual/10.0));
				break;
			}
		}
		if(pd!=NULL) 
		for(int s_pos=0; s_pos<pd->GetSnpCount(); s_pos++) {
			if(pd->GetSnp(s_pos)->GetPos() == (*snp_it)->GetPos()) {
				pall = pd->GetAllele(s_pos);
				break;
			}
		}
		double haprob = rd->GetHapProb();
		int hap = rd->GetHap();
		if(pall!='N')
			if((hap==pd->GetHap()&&all!=pall) || (hap!=pd->GetHap()&&all==pall))
				if((*snp_it)->GetKnown()==1)
					haprob = haprob;

		if(rd->GetSnpCount()==0||qualscore>=0.01) {
#ifdef DEBUG
cout << "For snp " << (*snp_it)->GetPos() << " skipping read " << rd->GetPos() << " for happrob calculation" << endl;
#endif
			// samtools support
			 continue;
		}
		if(rd->GetKnownCount()>0)
			known++;

#ifdef DEBUG
//cout << "Estimating snp " << (*snp_it)->GetPos() << ", qualscore: " << qualscore << " with read " << rd->GetPos() << " bearing haplotype " << hap << " with prob " << haprob << " and state " << state << endl;
#endif
//		if((*snp_it)->GetKnown()==0&&alt_ct>1) { continue;}
		if(all == ref) {
			probs[0] *= ((1.0-qualscore)*haprob); //(haprob); //Do nothing
			probs[1] *= (0.5*haprob);
			probs[2] *= (qualscore * (1.0 - haprob));
			ref_ct++;
		} else if(all == alt) {
				probs[0] *= (qualscore * (1.0 - haprob)); // (1.0 - haprob)
				probs[1] *= (0.5*haprob);
				probs[2] *= ((1.0 - qualscore)*haprob); //(haprob); //Do nothing
			//} else {
			//	probs[0] *= (haprob); //Do nothing
			//	probs[1] *= (0.5);
			//	probs[2] *= ((1.0-haprob) * (qualscore));
			//}
			alt_ct++;
		} else {
			probs[0] *= (errate);
			probs[1] *= (errate);
			probs[2] *= (errate);
		}
		real_count++;
#ifdef DEBUG
cout << "For snp " << (*snp_it)->GetPos() << " qualscore " << qualscore << " bearing allele " << all << ref << alt << " on read " << rd->GetPos() << " with haprob " << haprob << " for happrob calculation " << probs[0] << ":" << probs[1] << ":" << probs[2] << endl;
#endif
	}
	boost::math::binomial_distribution<double> distr(ref_ct+alt_ct,0.5);
	double val = boost::math::pdf(distr,alt_ct);
	//probs[1] *= val;
	//if(real_count==0||known==0) { // || ((*snp_it)->GetKnown()!=1&&(alt_ct==1&&ref_ct<=4||alt_ct<=2&&ref_ct<=8))) {
	if(real_count==0) { // || ((*snp_it)->GetKnown()!=1&&(alt_ct==1&&ref_ct<=4||alt_ct<=2&&ref_ct<=8))) {
		probs[0] = probs[1] = probs[2] = 1.0;
	} 
#ifdef DEBUG
cout << "For snp " << (*snp_it)->GetPos() << " happrob calculation " << probs[0] << ":" << probs[1] << ":" << probs[2] << endl;
#endif
}

void somaticHaplotypeProbability(vector<SNP*>::iterator snp_it, double probs[3])
{
	char ref = (*snp_it)->GetRef();
	char alt = (*snp_it)->GetAlt();
	int real_count = 0, ref_ct = 0, alt_ct = 0, ref1_ct = 0, alt1_ct = 0, ref2_ct = 0, alt2_ct = 0;
	double errate = 0.02;
	double prevelance = 1.0 - (double)(coverage)*0.01;
	double min_prevelance = 1.0 - (double)(coverage)*0.02;
	double threshold = 0.1 * (double)(coverage);
	double probs1[3], probs2[3];
	probs[0] = 1.0; probs[1] = 1.0; probs[2] = 1.0;
	probs1[0] = 1.0; probs1[1] = 1.0; probs1[2] = 1.0;
	probs2[0] = 1.0; probs2[1] = 1.0; probs2[2] = 1.0;

	// Here I assume that if an observed allele is not ref (or alt in case of BB),
	// then it belongs to the other haplotype.
	for(int count=0; count<(*snp_it)->GetReadCount(); count++) {
		char all, pall = 'N';
		double qual, qualscore;
		READ *rd = (*snp_it)->GetRead(count);
		READ *pd = (*snp_it)->GetRead(count-1);
		for(int s_pos=0; s_pos<rd->GetSnpCount(); s_pos++) {
			if(rd->GetSnp(s_pos)->GetPos() == (*snp_it)->GetPos()) {
				all = rd->GetAllele(s_pos);
				qual = rd->GetQualScore(s_pos)-33;
				qualscore = pow(10.0,-(qual/10.0));
				break;
			}
		}
		if(pd!=NULL) 
		for(int s_pos=0; s_pos<pd->GetSnpCount(); s_pos++) {
			if(pd->GetSnp(s_pos)->GetPos() == (*snp_it)->GetPos()) {
				pall = pd->GetAllele(s_pos);
				break;
			}
		}
		int hap = rd->GetHap();
		double haprob = rd->GetHapProb();

		if(pall!='N')
			if((hap==pd->GetHap()&&all!=pall) || (hap!=pd->GetHap()&&all==pall))
				if((*snp_it)->GetKnown()==1)
					haprob = haprob;

		if(rd->GetSnpCount()==0) {
#ifdef FULLDEBUG
cout << "For som " << (*snp_it)->GetPos() << " skipping read " << rd->GetPos() << " for happrob calculation" << endl;
#endif
			// samtools support
			 continue;
		}

#ifdef FULLDEBUG
cout << "Estimating som " << (*snp_it)->GetPos() << ", qualscore: " << qualscore << " with read " << rd->GetPos() << " bearing haplotype " << hap << " with prob " << haprob << endl;
#endif

//prevelance = 0.7;
		if(hap==1) {
			if(all == ref) {
				probs1[0] *= ((haprob)*(1.0-qualscore));
				probs1[1] *= (1.0 - (haprob + qualscore - 2.0*haprob*qualscore));
				probs1[2] *= ((1.0-haprob) * (qualscore));
				probs2[0] *= ((haprob)*(1.0-qualscore));
				probs2[1] *= ((2.0*haprob-1.0)*(prevelance*qualscore-2.0*prevelance*qualscore)+(1.0-haprob));
				probs2[2] *= ((1.0-haprob) * (qualscore));
				ref1_ct++;
			} else if(all == alt) {
				probs1[0] *= ((1.0-haprob) * (qualscore));
				probs1[1] *= ((haprob + qualscore - 2.0*haprob*qualscore));
				probs1[2] *= ((haprob)*(1.0-qualscore));
				probs2[0] *= ((1.0-haprob) * (qualscore));
				probs2[1] *= ((1.0-2.0*haprob)*(prevelance*qualscore-2.0*prevelance*qualscore)+haprob);
				probs2[2] *= ((haprob)*(1.0-qualscore));
				alt1_ct++;
			}
		} else if(hap==2) {
			if(all == ref) {
				probs1[0] *= ((haprob)*(1.0-qualscore));
				probs1[1] *= ((2.0*haprob-1.0)*(prevelance*qualscore-2.0*prevelance*qualscore)+(1.0-haprob));
				probs1[2] *= ((1.0-haprob) * (qualscore));
				probs2[0] *= ((haprob)*(1.0-qualscore));
				probs2[1] *= (1.0 - (haprob + qualscore - 2.0*haprob*qualscore));
				probs2[2] *= ((1.0-haprob) * (qualscore));
				ref2_ct++;
			} else if(all == alt) {
				probs1[0] *= ((1.0-haprob) * (qualscore));
				probs1[1] *= ((1.0-2.0*haprob)*(prevelance*qualscore-2.0*prevelance*qualscore)+haprob);
				probs1[2] *= ((haprob)*(1.0-qualscore));
				probs2[0] *= ((1.0-haprob) * (qualscore));
				probs2[1] *= ((haprob + qualscore - 2.0*haprob*qualscore));
				probs2[2] *= ((haprob)*(1.0-qualscore));
				alt2_ct++;
			}
		}
		real_count++;
#ifdef FULLDEBUG
cout << "For som " << (*snp_it)->GetPos() << " bearing allele " << all << ref << alt << " on read " << rd->GetPos() << " with haprob " << haprob << " for happrob calculation " << probs1[0] << ":" << probs1[1] << ":" << probs1[2] << " and " << probs2[0] << ":" << probs2[1] << ":" << probs2[2] << endl;
#endif
	}

	if(real_count > 0) {
		double ratio = 0.0;
		double ratio1 = ((double)ref1_ct)/((double)(ref1_ct+alt1_ct));
		double ratio2 = ((double)ref2_ct)/((double)(ref2_ct+alt2_ct));
		int ct_check1 = ref1_ct<=(int)(4*threshold) && alt1_ct==0 && ref2_ct>2*ref1_ct && alt2_ct<ref2_ct ? 1 : 0;
		int ct_check2 = ref2_ct<=(int)(4*threshold) && alt2_ct==0 && ref1_ct>2*ref2_ct && alt1_ct<ref1_ct ? 1 : 0;

		if(ref1_ct+alt1_ct==0) {
			probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
			ref_ct = ref2_ct; alt_ct = alt2_ct;
		} else if(ref2_ct+alt2_ct==0) {
			probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
			ref_ct = ref1_ct; alt_ct = alt1_ct;
		} else if(ratio1 > ratio2) {
			if(ct_check1) {
				probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
				ref_ct = ref1_ct; alt_ct = alt1_ct;
			} else {
				probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
				ref_ct = ref2_ct; alt_ct = alt2_ct;
			}
		} else if(ratio2 > ratio1) {
			if(ct_check2) {
				probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
				ref_ct = ref2_ct; alt_ct = alt2_ct;
			} else {
				probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
				ref_ct = ref1_ct; alt_ct = alt1_ct;
			}
		} else if(ref1_ct > ref2_ct) {
			probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
			ref_ct = ref2_ct; alt_ct = alt2_ct;
		} else if(ref2_ct > ref1_ct) {
			probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
			ref_ct = ref1_ct; alt_ct = alt1_ct;
		} else {
			probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
			ref_ct = ref1_ct; alt_ct = alt1_ct;
		}

		ratio = ((double)alt_ct)/((double)(ref_ct+alt_ct));
		threshold = (double)(ref_ct+alt_ct)/5.0;
		//if(ratio<prevelance && ref_ct+alt_ct > 0) {
#ifdef DEBUG
cout << "refct = " << ref_ct << "; altct = " << alt_ct << "; threshold = " << threshold << "; ratio = " << ratio << "; min_prevelance/2.0 = " << min_prevelance/2.0 << endl;
#endif

		//if(ref_ct+alt_ct>0 && (alt_ct>=threshold || ratio>=min_prevelance/2.0)) {
		if(ref_ct+alt_ct>0) { //&& (alt_ct>=2)) || ratio>=min_prevelance/2.0)) {

#ifdef DEBUG
cout << "Ratio = " << ratio << " prevelance = " << prevelance << " coverage " << coverage << " threshold = " << threshold << " and het prob = " << probs[1] << endl;
#endif
			boost::math::beta_distribution<double> distr(2.0*(ref_ct+alt_ct),ceil(threshold));
			double val = boost::math::pdf(distr,ratio)/coverage;

			//if(ratio=1.0)
				//val = 1.0;
			//probs[1] *= val;
#ifdef DEBUG
cout << "Val = " << val << ": New het prob = " << probs[1] << endl;
#endif
		} else {
			probs[1] = 0.0;
		}
			
	} else {
		probs[0] = probs[1] = probs[2] = 1.0;
	}
	if(probs[0]<0.0||probs[1]<0.0||probs[2]<0.0)
		cout << "Negative: " << probs[0] << ":" << probs[1] << ":" << probs[2] << endl;

#ifdef DEBUG
cout << "For som " << (*snp_it)->GetPos() << " happrob calculation " << ref1_ct << ";" << alt1_ct << ";" << ref2_ct << ";" << alt2_ct << ";" << probs1[0] << ":" << probs1[1] << ":" << probs1[2] << " and " << probs2[0] << ":" << probs2[1] << ":" << probs2[2] << endl;
#endif
}

/*
// Somatic mutation likelihood calculation does not depend on phase totally. A reference allele could still be observed on the alternate haplotype given prevelance < 1.0
void somaticHaplotypeProbability(vector<SNP*>::iterator snp_it, double probs[3])
{
	char ref = (*snp_it)->GetRef();
	char alt = (*snp_it)->GetAlt();
	int real_count = 0, ref_ct = 0, alt_ct = 0, ref1_ct = 0, alt1_ct = 0, ref2_ct = 0, alt2_ct = 0;
	double errate = 0.02;
	double prevelance = 1.0 - (double)(coverage)*0.01;
	double min_prevelance = 1.0 - (double)(coverage)*0.02;
	double threshold = 0.1 * (double)(coverage);
	double probs1[3], probs2[3];
prevelance = 0.6;
	probs[0] = 1.0; probs[1] = 1.0; probs[2] = 1.0;
	probs1[0] = 1.0; probs1[1] = 1.0; probs1[2] = 1.0;
	probs2[0] = 1.0; probs2[1] = 1.0; probs2[2] = 1.0;

	// Here I assume that if an observed allele is not ref (or alt in case of BB),
	// then it belongs to the other haplotype.
	for(int count=0; count<(*snp_it)->GetReadCount(); count++) {
		char all, pall = 'N';
		double qual, qualscore;
		READ *rd = (*snp_it)->GetRead(count);
		READ *pd = (*snp_it)->GetRead(count-1);
		for(int s_pos=0; s_pos<rd->GetSnpCount(); s_pos++) {
			if(rd->GetSnp(s_pos)->GetPos() == (*snp_it)->GetPos()) {
				all = rd->GetAllele(s_pos);
				qual = rd->GetQualScore(s_pos)-33;
				qualscore = pow(10.0,-(qual/10.0));
				break;
			}
		}
		if(pd!=NULL) 
		for(int s_pos=0; s_pos<pd->GetSnpCount(); s_pos++) {
			if(pd->GetSnp(s_pos)->GetPos() == (*snp_it)->GetPos()) {
				pall = pd->GetAllele(s_pos);
				break;
			}
		}
		int hap = rd->GetHap();
		double haprob = rd->GetHapProb();

		if(pall!='N')
			if((hap==pd->GetHap()&&all!=pall) || (hap!=pd->GetHap()&&all==pall))
				if((*snp_it)->GetKnown()==1)
					haprob = haprob;

		//if(rd->GetKnownCount()==0) {
		if(rd->GetSnpCount()==0) {
#ifdef FULLDEBUG
cout << "For som " << (*snp_it)->GetPos() << " skipping read " << rd->GetPos() << " for happrob calculation" << endl;
#endif
			// samtools support
			 continue;
		}

#ifdef FULLDEBUG
cout << "Estimating som " << (*snp_it)->GetPos() << ", qualscore: " << qualscore << " with read " << rd->GetPos() << " bearing haplotype " << hap << " with prob " << haprob << endl;
#endif

		if(hap==1) {
			if(all == ref) {
				probs1[0] *= (haprob);
		//		probs1[1] *= (0.5);
				probs1[2] *= ((haprob) * (qualscore));
				probs2[0] *= (haprob);
		//		probs2[1] *= 1;
				probs2[2] *= ((haprob) * (qualscore));
				ref1_ct++;
			} else if(all == alt) {
				probs1[0] *= ((haprob) * (qualscore));
		//		probs1[1] *= 1;
				probs1[2] *= (haprob);
				probs2[0] *= ((haprob) * (qualscore));
		//		probs2[1] *= (0.5);
				probs2[2] *= (haprob);
				alt1_ct++;
			}
		} else if(hap==2) {
			if(all == ref) {
				probs1[0] *= (haprob);
		//		probs1[1] *= 1;
				probs1[2] *= ((haprob) * (qualscore));
				probs2[0] *= (haprob);
		//		probs2[1] *= (0.5);
				probs2[2] *= ((haprob) * (qualscore));
				ref2_ct++;
			} else if(all == alt) {
				probs1[0] *= ((haprob) * (qualscore));
		//		probs1[1] *= (0.5);
				probs1[2] *= (haprob);
				probs2[0] *= ((haprob) * (qualscore));
		//		probs2[1] *= 1;
				probs2[2] *= (haprob);
				alt2_ct++;
			}
		}
		real_count++;
#ifdef FULLDEBUG
cout << "For som " << (*snp_it)->GetPos() << " bearing allele " << all << ref << alt << " on read " << rd->GetPos() << " with haprob " << haprob << " for happrob calculation " << probs1[0] << ":" << probs1[1] << ":" << probs1[2] << " and " << probs2[0] << ":" << probs2[1] << ":" << probs2[2] << endl;
#endif
	}

	if(real_count > 0) {
		double ratio = 0.0;
		double ratio1 = ((double)ref1_ct)/((double)(ref1_ct+alt1_ct));
		double ratio2 = ((double)ref2_ct)/((double)(ref2_ct+alt2_ct));
		int ct_check1 = ref1_ct<=(int)(4*threshold) && alt1_ct==0 && ref2_ct>2*ref1_ct && alt2_ct<ref2_ct ? 1 : 0;
		int ct_check2 = ref2_ct<=(int)(4*threshold) && alt2_ct==0 && ref1_ct>2*ref2_ct && alt1_ct<ref1_ct ? 1 : 0;

		if(ref1_ct+alt1_ct==0) {
			probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
			ref_ct = ref2_ct; alt_ct = alt2_ct;
		} else if(ref2_ct+alt2_ct==0) {
			probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
			ref_ct = ref1_ct; alt_ct = alt1_ct;
		} else if(ratio1 > ratio2) {
			if(ct_check1) {
				probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
				ref_ct = ref1_ct; alt_ct = alt1_ct;
			} else {
				probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
				ref_ct = ref2_ct; alt_ct = alt2_ct;
			}
		} else if(ratio2 > ratio1) {
			if(ct_check2) {
				probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
				ref_ct = ref2_ct; alt_ct = alt2_ct;
			} else {
				probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
				ref_ct = ref1_ct; alt_ct = alt1_ct;
			}
		} else if(ref1_ct > ref2_ct) {
			probs[0] = probs1[0]; probs[1] = probs1[1]; probs[2] = probs1[2];
			ref_ct = ref2_ct; alt_ct = alt2_ct;
		} else if(ref2_ct > ref1_ct) {
			probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
			ref_ct = ref1_ct; alt_ct = alt1_ct;
		} else {
			probs[0] = probs2[0]; probs[1] = probs2[1]; probs[2] = probs2[2];
			ref_ct = ref1_ct; alt_ct = alt1_ct;
		}

		ratio = ((double)alt_ct)/((double)(ref_ct+alt_ct));
		threshold = (double)(ref_ct+alt_ct)/5.0;
		//if(ratio<prevelance && ref_ct+alt_ct > 0) {
#ifdef DEBUG
cout << "refct = " << ref_ct << "; altct = " << alt_ct << "; threshold = " << threshold << "; ratio = " << ratio << "; min_prevelance/2.0 = " << min_prevelance/2.0 << endl;
#endif

		//if(ref_ct+alt_ct>0 && (alt_ct>=threshold || ratio>=min_prevelance/2.0)) {
		if(ref_ct+alt_ct>0) { //&& (alt_ct>=2)) || ratio>=min_prevelance/2.0)) {

#ifdef DEBUG
cout << "Ratio = " << ratio << " prevelance = " << prevelance << " coverage " << coverage << " threshold = " << threshold << " and het prob = " << probs[1] << endl;
#endif
			boost::math::beta_distribution<double> distr(2.0*(ref_ct+alt_ct),ceil(threshold));
			double val = boost::math::pdf(distr,ratio)/coverage;

			//if(ratio=1.0)
				//val = 1.0;
			probs[1] *= val;
#ifdef DEBUG
cout << "Val = " << val << ": New het prob = " << probs[1] << endl;
#endif
		} else {
			probs[1] = 0.0;
		}
			
	} else {
		probs[0] = probs[1] = probs[2] = 1.0;
	}
	if(probs[0]<0.0||probs[1]<0.0||probs[2]<0.0)
		cout << "Negative: " << probs[0] << ":" << probs[1] << ":" << probs[2] << endl;

#ifdef DEBUG
cout << "For som " << (*snp_it)->GetPos() << " happrob calculation " << ref1_ct << ";" << alt1_ct << ";" << ref2_ct << ";" << alt2_ct << ";" << probs1[0] << ":" << probs1[1] << ":" << probs1[2] << " and " << probs2[0] << ":" << probs2[1] << ":" << probs2[2] << endl;
#endif
}
*/

void genotypeProbability(vector<SNP*>::iterator snp_it, double probs[2])
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

void GetSnpList(SNP** known_snp_list, int *known_snp_count, int *known_index, int t)
{
	READ *prev_read = reads_list[t-1];
	int prev_read_snp_count = prev_read->GetSnpCount();
	SNP **prev_snp_list = prev_read->GetSnpList();
	READ *curr_read = reads_list[t];
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
			if(snp1->GetKnown()==1) {
				known_snp_list[*known_snp_count] = snp1;
				known_index[2*(*known_snp_count)] = it1;
				known_index[(2*(*known_snp_count))+1] = it2;
				known_snp_list[*known_snp_count]->IncrKnownOverlapCount();
				(*known_snp_count)++;
			}
			it1++;it2++;
		} else if(snp1->GetPos() < snp2->GetPos()) {
			it1++;
		} else if(snp1->GetPos() > snp2->GetPos()) {
			it2++;
		}
	}

}

void GetCommonSnpList(SNP**reads_snp_list, SNP** known_snp_list, int *common_snp_count, int *known_snp_count, int *common_index, int *known_index, int t)
{
	READ *prev_read = reads_list[t-1];
	int prev_read_snp_count = prev_read->GetSnpCount();
	SNP **prev_snp_list = prev_read->GetSnpList();
	READ *curr_read = reads_list[t];
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
			if(snp1->GetKnown()==1) {
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

double compute_new_emission(SNP **reads_snp_list, int count, int t, int *index, int hap)
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
	int pos = reads_snp_list[count]->GetPos();
	//char all1 = reads_list[t-1]->GetAllele(index[2*count]);
	//char all2 = reads_list[t]->GetAllele(index[2*count]);
	char all1 = GetPosAllele(t-1,pos);
	char all2 = GetPosAllele(t,pos);
	int obt = (ref==all1) ? ((ref==all2) ? 1 : 2) : ((ref==all2) ? 3 : 4);
	if(reads_snp_list[count]->GetKnown() == 1)
		snp_rate = known_snp_rate;
	gen_prior[0] = 1 - snp_rate - snp_rate*snp_rate;
	gen_prior[1] = snp_rate;
	gen_prior[2] = snp_rate*snp_rate;
	double good = (1-err_rate)*(1-err_rate);
	double bad = (1-err_rate)*err_rate;
	double ugly = err_rate*err_rate;

#ifdef DEBUG
cout << "\t" << all1 << "\t" << all2;
#endif
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

#ifdef DEBUG
cout << "\t" << gen_lik[i] << " " << obs_lik[i] << " " << gen_lik[i] * obs_lik[i];
#endif
		// REVISIT: Normalize here alone?
		prob += gen_lik[i] * obs_lik[i];
	}
	reads_snp_list[count]->addEmission(hap,prob);
#ifdef DEBUG
cout << "\t" << prob << "\t" << log(prob) << endl;
#endif
	//return gen_lik[1]*obs_lik[1];
	return prob;
}

