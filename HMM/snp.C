#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>

//#define DEBUG

using namespace std;

#include "snp.h"

//SNP::SNP(long snp_pos, char snp_ref, char snp_alt, int type, READ *read, vector<string> gl3, bool known_par)
SNP::SNP(string chr_s, long snp_pos, char snp_ref, char snp_alt, vector<string> gl3, int known_par, double qualscore)
{
	strcpy(chr,chr_s.c_str());
	known = known_par;
	ref = snp_ref;
	alt = snp_alt;
	refcount = altcount = errcount = 0;
	count = -1;
	posterior_count = 0;
	overlap_count = 1;
	known_overlap_count = 1;
	position = snp_pos;
	qual = qualscore;
	likelihood_ratio = -1;
	reads = new READ*[200];
	gl = new double[3];
	//posteriors = new double*[70];
	posterior = new double[3];

//	add_read(type, read);

	gl[0] = pow(10.0, -(atof(gl3[0].c_str()))/10.0);
	gl[1] = pow(10.0, -(atof(gl3[1].c_str()))/10.0);
	gl[2] = pow(10.0, -(atof(gl3[2].c_str()))/10.0);

	posterior[0] = gl[0];
	posterior[1] = gl[1];
	posterior[2] = gl[2];

	emission[1] = 0.0;
	emission[2] = 0.0;
}

SNP::~SNP()
{
//cout << "SNP Destructor called: " << position << endl;
	delete [] reads;
	delete [] gl;
//	delete[] posteriors;
	delete [] posterior;
}

char* SNP::GetChr()
{
	return chr;
}

long SNP::GetPos()
{
	return position;
}

void SNP::append(int type, READ *read)
{
	add_read(type, read);
}

void SNP::add_read(int type, READ *read)
{
	if(type==0)
		refcount++;
	else if(type==1)
		altcount++;
	else
		errcount++;
	reads[++count] = read;
}

void SNP::add_posteriors(double post[3])
{
	for(int i=0; i<3; i++) {
		posterior[i] = post[i];
	}
	posterior_count++;
}

void SNP::addEmission(int haplotype, double probability)
{
	emission[haplotype] = probability;
}

void SNP::assign_genotype(int gt, double genp)
{
	genotype = gt;
	genprob = genp;
}

void SNP::IncrOverlapCount()
{
	overlap_count++; 
}

void SNP::IncrKnownOverlapCount()
{
	known_overlap_count++; 
}

READ* SNP::GetRead(int pos)
{
	return reads[pos];
}

char SNP::GetRef()
{
	return ref;
}

char SNP::GetAlt()
{
	return alt;
}

int SNP::GetReadCount()
{
	return count + 1;
}

int SNP::GetRefCount()
{
	return refcount;
}

int SNP::GetAltCount()
{
	return altcount;
}

int SNP::GetErrCount()
{
	return errcount;
}

double* SNP::GetGenLik()
{
	return gl;
}

int SNP::GetKnown()

{
	return known;
}

double SNP::GetQualScore()
{
	return qual;
}

double* SNP::GetPosteriors()
{
	return posterior;
}

int SNP::GetEmissionState()
{
	if(emission[1]>emission[2])
		return 1;
	else if(emission[2]>emission[1])
		return 2;
	else {
		cout << "Error for snp " << position << " with equal emissions requested bearing value " << emission[1] << " " << emission[2] << endl;
	}
}

int SNP::GetOverlapCount()
{
	return overlap_count;
}

int SNP::GetKnownOverlapCount()
{
	return known_overlap_count;
}

/*
void SNP::CalculateLikelihoodRatio()
{
	double first, second;

	if(posterior_count>0) {
		for(int ct=0; ct<posterior_count; ct++) {
			posterior[0] += posteriors[ct][0];
			posterior[1] += posteriors[ct][1];
			posterior[2] += posteriors[ct][2];
		}
		posterior[0] /= posterior_count;
		posterior[1] /= posterior_count;
		posterior[2] /= posterior_count;
		first =  posterior[0]>posterior[1] ? (posterior[0]>posterior[2] ? posterior[0] : posterior[2]) : (posterior[1]>posterior[2] ? posterior[1] : posterior[2]);
//		second = posterior[0]>posterior[1] ? (posterior[1]>posterior[2] ? posterior[1] : posterior[2]) : (posterior[0]>posterior[2] ? posterior[0] : posterior[1]);
		if(posterior[0]==first)
			second = posterior[1] > posterior[2] ? posterior[1] : posterior[2];
		else if(posterior[1]==first)
			second = posterior[0] > posterior[2] ? posterior[0] : posterior[2];
		else
			second = posterior[0] > posterior[1] ? posterior[0] : posterior[1];
		likelihood_ratio = first/second;
#ifdef DEBUG
		cout << posterior[0]/posterior_count << "\t" << posterior[1]/posterior_count << "\t" << posterior[2]/posterior_count << "\t" << first/posterior_count << "\t" << second/posterior_count << "\t" << likelihood_ratio << endl;
#endif
	} else {
//		cout << "Snp " << position << " does not overlap among any reads" << endl;
	}
}

#ifdef DEBUG
void SNP::PrintPosterior()
{
	double posterior[3] = {0,0,0};

	if(posterior_count>0) {
		for(int ct=0; ct<posterior_count; ct++) {
			posterior[0] += posteriors[ct][0];
			posterior[1] += posteriors[ct][1];
			posterior[2] += posteriors[ct][2];
		}

		cout << known << "\t" << position << "\t" << count+1;
		cout << "\t" << refcount << "\t" << altcount << "\t" << errcount;
		cout << "\t" << posterior[0]/posterior_count << "\t" << posterior[1]/posterior_count << "\t" << posterior[2]/posterior_count << endl;
	} else {
		cout << "Snp " << position << " does not overlap among any reads" << endl;
	}
}
#endif

void SNP::PrintLR()
{
	CalculateLikelihoodRatio();
	cout << position << "\t" << posterior[0] << "\t" << posterior[1] << "\t" << posterior[2] << "\t" << likelihood_ratio << endl;
	//cout << known << "\t" << position << "\t" << likelihood_ratio << endl;
}
*/
