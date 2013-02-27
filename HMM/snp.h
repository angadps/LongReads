
#include<vector>

class READ;

class SNP {

public:

  SNP(void);
  //SNP(long snp_pos, char snp_ref, char snp_alt, int type, READ *read, vector<string>, bool known);
  SNP(string chr, long snp_pos, char snp_ref, char snp_alt, vector<string>, int known, double qual);

  ~SNP(void);

  char *GetChr();
  long GetPos();
  void append(int type, READ *read);
  READ *GetRead(int pos);
  char GetRef();
  char GetAlt();
  int GetReadCount();
  int GetRefCount();
  int GetAltCount();
  int GetErrCount();
  double* GetGenLik();
  int GetKnown();
  double GetQualScore();
  int GetOverlapCount();
  int GetKnownOverlapCount();
  int GetEmissionState();
  void add_posteriors(double posterior[3]);
  void assign_genotype(int gt, double genp);
  void addEmission(int haplotype, double probability);
  void IncrOverlapCount();
  void IncrKnownOverlapCount();
//  void PrintPosterior();
//  void PrintLR();
  double* GetPosteriors();

private:
  char chr[10];
  int known;
  char ref;
  char alt;
  int refcount;
  int altcount;
  int errcount;
  int count;
  int posterior_count;
  int genotype;
  int overlap_count;
  int known_overlap_count;
  long position;
  double qual;
  double likelihood_ratio;
  double genprob;
  READ **reads;
  double *gl;
//  double **posteriors;
  double *posterior;
  double emission[3];

void add_read(int type, READ *read);
//void CalculateLikelihoodRatio();

};

//===============================================================================
//===============================================================================
