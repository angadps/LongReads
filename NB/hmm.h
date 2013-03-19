
void NaiveBayes(int snp_start, int snp_end, int nbSymbols);
void GetCommonSnpList(SNP** reads_snp_list, SNP** known_snp_list, int *common_snp_count, int *known_snp_count, int *common_index, int *known_index, int t);
void GetSnpList(SNP** known_snp_list, int *known_snp_count, int *known_index, int t);
double compute_new_emission(SNP**reads_snp_list, int count, int t, int *index, int hap);
double compute_full_emission(SNP**snp_list, int count, int t, int *index, int hap);
void haplotypeProbability(vector<SNP*>::iterator snp_it, double happ[3]);
void somaticHaplotypeProbability(vector<SNP*>::iterator snp_it, double happ[3]);
void genotypeProbability(vector<SNP*>::iterator snp_it, double genp[2]);
void UpdateGenotypes(long start, long end);
void UpdateBaumWelchGenotypePosteriors(long start, long end);
void FindSomaticMutations(long start, long end);
char GetPosAllele(int t, long pos);

//===============================================================================
//===============================================================================
