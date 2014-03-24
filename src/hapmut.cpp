#include<stdio.h>
#include <assert.h>
#include <math.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<string>
#include <iomanip>
#include <fstream>
#include<map>
#include<vector>
#include<sstream>
#include<algorithm>
#include "api/BamReader.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;

#include "read.h"
#include "snp.h"
#include "hmm.h"


//#define DEBUG

int region = 0;
double errate = 0.005;
int MAX_READS = 10000;
int MAX_READ_LEN = 4000;
int INDEL_RANGE = 15;
long file_pos;

struct read_span {
	long start;
	long end;
};

struct DBSNP {
	int pos;
	char ref;
	char alt;
	int known;
};

vector<SNP*> snp_list;
vector<READ*> reads_list;
vector<struct DBSNP> known_snps;
vector<long int> known_soms;
RefVector refList;

read_span *read_info = new read_span[1000];

vector<string> &split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim))
		elems.push_back(item);
	return elems;
}

// Get a seg fault due to mem corruption error here.
// Fixed by setting \0 at end of line copied from file
vector<string> split(const string &s, char delim,int flag) {
	vector<string> elems;
	return split(s, delim, elems);
}

// Populates refList
// If single chromosome, runs that within specified region
// If region exceeds limits, resets (except start at 0)
int validate_region(const char *file_name, const char *stregion, int *start, int *end)
{
	int flag = 0;
	BamReader reader;
	if(!reader.Open(file_name)) {
		cerr << "Could not open BAM file" << file_name << endl;
		exit(1);
	}

	refList = reader.GetReferenceData();
	if(strcmp(stregion,"ALL")==0) {
		*start = -1;
		*end = -1;
		flag = 1;
	} else {
		string reg = stregion;
		vector<string> contig = split(reg,':',1);
		for(vector<RefData>::iterator rd = refList.begin(); rd!=refList.end(); rd++) {
			if((*rd).RefName==contig[0]) {
				if(contig.size()==1) {
					*start = -1;
					*end = -1;
					flag = 1;
				} else {
					vector<string> pos = split(contig[1],'-',1);
					if(pos.size()!=2||pos[0].empty()||pos[1].empty()) {
						flag = -1;
					} else if(atol(pos[0].c_str())&&atol(pos[1].c_str())) {
						*start = atol(pos[0].c_str());
						*end = atol(pos[1].c_str());
						if(*end>(*rd).RefLength) {
							*end = (*rd).RefLength;
						}
						flag = 1;
					} else {
						flag = -1;
					}
				}
				break;
			}
		}
		if(flag==0)
			flag = -1;
	}
	reader.Close();
	return flag;
}

void SetBamRegion(const char *samfile, BamReader &reader, string &refId, int start, int end)
{
	if(!reader.Open(samfile)) {
		cerr << "Could not open BAM file" << samfile << endl;
		exit(1);
	}

	int RefId = reader.GetReferenceID(refId);
	if(RefId<0) {
		cerr << "Reference sequence " << refId << " not found" << endl;
		exit(1);
	}

	BamRegion bregion(RefId,start,RefId,end);
	// REVISIT: SetRegion returns 0 irrespective of failure or success
	int result = reader.SetRegion(bregion);
	if(result) {
		cerr << "Buggy:Unable to locate region " << refId << ":" << start << "-" << end << endl;
		exit(1);
	}
	return;
}

// Sets file pointer to start of chromosome region specified
// DOES NOT rewind to start location and only moves pointer ahead.
void SetVcfRegion(FILE *file, string &chr)
{
	int flag = 0;
	while(true) {
		char str[10000], word[100];
		string dchr;
		fpos_t prev_pos;

		fgetpos(file,&prev_pos);
		fgets(str,sizeof(str),file);
		str[strlen(str)-1] = '\0';

		if(feof(file))
			break;
		if(str[0]=='#')
			continue;

		if(flag==0) {
			if((str[0]=='c'&&chr.c_str()[0]!='c') || (str[0]!='c'&&chr.c_str()[0]=='c')) {
				cerr << "chr namings in VCF files do not match" << endl;
				exit(1);
			}
			flag = 1;
		}
		int i = 0;
		while(str[i]!='\t')
			word[i] = str[i++];
		word[i] = '\0';
		if(strcmp(word,chr.c_str())!=0)
			continue;
		else {
			fsetpos(file,&prev_pos);
			break;
		}
	}
}

int get_indel_count(vector<CigarOp> &cigar, int limit, int *inend, int *instart, int *delend, int *delstart)
{
	//char c_cigar[1000];
	//strcpy(c_cigar,cigar.c_str());
	//c_cigar[cigar.length()] = '\0';
	int it = 0, len = 0, dlen = 0, mlen = 0, ilen = 0, dmlim1 = 0, dmlim3 = 0, lim1_flag = 0, lim2_flag = 0, dmid = 0, imid = 0;
	int limit2 = limit;
	int limit1 = limit - (INDEL_RANGE+1)/2;
	int limit3 = limit + (INDEL_RANGE-1)/2;

	//int cigar_len = strlen(c_cigar);
	//char *cigar_ptr = c_cigar+it;
	std::vector<CigarOp>::const_iterator cigarIter = cigar.begin();
	std::vector<CigarOp>::const_iterator cigarEnd = cigar.end();
	string cig;

	for ( ; cigarIter != cigarEnd; ++cigarIter) {
		const CigarOp& op = (*cigarIter);

		if(op.Type==Constants::BAM_CIGAR_SOFTCLIP_CHAR||op.Type==Constants::BAM_CIGAR_HARDCLIP_CHAR)
			continue;

		len = op.Length;
		if(op.Type==Constants::BAM_CIGAR_MATCH_CHAR) {
			mlen += len;
			len = 0;
			if(mlen+dlen-len>=limit1) {
				lim1_flag = 1;
				if(mlen+dlen-len>=limit2) {
					lim2_flag = 1;
					if(mlen+dlen-len>=limit3) {
						break;
					}
				}
			}
		} else if(op.Type==Constants::BAM_CIGAR_DEL_CHAR) {
			dlen += len;
			len = 0;
			if(mlen+dlen-len>=limit2&&lim2_flag==0)
				return 9999;
			if(mlen+dlen-len<=limit1)
				*delstart = (dlen);
			*delend = (dlen);
			if(mlen+dlen-len>limit1&&lim1_flag==0) {
				dmlim1 = (mlen+dlen-len)-limit1;
				*delstart = (dlen-dmlim1);
			}
			if(mlen+dlen-len<=limit2)
				dmid = dlen;
			if(mlen+dlen-len>=limit3) {
				dmlim3 = (mlen+dlen-len)-limit3;
				*delend = (dlen-dmlim3);
				break;
			}
		} else if(op.Type==Constants::BAM_CIGAR_INS_CHAR) {
			ilen += len;
			if(mlen+dlen-len+len<limit1)
				*instart = ilen;
			if(mlen+dlen-len+len<limit2)
				imid = ilen;
			len = 0;
			*inend = ilen;
		} else {
			cout << "Invalid character in cigar string " << op.Type  << endl;
		}
	}

	return dmid-imid;
}

int get_read_span(const char *file_name, string &chr, int start, int end)
{
	int iter = 0, flag = 0, span_iter = 0;
	string line, dchr;
	vector<string> read_line;
	char str[20000];

	BamReader reader;
	if(!reader.Open(file_name)) {
		cerr << "Could not open BAM file" << file_name << endl;
		exit(1);
	}

	int RefId = reader.GetReferenceID(chr);
	if(RefId<0) {
		cerr << "Reference sequence " << chr << " not found" << endl;
		exit(1);
	}

	BamRegion bregion(RefId,start,RefId,end);
	// REVISIT: SetRegion returns 0 irrespective of failure or success
	int result = reader.SetRegion(bregion);
	if(result) {
		cerr << "Buggy:Unable to locate region " << chr << ":" << start << "-" << end << endl;
		exit(1);
	}

	BamAlignment al;
	while(reader.GetNextAlignmentCore(al)) {
		if(flag==0) {
			read_info[span_iter].start = al.Position+1;
			flag = 1;
		}
		iter++;
		if(!(iter%MAX_READS)) {
			read_info[span_iter].end = al.Position+1;
			span_iter++;
			read_info[span_iter].start = al.Position+1;
		}
	}

	if(span_iter==0)
		printf("No reads in specified region to read %d:%s\n",chr.c_str(),file_name);
	else 
		read_info[span_iter].end = al.Position+1;

	reader.Close();
	return span_iter+1;
}

void read_known_snp_file(FILE *snp_file, string &chr, int start, int end)
{
	int flag = 0;
	string line;

	while(true) {
		char str[10000];
		string dchr;
		fpos_t prev_pos;

		fgetpos(snp_file,&prev_pos);
		fgets(str,sizeof(str),snp_file);
		str[strlen(str)-1] = '\0';

		if(feof(snp_file))
			break;
		if(str[0]=='#')
			continue;

		line = str;
		vector<string> snp_line = split(line, '\t',0);
		if(flag==0) {
			cout << "Reading dbsnp file..." << endl;
			flag = 1;
		}

		dchr = snp_line[0];
		//int prog = chrcmp(dchr,chr);
		//if(prog>0) {
		if(dchr!=chr) {
			cout << "Finished reading dbsnp..." << endl;
			fsetpos(snp_file,&prev_pos);
			break;
		} //else if(prog<0)
		//	continue;

		int ppos = atol(snp_line[1].c_str());
		int prog = (ppos>end) ? 1 : (ppos<start) ? -1 : 0;
		if(prog>0) {
			fsetpos(snp_file,&prev_pos);
			break;
		} else if(prog<0)
			continue;
		if(flag==1) {
			cout << "Chromosome match. Starting to read now..." << endl;
			flag = 2;
		}
		// Consider hets only. Others are not useful for haplotype calling
		if(snp_line[3]==snp_line[4] || snp_line[3].length()>1 || snp_line[4].length()>1)
			continue;
		//string snpit;
 		//if(snp_line[7].find("PH3;")==string::npos)
		//	snpit = snp_line[1] + "," + snp_line[3] + "," + snp_line[4] + ",1";
		//else
		//	snpit = snp_line[1] + "," + snp_line[3] + "," + snp_line[4] + ",2";

		struct DBSNP snpit;
		snpit.pos = ppos;
		snpit.ref = snp_line[3].c_str()[0];
		snpit.alt = snp_line[4].c_str()[0];
 		snpit.known = snp_line[7].find("PH3;")==string::npos ? 1 : 2;
		
		known_snps.insert(known_snps.end(),snpit);
	}
}

void read_snp_file(FILE *snp_file, string &chr, long snp_end, int start, int end)
{
	string line;

	//vector<string>::iterator vec_start = known_snps.begin();
	vector<struct DBSNP>::iterator vec_start = known_snps.begin();
	while(true) {
		char str[10000];
		fpos_t prev_pos;
		fgetpos(snp_file,&prev_pos);
		fgets(str,sizeof(str),snp_file);
		str[strlen(str)-1] = '\0';
		if(feof(snp_file))
			break;
		if(str[0]=='#')
			continue;
		line = str;
		vector<string> vcf_line = split(line, '\t',0);

		string dchr = vcf_line[0];
		//int prog = chrcmp(dchr,chr);
		//if(prog>0) {
		if(dchr!=chr) {
			fsetpos(snp_file,&prev_pos);
			break;
		} //else if(prog<0)
		//	continue;

		int ppos = atol(vcf_line[1].c_str());
		int prog = (ppos>end||ppos>snp_end) ? 1 : (ppos<start) ? -1 : 0;
		if(prog>0) {
			fsetpos(snp_file,&prev_pos);
			break;
		} else if(prog<0)
			continue;
		// Exclude non-biallelic snps and non-hets
		if(vcf_line[4].find(',')==string::npos) {
			int known = 0;

			// known = 0(error), 1(known), 2(HapMap)
			for(vector<struct DBSNP>::iterator known_snpit = vec_start; known_snpit != known_snps.end(); known_snpit++) {
				struct DBSNP known_line = *known_snpit;
				int pos_s = known_line.pos;
				if(pos_s>ppos)
					break;
				if(pos_s==ppos) {
					known = known_line.known;
					vec_start = known_snpit;
					break;
				}
			}
			// Changing quality score
			SNP *snp = new SNP(vcf_line[0], ppos, vcf_line[3].c_str()[0], vcf_line[4].c_str()[0], known, atof(vcf_line[5].c_str()));
			snp_list.insert(snp_list.end(), snp);
		}
	}
	cout << "SNP MAP size = " << snp_list.size() << endl;
}

int read_sam_files(BamReader &reader)
{
	int i, iter = 0, slen = 0;
	string line;
	BamAlignment al;

	vector<SNP*>::iterator snp_begin = snp_list.begin();
	while(reader.GetNextAlignmentCore(al) && iter<MAX_READS) {
		if(al.MapQuality<=0 || al.IsDuplicate() || al.IsFailedQC())
			continue;

		vector<int> clipSizes;
		vector<int> readPositions;
		vector<int> genomePositions;
		if(al.GetSoftClips(clipSizes, readPositions, genomePositions, false)) {
			int clipsum = 0;
			for(vector<int>::iterator clip=clipSizes.begin(); clip!=clipSizes.end(); clip++) {
				clipsum += (*clip);
			}
			slen = clipsum;
			if(10*clipsum>=3*al.Length)
				continue;
		} else {
			slen = 0;
		}

		al.BuildCharData();
		long start = al.Position+1;
		int length = al.Length;
		READ *read = new READ(start,length); // create new read object
		for(vector<SNP*>::iterator snp_it = snp_begin; snp_it != snp_list.end(); snp_it++) {
			long snp_pos = (*snp_it)->GetPos();
			int known_snp = (*snp_it)->GetKnown();
			if(snp_pos >= start) {
				if(snp_pos <= start+length-1) {
					bool inp,delp;
					char snp_ref = (*snp_it)->GetRef();
					char alt_ref = (*snp_it)->GetAlt();
					char allele;
					int indel_start=0, indel_end=0, inend=0, instart=0, delend=0, delstart=0;
					int qual = 0, type = 2;

					int ndels = get_indel_count(al.CigarData,snp_pos-start+1,&inend,&instart,&delend,&delstart);
					if(ndels==9999) {
						continue;
					} else {
						allele = al.QueryBases[snp_pos + slen - start -ndels];
						qual = al.Qualities[snp_pos + slen - start -ndels];
//cout << al.AlignedBases.length() << ":" << al.QueryBases.length() << ":" << ndels << endl;
						inp = inend>instart ? 1 : 0;
						delp = delend>delstart ? 1 : 0;
#ifdef DEBUG
cout << "Read:" << al.Position+1 << " snp_pos:" << snp_pos << " slen:" << slen << " start:" << start << " ndels:" << ndels << " allele:" << allele << " qual:" << qual << " inend:" << inend << " instart:" << instart << " inp:" << inp << " delend:" << delend << " delstart:" << delstart << " delp:" << delp << endl;
#endif
					}
					if(allele==snp_ref)
						type = 0;
					else if(allele==alt_ref)
						type = 1;

					read->addsnp(*snp_it,allele,qual,inp,delp,known_snp);
					(*snp_it)->append(type,read);
				} else {
					break;
				}
			} else {
				snp_begin = snp_it;
			}
		}
		reads_list.insert(reads_list.end(), read); // add it to the vector
		iter++;
	}

	cout << "number of reads = " << reads_list.size() << endl;
	cout << "number of snps = " << snp_list.size() << endl;
	return iter;
	//return reads_list.size();
}

void runNB(ofstream &stateOutput, ofstream &distanceOutput, long snp_start, long snp_end, int nbSymbols)
{
	char chr[10];

	NaiveBayes(snp_start,snp_end, nbSymbols);
	FindSomaticMutations(snp_start, snp_end);

	vector<SNP*>::iterator snp_list_end = snp_list.end();
	for(vector<SNP*>::iterator snp_it = snp_list.begin(); snp_it != snp_list_end; snp_it++) { // check for snps in vector
		if((*snp_it)->GetKnown())
			continue;
		int status = 0;
		long pos = (*snp_it)->GetPos();
		double som_post[3];
		if(pos<=snp_start)
			continue;
		else if(pos>snp_end)
			break;

		strcpy(chr,(*snp_it)->GetChr());
		chr[strlen((*snp_it)->GetChr())] = '\0';
		som_post[0] = (*snp_it)->GetSomaticPosteriors()[0];
		som_post[1] = (*snp_it)->GetSomaticPosteriors()[1];
		som_post[2] = (*snp_it)->GetSomaticPosteriors()[2];
		if(som_post[1]>=.75)
			status = 1;
		else
			status = 0;
		stateOutput << (*snp_it)->GetChr() << "\t" << (*snp_it)->GetPos() << "\t" << som_post[1] << "\t" << status << endl;
	}
	vector<READ*>::iterator reads_list_end = reads_list.end();
	for(vector<READ*>::iterator read_it = reads_list.begin()+reads_list.size()-(long)nbSymbols+1; read_it != reads_list_end; read_it++) {
		distanceOutput << chr << "\t" << (*read_it)->GetPos() << "\t" << (*read_it)->GetHap() << "\t" << (*read_it)->GetHapProb() << endl;
	}
	return;

/*
	int highHap = (*reads_list.begin())->GetHap();
	double newProb = (*reads_list.begin())->GetHapProb();
	if(highHap==1)
		newProb = newProb;
	else if(highHap==2)
		newProb = 1.0 - newProb;
	else
		cout << "Incorrect Haplotype for last read: " << newProb << endl;

	return newProb;
*/
}

void delete_objects(long snp_start, long snp_end)
{
	for(vector<SNP*>::iterator snp_it = snp_list.begin(); snp_it != snp_list.end();) {
		if((*snp_it)->GetPos() < snp_start) {
			cout << "Unexpected snp below lower limit. Had to be deleted earlier: " << (*snp_it)->GetPos() << endl;
		} else if((*snp_it)->GetPos() >= snp_end) {
			break;
		} else {
			delete (*snp_it);
			snp_list.erase(snp_list.begin());
		}
	}

	for(vector<READ*>::iterator read_it = reads_list.begin(); read_it != reads_list.end();) {
		if((*read_it)->GetPos() < snp_start) {
			delete (*read_it);
			reads_list.erase(reads_list.begin());
		} else if((*read_it)->GetPos()+(*read_it)->GetLen() >= snp_end) {
			break;
		} else {
			delete (*read_it);
			reads_list.erase(reads_list.begin());
		}
	}
}

int main(int argc, char **argv)
{
	int i;
	const char *samfile="";
	const char *vcffile="";
	const char *dbsnpfile="";
	const char *outbase="";
	const char *regionstr="";
	int sam_ck = 0, vcf_ck = 0, dbsnp_ck = 0, out_ck = 0, err_ck = 0, chr_ck = 0, max_reads_ck = 0, max_len_ck = 0, indel_range_ck = 0;
	if(argc>4) {
		int param = 1;
		for(param=1;param<argc;param++) {
			char *token[2];
			token[0] = strtok(argv[param], "=");
			if(token[0]==NULL||(token[1] = strtok(NULL, "=")) == NULL) {
				printf("Ignoring incorrect parameter #%d\n", param);
			} else {
				if(strcmp(token[0],"-sam_file")==0) {
					if(sam_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						sam_ck++;
						samfile = (const char *)token[1];
					}
				} else if(strcmp(token[0],"-vcf_file")==0) {
					if(vcf_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						vcf_ck++;
						vcffile = (const char *)token[1];
					}
				} else if(strcmp(token[0],"-dbsnp")==0) {
					if(dbsnp_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						dbsnp_ck++;
						dbsnpfile = (const char *)token[1];
					}
				} else if(strcmp(token[0],"-out_base")==0) {
					if(out_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						out_ck++;
						outbase = (const char *)token[1];
					}
				} else if(strcmp(token[0],"-err_rate")==0) {
					if(err_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						err_ck++;
						errate = atof((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-chr")==0) {
					if(chr_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						chr_ck++;
						regionstr = (const char *)token[1];
					}
				} else if(strcmp(token[0],"-max_reads")==0) {
					if(max_reads_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						max_reads_ck++;
						MAX_READS = atoi((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-max_read_len")==0) {
					if(max_len_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						max_len_ck++;
						MAX_READ_LEN = atoi((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-indel_range")==0) {
					if(indel_range_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						indel_range_ck++;
						INDEL_RANGE = atoi((const char *)token[1]);
					}
				} else {
					printf("Invalid parameter %s=%s. Ignoring..\n", token[0],token[1]);
				}
			}
		}
	} else {
		printf("Insufficient arguments\n");
		printf("Usage: hapmut -sam_file=<sam file> -vcf_file=<vcf file> -dbsnp=<dbsnp file> -out_base=<output base> [-err_rate=<error rate>] [-chr=<chr:start-end | ALL>] [-max_reads=<max reads at a time>] [-max_read_len=<maximum expected read length>] [-indel_range=<range of search for sequencing artifact>]\n");
		return 1;
	}
	if(sam_ck==0||vcf_ck==0||dbsnp_ck==0||out_ck==0) {
		printf("Insufficient arguments\n");
		printf("Usage: hapmut -sam_file=<sam file> -vcf_file=<vcf file> -dbsnp=<dbsnp file> -out_base=<output base> [-err_rate=<error rate>] [-chr=<chr:start-end | ALL>] [-max_reads=<max reads at a time>] [-max_read_len=<maximum expected read length>] [-indel_range=<range of search for sequencing artifact>]\n");
		return 1;
	}

	int iter=0, start = -1, end = -1;

	// All chromosome(s) from refList are to be run separately now
	if(validate_region(samfile,regionstr,&start,&end)<0) { //chrstrip(string(regionstr));
		printf("Incorrect chr value: %s\n",regionstr);
		return 1;
	}

	char stateOutputName[256];
	char distanceOutputName[256];

	sprintf(stateOutputName,"%s%s", outbase, ".sta");
	sprintf(distanceOutputName,"%s%s", outbase, ".obs");

	double prob1 = 0.99;

	FILE *dbsnp_file = fopen(dbsnpfile, "r");
	if(dbsnp_file==NULL) {
		printf("Cannot open dbsnp file %s\n",dbsnpfile);
		exit(1);
	}
	FILE *snp_file = fopen(vcffile, "r");
	if(snp_file==NULL) {
		printf("Cannot open vcf file %s\n",vcffile);
		exit(1);
	}
	ofstream stateOutput(stateOutputName);
	ofstream distanceOutput(distanceOutputName);

cout << "Entering loop" << endl;

	for(vector<RefData>::iterator refId = refList.begin(); refId!=refList.end(); refId++) {
		string chrname = (*refId).RefName;
		cout << "Reading span file; chromosome " << (*refId).RefName << endl;
		if(start==-1&&end==-1) {
			start = 1;
			end = (*refId).RefLength;
		}
		int span = get_read_span(samfile,(*refId).RefName,start,end);
		if(span==1)
			continue;

		cout << "Reading known snp file; chromosome " << (*refId).RefName << endl;
		SetVcfRegion(dbsnp_file,(*refId).RefName);
		read_known_snp_file(dbsnp_file,(*refId).RefName,start,end);

		BamReader reader;
		// Opens the Bam file and sets the specified region
		SetBamRegion(samfile,reader,(*refId).RefName,start,end);

		while(iter<span) {
			cout << "Reading vcf file; chromosome " << (*refId).RefName << endl;
			SetVcfRegion(snp_file,(*refId).RefName);
			read_snp_file(snp_file, (*refId).RefName, read_info[iter].end + MAX_READ_LEN, start, end); // Reads specified range or from pos 1
			cout << "Reading sam file; chromosome " << (*refId).RefName << endl;
			int read_ct = read_sam_files(reader); // Reads 10000 at a time
			if(iter==0) {
				distanceOutput << (*refId).RefName << "\t" << reads_list[0]->GetPos() << "\t" << reads_list[0]->GetHap() << "\t" << reads_list[0]->GetHapProb() << endl;
			} else {
				read_ct++;
			}
			// Runs haplotype calling on [start,end] start position range of reads
			// Runs mutation calling on (start,end]
			runNB(stateOutput, distanceOutput, read_info[iter].start, read_info[iter].end, read_ct);
			// Delete reads [start,end.end<end]
			// Delete snps [start,end)
			delete_objects(read_info[iter].start, read_info[iter].end);
			iter++;
		}
		reader.Close();
	}

	distanceOutput.close();
	stateOutput.close();
	fclose(snp_file);
	fclose(dbsnp_file);
	return 0;
}

/*
int chrstrip(string region)
{
	int index = 0;
	const char *regionstr = region.c_str();

	if(regionstr[0]=='c'&&regionstr[1]=='h'&&regionstr[2]=='r')
		index = 3;

	if(regionstr[index]=='X')
		return 23;
	else if(regionstr[index]=='Y')
		return 24;
	//else if(regionstr[index]=='M'&&regionstr[index+1]=='T')
	else if(regionstr[index]=='M')
		return 25;
	else if(atoi(regionstr+index)<1||atoi(regionstr+index)>22)
		return -1;
	else {
		//const char *strindex = regionstr+index;
		return atoi(regionstr+index);
	}
}

int chrcmp(string &dchr, string &chr)
{
	int p1 = chrstrip(dchr);
	int p2 = chrstrip(chr);

	if(p1>p2)
		return 1;
	else if(p1<p2)
		return -1;
	else
		return 0;
}
*/

