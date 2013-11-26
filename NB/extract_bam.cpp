#include "api/BamReader.h"
#include "api/BamWriter.h"
#include <stdlib.h>
#include <stdio.h>
#include<vector>

using namespace std;
using namespace BamTools;

//  g++ -g -Wall -I/home/singhann/tools/bamtools-master/include -L/home/singhann/tools/bamtools-master/lib extract_bam.cpp -lm -lbamtools -o extract_bam

//string inputFile("/home/singhann/files/30/Agilent_Run_30_35_HCC1143_BL-Normal.cleaned.bam.filtered.bam.chr1.bam.reheader.bam");
//string outputFile("/home/singhann/files/30/Agilent_Run_30_35_HCC1143_BL-Normal.cleaned.bam.filtered.bam.chr1.reheader.bam");


int main(int argc, char **argv) {

string inputFile(argv[1]);
//string outputFile(argv[2]);
BamReader reader;
if(!reader.Open(inputFile)) {
	cerr << "Could not open BAM file" << endl;
	return 1;
}

const SamHeader header = reader.GetHeader();
const RefVector references = reader.GetReferenceData();

BamWriter writer;
//if (!writer.Open(outputFile,header,references)) {
//	cerr << "Could not open output BAM file" << endl;
//	return 1;
//}

BamAlignment al;
while( reader.GetNextAlignmentCore(al)) {
	string str_type;
	int int_type;
	al.BuildCharData();
	std::vector<CigarOp>::const_iterator cigarIter = al.CigarData.begin();
	std::vector<CigarOp>::const_iterator cigarEnd = al.CigarData.end();
	string cig;
	std::vector<string> tags = al.GetTagNames();

	cout << "Name: " << al.Name << endl;
	cout << "Length: " << al.Length << endl;
	cout << "Position: " << al.Position << endl;
	cout << "Allele: " << endl;
	cout << "Quality: " << endl;
}
reader.Close();
return 0;
}
/*
	for ( ; cigarIter != cigarEnd; ++cigarIter) {
		const CigarOp& op = (*cigarIter);
		char ciglen[10];
		sprintf(ciglen,"%d",op.Length);
		
		switch (op.Type) {
			case Constants::BAM_CIGAR_DEL_CHAR :
				cig += ciglen;
				cig += "D";
				break;
			case Constants::BAM_CIGAR_MATCH_CHAR :
				cig += ciglen;
				cig += "M";
				break;
			case Constants::BAM_CIGAR_MISMATCH_CHAR :
				cig += ciglen;
				cig += "X";
				break;
			case Constants::BAM_CIGAR_REFSKIP_CHAR :
				cig += ciglen;
				cig += "N";
				break;
			case Constants::BAM_CIGAR_SEQMATCH_CHAR :
				cig += ciglen;
				cig += "=";
				break;
			case Constants::BAM_CIGAR_INS_CHAR :
				cig += ciglen;
				cig += "I";
				break;
			case Constants::BAM_CIGAR_SOFTCLIP_CHAR :
				cig += ciglen;
				cig += "S";
				break;
			case Constants::BAM_CIGAR_HARDCLIP_CHAR :
				cig += ciglen;
				cig += "H";
				break;
			case Constants::BAM_CIGAR_PAD_CHAR :
				cig += ciglen;
				cig += "P";
				break;
			default:
				cig += "*";
				break;
		}
	}
	if(cig.empty())
		cig += "*";
	string mateid;
	switch(al.MateRefID) {
		case 0:
			mateid = "=";
			break;
		case 22:
			mateid = "X";
			break;
		case 23:
			mateid = "Y";
			break;
		case 24:
			mateid = "MT";
			break;
		default:
			char mate[10];
			sprintf(mate,"%d",al.MateRefID+1);
			mateid = mate;
			break;
	}

	cout << al.Name.c_str() <<"\t" << al.AlignmentFlag << "\t1\t" << al.Position+1 << "\t" << al.MapQuality <<"\t" << cig << "\t" << mateid << "\t" << al.MatePosition+1 << "\t" << al.InsertSize << "\t" << al.QueryBases.c_str() << "\t" << al.Qualities.c_str();
	for(std::vector<string>::iterator it = tags.begin(); it!=tags.end();++it) {
		if((*it).compare("AM")==0||(*it).compare("AS")==0||(*it).compare("CM")==0||(*it).compare("CP")==0||(*it).compare("FI")==0||(*it).compare("H0")==0||(*it).compare("H1")==0||(*it).compare("H2")==0||(*it).compare("HI")==0||(*it).compare("IH")==0||(*it).compare("MQ")==0||(*it).compare("NH")==0||(*it).compare("NM")==0||(*it).compare("OP")==0||(*it).compare("PQ")==0||(*it).compare("SM")==0||(*it).compare("TC")==0||(*it).compare("UQ")==0||(*it).compare("X0")==0||(*it).compare("X1")==0||(*it).compare("XG")==0||(*it).compare("XM")==0||(*it).compare("XO")==0||(*it).compare("XN")==0||(*it).compare("XC")==0) {
			al.GetTag(*it,int_type);
			cout << "\t" << (*it).c_str() << ":i:" << int_type;
		} else if((*it).compare("BC")==0||(*it).compare("BQ")==0||(*it).compare("CC")==0||(*it).compare("CQ")==0||(*it).compare("CS")==0||(*it).compare("E2")==0||(*it).compare("FS")==0||(*it).compare("LB")==0||(*it).compare("MD")==0||(*it).compare("OQ")==0||(*it).compare("OC")==0||(*it).compare("PG")==0||(*it).compare("PU")==0||(*it).compare("Q2")==0||(*it).compare("R2")==0||(*it).compare("RG")==0||(*it).compare("U2")==0||(*it).compare("XA")==0) {
			al.GetTag(*it,str_type);
			cout << "\t" << (*it).c_str() << ":Z:" << str_type;
		} else if((*it).compare("XT")==0) {
			al.GetTag(*it,int_type);
			char char_type = int_type;
			cout << "\t" << (*it).c_str() << ":A:" << char_type;
		} else {
			cerr << "Invalid Tag Type " << (*it).c_str() << endl;
		}
	}
	cout << endl;
	//writer.SaveAlignment(al);
}

reader.Close();
writer.Close();

}
*/
