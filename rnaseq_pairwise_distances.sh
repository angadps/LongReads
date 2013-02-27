#!/bin/sh
#$ -cwd

#Script to get pairwise snp distances from an exome sequence sample

listf=$1

awk -v file="/ifs/scratch/c2b2/ys_lab/aps2157/shared/SeqCap_EZ_Exome_v2.bed" 'BEGIN{	i=1;
	while((getline line<file)>0) {
		split(line,pos,"\t");
		list[i,1]=pos[1];list[i,2]=pos[2];list[i,3]=pos[3];
		i++;
	}
	chr=list[1,1];start=list[1,2];end=list[1,3];
	ind=1;prev=0;
}
{
	dist=0;
	if("chr"$1!=chr) {
		while("chr"$1!=chr) {
			ind++; prev=0;
			chr=list[ind,1];start=list[ind,2];end=list[ind,3];
		}
	}
	if($2<end) {
		dist += $2-prev;
		prev=$2;
	} else {
		while($2>end) {
			ind++;
			if(prev>start)
				dist += end - prev + 1;
			else
				dist += end - start + 1;
			chr=list[ind,1];start=list[ind,2];end=list[ind,3];
		}
		dist += $2 - start;
		prev=$2;
	}
	print dist;
}' $listf

