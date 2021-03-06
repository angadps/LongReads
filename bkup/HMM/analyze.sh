#!/bin/sh
#$ -cwd

prefix=$1
thresh=`echo "(${SGE_TASK_ID})/100.0" | bc -l`

awk -v d=$thresh '{
	if($3==2) {
		if(($6+$7)>=d)
			trueab++;
		else
			abh++;
		if($9+$10>=d)
			strueab++;
		else
			sabh++;
	} else if($3==0) {
		if($6+$7<d)
			falseab++;
		else
			hab++;
		if($9+$10<d)
			sfalseab++;
		else
			shab++;
	}
}
END{
	absens = trueab/(trueab+abh);
	acc = trueab/(trueab+hab);
	spcab = falseab/(falseab+hab);
	sabsens = strueab/(strueab+sabh);
	sacc = strueab/(strueab+shab);
	sspcab = sfalseab/(sfalseab+shab);
	print 1-acc, "\t", 1-spcab, "\t", absens, "\t", d, "\t", 1-sacc, "\t", 1-sspcab, "\t", sabsens;
}' ${prefix}.sta > ${prefix}.b${SGE_TASK_ID}

awk -v d=$thresh '{
	if($13==1) {
		if($14>=d)
			truep++;
		else
			falsen++;
		if($9>=d)
			struep++;
		else
			sfalsen++;
	} else if($13==0) {
		if($14>=d)
			falsep++;
		else
			truen++;
		if($9>=d)
			sfalsep++;
		else
			struen++;
	}
}
END{
	sens = truep / (truep+falsen);
	acc = falsep / (truep+falsep);
	spc = truen / (truen+falsep);
	ssens = struep / (struep+sfalsen);
	sacc = struep / (struep+sfalsep);
	sspc = struen / (struen+sfalsep);
	print 1-acc, "\t", 1-spc, "\t", sens, "\t", d, "\t", 1-sacc, "\t", 1-sspc, "\t", ssens;
}' ${prefix}.sta > ${prefix}.c${SGE_TASK_ID}

# Code to obtain sensitivity to accuracy ratio for different cut-off for heterozygous snps

#awk -v d=$thresh '{ if($3==2) { if($6>=d) trueab++; else abh++; if($9>=d) strueab++; else sabh++; } else if($3==0) { if($6<d) falseab++; else hab++; if($9<d) sfalseab++; else shab++; } } END{absens = trueab/(trueab+abh); acc = trueab/(trueab+hab); spcab = falseab/(falseab+hab); sabsens = strueab/(strueab+sabh); sacc = strueab/(strueab+shab); sspcab = sfalseab/(sfalseab+shab); print 1-acc, "\t", 1-spcab, "\t", absens, "\t", d, "\t", 1-sacc, "\t", 1-sspcab, "\t", sabsens;}' ${prefix}.sta > ${prefix}.b${SGE_TASK_ID}

