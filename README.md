# haplotype_freq_filter

Read COPYRIGHT file for copyright information. 
Modification of code by K. Thornton to filter on individual haplotypes. 
Rather than filter on every SNP as in the original code, this filters only the FIRST SNP of an ms simulation, keeping haplotypes that meet frequency criteria based on the first SNP.i

## Usage:

* To filter based on major allele frequency use "-m freq";
* To filter based on derived allele frequency use "-d freq";
* Use -r to print the numbers of the haplotypes not filtered (i.e. kept) to standard error.
