#!usr/bin/env python3
import pysam 
import logging 
import sys

sample=sys.argv[1]
logging.basicConfig(format='%(message)s', level=logging.DEBUG)
bam_fn = '/working/joint_projects/Patch_CRC/clemencB/GandhiFL2.0/1.BWAKIT/'+sample+'.aln.bam'
with pysam.AlignmentFile(bam_fn,'rb') as bam : 

	hla_genes=[ ref for ref in bam.references if ref.startswith('HLA-')]
	alt_contigs=[ ref for ref in bam.references if ref.startswith('chr6_GL000')]+['chr6_KI270758v1_alt']

	## Retrieve the query names
	query_names = set()
	for read in bam.fetch('*'):
		query_names.add(read.query_name)		
	for ref in hla_genes + alt_contigs :
		#logging.debug("Getting reads for HLA gene: {}".format(ref))
		for read in bam.fetch(ref):
			#logging.debug("Found alignment {}".format(read))
			query_names.add(read.query_name)
		

	for read in bam.fetch(region='chr6:28510120-33480577'):
		query_names.add(read.query_name)

	## Create the new bam
	with pysam.AlignmentFile("HLAbam_"+sample+".aln.bam","wb",template=bam) as output:
		for read in bam.fetch('.'):
			if read.query_name in query_names:
				output.write(read)

