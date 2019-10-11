#!/usr/bin/env python3

import os
import sys
import sqlite3
import pandas as pd
from AddJunctionsToDatabase import connectToDB, commitAndClose

conn, cur = connectToDB()

MIN_READ = 5
MIN_NORM_READ = 0.05
MAX_GTEX_SEEN = 5

def countGTEX():
	cur.execute('select count(*) from SAMPLE_REF where type = 0;') # 0 = GTEX, 1 = PATIENT
	return cur.fetchone()[0]

def countPatients():
	cur.execute('select count(*) from SAMPLE_REF where type = 1;') # 0 = GTEX, 1 = PATIENT
	return cur.fetchone()[0]

def sampleSpecificJunctions(sample, min_read, min_norm_read):

	"""
	Generates a file with 
		1) junctions seen in a sample and not seen in any GTEx samples 
		2) with a read count equal to or greater than the specified minimum read count. 

	Note:
		1) Function does not discriminate against junctions seen in other patient samples.
		2) Function does not work for GTEx samples as the database query relies on n_gtex_seen being 0.

	The query provides information about read counts of only the one specified sample
	in the columns 'sample:read_count' and 'sample:norm_read_count', 

		Ex. 1:100-200 PATIENT2:344

	This occurs because we are joining on only one sample in the database as opposed
	to all.

	Args:
		sample, the name of the sample file you want to investigate, must include .bam extension
		min_read, the minimum number of reads a junction must have
		min_norm_read, the minimum normalized read count a junction must have or NULL

	    not reporting junctions with NULL norm_count (mostly NONE annotated)
	"""

	count = str(countGTEX())
	output = '_'.join([sample, 'specific', 'rc' + str(min_read), ('norm_rc' + str(min_norm_read)), 'n_gtex_' + count])

	df = pd.read_sql_query('''select gene_ref.gene,
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop) as pos,
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.total_patient_read_count as read_count,
		junction_counts.norm_read_count as norm_read_count,
		junction_ref.n_gtex_seen as n_gtex_seen,
		junction_ref.total_gtex_read_count as total_gtex_read_count
		from junction_counts, sample_ref, junction_ref, gene_ref
		where 
			sample_ref.sample_name = "{}" and
			junction_counts.read_count >= {} and
			junction_counts.norm_read_count >= {} and junction_counts.norm_read_count!='NULL' and
			junction_ref.n_gtex_seen <= {} and
			sample_ref.rowid=junction_counts.bam_id and
			junction_counts.junction_id = junction_ref.rowid and
			junction_ref.rowid = gene_ref.junction_id;'''.format(sample, min_read, min_norm_read, MAX_GTEX_SEEN), conn)

	df.to_csv(output, index=False, sep="\t")

def customSampleSpecificJunctions(sample, min_read, min_norm_read, max_n_gtex_seen, max_total_gtex_reads):

	"""
	Generates a text file using a query in which you can discover junctions specific to a sample
	with parameters for appearances in GTEx samples

	The query provides information about read counts of only the one specified sample
	in the columns 'sample:read_count' and 'sample:norm_read_count', 

		Ex. 1:100-200 PATIENT2:344

	This occurs because we are joining on only one sample in the database as opposed
	to all.

	Args:
		cur, a cursor to a connection to a database
		sample, the name of the sample file you want to investigate, must include .bam extension
		min_read, the minimum number of reads a junction must have
		min_norm_read, the minimum normalized read count a junction must have or NULL
		max_n_gtex_seen, the maximum number of gtex samples a junction can appear in
		max_total_gtex_reads, the maximum total read count for a junction in GTEx samples

	Returns:
	    None

	Raises:
	    None
	"""

	if not max_n_gtex_seen:
		max_n_gtex_seen = 0

	if not max_total_gtex_reads:
		max_total_gtex_reads = 0

	if not min_read:
		min_read = 0

	output = '_'.join([str(sample), ('rc' + str(min_read)), ('norm_rc' + str(min_norm_read)), ('maxGTEX' + str(max_n_gtex_seen)), ('maxGTEXrc' + str(max_total_gtex_reads))])

	# does not report events with norm_count==NULL
	df = pd.read_sql_query('''select group_concat(gene_ref.gene),
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop) as pos,
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.n_gtex_seen, 
		junction_ref.n_patients_seen,
		junction_ref.total_patient_read_count,
		junction_ref.total_gtex_read_count,
		junction_ref.total_read_count,
		group_concat(distinct junction_counts.read_count||':'||sample_ref.sample_name),
		group_concat(distinct junction_counts.norm_read_count||':'||sample_ref.sample_name)
		from 
		sample_ref inner join junction_counts on sample_ref.ROWID = junction_counts.bam_id 
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		inner join gene_ref on junction_ref.rowid = gene_ref.junction_id 
		where
		sample_ref.sample_name = "{}" and
		junction_counts.read_count >= {} and
		junction_ref.n_gtex_seen <= {} and
		junction_ref.total_gtex_read_count <= {} and
		junction_counts.norm_read_count >= {}
		group by 
		junction_ref.chromosome, junction_ref.start, junction_ref.stop;'''.format(sample, min_read, max_n_gtex_seen, max_total_gtex_reads, min_norm_read), conn)

	df.to_csv(output, index=False, sep="\t")

def sampleRatio(sample, junction_ratio):

	"""
	Dumps all junction information seen in all samples to a text file.

	The query provides information about read counts of a junction across all samples
	unlike the other queries in the columns 'sample:read_count' and 'sample:norm_read_count', 

		Ex. 1:100-200 GTEx1:20,GTEx3:211,PATIENT2:344

	This occurs because we are joining and grouping all sample names in the database as opposed
	to just one name.

	Args:
		sample, the name of the sample file you want to investigate, must include .bam extension
		junction_ratio, the minimum factor gain seen in the sample compared to the average GTEx sample

	Returns:
	    None

	Raises:
	    None
	"""

	output = '{}.junctions.{}x.control.tsv'.format(sample, junction_ratio)

	df = pd.read_sql_query('''select gene_ref.gene,
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop) as pos,
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation,
		junction_counts.read_count,
		junction_counts.norm_read_count,
		(CAST(junction_ref.total_patient_read_count as FLOAT) / junction_ref.n_patients_seen) as avg_patient_read_count,
		(CAST(junction_ref.total_gtex_read_count as FLOAT) / junction_ref.n_gtex_seen) as avg_gtex_read_count,
		(
			select MAX(j_counts.read_count) 
			from sample_ref s_ref, junction_ref j_ref, junction_counts j_counts 
			where s_ref.type=0 and 
			j_ref.chromosome=junction_ref.chromosome and
			j_ref.start=junction_ref.start and
			j_ref.stop=junction_ref.stop and
			s_ref.rowid=j_counts.bam_id and 
			j_counts.junction_id=j_ref.rowid
		) max_gtex_read_count
		from junction_counts, sample_ref, junction_ref, gene_ref
		where 
			sample_ref.sample_name = "{}" and
			junction_counts.read_count >= {} and
			junction_counts.norm_read_count >= {} and
			junction_counts.read_count >= {} * avg_gtex_read_count and
			junction_counts.norm_read_count!='NULL' and
			sample_ref.rowid=junction_counts.bam_id and
			junction_counts.junction_id = junction_ref.rowid and
			junction_ref.rowid = gene_ref.junction_id;'''.format(sample, MIN_READ, MIN_NORM_READ, junction_ratio), conn)

	df.to_csv(output, index=False, sep="\t")

def printSamplesInDB():

	"""
	Prints to stdout all samples in the database and their experiment type

	Args:
		cur, a cursor to a connection to a database

	Returns:
	    None

	Raises:
	    None
	"""

	cur.execute('''select sample_name,
	case
		when type = 0 then 'CONTROL'
		when type = 1 then 'PATIENT'
	END
	from SAMPLE_REF;''')

	for line in cur.fetchall():
		print('\t'.join(str(i) for i in line))

def printAllJunctions():

	"""
	Dumps all junction information seen in all samples to a text file.

	The query provides information about read counts of a junction across all samples
	unlike the other queries in the columns 'sample:read_count' and 'sample:norm_read_count', 

		Ex. 1:100-200 GTEx1:20,GTEx3:211,PATIENT2:344

	This occurs because we are joining and grouping all sample names in the database as opposed
	to just one name.

	Args:
		cur, a cursor to a connection to a database

	Returns:
	    None

	Raises:
	    None
	"""

	output = 'all_junctions_n_gtex_{}_n_patients_{}'.format(str(countGTEX()), str(countPatients()))

	df = pd.read_sql_query('''select group_concat(gene_ref.gene) as genes,
		(junction_ref.chromosome||':'||junction_ref.start||'-'||junction_ref.stop) as pos,
		case 
			when junction_ref.gencode_annotation = 0 then 'NONE'
			when junction_ref.gencode_annotation = 1 then 'START'
			when junction_ref.gencode_annotation = 2 then 'STOP'
			when junction_ref.gencode_annotation = 3 then 'BOTH'
			when junction_ref.gencode_annotation = 4 then 'EXON_SKIP'
		END as annotation, 
		junction_ref.n_gtex_seen, 
		junction_ref.n_patients_seen,
		junction_ref.total_patient_read_count,
		junction_ref.total_gtex_read_count,
		junction_ref.total_read_count,
		group_concat(distinct junction_counts.read_count||':'||sample_ref.sample_name) as sample_read_count,
		group_concat(distinct junction_counts.norm_read_count||':'||sample_ref.sample_name) as norm_sample_read_count
		from 
		sample_ref inner join junction_counts on sample_ref.ROWID = junction_counts.bam_id 
		inner join junction_ref on junction_counts.junction_id = junction_ref.rowid 
		inner join gene_ref on junction_ref.rowid = gene_ref.junction_id
		where junction_ref.total_read_count > 0
		group by 
		junction_ref.chromosome, junction_ref.start, junction_ref.stop;''', conn)

	df.to_csv(output, index=False, sep="\t")

if __name__=="__main__":
	
	if sys.argv[1] == '--printsamples':
		printSamplesInDB()
	elif sys.argv[1] == '--sample':
		print(sys.argv[4])
		sampleSpecificJunctions(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]))
	elif sys.argv[1] == '--sample_ratio':
		sampleRatio(sys.argv[2], float(sys.argv[3]))
	elif sys.argv[1] == '--custom_counts':
		customSampleSpecificJunctions(sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6])
	elif sys.argv[1] == '--all':
		printAllJunctions()
	else:
		print('Invalid option. Use one of the following:')
		print('--printsamples')
		print('--sample	[SAMPLE]	[MIN_READ]	[MIN_NORM_READ]')
		print('--sample_ratio [SAMPLE] [JUNCTION_RATIO]')
		print('--custom_counts [SAMPLE] [MIN_READ] [MIN_NORM_READ]	[MAX_N_GTEX_SEEN] [MAX_TOTAL_GTEX_READS]')
		print('--all')

	commitAndClose(conn)
