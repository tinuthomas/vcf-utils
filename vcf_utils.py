#!/use/bin/env python

import os 
import re
import sys 
import gzip 
import bz2
from collections import defaultdict

def return_file_handle(fname):
	"""
	Function open the file and returns the handler 
	"""
	try:
		if os.path.splitext(fname)[1] == ".gz":
			FH = gzip.open(fname, 'rb')
		elif os.path.splitext(fname)[1] == ".bz2":
			FH = bz2.BZ2File(fname, 'rb')
		else:
			FH = open(fname, 'rU')
	except Exception as error:
		sys.exit(error)
	return FH

def tab_delimit_vcf(fin):
	for line in fin :
		line = line.strip().split('\t')
        	if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
			data = line
		elif line[0] == '#CHROM':
			header = line
		elif line[0].startswith('##'):
			meta = line
			 

def VCFtoBED():
	"""
	Makes a BED file from a VCF
	"""
	try:
		vcffile = sys.argv[1]
	except:
		print __doc__
		sys.exit(-1)
	fin = return_file_handle(vciffile)
	
	for line in fin:
		line = line.strip().split('\t')
		if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
			print '\t'.join([line[0], str(int(line[1]) - 1), line[1]])
	fin.close()

def VCF_AF_filter():
	"""
	FunctionName :   VCF_AF_filter
	Author      :   Tinu Thomas
	Date        :   10/25/2013
	Purpose     :   Given a VCF file, filters and gives a VCF with variants greater or lesser than specified AF value by the user
	Usage       :   python AF_filter.py option VCF AF_value >OutputVCF
        	        Option      -   1 or 2
                 	   1) Filter variants with AF greater than AF specified by the user	 eg: python AF_filter.py Input.vcf  0.01 0
	                    2) Filter variants with AF value less than AF specified by the user	 eg: python AF_filter.py Input.vcf 0.01 0
        	            3) Filter variants within a range of AF values specified by the user eg: python AF_filter.py Input.vcf 0.01 0.05
	                AF_value    -   To be specified by the user(anything between 0 and 1)
        	        VCF         -   Input VCF file
			For Option 3   python AF_filter.py option VCF AF_lowerlimit AF_upperlimit >OutputVCF
	"""
	try:
	    option = float(sys.argv[1])
	    inVCF = sys.argv[2]
	    af_check = float(sys.argv[3])
	    upper_limit = float(sys.argv[4])
	except:
	    print __doc__
	    sys.exit(-1)

	if type(option) is not float or type(af_check) is not float:	# Checking is the correct option is entered
	    print "Option and AF must be numbers\nOption      -   1 or 2\n1) Filter variants with AF greater than AF specified by the user\n2) Filter variants with AF value less than AF specified by the user\nAF_value    -   To be specified by the user(anything between 0 and 1)\n"
	    sys.exit(-1)
	elif af_check < 0 or af_check > 1 :				# Checking if the permitable values are entered for Allele Frequency
	    print 'AF values should be only anything from 0 to 1'
	    sys.exit(-1)
	if float(option) == 1 :						# If option is 1, filters for variants with AF greater than the mentioned user entered
	    fin = return_file_handle(inVCF)
	    for line in fin:
	        line = line.strip('\n\r').split('\t')
        	if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
	            flag = 0
        	    for val in line[7].split(';'):
                	if '=' in val and val.split('=')[0] == 'AF' and ',' not in val.split('=')[1]  :
	                    if float(val.split('=')[1]) > af_check:
        	                flag = 1
                	        break
	                elif '=' in val and val.split('=')[0] == 'AF' and ',' in val.split('=')[1]:
        	            for af_val in val.split('=')[1].split(','):
                	        if float(af_val) > af_check:
                        	    flag = 1
	                            break
        	    if flag == 1 :
                	print '\t'.join(line)
	        else:
        	        print '\t'.join(line)
	    fin.close()
	elif float(option) == 2:				# If option is 2, filters for variants with AF less than the user mentioned AF
	    fin = return_file_handle(inVCF)
	    for line in fin:
        	line = line.strip('\n\r').split('\t')
	        if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
        	    flag = 0
	            for val in line[7].split(';'):
        	        if '=' in val and val.split('=')[0] == 'AF' and ',' not in val.split('=')[1]  :
                	    if float(val.split('=')[1]) < af_check:
                        	flag = 1
	                        break
        	        elif '=' in val and val.split('=')[0] == 'AF' and ',' in val.split('=')[1]:
                	    for af_val in val.split('=')[1].split(','):
                        	if float(af_val) < af_check:
	                            flag = 1
        	                    break
	            if flag == 1 :
        	        print '\t'.join(line)
	        else:
        	        print '\t'.join(line)
	    fin.close()
	elif float(option) == 3:				# If option is 3, filters for variants with AF in the range specified by the user
	    fin = return_file_handle(inVCF)
	    for line in fin:
	        line = line.strip('\n\r').split('\t')
        	if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
	            flag = 0
        	    for val in line[7].split(';'):
                	if '=' in val and val.split('=')[0] == 'AF' and ',' not in val.split('=')[1]  :
	                    if float(val.split('=')[1]) >= af_check and float(val.split('=')[1]) <= upper_limit :
        	                flag = 1
                	        break
	                elif '=' in val and val.split('=')[0] == 'AF' and ',' in val.split('=')[1]:
        	            for af_val in val.split('=')[1].split(','):
                	        if float(af_val) >= af_check and float(af_val) <= upper_limit:
                        	    flag = 1
	                            break
        	    if flag == 1 :
                	print '\t'.join(line)
	        else:
        	        print '\t'.join(line)
	else:
	    print "Please give the correct option 1 or 2\nOption      -   1 or 2\n1) Filter variants with AF greater than AF specified by the user\n2) Filter variants with AF value less than AF specified by the user\nAF_value    -   To be specified by the user\nVCF         -   Input VCF file"
	    sys.exit(-1)

def VCF_Info_filter():
	
	"""
	FunctionName :   VCF_Info_filer
	Author      :   Tinu Thomas
	Date        :   10/25/2013
	Purpose     :   Given a VCF file, filters and gives a VCF with variants greater or lesser than specified INFO field value by the user
	Usage       :   python  VCF_InfoField_Filter.py   Input.vcf  INFO_field   Value
        	        Option      -   1
                	    1) Filter variants with INFO field value less than the value specified by the user
			inVCF		-	Input VCF
			info_field	-	The INFO field in the VCF on which the filtering needs to be done
			value		-	Value of the INFO field by which the filtering needs to be done
	"""

	try:
		inVCF = sys.argv[1]
		info_field = sys.argv[2]
		value = float(sys.argv[3])
	except:
		print __doc__
		sys.exit(-1)
	info_pos = -9
	fin = return_file_handle(inVCF)
	for line in fin:
		line = line.strip('\n\r').split('\t')
		if line[0] == '#CHROM':		
			print '\t'.join(line)
		elif re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
			flag = 0
			for k,x in enumerate(line[7].split(';')):
				if x.split('=')[0] == info_field:
					if float(x.split('=')[1]) < value:
						flag = 1
						break
				
			if flag == 1 :
				print '\t'.join(line)
		else:
			print '\t'.join(line)


def VariantType_Filter():
	"""
	Given a VCF with VariantType INFO field, gives the statistics of SNPs and INDELs
	"""
	try:
		inVCF = sys.argv[1]
	except:
		print __doc__
		sys.exit(-1)
	ins_ct = 0
	del_ct = 0
	snp_ct = 0	
	var_ct = 0
	fin = return_file_handle(inVCF)
	for line in fin:
		line = line.strip().split('\t')
		if re.search(r'^chr(\d+|X|Y)|^(\d+|X|Y)',line[0]):
			var_ct += 1
			for val in line[7].split(';'):
				if val.split('=')[0] == 'VariantType' and val.split('=')[1].startswith('INSERTION'):
					ins_ct += 1
				elif val.split('=')[0] == 'VariantType' and val.split('=')[1].startswith('DELETION'):
					del_ct += 1
				elif val.split('=')[0] == 'VariantType' and val.split('=')[1] == 'SNP':
					snp_ct += 1
				
	fin.close()
	print "Variant_count ", var_ct
	print "SNP_count ", snp_ct
	print "INSERTION_count ", ins_ct


def VCF_SNP_filer():
	"""
	Given a VCF with or without any annotations, filters out bi-allelic SNPs
	"""

	file = sys.argv[1]
	fin = return_file_handle(file)
	for line in fin:
		line = line.strip().split('\t')
		if re.search(r'^(\d+|X|Y)', line[0]):
			if ',' not in line[4] and len(line[4]) == 1:
				print '\t'.join(line)
		else:	
			print '\t'.join(line)
	fin.close()	


def vcf_parser(input_vcf, format_field):
    """
    VCF main parsing program
	Program to convert VCF files to TEXT format

	Usage: python vcf_to_txt.py Input.vcf format_field >Input.vcf.txt 
	       format_field can only be GT, DP, GQ, PL, MQ0 or AD
    """
    if format_field != 'GT' and format_field != 'GQ' and format_field != 'DP' and format_field != 'PL' and format_field != 'AD' and format_field != 'MQ0':
        print 'This format field cannot be split'
        print __doc__
        sys.exit(-1)
    fin = return_file_handle(input_vcf)
    info_ct = format_ct = variant_ct = 0
    info_fields = {}
    info_header = []
    formatname_list = []
    for line in fin:
        line = line.strip('\n\r')
        if line.startswith('##INFO'):
            line = line.split(',')
            infoname = re.search(r'^##INFO=<ID=(.+)',line[0])
            infotype = re.search(r'^Type=(.*)',line[2])
            info_fields[infoname.group(1)] = infotype.group(1)        
        elif line.startswith('##FORMAT'):
            line = line.split(',')
            formatname = re.search(r'^##FORMAT=<ID=(.+)',line[0])
            formatname_list.append(formatname.group(1))
        elif line.startswith('#CHROM'):
            line = line.split('\t')
            if len(line) > 8 and format_field not in formatname_list:
                print "VCF doesnot contain ",format_field," format field"
                sys.exit(-1)
            for head in sorted(info_fields.keys()):
                info_header.append(head)
            if len(line) > 8:
                sample_names = []
                for gt in line[9:]:
                    sample_names.append(gt + '_' + format_field)
                #print '\t'.join( ['chr'+ ':'.join(line[0:2]),'VAR_TYPE' ]+ line[0:7] + ['minAF'] + info_header + sample_names)
                print '\t'.join( line[0:7] + ['chr'+ ':'.join(line[0:2]),'VAR_TYPE' ]+ ['minAF'] + info_header + sample_names)
            else: 
                #print '\t'.join( ['chr'+ ':'.join(line[0:2]), 'VAR_TYPE']+ line[0:7] + ['minAF'] + info_header)
                print '\t'.join( line[0:7] + ['chr'+ ':'.join(line[0:2]), 'VAR_TYPE']+ ['minAF'] + info_header)
        elif re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)',line):
            line = line.split('\t')
            info_values, AF_final = info_parser(info_header,info_fields,line[7])           
            if ',' in line[4]:
                var_type = 'MAV'
            else:
                if len(line[3]) > 1 or len(line[4]) > 1:
                    var_type = 'INDEL'
                elif len(line[3]) == 1 and len(line[4]) == 1:
                    var_type = 'SNP'
            if len(line) > 8 and format_field not in formatname_list:
                print "VCF doesnot contain ",format_field," format field"
                sys.exit(-1)
            elif len(line) > 8 and format_field in formatname_list:
                if line[0].startswith('chr'):
                    line[0] = line[0].replace('chr','')

                if format_field == 'GT':
                    GT = GT_parser(line[8], line[9:], len(sample_names))    
                    #print '\t'.join( ['chr'+ ':'.join(line[0:2]), var_type] + line[0:7] + [str(AF_final)] + info_values + GT )
                    print '\t'.join( line[0:7] + ['chr'+ ':'.join(line[0:2]), var_type] + [str(AF_final)] + info_values + GT )
                elif format_field == 'DP' or format_field == 'GQ' or format_field == 'PL' or format_field == 'AD' or format_field == 'MQ0':
                    processed_format = format_parser( line[8], line[9:], len(sample_names), format_field )
                    #print '\t'.join( ['chr'+ ':'.join(line[0:2]), var_type] + line[0:7] + [str(AF_final)]+ info_values + processed_format )
                    print '\t'.join( line[0:7] + ['chr'+ ':'.join(line[0:2]), var_type] + [str(AF_final)]+ info_values + processed_format )
                else:
                    print "Wrong format field name OR this field is not present in the VCF FORMAT Field"
                    sys.exit(-1)
    fin.close()


def info_parser(info_header, info_fields, info_col):
    """
    Function to parse INFO column
    This function will receive the list of INFO column names from the VCF header, i.e from 'info_fields dictionary' with key as INFO column name value as TYPE
    INFO column for a particular variant from the VCF
    """
    info = {}
    info_values = []
    for inf in sorted(info_fields.keys()):  #Initializing the dict values to None
        info[inf] = 'None'
    if info_col != '.':
        for col in info_col.split(';'):
            if info_fields[col.split('=')[0]] != 'Flag' :
                info[col.split('=')[0]] = col.split('=')[1]
            else: 
                info[col.split('=')[0]] = col.split('=')[0]
        for val in info_header:         #Getting the info column values in the same order as of info_header list in the list info_values 
            if val in info:
                info_values.append(info[val])
            else:
                print "Val Not Present"
    else:
        info_values =['.']*len(info_fields.keys())
    AF_list = []
    if ',' in info['AF']:
        AF_list = info['AF'].split(',')
        AF_final = min(AF_list)
    else:
        AF_final = info['AF']
    return info_values, AF_final

def GT_parser(format, Genotype_data,ln):
    """
    Gets the GT data for each variant and returns the GT as list
    """
    GT = []
    GT_pos = -9
    for k,x in enumerate(format.split(':')):
        if x == 'GT':
            GT_pos = k
    if GT_pos == -9:
        print "\n",'GT not available in the format field...',"\n"
        sys.exit(-1)
    for col in Genotype_data:
        GT.append(col.split(':')[GT_pos])
    return GT

def format_parser(format, Genotype_data,ln, format_field):
    """
    Parses the format_field specified by the user other than GT
    Had seen that in some variant format field may have just GT field and hence checking for length of FORMAT field, since GT would be the first genotype field and rest others comes after that
    """
    temp_format = []
    format_field_pos = -9
    if len(format.split(':')) > 1 :
        for k,x in enumerate(format.split(':')):
            if x == format_field:
                format_field_pos = k
        if format_field_pos == -9:
            print "\n",format_field, ' not available in the format field...', "\n"
            sys.exit(-1)
        for col in Genotype_data:
            if col != './.' and len(col.split(':')) >= format_field_pos + 1:
                temp_format.append(col.split(':')[format_field_pos])
            else:
                temp_format.append('.')
    else:
        temp_format = ['.']*ln

def vcf_add_ANNOVARannotations():

	"""
	FunctionName :   vcf_add_ANNOVARanotations
	Author      :   Tinu Thomas
	Date        :   05/12/2014
	Purpose     :   Given a VCF and VCF in text format with ANNOVAR annotations, adds the annotations on the VCF from the VCF in text format
	Usage       :   python add_ANNOVARannotations.py AnnovarOuput(*.hg19_mtuliannot)  Input.vcf
	"""

	try:
		vcftxt_annotations = sys.argv[1]
		vcf = sys.argv[2]
	except:
		print __doc__
		sys.exit(-1)

	variants = {}
	annot_type = {'Func.refGene':'String',
		'Gene.refGene':'String',
		'ExonicFunc.refGene':'String',
		'AAChange.refGene':'String',
		'Func.ensGene':'String',
		'Gene.ensGene':'String',
		'ExonicFunc.ensGene':'String',
		'AAChange.ensGene':'String',
		'snp137':'String',
		'snp137NonFlagged':'String',
		'cosmic68':'String',
		'popfreq_all':'Float',		
		'PopFreqMax':'Float',
		'1000G2012APR_ALL':'Float',
		'1000G2012APR_AFR':'Float',
		'1000G2012APR_AMR':'Float',
		'1000G2012APR_ASN':'Float',
		'1000G2012APR_EUR':'Float',
		'ESP6500si_AA':'Float',
		'ESP6500si_EA':'Float',
		'ESP6500si_ALL':'Float',
		'CG46':'Float',
		'clinvar_20140303':'String',
		'ljb23_all':'String',
		'cytoband':'String',
		'gwasCatalog':'String',
		'snp137':'String',
		'cosmic65':'String',
		'cg69':'String',
		'ljb2_all' : 'Float',
		'LJB2_SIFT':'Float',
		'LJB2_PolyPhen2_HDIV':'Float',
		'LJB2_PP2_HDIV_Pred':'String',
		'LJB2_PolyPhen2_HVAR':'Float',
		'LJB2_PolyPhen2_HVAR_Pred':'String',
		'LJB2_LRT':'Float',
		'LJB2_LRT_Pred':'String',
		'LJB2_MutationTaster':'Float',
		'LJB2_MutationTaster_Pred':'String',
		'LJB_MutationAssessor':'Float',
		'LJB_MutationAssessor_Pred':'String',
		'LJB2_FATHMM':'Float',
		'LJB2_GERP++':'Float',
		'LJB2_PhyloP':'Float',
		'LJB2_SiPhy':'Float',
		'LJB23_SIFT_score':'Float',
		'LJB23_SIFT_score_converted':'Float',
		'LJB23_SIFT_pred':'String',
		'LJB23_Polyphen2_HDIV_score':'Float',
		'LJB23_Polyphen2_HDIV_pred':'String',
		'LJB23_Polyphen2_HVAR_score':'Float',
		'LJB23_Polyphen2_HVAR_pred':'String',
		'LJB23_LRT_score':'Float',
		'LJB23_LRT_score_converted':'Float',
		'LJB23_LRT_pred':'String',
		'LJB23_MutationTaster_score':'Float',
		'LJB23_MutationTaster_score_converted':'Float',
		'LJB23_MutationTaster_pred':'String',
		'LJB23_MutationAssessor_score':'Float',
		'LJB23_MutationAssessor_score_converted':'Float',
		'LJB23_MutationAssessor_pred':'String',
		'LJB23_FATHMM_score':'Float',
		'LJB23_FATHMM_score_converted':'Float',
		'LJB23_FATHMM_pred':'String',
		'LJB23_RadialSVM_score':'Float',
		'LJB23_RadialSVM_score_converted':'Float',
		'LJB23_RadialSVM_pred':'String',
		'LJB23_LR_score':'Float',
		'LJB23_LR_pred':'String',
		'LJB23_GERP++':'Float',
		'LJB23_PhyloP':'Float',
		'LJB23_SiPhy':'Float',
		'cadd':'Float'
		}

	annot_vars = []
	chrom_pos = -9
	fin = return_file_handle(vcftxt_annotations)   #File with ANNOVAR annotations
	for line in fin:
		line = line.strip('\n\r').split('\t')
		if line[0] == 'Chr':	
			annot_index_start = line.index('Func.refGene')
			for k,x in enumerate(line[annot_index_start:]):
				if x != 'Otherinfo':
					annot = ( k + annot_index_start, x )
					annot_vars.append(annot)
				elif x == 'Otherinfo':
					chrom_pos = k + 5
		elif re.search(r'^chr(\d+|X|Y)|^(\d+|X|Y)', line[chrom_pos]):
			annot_list = []
			for vars in annot_vars:
				temp = line[vars[0]].replace(';',',')
				temp = temp.replace(' ','')
				if temp == '.':
					temp = '-999'
				annot_str = ';' + vars[1] + '=' + temp
				annot_list.append(annot_str)
			variants[ ':'.join(line[ chrom_pos:chrom_pos + 2 ] + [line[chrom_pos + 3 ]]) ] = ''.join(annot_list)
	n = 0
	p = 0
	fvcf = return_file_handle(vcf)
	for data in fvcf:
		data = data.strip('\n\r').split('\t')
		if re.search(r'^chr(\d+|X|Y)|^(\d+|X|Y)',data[0]):
			n += 1
			if ':'.join(data[0:2] +  [data[3]]) in variants :
				data[7] = data[7] + variants[':'.join(data[0:2] +  [data[3]])]	
				p += 1
				print '\t'.join(data)
			elif ':'.join(data[0:2] +  [data[3]]) not in variants:
				print '\t'.join(data)
		elif re.search(r'^#CHROM',data[0]):
			for y in annot_vars:
				print '##INFO=<ID=' + y[1] + ',Number=1,Type='+annot_type[y[1]]+',Description="ANNOVAR Annotation">'
			print '\t'.join(data)	
		else:
			print '\t'.join(data)
	fvcf.close()

def recalculateAF():
	"""
	Purpose :   Recalculates AF for ACAN run file and replaces old AF value with the new value. Else if AF is not found in the INFO field, calculates and adds it to the INFO field
	Usage   :   python /CommonDATA/Python_Scripts/recalculateAF.py  In.vcf  >Out.vcf
	"""

	try:
	    inVCF = sys.argv[1]
	except:
	    print __doc__
	    sys.exit(-1)
	info = {}
	fin = return_file_handle(inVCF)
	info_fields = {}
	for line in fin:
	    line = line.strip('\n\r')
	    if line.startswith('##INFO'):
        	print line
	        line = line.split(',')
	        infoname = re.search(r'^##INFO=<ID=(.+)',line[0])
        	infotype = re.search(r'^Type=(.*)',line[2])
	        info_fields[infoname.group(1)] = infotype.group(1)        
	    elif re.search(r'^chr(\d+|X|Y)|^(\d+|X|Y)',line):
        	if 'AC' in info_fields and 'AN' in info_fields:
	            line = line.split('\t')
        	    ALT = []
	            AC = []
        	    AF = []
	            if not re.search(r'AC=', line[7]):
        	        line[7] = line[7] + ';AC=0'
	            info_data = line[7].split(';')
        	    for k in info.keys():
                	info[k] = 'None'
	            for val in info_data:
        	        if val.split('=')[0] == 'AC' or val.split('=')[0] == 'AN':
                	    info[val.split('=')[0]] = val.split('=')[1]
	                    info[val.split('=')[0]] = val.split('=')[1]
        	            info[val.split('=')[0]] = val.split('=')[1]

	            #Checking for multiple alternate alleles in the vcf file and storing their names in the list named 'ALT'
        	    #Correspondingly storing multiple values of 'AC' in the list named AC
	            if ',' in line[4]:
        	        ALT = line[4].split(',')
                	AC = info['AC'].split(',')
	            else:
        	        ALT.append(line[4])
                	AC.append(info['AC'])
	            flag = 0
    	
        	    ln_ALT = len(ALT)
	            #Calculating the AF values and storing them in the list AF
        	    for alt_ct in range(0,ln_ALT):
                	if float(info['AN']) == 0:
	                    AF_value = 0
        	        else:
                	    AF_value = round(float(AC[alt_ct])/float(info['AN']) ,4)
	                AF.append(AF_value)
        	    AF_str = ','.join(str(e) for e in AF)
	            #if 'AF' in info_fields:
        	    line[7] = re.sub(r'AF=.*?;','AF='+AF_str+';',line[7])
	            #elif 'AF' not in info_fields and 'AC' in info:
        	    #    line[7] = line[7] + ';AF='+AF_str
	            print '\t'.join(line)
        	else:
	            print 'AC, AN or AF not present in the INFO field of the VCF'
        	    sys.exit(-1)
	    elif line.startswith('#CHROM'):
        	line = line.split('\t')
	        print '\t'.join(line)
	    else:
        	line = line.split('\t')
	        print '\t'.join(line)
	fin.close()


def RemoveACzerovariants():
	"""
	Purpose :   Given a VCF file would remove variants with AC = 0
	Usage   :   python /CommonDATA/Python_Scripts/RemoveACzerovariants.py Input.vcf >Output.vcf
	"""

	try:
	    inVCF = sys.argv[1]
	except:
	    print __doc__
	    sys.exit(-1)
	formatname_list = []
	fin = return_file_handle(inVCF)
	for line in fin:
	    line = line.strip('\n\r').split('\t')
	    if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
        	flag = 0
	        for val in line[7].split(';'):
        	    if '=' in val and val.split('=')[0] == 'AC' and ',' not in val.split('=')[1]  :
                	if float(val.split('=')[1]) != 0:
	                    flag = 1
        	            break
	            elif '=' in val and val.split('=')[0] == 'AC' and ',' in val.split('=')[1] :
        	        AC_values = list(set(val.split('=')[1].split(',')))
                	if AC_values != ['0']:
	                    flag = 1
        	if flag == 1 :
	            print '\t'.join(line)
	    else:
        	print '\t'.join(line)
	fin.close()



def CallRateFilter():
	"""
	Purpose : From the VCF calculates the callrate of variant per center and extracts only variants with call rate greater than user specified  genotyping call rate

	Usage   :   python GenotypingCallRateFilter.py  In.vcf  CallrateLimit  >Out.vcf
	"""
	try:
		inVCF = sys.argv[1]
		cr_limit = sys.argv[2]
	except:
		print __doc__
		sys.exit(-1)
	fin = return_file_handle(inVCF)
	for line in fin:
		line = line.strip().split('\t')
		msk_an = msk_cr = 0
		if re.search(r'^chr(\d+|X|Y)|^(\d+|X|Y)',line[0]):
			samples = len(line[9:])
			for val in line[7].split(';'):
				if val.split('=')[0] == 'AN':
					msk_an = float(val.split('=')[1])/2
					msk_cr = msk_an/float(samples)
			if (msk_cr >= float(cr_limit)):
				print '\t'.join(line)
			else:
				print '\t'.join(line)
	fin.close()

def CenterSpecific_ACANAF():
	"""
	This script would take AC, AN, AF fields from INFO column and adds center sepcific AC,AN, AF fields while 
	keeping original AC,AF, AN text intact. This would be useful when you want to merge VCFs with samples from different centers.
	From the AN field, calcualtes call rate for a variants and adds call rate to the info column	
	
	Usage 	: python CenterSpecific_ACANAF.py VCF CenterName >Output.vcf
	example : python CenterSpecific_ACANAF.py  MAYO_SimplexoSubset_GGA_VA_gatk2.6.4_Part2_07152013.vcf.gz Mayo >Out.vcf



	This script would take AC, AN, AF fields from INFO column and adds center sepcific AC,AN, AF fields while 
	keeping original AC,AF, AN text intact
	
	Usage 	: python CenterSpecific_ACANAF.py VCF CenterName
	example : python CenterSpecific_ACANAF.py  MAYO_SimplexoSubset_GGA_VA_gatk2.6.4_Part2_07152013.vcf.gz Mayo
	"""
	try:
		inVCF = sys.argv[1]
		CenterName = sys.argv[2]
	except:
		print CenterSpecific_ACANAF.__doc__
		sys.exit(-1)
	fin = return_file_handle(inVCF)
	for line in fin:
		line = line.strip().split('\t')
		if re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line[0]):
			flag_AF = 0
			flag_AC = 0
			flag_AN = 0
			common_AF_flag = 'N'
			gt_data = []
			samples = len(line[9:])
			for val in line[7].split(';'):
				if val.split('=')[0] == 'AF':
					cname_AF = CenterName + '_AF=' + val.split('=')[1] + ';'
					flag_AF = 1
					if ',' in val.split('=')[1] :
						if min(val.split('=')[1].split(',')) > 0.1 :
							common_AF_flag = 'Y'
					elif float(val.split('=')[1]) > 0.1 :
						common_AF_flag = 'Y'
					cname_common_flag = CenterName + '_common_AF=' + common_AF_flag + ';'
				elif val.split('=')[0] == 'AC':
					cname_AC = CenterName + '_AC=' + val.split('=')[1] + ';'
					flag_AC = 1
				elif val.split('=')[0] == 'AN':
					an = int(val.split('=')[1])/float(2)
					cr = an / float(samples)
					cname_AN = CenterName + '_AN=' + val.split('=')[1] + ';'
					flag_AN = 1
					cname_CR = CenterName + '_CR=' + str(round(float(cr),5)) + ';'
			if flag_AF == 1 and flag_AC == 1 and flag_AN == 1  :
				line[7] = cname_AC + cname_AF + cname_AN + cname_CR + cname_common_flag + line[7]
				print '\t'.join(line)
			else:
				print '\t'.join(line)


				for gt_dat in line[9:]:
					gt_data.append(gt_dat.split(':')[0])
				if all_same(gt_data):
					line[7] = line[7] + CenterName + '_AC=0;' + CenterName + '_AN=0;' + CenterName + '_AF=0'
					print '\t'.join(line)
				else:
					print "Missing AC, AN or AF values for the variant",'\t'.join(line[0:6])
					sys.exit(-1)					
		elif re.search(r'^#CHROM',line[0]):
			print '##INFO=<ID=' + CenterName + '_AC,Number=A,Type=Integer,Description="Center specific Allele count in genotypes, for each ALT allele, in the same order as listed">'
			print '##INFO=<ID=' + CenterName +'_AF,Number=A,Type=Float,Description="Center specific Allele Frequency, for each ALT allele, in the same order as listed">'
			print '##INFO=<ID=' + CenterName + '_AN,Number=1,Type=Integer,Description="Center specific Total number of alleles in called genotypes">'
			print '##INFO=<ID=' + CenterName + '_CR,Number=1,Type=Float,Description="Center specific call rate">'
			print '\t'.join(line)
		else:
			print '\t'.join(line)			
	fin.close()

def all_same(items) :
	return all(x == './.' for x in items)

def DP_Genotype_filter():
	"""
	Purpose :   Checks whether the DP in FORMAT field is lees than the specified DP_value by the user. IF yes, replaces the Genotype column with './.' for samples            less than the specified DP
	            Very Important to recalculate AC, AN and AF in INFO field after running this script
	Usage   :   python /CommonDATA/Python_Scripts/DP_Genotype_filter.py  In.vcf  DP_value  >Out.vcf 
	Example :   python /CommonDATA/Python_Scripts/DP_Genotype_filter.py  QUAL_recalibrated_snpEff.vcf  10  > DP_filtered.vcf
	"""

	try:
	    inVCF = sys.argv[1]
	    DP_value = sys.argv[2]
	except:
	    print __doc__
	    sys.exit(-1)
	formatname_list = []
	c= a =0
	fin = return_file_handle(inVCF)
	for line in fin:
	    line = line.strip('\n\r')
	    if line.startswith('##FORMAT'):
        	data = line
	        line = line.split(',')
	        formatname = re.search(r'^##FORMAT=<ID=(.+)',line[0])
        	formatname_list.append(formatname.group(1))
	        print '\t'.join(data.split('\t'))
	    elif re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line):
		c += 1
	        line = line.split('\t')
        	DP_pos = -9 
	        if len(line) <= 8:
        	    print "\n",'FORMAT field and Samples not found in the VCF ', line[0:2],"\n"
	            sys.exit(-1)
        	else:
	            for k,x in enumerate(line[8].split(':')):
        	        if x == 'DP':
                	    DP_pos = k
	        if DP_pos == -9:
		    print "\n",'DP not found in FORMAT field of the variant ',line[0:2],"\n"
	            sys.exit(-1)
        	else:
	            ln = len(sample_names) - 1
        	    for i in range(0, ln):
                	j = i + 9
	                if line[j] != '.' and line[j] != './.' and len(line[j].split(':')) >= DP_pos + 1:
        	            temp_DP = line[j].split(':')[DP_pos]
                	    if temp_DP != '.' and float(temp_DP) < float(DP_value):
                        	line[j] = './.'
	                else:
        	            line[j] = './.'
	            print '\t'.join(line)
	    elif line.startswith('#CHROM'):
        	if 'DP' not in formatname_list:
	            print "\n", 'DP not present in FORMAT field', "\n"
        	    sys.exit(-1)
	        else:
        	    line = line.split('\t')
	            sample_names = line[9:]
        	    print '\t'.join(line)
	    else:
        	line = line.split('\t')
	        print '\t'.join(line)
	fin.close()
    


def Filter_ChrPos_FromVCF(vcffile, PosvcfFile):
	"""
	Purpose	:	Given 2 VCF files, spits out variants in second VCF found in first.
	Usage   :   python Filter_ChrPos_FromVCF.py Variant.vcf Master.vcf >Output.vcf
	"""
	try:
		PosvcfFile = sys.argv[1]	
		vcffile = sys.argv[2]
	except:
		print Filter_ChrPos_FromVCF.__doc__
	        sys.exit(-1)
	fpos = return_file_handle(PosvcfFile)
	pos = {}
	pos_ct =0 
	for line in fpos:
		line = line.strip().split('\t')
		if re.search(r'^(\d+|X|Y|M)|^chr(\d+|X|Y|M)',line[0]):
			line[0] = line[0].replace("chr","")
			pos[":".join(line[0:2])]= 1			#Getting Gene names and storing them in the dictinary named 'pos'
	fpos.close()
	fin = return_file_handle(vcffile) #Opening vcf Family Distribution text file
	for data in fin:
		data = data.strip().split('\t')
		if data[0] == '#CHROM':
			print "\t".join(data)
		elif re.search(r'^(\d+|X|Y|M)|^chr(\d+|X|Y|M)',data[0]):	
			data[0] = data[0].replace('chr','')
			if ':'.join(data[0:2]) in pos:
				print "\t".join(data)
	fin.close()



def GQ_Genotype_filter():
	"""
	Purpose :   Checks whether the GQ in FORMAT field is lees than the specified GQ_value by the user. IF yes, replaces the Genotype column with './.' for samples            less than the specified GQ
            Very Important to recalculate AC, AN and AF in INFO field after running this script
	Usage   :   python /CommonDATA/Python_Scripts/GQ_Genotype_filter.py  In.vcf  GQ_value  >Out.vcf 
	Example :   python /CommonDATA/Python_Scripts/GQ_Genotype_filter.py  QUAL_recalibrated_snpEff.vcf  10  > GQ_filtered.vcf
	"""
	try:
		inVCF = sys.argv[1]
		GQ_value = sys.argv[2]
	except:
		print __doc__
		sys.exit(-1)
	formatname_list = []
	fin = return_file_handle(inVCF)
	for line in fin:
		line = line.strip('\n\r')
		if line.startswith('##FORMAT'):
			data = line
			line = line.split(',')
			formatname = re.search(r'^##FORMAT=<ID=(.+)',line[0])
			formatname_list.append(formatname.group(1))
			print '\t'.join(data.split('\t'))
		elif re.search(r'^(\d+|X|Y)|^chr(\d+|X|Y)', line):
			line = line.split('\t')
			GQ_pos = -9 
			if len(line) <= 8:
				print "\n",'FORMAT field and Samples not found in the VCF ', line[0:2],"\n"
				sys.exit(-1)
			else:
				for k,x in enumerate(line[8].split(':')):
					if x == 'GQ':
						GQ_pos = k
				if GQ_pos == -9:
				#print "\n",'GQ not found in FORMAT field of the variant ',line[0:2],"\n"
					print '\t'.join(line)
					#sys.exit(-1)
		"""
		else:
			ln = len(sample_names)
			for i in range(0, ln):
				j = i + 9
				if line[j] != '.' and line[j] != './.' and len(line[j].split(':')) >= GQ_pos + 1:
					temp_GQ = line[j].split(':')[GQ_pos]
			if temp_GQ != '.' and float(temp_GQ) < float(GQ_value):
				line[j] = './.'
			else:
				line[j] = './.'
			print '\t'.join(line)
		elif line.startswith('#CHROM'):
			if 'GQ' not in formatname_list:
				print "\n", 'GQ not present in FORMAT field', "\n"
				sys.exit(-1)
			else:
				line = line.split('\t')
				sample_names = line[9:]
				print '\t'.join(line)
			else:
				line = line.split('\t')
				print '\t'.join(line)
		"""
	fin.close()

def autosomal_dominance(pedfile, vcf):
	"""
	Date	:	07/11/2014
	Function :	autosomal_domimance
	Purpose	:	Given a VCF and PED file, gives the autosomal dominant alles in each family and also the autosomal dominant alleles per gene
	Usage	:	python scriptname.py PEDfile Input.vcf FamilyID 'Genotype to be considered for unknown affected status samples'
	"""	
	try:
		pedfile = sys.argv[1]
		inVCF = sys.argv[2]
		familyID = sys.argv[3]	
		unknown_geno = sys.argv[4]
	except:
		print autosomal_dominance.__doc__
		sys.exit(-1)

	ped_info = defaultdict(list)
	affected_status = {}
	affected_count = {}
	unaffected_status = {}
	unaffected_count = {}
	unknown_status = {}
	unknown_count = {}
	family = {}
	sample_ct = {}
	samples = defaultdict(list)
	affected_samples = []
	unaffected_samples = []
	unknown_samples = []
	fped = return_file_handle(pedfile)
	for ped in fped:
		if not ped.startswith('#'):
			ped = ped.strip().split('\t')
			ped_info[':'.join(ped[0:2])] = ped[2:]
			family[ped[1]] = ped[0]
			samples[ped[0]].append(ped[1])
			if ped[0] not in sample_ct :
				sample_ct[ped[0]] = 1
			else:
				sample_ct[ped[0]] = sample_ct[ped[0]] + 1
			if ped[0] not in affected_count and ped[5] == '2':	#Storing affected sample information
				affected_samples.append(ped[1])
				affected_count[ped[0]] = 1
			elif ped[0] in affected_count and ped[5] == '2':
				affected_samples.append(ped[1])
				affected_count[ped[0]] += 1
			if ped[0] not in unaffected_count and ped[5] == '1':
				unaffected_samples.append(ped[1])
				unaffected_count[ped[0]] = 1
			elif ped[0] in unaffected_count and ped[5] == '1':
				unaffected_samples.append(ped[1])
				unaffected_count[ped[0]] += 1
			if ped[0] not in unknown_count and ped[5] == '0':
				unknown_samples.append(ped[1])
				unknown_count[ped[0]] = 1
			elif ped[0] in unknown_count and ped[5] == '0':
				unknown_samples.append(ped[1])
				unknown_count[ped[0]] += 1

		#affected_status[] = ped[5]i
	#print 'family',family
	#print 'samples', samples
	#print 'ped_info', ped_info
	#print 'affected samples', affected_samples
	#print 'affected_count', affected_count
	#print 'unaffected_count', unaffected_count
	for sam in ped_info.items():
		if sam[1][3] == '2' :
			print sam[0].split(':')[1]
	sample_pos = {}
	fvcf = return_file_handle(inVCF)
	for line in fvcf:
		line = line.strip().split('\t')
		if re.search(r'^#CHROM',line[0]):
			for k,x in enumerate(line[9:]):
				sample_pos[x] = k
			samples = line[9:]	
			print '\t'.join(line)
		elif re.search(r'^chr(\d+|X|Y)|^(\d+|X|Y)',line[0]):
			sam_affect_ct = sam_unaffect_ct = sam_unknown_ct = 0
			for gt_pos, gt in enumerate(line[9:]):
				if  samples[gt_pos] in affected_samples and family[samples[gt_pos]] == familyID and gt.split(':')[0] == '0/1':
					sam_affect_ct += 1
				elif samples[gt_pos] in unaffected_samples and family[samples[gt_pos]] == familyID and gt.split(':')[0] == '0/0':
					sam_unaffect_ct += 1
				elif samples[gt_pos] in unknown_samples and family[samples[gt_pos]] == familyID and gt.split(':')[0] == unknown_geno:
					sam_unknown_ct += 1
			if sam_affect_ct == affected_count[familyID] and sam_unaffect_ct == unaffected_count[familyID] and sam_unknown_ct == unknown_count[familyID]:
				print '\t'.join(line)
		else:
			print '\t'.join(line)
	fvcf.close()
