## IMPORT
import xlsxwriter

## FUNCTIONS
def get_content(sample):
	"""Retrieves the content of hla.all"""
	content=[]
	with open(sample+'.hla.all','r') as hla_all:
		for line in hla_all.readlines():	
			l=line.rstrip().split('\t')
			content.append(l)
	return content

def organise_content(content):
	"""Creates a dictionnary where each HLA gene is the key linked to a list containing the different hla typing results"""
	cursor=content[0][0].split('*')
	cursor=cursor[0]
	HLA={'HLA-A':[]}
	for i in content:
		i[-3]=int(i[-3])
		i[-2]=int(i[-2])
		i[-1]=int(i[-1])
		hla=i[0].split('*')
		if cursor==hla[0]: 					# If the new result is still for the same gene
			HLA[hla[0]]=HLA[hla[0]]+[i] 	# Add the result to the list
		else:
			cursor=hla[0]		# If not, we create a new key 
			HLA[hla[0]]=[i] 	# And we initiate it as a new list
	return HLA


def check_rank(dic):
	"""Looks for a result with more exons typed than the first one"""
	HLA=[]
	for hla in dic:
		max_ex=dic[hla][0][-1]
		result=dic[hla][0]
		HLA=HLA+[[result]]
		found=0
		for r in dic[hla]:
			if found==0:
				if r[-1]>max_ex:
					possibility=r
					HLA[-1]=HLA[-1]+[possibility]
					found=1
	return HLA	# Return the list with either one element per gene (not correction to consider) or two elements (one result was typed using more exon than the first ranked result)
	
def relevant_correction(HLA,sample,rel):
	"""Summarizes per sample which gene has that correction possibility and stores the 2 options"""
	rel[sample]={}
	for i in HLA:
		if len(i)>1: # If there is a possibility of correction
			elem=i[0][0].split('*')
			elem=elem[0][4:]
			rel[sample][elem]=i # Create a dictionnary with the gene as the key and containing the two results
	return rel
	

def get_ID():
	"""Retrieves all samples ID"""
    ID=[]
    with open('list_ID.txt','r') as list_ID:
        for i in list_ID.readlines():
            ID=ID+[i.rstrip('\n')]
    return ID

def output_csv_mean(rel):
	"""Creates the output"""
	name="bwa_correction.xlsx"
	workbook = xlsxwriter.Workbook(name)
	worksheet = workbook.add_worksheet()
	bold = workbook.add_format({'bold': True})
	# Header
	row = 0
	col = 0
	worksheet.write(row, col,"Cases when there is another possibility of result where more exons were used to be typed")
	#1st column
	row=2
	for key in rel:
		col=0
		worksheet.write(row,col,key,bold)
		if rel[key]=={}:
			row+=1
		else:
			for gene in rel[key]:
				col=1
				worksheet.write(row,col,gene)
				for i in rel[key][gene]:
					col=2
					for j in i:
						worksheet.write(row,col,j)
						col+=1
					row+=1
	workbook.close()

## MAIN
ID=get_ID()
rel={}
for sample in ID:
	content=get_content(sample)
	HLA=organise_content(content)
	list_hla=check_rank(HLA)
	rel=relevant_correction(list_hla,sample,rel)
	output_csv_mean(rel)
	

