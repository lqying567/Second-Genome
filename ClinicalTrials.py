##Author: Lanqing Ying
import sys,os,argparse
import pandas as pd
from ClinTrials.ClinTrialsUtil import searchClinTrialDB
import requests
from bs4 import BeautifulSoup

p=argparse.ArgumentParser(description="Gereate actionfile for Onco result")
p.add_argument("-s", "--summary", type=str, help="summray output file")
p.add_argument("-m", "--mutation", type=str, help="somatic mutation txt file")
p.add_argument("-v", "--vcf", type=str, help="somatic mutation vcf file")
#p.add_argument("-o","--outFolder", type=str, help="output_file_folder", default="")
if len(sys.argv)<3:
	print sys.stderr, p.print_help()
	exit(0)
args=p.parse_args()

summary_file = open(args.summary,'r').readlines()         ## *_output.txt
somatic_mutation=open(args.mutation,'r').readlines()

cancer_type=['Breast Cancer']
cancer_type=[i.lower() for i in cancer_type]

knowledgeBase=open(path+'KnowledgeBaseGenes.txt','r').readlines()
knowledgeBaseGenes=[l.strip() for l in knowledgeBase]

suppressorGenes=open(path+'Tumor_suppressor_genes.txt','r').readlines()
suppressorGenes=[l.strip() for l in suppressorGenes]

FDA_info=pd.read_csv(path+'Drugs_FDA_info.txt',sep='\t')

therapy_db=pd.read_excel(knowledgeBaseFile, 'Targeted_Therapy')

drug_db=pd.read_excel(knowledgeBaseFile, 'DrugInfor')

sample=args.summary.split('/')[-1].split('_output')[0]
outputfile=args.summary.replace(sample+'_output.txt','actionfile_'+sample+'.txt')
#output = open(("actionfile_"+sample+'.txt'),'w')
output = open(outputfile,'w')

Temp_output=open(args.summary.replace("_output.txt", "_temp.txt"),'w')
output.write('Group\tGene\tResultColumn\tResultValue\tResultValueDetail\n')

state_to_code={'federated states of micronesia': 'fm', 'wyoming': 'wy', 'colorado': 'co', 'guam': 'gu', 'nebraska': 'ne', 'washington': 'wa', 'american samoa': 'as', 'alaska': 'ak', 'armed forces europe': 'ae', 'wisconsin': 'wi', 'nevada': 'nv', 'maine': 'me', 'north dakota': 'nd', 'mississippi': 'ms', 'south dakota': 'sd', 'palau': 'pw', 'new jersey': 'nj', 'new hampshire': 'nh', 'oklahoma': 'ok', 'delaware': 'de', 'minnesota': 'mn', 'north carolina': 'nc', 'illinois': 'il', 'new york': 'ny', 'arkansas': 'ar', 'west virginia': 'wv', 'puerto rico': 'pr', 'indiana': 'in', 'maryland': 'md', 'louisiana': 'la', 'idaho': 'id', 'armed forces americas': 'aa', 'iowa': 'ia', 'virgin islands': 'vi', 'arizona': 'az', 'michigan': 'mi', 'kansas': 'ks', 'utah': 'ut', 'virginia': 'va', 'oregon': 'or', 'armed forces africa': 'ae', 'connecticut': 'ct', 'montana': 'mt', 'california': 'ca', 'massachusetts': 'ma', 'armed forces canada': 'ae', 'rhode island': 'ri', 'armed forces middle east': 'ae', 'vermont': 'vt', 'georgia': 'ga', 'northern mariana islands': 'mp', 'tennessee': 'tn', 'florida': 'fl', 'hawaii': 'hi', 'kentucky': 'ky', 'pennsylvania': 'pa', 'district of columbia': 'dc', 'marshall islands': 'mh', 'texas': 'tx', 'armed forces pacific': 'ap', 'missouri': 'mo', 'south carolina': 'sc', 'ohio': 'oh', 'alabama': 'al', 'new mexico': 'nm'}
AA_code={'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}

def getTrialInfo(ID):
	zip_file=temp_folder+ID+'.zip'

	os.system('wget "https://clinicaltrials.gov/ct2/results?cond=&term=%s&cntry=&state=&city=&dist=&studyxml=true" -O %s -o /dev/null'%(ID,zip_file))

	os.system('unzip -o %s -d %s 2>/dev/null 1>/dev/null'%(zip_file,temp_folder))
	os.system("rm %s"%zip_file)

	file=temp_folder+ID+'.xml'
	lines=open(file,'r').readlines()
	t=[l for l in lines if '<brief_title>' in l][0]
	title=t.split('<brief_title>')[1].split('</brief_title>')[0]

	p=[l for l in lines if '<phase>' in l][0]
	phase=p.split('<phase>')[1].split('</phase>')[0]

	itv=[l.split('<intervention_name>')[1].split('</intervention_name>')[0] for l in lines if '<intervention_name>' in l]
	if len(itv)>5:
		itv=itv[:5]+['more']
	therapy=' + '.join(itv)

	state=[l.split('<state>')[1].split('</state>')[0].lower() for l in lines if '<state>' in l]
	state_uq=set(state)
	state_short=[state_to_code[s].upper() for s in state_uq if s in state_to_code.keys()]

	if len(state_uq)<=3:
		location=', '.join([x.capitalize() for x in state_uq])
	elif len(state_uq)<=10:
		location=', '.join(state_short)
	else:
		location=', '.join(state_short[0:5])+' + '+str(len(state_uq)-5)+' more'

	return [ID, therapy, title, phase, location]

pred_head=['dbNSFP_LRT_pred','dbNSFP_Polyphen2_HDIV_pred',
'dbNSFP_Polyphen2_HVAR_pred','dbNSFP_SIFT_pred',
'dbNSFP_MutationTaster_pred','dbNSFP_MutationAssessor_pred']

def getPred(chr,pos,file):
	vcf=open(file,'r').readlines()
	for line in vcf:
		if line.startswith('#'): continue
		v=line.split('\t')
		if v[0]==chr and int(v[1])==int(pos):
			count=0
			total=0
			info=v[7]
			for i in range(len(pred_head)):
				h=pred_head[i]+'='
				if h in info:
					total+=1
					pred=info.split(h)[1].split(';')[0].split(',')[0]
					if i<3 and pred=='D':
						count+=1
					if i>=3 and (pred!='L' and pred!='N'):
						count+=1
			if total==0:
				return 'Tolerated'
			if count/total*1.0>0.5:
				return 'Deleterious'
			else:
				return 'Tolerated'

def outputTherapy(group,num):
	if group==[]:
		return '%s\t\tMolecular_Abnormality\tNo medically actionable mutations were detected in this category.\t\n'%num
	else:
		line=''
		for item in group:
			line=line+'%s\t%s\tMolecular_Abnormality\t%s\t\n'%(num,item[0],item[1])
			line=line+'%s\t%s\tTherapies\t%s\t\n'%(num,item[0],item[2])
			line=line+'%s\t%s\tApproved_For\t%s\t\n'%(num,item[0],item[3])
			line=line+'%s\t%s\tLOE\t%s\t\n'%(num,item[0],item[4])
			line=line+'%s\t%s\tAllele_Freq\t%s\t\n'%(num,item[0],item[5])

		return line

def getTherapy(gene,mut,cancer,freq):
	mutation=mut.upper()
	mut_db=therapy_db[(therapy_db['Gene']==gene) & (therapy_db['Mutation']==mutation)]
	group1=[]
	group2=[]
	group3=[]

	if mut_db.empty:
		return [], [], []
	
	mut_db['Cancer_lower']=[c.lower() for c in mut_db['Cancer Type']]
	cancer_db=mut_db[mut_db['Cancer_lower'].isin(cancer)]
	if not cancer_db.empty:
		for i in range(len(cancer_db.index)):
			Therapy=cancer_db['Potential_benefit'].values[i]
			sensLOE=cancer_db['Sens_LOE'].values[i]
			#cancerType=cancer_db['Cancer Type'].values[i]
			cancerType='/'
			if Therapy in FDA_info['Drug'].values:
				cancerType=FDA_info[FDA_info['Drug']==Therapy]['Treatment'].values[0]

			if isinstance(sensLOE,float):
				sensLOE='/'

			Evidence=cancer_db['Sens_Evidence'].values[i]
			effect=cancer_db['Effect_Summary'].values[i]
			if isinstance(Therapy,float): continue
			row=[gene,mutation,Therapy,cancerType,sensLOE,freq,Evidence,effect]

			if sensLOE=='A' or sensLOE=='B':
				group1.append(row)

			if sensLOE=='C' or sensLOE=='D':
				group2.append(row)


			resTherapy=cancer_db['Lack_of_benefit'].values[i]
			resLOE=cancer_db['Res_LOE'].values[i]
			if resLOE!='A' and resLOE!='B': continue

			resEvidence=cancer_db['Res_Evidence'].values[i]
			if isinstance(resTherapy,float): continue
			row=[gene,mutation,resTherapy,cancerType,resLOE,freq,resEvidence]
			group3.append(row)

	non_cancer_db=mut_db[~mut_db['Cancer_lower'].isin(cancer)]
	if not non_cancer_db.empty:
		for i in range(len(non_cancer_db.index)):
			sensLOE=non_cancer_db['Sens_LOE'].values[i]
			Therapy=non_cancer_db['Potential_benefit'].values[i]
			if sensLOE!='A': continue
			#cancerType=non_cancer_db['Cancer Type'].values[i]
			cancerType='/'
			if isinstance(Therapy,float): continue
			if Therapy in FDA_info['Drug'].values:
				cancerType=FDA_info[FDA_info['Drug']==Therapy]['Treatment'].values[0]


			Evidence=non_cancer_db['Sens_Evidence'].values[i]
			effect=non_cancer_db['Effect_Summary'].values[i]

			row=[gene,mutation,Therapy,cancerType,'C',freq,Evidence,effect]
			group2.append(row)	

	def takeLOE(item):
		return item[4].replace('/','[')

	group1.sort(key=takeLOE)
	group2.sort(key=takeLOE)
	group3.sort(key=takeLOE)

	return group1, group2, group3


def sigVariants(group):
	for item in group:
		gm=item[0]+' '+item[1]
		if gm in geneMutation: continue
		geneMutation.append(gm)


def getDrug(group):
	line=''
	if group==[]:
		return line

	for item in group:
		therapy=item[2]
		if therapy in drug:
			continue
		else:
			drug.append(therapy)
			if therapy in drug_db['Drug.GenericName'].values:
				link=drug_db[drug_db['Drug.GenericName']==therapy]['FDA_Label_Link'].values[0]
				line=line+'7\t%s\tLink\t%s\t\n'%(therapy,link)

	return line


############ Generate Reference 
def generateRef(pubmedID):
	url='https://pubmed.ncbi.nlm.nih.gov/'+str(pubmedID)
	r = requests.get(url)

	if not r.history: 
		return ('')
		
	soup = BeautifulSoup(r.content, 'html.parser')

	title=soup.title.text.split('- PubMed')[0].strip()

	if title.startswith('No items found') or title.startswith('Page not found'):
		return ('')

	#jour=soup.find('div',class_='cit').text.split('doi:')[0].strip()
	source=soup.find('div',class_='article-source')
	jour1=source.find('button',class_="journal-actions-trigger trigger").text.strip()
	jour2=source.find('span',class_='cit').text.strip()
	jour=jour1+'. '+jour2

	#authors=soup.find('div',class_='auths').find_all('a')
	authors=soup.find('div',class_='authors-list').find_all('a',class_='full-name')

	auList=[au.text.strip() for au in authors]
	if len(auList)>5:
		auList=auList[:3]+['et al']

	auLine=', '.join(auList)

	ref=auLine+'. '+title+'. '+jour
	return (ref)

def getReferece(group):
	if group!=[]:
		for item in group:
			Evidence=item[6]
			if isinstance(Evidence,float):continue
			if str(Evidence).startswith('NCCN Guideline') or str(Evidence).startswith('FDA'):continue
			
			ref=generateRef(str(Evidence))
			if ref in citations or ref=='': continue
			citations.append(ref)


def getEffect(group):
	if group!=[]:
		for item in group:
			gene=item[0]
			if gene in geneList:continue
			effects.append([gene,item[7]])
			geneList.append(gene)




