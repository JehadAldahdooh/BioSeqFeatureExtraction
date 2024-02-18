

from __future__ import print_function 
import sys 
from Bio.SeqUtils import ProtParamData  # Local 
from Bio.SeqUtils import IsoelectricPoint  # Local 
from Bio.Seq import Seq 
from Bio import Entrez
from Bio.SeqUtils import GC, molecular_weight
#from Bio.Alphabet import IUPAC 
from Bio.Data import IUPACData 
from Bio.SeqUtils import molecular_weight 
import re
import math
import string




AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
_Hydrophobicity={"A":0.02,"R":-0.42,"N":-0.77,"D":-1.04,"C":0.77,"Q":-1.10,"E":-1.14,"G":-0.80,"H":0.26,"I":1.81,"L":1.14,"K":-0.41,"M":1.00,"F":1.35,"P":-0.09,"S":-0.97,"T":-0.77,"W":1.71,"Y":1.11,"V":1.13}
_AvFlexibility={"A":0.357,"R":0.529,"N":0.463,"D":0.511,"C":0.346,"Q":0.493,"E":0.497,"G":0.544,"H":0.323,"I":0.462,"L":0.365,"K":0.466,"M":0.295,"F":0.314,"P":0.509,"S":0.507,"T":0.444,"W":0.305,"Y":0.420,"V":0.386}
_Polarizability={"A":0.046,"R":0.291,"N":0.134,"D":0.105,"C":0.128,"Q":0.180,"E":0.151,"G":0.000,"H":0.230,"I":0.186,"L":0.186,"K":0.219,"M":0.221,"F":0.290,"P":0.131,"S":0.062,"T":0.108,"W":0.409,"Y":0.298,"V":0.140}
_FreeEnergy={"A":-0.368,"R":-1.03,"N":0.0,"D":2.06,"C":4.53,"Q":0.731,"E":1.77,"G":-0.525,"H":0.0,"I":0.791,"L":1.07,"K":0.0,"M":0.656,"F":1.06,"P":-2.24,"S":-0.524,"T":0.0,"W":1.60,"Y":4.91,"V":0.401}
_Steric={"A":0.52,"R":0.68,"N":0.76,"D":0.76,"C":0.62,"Q":0.68,"E":0.68,"G":0.00,"H":0.70,"I":1.02,"L":0.98,"K":0.68,"M":0.78,"F":0.70,"P":0.36,"S":0.53,"T":0.50,"W":0.70,"Y":0.70,"V":0.76}
_AAProperty=(_Hydrophobicity,_AvFlexibility,_Polarizability,_FreeEnergy,_Steric)
_AAPropertyName=('_Hydrophobicity','_AvFlexibility','_Polarizability','_FreeEnergy','_Steric')			

class ProteinAnalysis(object): 
    def __init__(self, prot_sequence, monoisotopic=False): 
        self.sequence = prot_sequence.upper()
        self.amino_acids_content = None
        self.amino_acids_percent = None
        self.length = len(self.sequence)
        self.monoisotopic = monoisotopic
 
    def count_amino_acids(self): 
        if self.amino_acids_content is None: 
            prot_dic = dict((k, 0) for k in IUPACData.protein_letters) 
            for aa in prot_dic: 
                prot_dic[aa] = self.sequence.count(aa)  
            self.amino_acids_content = prot_dic 
 
        return self.amino_acids_content 
 
    def get_amino_acids_percent(self): 
         if self.amino_acids_percent is None: 
             aa_counts = self.count_amino_acids() 
             percentages = {} 
             for aa in aa_counts: 
                 percentages[aa] = aa_counts[aa] / float(self.length) 
 
             self.amino_acids_percent = percentages 
 
         return self.amino_acids_percent 
 
    def molecular_weight(self): 
        return molecular_weight(self.sequence,seq_type="protein", monoisotopic=self.monoisotopic) 

    def aromaticity(self): 
        aromatic_aas = 'YWF' 
        aa_percentages = self.get_amino_acids_percent() 
        aromaticity = sum(aa_percentages[aa] for aa in aromatic_aas) 
        return aromaticity 

    def instability_index(self): 

        index = ProtParamData.DIWV 
        score = 0.0 
        for i in range(self.length - 1): 
            this, next = self.sequence[i:i + 2] 
            dipeptide_value = index[this][next] 
            score += dipeptide_value 

        return (10.0 / self.length) * score 
    
    def gravy(self): 
        total_gravy = sum(ProtParamData.kd[aa] for aa in self.sequence) 
        return total_gravy / self.length 
    
    def isoelectric_point(self): 
        aa_content = self.count_amino_acids() 
        ie_point = IsoelectricPoint.IsoelectricPoint(self.sequence, aa_content) 
        return ie_point.pi() 

##################################################################################################
def _mean(listvalue):
	"""
	The mean value of the list data.
	"""
	return sum(listvalue)/len(listvalue)
##################################################################################################
def _std(listvalue,ddof=1):
	"""
	The standard deviation of the list data.
	"""
	mean=_mean(listvalue)
	temp=[math.pow(i-mean,2) for i in listvalue]
	res=math.sqrt(sum(temp)/(len(listvalue)-ddof))
	return res
##################################################################################################
def NormalizeEachAAP(AAP):
	"""
	####################################################################################
	All of the amino acid indices are centralized and 
	
	standardized before the calculation.
	
	Usage:
	
	result=NormalizeEachAAP(AAP)
	
	Input: AAP is a dict form containing the properties of 20 amino acids.
	
	Output: result is the a dict form containing the normalized properties 
	
	of 20 amino acids.
	####################################################################################
	"""
	if len(AAP.values())!=20:
		print ('You can not input the correct number of properities of Amino acids!')
	else:
		Result={}
		for i,j in AAP.items():
			Result[i]=(j-_mean(AAP.values()))/_std(AAP.values(),ddof=0)

	return Result

####################################################################################
def CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,AAP,AAPName):
		
	"""
	####################################################################################
	you can use the function to compute MoreauBrotoAuto
	
	descriptors for different properties based on AADs.
	
	Usage:
	
	result=CalculateEachNormalizedMoreauBrotoAuto(protein,AAP,AAPName)
	
	Input: protein is a pure protein sequence.
	
	AAP is a dict form containing the properties of 20 amino acids (e.g., _AvFlexibility).
	
	AAPName is a string used for indicating the property (e.g., '_AvFlexibility'). 
	
	Output: result is a dict form containing 30 Normalized Moreau-Broto autocorrelation
	
	descriptors based on the given property.
	####################################################################################
	"""
		
	AAPdic=NormalizeEachAAP(AAP)

	Result={}
	for i in range(1,31):
		temp=0
		for j in range(len(ProteinSequence)-i):
			temp=temp+AAPdic[ProteinSequence[j]]*AAPdic[ProteinSequence[j+1]]
		if len(ProteinSequence)-i==0:
			Result[AAPName+str(i)]=round(temp/(len(ProteinSequence)),3)
		else:
			Result[AAPName+str(i)]=round(temp/(len(ProteinSequence)-i),3)

	return Result



def CalculateAAComposition(ProteinSequence):

	"""
	########################################################################
	Calculate the composition of Amino acids 
	
	for a given protein sequence.
	
	Usage:
	
	result=CalculateAAComposition(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing the composition of 
	
	20 amino acids.
	########################################################################
	"""
	LengthSequence=len(ProteinSequence)
	Result={}
	for i in AALetter:
		Result[i]=round(float(ProteinSequence.count(i))/LengthSequence*100,3)
	return Result
def CalculateDipeptideComposition(ProteinSequence):
	"""
	########################################################################
	Calculate the composition of dipeptidefor a given protein sequence.
	
	Usage: 
	
	result=CalculateDipeptideComposition(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing the composition of 
	
	400 dipeptides.
	########################################################################
	"""

	LengthSequence=len(ProteinSequence)
	Result={}
	for i in AALetter:
		for j in AALetter:
			Dipeptide=i+j
			Result[Dipeptide]=round(float(ProteinSequence.count(Dipeptide))/(LengthSequence-1)*100,2)
	return Result
def Getkmers():
	"""
	########################################################################
	Get the amino acid list of 3-mers. 
	
	Usage: 
	
	result=Getkmers()
	
	Output: result is a list form containing 8000 tri-peptides.
	
	########################################################################
	"""
	kmers=list()
	for i in AALetter:
		for j in AALetter:
			for k in AALetter:
				kmers.append(i+j+k)
	return kmers


def GetSpectrumDict(proteinsequence):
	"""
	########################################################################
	Calcualte the spectrum descriptors of 3-mers for a given protein.
	
	Usage: 
	
	result=GetSpectrumDict(protein)
	
	Input: protein is a pure protein sequence.
	
	Output: result is a dict form containing the composition values of 8000
	
	3-mers.
	########################################################################
	"""
	result={}
	kmers=Getkmers()
	for i in kmers:
		result[i]=len(re.findall(i,proteinsequence))
	return result
	
    
####################################################################################################################################
######                                 Main program starts from here                  #######################
####################################################################################################################################

import pandas as pd
import requests
from Bio import Entrez
Entrez.email = 'zia.rehman@helsinki.fi'

Uniprot_API='https://rest.uniprot.org/uniprotkb/'
# file=open('Proteins_features.csv',"r")
# proteins=file.readlines()   # one column file
# file.close()

proteins = ['P31152'] # list proteins here #

file=open('Protein_features_temp.csv','w')  

header="Uniprot_id;Pro_Sequence_len;"
for i in range (0,len(proteins)): 
    Uniprot = "temp"
    Uniprot=proteins[i].replace("\n","")
    request = requests.get(Uniprot_API+Uniprot)  # request  uniprot DB to extract protein sequence
    uniprot_API_json=request.json()

    if 'sequence' in uniprot_API_json:
        protein_seq=uniprot_API_json['sequence']['value']
    else: 
        handle = Entrez.efetch(db="protein", id=Uniprot, retmode="xml")  # request  NCBI DB to extract protein sequence if not found in UniProt
        records = Entrez.read(handle)
        seq_len=records[0]['GBSeq_length']
        prot_name=records[0]['GBSeq_definition'][0:].capitalize()
        protein_seq=str(records[0]["GBSeq_sequence"][0:]).upper()
    try:
        X = ProteinAnalysis(protein_seq) 
        seq_length=len(protein_seq)
    except:
        print("Something went wrong")
        continue

    Feature_AvFlexibility=CalculateEachNormalizedMoreauBrotoAuto(protein_seq,_AvFlexibility,'_AvFlexibility')
    try:
        write_line=Uniprot+";" + str(seq_length) + ";" + str(X.molecular_weight())+";" +str(X.aromaticity())+";" + str(X.gravy())+";" +str(X.instability_index())+";" +str(X.isoelectric_point())
        header=header+"Molecular_weight;Aromaticity;Gravy;Instability;Isoelectric_point"
    except:  continue  
    Feature_AvFlexibility=CalculateEachNormalizedMoreauBrotoAuto(protein_seq,_AvFlexibility,'_AvFlexibility')
    list_keys = [ k for k in Feature_AvFlexibility ]
    list_values = [ v for v in Feature_AvFlexibility. values()]    
    
    # Manually writing line by line to the output file
    if i==0:
        for k in range (len(list_keys)):
            header=header+';' + str(list_keys[k])   
           
    for k in range (len(list_values)):
        write_line=write_line+';' + str(list_values[k])
    Feature_Hydrophobicity=CalculateEachNormalizedMoreauBrotoAuto(protein_seq,_Hydrophobicity,'_Hydrophobicity')
    list_keys = [ k for k in Feature_Hydrophobicity ]
    list_values = [ v for v in Feature_Hydrophobicity. values()]
    if i==0:
        for k in range (len(list_keys)):
            header=header+';' + str(list_keys[k])   
    for k in range (len(list_values)):
        write_line=write_line+';' + str(list_values[k])
    Featur_Polarizabilitye=CalculateEachNormalizedMoreauBrotoAuto(protein_seq,_Polarizability,'_Polarizability')
    list_keys = [ k for k in Featur_Polarizabilitye ]
    list_values = [ v for v in Featur_Polarizabilitye. values()]
    if i==0:
        for k in range (len(list_keys)):
            header=header+';' + str(list_keys[k])   
    for k in range (len(list_values)):
        write_line=write_line+';' + str(list_values[k])
    Feature_FreeEnergy=CalculateEachNormalizedMoreauBrotoAuto(protein_seq,_FreeEnergy,'_FreeEnergy')
    list_keys = [ k for k in Feature_FreeEnergy ]
    list_values = [ v for v in Feature_FreeEnergy. values()]
    if i==0:
        for k in range (len(list_keys)):
            header=header+';' + str(list_keys[k])   
    for k in range (len(list_values)):
        write_line=write_line+';' + str(list_values[k])
    Feature_Steric=CalculateEachNormalizedMoreauBrotoAuto(protein_seq,_Steric,'_Steric')
    list_keys = [ k for k in Feature_Steric ]
    list_values = [ v for v in Feature_Steric. values()]
    if i==0:
        for k in range (len(list_keys)):
            header=header+';' + str(list_keys[k])   
    for k in range (len(list_values)):
        write_line=write_line+';' + str(list_values[k])
        
    Monopeptide_features=CalculateAAComposition(protein_seq)
    list_keys = [ k for k in Monopeptide_features ]
    list_values = [ v for v in Monopeptide_features. values() ]
    if i==0:
        for k in range (len(list_keys)):
            header=header+';' + str(list_keys[k])   
    for k in range (len(list_values)):
        write_line=write_line+';' + str(list_values[k])
    
    Dipeptide_features=CalculateDipeptideComposition(protein_seq)
    list_keys = [ k for k in Dipeptide_features ]
    list_values = [ v for v in Dipeptide_features. values() ]
    if i==0:
        for k in range (len(list_keys)):
            header=header+';' + str(list_keys[k])   
    for k in range (len(list_values)):
        write_line=write_line+';' + str(list_values[k])
        
       
    if i==0:
        file.write("%s \n" % header)
        file.write("%s \n" % write_line)
    else:
        file.write("%s \n" % write_line)
        
    print(i)
# for i in range (len(list_keys)):
#    print (list_keys[i],":",list_values[i])
    

#print("\nAmino_acids_percent:", X.get_amino_acids_percent()) 

file.close()
