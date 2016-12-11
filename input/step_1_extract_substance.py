import os
import sys
from Bio import Entrez
from Bio import Medline
import random
import numpy
import math

"""
This script can only download the part of "RN" in MEDLINE.
In the part of "RN", "0" is abandoned and in contrast "EC" and "Left" remain.
"""

def splitlist(inputlist, sublist_size):
	l = len(inputlist)
	outputlist = []
	if l <= sublist_size:
		outputlist.append(inputlist)
	else:
		s1 = l/int(sublist_size)
		s2 = l%int(sublist_size)
		for i in range(s1):
			outputlist.append(inputlist[i*sublist_size:(i+1)*sublist_size])
		if s2 != 0:
			outputlist.append(inputlist[s1*sublist_size:])
	return outputlist

def downloadMeshTerm(termToSearch):
	# Download resources (the part of "MH" in MEDLINE). Obviously, these paper will have "MeshTerm" as well.
	# Pay attention: all of the MeshTerm including 16 subcategories will be retained.
	# len(PMIDlist) == len(PMIDMeSHtermDic)
	Entrez.email = "howwhatppqq@163.com"
	batch_size = 10000 # batch_size could optional. It indicates the number of fetched documents each time via Entrez.efetch
	#Section1: Fetch the documents by searching PubMed 
	search_handle = Entrez.esearch(db="pubmed", term=termToSearch,retmax=500000)
	search_results = Entrez.read(search_handle)
	count_1 = int(search_results["Count"])
	idlist = search_results["IdList"]
	count_2 = len(idlist)
	sublist_size = 10000 #10000 is the max of len(search_results1["IdList"]). So this number should be constant
	PMIDlist = []
	PMIDMeSHtermDic = {}
	# The format of an element in PMIDMeSHtermDic is like:
	# PMIDMeSHtermDic["PMID"]=[[items of "MeshTerm"],[items of "MeshTerm" with '*']]
	# Here, at first we analyses all of the MeshTerm;
	# In the next step, we can classify the significant terms and perform association rule analysis towards specific class(e.g. compounds and protein).
	#outf = open("temp.txt","w")
	for subidlist in splitlist(idlist, sublist_size):
		fetch_handle = Entrez.efetch(db="pubmed", id=subidlist, rettype="medline", retmode="text")
		records = Medline.parse(fetch_handle)
		for record in records:
			if "MH" in record.keys():
				PMIDlist.append(record['PMID'])
				PMIDMeSHtermDic.update({record['PMID']:[[],[]]})
				for subrecord2 in record["MH"]:
					if "/" not in subrecord2:
						if "*" in subrecord2:
							newitem = ""
							for i in subrecord2:
								if i != "*":
									newitem = newitem + i
							itemTempWithStar = newitem
							itemTemp = newitem
						else:
							itemTemp = subrecord2
							itemTempWithStar = ""
					else:
						itemTemp = subrecord2.split("/")[0]
						if "*" in subrecord2:
							itemTempWithStar = subrecord2.split("/")[0]
						else:
							itemTempWithStar = ""
					if itemTemp not in PMIDMeSHtermDic[record['PMID']][0]:
						if "*" in itemTemp:
							newitem2 = ""
							for i in itemTemp:
								if i != "*":
									newitem2 = newitem2 + i
						else:
							newitem2 = itemTemp
						if newitem2 not in PMIDMeSHtermDic[record['PMID']][0]:
							PMIDMeSHtermDic[record['PMID']][0].append(newitem2)
					if (itemTempWithStar != "") and (itemTempWithStar not in PMIDMeSHtermDic[record['PMID']][1]):
						if "*" in itemTempWithStar:
							newitem3 = ""
							for k in itemTempWithStar:
								if k != "*":
									newitem2 = newitem2 + k
						else:
							newitem3 = itemTempWithStar
						if newitem3 not in PMIDMeSHtermDic[record['PMID']][1]:
							PMIDMeSHtermDic[record['PMID']][1].append(newitem3)
	return (PMIDlist,PMIDMeSHtermDic)

def downloadSubstance(termToSearch):
	# Download resources (the part of "RN" in MEDLINE, the paper must have "MeshTerm".  
	# If a article has "MH" and doesn't have "RN", it will also retained for further calculation) form PubMed using the search keywords
	# Pay attention: only the "Left" or "EC" in 'RN' will be retained , others (including "0" and No "RN") will be empty.
	# len(PMIDlist) == len(PMIDSubstanceDic)
	Entrez.email = "howwhatppqq@163.com"
	batch_size = 10000 # batch_size could optional. It indicates the number of fetched documents each time via Entrez.efetch
	#Section1: Fetch the documents by searching PubMed 
	search_handle = Entrez.esearch(db="pubmed", term=termToSearch,retmax=500000)
	search_results = Entrez.read(search_handle)
	count_1 = int(search_results["Count"])
	idlist = search_results["IdList"]
	count_2 = len(idlist)
	sublist_size = 10000 #10000 is the max of len(search_results1["IdList"]). So this number should be constant
	PMIDlist = []
	PMIDSubstanceDic = {}
	# The format of an element in PMIDSubstanceDic is like:
	# PMIDSubstanceDic["PMID"]=[[items of "Left"],[items if "EC"]]
	# Pay attention the sequence of lists in the list:
	# The first one belongs to "Left", and the second one belongs to "EC"
	# It will be extended when dealing with the MeshTerms and the abstracts.
	for subidlist in splitlist(idlist, sublist_size):
		fetch_handle = Entrez.efetch(db="pubmed", id=subidlist, rettype="medline", retmode="text")
		records = Medline.parse(fetch_handle)
		for record in records:
			if "MH" in record.keys():
				if "RN" in record.keys():
					categoryList = []
					for subrecord1 in record["RN"]:
						categoryList.append( subrecord1.split("(")[0][:-1])
					Tag = 0
					for item in categoryList:
						if item != "0":
							Tag = 1
							break
					if Tag == 1:
						PMIDlist.append(record['PMID'])
						PMIDSubstanceDic.update({record['PMID']:[[],[]]})
						for subrecord2 in record["RN"]:
							contentTemp = subrecord2.split("(")
							name = "".join(contentTemp[1:])[:-1]
							if contentTemp[0][:-1] == "0":
								pass
							elif "EC " in contentTemp[0][:-1]:
								if len(name) != 0:
									PMIDSubstanceDic[record['PMID']][1].append(name)
							else:
								if len(name) != 0:
									PMIDSubstanceDic[record['PMID']][0].append(name)
					else:
						PMIDlist.append(record['PMID'])
						PMIDSubstanceDic.update({record['PMID']:[[],[]]})
				else:
					PMIDlist.append(record['PMID'])
					PMIDSubstanceDic.update({record['PMID']:[[],[]]})
	return (PMIDlist,PMIDSubstanceDic)

def main(arg = sys.argv):
	if len(arg) != 3:
		print >> sys.stderr,"*.py inputFile1(including the name list of drugs) inputFile2(para.txt)"
		sys.exit(1)
	else:
		inf = open(arg[1],"r")
		queryList = []
		for line in inf:
			queryList.append(line.strip())
		inf.close()
		inf2 = open(arg[2],"r")
		for line in inf2:
			content = line.strip().split("\t")
			if content[0] == "endDate":
				endDate = content[1]
				break
			else:
				pass
		inf2.close()
		timeInterval = "(0001/01/01[PDAT] : " + endDate + "[PDAT])"
		for drugName in queryList:
			#searchTerm = drugName.lower()
			searchTerm_temp = drugName.lower()
			searchTerm = drugName

			#outf5 = open(searchTerm_temp+"_PMID_MeShterms.txt","w")
			outf5 = open(searchTerm_temp+"_PMID_Substances.txt","w")
			keyword = drugName
			searchTerm_1 = searchTerm + " AND drug interaction[MeSH Terms] AND "+timeInterval
			#searchTerm_1 = "((((((((Drug Interactions[MeSH Terms]) OR Drug Agonism[MeSH Terms]) OR Drug Partial Agonism[MeSH Terms]) OR Drug Antagonism[MeSH Terms]) OR Drug Inverse Agonism[MeSH Terms]) OR Drug Synergism[MeSH Terms]) OR Food-Drug Interactions[MeSH Terms]) OR Herb-Drug Interactions[MeSH Terms]) AND " + searchTerm+"[MeSH Terms] AND "+timeInterval
			searchTerm_2 = searchTerm+"[MeSH Terms] AND "+timeInterval
			searchTerm_2 = searchTerm+timeInterval
			print "Download %s"%(searchTerm)
			#result_1 = downloadMeshTerm(searchTerm_1)
			result_1 = downloadSubstance(searchTerm_1)
			print "First Part is Done"
			result_1_PMIDList = result_1[0]
			result_1_dict = result_1[1]
			
			outf5.write("PMID\tSubstances\n")
			for pmidTemp in result_1_PMIDList:
				for submeshTerms in result_1_dict[pmidTemp][0]:
					outf5.write(pmidTemp+"\t"+submeshTerms+"\n")
				for submeshTerms in result_1_dict[pmidTemp][1]:
					outf5.write(pmidTemp+"\t"+submeshTerms+"\n")
			outf5.close()
			
			#result_2 = downloadMeshTerm(searchTerm_2)
			result_2 = downloadSubstance(searchTerm_2)
			print "Second Part is Done"

			result_2_dict = result_2[1]
			with open(searchTerm_temp+"_result_1_compounds.txt", 'w') as f1, open(searchTerm_temp+"_result_1_proteins.txt", 'w') as f2:
				for key, value in result_1_dict.iteritems():
					s = ''
					while drugName in value[0]:
						value[0].remove(drugName)
					for x in value[0]:
						s = s + x + '~'
					if (s!='') and (s!='~'):
						f1.write('%s\n' % s)
					s = ''
					while drugName in value[1]:
						value[1].remove(drugName)
					for x in value[1]:
						s = s + x + '~'
					if (s!='') and (s!='~'):
						f2.write('%s\n' % s)
			with open(searchTerm_temp+"_result_2_compounds.txt", 'w') as f1, open(searchTerm_temp+"_result_2_proteins.txt", 'w') as f2:
				for key, value in result_2_dict.iteritems():
					s = ''
					if key not in result_1_PMIDList:
						while drugName in value[0]:
							value[0].remove(drugName)
						# cnt = 0
						for x in value[0]:
							if (x!=''):
								s = s + x + '~'
								# cnt += 1
						if (s!='') and (s!='~'): # and (cnt >=1):
							f1.write("%s\n" % s)
						s = ''
						while drugName in value[1]:
							value[1].remove(drugName)
						# cnt = 0
						for x in value[1]:
							if (x!=''):
								s = s + x + '~'
								# cnt += 1
						if (s!='') and (s!='~'): # and (cnt >=1):
							f2.write("%s\n" % s)
			
if __name__ == "__main__":
	main()
