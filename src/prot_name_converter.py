'''
Created on 11 Dec 2018

@author: mate

parse human-mouse gene table and use it to map mouse IDs to human ones
'''

def file_importer(relPath, methodS = "r", encodeS = None):
  """open a file specified by the relPath variable, 
  even if its in another directory"""
  import os
  scriptDir = os.path.dirname(os.path.realpath(__file__)) # absolute dir this script is in
  absFilePath = os.path.join(scriptDir, relPath)
  # print(encodingS)
  inpF = open(absFilePath, methodS, encoding = encodeS)
  return inpF


def homologue_parser():
  """Open file from http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt and extract human-mouse gene pairs. 
  Return in a dict with human gene names as keys and mouse gene names as values."""
  
  relPath = "data/HOM_MouseHumanSequence.rpt"
  
  inpF = file_importer(relPath)
  
  fileDict = {}
  compDict = {}
  
  next(inpF)
  
  # populate filedict with homology group IDs as keys, and values as a list of lists. the first element is the list of human gene names, second is the mouse gene names
  for inpLine in inpF:
    inpList = inpLine.split("\t")
    
    
    if inpList[0] in fileDict:
      if inpList[2] == "9606": fileDict[inpList[0]][0].append(inpList[3])
      elif inpList[2] == "10090": fileDict[inpList[0]][1].append(inpList[3])  
      else:
        print(inpList)
        raise ValueError
    
    else:
      if inpList[2] == "9606": fileDict[inpList[0]] = [[inpList[3]],[]] 
      elif inpList[2] == "10090": fileDict[inpList[0]] = [[],[inpList[3]]]
      else:
        print(inpList)
        raise ValueError      
      
  
  inpF.close()
  # print(fileDict)
  
  for keyS in fileDict:
    curVal = fileDict[keyS]
    for keyI in curVal[0]:
      if keyI == "": continue
      if keyI in compDict:
        if compDict[keyI] == []: del compDict[keyI]
        else:
          print(curVal)
          print(compDict[keyI])
          raise ValueError
      compDict[keyI] = curVal[1]
      
  for keyS, valueS in list(compDict.items()):
    if valueS == []: del compDict[keyS]
    
  return compDict
    
#   for keyS, valueS in compDict.items():
#     print(keyS, valueS) 

def gene_list_reader():
  """read in .csv file with entrez gene IDs from Caludio's lab and return a list of those IDs"""
  
  relPath = "data/genes_met_modelling_human.csv"
  
  geneL = []
  with file_importer(relPath, encodeS = "utf-8-sig") as inpF:
    for inpLine in inpF:
      inpI = inpLine.strip("\n'").split(".")[0]
      if inpI not in geneL: geneL.append(inpI)
      
  return geneL
    
def prot_id_converter(protList, orgnID = "10090", inpDB = "uniprotaccession", outDB = "genbankproteingi"):
  """copied from tools.py in the ed project.
  
  take in a list of uniprot entry names 
  and convert them to protein IDs using biodbnet.
  return a list of protein IDs. 
  bioDBnet often takes time to load so this might take several minutes to complete.
  
  For orgnID use 10090 for mouse and 9606 for human.
  
  Human gene symbols are all caps: VAV1
  Mouse gene symbols are first caps: Vav1
  
  inpDB is input database type, use "uniprotaccession" for uniprot entry IDs or "genesymbol" for gene symbols
  outDB is output database type, use "geneid" for Gene ID or "genbankproteinaccession" or "genbankproteingi" or "refseqproteingi" for those
 
  """
  import urllib.request, json
  urlStr = "http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=" + inpDB + "&inputValues=" + ",".join(protList) + "&outputs=" + outDB + "&taxonId=" + orgnID  
  print("connecting to biodbnet. This might take a while...")
  uParsed = urllib.request.urlopen(urlStr)  
  print("connection successful")
  responseJson = uParsed.read()
  # print(responseJson)
  parsedJson = json.loads(responseJson.decode('utf-8'))
  # print parsedJson
  # parsedJson = [{u'Gene ID': u'54196', u'InputValue': u'Q8CCS6'}, {u'Gene ID': u'99982', u'InputValue': u'Q6ZQ88'}]
  # parsedJson = [{u'GenBank Protein Accession': u'BAC27741//Q8CCS6//EDL36322//AAH55866//NP_062275//XP_006519335//AAC00210////EDL36323', u'InputValue': u'Q8CCS6'}, {u'GenBank Protein Accession': u'AAH19417//XP_006539394//XP_006539393//NP_598633//AAH59885//CBY79415//CBY88367////XP_006539392//EDL29935//Q6ZQ88//BAC97980', u'InputValue': u'Q6ZQ88'}]
  # parsedJson = [{u'GenBank Protein GI': u'9506945//148704376//46396417//26328001//33585929//2351846////148704375//568988212', u'InputValue': u'NP_062275'}, {u'GenBank Protein GI': u'51315882//18044445//224994233////568932208//317440660//37589595//568932212//315003691//148697988//37360004//568932210', u'InputValue': u'NP_598633'}]
  # parsedJson = [{u'RefSeq Protein GI': u'//6005942', u'InputValue': u'P55072'}, {u'RefSeq Protein GI': u'530368795////46488944//767919614//578804849//31455611', u'InputValue': u'P43403'}, {u'RefSeq Protein GI': u'7108367//384551646//768003854//530425424//384551649', u'InputValue': u'P15498'}, {u'RefSeq Protein GI': u'767904317//112789546//767904319////112789548//767904315', u'InputValue': u'P06239'}, {u'RefSeq Protein GI': u'4502671//', u'InputValue': u'P07766'}, {u'RefSeq Protein GI': u'7108367//384551646//768003854//530425424//384551649', u'InputValue': u'P15498'}, {u'RefSeq Protein GI': u'767904317//112789546//767904319////112789548//767904315', u'InputValue': u'P06239'}, {u'RefSeq Protein GI': u'767910875//37595565//4557431////767910873', u'InputValue': u'P20963'}, {u'RefSeq Protein GI': u'4502671//', u'InputValue': u'P07766'}, {u'RefSeq Protein GI': u'530368795////46488944//767919614//578804849//31455611', u'InputValue': u'P43403'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'767910875//37595565//4557431////767910873', u'InputValue': u'P20963'}, {u'RefSeq Protein GI': u'767910875//37595565//4557431////767910873', u'InputValue': u'P20963'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'P32577'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'//768033853//4507909', u'InputValue': u'P42768'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'767985679//767985664////767985659//767985662//767985670//20149528//767985683//578827539//767985657//767985674//767985677//767985681//767985668', u'InputValue': u'O43586'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'767985679//767985664////767985659//767985662//767985670//20149528//767985683//578827539//767985657//767985674//767985677//767985681//767985668', u'InputValue': u'O43586'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'//29725609//41327734//41327732//41327736', u'InputValue': u'P00533'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'744066863//16753212', u'InputValue': u'O75563'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'EBI-9974954'}, {u'RefSeq Protein GI': u'744066863//16753212', u'InputValue': u'O75563'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'Q3UND0'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'P97814'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'Q99M15'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}]
  # parsedJson = [{u'UniProt Accession': u'Q3TIM2//Q3TFH9//Q3TXN9//Q01853//Q6PI18//Q8BSR6//Q8CEG4', u'InputValue': u'Vcp'}, {u'UniProt Accession': u'P43404//P97455//Q80VV2//Q8CHJ3', u'InputValue': u'Zap70'}, {u'UniProt Accession': u'P27870//Q8BTV7', u'InputValue': u'Vav1'}, {u'UniProt Accession': u'Q91X65//Q62320//E9Q696//Q61794//P06240//Q61795', u'InputValue': u'Lck'}, {u'UniProt Accession': u'P22646', u'InputValue': u'Cd3e'}, {u'UniProt Accession': u'Q9D3G3//P29020//P24161', u'InputValue': u'Cd247'}, {u'UniProt Accession': u'Q03143//Q8VCW1//P41241//Q80WU4', u'InputValue': u'Csk'}, {u'UniProt Accession': u'Q64260//Q6GU11//P15116', u'InputValue': u'Cdh2'}, {u'UniProt Accession': u'P70424//Q61525//Q6ZPE0', u'InputValue': u'Erbb2'}, {u'UniProt Accession': u'E9QPE2//P05622', u'InputValue': u'Pdgfrb'}, {u'UniProt Accession': u'Q3UFB7', u'InputValue': u'Ntrk1'}, {u'UniProt Accession': u'A6H6U3//Q9Z2A0//Q9R1D8//Q9R215', u'InputValue': u'Pdpk1'}, {u'UniProt Accession': u'P16590//Q920Z3//Q9R264//P16882//Q61653//Q6DI66//Q80W86//Q8R1M5', u'InputValue': u'Ghr'}, {u'UniProt Accession': u'P70315', u'InputValue': u'Was'}, {u'UniProt Accession': u'Q60631//Q61240', u'InputValue': u'Grb2'}, {u'UniProt Accession': u'Q9EP98//Q01279', u'InputValue': u'Egfr'}, {u'UniProt Accession': u'Q8BK74//Q3UND0//Q9Z2K4', u'InputValue': u'Skap2'}, {u'UniProt Accession': u'Q4V9R4//P97814', u'InputValue': u'Pstpip1'}, {u'UniProt Accession': u'Q99M15//Q6GTF6//Q9Z189', u'InputValue': u'Pstpip2'}, {u'UniProt Accession': u'P08032//B2RWX6//P97502', u'InputValue': u'Spta1'}, {u'UniProt Accession': u'Q3U527//Q8CEA1//P22682', u'InputValue': u'Cbl'}]
  # parsedJson = [{u'KEGG Gene ID': u'mmu:16656', u'InputValue': u'A2A884'}, {u'KEGG Gene ID': u'mmu:11737', u'InputValue': u'O35381'}, {u'KEGG Gene ID': u'mmu:19763', u'InputValue': u'O35730'}]
  # parsedJson = [{u'InputValue': u'Q5SWU9', u'Gene Symbol': u'Acaca'}, {u'InputValue': u'Q8VDD5', u'Gene Symbol': u'Myh9'}, {u'InputValue': u'Q3T9S7', u'Gene Symbol': u'Pcx'}, {u'InputValue': u'B2RRX1', u'Gene Symbol': u'Actb'}, {u'InputValue': u'Q71LX8', u'Gene Symbol': u'Hsp90ab1'}, {u'InputValue': u'B2RRE2', u'Gene Symbol': u'Myo18a'}, {u'InputValue': u'Q3U2W2', u'Gene Symbol': u'Mybbp1a'}, {u'InputValue': u'Q3TII3', u'Gene Symbol': u'Eef1a1'}, {u'InputValue': u'P99024', u'Gene Symbol': u'Tubb5'}, {u'InputValue': u'E9QAS3', u'Gene Symbol': u'Ptpn22'}, {u'InputValue': u'Q99MR8', u'Gene Symbol': u'Mccc1'}, {u'InputValue': u'Q3THE2', u'Gene Symbol': u'Myl12b'}, {u'InputValue': u'D3YZ62', u'Gene Symbol': u'Myo5a'}, {u'InputValue': u'Q3UGC8', u'Gene Symbol': u'Pcca'}, {u'InputValue': u'Q6S385', u'Gene Symbol': u'Plec'}, {u'InputValue': u'B2RTP7', u'Gene Symbol': u'Krt2'}, {u'InputValue': u'B1AQ77', u'Gene Symbol': u'Krt15'}, {u'InputValue': u'D3Z6I8', u'Gene Symbol': u'Tpm3'}, {u'InputValue': u'B2RTM0', u'Gene Symbol': u'Hist2h4'}, {u'InputValue': u'Q8K0Z5', u'Gene Symbol': u'Tpm3'}, {u'InputValue': u'Q3TNH0', u'Gene Symbol': u'Tmpo'}, {u'InputValue': u'Q3TIG9', u'Gene Symbol': u'Myl6'}, {u'InputValue': u'D2KHZ9', u'Gene Symbol': u'GAPDH'}, {u'InputValue': u'Q6P5D8', u'Gene Symbol': u'Smchd1'}, {u'InputValue': u'Q4FZG4', u'Gene Symbol': u'Flna'}, {u'InputValue': u'F1DGF6', u'Gene Symbol': u'Prkcd'}, {u'InputValue': u'Q3TFG3', u'Gene Symbol': u'Eif4a1'}, {u'InputValue': u'B2RPX1', u'Gene Symbol': u'Iqcd'}, {u'InputValue': u'Q8BQ35', u'Gene Symbol': u'Sptbn1'}, {u'InputValue': u'E0CZ27', u'Gene Symbol': u'H3f3a'}, {u'InputValue': u'Q9CR57', u'Gene Symbol': u'Rpl14'}, {u'InputValue': u'Q0VG47', u'Gene Symbol': u'Hnrnpa3'}, {u'InputValue': u'Q8C553', u'Gene Symbol': u'Lmnb1'}, {u'InputValue': u'Q3T9U9', u'Gene Symbol': u'Rpl3'}, {u'InputValue': u'Q3KQJ4', u'Gene Symbol': u'Hspa8'}, {u'InputValue': u'Q3U7D2', u'Gene Symbol': u'Rpl15'}, {u'InputValue': u'A0PJE6', u'Gene Symbol': u'Pccb'}, {u'InputValue': u'Q68FG3', u'Gene Symbol': u'Spty2d1'}, {u'InputValue': u'Q0VB76', u'Gene Symbol': u'Gzmc'}, {u'InputValue': u'Q32P04', u'Gene Symbol': u'Krt5'}, {u'InputValue': u'D3Z6F5', u'Gene Symbol': u'Atp5a1'}, {u'InputValue': u'Q3U0I3', u'Gene Symbol': u'Cct3'}, {u'InputValue': u'Q3TJZ1', u'Gene Symbol': u'Eef2'}, {u'InputValue': u'Q3UI57', u'Gene Symbol': u'Mcm3'}]

  # 2018:
  # b'[\n    {\n        "InputValue": "100",\n        "Gene Symbol": "ADA"\n    },\n    {\n        "InputValue": "10005",\n        "Gene Symbol": "ACOT8"\n    }

  if "Gene Symbol" in parsedJson[0]:
    return gene_symbol_wrangler(parsedJson)
  
  if "Gene ID" in parsedJson[0]:
    return gene_ID_wrangler(parsedJson)

  if "RefSeq Protein GI" in parsedJson[0]:
    return refseq_gi_wrangler(parsedJson)
  
                         
  else: 
    print("was expecting gene symbols but got something else:")
    print(parsedJson)
    raise ValueError

def gene_symbol_wrangler(inpAcc):
  """rewritten in 2018, not the same as in tools.py
  
  to be called by prot_id_converter
  Take in a loaded json file which has gene symbols as results.
  Return a list with with the gene names from the query"""
  
  print("processing gene symbols")
  
  resD = {}
  
  for convI in inpAcc:
    keyI = convI["InputValue"]
    valueI = convI["Gene Symbol"]
    resD[keyI] = valueI

  return resD

def gene_ID_wrangler(inpAcc):
  """rewritten in 2018, not the same as in tools.py
  
  to be called by prot_id_converter
  Take in a loaded json file which has gene symbols as results.
  Return a list with with the gene names from the query"""
  
  print("processing gene symbols")
  
  resD = {}
  
  for convI in inpAcc:
    keyI = convI["InputValue"]
    valueI = convI["Gene ID"]
    resD[keyI] = valueI

  return resD

def refseq_gi_wrangler(inpAccessions):
  """to be called by prot_id_converter.
  Take in a loaded json entry from biodbnet db2db which has refseq protein accession ID results. 
  pick the smallest GI number for each entry (this should be a nice, full length sequence) and return a list of protein GIs. 
  These are unique identifiers (UIDs) that can be used by entrez.epost, entrez.efetch and such things.
  """
  print("processing GenBank Protein GIs")
  resD = {}
  for inpAccD in inpAccessions:
    queryL = inpAccD["RefSeq Protein GI"].split("//")
    minQ = 999999999999999 # this should be big enough
    for queryI in queryL:
      if queryI == "" or queryI == "-": continue
      # curQuery = int(queryI.encode("ascii","ignore"))
      curQuery = int(queryI)
      if curQuery < minQ:
        minQ = curQuery
    if minQ == 999999999999999: 
      print("GI not found in this query:")
      print(inpAccD)
    else:
      resD[inpAccD["InputValue"]] = minQ   

  return resD 

def gene_list_converter():
  """take in a list of human gene IDs and return a dict of gene ID: gene name key-value pairs"""
  geneL = gene_list_reader()
  humD = prot_id_converter(geneL, "9606", "geneid", "genesymbol")
  return humD


def prot_sequence_finder(protL):
  """take in a gene symbol list and return a dict with names as keys and protein sequence strings as values """
  
  idDict = prot_id_converter(protL, "9606", inpDB = "genesymbol",outDB="refseqproteingi")
  seqD = prot_entrez_fetch(idDict, retM="gb", retT="fasta")
  
  protD = {}
  
  for keyS, valueS in idDict.items():
    protD[keyS] = seqD[valueS]
  
  return protD
  
  
def entrez_fasta_parser(handleFasta):
  """Split and organise a single efetch query which is a string containing fasta sequences
  into a list of fasta sequences. Return the list."""
  fullList = handleFasta.read().split("\n") 
  resL = []
  seqFlag = False
  for fullLine in fullList:
    if fullLine == "":
      seqFlag = False
      continue
    elif fullLine[0] == ">":
      resL.append(fullLine + "\n")
      seqFlag = True
    elif seqFlag:
      resL[-1] += fullLine     
  return resL

def prot_entrez_fetch(proteinDict, retM="text", retT="fasta"):
  """take in a list of protein GI accessions and return their corresponding fasta sequences as a list using Entrez.
  Each returned list item is the fasta header + new line + sequence.
  Do not give it more than 200 GIs at once as requested by entrez.
  
  GI accessions are deprecated. need to use something else.
  """
  from Bio import Entrez
  Entrez.email ="mate.ravasz@ed.ac.uk"
  for i in proteinDict.values():
    try:
      int(i) # test if really a list of UIDs
    except ValueError:
      print("was expecting UIDs like \"12345678\", but got something else instead:")
      print(i)
      raise
  inpList = list(proteinDict.items())
  proteinIntList = []
  for k in inpList:
    proteinIntList.append(k[1])
    
  proteinList = list(map(str, proteinIntList))
  # print proteinList
  
  print("connecting to Entrez...")
  requestR = Entrez.epost("protein",id=",".join(proteinList)) # send all UIDs as a single query to entrez. 
  resultO = Entrez.read(requestR)
  webEnv = resultO["WebEnv"]
  queryKey = resultO["QueryKey"]
  handleO = Entrez.efetch(db="protein",retmode=retM, rettype=retT, webenv=webEnv, query_key=queryKey) # retrieve all results in batch
  print("connection successful")
  if retT == "fasta":
    fastaL = entrez_fasta_parser(handleO)
    protD = {}
    for j in range(len(proteinIntList)):
      protD[proteinIntList[j]] = fastaL[j].split("\n")[1]
    return protD
    
  elif retM == "text" and retT == "gp": # use "gp" for genpept flatfile format, and "gb" for genbank flatfile for genes
    return handleO.read()
  else:
    print("this data format was not expected:")
    print("retmode: ", retM)
    print("rettype: ", retT)
    raise ValueError


def blaster(protSeq, orgnID):
  """take in an amino acid sequence and return the best matching gene name from the organism defined"""


def species_converter(): 
  """convert human gene names to mouse gene names.
  make a dict with human gene ID as key, and a dict as value
  where the keys are human gene names, and values are mouse gene names
  write these out to a csv file, with first column being the human gene ID
  second column the human gene symbol, third column the mouse gene symbol
  fourth column the mouse gene ID"""
  
  geneD = gene_list_converter()
  convertD = homologue_parser()
#   k = 0
#   for i,j in geneD.items():
#     print(i,j)
#     k += 1
#     if k == 10000: break
#   print("------------------------------------------------")
#   k = 0
#   for i,j in convertD.items():
#     print(i,j)
#     k += 1
#     if k == 100: break
  
  procD = {}
  
  for geneK, geneV in geneD.items():
    if geneV == "-": # handle missing gene names
      procD[geneK] = {"-": ["-"]}
      continue
    
    if geneV in convertD:
      if geneK in procD:
        print(geneK,geneV,procD[geneK])
        raise ValueError
      
      else: procD[geneK] = {geneV: convertD[geneV]}
      
    else: procD[geneK] = {geneV: ["-"]}
    
  mouseGeneL = []
  for keyS in procD.keys():
    for valueL  in procD[keyS].values(): # @UnusedVariable
      if valueL == ["-"]: continue
      for valueI in valueL:
        if valueI in mouseGeneL: continue
        mouseGeneL.append(valueI)
          
  mouseD = prot_id_converter(mouseGeneL, "10090", "genesymbol", "geneid")
  
  with file_importer("data/converted_gene_list.csv", "w") as outF:
    for keyS in procD.keys():
      outF.write(keyS + ",")
      for keyN, valueN in procD[keyS].items():
        outF.write(keyN + ",")
        if keyN == "-":
          outF.write("-,-\n")
          continue
        valL = []
        for valueI in valueN:
          if valueI in mouseD: valL.append(mouseD[valueI])
          else: valL.append("-")
          
          if valueI is valueN[-1]: outF.write(valueI + ",")
          else: outF.write(valueI + ";")
          
        for valI in valL:
          if valI is valL[-1]: 
            if "//" in valI: 
              valIL = valI.split("//")
              for valILI in valIL:
                if valILI is valIL[-1]: outF.write(valILI + "\n")
                else: outF.write(valILI + ";")
            else: outF.write(valI + "\n")
            
          else: 
            if "//" in valI: 
              valIL = valI.split("//")
              for valILI in valIL:
                outF.write(valILI + ";")
            else: outF.write(valI + ";")
  print("file written")
        
      
  
  
print(prot_sequence_finder(["COX7C","SLC2A14","CYP2D6","CYP3A4","LDHAL6A"]))
