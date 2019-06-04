'''
Created on 1 Feb 2019

@author: mate

- take in uniprot IDs from processed proteinGroups file.
- convert mouse uniprot IDs to mouse gene symbols
- convert mouse gene symbols to human gene symbols
- convert mouse gene symbols with no human counterpart by either looking fo the same gene name in humans, or blast-ing the sequence for a human homologue
- convert human gene names to HGNC IDs

some code here copied over from prot_name_converter in ed repo
'''


def file_importer(relPath, methodS = "r", encodeS = None):
  """open a file specified by the relPath variable, 
  even if its in another directory"""
  import os
  scriptDir = os.path.dirname(os.path.realpath(__file__)) # absolute dir this script is in
  absFilePath = os.path.join(scriptDir, relPath)
  # print(encodingS)
  if methodS == "w":
    filePartL = absFilePath.split(".")
    if len(filePartL) > 1:
      i = 1
      filePartL = absFilePath.split(".")
      
      while os.path.exists(filePartL[0] + "-" + str(i) + "." + filePartL[1]):
        i += 1
      
      absFilePath = filePartL[0] + "-" + str(i) + "." + filePartL[1]

    else:
      print("no extension given in file name. please add an extension to your file, like .txt or something, as this will not do: " + relPath)
  inpF = open(absFilePath, methodS, encoding = encodeS)
  return inpF


def homologue_parser():
  """Open file from http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt and extract human-mouse gene pairs. 
  Return in a dict with human gene names as keys and a list of? mouse gene names as values."""
  
  relPath = "data/HOM_MouseHumanSequence.rpt"
  print("Creating homologue dict...")
  
  inpF = file_importer(relPath)
  
  fileDict = {}
  compDict = {}
  
  next(inpF)
  
  # populate filedict with homology group IDs as keys, and values as a list of lists. the first element is the list of mouse gene names, second is the human gene names
  for inpLine in inpF:
    inpList = inpLine.split("\t")
    
    
    if inpList[0] in fileDict:
      if inpList[2] == "9606": fileDict[inpList[0]][1].append(inpList[3])
      elif inpList[2] == "10090": fileDict[inpList[0]][0].append(inpList[3])  
      else:
        print(inpList)
        raise ValueError
    
    else:
      if inpList[2] == "9606": fileDict[inpList[0]] = [[],[inpList[3]]] 
      elif inpList[2] == "10090": fileDict[inpList[0]] = [[inpList[3]],[]]
      else:
        print(inpList)
        raise ValueError      
      
  
  inpF.close()
  
  for keyS in fileDict:
    curVal = fileDict[keyS]
    for keyI in curVal[0]:
      if keyI == "": continue
      if keyI in compDict:
        if compDict[keyI] == []: 
          compDict[keyI] = curVal[1]
        else: 
          for curI in curVal[1]:
            if curI not in compDict[keyI]: compDict[keyI].append(curI)
      else: compDict[keyI] = curVal[1]
      
  for keyS, valueS in list(compDict.items()):
    if valueS == []: del compDict[keyS]
    
  print("Done.\n")
  return compDict
    
#   for keyS, valueS in compDict.items():
#     print(keyS, valueS) 

def gene_list_reader(splitLines = True):
  """read in a processed proteingroups.csv file and return the gene names as a list"""
  
  relPath = "data/24H_T4_original_proteinGoups_avg.csv"
  
  geneL = []
  with file_importer(relPath) as inpF:
    next(inpF)
    for inpLine in inpF:
      inpI = inpLine.rstrip("\n").split(",")[0]
      if splitLines and ";" in inpI:
        inpL = inpI.split(";")
        for inpLI in inpL:
          if "-" in inpLI or "__" in inpI: continue
          if inpLI not in geneL: geneL.append(inpLI)
      elif inpI not in geneL: 
        if "__" in inpI: continue
        geneL.append(inpI.rstrip("\n"))
      
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
  
  def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
        
  finD = {}
  for chunkL in chunks(protList,400):
  
    urlStr = "http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=" + inpDB + "&inputValues=" + ",".join(chunkL) + "&outputs=" + outDB + "&taxonId=" + orgnID  
    # print("connecting to biodbnet. This might take a while...")
    print(".", end = "")
    uParsed = urllib.request.urlopen(urlStr)  
    # print("connection successful")
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
      finD.update(gene_symbol_wrangler(parsedJson))
    
    elif "Gene ID" in parsedJson[0]:
      finD.update(gene_ID_wrangler(parsedJson))
  
    elif "RefSeq Protein GI" in parsedJson[0]:
      finD.update(refseq_gi_wrangler(parsedJson))
      
    elif "UniProt Accession" in parsedJson[0]:
      finD.update(uniprot_acc_wrangler(parsedJson))
      
    elif("HGNC ID") in parsedJson[0]:
      finD.update(hgnc_id_wrangler(parsedJson))
  
                         
    else: 
      print("was expecting gene symbols but got something else:")
      print(parsedJson)
      raise ValueError
    
  return finD

def hgnc_id_wrangler(inpAcc):
  
  resD = {}
  
  for convI in inpAcc:
    keyI = convI["InputValue"]
    valueL = convI["HGNC ID"].split("//")
    resD[keyI] = []
    for valueI in valueL: resD[keyI].append(valueI)

  return resD


def uniprot_acc_wrangler(inpAcc):
  
  resD = {}
  
  for convI in inpAcc:
    keyI = convI["InputValue"]
    valueL = convI["UniProt Accession"].split("//")
    resD[keyI] = []
    for valueI in valueL: resD[keyI].append(valueI)

  return resD

def gene_symbol_wrangler(inpAcc):
  """rewritten in 2018, not the same as in tools.py
  
  to be called by prot_id_converter
  Take in a loaded json file which has gene symbols as results.
  Return a list with with the gene names from the query"""
  
  # print("processing gene symbols")
  
  resD = {}
  
  for convI in inpAcc:
    keyI = convI["InputValue"]
    valueL = convI["Gene Symbol"].split("//")
    resD[keyI] = []
    for valueI in valueL: resD[keyI].append(valueI)

  return resD

def gene_ID_wrangler(inpAcc):
  """rewritten in 2018, not the same as in tools.py
  
  to be called by prot_id_converter
  Take in a loaded json file which has gene symbols as results.
  Return a list with with the gene names from the query"""
  
  # print("processing gene symbols")
  
  resD = {}
  
  for convI in inpAcc:
    keyI = convI["InputValue"]
    valueL = convI["Gene ID"].split("//")
    resD[keyI] = []
    for valueI in valueL: resD[keyI].append(valueI)

  return resD

def refseq_gi_wrangler(inpAccessions):
  """to be called by prot_id_converter.
  Take in a loaded json entry from biodbnet db2db which has refseq protein accession ID results. 
  pick the smallest GI number for each entry (this should be a nice, full length sequence) and return a list of protein GIs. 
  These are unique identifiers (UIDs) that can be used by entrez.epost, entrez.efetch and such things.
  """
  # print("processing GenBank Protein GIs")
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
  """take in a list of mouse uniprot accessions and return a dict of accession: gene name key-value pairs"""
  geneL = gene_list_reader()

  print("converting uniprot IDs to gene names...")
  humD = prot_id_converter(geneL, "10090", "uniprotaccession", "genesymbol")
    
  print("\nDone.\n")
  return humD


def prot_sequence_finder(protL):
  """take in a gene symbol list and return a dict with names as keys and protein sequence strings as values """
  print("\nconverting missing mouse gene names to refseq identifiers...")
  idDict = prot_id_converter(protL, "10090", inpDB = "genesymbol",outDB="refseqproteingi")
  print("\nDone.\n")
  print("fetching {} sequences...".format(len(idDict)))
  seqD = prot_entrez_fetch(idDict, retM="gb", retT="fasta")
  
  protD = {}
  
  for keyS, valueS in idDict.items():
    protD[keyS] = seqD[valueS]
  
  print("\nDone.\n")
  
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


def blaster(protSeq, orgnID = "Homo sapiens"):
  """take in an amino acid sequence and return the best matching gene name from the organism defined"""
  
  from Bio.Blast.NCBIWWW import qblast
  from Bio.Blast import NCBIXML
  from sys import exit
  
  print("\nconnecting to BLAST server. this will take some time...")
  i = 1
  while i < 4: # BLAST sometimes returns empty results. if so, try once more, it happens quite rarely and resending the query seems to fix it.
    print("attempt number " + str(i))
    i += 1
    resX = qblast("blastp","refseq_protein", protSeq, entrez_query = orgnID + "[organism]", descriptions = 100, alignments = 100)
    resO = NCBIXML.read(resX)
    if resO.descriptions != []: break 
  if resO.descriptions == []: 
    print("connection unsuccessful. The BLAST server is acting up. Try again later.")
    exit(0)
  
  else: print("connection successful")
    
  print(resO.descriptions[0])
  descO = resO.descriptions[0]
  if descO.e < 1e-137: # set identification threshold here. 0.01 still returns hits in most cases, 1e-140 is based on ptpn22 similarity between mouse and human
    descID = descO.title.split("|")[3]
    if "." in descID: return descID.split(".")[0]
    else: return descID
      
  else: return "-"
    

def blast_backup_parser(fName):
  """blaster tends to fail often leading to the whole script crashing when only part of the job is finished.
  This script helps to pick up the pieces and read in the results before the crash, greatly shortening run times."""
  
  blastD = {}
  with file_importer(fName) as blastF:
    for blastLine in blastF:
      blastL = blastLine.rstrip().split(",")
      blastD[blastL[0]] = blastL[1]
  
  return blastD
    

def species_converter(): 
  """convert human gene names to mouse gene names.
  make a dict with human gene ID as key, and a dict as value
  where the keys are human gene names, and values are mouse gene names
  write these out to a csv file, with first column being the human gene ID
  second column the human gene symbol, third column the mouse gene symbol
  fourth column the mouse gene ID"""
  
  geneD = gene_list_converter() # dict with keys as mouse uniprot IDS, and values as mouse gene symbols

  convertD = homologue_parser() # dict with human gene symbols as keys, and mouse gene symbols as values

  
  procD = {} # key is human gene ID, value is a dict with human gene symbol as key, mouse gene symbol as value
  missL = []
  missUpperL = []
  missIDL = []
  hitCount = 0
  
  print("mapping mouse gene names to human homologues...")
  
  for geneK, geneV in geneD.items(): # geneK: mouse uniprot accession, geneV = mouse gene symbols
    if geneV == ["-"]: # handle missing gene symbols
      procD[geneK] = {"-": ["-"]}
      continue
    
    for geneVI in geneV:
      if type(geneVI) == type(["hi"]): print(geneVI)
      if geneVI in convertD:
        if geneK in procD:
          if geneVI in procD[geneK]:
            print(geneK,geneVI,procD[geneK])
            raise ValueError
          else:
            procD[geneK][geneVI] = convertD[geneVI]
            # print(convertD[geneVI])
      
        else: 
          procD[geneK] = {geneVI: convertD[geneVI]} # mouse uniprot accession: mouse gene symbol: human gene symbol
          # print(procD[geneK])
          hitCount += 1
        
      
      else: 
        if geneVI in missL: print(geneVI)
        missL.append(geneVI) # collect gene symbols for which no human homologue was found
        missUpperL.append(geneVI.upper())
        missIDL.append(geneK) # collect matching uniprot accessions too
        procD[geneK] = {geneVI: ["-"]} # store entries as missing for now
      
  print(str(hitCount) + " mouse genes mapped to human names.")
  print(len(procD))
  print("Attempting to map unmapped mouse genes to humans via gene name similarity...")
  
  missNameD = prot_id_converter(missUpperL, "9606", "genesymbol", "uniprotaccession")
  
  removeL = []
  hopelessL = []
  
  countNameNum = 0
  for keyName, valueName in missNameD.items(): # keys are human gene symbols, values are human uniprot accessions. This is only really needed to test if the keys are valid gene names in humans
    if valueName != ["-"]:
      missPos = missUpperL.index(keyName) 
      procD[missIDL[missPos]] = {missL[missPos]: [keyName]}
      removeL.append(missPos)
      # print(keyName, valueName)
      countNameNum += 1
  

  hopeIDL = []
  for j in range(len(missL)):
    if j not in removeL: 
      hopelessL.append(missL[j])
      hopeIDL.append(missIDL[j])

      
  missIDL = hopeIDL

      
  print("\nDone, {} hits found\n".format(countNameNum))
  print("Attempting to map unmapped mouse genes to humans via sequence similarity...")
  
  missSeqD = prot_sequence_finder(hopelessL) # prepare a dict with keys as missing mouse gene symbols and values as their sequences in mice where applicable
  missNameL = []
  missGIL = []
  # print(hopelessL)
  
  relPath = "data/blast_output.txt"  
  blastBackupD = blast_backup_parser(relPath)
   
  for resumeI in blastBackupD.keys(): # find previous version of blast run and use file contents to shorten search
    if resumeI in missSeqD: del missSeqD[resumeI]
    missNameL.append(resumeI)
    missGIL.append(blastBackupD[resumeI])
  
  outF = file_importer(relPath, "w")
  print("BLASTing " + str(len(missSeqD)) + " sequences...")
  for keyS, valueS in missSeqD.items():
    resI = blaster(valueS)
    missGIL.append(resI) # blast sequences and get their human refseq protein GI. this step will take a lot of time. this list will contain human protein genbank accessions
    missNameL.append(keyS) # the matching mouse gene symbols
    outF.write(keyS + "," + resI + "\n")
    
  outF.close()
  print("Creating homologue dict...")
 
    
  missSymbolD = prot_id_converter(missGIL, "9606", "genbankproteinaccession", "genesymbol") # convert protein GIs to gene symbols. keys are human protein GIs, values are human gene symbols
  # print(missSymbolD)

  print("\nprocessing converted names...")
  print("mouse gene names not converted to human:")
  for i in range(len(hopelessL)):
    if hopelessL[i] in missNameL: # human gene names vs human protein IDs. not good.
      # print(hopelessL[i])
      procD[missIDL[i]] = {hopelessL[i]: missSymbolD[missGIL[missNameL.index(hopelessL[i])]]}  # put the newly found gene names into the results dict
      # print([missSymbolD[missGIL[missNameL.index(hopelessL[i])]]])
    else:
      print(hopelessL[i])
    
  mouseGeneL = []
  for keyS in procD.keys():
    for valueL  in procD[keyS].values(): 
      if valueL == ["-"]: continue
      for valueI in valueL:
        if valueI in mouseGeneL: continue
        mouseGeneL.append(valueI)
  
  print("converting human gene symbols to HGNC IDs...")        
  mouseD = prot_id_converter(mouseGeneL, "9606", "genesymbol", "hgncid")
  print("\nDone.")
  

  keyList = gene_list_reader(False)
  
  print("\nwriting results to a file...")
  with file_importer("data/converted_gene_list.csv", "w") as outF:

    for keyItem in keyList:
      outF.write(keyItem + ",") # write out uniprot IDs
      keyItemList = keyItem.split(";")
      currList = []
      for keyItemListItem in keyItemList:
        if keyItemListItem in procD: 
          currList.append(procD[keyItemListItem])
        else: currList.append("-")

      currKL = []
      currVL = []
  
      for currI in currList:
        if currI == "-":
          currKL.append("-")
          currVL.append(["-"])
        
        else:  
          for currK, currV in currI.items():
            currKL.append(currK)
            currVL.append(currV)
        
      for i in range(len(currKL)):# write out mouse gene symbols
        currKI = currKL[i]
        if i == len(currKL) - 1: outF.write(currKI + ",") 
        else: outF.write(currKI + ";")
        
      valL = []    
      nameFlag = False
      for l in range(len(currVL)): # write out human gene symbols
        currVIL = currVL[l]
        if l == len(currVL) - 1: nameFlag = True
        for m in range(len(currVIL)):
          currVI = currVIL[m]

          if currVI in mouseD: valL.append(mouseD[currVI])
          else: valL.append(["-"])
            
          if m == len(currVIL)-1 and nameFlag: outF.write(currVI + ",")
          else: outF.write(currVI + ";")    
      
      endFlag = False
      for k in range(len(valL)): # write out HGNC IDs
        valI = valL[k]
        if k == len(valL) - 1: endFlag = True
        for j in range(len(valI)):
          valII = valI[j]

          if j == len(valI) - 1 and endFlag: outF.write(valII + "\n")

          else: outF.write(valII.rstrip("\n") + ";")            
      

  print("file written")
        
      
species_converter()
  
# blaster("MEFHNGGHVSGIGGFLVSLTSRMKPHTLAVTPALIFAITVATIGSFQFGYNTGVINAPETIIKEFINKTLTDKANAPPSEVLLTNLWSLSVAIFSVGGMIGSFSVGLFVNRFGRRNSMLIVNLLAATGGCLMGLCKIAESVEMLILGRLVIGLFCGLCTGFVPMYIGEISPTALRGAFGTLNQLGIVIGILVAQIFGLELILGSEELWPVLLGFTILPAILQSAALPCCPESPRFLLINRKKEENATRILQRLWGTQDVSQDIQEMKDESARMSQEKQVTVLELFRVSSYRQPIIISIVLQLSQQLSGINAVFYYSTGIFKDAGVQQPIYATISAGVVNTIFTLLSLFLVERAGRRTLHMIGLGGMAFCSTLMTVSLLLKNHYNGMSFVCIGAILVFVACFEIGPGPIPWFIVAELFSQGPRPAAMAVAGCSNWTSNFLVGLLFPSAAYYLGAYVFIIFTGFLITFLAFTFFKVPETRGRTFEDITRAFEGQAHGADRSGKDGVMGMNSIEPAKETTTNV")
  
# print(prot_sequence_finder(["COX7C","SLC2A14","CYP2D6","CYP3A4","LDHAL6A"]))

