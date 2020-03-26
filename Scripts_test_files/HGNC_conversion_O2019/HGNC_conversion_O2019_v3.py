# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:01:38 2019

@author: Admin
"""

def readfile(filetoread, p=None):
    if p is None:
        listoflines = []
    openfile = open(filetoread)
    outputpar = openfile.readlines()
    listoflines = [linein for linein in outputpar if linein.split() != [] and linein.split() != '' ]    
    openfile.close()
    listoflines[-1] = listoflines[-1] + '\n'
    return listoflines

#Writing files: without and with space
def writefilewospace(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w") as f: f.writelines("".join(datatowrite))

def writefile(filenametowrite, datatowrite):
    with open("%s.txt" %filenametowrite, "w") as f: f.writelines("\n".join(datatowrite))

def writefilewithspace(filenametowrite, datatowrite):
    datatowritewithspace = [item +'\n' for item in datatowrite]
    with open("%s.txt" %filenametowrite, "wb+") as f: f.writelines("".join(datatowritewithspace))
    
def HGNCextraction(HGNCfname,threecolfname):
    HGNCfile = readfile('%s.txt'%HGNCfname)
    FiltHGNCfile = [item2 for item2 in HGNCfile if 'Entry Withdrawn' not in item2]
    FiltHGNCfile2 = [item22 for item22 in FiltHGNCfile if 'symbol withdrawn' not in item22]
    firstline = HGNCfile[0].split('\t')
    colind1 = [ind1 for ind1,val1 in enumerate(firstline) if 'symbol' in val1 or 'Synonyms' in val1]
    threecols =[item4.split('\t')[colind1[0]]+'\t'+item4.split('\t')[colind1[1]] +'\t'+item4.split('\t')[colind1[2]] for item4 in FiltHGNCfile2[1:]]
    writefile(threecolfname, threecols)
    return firstline



#PPI conversion
    #HPRD
    
def RefListCreation(fnametoread, Reflist, lreset = None ):
    if lreset is None:
       list1,splitlist,lista,list2part,listaa,listad2,totlist,sortedtotlist = ([] for i in range(8))
    Processed_PPI_file = readfile("%s.txt" %fnametoread)
    splitlist = [item.split() for item in Processed_PPI_file]
# Part 1: If len of item in the HGNC file is 1, create a tuple of the name. It is already an approved symbol, with no alternative names
    list1 = [lst[0]+'\t'+ lst[0] for lst in splitlist if len(lst) ==1]  #P1 of the Ref set
# Part 2: If len of item is 2, create a tuple for approved symbol and alternative names.
    list2part =[item2 for item2 in splitlist if len(item2) >= 2]
    for lst1 in list2part:   # i.e. If multiple old names present
            for gene in lst1: 
                gene2 = gene.split(',')#Create a tuple of approved symbol and every alternate symbol
                lista = lst1[0]+'\t'+ gene2[0]  # Tuple of Approved symbol(lst[0]) +\t + Alias symbol(gene)
                listaa.append(lista)

    totlist = list1 + listaa 
    sortedtotlist = sorted(totlist)
    writefile(Reflist,sortedtotlist)


def RemSingInts(PPIfname,Processedfname, lreset = None):
    if lreset is None:
        intlist,removedInteraction, twocol =([] for i in range(3))
    PPIfileread = readfile('%s.txt'%PPIfname) #Extracted BioGRID PPI file with 2 columns    
    for line11 in PPIfileread:
        inp2 = line11.split()
        if '-' in inp2[0] and len(inp2[0])==1:
            removedInteraction.append(line11)
        elif '-' in inp2[3] and len(inp2[3])==1:
            removedInteraction.append(line11)
            print("In second case")
        else:
            intlist.append(line11)
            intlist2 = inp2[0]+'\t'+ inp2[3]
            twocol.append(intlist2)
    writefile(Processedfname, twocol)


def ConvertPPI(HGNCf, PPI2col,ConvPPIfname, lreset = None):
    if lreset is None:
        ListofInts, Missing = ([] for ij in range(2))
    HGNCfile2 = readfile('%s.txt' %HGNCf)  #HGNC reference file containing tuples
    PPIf = readfile('%s.txt'%PPI2col) #Extracted BioGRID PPI file with 2 columns    
    HGNCL1 = [it3.split()[0] for it3 in HGNCfile2]   #Split the HGNC tuples into two separate lists, for searching convenience
    HGNCL2 = [it4.split()[1] for it4 in HGNCfile2]
    for line in PPIf:
        CorrL1 = ''
        CorrL2 = ''
        linespt = line.split()      # Split to check for each interactor
        if linespt[0] in HGNCL1:    # Quick search to check if gene1 in 1st column
            for ind, val in enumerate(HGNCL1):
                if linespt[0] == val:    #Ensure the name of the gene is an exact match, not a subset
                    CorrL1 = linespt[0]  # If match, then keep the original symbol
            
        elif linespt[0] in HGNCL2:      #If gene in list2[alternate symbol]
            for ind2, val2 in enumerate(HGNCL2):
                if linespt[0] == val2:  #Ensure the name of the gene is an exact match, not a subset
                    CorrL1 = HGNCL1[ind2] #Corr symbol is the symbol in the corresponding list. 
        else:
            pass
    
            
        if linespt[1] in HGNCL1: #Check for the second interactor
            for ind3, val3 in enumerate(HGNCL1):
                if linespt[1] == val3:
                    CorrL2 = linespt[1]
        elif linespt[1] in HGNCL2:
            for ind4, val4 in enumerate(HGNCL2):
                if linespt[1] == val4:
                    CorrL2 = HGNCL1[ind4]
        else:
            pass
    
        if CorrL1 != '' and CorrL2 != '':
            interaction = CorrL1 + '\t' + CorrL2
            ListofInts.append(interaction)


    setofints = set(ListofInts)
    lstofints = list(setofints)
    writefile(ConvPPIfname, lstofints)

def MDGeneListHGNCconv(HGNCReffname,MDGenelist,ConvGLfname,MissGL,lmn = None):
    if lmn is None:
        allnamelist, missinglist = ([] for i in range(2))
    RefFile = readfile('%s.txt'%HGNCReffname)  #List of correct names
    RawFile2= readfile('%s.txt'%MDGenelist) #List to be corrected
    RawFile = sorted(RawFile2)


    for lineinRawFile in RawFile:
        lineinRawFile2 = lineinRawFile.split()
        for lineinRefFile in RefFile:
            lineinRefFile2 = lineinRefFile.split()
            if lineinRawFile2[0] == lineinRefFile2[0]: #If the symbol is already an approved symbol, keep it
                allnamelist.append(lineinRefFile2[0])
            elif lineinRawFile2[0] == lineinRefFile2[1]: #If the symbol in the list matches an old symbol, replace it with the approved one in the new list
                allnamelist.append(lineinRefFile2[0])
    
    nrlist = set(allnamelist)
    namelist = list(nrlist)
    writefile(ConvGLfname,nrlist)
    misslist = [item10 for item10 in RawFile if item10.split()[0] not in namelist]    
    writefile(MissGL,misslist)


HGNCextraction('HGNC_downloaded_Test','HGNC_completeset_processed_Test')
RefListCreation('HGNC_completeset_processed_Test','HGNC_RefList_Test')    
RemSingInts('HPRD_Test','Interactions_SIR_Test')
ConvertPPI('HGNC_RefList_Test2', 'Interactions_SIR_Test','HPRD_HGNC_Test')
MDGeneListHGNCconv('HGNC_RefList_Test2','MD_genes_test','MD_HGNC_Test', 'Missing_TEST')
