import numpy as np
import time
import pickle
import os
import sys
import parmap
import subprocess as sp
import random
from multiprocessing import Manager
import scipy.stats
import math
from Bio import pairwise2
import copy

def BinSaving(JobPath, PythonData, FileName):#c     
    with open(JobPath + FileName+'.Pybin', 'wb') as fb: pickle.dump(PythonData, fb)

def BinLoading(JobPath, FileName):#c      
    with open(JobPath + FileName+'.Pybin', 'rb') as rb: return pickle.load(rb)

def Pre1_DepthDicConstruction(BamPath): # non bed input, non bed output
    os.system("samtools depth -a {0} > {0}.1bpDepth4Bias".format(BamPath))
    with open(BamPath + '.1bpDepth4Bias', 'r') as fr: DL = fr.read().split('\n')
    while DL.count('') > 0: DL.remove('')
    DepthDic = {} # { Contig: {Pos : depth}}
    for i in DL:
        DT = i.split('\t')
        Contig, Pos, Depth = DT[0], int(DT[1]), int(DT[2])
        if Contig in DepthDic: DepthDic[Contig][Pos] = Depth
        else: DepthDic[Contig] = {Pos:Depth}
    return DepthDic


def OneCopyGetFasta(FastaPath, BedPath): # 150kmer extraction fail -> 20kmer + 1cp merge -> 150< region to 150size binning
    print("One Copy Get")
    os.system("awk -F \"\t\" '($3-$2)+20 >= 150 {{print $1\"\t\"$2\"\t\"$3+19}}' {0} > {0}.150cutPlus20.bed".format(BedPath))
    print("Get Fasta from bedtools")
    os.system("bedtools getfasta -fi {0} -bed {1}.150cutPlus20.bed -tab -fo {1}.150cutPlus20_getfasta.bed".format(FastaPath, BedPath))

def OneCopyBed2KmerDic(BedPath, FaiPath): # KmerDic = {Chr : {Pos:[Dat]}}
    BedDat = sp.getoutput("cat {0}".format(BedPath + '.150cutPlus20_getfasta.bed')).split('\n')
    FaiDat = sp.getoutput("cat {0}".format(FaiPath)).split('\n')
    if BedDat.count('') >0: BedDat.remove('')
    if FaiDat.count('') >0: FaiDat.remove('')
    KmerDic = {}
    #### KmerDic Init. ####
    for i in FaiDat:
        DT = i.split('\t')
        Contig, Len = DT[0], int(DT[1])
        
        Index, KmerDic[Contig] = 0, {} # 1lv INIT.
        while Index < Len: 
            KmerDic[Contig][Index] = {'Seq':'Non1cp', 'GC':'Non1cp', 'AT':'Non1cp', 'AT15':'Non1cp'} # 2lv INIT.
            Index += 1
    
    #### 150bp Seq Profiling ####
    Step, DataLen = 0, len(BedDat)
    for i in BedDat:
        if (DataLen - Step) % 1000 == 0: print(DataLen - Step)
        DT = i.split('\t')
        Pos, Seq = DT[0], DT[1].upper() # DT[1] hot fix 20220725
        Contig, Start, End = Pos.split(':')[0], int(Pos.split(':')[1].split('-')[0]), int(Pos.split('-')[1])
        Index = 0
        while Index + 150 <= End-Start:
            SubSeq = Seq[Index:Index + 150]
            GCProp = round( ((SubSeq.count("G") + SubSeq.count("C")) / len(SubSeq)), 3) * 100
            ATProp = round( ((SubSeq.count("A") + SubSeq.count("T")) / len(SubSeq)), 3) * 100
            AT15 = SubSeq.count("ATATATATATATATATATATATATATATAT")
            KmerDic[Contig][Start+Index] = {'Seq':SubSeq, 'GC':GCProp, 'AT':ATProp, 'AT15':AT15}
            Index += 1
        Step += 1
    return KmerDic
    #### 


#def TotalProfiling(KstarTabPath, KmerDic, DepthDic1BP):
#    os.system("awk -F \"\t\" '{{GCPer = (count("G") + count("C")) / 150}} {{ATPer = (count("A") + count("T")) / 150}} {{AT15 = count("ATATATATATATATATATATATATATATAT")}} {{print $5\"\t\"$6+$8\"\t\"GCPer\"\t\"ATPer\"\t\"AT15}}' {0} > {0}_150Mer".format(MerquryTabPath))
#    os.system("paste {0} {1}_150Mer > {0}_BiasProfiling.txt".format(KstarTabPath, MerquryTabPath))



def MeanDepthCal(KstarDat, DepthDic1BP, KmerDic):
    DL = sp.getoutput('cat {0}'.format(KstarDat)).split('\n')
    if DL.count('') > 0: DL.remove('')
    f0 = open(KstarDat + '_TotalProfile_150.txt', 'w')
    Step, DataLen = 0, len(DL)
    for i in DL:
        if (DataLen - Step) % 1000 == 0: print(DataLen - Step)
        Step += 1
        DT = i.split('\t')
        Contig, Start, KStar = DT[0], int(DT[1]), float(DT[4])
        if Contig in KmerDic: SubSeq, GCProp, ATProp, AT15 = KmerDic[Contig][Start]['Seq'], KmerDic[Contig][Start]['GC'], KmerDic[Contig][Start]['AT'], KmerDic[Contig][Start]['AT15']
        else: continue
        if KStar == 0 and SubSeq != 'Non1cp': pass  # Moved from below
        else: continue                              # Moved from below
        #### Mean Depth Cal ####
        Depth, Index = 0, 0
        while Index < 150 - 1:
            if Contig in DepthDic1BP:
                if Start + Index + 1 in DepthDic1BP[Contig]:
                    Depth += DepthDic1BP[Contig][Start + Index + 1]
                    Index += 1
                else: 
                    print(Contig + '\t' + str(Start))
                    Depth = 'ND'
                    continue
            else: 
                print(Contig)
                Depth = 'ND'
                continue
        if Depth == 'ND': MeanDepth = 'ND'
        else: MeanDepth = round(Depth / 150, 1)
        
        #######################
        #### Single Region ####  >>> prototype now. MQ 60 reads set check will be implemented. >>> non-data writing for non-single, erroneous region
        #if KStar == 0 and SubSeq != 'Non1cp': 
        #Type = "NE" # non-error
        Type = "NE"
        f0.write(i + '\t' + '\t'.join(list(map(str, [GCProp, ATProp, AT15, MeanDepth, Type, SubSeq]))) + '\n')
        #else: 
        #    Type = "PE" # potential error
        #    #f0.write(i + '\t' + '\t'.join(list(map(str, [GCProp, ATProp, AT15, MeanDepth, Type, '-']))) + '\n') # Space Save
    DL = 0

def LocusDepthBiasDic(TotalProfilingPath):
    DL = sp.getoutput('cat {0}'.format(TotalProfilingPath)).split('\n')
    if DL.count('') >0: DL.remove('')
    LocusDic, GCDic, ATDic, AT15Dic = {}, {}, {}, {}
    Step, DataLen = 0, len(DL)
    for i in DL:
        if (DataLen - Step) % 1000 == 0: print(DataLen - Step)
        Step += 1
        DT = i.split('\t')
        Chr, Start, GC, AT, AT15, MeanDepth = DT[0], DT[1], round(float(DT[5]),1), round(float(DT[6]),1), DT[7], DT[8]
        Locus = Chr + ':' + str(Start) + '-' + str(int(Start) + 150)
        LocusDic[Locus] = [MeanDepth, GC, AT, AT15]
        if GC in GCDic: GCDic[GC].append(MeanDepth)
        else: GCDic[GC] = [MeanDepth]
        if AT in ATDic: ATDic[AT].append(MeanDepth)
        else: ATDic[AT] = [MeanDepth]
        if AT15 in AT15Dic: AT15Dic[AT15].append(MeanDepth)
        else: AT15Dic[AT15] = [MeanDepth]
    DL = 0
    return LocusDic, GCDic, ATDic, AT15Dic


def Dic2Txt(BiasDic, TxtPath):
    f0 = open(TxtPath, 'w')
    BiasTup = sorted(BiasDic.items())
    for i in BiasTup:
        GC, GCDat = i[0], list(map(float, i[1]))
        Mean, Var = round(sum(GCDat) / len(GCDat), 1), round(np.std(GCDat), 2)
        f0.write(str(GC) + '\t'+ str(Mean) + '\t'+ str(Var) + '\t' + str(len(GCDat)) + '\n')
    f0.close()


def Dic2IntDic(BiasDic):
    NewBiasDic = {}
    for i in BiasDic:
        NewKey, Dat = round(i), BiasDic[i]
        if NewKey in NewBiasDic: NewBiasDic[NewKey] += Dat
        else: NewBiasDic[NewKey] = Dat
    return NewBiasDic


if __name__=="__main__": 
    ############ OPERATION #############
    OneCPBedPath, FaiPath, KstarDatPath = 'bTaeGut2_CCS.ass.1cp.bed', 'bTaeGut2_CCS_Diploid.fasta.fai', 'bTaeGut_CCS_merfin.dump'
    BamPath, FastaPath = 'bTaeGut2_CCS_10X.bam', 'bTaeGut2_CCS.fasta.masked'
    DepthDic = Pre1_DepthDicConstruction(BamPath) # previous 1bp bam file of above was deleted.
    BinSaving('./', DepthDic, '1BPDepthDic4Bias')
    OneCopyGetFasta(FastaPath, OneCPBedPath)
    KmerDic = OneCopyBed2KmerDic(OneCPBedPath, FaiPath)
    BinSaving('./', KmerDic, 'KmerMetaDic')

    DepthDic = BinLoading('./', '1BPDepthDic4Bias')
    KmerDic = BinLoading('./', 'KmerMetaDic')
    MeanDepthCal(KstarDatPath, DepthDic, KmerDic)

    ResultDic = LocusDepthBiasDic('bTaeGut_CCS_merfin.dump_TotalProfile_150.txt')
    BinSaving('./', ResultDic[0], 'LocusDic')
    BinSaving('./', ResultDic[1], 'GCDic')
    BinSaving('./', ResultDic[2], 'ATDic')
    BinSaving('./', ResultDic[3], 'AT15Dic')
    GCDic = BinLoading('./', "GCDic")
    GCDic = Dic2IntDic(GCDic)
    Dic2Txt(GCDic, "GCBiasHistogram.txt")

