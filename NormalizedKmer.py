import numpy as np
import time
import pickle
import os
import sys
import parmap
import subprocess as sp
from numpy import random
from multiprocessing import Manager
import scipy.stats
import math
#from Bio import pairwise2
import copy
import gzip
from Bio.Seq import Seq as BioSeq

## genome must have expansion 'fasta'


def BinSaving(JobPath, PythonData, FileName):#c 
    with open(JobPath + FileName+'.Pybin', 'wb') as fb: pickle.dump(PythonData, fb)

def BinLoading(JobPath, FileName):#c  
    with open(JobPath + FileName+'.Pybin', 'rb') as rb: return pickle.load(rb)

def Base2Number(Base):
    Transfer = {'A':'1', 'C':'2', 'G':'3', 'T':'4'}
    Number = ''
    for i in Base: Number += Transfer[i]
    return int(Number)

def LoadGCBiasHist(GCHistPath):
    Dat = sp.getoutput("cat {0}".format(GCHistPath)).split('\n')
    if Dat.count('') > 0: Dat.remove('')
    GCDic = {} # GC:weight (weight = MeanDepth / Mode)
    #ModeDepth, ModeFreq = 0, 0
    # Find Mode GC
    MaxDepth = 0
    for i in Dat:
        DT = i.split('\t')
        GC, Depth, Var, Freq = int(DT[0]), float(DT[1]), float(DT[2]), int(DT[3])
        if MaxDepth < Depth: MaxDepth = Depth
        #if ModeFreq < Freq: ModeDepth, ModeFreq = Depth, Freq
    for i in Dat:
        DT = i.split('\t')
        GC, Depth, Var, Freq = int(DT[0]), float(DT[1]), float(DT[2]), int(DT[3])
        GCDic[GC] = round(1 / (Depth / MaxDepth), 3) #ModeDepth), 3)
    UpperValue = GCDic[max(list(GCDic.keys()))]
    LowerValue = GCDic[min(list(GCDic.keys()))]
    for i in list(range(max(list(GCDic.keys())) + 1, 100 + 1)):
        GCDic[i] = UpperValue
    for i in list(range(0, min(list(GCDic.keys())))):
        GCDic[i] = LowerValue
    return GCDic


def Pre1_DepthDicConstruction(BamPath): # non bed input, non bed output
    os.system("samtools depth -aa {0} > {0}.1bpDepth4Bias".format(BamPath))
    with open(BamPath + '.1bpDepth4Bias', 'r') as f0: 
        Step = 0
        ContigDic, PContig = {}, 'ND' #{Pos:Depth}
        os.system("mkdir DepthDic")
        while True:
            Step += 1
            if Step % 100000 == 0: print(Step)
            DL = f0.readline().split('\n')[0]
            if not DL: break
            if DL == '': break
            DT = DL.split('\t')
            Contig, Pos, Depth = DT[0], int(DT[1]), int(DT[2])
            if Contig != PContig: 
                BinSaving('./DepthDic/', ContigDic, PContig)
                ContigDic = {Pos : Depth} # Second Init.
            else: 
                ContigDic[Pos] = Depth
            PContig = Contig
        BinSaving('./DepthDic/', ContigDic, PContig)


def Pre2_MerquryMultiplicities(AsmPath1, AsmPath2, CuttedMergedRead1, MergedRead2, KLength):
    os.system("meryl k={0} count output Merged_R1.meryl {1}".format(KLength, CuttedMergedRead1))
    os.system("gzip -d {0}.gz".format(MergedRead2))
    os.system("meryl k={0} count output Merged_R2.meryl {1}".format(KLength, MergedRead2))
    os.system("meryl union-sum output 10XRead.meryl Merged_R*.meryl")
    ##os.system("rm -r Merged_R1.meryl")
    ##os.system("rm -r Merged_R2.meryl")
    os.system("cat {0} | awk -F \"\t\" 'index($1,\">\") == 1 {{print $1\"_alt\"}} index($1, \">\") == 0 {{print $1}}' > {0}_header".format(AsmPath2))
    os.system("cat {0} {1} > Diploid.fasta".format(AsmPath1, AsmPath2 + '_header'))
    AsmPath = 'Diploid.fasta'
    os.system("mkdir merqury")
    ###### "export MERQURY=$PWD" ########
    ###### "export PATH=<Path>/meryl-1.0/Linux-amd64/bin:$PATH" #######
    os.chdir("./merqury/")
    os.system("ln -s ../10XRead.meryl")
    os.system("~/Program/merqury/merqury.sh 10XRead.meryl ../{0} 10Xmerqury > 10Xmerqury.log".format(AsmPath))
    os.system("meryl-lookup -dump -memory 100 -sequence ../{0} -mers {1}.meryl > 10Xmulti.ass.bed".format(AsmPath, AsmPath.split('.fasta')[0]))
    os.system("meryl-lookup -dump -memory 100 -sequence ../{0} -mers 10XRead.meryl > 10Xmulti.read.bed".format(AsmPath))
    os.chdir("../")
    os.system("paste ./merqury/10Xmulti.read.bed ./merqury/10Xmulti.ass.bed > ./merqury/10Xmulti.merged.bed") # merge
    os.system("rm ./merqury/10Xmulti.ass.bed")
    os.system("rm ./merqury/10Xmulti.read.bed")
    os.system("mv ./merqury/10Xmulti.merged.bed ./")


def P3_DoGenomeScope(ReadHistPath, KLen, Name4Write):
    os.system("~/Program/genomescope2.0/genomescope.R -i {0} -o ./{2}/ -k {1} --fitted_hist".format(ReadHistPath, KLen, Name4Write))
    HapPeak = float(sp.getoutput("cat ./{0}/model.txt | awk '$1 == \"kmercov\" {{print $2}}'".format(Name4Write)))
    Threshold = int(sp.getoutput("cat ./{0}/lookup_table.txt | wc -l".format(Name4Write))) + 1
    LookupTablePath = './{0}/lookup_table.txt'.format(Name4Write)
    return HapPeak, Threshold, LookupTablePath

def P3_LookupTable2Dic(LookupTablePath):
    CNDic = {0:0}
    DL = sp.getoutput('cat {0}'.format(LookupTablePath)).split('\n')
    if DL.count('') > 0: DL.remove('')
    Multi = 1
    for i in DL:
        DT = i.split(',')
        CN = int(DT[0])
        CNDic[Multi] = CN
        Multi += 1
    return CNDic

def P3_KStarCal(RCN, ACN, LookupDic, HapPeak, Cutoff): # Cutoff = 178 in bTaeGut2
    if RCN < Cutoff: ExpCN = LookupDic[RCN]
    else: ExpCN = int(round(RCN / HapPeak, 0)) # need to check using hap peak -> OK
    if ACN == 0 and ExpCN == 0: KStar = 'false sequence in Ref'
    elif ACN == 0: KStar = 'inf'
    elif ExpCN == 0: KStar = '-inf'
    else: KStar = round((ExpCN - ACN) / min(ExpCN, ACN), 2)
    return KStar, ExpCN

def Pre3_KStarCal(LookupTablePath, HapPeak, Cutoff):
    f0 = open('10Xmulti.merged.bed', 'r')
    f1 = open('10Xmulti.merged.bed' + ".KStar", 'w')
    LTDic = P3_LookupTable2Dic(LookupTablePath)
    Step = 0
    while True:
        Step += 1
        if Step % 100000 == 0: print(Step)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split('\t')
        Chr, Start, Kmer1, Kmer2, RCN, ACN = DT[0], DT[2], DT[4], DT[6], int(DT[5]) + int(DT[7]), int(DT[13]) + int(DT[15])
        IKmer1, IKmer2 = Base2Number(Kmer1), Base2Number(Kmer2)
        IKmer = min(IKmer1, IKmer2)
        KStarsCal = P3_KStarCal(RCN, ACN, LTDic, HapPeak, Cutoff)
        KStar, ExpACN = KStarsCal[0], KStarsCal[1]
        f1.write("\t".join([Chr, Start, str(int(Start) + 1), str(IKmer), str(ACN), str(RCN), str(ExpACN), str(KStar)]) + '\n') # > will be used for further part

def Pre4_Pri1Cp(AsmPath, KLength):
    os.chdir("./merqury/")
    os.system("meryl k={0} count ../{1} output {2}.meryl".format(KLength, AsmPath, AsmPath.split('.fasta')[0]))
    os.system("meryl-lookup -dump -memory 100 -sequence ../{0} -mers {1}.meryl > Pri.ass.bed".format(AsmPath, AsmPath.split('.fasta')[0]))
    os.system("awk -F \"\t\" '{{print $6+$8}}' Pri.ass.bed > Pri.ass.bed.cp")
    os.chdir("../")
    os.system("mv ./merqury/Pri.ass.bed.cp ./")
    
    #[hotfix 20221219]
    #os.system("awk -F \"\t\" '$6 + $8 == 1 {{print $1\"\t\"$3\"\t\"$3+1\"\t\"$6+$8}}' Pri.ass.bed > Pri.ass.bed.1cp")
    #os.chdir("../")
    #os.system("mv ./merqury/Pri.ass.bed.1cp ./")


def Main1_Get1CP0KStar(FlankingSize, ZeroPropCutoff):
    SizeLimit = FlankingSize * 2 + 150
    #Get Pri 1cp
    ##[hotfix] os.system("paste 10Xmulti.merged.bed.KStar Pri.ass.bed.1cp > 10Xmulti.merged.bed.KStar.Pri1cp") # needto intersect 근데시간너무많이걸려서 paste로바꾼거. 따라서 위 Pre4  수정 필요.
    os.system("paste 10Xmulti.merged.bed.KStar Pri.ass.bed.cp > 10Xmulti.merged.bed.KStar.Pricp") # Pri1cp inf is just merged
    os.system("awk -F \"\t\" '$5 == 2 && $9 == 1 {{print $1\"\t\"$2\"\t\"$3\"\t\"$8}}' 10Xmulti.merged.bed.KStar.Pricp > 10Xmulti.merged.bed.KStar.Pri1cp.Di2cp") # $12 -> $9
    #1cp merge
    os.system("~/Program/bedtools2/bin/bedtools merge -i 10Xmulti.merged.bed.KStar.Pri1cp.Di2cp -d 1 -c 4 -o collapse -delim \"|\" > 10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge") #
    os.system("awk -F \"\t\" '$3-$2 >= {0} {{print $1\"\t\"$2 + {1}\"\t\"$3 - {1}\"\t\"$4}}'  10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge >  10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered ".format(SizeLimit, FlankingSize))
    f0 = open('10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered', 'r')
    f1 = open('10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered.KStarFiltered', 'w')
    Step = 0
    while True:
        Step += 1
        if Step % 100000 == 0: print(Step)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split('\t')
        Chr, Start, End, KStars = DT[0], DT[1], DT[2], list(map(float, DT[3].split('|')))
        #print(KStars)
        if Chr.find('_alt') > -1: continue # alt loci
        ErrProp =  len(list(filter(lambda x: x>0, KStars))) + KStars.count(float('-inf'))  / len(KStars) # GC bais can derive negative KStar by lower ExpRCN
        if ErrProp < ZeroPropCutoff: f1.write(DL + '\t' + str(round(ErrProp, 2)) +'\n') # merged bed (from 1bp) after filtering with least size and zero proportion of Kmer

    f0.close()
    f1.close()


def M2_DepthDicLoad(DepthDic1BP, Start, Contig):
    Depth, Index = 0, 0
    while Index < 150 - 1:
        if Start + Index + 1 in DepthDic1BP:
            Depth += DepthDic1BP[Start + Index + 1]
            Index += 1
        else:
            print(Contig + '\t' + str(Start))
            Depth = 'ND'
            continue
    if Depth == 'ND': MeanDepth = 'ND'
    else: MeanDepth = round(Depth / 150, 1)
    return MeanDepth


def Main2_GetGCPropAndDepth(FastaPath, KSize): # 150kmer extraction fail -> 20kmer + 1cp merge -> 150< region to 150size binning
    print("One Copy Get")
    os.system("awk -F \"\t\" '{{print $1\"\t\"$2\"\t\"$3+{0}}}' 10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered.KStarFiltered > 10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered.KStarFiltered.PlusK.bed".format(KSize-1)) # for expansion to terminal kmer region
    print("Get Fasta from bedtools") # for 150Kmer construction
    os.system("bedtools getfasta -fi {0} -bed 10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered.KStarFiltered.PlusK.bed -tab -fo 10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered.KStarFiltered.PlusK_getfasta.bed".format(FastaPath))
    BedDat = sp.getoutput("cat {0}".format('10Xmulti.merged.bed.KStar.Pri1cp.Di2cp.bedmerge.SizeFiltered.KStarFiltered.PlusK_getfasta.bed')).split('\n')
    if BedDat.count('') >0: BedDat.remove('')
    
    f0 = open('FinalFiltered_150Profile.bed', 'w')
    Step, DataLen, LoadedContig = 0, len(BedDat), "INIT"
    for i in BedDat:
        if (DataLen - Step) % 1000 == 0: print(DataLen - Step)
        DT = i.split('\t')
        Pos, Seq = DT[0], DT[1].upper() # DT[1] hot fix 20220725
        Contig, Start, End, = Pos.split(':')[0], int(Pos.split(':')[1].split('-')[0]), int(Pos.split('-')[1])
        if Contig == LoadedContig: pass
        else: 
            DepthDic = BinLoading('./DepthDic/', Contig)
            print("DepthDic Loading: {0}".format(Contig))
        Index = 0
        while Index + 150 <= End-Start: # 150bp window calculation of GC and GA
            SubSeq = Seq[Index:Index + 150]
            GCProp = round( ((SubSeq.count("G") + SubSeq.count("C")) / 150), 3) * 100
            GAProp = round( ((SubSeq.count("G") + SubSeq.count("A")) / 150), 3) * 100
            AT15 = SubSeq.count("ATATATATATATATATATATATATATATAT")
            MeanDepth = M2_DepthDicLoad(DepthDic, Start + Index, Contig)
            f0.write('\t'.join(list(map(str, [Contig, Start + Index, SubSeq, GCProp, GAProp, AT15, MeanDepth]))) + '\n')

            Index += 1

        LoadedContig = Contig
        Step += 1

def Main3_LocusDepthBiasDic():
    DL = sp.getoutput('cat FinalFiltered_150Profile.bed').split('\n')
    if DL.count('') >0: DL.remove('')
    LocusDic, GCDic, GADic, AT15Dic = {}, {}, {}, {}
    Step, DataLen = 0, len(DL)
    for i in DL:
        if (DataLen - Step) % 1000 == 0: print(DataLen - Step)
        Step += 1
        DT = i.split('\t')
        Chr, Start, GC, GA, AT15, MeanDepth = DT[0], DT[1], round(float(DT[3]),1), round(float(DT[4]),1), DT[5], DT[6]
        Locus = Chr + ':' + str(Start) + '-' + str(int(Start) + 150)
        LocusDic[Locus] = [MeanDepth, GC, GA, AT15]
        if GC in GCDic: GCDic[GC].append(MeanDepth)
        else: GCDic[GC] = [MeanDepth]
        if GA in GADic: GADic[GA].append(MeanDepth)
        else: GADic[GA] = [MeanDepth]
        if AT15 in AT15Dic: AT15Dic[AT15].append(MeanDepth)
        else: AT15Dic[AT15] = [MeanDepth]
    DL = 0
    return LocusDic, GCDic, GADic, AT15Dic


def Main4_Dic2Txt(BiasDic, TxtPath):
    BiasDic = M4_Dic2IntDic(BiasDic)
    f0 = open(TxtPath, 'w')
    BiasTup = sorted(BiasDic.items())
    for i in BiasTup:
        GC, GCDat = i[0], list(map(float, i[1]))
        Mean, Var = round(sum(GCDat) / len(GCDat), 1), round(np.std(GCDat), 2)
        f0.write(str(GC) + '\t'+ str(Mean) + '\t'+ str(Var) + '\t' + str(len(GCDat)) + '\n')
    f0.close()

def M4_Dic2IntDic(BiasDic):
    NewBiasDic = {}
    for i in BiasDic:
        NewKey, Dat = round(i), BiasDic[i]
        if NewKey in NewBiasDic: NewBiasDic[NewKey] += Dat
        else: NewBiasDic[NewKey] = Dat
    return NewBiasDic


def M5_Read2Kmer4GC(Seq, KSize, fw): # {Kmer:{'Ct':Nubmer of Count, 'GC':GC percent List ~ [30, 25, 30..], 'AT':AT percent, 'AT15':[1,2,1,3,1..]}, Kmer2:count2...}
    n_Kmers = len(Seq) - KSize + 1
    GCProp = round( float(str(((Seq.count("G") + Seq.count("C")) / len(Seq)))[0:5]), 3) * 100   # GC prop is calculted for the whole sequence imported 
    GAProp = round( float(str(((Seq.count("G") + Seq.count("A")) / len(Seq)))[0:5]), 3) * 100
    AT15 = Seq.count("ATATATATATATATATATATATATATATAT")
    GCProp, GAProp = int(round(GCProp, 0)), int(round(GAProp, 0))
    for i in range(n_Kmers):
        Kmer = Seq[i:i + KSize]
        Kmer1, Kmer2 = Kmer, str(BioSeq(Kmer).reverse_complement())    
        IKmer1, IKmer2 = Base2Number(Kmer1), Base2Number(Kmer2)
        IKmer = min(IKmer1, IKmer2)
        fw.write(str(IKmer) + "\t" + "\t".join(list(map(str, [GCProp, GAProp]))) + '\n')
    

def M5_FastqSpliter(Fastq, DirPath, Cores):
    HalfCores = divmod(Cores, 2)[0]
    print("HalfCores: {0}".format(HalfCores))
    DatLen = int(sp.getoutput("cat {0} | wc -l".format(Fastq)))
    Share = divmod(DatLen, HalfCores)[0] + divmod(DatLen, HalfCores)[0] %4
    print(str(DatLen) + '\t' + str(Share) + '\t' + str(DatLen/Share) + '\t' + str((DataLen / (DatLen/Share)) / 4))
    os.system("split -l {0} {1} {1}_fastqsplit".format(Share, Fastq))
    os.system("mv *_fastqsplit* {0}".format(DirPath))


def Main5_Read2Kmer(MergedAdapterRemoved10XFastq, Merged10XFastq, KSize, Cores):
    os.system("mkdir readsplit")
    M5_FastqSpliter(MergedAdapterRemoved10XFastq, './readsplit', Cores)
    M5_FastqSpliter(Merged10XFastq, './readsplit', Cores)
    file_list = os.listdir('./readsplit')
    print(file_list) #non-np.array_split. FastqSpliter will be split he fastq as much as the number of Cores u input.
    parmap.map(M_M5_Read2Kmer, file_list, './readsplit/', KSize, pm_pbar = True)


def M_M5_Read2Kmer(File, SuffixPath, KSize):
    f0 = open(SuffixPath + File, 'r')
    fw = open(SuffixPath + File + ".rawcount", 'w')
    Line = 4
    NumRead, NSeq = 0, 0
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        if Line % 4 == 2: pass
        else: continue
        if DL.find("N") > -1: # need to check 0 or -1 
            NumRead += 1
            NSeq += 1
            continue
        else: NumRead += 1
        Seq = DL
        M5_Read2Kmer4GC(Seq, KSize, fw)
    print(NumRead)
    print(NSeq)

######################################
############# SUPPRESSED #############
######################################

def Main6_KmerCounting_IndividualMerge(BiasDic): # [Notused now]
    file_list = os.listdir("./readsplit/")
    for i in file_list:
        if i.find("rawcount") > -1:
            if i.find("sorted") > -1: continue
            print(i)
            Prefix = i
            M6_KmerCounting_Multiproc(BiasDic, Prefix)

def Main6_KmerCounting_PrefixMerge(BiasDic, Prefix): # 1. merge the files by sorting. 2. line by line merge for 1count kmer
    #os.system("cat ./readsplit/{0}*.rawcount > {0}.intermerged.R1R2".format(Prefix)) [old]
    #os.system("sort -S 90% -T ./ --parallel=20 -nk1 {0}.intermerged.R1R2 > {0}.intermerged.R1R2.sorted".format(Prefix)) [old]
    os.system("sort -S 85% -T ./ --parallel=20 -nk1 ./readsplit/{0}*.rawcount > {0}.intermerged.R1R2.sorted".format(Prefix))
    f0 = open('{0}.intermerged.R1R2.sorted'.format(Prefix), 'r')
    f1 = open('{0}.intermerged.R1R2.sorted.merged'.format(Prefix), 'w')
    Line, PKmer, Temp = 0, "INIT", ['', 0, [], [], 0, 0] # kmer, Count, GCWeight, GAWeight,  mean GCWeight, mean GAWeight
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        #IKmer, Count, GCProp, GAProp, GCWeight, GAWeight = int(DT[0]), int(DT[1]), int(DT[2]), int(DT[3]), float(DT[4]), float(DT[5]) # Prop will be loaded as int for BiasDic loading
        IKmer, GCProp, GAProp = int(DT[0]), int(DT[1]), int(DT[2])
        GCWeight, GAWeight = round(BiasDic[GCProp], 2), round(BiasDic[GAProp], 2)
        if IKmer == PKmer:
            Temp[1] += 1
            Temp[2].append(GCProp)
            Temp[3].append(GAProp)
            Temp[4] += GCWeight
            Temp[5] += GAWeight
        else:
            if PKmer == "INIT": pass
            else:
                Temp[2], Temp[3] = round(sum(Temp[2]) / Temp[1], 2), round(sum(Temp[3]) / Temp[1], 2)
                f1.write("\t".join(list(map(str, Temp))) + '\n') # Writing
            Temp = [IKmer, 1, [GCProp], [GAProp], GCWeight, GAWeight] # Init           
        PKmer = IKmer
    Temp[2], Temp[3] = round(sum(Temp[2]) / Temp[1], 2), round(sum(Temp[3]) / Temp[1], 2)
    f1.write("\t".join(list(map(str, Temp))) + '\n') # Final writing

###############################################
###############################################
###############################################

def Main6_KmerCounting_PrefixMerge_MultiProc(BiasDic, Prefix): # 1. merge the files by sorting. 2. line by line merge for 1count kmer
    os.system("sort -S 50% -T ./ --parallel=20 -nk1 ./readsplit/{0}*.rawcount > {0}.intermerged.R1R2.sorted".format(Prefix))
    Prefix2 = "{0}.intermerged.R1R2".format(Prefix)
    M6_KmerCounting_Multiproc(BiasDic, Prefix2)    


def M6_3lvSpliter(TargetPath, Prefix):
    f0 = open(TargetPath, 'r')
    os.system("mkdir ./{0}_kmersplit/".format(Prefix))
    f1 = open('./{0}_kmersplit/111.rawcount.split'.format(Prefix), 'w') 
    PIndex, Line = 'ND', 0
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        IKmer = DT[0]
        Index = IKmer[0:3]
        if Index == PIndex: f1.write(DL + '\n')
        else:
            f1.close()
            f1 = open('./{0}_kmersplit/{1}.rawcount.split'.format(Prefix, Index), 'w')
            f1.write(DL + '\n')
            PIndex = Index
    f1.close()


def M6_KmerCounting_Multiproc(BiasDic, Prefix): # Num cores used in multiproc is fixed as 16 by 2lv indexing ( 4*4 )
    M6_3lvSpliter('./{0}.sorted'.format(Prefix), Prefix)
    file_list = os.listdir('./{0}_kmersplit/'.format(Prefix)) # Number of file should be on the cores is fixed as 16
    print(file_list)
    parmap.map(M_M6_KmerCounting, file_list, BiasDic, Prefix, pm_pbar = True)
    os.system("cat ./{0}_kmersplit/*.KmerMerged > {0}.sorted.merged".format(Prefix))


def M_M6_KmerCounting(FileName, BiasDic, Prefix): # sorted splited raw count by suffix
    print(FileName)
    f0 = open('./{0}_kmersplit/'.format(Prefix) + FileName, 'r')
    f1 = open('./{0}_kmersplit/'.format(Prefix) + FileName + '.KmerMerged', 'w')
    Line, PKmer, Temp, TempStr, TempLen = 0, "INIT", ['', 0, [], [], 0, 0], '', 0 # kmer, Count, GCWeight, GAWeight,  mean GCWeight, mean GAWeight
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        IKmer, GCProp, GAProp = int(DT[0]), int(DT[1]), int(DT[2])
        GCWeight, GAWeight = round(BiasDic[GCProp], 2), round(BiasDic[GAProp], 2)
        if IKmer == PKmer:
            Temp[1] += 1
            Temp[2].append(GCProp)
            Temp[3].append(GAProp)
            Temp[4] += GCWeight
            Temp[5] += GAWeight
        else:
            if PKmer == "INIT": pass
            else:
                Temp[2], Temp[3] = round(sum(Temp[2]) / Temp[1], 2), round(sum(Temp[3]) / Temp[1], 2)
                TempStr += "\t".join(list(map(str, Temp))) + '\n'
                TempLen += 1
                if TempLen == 100000:
                    f1.write(TempStr) # Writing
                    TempStr, TempLen = '', 0 # INIT TempData

            Temp = [IKmer, 1, [GCProp], [GAProp], GCWeight, GAWeight] # Init
        PKmer = IKmer
    Temp[2], Temp[3] = round(sum(Temp[2]) / Temp[1], 2), round(sum(Temp[3]) / Temp[1], 2)
    TempStr += "\t".join(list(map(str, Temp))) + '\n' # Final Add to temp
    f1.write(TempStr) # Final writing


def M6_KmerCounting_Multiproc_Final(BiasDic, Prefix): # Num cores used in multiproc is fixed as 16 by 2lv indexing ( 4*4 )
    M6_3lvSpliter('./{0}.sorted'.format(Prefix), Prefix)
    file_list = os.listdir('./{0}_kmersplit/'.format(Prefix)) # Number of file should be on the cores is fixed as 16
    file_list.sort()
    print(file_list)
    parmap.map(M_M6_KmerCounting_Final, file_list, BiasDic, Prefix, pm_pbar = True)
    os.system("cat ./{0}_kmersplit/*.KmerMerged > {0}.sorted.merged".format(Prefix))


def M_M6_KmerCounting_Final(FileName, BiasDic, Prefix): # sorted splited raw count by suffix
    print(FileName)
    f0 = open('./{0}_kmersplit/'.format(Prefix) + FileName, 'r')
    f1 = open('./{0}_kmersplit/'.format(Prefix) + FileName + '.KmerMerged', 'w')
    Line, PKmer, Temp, TempStr, TempLen = 0, "INIT", ['', 0, 0, 0, 0, 0], '', 0 # kmer, Count, GCWeight, GAWeight,  mean GCWeight, mean GAWeight
    while True:
        Line += 1
        #if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        #IKmer, Count, GCProp, GAProp, GCWeight, GAWeight = int(DT[0]), int(DT[1]), float(DT[2]), float(DT[3]), int(round(float(DT[4]), 0)), int(round(float(DT[5]), 0)) # Prop will be loaded as int for BiasDic loading
        IKmer, Count, GCProp, GAProp, GCWeight, GAWeight = int(DT[0]), int(DT[1]), float(DT[2]), 0, int(round(float(DT[4]), 0)), 0
        if IKmer == PKmer:
            Temp[1] += Count
            Temp[2] = round(((Temp[2] * (Temp[1]-Count)) + GCProp * Count) / Temp[1], 2) # Mean GC Proportion Calculation by update by line
            #Temp[3] = round(((Temp[3] * (Temp[1]-Count)) + GAProp * Count) / Temp[1], 2)
            Temp[4] += GCWeight
            #Temp[5] += GAWeight
        else:
            if PKmer == "INIT": pass
            else:
                TempStr += "\t".join(list(map(str, Temp))) + '\n'
                TempLen += 1
                if TempLen == 1000000:
                    f1.write(TempStr) # Writing
                    TempStr, TempLen = '', 0 # INIT TempData

            Temp = [IKmer, Count, GCProp, GAProp, GCWeight, GAWeight] # Init
        PKmer = IKmer
    TempStr += "\t".join(list(map(str, Temp))) + '\n' # Final Add to temp
    f1.write(TempStr) # Final writing


def Main6_KmerCounting_FinalMerge(BiasDic): # 1. sort intermerged 2 files, final merge by merged counts
    os.system("sort -S 80% -T ./ --parallel=20 -nk1 *.intermerged.R1R2.sorted.merged > Finalintermerged.sorted")
    Prefix = 'Finalintermerged'
    M6_KmerCounting_Multiproc_Final(BiasDic, Prefix)
    '''
    f0 = open('Finalintermerged.sorted', 'r')
    f1 = open('Finalintermerged.sorted.merged', 'w')
    Line, PKmer, Temp = 0, "INIT", ['', 0, 0, 0, 0, 0] # kmer, Count, GCWeight, GAWeight,  mean GCWeight, mean GAWeight
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        IKmer, Count, GCProp, GAProp, GCWeight, GAWeight = int(DT[0]), int(DT[1]), float(DT[2]), float(DT[3]), float(DT[4]), float(DT[5]) # Prop will be loaded as int for BiasDic loading
        if IKmer == PKmer:
            Temp[1] += Count
            Temp[2] = round(((Temp[2] * Temp[1]-Count) + GCProp * Count) / Temp[1], 2) # Mean GC Proportion Calculation by update by line
            Temp[3] = round(((Temp[3] * Temp[1]-Count) + GAProp * Count) / Temp[1], 2)
            Temp[4] += GCWeight
            Temp[5] += GAWeight
        else:
            if PKmer == "INIT": pass
            else: 
                f1.write("\t".join(list(map(str, Temp))) + '\n') # Writing
            Temp = [IKmer, Count, GCProp, GAProp, GCWeight, GAWeight] # Init   
        PKmer = IKmer
    f1.write("\t".join(list(map(str, Temp))) + '\n') # Final writin
    '''


def Main7_NKmerSync(AsmBedPath, FinalKmerPath):
    os.system("sort -S 75% -T ./ --parallel=20 -nk4 {0} > {0}.sorted".format(AsmBedPath))
    f0 = open(AsmBedPath + '.sorted' , 'r')
    f1 = open(AsmBedPath + '.sorted.NKmer', 'w')
    f2 = open(FinalKmerPath, 'r')
    Line, CKmer, CKmerDat = 0, 0, "ND"
    print(f0.readline())
    while True:
        DL = f0.readline().split('\n')[0]
        IKmer = int(DL.split('\t')[3])
        if not DL: break
        while CKmer < IKmer: 
            Line += 1
            if Line % 100000 == 0: print(Line)
            KDL = f2.readline()
            CKmer = int(KDL[0:20])
        if CKmer == IKmer: f1.write(DL + '\t' + KDL[20:-2] + '\n')
        else: f1.write(DL + "\tND in KmerDB" + '\n')


############################## SUPPRESSED #############################
#######################################################################

def M_M7_Count2DicGC(KmerCountPath):
    f0 = open(KmerCountPath, 'r')
    Line = 0
    KmerDic = {} # INIT, Kmer:(count, GCProp, GAProp, GCWeight, GAWeight)
    while True:
        Line += 1
        if Line % 100000 == 0: print(KmerCountPath +': ' +str(Line))
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        KmerDic[int(DT[0])] = (round(float(DT[2]), 2), int(round(float(DT[4]), 0)))
        #IKmer, RC, GCProp, GCWeight = int(DT[0]), int(DT[1]), round(float(DT[2]), 2), int(round(float(DT[4]), 0))
        #KmerDic[IKmer] = (GCProp, GCWeight) # Add Data
    BinSaving('./IKmerDic/', KmerDic, KmerCountPath)


def Main7_Count2DicGC(KmerCountPath):
    f0 = open(KmerCountPath, 'r')
    Line, FileIndex = 0, 11
    KmerDic = {} # INIT, Kmer:(count, GCProp, GAProp, GCWeight, GAWeight)
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        if Line> 1: continue
        #IKmer, RC, GCProp, GCWeight = int(DT[0]), int(DT[1]), round(float(DT[2]), 2), int(round(float(DT[4]), 0))
        #Dat = (RC, GCProp, GCWeight)
        #Prefixer = int(str(IKmer)[0:2])
        Prefixer = int(str(int(DT[0]))[0:2])
        if Prefixer == FileIndex:
            KmerDic[int(DT[0])] = (round(float(DT[2]), 2), int(round(float(DT[4]), 0)))
            #KmerDic[IKmer] = Dat # Add Data
        else:
            BinSaving('./IKmerDic/', KmerDic, str(FileIndex)) # Save Data
            FileIndex = Prefixer #  FileIndex Change
            KmerDic = {int(DT[0]) : (round(float(DT[2]), 2), int(round(float(DT[4]), 0)))}
            #KmerDic = {IKmer:Dat} # INIT continuous
    BinSaving('./IKmerDic/', KmerDic, str(FileIndex)) # Final Save

##############################################################################################
##############################################################################################


def M8_GetNormalizedLookupTable(WeightCountPosition, KLength): # GC = 5, GA = 6 (awk position)
    os.system("awk -F \"\t\" 'BEGIN{{ arr[0] = 0 }} {{ $5 = int($5 + 0.5) }} {{ if($5 in arr){{ arr[$5] += 1 }} else{{ arr[$5] = 1 }} }} END{{ for (i in arr){{ print i\"\t\"arr[i]}} }}' Finalintermerged.sorted.merged | sort -n -k1 > Finalintermerged.sorted.merged.hist".format(WeightCountPosition)
    GS = P3_DoGenomeScope("Finalintermerged.sorted.merged.hist", KLength, "Final_GS")
    HapPeak, Threshold, LookupTablePath = GS[0], GS[1], GS[2]
    return HapPeak, Threshold, LookupTablePath


def Main8_FinalMerge(BiasDic, LookupTablePath1, LookupTablePath2, GCHapPeak, GAHapPeak, GCCutoff, GACutoff): # multiplicities.bed + KmerCounting or AsmCp + KmerCounting 
    LTDic = P3_LookupTable2Dic(LookupTablePath1)
    ##############
    f0 = open('10Xmulti.merged.bed.KStar.sorted.NKmer', 'r')
    f1 = open('10Xmulti.merged.bed.KStar.sorted.NKmer.NKStar', 'w')
    Line = 0
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        if '' in DT: DT.remove('')
        else:
            if DT[8] == 'ND in KmerDB': 
                f1.write(DL + '\n')
                continue
        ACN, Count, NCount, GCRatio, GARatio, NCountGC, NCountGA = int(DT[4]), DT[5], DT[8], DT[9], 0, round(float(DT[11]), 0), 0
        GCKStarsCal = P3_KStarCal(NCountGC, ACN, LTDic, GCHapPeak, GCCutoff)
        GCKStar, GCExpACN = GCKStarsCal[0], GCKStarsCal[1]
        #GAKStarsCal = P3_KStarCal(NCountGA, ACN, LTDic, GAHapPeak, GACutoff)
        #GAKStar, GAExpACN = GAKStarsCal[0], GAKStarsCal[1]
        GAKStar, GAExpACN = 0, 0
        f1.write(DL + "\t" +  str(GCExpACN) + '\t' + str(GCKStar) + '\t' + str(GAExpACN) + '\t' + str(GAKStar) + '\n')

    #os.system("sort -m -t -k1,1 -nk2,2 ?? 10Xmulti.merged.bed.KStar.sorted.NKmer > 10Xmulti.merged.bed.KStar.sorted.NKmer.sorted")

def Post1_AsmKmertoGCKmerDic():
    GCDic = {} # GCPercent:{Kmer:META} 
    f0 = open('10Xmulti.merged.bed.KStar.scaffold1.sorted.NKmer.NKStar', 'r')
    Line = 0
    while True:
        Line += 1
        if Line % 100000 == 0: print(Line)
        DL = f0.readline().split('\n')[0]
        if not DL: break
        if DL == '': break
        DT = DL.split("\t")
        if '' in DT: DT.remove('')
        else:
            if DT[8] == 'ND in KmerDB': continue
        IKmer, GCRatio, Count, NCountGC, KStar, NKStar = int(DT[3]), int(math.floor(float(DT[9]))), int(DT[5]), int(round(float(DT[11]),0)), DT[7], DT[14]
        Meta = [Count, NCountGC, KStar, NKStar]
        if GCRatio in GCDic:
            if IKmer in GCDic[GCRatio]: pass
            else: GCDic[GCRatio][IKmer] = Meta
        else: GCDic[GCRatio] = {IKmer:Meta}
    BinSaving('./', GCDic, "GCKmerDic")


def Post2_RandomSampling(SampleSize):
    GCDic = BinLoading('./', "GCKmerDic")
    print("Load Done")
    f0 = open('RandomSampledKmerDat.txt', 'w')
    for GCRatio in GCDic:
        IKmerList = list(GCDic[GCRatio].keys())
        RanList = random.choice(IKmerList, size=SampleSize, replace=True)
        print("Sampling Done")
        for RanKmer in RanList:
            DL = GCDic[GCRatio][RanKmer]
            f0.write(str(GCRatio) + '\t' + "\t".join(list(map(str, DL))) + '\n')




# bTG2 #
BamPath, AsmPath1, AsmPath2, CuttedMergedRead1, MergedRead2, KLength = 'bTaeGut2_trio.rebinned.hap1.s2.fasta.bam', 'bTaeGut2_trio.rebinned.hap1.s2.fasta', 'bTaeGut2_trio.rebinned.hap2.s2.fasta', 'Merged_10X_R1_rmadt.fastq', 'Merged_10X_R2_.fastq', 20
FlankingSize, ZeroPropCutoff = 50, 0.5
Cores = 20
WeightCountPosition = 5
Pre1_DepthDicConstruction(BamPath)
Pre2_MerquryMultiplicities(AsmPath1, AsmPath2, CuttedMergedRead1, MergedRead2, KLength)

RunGS = P3_DoGenomeScope('./merqury/10XRead.hist', KLength, 'bTG2')
HapPeak, Cutoff, LookupTablePath = RunGS[0], RunGS[1], RunGS[2]
print(HapPeak)
print(Cutoff)
print(LookupTablePath)
Pre3_KStarCal(LookupTablePath, HapPeak, Cutoff)
Pre4_Pri1Cp(AsmPath1, KLength)
Main1_Get1CP0KStar(FlankingSize, ZeroPropCutoff)
Main2_GetGCPropAndDepth(AsmPath1, KLength)
BiasDic = Main3_LocusDepthBiasDic()
Main4_Dic2Txt(BiasDic[1], "GCBiasDic.txt")
#Main4_Dic2Txt(BiasDic[2], "GABiasDic.txt")

#### [read parsing] ####
Main5_Read2Kmer(CuttedMergedRead1, MergedRead2, KLength, Cores)

BiasDic = LoadGCBiasHist("GCBiasDic.txt")
BinSaving('./', BiasDic, "BiasDic")
Main6_KmerCounting_PrefixMerge(BiasDic, "Merged_10X_R1")  #NTC
Main6_KmerCounting_PrefixMerge_MultiProc(BiasDic, "Merged_10X_R2")
Main6_KmerCounting_FinalMerge(BiasDic)
Main7_3lvCount2DicGC("Finalintermerged")

Main7_NKmerSync("10Xmulti.merged.bed.KStar", "Finalintermerged.sorted.merged")

WeightGS = M8_GetNormalizedLookupTable(WeightCountPosition, KLength)
HapPeak, Cutoff, LookupTablePath = WeightGS[0], WeightGS[1], WeightGS[2] 
print(HapPeak)
print(Cutoff)
print(LookupTablePath)
Main8_FinalMerge(BiasDic, LookupTablePath, LookupTablePath, HapPeak, HapPeak, Cutoff, Cutoff)

Post1_AsmKmertoGCKmerDic()
Post2_RandomSampling(10000)

