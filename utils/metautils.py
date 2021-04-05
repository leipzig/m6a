import pandas
import os
import sys
import json
import cloudstorage as gcs



#different readlengths for trimmomatic outputs
READLENGTHS = [36,75]

#from utils import metautils
#metautils.srpMeta('SRP012098')
#study_accession            experiment_accession            sample_accession            run_accession            experiment_title            experiment_attribute            taxon_id            library_selection            library_layout            library_strategy            library_source            library_name            bases            spots            adapter_spec            avg_read_length
#      SRP141397                      SRX3980107                  SRS3205258                SRR7049033                S78_shotgun                                                   0                       RANDOM                  PAIRED -                         WGS               METAGENOMIC             S78_shotgun        229716627          1009555                                 227.54245880610765

	
#study_accession	experiment_accession	experiment_title	     experiment_desc	      organism_taxid 	organism_name	library_name	library_strategy	library_source	library_selection	library_layout	sample_accession	sample_title	instrument	instrument_model	instrument_model_desc	total_spots	total_size	run_accession	run_total_spots	run_total_bases
#SRP012098          SRX136624               GSM908336: HepG2_IFN_Input; Homo sapiens; RNA-Seq GSM908336: HepG2_IFN_Input; Homo sapiens; RNA-Seq 9606 Homo sapiens GSM908336: HepG2_IFN_Input RNA-Seq TRANSCRIPTOMIC cDNA SINGLE SRS308332	Illumina Genome Analyzer IIx Illumina Genome Analyzer IIx ILLUMINA 30351576 610197365 SRR456549 30351576 1092656736
#read in an SPR metadata

#python -c "from utils import metautils;srp=metautils.srpMeta('SRP012098');print(srp.getFilesByTreatment('IFN'));"
class srpMeta():
    def __init__(self,SRP,gstype="gs"):
        self.srp=SRP
        self.st = pandas.read_csv("metadata/"+SRP+".metadata",sep="\t")
        #library name GSM908336: HepG2_IFN_Input
        self.st[['GSM','celltype','treatment','casecontrol']] = pandas.read_csv("metadata/"+SRP+".metadata",sep="\t")["library_name"].str.replace(' ','',regex=False).str.split('[_:]',expand=True)
        self.st['run']=self.st['celltype']+'_'+self.st['treatment']+'_'+self.st['casecontrol']
        self.st.columns = self.st.columns.str.replace('.', '',regex=False)
        self.st.columns = self.st.columns.str.replace(' ', '',regex=False)
        self.process(gstype)
        #['NA06984.1.M_111124_4_1.fastq.gz', 'NA06984.1.M_111124_4_2.fastq.gz', 
    
    def process(self,gstype='gs',compressed=True):
        #create the gs specific paths
        #gs://truwl-quertermous/SRR7058289/102901.1.fastq.gz
        #gs://truwl-quertermous/SRR7058289/102901.2.fastq.gz
        #from utils.metautils import srpMeta
        #srp=srpMeta('SRP142360').process()
        if gstype=='https':
            prefix="https://storage.googleapis.com/truwl-dominissini"
        elif gstype=='gs':
            prefix="gs://truwl-dominissini"
        elif gstype=='local':
            prefix=self.srp
        else:
            raise ValueError("need a gstype of https or gs or local")
        if compressed:
            suffix='fastq.gz'
        else:
            suffix='fastq'
        self.st['sequence_file'] = self.st.apply(lambda x: "{0}/{1}/{2}_{3}_{4}.{5}".format(prefix,x['run_accession'],x['celltype'],x['treatment'],x['casecontrol'],suffix), axis=1)


    def getSRRs(self):
        return(self.st['run_accession'].tolist())
    
    def getAllFiles(self,flatten=True,compress=True,gstype='gs'):
        """
        flatten - all files go in top directory vs SRP/EXP/SRR
        compress - .gz suffix
        """
        files=[]
        filesofinterest=self.st.copy()
        files=list(filesofinterest['sequence_file'].astype(str))
        return(files)
    
    def getTreatments(self):
        return(list(set(self.st['treatment'].astype(str))))
    
    def getFilesByTreatment(self,treatment,gstype='local',flatten=True,compress=True):
        """
        flatten - all files go in top directory vs SRP/EXP/SRR
        compress - .gz suffix
        """
        self.process(gstype=gstype)
        files=[]
        filesofinterest=self.st[self.st['treatment']==treatment].copy()
        types=list(filesofinterest['casecontrol'].astype(str))
        files=list(filesofinterest['sequence_file'].astype(str))
        res = dict(zip(types, files))
        return(res)
    
    def getRunsByTreatment(self,treatment,gstype='local',flatten=True,compress=True):
        """
        flatten - all files go in top directory vs SRP/EXP/SRR
        compress - .gz suffix
        """
        self.process(gstype=gstype)
        files=[]
        filesofinterest=self.st[self.st['treatment']==treatment].copy()
        types=list(filesofinterest['casecontrol'].astype(str))
        files=list(filesofinterest['run'].astype(str))
        res = dict(zip(types, files))
        return(res)
        
    def getControlRun(self,treatment):
        run=self.st.loc[(self.st['treatment']==treatment) & (self.st['casecontrol']=='Input'),'run']
        return(run.to_string(index=False))

    def getReplicatesofMate(self,strategy,cell,mate):
        pair = "truwl_pair{}".format(mate)
        foi=self.st.loc[(self.st['library_strategy']==strategy) & (self.st['library_source']==cell), pair]
        files = list(foi.astype(str))
        return(files)
    
    def printmeripseqManifest(self,gstype='local',compress=False):
        """
        Manifest for https://github.com/kingzhuky/meripseqpipe
        """
        #Sample_ID,input_FileName,ip_FileName,Group
        #HS,HepG2_HS_Input,HepG2_HS_IP,HS
        #IFN,HepG2_IFN_Input,HepG2_IFN_IP,IFN
        #HGF,HepG2_HGF_Input,HepG2_HGF_IP,HGF
        #UV,HepG2_UV_Input,HepG2_UV_IP,UV
        print("Sample_ID,input_FileName,ip_FileName,Group")
        for treatment in self.getTreatments():
            tdict=self.getFilesByTreatment(treatment,gstype=gstype,compress=compress)
            print("{},{},{},{}".format(treatment,tdict['Input'],tdict['IP'],treatment))

    def printmeripseqConfig(self,gstype='local',compress=False):
        for treatment in self.getTreatments():
            tdict=self.getRunsByTreatment(treatment,gstype=gstype,compress=compress)
            print("['{}.input',['{}']],".format(treatment,tdict['Input']))
            print("['{}.ip',['{}']],".format(treatment,tdict['IP']))

    def printchipnfManifest(self):
        """
        Mainfest for https://github.com/guigolab/chip-nf
        """
        # sample1    sample1_run1     /path/to/sample1_run1.fastq.gz    -           H3
        # sample1    sample1_run2     /path/to/sample1_run2.fastq.gz    -           H3
        # sample1    sample1_run3     /path/to/sample1_run3.fastq.gz    -           H3
        # sample1    sample1_run4     /path/to/sample1_run4.fastq.gz    -           H3
        # sample2    sample2_run1     /path/to/sample2_run1.fastq.gz    control1    H3K4me2
        # control1   control1_run1    /path/to/control1_run1.fastq.gz   control1    input
        self.process(gstype='local')
        self.st['controlrun']=self.st.apply(lambda x: self.getControlRun(x['treatment']), axis=1)
        #self.st.apply(lambda x: "{0}\t{1}\t{2}\t{3}\t{4}".format(x['study_accession'],x['run'],x['sequence_file'],self.getControlRun(x['treatment']),x['casecontrol']), axis=1)
        print(self.st[['run_accession','run','sequence_file','controlrun','casecontrol']].to_csv(header=False,index=False,sep="\t"))

    def printnfcoreManifest(self):
        """
        Manifest for https://github.com/nf-core/chipseq
        """
        #group,replicate,fastq_1,fastq_2,antibody,control
        #WT_BCATENIN_IP,1,BLA203A1_S27_L006_R1_001.fastq.gz,,BCATENIN,WT_INPUT
        #WT_BCATENIN_IP,2,BLA203A25_S16_L002_R1_001.fastq.gz,,BCATENIN,WT_INPUT
        #WT_BCATENIN_IP,3,BLA203A49_S40_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
        #WT_INPUT,1,BLA203A6_S32_L006_R1_001.fastq.gz,,,
        #WT_INPUT,2,BLA203A30_S21_L002_R1_001.fastq.gz,,,
        #WT_INPUT,3,BLA203A31_S21_L003_R1_001.fastq.gz,,,
        self.process(gstype='local')
        self.st['controlrun']=self.st.apply(lambda x: self.getControlRun(x['treatment']) if x['casecontrol']!='Input' else '', axis=1)
        self.st['replicate']=1
        self.st['empty']=''
        print("group,replicate,fastq_1,fastq_2,antibody,control")
        print(self.st[['run','replicate','sequence_file','empty','casecontrol','controlrun']].to_csv(header=False,index=False))
        
    def getATACManifest(self):
        cells = list(set(self.st.loc[self.st['library_strategy']=='ATAC-seq', 'library_source'].astype(int)))
        manifest = {}
        manifest["atac.pipeline_type"] = "atac"
        manifest["atac.genome_tsv"] = "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv"
        manifest["atac.paired_end"] = "true"
        manifest["atac.auto_detect_adapter"] = "true"
        manifest["atac.enable_xcor"] = "true"
        manifest["atac.title"] = "SRP142360ATAC"
        manifest["atac.description"] = "ATAC-seq on human coronary artery smooth muscle cells (HCASMC)"
        for i in range(len(cells)):
            for pair in [1,2]:
                reps=self.getReplicatesofMate('ATAC-seq',cells[i],pair)
                manifest["atac.fastqs_rep{}_R{}".format(i+1,pair)] = reps
        print(json.dumps(manifest, indent=4))
#         {
#     "atac.pipeline_type" : "atac",
#     "atac.genome_tsv" : "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv",
#     "atac.fastqs_rep1_R1" : [
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF341MYG.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF106QGY.subsampled.400.fastq.gz"
#     ],
#     "atac.fastqs_rep1_R2" : [
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF248EJF.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF368TYI.subsampled.400.fastq.gz"
#     ],
#     "atac.fastqs_rep2_R1" : [
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF641SFZ.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF751XTV.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF927LSG.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF859BDM.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF193RRC.subsampled.400.fastq.gz",
#         "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF366DFI.subsampled.400.fastq.gz"
#     ],
#     "atac.fastqs_rep2_R2" : [
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF031ARQ.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF590SYZ.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF734PEQ.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF007USV.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF886FSC.subsampled.400.fastq.gz",
#          "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF573UXK.subsampled.400.fastq.gz"
#     ],
#     "atac.paired_end" : true,
#     "atac.auto_detect_adapter" : true,
#     "atac.enable_xcor" : true,
#     "atac.title" : "SRP142360ATAC",
#     "atac.description" : "ATAC-seq on primary keratinocytes in day 0.0 of differentiation"
# }
    
    def getRNAManifest(self):
        cells = list(set(self.st.loc[self.st['library_strategy']=='RNA-Seq', 'library_source'].astype(int)))
        manifest = {}
        manifest["rna.pipeline_type"] = "atac"
        manifest["rna.genome_tsv"] = "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv"
        manifest["rna.paired_end"] = "true"
        manifest["rna.auto_detect_adapter"] = "true"
        manifest["rna.enable_xcor"] = "true"
        manifest["rna.title"] = "SRP142360RNA"
        manifest["rna.description"] = "RNA-seq on human coronary artery smooth muscle cells (HCASMC)"
        for i in range(len(cells)):
            for pair in [1,2]:
                reps=self.getReplicatesofMate('RNA-Seq',cells[i],pair)
                manifest["rna.fastqs_rep{}_R{}".format(i+1,pair)] = reps
        print(json.dumps(manifest, indent=4))

    def getAllFiles(self):
        getFilesFromRunList(self.populationRuns('ALL'))

    #>>> metautils.calculateExpected('ALL')
    #Expecting 413 fastq and 253 bams
    def calculateExpected(self,population):
        fastq = self.getFilesFromRunList(self.populationRuns(population))
        bams = self.populationRuns(population)
        print("Expecting {0} fastq and {1} bams".format(len(fastq),len(bams)))
    
    #GBR FIN ALL YRI TSI 
    #>>> metautils.populationRuns('ALL')
    #['NA12286.2.MI_120126_5', 'NA10851.4.M_120208_1',
    def populationRuns(self,population):
        return(self.getSRRs())

    #Check to see if paired-end sequencing was used in each sample so appropriate stAR parameters are called
    def isPaired(self,run):
        if self.st.loc[self.st['Assay Name']==run]['Comment[LIBRARY_LAYOUT]'].all()=='PAIRED':
            return(True)
    
    #Either run star with paired or unpaired options
    def starProgram(self,run):
        if self.isPaired(run):
            return("star_paired.sh")
        else:
            return("star_unpaired.sh")
    
    #Run trimmomatic bash script with paired or unpaired options
    def trimProgram(self,run):
        if self.isPaired(run):
            return("trimmomatic_paired.sh")
        else:
            return("trimmomatic_unpaired.sh")
            
    #pairedOrSinglefastqInput('NA06984.1.M_111124_4')
    #['NA06984.1.M_111124_4_1.fastq.gz', 'NA06984.1.M_111124_4_2.fastq.gz']
    #pairedOrSinglefastqInput('NA06985.1.MI_120104_3')
    #['NA06985.1.MI_120104_3_1.fastq.gz']
    
    #set file names for trimmomatic based on whether a trimmed input or output is expected
    def pairedOrSinglefastqInput(self,run,trimmed=False,readlength=None):
        if trimmed:
            extName=".trimmed_{0}.fastq.gz".format(readlength)
        else:
            extName=".fastq.gz"
        if self.isPaired(run):
            return([run+"_1"+extName,run+"_2"+extName])
        else:
            return([run+"_1"+extName])
    

    def getPathFromName(self,name):
        return(ALL_LOOKUP[name])
    
    #ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188404/ERR188404_1.fastq.gz

    
    def getType(self,runName):
        return(self.getProperty(runName, 'Platform'))
    
    
    def getProperty(self, runName, column):
        stseries = self.st.loc[self.st['Run'] == runName, column]
        if len(self.stseries) == 1:
            return(self.stseries.to_self.string(index=False))  # pandas is ridiculous
        elif len(self.stseries) > 1:
            raise ValueError("Not expecting to match multiple run names, just 1")
        else:
            raise ValueError("Can't find that run {0}".format(runName))
    
    
    def getMemory(self, runName):
        return(int(getProperty(runName, 'size_MB')))
    
    #Up to 128 letters (uppercase and lowercase), numbers, hyphens, and underscores are allowed
    def cleanJobname(self, name):
        return(name.replace('.','_'))
    
    def getECS(self, runName, units, program):
        if program == 'minimap':
            mb = 32000
        elif program == 'star':
            superIntensive = ['NA12716.7.M_120219_6', 'NA12775.7.M_120219_5',
                              'NA12814.7.M_120219_3', 'NA11831.7.M_120219_2',
                              'NA11994.2.M_111216_7', 'NA11893.7.M_120219_3']
            if runName in superIntensive:
                mb = 64000
            else:
                mb = 48000
        elif program == 'IsoModule':
            mb = 192000
        elif program == 'samtoolsindex':
            mb = 16000
        elif program == 'samtoolsmerge':
            mb = 16000
        elif program == 'samtoolssubsample':
            mb = 16000
        elif program == 'fastx':
            mb = 8000
        elif program == 'trimmomatic':
            mb = 8000
        elif program == 'rmats':
            mb = 16000
        else:
            raise ValueError
        if units == 'bytes':
            return(1048576*mb)
        elif units == 'mb':
            return(mb)
    # 2  # number of sample
    # 3  # number of replicates in sample #1
    # sam1_rep1.bam
    # sam1_rep2.bam
    # sam1_rep3.bam
    # 3  # number of replicates in sample #2
    # sam2_rep1.bam
    # sam2_rep2.bam
    # sam2_rep3.bam
    
    #Used for a two way comparison between various samples
    def twoSampleComparisonManifest(self, biosamp1, biosamp2, filename):
        text_file = open(filename, "w")
        text_file.write("2\n")  # two-way comparison
        for biosamp in [biosamp1, biosamp2]:
            run = getRunBamsFromBioSample(biosamp, include_bai=False)
            text_file.write("{0}\n".format(len(run)))
            for replicate in run:
                text_file.write("{0}\n".format(replicate))
    
        #Used for a two way comparison between various samples
    def twoTreatmentComparisonManifest(self, treatment1, treatment2, filename):
        text_file = open(filename, "w")
        text_file.write("2\n")  # two-way comparison
        for treatment in [treatment1, treatment2]:
            run = getRunBamsFromTreatment(treatment, include_bai=False)
            text_file.write("{0}\n".format(len(run)))
            for replicate in run:
                text_file.write("{0}\n".format(replicate))
    # "panorama-clk-repro/SRP091981/
    
    #>>> metautils.getRunBamsFromBioSample('NA12778')
    #['NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam.bai', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam', 'NA12778.1.M_111124_4.Aligned.sortedByCoord.out.bam.bai', 'NA12778.1.MI_120104_3.Aligned.sortedByCoord.out.bam', 'NA12778.1.MI_120104_3.Aligned.sortedByCoord.out.bam.bai']
    #Get a List of bams originating from the same biological source with muleiples removed (useful for the merge step, where multiple runs must be collapsed down to one)
    def getRunBamsFromBioSample(self, biosamp, readlength, include_s3=None, include_bai=True):
    #    if include_bai:
    #        exts = ['bam', 'bam.bai']
    #    else:
        exts = ['bam']
        runs = self.getRunsFromBioSample(biosamp)
        if include_s3:
            bams = ["{0}/{1}.{3}_trimmed.Aligned.sortedByCoord.out.{2}".format(
                include_s3, replicate, ext, readlength) for replicate in runs for ext in exts]
        else:
            bams = ["{0}.trimmed_{2}.Aligned.sortedByCoord.out.{1}".format(
                replicate, ext, readlength) for replicate in runs for ext in exts]
        return(bams)
    
    def getRunBamsFromTreatment(self, treatment, readlength, include_s3=None, include_bai=True):
    #    if include_bai:
    #        exts = ['bam', 'bam.bai']
    #    else:
        exts = ['bam']
        runs = self.getRunsFromTreatment(treatment)
        if include_s3:
            bams = ["{0}/{1}.{3}_trimmed.Aligned.sortedByCoord.out.{2}".format(
                include_s3, replicate, ext, readlength) for replicate in runs for ext in exts]
        else:
            bams = ["{0}.trimmed_{2}.Aligned.sortedByCoord.out.{1}".format(
                replicate, ext, readlength) for replicate in runs for ext in exts]
        return(bams)
    #>>> metautils.getRunsFromBioSample('NA12778')
    #['NA12778.1.M_111124_4', 'NA12778.1.M_111124_4', 'NA12778.1.MI_120104_3']
    #Get a List of runs from the same source with muleiples removed (useful for the merge step, where multiple runs must be collapsed down to one)
    def getRunsFromBioSample(self, biosample,include_bai = True,allowSingle=True):
        return(List(set(self.st.loc[(self.st['Source Name'] == biosample) & (self.st['Comment[Quality Control passed]'] == 1)]['Assay Name'].toList())))
    
    def getRunsFromTreatment(self, treatment,include_bai = True,allowSingle=True):
        return(List(set(self.st.loc[(self.st['treatment'] == treatment) & (self.st['Comment[Quality Control passed]'] == 1)]['Assay Name'].toList())))

    #metautils.populationBiosamples('ALL')
    #['NA12778', 'NA12045', 'NA12144', ...
    #Returns either all runs that pass the QC check for a specific population or all populations depending on arguments
    def populationBiosamples(self,population):
        if(population=='ALL'):
            return(list(set(self.st.loc[(self.st['Comment[Quality Control passed]'] == 1)]['Source Name'])))
        else:
            return(list(set(self.st.loc[(self.st['Characteristics[population]'] == population) & (self.st['Comment[Quality Control passed]'] == 1)]['Source Name'])))
    
    #This returns all sequencing runs over 70 nt (longruns) looks through the population being considered appends them to the sample List using the mergedbam file naming convention, and then appends all runs from that population with a 36nt trimmed suffix
    def getListOfLengthOutputs(self,population):
        x = 0
        biosamps = self.populationBiosamples(population)
        biosamp_lengthList = []
        longruns = list(set(self.st.loc[self.st['Comment[SEQUENCE_LENGTH]'].astype(int)<70]['Source Name'].toList()))
        for item in longruns:
            if (item in biosamps) and (item != 'NA07000'):
                biosamp_lengthList.append("trimmed_75nt/"+item +"/"+item + ".rmats")
        for item in biosamps:
            if item != 'NA07000':
                biosamp_lengthList.append("trimmed_75nt/"+item +"/"+item + ".rmats")
        return(biosamp_lengthList)

def GCSExists(gcs_file):
    '''
    True if file exists; pass complete /bucket/file
    '''
    try:
        file = gcs.open(gcs_file,'r')
        file.close()
        status = True
    except:
        status = False
    return status

if __name__== "__main__":
    srpMeta(sys.argv[1]).get16SFiles()
    
