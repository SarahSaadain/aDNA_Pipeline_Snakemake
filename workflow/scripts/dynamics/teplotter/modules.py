from io import StringIO, TextIOBase
from collections import defaultdict
import logging
import re


class FileWriter:

    def __init__(self,outfile):
        self.outfile=outfile
        self.file_handle=None
        self.should_close=False
        if outfile is not None:
            self.file_handle=open(outfile,"w")
            self.should_close=True
        
    def write(self,towrite:str):
        if self.file_handle is not None:
            self.file_handle.write(towrite+"\n")
        else:
            print(towrite)
    
    def __exit__(self):
            if self.should_close and self.file_handle is not None:
                self.file_handle.close()
                self.file_handle = None

class SequenceEntryReader:
    """
    Simple iterator over SeqEntry file that yields one record at a time.

    
    Usage:
        for se in SeqEntry("seqentryfile.se"):
            ...
    """
    
    def __init__(self, file):
        self.file = file
        self._file = None
        self._should_close = False

    def __iter__(self):
        self._open_file()
        self._activeSeq = None
        self._activeLines = []
        return self

    def __next__(self):
        while True:
            line = self._file.readline()
            if not line:
                # End of file — yield last record if any
                if self._activeSeq:
                    se = SequenceEntry.parse(self._activeLines)
                    self._activeSeq = None
                    self._activeLines = []
                    return se
                raise StopIteration

            line = line.rstrip('\n\r')
            line = line.strip()
            seqName=line.split("\t")[0]
            

            # first record; initialize
            if self._activeSeq is None:
                self._activeSeq=seqName
                self._activeLines=[line]
            # new record; safe and start new one
            elif seqName!=self._activeSeq:
                # New record starts
                # Yield previous record
                se = SequenceEntry.parse(self._activeLines)
                self._activeSeq=seqName
                self._activeLines = [line]
                return se
            elif line:  # skip empty lines
                self._activeLines.append(line)

    def _open_file(self):
        if self._file is not None:
            return
        
        if hasattr(self.file, 'readline'):
            # already a file object
            self._file = self.file
            self._should_close = False
        else:
            # assume it's a path
            path = self.file
            if path.endswith(('.gz', '.gzip')):
                import gzip
                self._file = gzip.open(path, 'rt')
            else:
                self._file = open(path, 'r')
            self._should_close = True

    

    def close(self):
        if self._should_close and self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        self._open_file()
        return self

    def __exit__(self): 
        self.close()


def is_snp(refc,hash,cov,minc,minfreq):
    # refc = A
    # hash =
    if cov==0:
        return False
    if 'A'!= refc:
        ac=hash['A']
        af=float(ac)/float(cov)
        if ac>=minc and af>=minfreq:
            return True
    if 'T' != refc:
        tc=hash['T']
        tf=float(tc)/float(cov)
        if tc>=minc and tf >=minfreq:
            return True
    if 'C'!= refc:
        cc=hash['C']
        cf=float(cc)/float(cov)
        if cc>=minc and cf>=minfreq:
            return True
    if 'G' != refc:
        gc=hash['G']
        gf=float(gc)/float(cov)
        if gc>=minc and gf >=minfreq:
            return True
    return False


class Indel:
    @classmethod
    def parse(cls,e):
        if len(e)!=5:
            raise Exception(f"Cannot parse Indel {e}")
        ref,type,pos,length,count=e[0],e[1],int(e[2]),int(e[3]),float(e[4])
        return Indel(ref,type,pos,length,count)
 
    def __init__(self,ref:str,type:str,pos:int,length:int,count):
        self.ref=ref
        self.type=type # ins or del
        self.pos=pos
        self.length=length
        self.count=count
    
    def __str__(self):
        # ref, type, pos, length, count
        tp=[self.ref, self.type, f"{self.pos}",f"{self.length}",f"{self.count:.2f}"]
        tp="\t".join(tp)
        return tp
    
    def normalize(self, normfactor:float):
        ni=Indel(self.ref,self.type,self.pos,self.length,float(self.count)/normfactor)
        return ni


class SNP:
    @classmethod
    def parse(cls,e):
        if len(e)!=8:
            raise Exception(f"Cannot parse SNP {e}")
        ref,type,pos,refc,ac,tc,cc,gc =e[0],e[1],int(e[2]),e[3],float(e[4]),float(e[5]),float(e[6]),float(e[7])
        return SNP(ref,pos,refc,ac,tc,cc,gc)


    def __init__(self,ref:str,pos:int,refc:str,ac,tc,cc,gc):
        self.ref=ref
        self.pos=pos
        self.refc=refc
        self.ac=ac
        self.tc=tc
        self.gc=gc
        self.cc=cc
    
    def __str__(self):
        # ref, 'snp', pos, refc, ac tc cc gc
        tp=[self.ref, "snp",  f"{self.pos}", self.refc,f"{self.ac:.2f}",f"{self.tc:.2f}",f"{self.cc:.2f}",f"{self.gc:.2f}"] 
        tp="\t".join(tp)
        return tp
    
    def normalize(self,normFactor:float):
        acn=float(self.ac)/normFactor
        tcn=float(self.tc)/normFactor
        ccn=float(self.cc)/normFactor
        gcn=float(self.gc)/normFactor
        ns=SNP(self.ref,self.pos,self.refc,acn,tcn,ccn,gcn)
        return ns


class SequenceEntry:

    @classmethod 
    def parse(cls,lines):
        activeName=None
        covar=None
        ambcovar=None
        snplist=[]
        indellist=[]

        for l in lines:
            tmp=l.split("\t")
            sn=tmp[0]
            if activeName is None:
                activeName=sn
            assert sn == activeName
            feature=tmp[1]
            if feature=="ins" or feature == "del":
                indel=Indel.parse(tmp)
                indellist.append(indel)
            elif feature == "snp":
                snp=SNP.parse(tmp)
                snplist.append(snp)
            elif feature =="cov":
                if covar is not None:
                    raise Exception(f"two coverage arrays for sequence {sn}")
                covar = [float(x) for x in tmp[2].split()]
            elif feature =="ambcov":
                if ambcovar is not None:
                    raise Exception(f"two amb coverage arrays for sequence {sn}")
                ambcovar = [float(x) for x in tmp[2].split()]
            else:
                raise Exception(f"Unknown feature {feature}")
        if covar is None:
            raise Exception(f"No coverage for {activeName}")
        if ambcovar is None:
            raise Exception(f"No ambiguous coverage for {activeName}")
        return SequenceEntry(activeName,covar,ambcovar,snplist,indellist)
    
    def __init__(self,sequence_name:str,coverage,ambiguous_coverage,snp_list,indel_list):
        self.sequence_name=sequence_name
        self.coverage=coverage
        self.ambiguous_coverage=ambiguous_coverage
        self.snp_list=snp_list
        self.indel_list=indel_list
    
    def __str__(self):
        # cov
        tmp=" ".join([f"{i:.2f}" for i in self.coverage])
        tpcov="\t".join([self.sequence_name,"cov",tmp])
        #ambcov
        tmp=" ".join([f"{i:.2f}"  for i in self.ambiguous_coverage])
        tpambcov="\t".join([self.sequence_name,"ambcov",tmp])
        tp=[tpcov,tpambcov]
        for s in self.snp_list:
            tp.append(str(s))
        for id in self.indel_list:
            tp.append(str(id))
        topr="\n".join(tp)
        return topr
    
    def normalize(self,normfactor:float):
        cov=[float(i)/normfactor for i in self.coverage]
        ambcov=[float(i)/normfactor for i in self.ambiguous_coverage]
        snplist=[]
        for s in self.snp_list:
            snplist.append(s.normalize(normfactor))
        indellist=[]
        for i in self.indel_list:
            indellist.append(i.normalize(normfactor))
        return SequenceEntry(self.sequence_name,cov,ambcov,snplist,indellist)


class SequenceEntryBuilder:
    def __init__(self,ref_sequence_string:str,ref_sequence_name:str,min_mapping_quality:int):
        self.ref_sequence_string=ref_sequence_string
        self.ref_sequence_length=len(ref_sequence_string)
        self.ref_sequence_name=ref_sequence_name

        # following variables are filled by add_read and used to create a SequenceEntry by to_SequenceEntry
        
        self.snp_overview_list=[{'A':0,'T':0,'C':0,'G':0} for i in list(ref_sequence_string)]
        self.coverage_list=[0 for i in list(ref_sequence_string)]
        self.ambiguous_coverage=[0 for i in list(ref_sequence_string)]
        self.insertions_list=[]
        self.deletions_list=[]
        self.minmapq=min_mapping_quality

    def __parse_cigar(self,cigar: str) -> list[tuple[str, int]]:
        """Parse CIGAR string into list of (op, length) tuples."""
        ops = []
        i = 0
        num = ""
        while i < len(cigar):
            if cigar[i].isdigit():
                num += cigar[i]
            else:
                if num:
                    ops.append((cigar[i], int(num)))
                    num = ""
            i += 1
        return ops
    
    def __add_coverage(self,alignment_start_pos: int,cigar_list,mapq:int):

        reference_position = alignment_start_pos-1
        query_position = 0
        ### ref     ATTTAAACCCC---AAAA
        ### que.    ATTT---CCCCTTTAAAA
        ###             D      I
        for cigar_operator, length in cigar_list:
            if cigar_operator in ('H', 'S'):  # Hard/soft clip: consumes query only; does not add to coverage
                query_position += length
            elif cigar_operator == 'I':  # Insertion: consumes query only; does not add to coverage
                query_position += length
            elif cigar_operator in('D','N'):  # Deletion: consumes reference only; does add to coverage
                for i in range(length):
                    p = reference_position + i

                    if p >= self.ref_sequence_length:
                        break
                    
                    self.coverage_list[p]+=0 # should i add to coverage -> at moment no

                    if mapq<self.minmapq:
                        self.ambiguous_coverage[p]+=0 # should i add to coverage -> at moment no
                
                reference_position += length

            elif cigar_operator in ('M', '=', 'X'):  # Match/mismatch: consumes both; adds coverage
                for i in range(length):
                    p=reference_position+i
                    
                    if p>=self.ref_sequence_length:
                        break

                    self.coverage_list[p]+=1

                    if mapq<self.minmapq:
                        self.ambiguous_coverage[p]+=1

                reference_position += length
                query_position += length

    def __add_indels(self,alignment_start_pos:int, cigar_list):
        reference_position = alignment_start_pos-1
        query_position = 0
        
        ### ref     ATTTAAACCCC---AAAA
        ### que.    ATTT---CCCCTTTAAAA
        ###             D      I
        for cigar_operator, length in cigar_list:
            if cigar_operator in ('H', 'S'):  # Hard/soft clip: consumes query only; does not add to coverage
                query_position += length
            elif cigar_operator == 'I':  # Insertion: consumes query only; does not add to coverage
                self.insertions_list.append((reference_position,length))
                query_position += length
            elif cigar_operator in('D','N'):  # Deletion: consumes reference only; does add to coverage
                self.deletions_list.append((reference_position,length))
                reference_position += length
            elif cigar_operator in ('M', '=', 'X'):  # Match/mismatch: consumes both
                reference_position += length
                query_position += length

    def __add_snps(self,alignment_start_pos:int,cigar_list,sequence:str):
        
        ### ref     ATTTAAACCCC---AAAA
        ### que.    ATTT---CCCCTTTAAAA
        reference_position=alignment_start_pos-1
        query_position = 0

        for cigar_operator, length in cigar_list:
            if cigar_operator in ('H', 'S'):  # Hard/soft clip: consumes query only
                query_position += length
            elif cigar_operator == 'I':  # Insertion: consumes query only
                query_position += length
            elif cigar_operator in('D','N'):  # Deletion: consumes reference only
                reference_position += length
            elif cigar_operator in ('M', '=', 'X'):  # Match/mismatch: consumes both
                for i in range(length):
                    base = sequence[query_position + i]
                    if base.upper() in 'ATCG':
                        p=reference_position+i
                        if p>=self.ref_sequence_length:
                            break
                        self.snp_overview_list[p][base]+=1
                reference_position += length
                query_position += length
            # Ignore N (skipped reference), P (padding) if present

    
    def add_read(self,alignment_start_pos:int,cigar:str,mapping_quality:int,sequence:str):

        cigar_list = self.__parse_cigar(cigar)
        self.__add_coverage(alignment_start_pos, cigar_list, mapping_quality) # increase coverage; only cigar and mapquality considered
        self.__add_indels(alignment_start_pos, cigar_list)        # add indels; only cigar considered; mapq ignored
        self.__add_snps(alignment_start_pos, cigar_list, sequence)      # add snps; only cigar considered; mapq ignored
    
    def to_SequenceEntry(self,snp_min_count, snp_min_frequency, indel_min_count, indel_min_frequency):
        
        snp_list=[]
        for i,snp in enumerate(self.snp_overview_list):
            refc = self.ref_sequence_string[i]
            cov = self.coverage_list[i]
            if is_snp(refc,snp,cov,snp_min_count,snp_min_frequency):
                snp_entry=SNP(self.ref_sequence_name,i+1,refc,snp['A'],snp['T'],snp['C'],snp['G']) # one based snp position
                snp_list.append(snp_entry)
        
        indel_list=[]
        # INSERTIONS
        tmp=defaultdict(int)
        for ins in self.insertions_list:
            tmp[ins]+=1

        for ins, count in tmp.items():
            pos=ins[0]-1 # position in ins is 1-based but coverage is 0-based
            cov=self.coverage_list[pos]
            insfreq=float(count)/float(cov)
            if count>=indel_min_count and insfreq>=indel_min_frequency:
                id=Indel(self.ref_sequence_name,"ins",ins[0],ins[1],count)
                indel_list.append(id)

        # DELETIONS; kept separate on purpose; in case I want to treat them differentially later
        tmp=defaultdict(int)
        for de in self.deletions_list:
            tmp[de]+=1

        for de,count in tmp.items():
            pos=de[0]-1 # position in ins is 1-based but coverage is 0-based
            cov=self.coverage_list[pos]
            defreq=float(count)/float(cov)
            if count>=indel_min_count and defreq>=indel_min_frequency:
                id=Indel(self.ref_sequence_name,"del",de[0],de[1],count)
                indel_list.append(id)

        se=SequenceEntry(self.ref_sequence_name,self.coverage_list,self.ambiguous_coverage,snp_list,indel_list)
        return se

class NormFactor:

    def _get_coverage_triplet(cov:list,qlen:int):
        if qlen==0:
            mcov=float(sum(cov))/float(len(cov))
            return [mcov,None,None]
        
        cov.sort()
        first=  cov[:qlen]
        middle=  cov[qlen:-qlen]
        last=  cov[-qlen:]
        mfirst=float(sum(first))/float(len(first))
        mmiddle=float(sum(middle))/float(len(middle))
        mlast=float(sum(last))/float(len(last))
        return [mmiddle,mfirst,mlast]



    @classmethod
    def get_coverage_stat(cls, se, minDistance:int, quantile:int):
        assert quantile<50 and quantile>=0
        assert minDistance >=0
        cov=se.cov
        ambcov=se.ambcov
        if minDistance>0:
                cov=cov[minDistance:-minDistance]
                ambcov=ambcov[minDistance:-minDistance]
        assert len(cov)==len(ambcov)
        if len(cov)==0:
            return [None,]*6
        qfrac=float(quantile)/100.0
        qlen=int(len(cov)*qfrac)
        covtrip=NormFactor._get_coverage_triplet(cov,qlen)
        ambcovtrip=NormFactor._get_coverage_triplet(ambcov,qlen)
        toret=[]
        toret.extend(covtrip)
        toret.extend(ambcovtrip)
        return toret


    @classmethod
    def compute_normalization_factor_for_file(cls, filename:str, scg_suffix:str, min_end_distance: int, quanitle:int):
        # compute the normalization factor from a seq-entry file (seq overview file so-file)
        scgs=[]

        for se in SequenceEntryReader(filename):
            if se.sequence_name.endswith(scg_suffix):
                scgs.append(se)

        logging.debug(f"Found {len(scgs)} single copy genes with suffix '{scg_suffix}' in file {filename} for normalization.")

        if len(scgs)==0:
            logging.error(f"No single copy genes found with suffix '{scg_suffix}' in file {filename}. Cannot compute normalization factor.")
            raise Exception("No single copy genes found for normalization")
        
        normfactor = NormFactor.compute_normalization_factor_for_sequence_entries(scgs,min_end_distance,quanitle)
        return normfactor

    @classmethod
    def compute_normalization_factor_for_sequence_entries(cls, sequence_entities: list[SequenceEntry], ignore_read_ends_base_count:int,quantile:int) -> float:
        assert quantile<50 and quantile>=0
        assert ignore_read_ends_base_count >=0

        logging.debug(f"Computing normalization factor using {len(sequence_entities)} single copy genes. Ignoring ends of single copy genes by {ignore_read_ends_base_count} bases and excluding {quantile}% of the most extreme coverage values based on quantiles.")

        # compute normalizatino factor for seq-entries
        coverages_list=[]
        for sequence_entity in sequence_entities:
            # ignore the ends of the entries
            if len(sequence_entity.coverage) <= 2 * ignore_read_ends_base_count:
                logging.debug(f"Skipping single copy gene '{sequence_entity.sequence_name}' for normalization factor computation because its length ({len(sequence_entity.coverage)}) is less than or equal to twice the minimum end distance ({2 * ignore_read_ends_base_count}).")
                continue

            total_coverage = sequence_entity.coverage
            if ignore_read_ends_base_count > 0:
                # exclude the ends of the scgs
                total_coverage = total_coverage[ignore_read_ends_base_count:-ignore_read_ends_base_count]

            coverages_list.extend(total_coverage)

        # finaly exclude the quantiles of the largest and smallest coverages        
        coverages_list.sort()
       
        if quantile>0:
            quantile_fraction=float(quantile)/100.0
            quantile_drop_count = int(len(coverages_list) * quantile_fraction)

            logging.debug(f"Excluding {quantile_drop_count} coverage values from the lower and the upper end for normalization factor computation based on quantile {quantile}%.")

            coverages_list = coverages_list[quantile_drop_count:-quantile_drop_count]

        if len(coverages_list)==0:
            raise Exception("Unable to normalize; no valid coverage for a single copy gene")

        mean = float(sum(coverages_list))/float(len(coverages_list))
        return mean

def load_fasta(fafile):
    entries = {}
    current_header = None
    current_sequence = []
    fh=None
    if  isinstance(fafile, TextIOBase):
        fh=fafile
    else:
        fh=open(fafile,'r')
    
    for line in fh:
        line = line.rstrip()  # remove trailing \n
        
        if line.startswith('>'):
            # Save previous entry if exists
            if current_header is not None:
                seq = ''.join(current_sequence)
                entries[current_header]=seq
            
            # Start new entry; get rid of >
            current_header = line[1:]
            # split and get rid of anything after whitespace
            if re.search(r'\s', current_header):
                current_header=re.split(r'\s+', current_header)[0]
            current_sequence = []
        elif line and current_header is not None:
            # Add sequence line (skip empty lines)
            current_sequence.append(line)
    
    # Don't forget the last entry!
    if current_header is not None:
        seq = ''.join(current_sequence)
        entries[current_header]=seq
    fh.close()
    
    return entries

def test_computeNormalization():
    # todo add tests for quantiles
    ses=[SequenceEntry("t",[10,]*10,[],[],[]),SequenceEntry("t",[2,]*10,[],[],[])]
    nf=NormFactor.compute_normalization_factor_for_sequence_entries(ses,0,0)
    assert nf==6, "test1"

    ses=[SequenceEntry("t",[1,10,10,10,10,10,1],[],[],[]),SequenceEntry("t",[1,2,2,2,2,2,1],[],[],[])]
    nf=NormFactor.compute_normalization_factor_for_sequence_entries(ses,0,0)
    assert nf<5, "test2"
    nf=NormFactor.compute_normalization_factor_for_sequence_entries(ses,1,0)
    assert nf==6, "test3"

    print("Quick test computation of normalization factor passed ✓")

def test_covstat():
    se=SequenceEntry("t",[0,2,2,2,2,2,1,2,3,2,2,0],[99,5,5,5,6,4,5,5,5,5,5,99],[],[])
    cs=NormFactor.get_coverage_stat(se,1,10)

   

    assert cs[0]==2
    assert cs[1]==1
    assert cs[2]==3
    assert cs[3]==5
    assert cs[4]==4
    assert cs[5]==6
    print("Quick test computations of coverage statistic passed ✓")

def test_normalize():
    s=SNP("chr1",1,"A",5,6,7,1)
    sn=s.normalize(2.0)
    assert sn.ref=="chr1"
    assert sn.pos== 1
    assert sn.refc=="A"
    assert sn.ac==2.5
    assert sn.tc==3
    assert sn.cc==3.5
    assert sn.gc==0.5

    print("Quick test of SNP normalization PASSED ✓")
    id=Indel("chr2", "ins",5,2,11)
    idn=id.normalize(2.0)
    assert idn.ref == "chr2"
    assert idn.type=="ins"
    assert idn.pos==5
    assert idn.count==5.5
    assert idn.length==2
    print("Quick test of Insertion normalization PASSED ✓")
    deli=Indel("chr3", "del",5,2,20)
    de=deli.normalize(5.0)
    assert de.ref == "chr3"
    assert de.type=="del"
    assert de.pos==5
    assert de.count==4
    assert de.length==2
    print("Quick test of Deletion normalization PASSED ✓")

    id=Indel("chr2", "ins",5,2,11)
    deli=Indel("chr3", "del",5,2,20)
    s=SNP("chr1",1,"A",5,6,7,1)
    se=SequenceEntry("te1",[5,6,6,4,2],[2,3,4,6,1],[s],[id,deli])
    sn=se.normalize(2)
    assert sn.coverage[0]==2.5
    assert sn.coverage[1]==3
    assert sn.coverage[4]==1
    assert sn.ambiguous_coverage[0]==1
    assert sn.ambiguous_coverage[1]==1.5
    assert sn.ambiguous_coverage[4]==0.5
    assert sn.ambiguous_coverage[3]==3
    assert sn.snp_list[0].ac==2.5
    assert sn.indel_list[0].count==5.5
    assert sn.indel_list[1].count==10
    print("Quick test of SeqEntry normalization PASSED ✓")




def test_Seq_Builder_add():
    sb=SequenceEntryBuilder("AAATTTCCCGGG","hans",5)
    sb.add_read(1,"3M",4,"AAA")
    sb.add_read(1,"3M",5,"TTT")

    assert sb.coverage_list[0]==2
    assert sb.ambiguous_coverage[0]==1
    assert sb.coverage_list[1]==2
    assert sb.ambiguous_coverage[1]==1
    assert sb.coverage_list[2]==2
    assert sb.ambiguous_coverage[2]==1
    assert sb.coverage_list[3]==0
    assert sb.ambiguous_coverage[3]==0
    assert sb.snp_overview_list[0]['A']==1
    assert sb.snp_overview_list[0]['T']==1
    

    # AAATTT---CCCGGG
    # 123456---789012
    #    TTTAAACCC
    sb.add_read(4,"3=3I3X",5,"TTTAAACCC")
    assert sb.coverage_list[3]==1
    assert sb.coverage_list[4]==1
    assert sb.coverage_list[5]==1
    assert sb.coverage_list[6]==1
    assert sb.coverage_list[7]==1
    assert sb.coverage_list[8]==1
    assert sb.coverage_list[9]==0
    assert sb.snp_overview_list[3]['T']==1
    assert sb.snp_overview_list[6]['A']==0
    assert sb.snp_overview_list[6]['C']==1
    assert sb.insertions_list[0]==(6,3),  f"got {sb.insertions_list[0]}"


    # AAATTTCCCGGG
    # 123456789012
    #    TTT---AAA
    sb.add_read(4,"3M3D3M",5,"TTTAAA")
    assert sb.coverage_list[3]==2
    assert sb.coverage_list[4]==2
    assert sb.coverage_list[5]==2
    assert sb.coverage_list[6]==1
    assert sb.coverage_list[7]==1
    assert sb.coverage_list[8]==1
    assert sb.coverage_list[9]==1
    assert sb.coverage_list[10]==1
    assert sb.coverage_list[11]==1
    assert sb.deletions_list[0]==(6,3), f"got {sb.deletions_list[0]}"

    sb.add_read(12,"3M",5,"TTT")


    print("Quick test of SeqBuilder add PASSED ✓")


def test_Seq_Builder_init():
    sb=SequenceEntryBuilder("AAATTTCCCGGG","hans",5)
    assert sb.ref_sequence_string == "AAATTTCCCGGG",           f"sequence"
    assert sb.ref_sequence_name == "hans",            f"seqname"
    assert sb.minmapq ==5,                  f"minmapq"
    assert len(sb.coverage_list) == 12,             f"length of covar"
    assert len(sb.ambiguous_coverage) == 12,           f"length of ambcovar"
    assert len(sb.snp_overview_list) == 12,           f"length of ambcovar"
    assert len(sb.insertions_list) == 0,           f"length of ambcovar"
    assert len(sb.deletions_list) == 0,           f"length of ambcovar"
    print("Quick test of SeqBuilder __init__ PASSED ✓")

def test_fasta_loader():
    from io import StringIO
    
    test_content = """>seq1 some description
ACGTACGT
GCTA
>seq2
NNNNNNNNNN
>seq3 empty sequence

>seq4
ATGCATGCATGC
"""

    result = load_fasta(StringIO(test_content))
    
    expected = {
        "seq1": "ACGTACGTGCTA",
        "seq2": "NNNNNNNNNN",
        "seq3": "",
        "seq4": "ATGCATGCATGC"
    }
    
    assert len(result) == 4,           f"Expected 4 sequences, got {len(result)}"
    assert "seq3" in result,           "Missing empty sequence entry"
    assert result["seq3"] == "",       "Empty sequence should be empty string"
    assert result == expected,         "Dictionary content doesn't match expected"
    
    print("Quick test of fasta_loader PASSED ✓")

if __name__ == "__main__":
    test_fasta_loader()
    test_Seq_Builder_init()
    test_Seq_Builder_add()
    test_normalize()
    test_computeNormalization()
    test_covstat()