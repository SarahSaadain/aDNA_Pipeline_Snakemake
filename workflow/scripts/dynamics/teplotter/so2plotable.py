#!/usr/bin/env python
import argparse
import logging
from modules import SequenceEntry, SequenceEntryReader, FileWriter, NormFactor
import os

padto=9
def format_col(topr:list):
    if padto==0:
        return "\t".join(topr)
    tf=[]
    for i in range(padto):
        if(i< len(topr)):
            tf.append(str(topr[i]))
        else:
            tf.append("")
    return "\t".join(tf)

def prepareCoveragForPrint(set:list, sampleid:str,covtype:str):
    tmp=[]
    for i,c in enumerate(set):
        # seqname, sampleid, cov, pos, count
        t=[se.sequence_name,sampleid,covtype,str(i+1),str(c)]
        tmp.append(t)
       
    first,last=tmp[0],tmp[-1]
    newfirst=[first[0],first[1],first[2],first[3],"0.0"]
    newlast=[last[0],last[1],last[2],last[3],"0.0"]
    tmp.insert(0,newfirst)
    tmp.append(newlast)

    topr=[]
    for i in tmp:
        topr.append(format_col(i))
    return topr

def prepareForPrint(se:SequenceEntry, sampleid:str):
    lines=[]
    covt=prepareCoveragForPrint(se.coverage,sampleid,"cov")
    ambcovt=prepareCoveragForPrint(se.ambiguous_coverage,sampleid,"ambcov")
    lines.extend(covt)
    lines.extend(ambcovt)

    for s in se.snp_list:
        # seqname, sampleid, snp, pos, refc, ac, tc, cc, gc
        # SNP(ref,pos,refc,ac,tc,cc,gc)
        a={"A":s.ac,"T":s.tc,"C":s.cc,"G":s.gc}
        for base,count in a.items():
            if count ==0 or base==s.refc:
                continue
            tmp=[se.sequence_name,sampleid,"snp",str(s.pos), s.refc,base,str(count)]
            lines.append(format_col(tmp))
    
    for i in se.indel_list:
        if i.type=="ins":
            # seqname, sampleid, del, pos, length, count
            tmp=[se.sequence_name,sampleid,"ins",str(i.pos),str(i.length),str(i.count)]
            lines.append(format_col(tmp))
            # ref:str,type:str,pos:int,length:int,count

        elif i.type=="del":
            # seqname, sampleid, ins, startpos, endpos, startcov,endcov, count
                # AAATTTCCCGGG
                # 123456789012
                #    TTT---AAA
                # pos = 6 and len = 3
                # bow from 6 to 10 (actual 0-based coverages are 5 and 9)
            startpos=i.pos
            endpos=startpos+i.length+1
            startcov=se.coverage[startpos-1]
            endcov=se.coverage[endpos-1]
            tmp=[se.sequence_name,sampleid,"del",str(startpos),str(endpos),str(startcov),str(endcov),str(i.count)]
            lines.append(format_col(tmp))

        else:
            raise Exception(f"invalid type{i.type}")
        
    tr="\n".join(lines)
    return tr

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser(description="""           
normalizes the coverage for seqentries
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Sarah Saadain
""")
parser.add_argument('--so', type=str, default=None,dest="so", required=True, help="A sequence overview (*.so) file")
parser.add_argument("--seq-ids", type=str, required=False, dest="seqids", default="ALL", help="IDs of the entries that should be plotted; separated by comma; can also be 'ALL'")
parser.add_argument("--sample-id", type=str, required=False, dest="sampleid", default="x", help="the ID of current sample")
parser.add_argument("--prefix", type=str, required=False, dest="prefix", default="", help="the prefix for the output file")
parser.add_argument("--outdir", type=str, required=False, dest="outputdir", default=None, help="the output directory; a plotable will be written for each fasta entry")
parser.add_argument("--outfile", type=str, required=False, dest="outfile", default=None, help="output file in plotable format;")
parser.add_argument("--log-level", type=str, required=False, dest="loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)")

args = parser.parse_args()
logging.getLogger().setLevel(args.loglevel)

# check for invalid parameter combinations
if args.outfile is not None and args.outputdir is not None:
    raise Exception("invalid parameters; either provide output-dir or output-file; not both")

# prefix only with output-dir
if args.prefix != "" and args.outputdir is None:
    raise Exception("invalid parameters; prefix only works if output-dir is provided")

# initialize writer
writer=FileWriter(args.outfile)

#if no output file is provided, don't write log to screen, otherwise it will mess up the output
if args.outfile is None:
    logging.getLogger().setLevel("ERROR")

seqset=None
if "," in args.seqids:
    seqset=set(args.seqids.split(","))
else:
    seqset=set([args.seqids])

is_print_all_requested = ( args.seqids.upper() == "ALL" )

# create output directory if it does not exist
if args.outputdir is not None:
    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

prefix = args.prefix

for se in SequenceEntryReader(args.so):
    if is_print_all_requested or se.sequence_name in seqset:

        tp = prepareForPrint(se, args.sampleid)
        
        if args.outputdir is not None:
            filename=se.sequence_name
            filename=filename.replace("/","_")
            filename=filename.replace(" ","_")
            filename=prefix + filename + ".plotable"
            full_path = os.path.join(args.outputdir, filename)
            with open(full_path, "w") as f:
                f.write(tp)
        else:
            writer.write(tp)