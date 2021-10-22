from . import __path__
import sys
import subprocess
from glob import glob
import numpy as np
import time
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn,TimeRemainingColumn
from rich.console import Console
from rich.table import Table
from multiprocessing import Pool, Process, Manager
import os
from pandas.errors import EmptyDataError
import pandas as pd
import anndata
import argparse
from argparse import RawTextHelpFormatter
from functools import partial
from shutil import which

import logging
from rich.logging import RichHandler
FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")

fqcounter = os.path.join(__path__[0], "fqcount")

help = 'Full bulk pipeline, from fastq to adata count matrix!\n'\
    'Performs the following: fastq -STAR-> bam -featureCounts-> anndata.h5ad'

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument("--fq_path", help="Path for input fastq files (relative, default: fastq).")
#parser.add_argument("--bam_path", help="Path for STAS (default: fastq).")
parser.add_argument("--star_ref", help="STAR index path.")
parser.add_argument("--gtf", help="GTF file path for featureCounts.")
parser.add_argument("--n_threads","-n", help="number of threads per fastq file, for both STAR and featureCounts.")
parser.add_argument("--adata_out", help="Path for the adata output (relative, default: adata_bulk_star.h5ad).")



def ref_loader(star_ref):
    loadref="STAR --genomeLoad LoadAndExit --genomeDir %s" % star_ref
    proc=subprocess.Popen(loadref.split(),
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    proc.wait()
    
def ref_remover(star_ref):
    loadref="STAR --genomeLoad Remove --genomeDir %s" % star_ref
    proc=subprocess.Popen(loadref.split(),
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    proc.wait()

def runcom(c):
    out=subprocess.check_output(c.split())
    return int(out)
   
def run_star(n_threads,star_ref,fq_path,sample):
    fqs=" ".join(glob(os.path.join(fq_path,sample+"*.gz")))
    
    runstar="STAR --runThreadN %s --limitBAMsortRAM 10000000000 --genomeLoad LoadAndKeep --genomeDir %s --readFilesIn %s --outFileNamePrefix aligned/%s_ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard" %(n_threads,star_ref,fqs,sample)
    proc=subprocess.Popen(runstar.split(),
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    proc.wait()
    return proc.returncode,proc.communicate()[1]


def gather_starlogs(f):
    a_file = open(f, "r")
    lines = a_file.readlines()
    if len(lines)<=2:
        return 0
    if len(lines)==3 & lines[-1] == 'ALL DONE!\n':
        return 0
    else:
        l = 2 if lines[-1] == 'ALL DONE!\n' else 1
        return int(lines[-l].split("    ")[2])

def run_fc(n_threads,gtf,bam):
    name=bam.split("/")[1].split("_")[0]
    runfc="featureCounts -T %s -a %s -g gene_name -o fc/%s.txt %s" %(n_threads,gtf,name,bam)
    proc=subprocess.Popen(runfc.split(),
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    proc.wait()

def run_fc_par(n_threads,gtf):
    pool = Pool()
    func = partial(run_fc, n_threads,gtf)
    pool.map(func, glob("aligned/*.bam"))
    

def fc2adata(p):
    try:
        df = pd.read_csv(p, sep='\t', comment='#',usecols=[0,6],index_col=0)
        df.columns=[p.split(".")[0]]
        return df
    except EmptyDataError:
        pass
    
def fc2adata_par(out,samples):
    pool = Pool()
    fcs = pool.map(fc2adata, glob("fc/*.txt"))
    allcounts = pd.concat(fcs,axis=1)
    adata = anndata.AnnData(allcounts.T)
    adata.obs_names = samples
    adata.write_h5ad(out)

def main():
    console = Console()
    console.print("bulktools 0.1",style="bold")
    args = parser.parse_args()
    manager = Manager()
    return_dict = manager.dict()
    
    if which("STAR") is None:
        log.error("STAR not installed!")
        sys.exit(1)
        
    if which("featureCounts") is None:
        log.error("subread not installed!")
        sys.exit(1)
    
    if args.star_ref is not None:
        star_ref = args.star_ref
    else:
        log.error("STAR ref missing!")
        parser.print_help(sys.stderr)
        sys.exit(1)
    if args.gtf is not None:
        gtf = args.gtf
    else:
        log.error("GTF missing!")
        parser.print_help(sys.stderr)
        sys.exit(1)

    n_threads = int(args.n_threads) if args.n_threads is not None else 1
    adata_out = args.adata_out if args.adata_out is not None else "adata_bulk_star.h5ad"
    fq_path = args.fq_path if args.fq_path is not None else "fastq"
    cpuCount = os.cpu_count()
    fastqs=glob(os.path.join(fq_path,"*.gz"))
    
    if len(fastqs) == 0:
        log.error("No fastq detected! (path: %s)" %fq_path)
        sys.exit(1)
        
    samples = np.unique([os.path.basename(f).split("_")[0] for f in fastqs])
    n_samples = len(samples)
    
    console.print("Pipeline parameters:", style="bold")
    console.print("fastq files path: %s" %os.path.abspath(fq_path))
    console.print("STAR reference path: %s" %star_ref)
    console.print("GTF file path: %s" %gtf)
    console.print("Output adata path: %s" %os.path.abspath(adata_out))
    
    console.print("total CPU threads: %s" %cpuCount)
    if cpuCount < n_samples*int(n_threads):
        console.print("WARNING: current settings will use more threads than available, reducing!",
                      style="bold orange3")
    while cpuCount < n_samples*n_threads:
        n_threads = n_threads-1
    console.print("Number of threads per sample: %s" %n_threads)
    console.print("With %s samples, peak utilization will use %s threads" %(n_samples,n_samples*n_threads))
    
    def countreads(fastqs):
        countfq=[fqcounter+" %s" %f for f in fastqs]
        pool = Pool()
        outs=pool.map(runcom, countfq)
        df=pd.Series(outs,index=[os.path.basename(f).split("_")[0] for f in fastqs])
        return_dict[0] = df.groupby(df.index).sum().to_dict()
    
    
    with console.status("[bold green]Counting number of reads...") as status:
        return_dict = manager.dict()
        console.log("Loading reference")
        p1 = Process(target=ref_loader,args=(star_ref,))
        p1.start()
        p2 = Process(target=countreads,args=(fastqs,))
        p2.start()
        while p1.is_alive():
            time.sleep(.2)
            p1.join()
        console.log("Reference loaded")

        p2.join()    
        
    table = Table(show_header=True, header_style="bold magenta",show_lines=True) 
    table.add_column("Samples", style="bold")
    table.add_column("Input files")
    table.add_column("Number of reads", justify="right")
    
    dct_nreads = return_dict.values()[0]
    
    for sample, nreads in dct_nreads.items():
        table.add_row(
            sample, "\n".join(glob(os.path.join(fq_path,sample+"*.gz"))), str(nreads)
        )
    
    console.print(table)
    
    console.log("Aligning reads")
    
    def run_star_par(n_threads,star_ref,fq_path,samples):
        pool = Pool()
        func = partial(run_star, n_threads,star_ref,fq_path)
        outs = pool.map(func, samples)
        return_dict["run_star"] = dict(zip(samples,outs))

    with Progress(SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeRemainingColumn(),
            console=console,
            transient=True,
            speed_estimate_period=120) as progress:
        
        dct_task = {}
        for sample,nreads in dct_nreads.items():
            dct_task[sample] = progress.add_task("[green]Aligning reads for sample %s..."%sample, 
                                            total=nreads)
        #return_dict = manager.dict()
        p3 = Process(target=run_star_par,args=(n_threads,star_ref,fq_path,samples))
        p3.start()

        #while len(glob("aligned/*.progress.out")) != len(glob(os.path.join(fq_path,"*.gz"))):
        #    time.sleep(1)
        
        while "run_star" not in return_dict:
            time.sleep(1)

        prog = True
        err = False
        errmsg=[]
        while prog:                      
            for sample,task in dct_task.items():
                if return_dict["run_star"][sample][0]!=0:
                    progress.update(task,visible=False)
                    err = True
                    errmsg.append(return_dict["run_star"][sample][1].decode('utf-8'))
                elif len(glob("aligned/*.progress.out")) != len(samples):
                    progress.update(task,completed=gather_starlogs("aligned/%s_Log.progress.out" %sample))
                    if len(glob("aligned/%s_Log.final.out" %sample))==1:
                        progress.update(task,completed=dct_nreads[sample])
            
            if err:
                time.sleep(.1)
                for e in errmsg:
                    log.error(e)
                log.info("exiting...")
                p3b = Process(target=ref_remover,args=(star_ref,))
                p3b.start()
                p3b.join()
                sys.exit(0)
            
            if len(glob("aligned/*final*")) == len(samples):
                prog = False
            progress.refresh()
            time.sleep(1)
        p3.join()
    p3b = Process(target=ref_remover,args=(star_ref,))
    p3b.start()
    console.log("Alignment done!") 

    if os.path.isdir("fc")==False:
        os.mkdir("fc")

    with console.status("[bold green]Counting features from alignments...") as status:
        p4 = Process(target=run_fc_par,args=(n_threads,gtf,))
        p4.start()
        p4.join()

    console.log("Counting done!")

    with console.status("[bold green]Merging results into one adata...") as status:
        p5 = Process(target=fc2adata_par,args=(adata_out,samples,))
        p5.start()
        p5.join()

    console.log("Merging done!")
