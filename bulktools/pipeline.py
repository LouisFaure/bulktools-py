from . import __path__, __version__
version = __version__
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
import concurrent.futures

import logging
from rich.logging import RichHandler
FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")

fqcounter = os.path.join(__path__[0], "fqcount")

manager = Manager()
return_dict = manager.dict()

help = 'Full bulk pipeline, from fastq to adata count matrix!\n'\
    'Performs the following: fastq -STAR-> bam -featureCounts-> anndata.h5ad'

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument("--fq_path","-f", help="Path for input fastq files (relative, default: fastq).")
parser.add_argument("--bam_path","-b", help="Path for aligned BAMs (default: aligned).")
parser.add_argument("--star_ref","-s", help="STAR index path.")
parser.add_argument("--gtf","-g", help="GTF file path for featureCounts.")
parser.add_argument("--n_threads","-n", help="number of threads per fastq file, for both STAR and featureCounts.")
parser.add_argument("--mem","-m", help="how much Gb to pass to --limitBAMsortRAM for STAR alignment (in Gb, default 10).")
parser.add_argument("--adata_out","-o", help="Path for the adata output (relative, default: adata_bulk_star.h5ad).")

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
   
def run_star(n_threads,mem,star_ref,fq_path,bam_path,sample):
    fqs=" ".join(glob(os.path.join(fq_path,sample+"*.gz")))
    
    runstar=f"STAR --limitBAMsortRAM {mem*1024**3} --genomeLoad LoadAndKeep --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --runThreadN {n_threads} --genomeDir {star_ref} --readFilesIn {fqs} --outFileNamePrefix {bam_path}/{sample}_"
    proc=subprocess.Popen(runstar.split(),
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    while proc.returncode == None:
        proc.poll()
        for keys,values in return_dict.items():
            if len(values)==3:
                if values[0]!=0:
                    proc.kill()
        time.sleep(1)
            
    else:   
        proc.wait()
        return_dict[sample]=[proc.returncode,proc.communicate()[1],runstar]
        
def run_star_par(n_threads,mem,star_ref,fq_path,bam_path,samples):
    pool = Pool()
    func = partial(run_star, n_threads,mem,star_ref,fq_path,bam_path)
    pool.map(func, samples)


def gather_starlogs(f):
    a_file = open(f, "r")
    lines = a_file.readlines()
    if len(lines)<=2:
        return 0
    if (len(lines)==3) & (lines[-1] == 'ALL DONE!\n'):
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

def run_fc_par(n_threads,gtf,bam_path):
    pool = Pool()
    func = partial(run_fc, n_threads,gtf)
    pool.map(func, glob(f"{bam_path}/*.bam"))
    

def fc2adata(p):
    try:
        df = pd.read_csv(p, sep='\t', comment='#',usecols=[0,6],index_col=0)
        name = p.split(".")[0]
        df.columns = [name.split("/")[-1]]
        return df
    except EmptyDataError:
        pass
    
def fc2adata_par(out,samples):
    pool = Pool()
    fcs = pool.map(fc2adata, glob("fc/*.txt"))
    allcounts = pd.concat(fcs,axis=1)
    adata = anndata.AnnData(allcounts.T)
    adata.write_h5ad(out)

def main():
    console = Console(record=True)
    console.print("bulktools %s" %version,style="bold")
    args = parser.parse_args()
    
    
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

    adata_out = args.adata_out if args.adata_out is not None else "adata_bulk_star.h5ad"
    fq_path = args.fq_path if args.fq_path is not None else "fastq"
    bam_path = args.bam_path if args.bam_path is not None else "aligned"
    mem = args.mem if args.mem is not None else 10
    cpuCount = os.cpu_count()
    fastqs=glob(os.path.join(fq_path,"*.gz"))
    
    if os.path.isdir(bam_path):
        if os.listdir(bam_path):
            log.error("Bam path is not empty! (path: %s)" %os.path.abspath(bam_path))
            sys.exit(1)
    
    if len(fastqs) == 0:
        log.error("No fastq detected! (path: %s)" %os.path.abspath(fq_path))
        sys.exit(1)
        
    samples = np.unique([os.path.basename(f).split("_")[0] for f in fastqs])
    n_samples = len(samples)
    
    n_threads = int(args.n_threads) if args.n_threads is not None else int(np.floor(cpuCount/n_samples))
    
    console.print("Pipeline parameters:", style="bold")
    console.print("fastq files path: %s" %os.path.abspath(fq_path))
    console.print("bam files path: %s" %os.path.abspath(bam_path))
    console.print("STAR reference path: %s" %os.path.abspath(star_ref))
    console.print("GTF file path: %s" %os.path.abspath(gtf))
    console.print("Output adata path: %s" %os.path.abspath(adata_out))
    
    console.print("total CPU threads: %s" %cpuCount)
    if cpuCount < n_samples*int(n_threads):
        console.print("WARNING: current settings will use more threads than available cores, reducing!",
                      style="bold orange3")

    while (cpuCount < n_samples*n_threads) & (n_threads!=1):
        n_threads = n_threads-1
    console.print("Number of threads per sample: %s" %n_threads)
    console.print("With %s samples, peak utilization will use %s threads" %(n_samples,n_samples*n_threads))
    
    def countreads(fastqs):
        countfq=[fqcounter+" %s" %f for f in fastqs]
        pool = Pool()
        outs=pool.map(runcom, countfq)
        df=pd.Series(outs,index=[os.path.basename(f).split("_")[0] for f in fastqs])
        return_dict[0] = df.groupby(df.index).sum().to_dict()
    
    
    with Progress(SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
            transient=True) as status:
        
        dct_task = {}
        dct_task["counting"] = status.add_task("[bold green]Counting number of reads..")
        p1 = Process(target=countreads,args=(fastqs,))
        p1.start()
        dct_task["loadindex"] = status.add_task("[bold green]Loading index...")
        p2 = Process(target=ref_loader,args=(star_ref,))
        p2.start()
        counting = True
        loading = True
        first = []
        while (p1.is_alive()) | (p2.is_alive()):
            if (not p1.is_alive()) & (counting):
                status.update(dct_task["counting"],visible=False)
                p1.join()
                console.log("Counting done")
                counting = False
                first.append("counting")
            if (not p2.is_alive()) & (loading):
                status.update(dct_task["loadindex"],visible=False)
                p2.join()
                console.log("Reference loaded")
                loading = False
                first.append("loading")

    console.log("Counting done" if first[0]=="loading" else "Reference loaded")
        
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
    
    allreads=list(dct_nreads.values())
    time.sleep(1)
    console.log("Aligning reads")
    
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
        
        p3 = Process(target=run_star_par,args=(n_threads,mem,star_ref,fq_path,bam_path,samples))
        p3.start()
        
        while not os.path.isdir(bam_path):
            time.sleep(1)

        prog = True
        err = False
        while prog:
            for sample,task in dct_task.items():
                if sample in return_dict:
                    if return_dict[sample][0]!=0:
                        err = True
                        errmsg=[sample,return_dict[sample][1].decode('utf-8'),return_dict[sample][2]]
                elif len(glob(f"{bam_path}/*.progress.out")) == n_samples:
                    progress.update(task,completed=gather_starlogs(f"{bam_path}/{sample}_Log.progress.out"))
                    if len(glob(f"{bam_path}/{sample}_Log.final.out"))==1:
                        progress.update(task,completed=dct_nreads[sample])
            
            if err:
                for sample,task in dct_task.items():
                    progress.update(task,visible=False)
                time.sleep(.1)
                log.error(f"Sample {errmsg[0]} experienced the following error:")
                print(errmsg[2])
                print(errmsg[1])
                time.sleep(2)
                p3.kill()
                log.info("exiting...")
                p3b = Process(target=ref_remover,args=(star_ref,))
                p3b.start()
                p3b.join()    
                time.sleep(.1)              
                sys.exit(1)
            
            if len(glob(f"{bam_path}/*final*")) == len(samples):
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
        p4 = Process(target=run_fc_par,args=(n_threads,gtf,bam_path,))
        p4.start()
        p4.join()

    console.log("Counting done!")

    with console.status("[bold green]Merging results into one adata...") as status:
        p5 = Process(target=fc2adata_par,args=(adata_out,samples,))
        p5.start()
        p5.join()

    console.log("Merging done!")
    console.save_text("bt_run.log")
