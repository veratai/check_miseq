"""
Calculate the frequency of reads mapping to a library
of reference genomes from a MiSeqsample (paired FASTQ files).
"""
import sys
import os
import argparse
import subprocess
import re
from glob import glob
from datetime import datetime

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token


def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    IS_FIRST_SEGMENT = 0x40
    return (int(flag) & IS_FIRST_SEGMENT) != 0

def count_file_lines(path):
    """ Run the wc command to count lines in a file, as shown here:
    https://gist.github.com/zed/0ac760859e614cd03652
    """
    wc_output = subprocess.check_output(['wc', '-l', path])
    return int(wc_output.split()[0])


def do_map(fastq1, fastq2, refpath, bowtie_threads, min_match_len, min_mapq, min_score):
    """
    Process SAM output as it is streamed from bowtie2 to assign
    short reads to reference sequences.
    """
    # get size of FASTQ
    nrecords = count_file_lines(fastq1) / 2
    
    # stream output from bowtie2
    bowtie_args = ['bowtie2',
                   '--quiet',
                   '-x', refpath,
                   '-1', fastq1,
                   '-2', fastq2,
                   '--no-unal',
                   '--local',
                   '-p', str(bowtie_threads)]
                   
    p = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)
    
    # collect SAM output by refname
    counts = {}
    n_short = 0
    n_hybrid = 0
    n_mapq = 0
    n_low = 0
    progress = 0
    total_count = 0
    
    with p.stdout:
        for line in p.stdout:
            if line.startswith('@'):
                # skip header line
                continue

            items = line.split('\t')
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = items[:11]
            progress += 1
            if progress % 5000 == 0:
                print('[%s] (%d/%d) mapped %d (%d/%d/%d/%d)' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    progress, nrecords, 2*total_count, n_mapq, n_hybrid, n_short, n_low))
                
            taxon = rname.replace('-', '_').split('_')[0]  # 'HCV-1a' -> 'HCV_1a' -> 'HCV'
            
            # ignore second reads
            if not is_first_read(flag):
                continue
            
            # discard reads with low map quality
            q = int(mapq)
            if q < min_mapq:
                n_mapq += 1
                continue

            # ignore reads whose mate mapped to a different taxon
            if rnext != '=' and taxon != rnext.split('-')[0]:
                n_hybrid += 1
                continue
            
            # filter out reads based on match length - e.g., primers
            tokens = cigar_re.findall(cigar)
            match_len = sum([int(token.strip('M')) for token in tokens if token.endswith('M')])
            if match_len < min_match_len:
                n_short += 1
                continue

            # filter out reads with low alignment score
            optionals = dict([item.split(':i:') for item in items[11:] if ':i:' in item])
            if 'AS' not in optionals or int(optionals['AS']) < min_score:
                n_low += 1
                continue

            # update counts
            if taxon not in counts:
                counts.update({taxon: 0})
            counts[taxon] += 1
            total_count += 1
            
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, bowtie_args)

    # output results
    keys = counts.keys()
    keys.sort()

    return counts, n_short, n_mapq, n_hybrid



def main():
    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')
    
    parser.add_argument('-path', help='<input> folder containing FASTQ files')
    parser.add_argument('-pathlist', help='<input> file containing list of folders')
    parser.add_argument('-output', help='<output CSV> file to write results')
    parser.add_argument('-log', help='path to logfile', default='checkMiSeq.batch.log')
    parser.add_argument('-x', help='path to bowtie2 index (.bt2)', default='/home/veratai/RefSeqs/checkMiSeq_refseqs')
    parser.add_argument('-p', type=int, help='number of bowtie2 threads', default=6)
    parser.add_argument('-minlen', type=int, help='minimum match length (CIGAR M)', default=100)
    parser.add_argument('-minq', type=int, help='minimum mapping quality (MAPQ)', default=0)
    parser.add_argument('-mins', type=int, help='minimum alignment score', default=0)

    args = parser.parse_args()
    
    paths = []
    if args.pathlist:
        handle = open(args.pathlist, 'rU')
        for line in handle:
            paths.append(line.strip('\n'))
        handle.close()
    elif args.path:
        paths.append(args.path)
    else:
        parser.print_help()
        sys.exit()
    
    # check that we have access to bowtie2
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')
    
    
    
    complete = {}
    if os.path.exists(args.output):
        handle = open(args.output, 'rU')
        header = handle.next()
        for line in handle:
            runname, sample, snum = line.strip('\n').split(',')[:3]
            complete.update({(runname, sample, snum): 1})
        handle.close()
        # re-open file for appending
        handle = open(args.output, 'a')
    else:
        # new output file
        handle = open(args.output, 'w')
        handle.write('runname,sample,snum,taxon,count,total\n')
    
    for path in paths:
        # check that the inputs exist
        if not os.path.exists(path):
            print('No folder found at %s', path)
            sys.exit(1)

        if not path.endswith('/'):
            path += '/'

        files = glob(path + '*_R1_001.fastq')
        if len(files) == 0:
            print 'ERROR: No FASTQ R1 files found at', path
            sys.exit()
   
        runname = path.split('/')[4]

        for f1 in files:
            f2 = f1.replace('_R1_', '_R2_')
            filename = os.path.basename(f1)
            if filename.startswith('Undetermined'):
                continue
            
            sample, snum = filename.split('_')[:2]
            
            if (runname, sample, snum) in complete:
                # skip previously processed sample
                continue
                
            logfile = open(args.log, 'a')
            logfile.write('[%s] start processing %s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), filename))
            logfile.close()

            counts, n_short, n_mapq, n_hybrid = do_map(f1, f2, args.x, args.p, args.minlen, args.minq, args.mins)
            n_discard = n_short + n_mapq + n_hybrid
            total_count = sum(counts.values()) + n_discard

            for taxon, count in counts.iteritems():
                handle.write('%s,%s,%s,%s,%d,%d\n' % (runname, sample, snum, taxon, count, total_count))

            # record number of reads that failed to map
            handle.write('%s,%s,%s,not_referenced,%d,%d\n' % (runname, sample, snum, n_discard, total_count))
            handle.flush()  # clear write buffer
            
            logfile = open(args.log, 'a')
            logfile.write('[%s] end processing %s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), filename))
            logfile.close()
            
    handle.close()
    logfile.close()


if __name__ == '__main__':
    main()
