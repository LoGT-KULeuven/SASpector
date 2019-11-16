from Bio.Blast.Applications import NcbiblastnCommandline as blast
import argparse

def blast2Unmap(args):
    cline = blast(cmd = 'blastn', out = '{prefix}_blast2Unmap.tsv'.format(prefix = args.prefix), outfmt = '6 qseqid qlen qstart qend sseqid sstart send pident evalue qcovs sstrand',
    task = 'megablast', query = args.unmapped, subject = args.contigs)
    stdout, stderr = cline()

def blast2Conf(args):
    cline = blast(cmd = 'blastn', out = '{prefix}_blast2Conf.tsv'.format(prefix = args.prefix), outfmt = '6 qseqid qlen qstart qend sseqid sstart send pident evalue qcovs sstrand',
    task = 'megablast', query = args.conflict, subject = args.reference)
    stdout, stderr = cline()

def main():
    parser = argparse.ArgumentParser(prog = 'SASpector - Checker', description = 'Multiple alignment of unmapped regions and conflict regions to detect false positives in the genome mapping')
    subparsers = parser.add_subparsers()

    parser_blast2Unmap = subparsers.add_parser('blast2Unmap', help = 'Alignment of unmapped regions', description = 'BLAST for the unmapped regions with the contigs of the short-read assembly')
    parser_blast2Unmap.add_argument('unmapped', type = str)
    parser_blast2Unmap.add_argument('contigs', type = str)
    parser_blast2Unmap.set_defaults(function = blast2Unmap)

    parser_blast2Conf = subparsers.add_parser('blast2Conf', help = 'Alignment of conflict contigs', description = 'BLAST for the conflict regions with the reference genome')
    parser_blast2Conf.add_argument('conflict', type = str)
    parser_blast2Conf.add_argument('reference', type = str)
    parser_blast2Conf.set_defaults(function = blast2Conf)

    parser.add_argument('prefix', type = str)

    args = parser.parse_args()
    args.function(args)

if __name__ == '__main__':
    main()
