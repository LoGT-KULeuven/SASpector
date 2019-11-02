from __future__ import print_function
import subprocess
import shlex
import argparse

def mauve(reference, contigs, prefix):
    cmd = 'progressiveMauve {reference} {contigs} --output={prefix}.alignment --backbone-output={prefix}.backbone'.format(
    reference = reference, contigs = contigs, prefix = prefix)
    process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip().decode())
        rc = process.poll()
    return rc


def main():
    parser = argparse.ArgumentParser(prog = 'What the *** is wrong with my Illumina Assembly? -- Mapper', description = 'Genome mapper with progressiveMauve for Hybrid and De novo assemblies')
    parser.add_argument('reference', type = str, help = 'Hybrid assembly FASTA file as reference genome')
    parser.add_argument('contigs', type = str, help = 'De novo assembly FASTA file as contigs/draft genome')
    parser.add_argument('prefix', type = str, help = 'Genome ID')
    args = parser.parse_args()

    mauve(args.reference, args.contigs, args.prefix)

if __name__ == '__main__':
    main()
