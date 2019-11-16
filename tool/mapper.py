import subprocess
import shlex
import os
import shutil


def mauve(reference, contigs, prefix, out):
    cmd = 'progressiveMauve {reference} {contigs} --output={prefix}.alignment --backbone-output={prefix}.backbone'.format(
    reference = reference, contigs = contigs, prefix = prefix, out = out)
    process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
    while process.poll() is None:
        l = process.stdout.readline()
        print(l.decode())
    print(process.stdout.read().decode())
    
    newdir = 'alignment'
    os.makedirs(os.path.join(out,newdir))
    shutil.move('{prefix}.alignment'.format(prefix = prefix), '{out}/alignment/{prefix}.alignment'.format(prefix = prefix, out = out))
    shutil.move('{prefix}.alignment.bbcols'.format(prefix = prefix), '{out}/alignment/{prefix}.alignment.bbcols'.format(out = out, prefix = prefix))
    shutil.move('{prefix}.backbone'.format(prefix = prefix), '{out}/alignment/{prefix}.backbone'.format(out = out, prefix = prefix))



