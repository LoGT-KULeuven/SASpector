import subprocess
import shlex
import os
import shutil
import progressbar 



def mauve(reference, contigs, prefix, out):
    bar = progressbar.ProgressBar(widgets = ['Aligning: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        cmd = 'progressiveMauve {reference} {contigs} --output={prefix}.alignment --backbone-output={prefix}.backbone'.format(
        reference = reference, contigs = contigs, prefix = prefix, out = out)
        process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
        while process.poll() is None:
            l = process.stdout.readline()
        #time.sleep(0.09)
    
    newdir = 'alignment'
    os.makedirs(os.path.join(out,newdir))
    shutil.move('{prefix}.alignment'.format(prefix = prefix), '{out}/alignment/{prefix}.alignment'.format(prefix = prefix, out = out))
    shutil.move('{prefix}.alignment.bbcols'.format(prefix = prefix), '{out}/alignment/{prefix}.alignment.bbcols'.format(out = out, prefix = prefix))
    shutil.move('{prefix}.backbone'.format(prefix = prefix), '{out}/alignment/{prefix}.backbone'.format(out = out, prefix = prefix))
    shutil.move('{reference}.sslist'.format(reference = reference), '{out}/alignment/{reference}.sslist'.format(out = out, reference = reference))
    shutil.move('{contigs}.sslist'.format(contigs = contigs), '{out}/alignment/{contigs}.sslist'.format(out = out, contigs = contigs))



