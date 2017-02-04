#!/usr/bin/python

import os
import argparse
import traceback
import sys
import subprocess
from argparse import RawTextHelpFormatter

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= "Wrapper script: runs image data conversion, organization, and image processing",
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('sourcedir', help="Path to DICOM source directory, where directories of raw image files are located")
    parser.add_argument('outputdir', help="Output directory -- must exist prior to running")
    parser.add_argument('-subjID', dest='subjID', help="Subject ID/prefix (must not contain hyphens or underscores", required=True)
    parser.add_argument('-dataset', dest='dataset', help="Data set name", required=True)
    parser.add_argument('--license_key', dest='license_key', help="FreeSurfer license key", required=True)
    parser.add_argument('--n_cpus', help='Number of CPUs/cores available to use.',
                        default=1, type=int)
    parser.add_argument('--stages', help='Which stages to run. Space separated list.',
                        nargs="+", choices=['PreFreeSurfer', 'FreeSurfer',
                                            'PostFreeSurfer', 'fMRIVolume',
                                            'fMRISurface', 'DiffusionPreprocessing','TaskfMRIAnalysis'],
                        default=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer',
                                 'fMRIVolume', 'fMRISurface',
                                 'DiffusionPreprocessing', 'TaskfMRIAnalysis'])

    args = parser.parse_args()

    try:
        bidsconv = "python /bidsconversion/bin/run.py " + args.sourcedir + ' ' + args.outputdir + ' /dcm2niix/bin/dcm2niibatch ' + '-subj ' + args.subjID + ' -dataset ' + args.dataset
        subprocess.call(bidsconv, shell=True)

        os.mkdir(args.outputdir+ '/' + args.subjID+'_output')
        runpy = "/run.py " + args.outputdir+'/'+args.subjID+'_bids/'+args.dataset + ' '+ args.outputdir+'/'+args.subjID+'_output '+" participant --n_cpus " + str(args.n_cpus) + ' ' + '--license_key ' + args.license_key + ' --stages ' + ' '.join(args.stages)
        subprocess.call(runpy, shell=True)

    except:
        print traceback.print_exc(file=sys.stdout)