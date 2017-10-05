#!/usr/bin/python
"""
Wrapper script to convert and organize DICOM raw files into BIDS-compatible files.
Also launches run.py.

"""

import os, shutil
import sys
import argparse
import traceback
import sys
import subprocess
from argparse import RawTextHelpFormatter
sys.path.append('/hcpbin')
from psychopy2evs import gen_onsets

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= "Wrapper script: runs image data conversion, organization, and image processing",
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('sourcedir', help="Path to DICOM source directory, where directories of raw image files are located."
                                          "If you already have BIDS-compatible dataset, set this positional argument to"
                                          "'skip'.")
    parser.add_argument('outputdir', help="Output directory -- must exist prior to running")
    parser.add_argument('-subjID', dest='subjID', help="Subject ID/prefix (must not contain hyphens or underscores", required=True)
    parser.add_argument('-dataset', dest='dataset', help="Data set name", required=False)
    parser.add_argument('--license_key', dest='license_key', help="FreeSurfer license key", required=True)
    parser.add_argument('--n_cpus', help='Number of CPUs/cores available to use.',
                        default=1, type=int)
    parser.add_argument('--stages', help='Which stages to run. Space separated list.',
                        nargs="+", choices=['PreFreeSurfer', 'FreeSurfer',
                                            'PostFreeSurfer', 'fMRIVolume',
                                            'fMRISurface', 'DiffusionPreprocessing', 'melodic', 'fix',
                                            'generateLevel1fsf','TaskfMRIAnalysis'],
                        default=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer',
                                 'fMRIVolume', 'fMRISurface',
                                 'DiffusionPreprocessing','melodic', 'generateLevel1fsf','TaskfMRIAnalysis'])
    parser.add_argument('-EV', dest='EV', help='Parent directory which contains folders with Psychopy output.\n'
                                               'This will generate EV files compatible for FSL task fmri analysis', required=False)
    parser.add_argument('--fslEV', help='Indicates that the EV folder path given contains FSL-compatible EV files. \n'
                                        'No conversion is required.', action='store_true', required=False)
    parser.add_argument('--bids', dest='bidsdir', help='If you have already converted your data according to BIDS specs,\n'
                                                       'then you can define the path to the converted files using this flag.',
                        required=False)

    args = parser.parse_args()

    try:
        # generate EV files
        if args.EV is not None:
            if not args.fslEV:
                gen_onsets(args.subjID, args.EV, args.outputdir)
            # psychopy = "python /psychopy2evs.py -l " + args.EV + " -o " + args.outputdir + " -s " + args.subjID
            # subprocess.call(psychopy, shell=True)

        #TODO: take out conversions?
        if args.dataset:
            dataset = " -dataset {0}".format(args.dataset)
        else:
            dataset = ""
        if not args.bidsdir and not os.path.exists(args.outputdir + '/' + args.subjID + '_bids' ):
            bidsconv = "python /bidsconversion/bin/run.py " + \
                       args.sourcedir + ' ' + \
                       args.outputdir + \
                       ' -subj ' + args.subjID + \
                       dataset
            subprocess.call(bidsconv, shell=True)

        outputFolder = args.subjID+'_output'
        if not os.path.exists(os.path.join(args.outputdir, outputFolder)):
            os.mkdir(os.path.join(args.outputdir, outputFolder))

        # move previously created log files into logs directory
        logsfpath = os.path.join(args.outputdir, outputFolder, "logs")
        try:
            if not os.path.exists(logsfpath):
                os.mkdir(logsfpath)
                shutil.move(os.path.join(args.outputdir, "bids_conversion_logs"), os.path.join(logsfpath, "bids_conversion_logs"))
        except:
            print("BIDS conversion logs not found. Not a problem! Starting processing...")

        if not args.bidsdir and args.dataset:
            bidsDataSet = os.path.join(args.outputdir, args.subjID+'_bids/'+args.dataset)
        if not args.bidsdir and not args.dataset:
            bidsDataSet = os.path.join(args.outputdir, args.subjID + '_bids/')
        else:
            bidsDataSet = args.bidsdir
        runpy = "/hcpbin/run.py " + \
                bidsDataSet + ' ' + \
                os.path.join(args.outputdir, outputFolder) + \
                " participant --n_cpus " + str(args.n_cpus) + ' ' + \
                '--license_key ' + args.license_key + \
                ' --stages ' + ' '.join(args.stages)
        subprocess.call(runpy, shell=True)

    except:
        print traceback.print_exc(file=sys.stdout)