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
                                            'fMRISurface', 'DiffusionPreprocessing',
                                            'generateLevel1fsf','TaskfMRIAnalysis'],
                        default=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer',
                                 'fMRIVolume', 'fMRISurface',
                                 'DiffusionPreprocessing', 'generateLevel1fsf','TaskfMRIAnalysis'])
    parser.add_argument('-EV', dest='EV', help='Parent directory which contains folders with Psychopy output.\n'
                                               'This will generate EV files compatible for FSL task fmri analysis', required=False)

    args = parser.parse_args()

    try:
        # generate EV files
        if args.EV is not None:
            gen_onsets(args.subjID, args.EV, args.outputdir)
            # psychopy = "python /psychopy2evs.py -l " + args.EV + " -o " + args.outputdir + " -s " + args.subjID
            # subprocess.call(psychopy, shell=True)

        if not os.path.exists(args.outputdir + '/' + args.subjID + '_bids' ):
            bidsconv = "python /bidsconversion/bin/run.py " + args.sourcedir + ' ' + args.outputdir + ' /dcm2niix/build/bin/dcm2niibatch ' + '-subj ' + args.subjID + ' -dataset ' + args.dataset
            subprocess.call(bidsconv, shell=True)

        if not os.path.exists(args.outputdir+ '/' + args.subjID+'_output'):
            os.mkdir(args.outputdir+ '/' + args.subjID+'_output')

        # move previously created log files into logs directory
        logsfpath = os.path.join(args.outputdir, args.subjID+'_output', "logs")
        if not os.path.exists(logsfpath):
            os.mkdir(logsfpath)
            shutil.move(os.path.join(args.outputdir, "bids_conversion_logs"), os.path.join(logsfpath, "bids_conversion_logs"))

        runpy = "/run.py " + args.outputdir+'/'+args.subjID+'_bids/'+args.dataset + ' '+ args.outputdir+'/'+args.subjID+'_output '+" participant --n_cpus " + str(args.n_cpus) + ' ' + '--license_key ' + args.license_key + ' --stages ' + ' '.join(args.stages)
        subprocess.call(runpy, shell=True)

    except:
        print traceback.print_exc(file=sys.stdout)