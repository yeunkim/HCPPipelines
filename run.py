#!/usr/bin/python
"""
Calls HCP scripts to run the HCP Pipelines.
"""

from __future__ import print_function
import argparse
import os
import shutil
import nibabel
from glob import glob
from subprocess import Popen, PIPE, check_output
from shutil import rmtree
import subprocess
from bids.grabbids import BIDSLayout
from functools import partial
from collections import OrderedDict
import time
from multiprocessing import Process, Pool
import logging
import datetime

ms = lambda: int(round(time.time() * 1000))

def run(command, env={}, cwd=None, stage='', filename='', subject=''):
    merged_env = os.environ
    merged_env.update(env)
    merged_env.pop("DEBUG", None)
    # suffix = ms()
    # logfn = stage + '_' + str(suffix) + '.log'
    logfn = subject + '_' + stage + filename + '.log'
    logpath = os.path.join(str(cwd),'logs', logfn)
    logfile = open(logpath, 'w')
    process = Popen(command, stdout=PIPE, stderr=subprocess.STDOUT,
                    shell=True, env=merged_env, cwd=cwd)

    for line in process.stdout:
        logfile.write(line)

    while True:
        line = process.stdout.readline()
        line = str(line)[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)

grayordinatesres = "2" # This is currently the only option for which the is an atlas
lowresmesh = 32

def run_pre_freesurfer(**args):
    args.update(os.environ)
    args["t1"] = "@".join(t1ws)
    args["t2"] = "@".join(t2ws)

    cmd = '{HCPPIPEDIR}/PreFreeSurfer/PreFreeSurferPipeline.sh ' + \
    '--path="{path}" ' + \
    '--subject="{subject}" ' + \
    '--t1="{t1}" ' + \
    '--t2="{t2}" ' + \
    '--t1template="{HCPPIPEDIR_Templates}/MNI152_T1_{t1_template_res:.1f}mm.nii.gz" ' + \
    '--t1templatebrain="{HCPPIPEDIR_Templates}/MNI152_T1_{t1_template_res:.1f}mm_brain.nii.gz" ' + \
    '--t1template2mm="{HCPPIPEDIR_Templates}/MNI152_T1_2mm.nii.gz" ' + \
    '--t2template="{HCPPIPEDIR_Templates}/MNI152_T2_{t2_template_res:.1f}mm.nii.gz" ' + \
    '--t2templatebrain="{HCPPIPEDIR_Templates}/MNI152_T2_{t2_template_res:.1f}mm_brain.nii.gz" ' + \
    '--t2template2mm="{HCPPIPEDIR_Templates}/MNI152_T2_2mm.nii.gz" ' + \
    '--templatemask="{HCPPIPEDIR_Templates}/MNI152_T1_{t1_template_res:.1f}mm_brain_mask.nii.gz" ' + \
    '--template2mmmask="{HCPPIPEDIR_Templates}/MNI152_T1_2mm_brain_mask_dil.nii.gz" ' + \
    '--brainsize="150" ' + \
    '--fnirtconfig="{HCPPIPEDIR_Config}/T1_2_MNI152_2mm.cnf" ' + \
    '--fmapmag="{fmapmag}" ' + \
    '--fmapphase="{fmapphase}" ' + \
    '--fmapgeneralelectric="NONE" ' + \
    '--echodiff="{echodiff}" ' + \
    '--SEPhaseNeg="{SEPhaseNeg}" ' + \
    '--SEPhasePos="{SEPhasePos}" ' + \
    '--echospacing="{echospacing}" ' + \
    '--seunwarpdir="{seunwarpdir}" ' + \
    '--t1samplespacing="{t1samplespacing}" ' + \
    '--t2samplespacing="{t2samplespacing}" ' + \
    '--unwarpdir="{unwarpdir}" ' + \
    '--gdcoeffs="NONE" ' + \
    '--avgrdcmethod={avgrdcmethod} ' + \
    '--topupconfig="{HCPPIPEDIR_Config}/b02b0.cnf" ' + \
    '--printcom=""'
    cmd = cmd.format(**args)
    logging.info(" {0} : Running PreFreeSurfer".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    t = time.time()
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])}, stage="PreFreeSurfer", subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished running PreFreeSurfer. Time duration: {0} minutes".format(str(elapsed)))
    os.sys.stdout.write("\nElapsed time for PreFreeSurfer is " + str(elapsed) + " minutes. \n")

def run_freesurfer(**args):
    args.update(os.environ)
    args["subjectDIR"] = os.path.join(args["path"], args["subject"], "T1w")
    cmd = '{HCPPIPEDIR}/FreeSurfer/FreeSurferPipeline.sh ' + \
      '--subject="{subject}" ' + \
      '--subjectDIR="{subjectDIR}" ' + \
      '--t1="{path}/{subject}/T1w/T1w_acpc_dc_restore.nii.gz" ' + \
      '--t1brain="{path}/{subject}/T1w/T1w_acpc_dc_restore_brain.nii.gz" ' + \
      '--t2="{path}/{subject}/T1w/T2w_acpc_dc_restore.nii.gz" ' + \
      '--printcom=""'
    cmd = cmd.format(**args)

    if not os.path.exists(os.path.join(args["subjectDIR"], "fsaverage")):
        shutil.copytree(os.path.join(os.environ["SUBJECTS_DIR"], "fsaverage"),
                        os.path.join(args["subjectDIR"], "fsaverage"))
    if not os.path.exists(os.path.join(args["subjectDIR"], "lh.EC_average")):
        shutil.copytree(os.path.join(os.environ["SUBJECTS_DIR"], "lh.EC_average"),
                        os.path.join(args["subjectDIR"], "lh.EC_average"))
    if not os.path.exists(os.path.join(args["subjectDIR"], "rh.EC_average")):
        shutil.copytree(os.path.join(os.environ["SUBJECTS_DIR"], "rh.EC_average"),
                        os.path.join(args["subjectDIR"], "rh.EC_average"))
    logging.info(" {0} : Running FreeSurfer".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    t = time.time()
    run(cmd, cwd=args["path"], env={"NSLOTS": str(args["n_cpus"]),
                                    "OMP_NUM_THREADS": str(args["n_cpus"])}, stage="FreeSurfer",subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished running FreeSurfer. Time duration: {0} minutes".format(str(elapsed)))
    os.sys.stdout.write("\nElapsed time for FreeSurfer is " + str(elapsed) + " minutes. \n")

def run_post_freesurfer(**args):
    args.update(os.environ)
    cmd = '{HCPPIPEDIR}/PostFreeSurfer/PostFreeSurferPipeline.sh ' + \
      '--path="{path}" ' + \
      '--subject="{subject}" ' + \
      '--surfatlasdir="{HCPPIPEDIR_Templates}/standard_mesh_atlases" ' + \
      '--grayordinatesdir="{HCPPIPEDIR_Templates}/91282_Greyordinates" ' + \
      '--grayordinatesres="{grayordinatesres:s}" ' + \
      '--hiresmesh="164" ' + \
      '--lowresmesh="{lowresmesh:d}" ' + \
      '--subcortgraylabels="{HCPPIPEDIR_Config}/FreeSurferSubcorticalLabelTableLut.txt" ' + \
      '--freesurferlabels="{HCPPIPEDIR_Config}/FreeSurferAllLut.txt" ' + \
      '--refmyelinmaps="{HCPPIPEDIR_Templates}/standard_mesh_atlases/Conte69.MyelinMap_BC.164k_fs_LR.dscalar.nii" ' + \
      '--regname="FS" ' + \
      '--printcom=""'
    cmd = cmd.format(**args)
    logging.info(" {0} : Running PostFreeSurfer".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    t = time.time()
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])}, stage="PostFreeSurfer",subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished running PostFreeSurfer. Time duration: {0} minutes".format(str(elapsed)))
    os.sys.stdout.write("\nElapsed time for PostFreeSurfer is " + str(elapsed) + " minutes. \n")

def run_generic_fMRI_volume_processsing(**args):
    args.update(os.environ)
    cmd = '{HCPPIPEDIR}/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh ' + \
      '--path={path} ' + \
      '--subject={subject} ' + \
      '--fmriname={fmriname} ' + \
      '--fmritcs={fmritcs} ' + \
      '--fmriscout={fmriscout} ' + \
      '--SEPhaseNeg={SEPhaseNeg} ' + \
      '--SEPhasePos={SEPhasePos} ' + \
      '--fmapmag="NONE" ' + \
      '--fmapphase="NONE" ' + \
      '--fmapgeneralelectric="NONE" ' + \
      '--echospacing={echospacing} ' + \
      '--echodiff="NONE" ' + \
      '--unwarpdir={unwarpdir} ' + \
      '--fmrires={fmrires:s} ' + \
      '--dcmethod={dcmethod} ' + \
      '--gdcoeffs="NONE" ' + \
      '--topupconfig={HCPPIPEDIR_Config}/b02b0.cnf ' + \
      '--printcom="" ' + \
      '--biascorrection={biascorrection} ' + \
      '--mctype="MCFLIRT"'
    cmd = cmd.format(**args)
    logging.info(" {0} : Running fMRIVolume".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    t = time.time()
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])}, stage="fMRIVolume", filename='_{0}'.format(args["fmriname"]),subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished running fMRIVolume. Time duration: {0} minutes".format(str(elapsed)))
    os.sys.stdout.write("\nElapsed time for fMRIVolume is " + str(elapsed) + " minutes. \n")

def run_generic_fMRI_surface_processsing(**args):
    # print(args)
    args.update(os.environ)
    cmd = '{HCPPIPEDIR}/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh ' + \
      '--path={path} ' + \
      '--subject={subject} ' + \
      '--fmriname={fmriname} ' + \
      '--lowresmesh="{lowresmesh:d}" ' + \
      '--fmrires={fmrires:s} ' + \
      '--smoothingFWHM={fmrires:s} ' + \
      '--grayordinatesres="{grayordinatesres:s}" ' + \
      '--regname="FS"'
    cmd = cmd.format(**args)
    logging.info(" {0} : Running fMRISurface".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    t = time.time()
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])}, stage="fMRISurface", filename='_{0}'.format(args["fmriname"]),subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished running fMRISurface. Time duration: {0} minutes".format(str(elapsed)))
    os.sys.stdout.write("\nElapsed time for fMRISurface is " + str(elapsed) + " minutes. \n")

def run_diffusion_processsing(**args):
    # print(args)
    args.update(os.environ)
    cmd = '{HCPPIPEDIR}/DiffusionPreprocessing/DiffPreprocPipeline.sh ' + \
      '--posData="{posData}" ' +\
      '--negData="{negData}" ' + \
      '--path="{path}" ' + \
      '--subject="{subject}" ' + \
      '--echospacing="{echospacing}" '+ \
      '--PEdir={PEdir} ' + \
      '--gdcoeffs="NONE" ' + \
      '--dwiname="{dwiname}" ' + \
      '--printcom=""'
    cmd = cmd.format(**args)
    t = time.time()
    logging.info(" {0} : Running DiffusionPreprocessing".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])}, stage='DiffusionPreprocessing', filename='_{0}'.format(args["dwiname"]),subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished running DiffusionPreprocessing. Time duration: {0} minutes".format(str(elapsed)))
    os.sys.stdout.write("\nElapsed time for DiffusionPreprocessing is " + str(elapsed) + " minutes. \n")

def generate_level1_fsf(**args):
    # print(args)
    args.update(os.environ)
    cmd = '{HCPPIPEDIR}/Examples/Scripts/generate_level_1_fsf_dev.sh ' + \
        '--studyfolder={studyfolder} ' + \
        '--subject={subject} ' + \
        '--taskname={taskname} ' + \
        '--templatedir={HCPPIPEDIR}/{templatedir} ' + \
        '--outdir={outdir} ' + \
        '--dir={dir}'
    cmd = cmd.format(**args)
    t = time.time()
    logging.info(" {0} : Generating level 1 FSF file".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    # print('\n', cmd, '\n')
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])},subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished generating level 1 FSF file. Time duration: {0} minutes".format(str(elapsed)))

def generate_level2_fsf(**args):
    # print(args)
    args.update(os.environ)
    cmd = '{HCPPIPEDIR}/Examples/Scripts/generate_level_2_fsf_dev.sh ' + \
        '--studyfolder={studyfolder} ' + \
        '--subject={subject} ' + \
        '--taskname={taskname} ' + \
        '--templatedir={HCPPIPEDIR}/{templatedir} ' + \
        '--outdir={outdir} '
    cmd = cmd.format(**args)
    t = time.time()
    logging.info(" {0} : Generating level 2 FSF file".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    # print('\n', cmd, '\n')
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])},subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished generating level 2 FSF file. Time duration: {0} minutes".format(str(elapsed)))

def run_task_fmri_analysis(**args):
    print(args)
    args.update(os.environ)
    cmd = '{HCPPIPEDIR}/TaskfMRIAnalysis/TaskfMRIAnalysis.sh ' + \
        '--path={path} ' + \
        '--subject={subject} ' + \
        '--lvl1tasks={lvl1tasks} ' + \
        '--lvl1fsfs={lvl1fsfs} ' + \
        '--lvl2task={lvl2task} ' + \
        '--lvl2fsf={lvl2fsf} ' + \
        '--lowresmesh={lowresmesh} ' + \
        '--grayordinatesres="{grayordinatesres:s}" ' + \
        '--confound={confound} ' + \
        '--finalsmoothingFWHM={finalsmoothingFWHM} ' + \
        '--temporalfilter={temporalfilter} ' + \
        '--vba={vba} ' + \
        '--regname={regname} ' + \
        '--parcellation={parcellation} ' + \
        '--parcellationfile={parcellationfile} ' + \
        '--printcom=""'
    cmd = cmd.format(**args)
    t = time.time()
    logging.info(" {0} : Running TaskfMRIAnalysis".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
    logging.info(cmd)
    # print('\n', cmd, '\n')
    run(cmd, cwd=args["path"], env={"OMP_NUM_THREADS": str(args["n_cpus"])}, stage='TaskfMRIAnalysis', filename='_{0}'.format(args["lvl1tasks"]),subject=args["subject"])
    elapsed = time.time() - t
    elapsed = elapsed / 60
    logging.info("Finished running TaskfMRIAnalysis. Time duration: {0} minutes".format(str(elapsed)))

def func_stages(stages_dict):
    for stage, stage_func in stages_dict.iteritems():
        if stage in args.stages:
            stage_func()

__version__ = open('/version').read()

parser = argparse.ArgumentParser(description='HCP Pipelines BIDS App (T1w, T2w, fMRI)')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
                    'should be stored. If you are running group level analysis '
                    'this folder should be prepopulated with the results of the'
                    'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                    'Multiple participant level analyses can be run independently '
                    '(in parallel) using the same output_dir.',
                    choices=['participant'])
parser.add_argument('--participant_label', help='The label of the participant that should be analyzed. The label '
                   'corresponds to sub-<participant_label> from the BIDS spec '
                   '(so it does not include "sub-"). If this parameter is not '
                   'provided all subjects should be analyzed. Multiple '
                   'participants can be specified with a space separated list.',
                   nargs="+")
parser.add_argument('--n_cpus', help='Number of CPUs/cores available to use.',
                   default=1, type=int)
parser.add_argument('--stages', help='Which stages to run. Space separated list.',
                   nargs="+", choices=['PreFreeSurfer', 'FreeSurfer',
                                       'PostFreeSurfer', 'fMRIVolume',
                                       'fMRISurface', 'DiffusionPreprocessing', 'TaskfMRIAnalysis'],
                   default=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer',
                            'fMRIVolume', 'fMRISurface',
                            'DiffusionPreprocessing', 'TaskfMRIAnalysis'])
parser.add_argument('--license_key', help='FreeSurfer license key - letters and numbers after "*" in the email you received after registration. To register (for free) visit https://surfer.nmr.mgh.harvard.edu/registration.html',
                    required=True)
parser.add_argument('-v', '--version', action='version',
                    version='HCP Pielines BIDS App version {}'.format(__version__))

args = parser.parse_args()

starttime = time.time()
# suffix = ms()
## make logs directory
if not os.path.exists(os.path.join(args.output_dir, "logs")):
    os.makedirs(os.path.join(args.output_dir, "logs"))

subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
logging.basicConfig(filename=os.path.join(args.output_dir, 'logs','sub-{0}_HCPPipelines.log'.format(subjects_to_analyze[0])),
                    level=logging.DEBUG, format='%(levelname)s:%(message)s', filemode='w')
logging.info(" {0} : Starting the HCP processing pipeline...".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
# subprocess.call("python genStatusFile.py {0} /public_html".format(subjects_to_analyze[0]), shell=True)

logging.info(" {0} : Running bids-validator".format(datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Z %Y")))
#
#
t = time.time()
run("bids-validator " + args.bids_dir, cwd=args.output_dir, stage="bids-validator", subject='sub-{0}'.format(subjects_to_analyze[0]))
elapsed = time.time() - t
elapsed = elapsed / 60
logging.info("Finished running bids-validator. Time duration: {0} minutes".format(str(elapsed)))

layout = BIDSLayout(args.bids_dir)
subjects_to_analyze = []
# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# running participant level
if args.analysis_level == "participant":
    # find all T1s and skullstrip them
    for subject_label in subjects_to_analyze:
        logging.info("Running subject: {0}".format(subjects_to_analyze[0]))
        boldsnumruns = 0
        dwinumruns = 0
        dirnums = []

        numruns = set(layout.get(target='run', return_type='id',
                                 subject=subject_label, type='dwi',
                                 extensions=["nii.gz", "nii"]))

        for session in numruns:
            acqs = set(layout.get(target='acquisition', return_type='id',
                                  subject=subject_label, type='dwi',
                                  extensions=["nii.gz", "nii"]))
            for acq in acqs:
                y = [int(s) for s in acq[0:len(acq)] if s.isdigit()]
                y = [str(s) for s in y]
                y = ''.join(y)
                dirnums.append(y)

        dwinumruns = len(numruns) * len(set(dirnums))
        logging.info("Total number of DiffusionPreprocessing processes: {0}".format(str(dwinumruns)))

        bolds = [f.filename for f in layout.get(subject=subject_label,
                                                type='bold',
                                                extensions=["nii.gz", "nii"])]

        boldsnumruns = len(bolds)
        logging.info("Total number of fMRI processes: {0}".format(str(boldsnumruns)))



        t1ws = [f.filename for f in layout.get(subject=subject_label,
                                               type='T1w',
                                               extensions=["nii.gz", "nii"])]
        t2ws = [f.filename for f in layout.get(subject=subject_label,
                                               type='T2w',
                                               extensions=["nii.gz", "nii"])]
        assert (len(t1ws) > 0), "No T1w files found for subject %s!"%subject_label
        assert (len(t2ws) > 0), "No T2w files found for subject %s!"%subject_label

        available_resolutions = [0.7, 0.8, 1.0]
        t1_zooms = nibabel.load(t1ws[0]).get_header().get_zooms()
        t1_res = float(min(t1_zooms[:3]))
        t1_template_res = min(available_resolutions, key=lambda x:abs(x-t1_res))
        t2_zooms = nibabel.load(t2ws[0]).get_header().get_zooms()
        t2_res = float(min(t2_zooms[:3]))
        t2_template_res = min(available_resolutions, key=lambda x:abs(x-t2_res))

        fieldmap_set = layout.get_fieldmap(t1ws[0])
        fmap_args = {"fmapmag": "NONE",
                     "fmapphase": "NONE",
                     "echodiff": "NONE",
                     "t1samplespacing": "NONE",
                     "t2samplespacing": "NONE",
                     "unwarpdir": "NONE",
                     "avgrdcmethod": "NONE",
                     "SEPhaseNeg": "NONE",
                     "SEPhasePos": "NONE",
                     "echospacing": "NONE",
                     "seunwarpdir": "NONE"}

        if fieldmap_set:
            t1_spacing = layout.get_metadata(t1ws[0])["RealDwellTime"]
            t2_spacing = layout.get_metadata(t2ws[0])["RealDwellTime"]

            unwarpdir = layout.get_metadata(t1ws[0])["PhaseEncodingDirection"]
            unwarpdir = unwarpdir.replace("i","x").replace("j", "y").replace("k", "z")
            if len(unwarpdir) == 2:
                unwarpdir = "-" + unwarpdir[0]

            fmap_args.update({"t1samplespacing": "%.8f"%t1_spacing,
                              "t2samplespacing": "%.8f"%t2_spacing,
                              "unwarpdir": unwarpdir})

            if fieldmap_set["type"] == "phasediff":
                merged_file = "%s/tmp/%s/magfile.nii.gz"%(args.output_dir, subject_label)
                run("mkdir -p %s/tmp/%s/ && fslmerge -t %s %s %s"%(args.output_dir,
                subject_label,
                merged_file,
                fieldmap_set["magnitude1"],
                fieldmap_set["magnitude2"]))

                phasediff_metadata = layout.get_metadata(fieldmap_set["phasediff"])
                te_diff = phasediff_metadata["EchoTime2"] - phasediff_metadata["EchoTime1"]
                # HCP expects TE in miliseconds
                te_diff = te_diff*1000.0

                fmap_args.update({"fmapmag": merged_file,
                                  "fmapphase": fieldmap_set["phasediff"],
                                  "echodiff": "%.6f"%te_diff,
                                  "avgrdcmethod": "SiemensFieldMap"})
            elif fieldmap_set["type"] == "epi":
                SEPhaseNeg = None
                SEPhasePos = None
                for fieldmap in fieldmap_set["epi"]:
                    enc_dir = layout.get_metadata(fieldmap)["PhaseEncodingDirection"]
                    if "-" in enc_dir:
                        SEPhaseNeg = fieldmap
                    else:
                        SEPhasePos = fieldmap

                seunwarpdir = layout.get_metadata(fieldmap_set["epi"][0])["PhaseEncodingDirection"]
                seunwarpdir = seunwarpdir.replace("-", "").replace("i","x").replace("j", "y").replace("k", "z")

                #TODO check consistency of echo spacing instead of assuming it's all the same
                if "EffectiveEchoSpacing" in layout.get_metadata(fieldmap_set["epi"][0]):
                    echospacing = layout.get_metadata(fieldmap_set["epi"][0])["EffectiveEchoSpacing"]
                elif "TotalReadoutTime" in layout.get_metadata(fieldmap_set["epi"][0]):
                    # HCP Pipelines do not allow users to specify total readout time directly
                    # Hence we need to reverse the calculations to provide echo spacing that would
                    # result in the right total read out total read out time
                    # see https://github.com/Washington-University/Pipelines/blob/master/global/scripts/TopupPreprocessingAll.sh#L202
                    print("BIDS App wrapper: Did not find EffectiveEchoSpacing, calculating it from TotalReadoutTime")
                    # TotalReadoutTime = EffectiveEchoSpacing * (len(PhaseEncodingDirection) - 1)
                    total_readout_time = layout.get_metadata(fieldmap_set["epi"][0])["TotalReadoutTime"]
                    phase_len = nibabel.load(fieldmap_set["epi"][0]).shape[{"x": 0, "y": 1}[seunwarpdir]]
                    echospacing = TotalReadoutTime / float(phase_len - 1)
                else:
                    raise RuntimeError("EffectiveEchoSpacing or TotalReadoutTime defined for the fieldmap intended for T1w image. Please fix your BIDS dataset.")

                fmap_args.update({"SEPhaseNeg": SEPhaseNeg,
                                  "SEPhasePos": SEPhasePos,
                                  "echospacing": "%.6f"%echospacing,
                                  "seunwarpdir": seunwarpdir,
                                  "avgrdcmethod": "TOPUP"})
        #TODO add support for GE fieldmaps

        struct_stages_dict = OrderedDict([("PreFreeSurfer", partial(run_pre_freesurfer,
                                                path=args.output_dir,
                                                subject="sub-%s"%subject_label,
                                                t1ws=t1ws,
                                                t2ws=t2ws,
                                                n_cpus=args.n_cpus,
                                                t1_template_res=t1_template_res,
                                                t2_template_res=t2_template_res,
                                                **fmap_args)),
                       ("FreeSurfer", partial(run_freesurfer,
                                             path=args.output_dir,
                                             subject="sub-%s"%subject_label,
                                             n_cpus=args.n_cpus)),
                       ("PostFreeSurfer", partial(run_post_freesurfer,
                                                 path=args.output_dir,
                                                 subject="sub-%s"%subject_label,
                                                 grayordinatesres=grayordinatesres,
                                                 lowresmesh=lowresmesh,
                                                 n_cpus=args.n_cpus))
                       ])
        for stage, stage_func in struct_stages_dict.iteritems():
            if stage in args.stages:
                stage_func()

        # dwis = layout.get(subject=subject_label, type='dwi',
        #                                          extensions=["nii.gz", "nii"])

        # print(dwis)

        posData = []
        negData = []
        PEdir = "None"
        dwiname = "Diffusion"
        dirnums = []


        #TODO: remove redundant code after testing
        numruns = set(layout.get(target='run', return_type='id',
                                 subject=subject_label, type='dwi',
                                 extensions=["nii.gz", "nii"]))

        for session in numruns:
            acqs = set(layout.get(target='acquisition', return_type='id',
                                  subject=subject_label, type='dwi',
                                  extensions=["nii.gz", "nii"]))
            for acq in acqs:
                y = [int(s) for s in acq[0:len(acq)] if s.isdigit()]
                y = [str(s) for s in y]
                y = ''.join(y)
                dirnums.append(y)

            dwinumruns = len(numruns) * len(dirnums)

            for dirnum in set(dirnums):
                dwiname = "Diffusion" + "_dir-" + dirnum + "_" + session + "_corr"

                diracqs = [x for x in acqs if dirnum in x]
                if "AP" or "PA" in diracqs:
                    PEdir = 2
                elif "LR" or "RL" in diracqs:
                    PEdir = 1
                else:
                    RuntimeError("Acquisition direction not specified on dwi file")
                pos = "EMPTY"
                neg = "EMPTY"
                gdcoeffs = "None"

                for d in set(diracqs):
                    dwis = layout.get(subject=subject_label,
                                      type='dwi', acquisition=d, run=session,
                                      extensions=["nii.gz", "nii"])
                    for dwi in dwis:
                        dwi = dwi.filename
                        if "-" in layout.get_metadata(dwi)["PhaseEncodingDirection"]:
                            neg = dwi
                            # negData.append(neg)
                        else:
                            pos = dwi
                            # posData.append(pos)

                # assert len(dwis) <= 2
                # for dwi in dwis:
                #     dwi = dwi.filename
                #     if "-" in layout.get_metadata(dwi)["PhaseEncodingDirection"]:
                #         neg = dwi
                #         # negData.append(neg)
                #     else:
                #         pos = dwi
                #         # posData.append(pos)

                echospacing = layout.get_metadata(pos)["EffectiveEchoSpacing"] * 1000
                dwi_stage_dict = OrderedDict([("DiffusionPreprocessing", partial(run_diffusion_processsing,
                                                                                 posData=pos,
                                                                                 negData=neg,
                                                                                 path=args.output_dir,
                                                                                 subject="sub-%s" % subject_label,
                                                                                 echospacing=echospacing,
                                                                                 PEdir=PEdir,
                                                                                 gdcoeffs="NONE",
                                                                                 dwiname=dwiname,
                                                                                 n_cpus=args.n_cpus))])
                for stage, stage_func in dwi_stage_dict.iteritems():
                    if stage in args.stages:
                        # stage_func()
                        try:
                            Process(target=stage_func).start()
                        except:
                            logging.error("{0} stage ended with error. Please check.".format(stage))

        # logging.info("Total number of DiffusionPreprocessing processes: {0}".format(str(dwinumruns)))

        # TODO: remove redundant code after testing
        bolds = [f.filename for f in layout.get(subject=subject_label,
                                                type='bold',
                                                extensions=["nii.gz", "nii"])]

        boldsnumruns = len(bolds)

        # logging.info("Total number of fMRI processes: {0}".format(str(boldsnumruns)))
        # totalnumofprocs = 3 + dwinumruns + boldsnumruns


        for fmritcs in bolds:
            fmriname = "_".join(fmritcs.split("sub-")[-1].split("_")[1:-1]).split(".")[0]
            assert fmriname

            fmriscout = fmritcs.replace("_bold", "_sbref")
            if not os.path.exists(fmriscout):
                fmriscout = "NONE"

            fieldmap_set = layout.get_fieldmap(fmritcs)
            # print(fieldmap_set)
            if fieldmap_set and fieldmap_set["type"] == "epi":
                SEPhaseNeg = None
                SEPhasePos = None
                for fieldmap in fieldmap_set["epi"]:
                    enc_dir = layout.get_metadata(fieldmap)["PhaseEncodingDirection"]
                    if "-" in enc_dir:
                        SEPhaseNeg = fieldmap
                    else:
                        SEPhasePos = fieldmap
                echospacing = layout.get_metadata(fmritcs)["EffectiveEchoSpacing"]
                unwarpdir = layout.get_metadata(fmritcs)["PhaseEncodingDirection"]
                unwarpdir = unwarpdir.replace("i","x").replace("j", "y").replace("k", "z")
                if len(unwarpdir) == 2:
                    unwarpdir = "-" + unwarpdir[0]
                dcmethod = "TOPUP"
                biascorrection = "SEBASED"
            else:
                SEPhaseNeg = "NONE"
                SEPhasePos = "NONE"
                echospacing = "NONE"
                unwarpdir = "NONE"
                dcmethod = "NONE"
                biascorrection = "NONE"

            zooms = nibabel.load(fmritcs).get_header().get_zooms()
            fmrires = float(min(zooms[:3]))
            fmrires = "2" # https://github.com/Washington-University/Pipelines/blob/637b35f73697b77dcb1d529902fc55f431f03af7/fMRISurface/scripts/SubcorticalProcessing.sh#L43
            # While running '/usr/bin/wb_command -cifti-create-dense-timeseries /scratch/users/chrisgor/hcp_output2/sub-100307/MNINonLinear/Results/EMOTION/EMOTION_temp_subject.dtseries.nii -volume /scratch/users/chrisgor/hcp_output2/sub-100307/MNINonLinear/Results/EMOTION/EMOTION.nii.gz /scratch/users/chrisgor/hcp_output2/sub-100307/MNINonLinear/ROIs/ROIs.2.nii.gz':
            # ERROR: label volume has a different volume space than data volume


            func_stages_dict= OrderedDict([("fMRIVolume", partial(run_generic_fMRI_volume_processsing,
                                                      path=args.output_dir,
                                                      subject="sub-%s"%subject_label,
                                                      fmriname=fmriname,
                                                      fmritcs=fmritcs,
                                                      fmriscout=fmriscout,
                                                      SEPhaseNeg=SEPhaseNeg,
                                                      SEPhasePos=SEPhasePos,
                                                      echospacing=echospacing,
                                                      unwarpdir=unwarpdir,
                                                      fmrires=fmrires,
                                                      dcmethod=dcmethod,
                                                      biascorrection=biascorrection,
                                                      n_cpus=args.n_cpus)),
                                ("fMRISurface", partial(run_generic_fMRI_surface_processsing,
                                                       path=args.output_dir,
                                                       subject="sub-%s"%subject_label,
                                                       fmriname=fmriname,
                                                       fmrires=fmrires,
                                                       n_cpus=args.n_cpus,
                                                       grayordinatesres=grayordinatesres,
                                                       lowresmesh=lowresmesh))
                                ])
            for stage, stage_func in func_stages_dict.iteritems():
                if stage in args.stages:
                    try:
                        stage_func()
                    except:
                        logging.error("{0} stage ended with error. Please check. Continuing...".format(stage))
        # Process(target= func_stages, args=(func_stages_dict, )).start()

        # task fmris
        tasks =[]
        TemplateDir='/Examples/fsf_templates/'
        FinalSmoothingFWHM="2"
        LowResMesh="32"
        GrayOrdinatesResolution="2"
        OriginalSmoothingFWHM="2"
        Confound="None"
        TemporalFilter="200"
        VolumeBasedProcessing="NO"
        RegNames="NONE"
        ParcellationList="NONE"
        ParcellationFileList="NONE"

        # numruns = set(layout.get(target='run', return_type='id',
        #                          subject=subject_label, type='bold',
        #                          extensions=["nii.gz", "nii"]))
        #
        # for session in numruns:
        #     tfmris = [f.filename for f in layout.get(subject=subject_label,
        #                                             type='bold', run=session,
        #                                             extensions=["nii.gz", "nii"])]
        #     for f in tfmris:
        #         taskname = f.split("task-")[1].split("_")[0]
        #         tasks.append(taskname)
        #
        #     ##check that there are two runs of each task (i.e. AP, PA)
        #     for i in range(0,len(tasks)):
        #         tasknumcheck = layout.get(task=tasks[i], subject=subject_label, type="bold",
        #                                   run=session, extensions=["nii.gz", "nii"])
        #         print(tasknumcheck)
        #         # assert len(tasknumcheck) >= 2
        #
        #     for task in set(tasks):
        #         if task is not "rest":
        #             acqs = set(layout.get(target='acquisition', return_type='id', task=task,
        #                                   subject=subject_label, type='bold', run=session,
        #                                   extensions=["nii.gz", "nii"]))
        #             # assert len(acqs) >= 2
        #             acq = list(acqs)
        #
        #             taskone = 'task-{0}_acq-{1}_{2}'.format(task, acq[0],session)
        #             tasktwo = 'task-{0}_acq-{1}_{2}'.format(task, acq[1], session)
        #             OutDir = '{0}/sub-{1}/MNINonLinear/Results/'.format(args.output_dir, subject_label)
        #
        #             LevelOneTasks= 'task-{0}_acq-{1}_{2}@task-{3}_acq-{4}_{5}'.format(task,acq[0],session,task,acq[1],session)
        #             LevelOneFSFs= 'task-{0}_acq-{1}_{2}@task-{3}_acq-{4}_{5}'.format(task,acq[0],session,task,acq[1],session)
        #             LevelTwoTask= 'task-{0}_{1}'.format(task, session)
        #             LevelTwoFSF= 'task-{0}_{1}'.format(task, session)
        #
        #             func_stages_dict = OrderedDict([('generateLevel1fsf', partial(generate_level1_fsf,
        #                                                                           path=args.output_dir,
        #                                                                           studyfolder=args.output_dir,
        #                                                                           subject="sub-%s" % subject_label,
        #                                                                           taskname=taskone,
        #                                                                           templatedir=TemplateDir,
        #                                                                           outdir='{0}/{1}'.format(OutDir,taskone),
        #                                                                           dir=acq[0],
        #                                                                           n_cpus=args.n_cpus
        #                                                                           )),
        #                                             ('generateLevel1fsfb', partial(generate_level1_fsf,
        #                                                                           path=args.output_dir,
        #                                                                           studyfolder=args.output_dir,
        #                                                                           subject="sub-%s" % subject_label,
        #                                                                           taskname=tasktwo,
        #                                                                           templatedir=TemplateDir,
        #                                                                           outdir='{0}/{1}'.format(OutDir,tasktwo),
        #                                                                           dir=acq[1],
        #                                                                           n_cpus=args.n_cpus
        #                                                                           )),
        #                                             ('generateLevel2fsf', partial(generate_level2_fsf,
        #                                                                           path=args.output_dir,
        #                                                                           studyfolder=args.output_dir,
        #                                                                           subject="sub-%s" % subject_label,
        #                                                                           taskname='task-{0}_{1}'.format(task,session),
        #                                                                           templatedir=TemplateDir,
        #                                                                           outdir='{0}/task-{1}_{2}'.format(OutDir,task,session),
        #                                                                           folderone=taskone,
        #                                                                           foldertwo=tasktwo,
        #                                                                           n_cpus=args.n_cpus
        #                                                                           )),
        #                                             ('TaskfMRIAnalysis', partial(run_task_fmri_analysis,
        #                                                                           path=args.output_dir,
        #                                                                           subject="sub-%s" % subject_label,
        #                                                                           lvl1tasks=LevelOneTasks,
        #                                                                           lvl1fsfs=LevelOneFSFs,
        #                                                                           lvl2task=LevelTwoTask,
        #                                                                           lvl2fsf=LevelTwoFSF,
        #                                                                           lowresmesh=LowResMesh,
        #                                                                           grayordinatesres=grayordinatesres,
        #                                                                           origsmoothingFWHM=OriginalSmoothingFWHM,
        #                                                                           confound=Confound,
        #                                                                           finalsmoothingFWHM=FinalSmoothingFWHM,
        #                                                                           temporalfilter=TemporalFilter,
        #                                                                           vba=VolumeBasedProcessing,
        #                                                                           regname=RegNames,
        #                                                                           parcellation=ParcellationList,
        #                                                                           parcellationfile=ParcellationFileList,
        #                                                                           n_cpus=args.n_cpus
        #                                                                           ))])
        #             for stage, stage_func in func_stages_dict.iteritems():
        #                 if stage in args.stages:
        #                     stage_func()

    elapsed = time.time() - starttime
    elapsed = elapsed / 60
    logging.info("HCP Processing completed. Time duration: {0} minutes".format(str(elapsed)))



