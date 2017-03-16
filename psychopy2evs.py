# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

Edited from "gen_banda_onsets.py" created by Mathias Goncalves @ MGH:
    Generates bold onsets for facematching, gambling, and
    conflict tasks. Logfile directory (-l) and output
    directory (-o) have already been set for specified
    locations within this project, but can be changed
    with their appropriate flags if needed.
    
Edited by Amber Leaver, 02/27/2017 to remove code for irrelevant tasks and to add new code for CARIT task
Edited by Yeun Kim 03/15/17.

"""

import os
from glob import glob
from csv import reader, writer
from collections import OrderedDict
import argparse
import pandas as pd

def write_onset(conds, outdir, duration):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for idx,cond in enumerate(conds):
        with open(os.path.join(outdir,'cond%03d.txt'%(idx+1)),'w') as csvfile:
            csvwriter=writer(csvfile, delimiter='\t')
            for onset in conds[cond]:
                csvwriter.writerow(['%.1f'%onset, '%.1f'%duration, '%.1f'%1])
    return


def get_face(logdir, run, outdir=None):
    log = [x for x in os.listdir(logdir) if x.endswith('.csv') and 'Scanner' in x and '_AB_FaceMatching' in x][-1]
    df = pd.read_csv(os.path.join(logdir, log))
    block_dur = 18.0
    objs = []
    happ = []
    neut = []
    fear = []
    fix = [] 
    current_cond = ''
    time = 0
    # unsure about fixation accuracy
    for i,row in df.iterrows():
        try:
            check_cond = row['Condition']
            check_run = row['Run']
            if 'fixation' in check_cond and run in check_run:
                fix.append(time)
                time += block_dur
            elif 'Happy' in check_cond and current_cond is not check_cond and run in check_run:
                happ.append(time)
                time += block_dur
                current_cond = check_cond
            elif 'Neutral' in check_cond and current_cond is not check_cond and run in check_run:
                neut.append(time)
                time += block_dur
                current_cond = check_cond
            elif 'Fearful' in check_cond and current_cond is not check_cond and run in check_run:
                fear.append(time)
                time += block_dur
                current_cond = check_cond
            elif 'Object' in check_cond and current_cond is not check_cond and run in check_run:
                objs.append(time)
                time += block_dur
                current_cond = check_cond
            else:
                continue
        except TypeError:
            #skipping blank lines
            continue        
    conds = OrderedDict([('fear', fear), ('happy', happ), ('neutral', neut), ('object',objs)])
    if outdir:
        write_onset(conds, outdir, block_dur)
    return conds

def get_carit(logdir, run, outdir=None):
    log = [x for x in os.listdir(logdir) if x.endswith('.csv') and 'wide' in x and 'CARIT' in x][-1]
    df = pd.read_csv(os.path.join(logdir, log))
    hit = []
    miss = []
    cr = []
    fa = []
    countdown = [] #let's not model this for now
    time = 0
    
    for i,row in df.iterrows():
        try:
            check_cond = row['corrRespTrialType']
            check_run = row['nRuns']
            time = row['shapeStartTime']
            if 'Hit' in check_cond:
                hit.append(time)
            elif 'Miss' in check_cond:
                miss.append(time)
            elif 'corReject' in check_cond:
                cr.append(time)
            elif 'falseAlarm' in check_cond:
                fa.append(time)
            else:
                continue
        except TypeError:
            #skipping blank lines
            continue
    conds = OrderedDict([('Hit', hit), ('Miss', miss), ('corReject', cr), ('falseAlarm', fa)])
    if outdir:
        write_onset(conds, outdir, 0.6)
    return conds


def gen_onsets(subjs, logdir, outdir):
    # need error catching
    for subj in subjs:
        print('Processing %s...'%subj)
        subjlogs = os.path.join(logdir, '{}'.format(subj))
        subjout = os.path.join(outdir, 'EVs')
        for i,run in enumerate(['A', 'B'],1):
            try:
                #get_face(subjlogs, run, os.path.join(subjout, 'task001_run%03d'%i))
                get_face(subjlogs, run, os.path.join(subjout, 'task-face_run-%02d'%i))
            except OSError:
                print('Error: Unable to find face-matching run %d'%i)
        for i,run in enumerate(['1'],1):
            try:
                #get_gamble(subjlogs, run, os.path.join(subjout, 'task002_run%03d'%i))
                get_carit(subjlogs, run, os.path.join(subjout, 'task-carit_run-%02d'%i))
            except OSError:
                    print('Error: Unable to find CARIT run %d'%i)
    
    
        
if __name__ == '__main__':
    docstr = '\n'.join((__doc__,
"""
           Example usage:

           python psychopy2evs.py -s s000701 s000702
"""))
    defstr = ' (default %(default)s)'
    parser = argparse.ArgumentParser(prog='psychopy2evs.py',
                                     description=docstr)
    parser.add_argument('-s', dest='subjects', default=[],
                        type=str, nargs='+',
                        help="Subject/session IDs (e.g. s000701)", required=True)
    parser.add_argument('-l', dest='logdir', type=str,
                        help="Onset files directory")
    parser.add_argument('-o', dest='outdir',
                        help="Output directory")
    parser.add_argument('-m', dest='model', default=1, type=int,
                        help="Model index" + defstr)
    args = parser.parse_args()
    
    gen_onsets(args.subjects, args.logdir, args.outdir)
