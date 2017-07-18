#!/usr/bin/python
"""
Parses the output HCPProcessing.log file and generates json files indicating the status of the processes.
"""

from __future__ import unicode_literals, print_function
import os, sys, json, time
import argparse

global DERIV_BASE_DIR
global PUBLIC
global SUBJ

DERIV_BASE_DIR = ""
SUBJ = ""
LOGFN = "{0}_HCPPipelines.log".format(SUBJ)
STRUCT = 3
DWI = 0
FMRI = 0
ALL_DONE = False

FIN = False

def createStatusPath(subject):
    return os.path.join(DERIV_BASE_DIR, 'output', 'logs', 'sub-{0}_HCPPipelines.log'.format(subject))

def createHCPStatePath():
    return os.path.join(PUBLIC, "hcp_pipelines_state.json")



def generateJSON(processingComplete):
    j = {}

    if processingComplete:
        j['status'] = "All processing complete"
    else:
        j['status'] = "Processing in progress"

    j['subject'] = SUBJ

    currentJson = {}
    currentJson['name'] = SUBJ
    seenNotDone = False

    with open(createStatusPath(SUBJ), 'rb') as f:
        f.seek(-2, 2)
        last = f.readlines()[-1].decode()

    if "HCP Processing completed" in last:
        seenNotDone = True

    currentJson['state'] = 0

    if not seenNotDone:
        global ALL_DONE
        ALL_DONE = True
    return json.dumps(j)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Reads and parses log files to create status text files.\n")
    # parser.add_argument('logfile', help='log file to read status from\n')
    parser.add_argument('subject', help='subject ID')
    parser.add_argument('publichtml', help='public_html is the path to save hcppipeline_state.json file\n')

    args = parser.parse_args()

    DERIV_BASE_DIR = os.path.abspath(os.path.dirname(LOGFN))
    PUBLIC = os.path.abspath(args.publichtml)
    SUBJ = args.subject

    while not ALL_DONE:
        jsonToWrite = generateJSON(False)
        hcpState = open(createHCPStatePath(), 'w')
        hcpState.write(jsonToWrite)
        hcpState.close()
        time.sleep(1)

    jsonToWrite = generateJSON(True)
    hcpState = open(createHCPStatePath(), 'w')
    hcpState.write(jsonToWrite)
    hcpState.close()

    exit(0)