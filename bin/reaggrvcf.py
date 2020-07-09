#!/usr/bin/env python3
import sys
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='extract annotation for specific position')
    parser.add_argument('--vcfi',type=str,required=True, help="file contains chro and list of files with annotation")
    parser.add_argument('--vcfout', type=str,help="output")
    args = parser.parse_args()
    return args

args = parseArguments()


