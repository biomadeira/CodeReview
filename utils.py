#!/usr/bin/env python2.7

import os
import sys
import requests
import time
from datetime import datetime
from collections import OrderedDict


def current_time():
    """
    Gets current time to load into log files.

    @return: returns formatted time
    """

    date = datetime.now()
    year = date.strftime("%Y")
    month = date.strftime("%m")
    day = date.strftime("%d")
    hour = date.strftime("%H")
    minute = date.strftime("%M")
    second = date.strftime("%S")
    curr_time = "%s/%s/%s %s:%s:%s" % (day, month, year, hour, minute, second)

    return curr_time


def flash(message):
    """
    Flashes a message out.

    @param message: input message str()
    """

    message = str(message)

    print(message)
    sys.stdout.flush()
    return


def create_directory(directory):
    """
    Creates a directory structure if it does not exist.

    @param directory: output directory name
    """

    path = os.getcwd() + "/"
    directory = "%s%s" % (path, directory)
    if not os.path.exists(directory):
        os.makedirs(directory)

    return


def load_lines(input_path, verbose=True):
    """
    Parses and loads lines from an input file.

    @param input_path: path to the fasta file
    @param verbose: Boolean
    @return: returns fasta object
    """

    if verbose:
        flash("Loading lines from %s..." % input_path)

    with open(input_path) as readlines:
        lines = [line for line in readlines]

    return lines


def write_log(message, output_file):
    """
    Appends a message to a log file.

    @param message: message to be printed
    @param output_file: path to the output file
    """

    with open(output_file, "a") as outlog:
        outlog.write(message + "\n")

    return


def request_info_url(identifier, url, lines=False, verbose=True):
    """
    Gets formatted content from the provided url.

    If it fails, tries to access it a few more times to make sure it
    was not a temporary error.

    @param identifier: identifier used for error handling
    @param url: input URL
    @param lines: Boolean outputs either a list or a string
    @param verbose: Boolean
    @return: returns a data object from *requests*
    """

    info = ""
    req = requests.get(url)
    if verbose:
        print req.status_code, url

    if req.status_code == 200:
        if lines:
            info = [line for line in req.iter_lines()]
        else:
            info = req.text
    elif req.status_code == 400 or req.status_code == 404 or req.status_code == 502:
        if lines:
            info = []
        else:
            info = ""
    elif req.status_code == 429 or req.status_code == 503 or req.status_code == 504:
        time.sleep(0.5)
        request_info_url(identifier, url, lines=lines, verbose=verbose)
    else:
        status = req.status_code
        message = "%s\tError %s: Could not download the data from %s at this time..." % (identifier, status, url)
        print(message)
        path = os.getcwd() + "/"
        output_file = "e_url.log"
        with open(path + output_file, "a") as outlog:
            outlog.write(message + "\n")

    return info


def uniprot_summary_mapping(entities_dict, residues_dict, form="plain"):
    """
    Gets mappings of uniprot/domain/ptms/variants/features.

    @param entities_dict: Dictionary of dictionaries containing information
    per residue.
    @param residues_dict: Dictionary of dictionaries containing information
    per residue.
    @param form: output format
    @return: returns a dictionary of 'alignments'
    """

    assert len(entities_dict) == len(residues_dict)

    # get feature entries that appear in the residues entries
    # get a 'alignment' view of residues, ss and domains
    # all under a dictionary structure
    maps = OrderedDict()

    # RES all
    fet_align = ""
    for res in residues_dict:
        try:
            fet_align += residues_dict[res]["UNIPROT_NAME"]
        except:
            fet_align += "-"
    maps["UNIPROT_NAME"] = fet_align

    # Variants
    fet_align = ""
    for res in entities_dict:
        if entities_dict[res]["VARIANTS"] != []:
            fet_align += "V"
        else:
            fet_align += "-"
    maps["VARIANTS"] = fet_align

    # Mutations
    fet_align = ""
    for res in entities_dict:
        if entities_dict[res]["MUTATIONS"] != []:
            fet_align += "M"
        else:
            fet_align += "-"
    maps["MUTATIONS"] = fet_align

    if form is "plain":
        print(maps)

    return maps


if __name__ == "__main__":
    # testing routines
    pass