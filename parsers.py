#!/usr/bin/env python2.7

from utils import flash


def parse_information_from_uniprot(fasta_object, form="plain", verbose=True):
    """
    Gets fasta SEQ from uniprot with header as optional

    @param fasta_object: UniProt fasta lines
    @param form: output format
    @param verbose: Boolean
    @return: returns sequence, name, gene and species
    """

    if verbose:
        flash("Parsing Features from UniProt...")

    sequence = ""
    name = ""
    gene = ""
    species = ""
    fasta = ""

    for count, line in enumerate(fasta_object):
        if count is 0:
            header = line.rstrip("\r\n")
            # example header: it includes the '>' (fasta format)
            """
            >sp|P04217|A1BG_HUMAN Alpha-1B-glycoprotein OS=Homo sapiens GN=A1BG PE=1 SV=4
            """
            try:
                name = (header.split("|")[2]).split()[0]
            except:
                name = "-"
            try:
                gene = (header.split("GN=")[1]).split()[0]
            except:
                gene = "-"
            try:
                species = ("_".join((header.split("OS=")[1]).split()[0:2])).lower()
            except:
                species = "-"
        else:
            sequence += line.rstrip("\r\n")
        fasta += line

    if form is "plain":
        print(sequence, name, gene, species)

    return sequence, name, gene, species


def parse_ensembl_from_uniprot(identifier, txt_object, form="plain", verbose=True):
    """
    Gets Ensembl ids for Gene, Transcript and Protein.
    Converts UniProt id into respective Ensembl ids.

    @param identifier: UniProt ID
    @param txt_object: UniProt txt lines object
    @param form: output format
    @param verbose: Boolean
    @return: returns ensembl gene, transcript and protein identifiers
    """

    if verbose:
        flash("Fetching Ensembl IDs from UniProt...")

    # example lines
    """
    DR   Ensembl; ENSOCUT00000003872; ENSOCUP00000003357; ENSOCUG00000003871.

    DR   Ensembl; ENST00000269305; ENSP00000269305; ENSG00000141510. [P04637-1]
    DR   Ensembl; ENST00000420246; ENSP00000391127; ENSG00000141510. [P04637-2]
    DR   Ensembl; ENST00000445888; ENSP00000391478; ENSG00000141510. [P04637-1]
    DR   Ensembl; ENST00000455263; ENSP00000398846; ENSG00000141510. [P04637-3]

    DR   Ensembl; ENST00000553106; ENSP00000448059; ENSG00000171759.
    """

    ensemblg = []
    ensemblt = []
    ensemblp = []

    for line in txt_object:
        if "DR   Ensembl;" in line[0:14]:
            skip = False
            line = line.rstrip("\r\n")
            line = line.split(";")
            enst = line[1].strip()
            ensp = line[2].strip()
            ensg = line[3].strip()

            if "[" in ensg:
                ensg = ensg.split()
                uni = ensg[1][1:-1]
                if "%s" % identifier == uni:
                    ensg = ensg[0].rstrip(".")
                elif "%s-1" % identifier == uni:
                    ensg = ensg[0].rstrip(".")
                else:
                    skip = True
            else:
                ensg = ensg.rstrip(".")

            if not skip:
                if ensg not in ensemblg:
                    ensemblg.append(ensg)
                if enst not in ensemblt:
                    ensemblt.append(enst)
                if ensp not in ensemblp:
                    ensemblp.append(ensp)

    if form is "plain":
        print(ensemblg, ensemblt, ensemblp)

    return ensemblg, ensemblt, ensemblp


if __name__ == "__main__":
    # testing routines
    pass