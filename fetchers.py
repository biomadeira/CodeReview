#!/usr/bin/env python2.7

import os
import json
from collections import OrderedDict

from utils import flash
from utils import write_log
from utils import request_info_url

from library import aa_symbols_ext
from library import aa_symbols_rev_ext
from library import aa_physicochemical_full


def fetch_variants_from_ensembl_rest(identifier, sequence, species, ensemblg, ensemblt, ensemblp,
                                     method="ENSEMBL", full=False, form="plain", verbose=True):
    """
    Fetchs variants using the Ensembl REST API. As of July 2014, Ensembl
    variants include 1000 Genomes Project, dbSNP, HAPMAP and COSMIC.
    Gets variation specific info on populations, genotypes and phenotypes.

    More info in
    http://www.ensembl.org/info/genome/variation/sources_documentation.html

    @param identifier: UNIPROT identifier
    @param sequence: UNIPROT sequence
    @param species: UNIPROT species
    @param ensemblg: list of ENSEMBL gene identifiers
    @param ensemblt: list of ENSEMBL transcript identifiers
    @param ensemblp: list of ENSEMBL protein identifiers
    @param method: ENSEMBL query of transcript variation "ENSEMBL" or somatic var "COSMIC"
    @param full: Boolean (full variant information with
        populations, genotypes, phenotypes, etc.)
    @param form: output format
    @param verbose: Boolean
    @return: returns a list of variants
    """

    if verbose:
        flash("Fetching Variants from Ensembl using REST API...")

    # example lines to parse
    """
    [{
        "end": 0,
        "feature_type": "somatic_transcript_variation",
        "start": 0,
        "type": "intron_variant",
        "codons": "",
        "polyphen": null,
        "seq_region_name": "ENSP00000470284",
        "sift": null,
        "residues": "",
        "allele": "G/A",
        "translation": "ENSP00000470284",
        "minor_allele_frequency": null,
        "id": "COSM1393501"
    },
    {
        "end": 349,
        "feature_type": "somatic_transcript_variation",
        "start": 349,
        "type": "missense_variant",
        "codons": "gGg/gTg",
        "polyphen": 1,
        "seq_region_name": "ENSP00000470284",
        "sift": 0,
        "residues": "G/V",
        "allele": "G/T",
        "translation": "ENSP00000470284",
        "minor_allele_frequency": null,
        "id": "COSM48634"
    }]"""

    # example lines to parse
    """
    {"mappings": [
        {
            "end": 39664264,
            "start": 39664264,
            "coord_system": "chromosome",
            "allele_string": "G/A",
            "seq_region_name": "19",
            "location": "19:39664264-39664264",
            "strand": 1
        }
    ],
    "populations": [
        {
            "allele": "G",
            "frequency": "0.99988",
            "allele_count": "8365",
            "population": "ESP6500:European_American"
        },
        {
            "allele": "A",
            "frequency": "0.000119531",
            "allele_count": "1",
            "population": "ESP6500:European_American"
        },
    ],
    "most_severe_consequence": "Missense variant",
    "var_class": "SNP",
    "evidence": [
        "ESP"
    ],
    "source": "Variants (including SNPs and indels) imported from dbSNP",
    "ambiguity": "R",
    "MAF": null,
    "ancestral_allele": "G",
    "name": "rs371973365",
    "synonyms": [],
    "genotypes": [],
    "phenotypes": [],
    }
    """

    variants = []

    ensg_list = ensemblg
    enst_list = ensemblt
    ensp_list = ensemblp

    for ensemblt, ensemblp in zip(enst_list, ensp_list):

        if verbose:
            flash("Ensembl Protein %s..." % ensemblp)
        # first compares the uniprot sequence to the ensemblp sequence
        url = "http://rest.ensembl.org/sequence/id/%s" % ensemblp
        url += "?content-type=text/plain;type=protein"
        read = request_info_url(identifier, url, lines=True, verbose=verbose)

        if read != []:
            ensembl_seq = read[0].rstrip("\r\n")
            if len(sequence) == len(ensembl_seq):

                # get ensembl gene ID: important if there are >1 gene ids
                if len(ensg_list) == 1:
                    ensemblg = ensg_list[0]
                else:
                    url = "http://rest.ensembl.org/overlap/id/%s" % ensemblt
                    url += "?feature=transcript;content-type=application/json"
                    read = request_info_url(identifier, url, lines=False, verbose=verbose)

                    ensemblg = "-"
                    if read != "":
                        info = json.loads(read)
                        for entry in info:
                            if entry["id"] == ensemblt:
                                ensemblg = entry["Parent"]

                # gets the variants: either ENSEMBL (transcript variants)
                # or COSMIC (somatic variants)
                url = "http://rest.ensembl.org/overlap/translation/%s" % ensemblp
                if method == "ENSEMBL":
                    url += "?feature=transcript_variation;"
                elif method == "COSMIC":
                    url += "?feature=somatic_transcript_variation;"
                else:
                    url += "?feature=transcript_variation;feature=somatic_transcript_variation;"
                url += "content-type=application/json"
                read = request_info_url(identifier, url, lines=False, verbose=verbose)

                if read != "":
                    # gets result and parse arguments
                    info = json.loads(read)
                    for entry in info:
                        vres = entry["residues"]
                        vtype = entry["type"]
                        vid = entry["id"]
                        vsite = int(entry["start"])

                        # filtering out synonymous variants and other more complicated frameshift variants
                        if vres != "" and "/" in vres and vtype != "synonymous_variant" and len(vres) == 3 and \
                            (vid[0:2] == "rs" or vid[0:4] == "COSM") and vsite < len(sequence):

                            # loading variation dictionary
                            variant = OrderedDict()

                            vres1 = aa_symbols_rev_ext[vres.split("/")[0]]
                            vres2 = aa_symbols_rev_ext[vres.split("/")[1]]

                            variant["VARIATION"] = "p.%s%s%s" % (vres1, vsite, vres2)
                            variant["SITE"] = str(vsite)

                            variant["RES1"] = vres1
                            variant["RES2"] = vres2
                            try:
                                variant["RES1_PROP"] = aa_physicochemical_full[vres1]
                            except:
                                 variant["RES1_PROP"] = []
                            try:
                                variant["RES2_PROP"] = aa_physicochemical_full[vres2]
                            except:
                                 variant["RES2_PROP"] = []
                            variant["SOURCE"] = vid
                            variant["ENSEMBL_GENE"] = ensemblg
                            variant["ENSEMBL_TRANSCRIPT"] = ensemblt
                            variant["ENSEMBL_PROTEIN"] = ensemblp
                            variant["TYPE"] = " ".join(entry["type"].split("_"))
                            variant["FEATURE_TYPE"] = entry["feature_type"]
                            variant["CODONS"] = entry["codons"]
                            variant["ALLELE"] = entry["allele"]
                            alle_freq = entry["minor_allele_frequency"]
                            if alle_freq is None:
                                variant["ALLELE_FREQUENCY"] = "-"
                            else:
                                variant["ALLELE_FREQUENCY"] = alle_freq

                            variant["LOCATION"] = "-"
                            variant["CHROMOSSOME"] = "-"
                            variant["TRAIT"] = "-"
                            variant["TRAIT_DB"] = "-"

                            if full:
                                # get id specific info on populations, genotypes and phenotypes
                                url = "http://rest.ensembl.org/variation/%s/%s" % (species, vid)
                                url += "?pops=1;phenotypes=1;genotypes=1;"
                                url += "content-type=application/json"
                                read = request_info_url(identifier, url, lines=False, verbose=verbose)

                                if read != "":
                                    # get result and parse arguments
                                    entry2 = json.loads(read)
                                    variant["LOCATION"] = entry2["mappings"][0]["location"]
                                    variant["CHROMOSSOME"] = entry2["mappings"][0]["assembly_name"]

                                    try:
                                        variant["TRAIT"] = entry2["phenotypes"][0]["trait"]
                                    except:
                                        variant["TRAIT"] = "-"
                                    try:
                                        variant["TRAIT_DB"] = entry2["phenotypes"][0]["source"]
                                    except:
                                        variant["TRAIT_DB"] = "-"

                                    variant["MAPPINGS"] = entry2["mappings"]
                                    variant["GENOTYPES"] = entry2["genotypes"]
                                    variant["PHENOTYPES"] = entry2["phenotypes"]
                                    variant["SYNONYMS"] = entry2["synonyms"]
                                    variant["POPULATIONS"] = entry2["populations"]
                                    variant["EVIDENCE"] = entry2["evidence"]
                                    variant["CONSEQUENCE"] = entry2["most_severe_consequence"]
                                    variant["SOURCE_DB"] = entry2["source"]

                            else:
                                # get only essential information for each id
                                url = "http://rest.ensembl.org/variation/%s/%s" % (species, vid)
                                url += "?phenotypes=1;"
                                url += "content-type=application/json"
                                read = request_info_url(identifier, url, lines=False, verbose=verbose)

                                if read != "":
                                    # get result and parse arguments
                                    entry2 = json.loads(read)

                                    variant["LOCATION"] = entry2["mappings"][0]["location"]
                                    variant["CHROMOSSOME"] = entry2["mappings"][0]["assembly_name"]
                                    try:
                                        variant["TRAIT"] = entry2["phenotypes"][0]["trait"]
                                    except:
                                        variant["TRAIT"] = "-"
                                    try:
                                        variant["TRAIT_DB"] = entry2["phenotypes"][0]["source"]
                                    except:
                                        variant["TRAIT_DB"] = "-"

                            # test for mapping between ensembl and uniprot sequences
                            try:
                                # compares the residue that is supposed to match the variation entry with
                                # the one observed in uniprot sequence

                                if vres1 != "---" and vres1 != "***":
                                    assert aa_symbols_ext[vres1] == sequence[vsite - 1]

                                if variant not in variants:
                                    variants.append(variant)

                            except:
                                message = "%s\tWarning: %s in sequence position %s, does not match the %s for %s" % \
                                            (identifier, aa_symbols_ext[vres1], vsite, sequence[vsite - 1], ensemblp)
                                print(message)
                                path = os.getcwd() + "/"
                                output_file = "error_variants.log"
                                write_log(message, path + output_file)


    if form is "plain":
        print(variants)
    elif form is "json":
        print(json.dumps(variants, sort_keys=False, indent=4))

    return variants


if __name__ == "__main__":
    # testing routines
    pass