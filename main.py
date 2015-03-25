#!/usr/bin/env python2.7

import os
import json
from collections import OrderedDict

from parsers import parse_information_from_uniprot
from parsers import parse_ensembl_from_uniprot
from fetchers import fetch_variants_from_ensembl_rest

from utils import flash
from utils import current_time
from utils import create_directory
from utils import write_log
from utils import request_info_url
from utils import load_lines
from utils import uniprot_summary_mapping
from library import ensembl_species


class CoreSEQUENCE(object):
    """
    Loads information based on a UNIPROT identifier.
    This includes path to the local files, web fetched
    data and SIFTS features.
    """

    def __init__(self, identifier=None, db=True, form="plain", verbose=True):
        """
        Init method. Quickly loads uniprot data from fasta and txt.

        @param identifier: UNIPROT identifier
        @param form: output format
        @param verbose: Boolean
        @return: returns sequence, name, gene, species, fasta
                ensemblg, ensemblt, ensemblp, txt, domains,
               variants, ptms, features, etc.
        """

        self.identifier = identifier
        self.form = form
        self.verbose = verbose
        if self.identifier is not None:
            if self.verbose:
                flash("Loading UNIPROT ID %s..." % self.identifier)

            # assert identifier length
            try:
                assert len(identifier) is 6
            except:
                message = "Error: UniProt identifiers are usually 6-character long..."
                if self.verbose:
                    flash(message)
                raise AssertionError(message)

            # assert current work directory
            cwd_path = os.getcwd()
            self.uniprot_path_uniprot = cwd_path + "/Data"
            create_directory("Data")
            self.uniprot_fullpath_fasta = "%s/%s.%s" % (self.uniprot_path_uniprot, self.identifier, "fasta")
            if os.path.isfile(self.uniprot_fullpath_fasta):
                try:
                    self.uniprot_fasta = True
                    self.uniprot_fasta_object = load_lines(self.uniprot_fullpath_fasta, verbose=self.verbose)
                except:
                    self.uniprot_fasta = False
            else:
                # load source information into the object so it doesn't need be accessed all the time
                url = 'http://www.uniprot.org/uniprot/%s.fasta' % identifier
                self.uniprot_fasta_object = request_info_url(identifier, url, lines=True, verbose=self.verbose)

                if self.uniprot_fasta_object != []:
                    self.uniprot_fasta = True
                    with open(self.uniprot_fullpath_fasta, "w") as output:
                        for line in self.uniprot_fasta_object:
                            output.write(line + "\n")
                else:
                    self.uniprot_fasta = False
                    message = "Warning: %s.fasta not available for download." % identifier
                    print(message)
                    path = os.getcwd() + "/"
                    output_file = "error_uniprot.log"
                    write_log("%s\t%s" % (current_time(), message), path + output_file)

            self.uniprot_fullpath_txt = "%s/%s.%s" % (self.uniprot_path_uniprot, self.identifier, "txt")
            if os.path.isfile(self.uniprot_fullpath_txt):
                try:
                    self.uniprot_txt = True
                    self.uniprot_txt_object = load_lines(self.uniprot_fullpath_txt, verbose=self.verbose)
                except:
                    self.uniprot_txt = False
            else:
                url = 'http://www.uniprot.org/uniprot/%s.txt' % identifier
                self.uniprot_txt_object = request_info_url(identifier, url, lines=True, verbose=self.verbose)

                if self.uniprot_txt_object != []:
                    self.uniprot_txt = True
                    with open(self.uniprot_fullpath_txt, "w") as output:
                        for line in self.uniprot_txt_object:
                            output.write(line + "\n")
                else:
                    self.uniprot_txt = False
                    message = "Warning: %s.txt not available for download." % identifier
                    print(message)
                    path = os.getcwd() + "/"
                    output_file = "error_uniprot.log"
                    write_log("%s\t%s" % (current_time(), message), path + output_file)

            # additional self variables to updated once the methods are launched
            self.sequence = None
            self.name = None
            self.gene = None
            self.species = None
            self.ensemblg = None
            self.ensemblt = None
            self.ensemblp = None

            self.information = None
            self.variants = None
            self.mutations = None
            self.entities = None
            self.residues = None
            self.summary = None

            # trying to load data from DB
            if db:
                # I removed this bit from here to make the code shorter
                pass

        return

    def load_identifier(self, identifier, db, form, verbose):
        """Initiates the class with a UniProt identifier"""

        return self.__init__(identifier, db, form, verbose)

    def get_information(self):
        """
        Gets a list of features for that UniProt identifier
        from UNIPROT and ENSEMBL.
        """

        if self.verbose:
            flash("Getting Information...")

        feat = OrderedDict()
        if self.uniprot_fasta and self.uniprot_txt:
            self.sequence, self.name, self.gene, self.species = parse_information_from_uniprot(self.uniprot_fasta_object,
                                                                                               form=self.form,
                                                                                               verbose=self.verbose)
            self.ensemblg, self.ensemblt, self.ensemblp = parse_ensembl_from_uniprot(self.identifier,
                                                                                     self.uniprot_txt_object,
                                                                                     form=self.form,
                                                                                     verbose=self.verbose)

            feat["NAME"] = self.name
            feat["SEQUENCE"] = self.sequence
            feat["GENE"] = self.gene
            feat["SPECIES"] = self.species
            feat["ENSEMBL_GENE"] = self.ensemblg
            feat["ENSEMBL_TRANSCRIPT"] = self.ensemblt
            feat["ENSEMBL_PROTEIN"] = self.ensemblp

        self.information = feat

        return self.information

    def get_variants(self):
        """
        Gets a list of variants for that UniProt identifier.
        Variants are from UNIPROT and ENSEMBL (as of July 2014).
        """

        if self.verbose:
            flash("Getting Variants Information (Ensembl)...")

        try:
            assert isinstance(self.information, dict)
        except:
            self.get_information()

        self.variants = []
        if self.uniprot_fasta and self.uniprot_txt:
            if self.species in ensembl_species:
                self.variants = fetch_variants_from_ensembl_rest(self.identifier, self.sequence, self.species,
                                                                 self.ensemblg, self.ensemblt, self.ensemblp,
                                                                 method="ENSEMBL", full=False,
                                                                 form=self.form, verbose=self.verbose)
        return self.variants

    def get_mutations(self):
        """
        Gets a list of variants for that Uniprot identifier.
        Variants are from UNIPROT and ENSEMBL (as of July 2014).
        """

        if self.verbose:
            flash("Getting Mutations Information (Ensembl)...")

        try:
            assert isinstance(self.information, dict)
        except:
            self.get_information()

        self.mutations = []
        if self.uniprot_fasta and self.uniprot_txt:
            if self.species in ensembl_species:

                self.mutations = fetch_variants_from_ensembl_rest(self.identifier, self.sequence, self.species,
                                                                  self.ensemblg, self.ensemblt, self.ensemblp,
                                                                  method="COSMIC", full=False,
                                                                  form=self.form, verbose=self.verbose)

        return self.mutations

    def get_entities(self):
        """
        Gets Residues information similar to get_residues in CoreSIFTS.
        """

        if self.verbose:
            flash("Getting Entities...")

        try:
            assert isinstance(self.variants, list)
        except:
            self.get_variants()

        try:
            assert isinstance(self.mutations, list)
        except:
            self.get_mutations()

        # got over each residue in the sequence
        res_entries = OrderedDict()
        if self.uniprot_fasta and self.uniprot_txt:
            for i in range(len(self.sequence)):
                site = str(i + 1)
                res = OrderedDict()

                var_list = []
                for var in self.variants:
                    if var["SITE"] == site:
                        var_list.append(var["SOURCE"])
                res["VARIANTS"] = var_list

                var_list = []
                for var in self.mutations:
                    if var["SITE"] == site:
                        var_list.append(var["SOURCE"])
                res["MUTATIONS"] = var_list

                res_entries[str(i + 1)] = res

        self.entities = res_entries

        return self.entities

    def get_residues(self):
        """
        Gets Residues information similar to get_residues in CoreSIFTS.
        """

        if self.verbose:
            flash("Getting Residues...")

        try:
            assert isinstance(self.variants, list)
        except:
            self.get_variants()

        try:
            assert isinstance(self.mutations, list)
        except:
            self.get_mutations()

        # get a dictionary of var and site
        variants = OrderedDict()
        for var in self.variants:
            source = var["SOURCE"]
            variation = var["VARIATION"]
            site = var["SITE"]
            variants["%s:%s" % (source, variation)] = [site]

        mutations = OrderedDict()
        for var in self.mutations:
            source = var["SOURCE"]
            variation = var["VARIATION"]
            site = var["SITE"]
            mutations["%s:%s" % (source, variation)] = [site]

        # got over each residue in the sequence
        res_entries = OrderedDict()
        if self.uniprot_fasta and self.uniprot_txt:
            for i in range(len(self.sequence)):
                site = str(i + 1)
                res = OrderedDict()
                res["UNIPROT_ID"] = site
                res["UNIPROT_NAME"] = self.sequence[i]
                res["UNIPROT_ACC"] = self.identifier

                for var in variants:
                    sites = variants[var]
                    if site in sites:
                        res[var] = self.sequence[i]

                for var in mutations:
                    sites = mutations[var]
                    if site in sites:
                        res[var] = self.sequence[i]

                res_entries[str(i + 1)] = res

        self.residues = res_entries

        return self.residues

    def get_summary(self):
        """
        Gets a summary view similar to get_summary in CoreSIFTS.
        """

        if self.verbose:
            flash("Getting Summary...")

        try:
            assert isinstance(self.entities, dict)
        except:
            self.get_entities()

        try:
            assert isinstance(self.residues, dict)
        except:
            self.get_residues()

        self.summary = OrderedDict()
        if self.uniprot_fasta and self.uniprot_txt:
            self.summary = uniprot_summary_mapping(self.entities, self.residues, form=self.form)

        return self.summary


def main_handler(identifier, form="plain", verbose=True):
    """
    Calls program pipelines according to the input argument.

    @param identifier: UNIPROT identifier
    @param form: output format
    @param verbose: Boolean
    @return: returns a dictionary with variation information
    """

    information = OrderedDict()

    sequence = CoreSEQUENCE(identifier, db=False, form="", verbose=verbose)
    info = sequence.get_information()
    variants = sequence.get_variants()
    mutations = sequence.get_mutations()
    residues = sequence.get_residues()
    summary = sequence.get_summary()

    information["INFORMATION"] = info
    information["VARIANTS"] = variants
    information["MUTATIONS"] = mutations
    information["RESIDUES"] = residues
    information["SUMMARY"] = summary

    if form is "plain":
        print(information)
    elif form is "json":
        print(json.dumps(information, sort_keys=False, indent=4))

    return information


def main():
    """
    Simple main option parser. It does not use OptionParser and ArgumentParser,
    because I am using a very odd load of arguments and joint strings.
    """

    import argparse
    parser = argparse.ArgumentParser(description='Gets variants for input UniProt ID(s) (example P00439).')

    parser.add_argument('-i', metavar='N', type=str, nargs='+',
                        dest='entries', help='input UniProt ID(s)')
    parser.add_argument('-v', '-verbose', dest='verbose', default=False,
                        help='turns verbosity on', action='store_true')

    args = parser.parse_args()
    if isinstance(args.entries, list):
        for identifier in args.entries:
            main_handler(identifier, form="json", verbose=args.verbose)
    else:
        print "...No input provided! Check the program help with 'python main.py -h'"


if __name__ == "__main__":
    main()