#!/usr/bin/env python3

# Copyright 2016-2026, Robin Kobus (github.com/funatiq)
#
# This file is part of the MetaCache taxonomic sequence classification tool.
#
# MetaCache is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MetaCache is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MetaCache.  If not, see <http://www.gnu.org/licenses/>.


# - query datatabase with options '-abundances' and '-abundance-per species':
#       ./metacache query dbname read_files -out classification.txt -abundances abund.txt -abundance-per species
# - then run this script on 'abund.txt':
#       ./krona-from-abundances abund.txt 


import os, argparse, subprocess

# parse arguments
parser = argparse.ArgumentParser(description="Generate Krona diagrams from MetaCache abundance files. Requires Krona to be installed. Query MetaCache database with options '-abundances <abundance_file> -abundance-per species'.")
parser.add_argument('abundance_file',type=argparse.FileType('r'),
	help='metacache abundance file')

args=parser.parse_args()

abundance_fname = os.path.splitext(args.abundance_file.name)[0]

abundance_orig_fname = abundance_fname + "_orig.txt"
abundance_est_fname = abundance_fname + "_est.txt"

abundance_orig_file = open(abundance_orig_fname, 'w')
abundance_est_file = open(abundance_est_fname, 'w')

stage = 0

for line in args.abundance_file:
	if stage == 0:
		abundance_orig_file.write(line)
		if line[0] != '#':
			stage = 1
	elif stage == 1:
		if line[0] != '#':
			abundance_orig_file.write(line)
		else:
			stage = 2
	if stage == 2:
		abundance_est_file.write(line)

abundance_orig_file.close()
abundance_est_file.close()

krona_fname = abundance_fname + ".krona.html"

subprocess.run(["ktImportTaxonomy", "-t", "3", "-s", "0", "-m", "5", "-o", krona_fname, abundance_est_fname])
