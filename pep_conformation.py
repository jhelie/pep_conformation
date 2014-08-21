#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb="0.0.1"
parser = argparse.ArgumentParser(prog='pep_conformation', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
***********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/pep_conformation
DOI: 
***********************************************

[ DESCRIPTION ]

This script identities the conformational state of a peptide in a flat membrane.

It works by comparing the position of pre-defined residues (peptide dependent) with
respect to the middle bilayer plane.

This script is intended for small SA systems containing a planar bilayer (the z 
coordinate is used) and only one peptide.

The output is coded as follows:
- penetratin: 0 surfacic, 1 TM
- transportan: 0 surfacic, 1 TM, 2 TM*, 3 U-shape

In case the --model option is used with --transportan, surfacic and TM (0 and 1) are the
only two possible conformations.

[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - numpy
 - scipy

[ NOTES ]

Leaflet finding method:
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use a 15 Angstrom cutoff value:
     -> '--leaflet 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflet large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the gro file supplied in
    order to get a meaningful outcome.
  
[ USAGE ]
	
Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		1	: process every t-frames

Peptide selection
-----------------------------------------------------
--penetratin		: to select penetratin
--transportan		: to select transportan
--model			: to specify helical model (only useful for transportan)

Leaflet identification  
-----------------------------------------------------
--leaflets	optimise: leaflet identification technique
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[10000000000000], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[1], type=int, help=argparse.SUPPRESS)

#peptide selection
parser.add_argument('--penetratin', dest='penetratin', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--transportan', dest='transportan', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--model', dest='model', action='store_true', help=argparse.SUPPRESS)

#lipids identification options
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

#parse user inputs
#-----------------
args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start = args.t_start[0]
args.t_end = args.t_end[0]
args.frames_dt = args.frames_dt[0]

#lipids identification options
args.cutoff_leaflet = args.cutoff_leaflet[0]

#process options
#---------------

global lipids_ff_nb
lipids_ff_nb = 0
	
#leaflet identification
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should be a number or 'large', see note."
		sys.exit(1)

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy as sp
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================

if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)

if args.xtcfilename == "no":
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

if args.penetratin and args.transportan:
	print "Error: you can't specify both --penetratin and --transportan."
	sys.exit(1)
if not args.penetratin and not args.transportan:
	print "Error: you need to specify either --penetratin or --transportan."
	sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "pep_conf_" + args.grofilename[:-4]
	else:
		args.output_folder = "pep_conf_" + args.xtcfilename[:-4]

if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	os.mkdir(args.output_folder)
	
	#create log
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/pep_conformation.log'
	output_log=open(filename_log, 'w')		
	output_log.write("[ff_contacts_sizes v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python pep_conformation.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def set_lipids_beads():

	global leaflet_sele_string
	leaflet_sele_string = "name PO4 or name PO3 or name B1A"

	return
def load_MDA_universe():												#DONE
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global nb_frames_to_process
	global f_start
	f_start = 0
	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		frames_to_write = [True]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		U = Universe(args.grofilename, args.xtcfilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes
		U.trajectory.rewind()
		#sanity check
		if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
			print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorter than the starting stime specified (" + str(args.t_start) + "ns)."
			sys.exit(1)
		if U.trajectory.numframes < args.frames_dt:
			print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

		#create list of index of frames to process
		if args.t_start > 0:
			for ts in U.trajectory:
				progress = '\r -skipping frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
				sys.stdout.flush()
				sys.stdout.write(progress)
				if ts.time/float(1000) > args.t_start:
					f_start = ts.frame-1
					break
			print ''
		if (nb_frames_xtc - f_start)%args.frames_dt == 0:
			tmp_offset = 0
		else:
			tmp_offset = 1
		frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(nb_frames_xtc - f_start)//args.frames_dt+tmp_offset))
		nb_frames_to_process = len(frames_to_process)
				
	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	return
def identify_proteins():	
	print "\nIdentifying proteins..."
	
	#import modules
	global nx
	import networkx as nx

	#declare variables
	global proteins_nb
	global proteins_sele
	global proteins_sele_string
	global proteins_sele_string_VMD
	global proteins_boundaries
	global proteins_nb_atoms
	global nb_atom_per_protein	
	global b1
	global b2	
	global b3
	proteins_nb = 0
	proteins_sele = {}
	proteins_sele_string = {}
	proteins_sele_string_VMD = {}
	proteins_boundaries = {}
	
	#check for protein presence
	if U.selectAtoms("protein").numberOfAtoms() == 0:
		print "Error: no protein detected."
		sys.exit(1)
	
	#declare local variables
	proteins_ca_nb = {}
	proteins_ca_nmax = 0
	proteins_ca_group = {}
	proteins_boundaries = {}

	#retrieve 1st atom info
	proteins_sele["all"] = U.selectAtoms("protein")
	proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
	prec_resnum = proteins_sele["all"][0].resnum
	prec_segid = proteins_sele["all"][0].segid
	prec_atnum = proteins_sele["all"][0].number+1
	prev_atnum = proteins_sele["all"][0].number+1						#atom corresponding to the beginning of the current protein
	#browse following atoms
	for a in proteins_sele["all"][1:]:
		delta_res = a.resnum-prec_resnum
		delta_atm = a.number+1-prec_atnum
		if delta_res < 0 or a.segid != prec_segid or delta_atm > 1:
			proteins_boundaries[proteins_nb] = [prev_atnum,prec_atnum]
			proteins_nb += 1
			prev_atnum = a.number + 1
		prec_resnum = a.resnum
		prec_atnum = a.number + 1
		prec_segid = a.segid		
	#add last protein section
	if prev_atnum < proteins_sele["all"][proteins_nb_atoms-1].number:
		proteins_boundaries[proteins_nb] = [prev_atnum,proteins_sele["all"][proteins_nb_atoms-1].number+1]
		proteins_nb += 1
	
	#display results
	print " -protein found:", proteins_nb

	#create protein selections and save into a txt file
	filename_sele=os.getcwd() + '/' + str(args.output_folder) + '/proteins.sele'
	for p_index in range(0, proteins_nb):
		progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
		sys.stdout.flush()
		sys.stdout.write(progress)
		proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
		proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
		proteins_sele[p_index] = U.selectAtoms(proteins_sele_string[p_index])
	nb_atom_per_protein = proteins_sele[0].numberOfAtoms()

	#sanity check
	if proteins_nb > 1:
		print "Error: more than 1 peptide detected."
		sys.exit(1)
	if args.penetratin and nb_atom_per_protein != 47:
		print "Error: --penetratin selected but " + str(nb_atom_per_protein) + " atoms found instead of 47."
		sys.exit(1)	
	if args.transportan and nb_atom_per_protein != 55:
		print "Error: --transportan selected but " + str(nb_atom_per_protein) + " atoms found instead of 55."
		sys.exit(1)

	#create beads selection
	if args.transportan:
		b1 = proteins_sele[0].selectAtoms("bynum 1")
		b2 = proteins_sele[0].selectAtoms("bynum 31")
		b3 = proteins_sele[0].selectAtoms("bynum 49")
	if args.penetratin:
		b1 = proteins_sele[0].selectAtoms("bynum 1")
		b2 = proteins_sele[0].selectAtoms("bynum 45")

	print ""

	return
def identify_leaflets():
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
		else:
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if np.shape(L.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			leaflet_sele["upper"]["all species"] = L.group(0)
			leaflet_sele["lower"]["all species"] = L.group(1)
		else:
			leaflet_sele["upper"]["all species"] = L.group(1)
			leaflet_sele["lower"]["all species"] = L.group(0)
		leaflet_sele["both"]["all species"] = leaflet_sele["lower"]["all species"] + leaflet_sele["upper"]["all species"]
		if np.shape(L.groups())[0] == 2:
			print " -found 2 leaflets: ", leaflet_sele["upper"]["all species"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"]["all species"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, np.shape(L.groups())[0]):
				other_lipids += L.group(g).numberOfResidues()
			print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"]["all species"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"]["all species"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	#use cof:
	else:
		leaflet_sele["both"]["all species"] = U.selectAtoms(leaflet_sele_string)
		tmp_lipids_avg_z = leaflet_sele["both"]["all species"].centerOfGeometry()[2]
		leaflet_sele["upper"]["all species"] = leaflet_sele["both"]["all species"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		leaflet_sele["lower"]["all species"] = leaflet_sele["both"]["all species"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))
		print " -found 2 leaflets: ", leaflet_sele["upper"]["all species"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"]["all species"].numberOfResidues(), '(lower) lipids'
	
	#store full selections
	for l in ["lower","upper","both"]:
		leaflet_sele_atoms[l]["all species"] = leaflet_sele[l]["all species"].residues.atoms
	
	return

#=========================================================================================
# data structures
#=========================================================================================

def data_struct_time():													

	global frames_nb
	global frames_time
	frames_nb = np.zeros(nb_frames_to_process)
	frames_time = np.zeros(nb_frames_to_process)

	return
def data_pep_conf():

	global pep_conf
	pep_conf = {}
			
	return

#=========================================================================================
# core functions
#=========================================================================================

def detect_pep_conf(f_time):

	#find middle membrane plane
	tmp_zu = leaflet_sele["upper"]["all species"].centerOfGeometry()[2]
	tmp_zl = leaflet_sele["lower"]["all species"].centerOfGeometry()[2]
	z_median = tmp_zl + (tmp_zu - tmp_zl)/float(2)

	#case: penetratin
	#----------------
	if args.penetratin:
		b1_pos = b1.centroid()[2] - z_median
		b2_pos = b2.centroid()[2] - z_median
		if np.sign(b1_pos) == np.sign(b2_pos):
			pep_conf[f_time] = 0
		else:
			pep_conf[f_time] = 1
	
	#case: transportan
	#-----------------
	if args.transportan:
		if args.model:
			b1_pos = b1.centroid()[2] - z_median
			b3_pos = b3.centroid()[2] - z_median
			if np.sign(b1_pos) == np.sign(b3_pos):
				pep_conf[f_time] = 0
			else:
				pep_conf[f_time] = 1
		else:
			b1_pos = b1.centroid()[2] - z_median
			b2_pos = b2.centroid()[2] - z_median
			b3_pos = b3.centroid()[2] - z_median
			if np.sign(b1_pos) == np.sign(b2_pos) and np.sign(b1_pos) == np.sign(b3_pos):
				pep_conf[f_time] = 0
			elif np.sign(b1_pos) == np.sign(b2_pos):
				pep_conf[f_time] = 1
			elif np.sign(b2_pos) == np.sign(b3_pos):
				pep_conf[f_time] = 2
			elif np.sign(b1_pos) == np.sign(b3_pos):
				pep_conf[f_time] = 3

	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():
	
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/' + str(args.output_folder) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("# [pep_conformation v" + str(version_nb) + "]\n")
	output_xvg.write("# The data is organised as follows:\n")
	if args.penetratin:
		output_xvg.write("#  - 0: surfacic\n")
		output_xvg.write("#  - 1: TM\n")
		output_xvg.write("@ title \"Evolution of penetratin conformation\"\n")
	else:
		output_xvg.write("#  - 0: surfacic\n")
		output_xvg.write("#  - 1: TM\n")
		output_xvg.write("#  - 2: TM*\n")
		output_xvg.write("#  - 3: U-shape\n")
		output_xvg.write("@ title \"Evolution of transportan conformation\"\n")
	output_xvg.write("@ xaxis label \"time\"\n")
	output_xvg.write("@ yaxis label \"confomration code\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 1\n")
	output_xvg.write("@ s0 legend \"conformation\"\n")
	for f_index in range(0,len(frames_time)):
		f_time = frames_time[f_index]
		results = str(f_time) + "	" + str(pep_conf[f_time])
		output_xvg.write(results + "\n")
	output_xvg.close()

	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
#process inputs
#=========================================================================================
#data loading
set_lipids_beads()
load_MDA_universe()
identify_proteins()
identify_leaflets()

#create data structures
data_struct_time()
data_pep_conf()

#=========================================================================================
# generate data
#=========================================================================================
print "\nDetecting peptide conformation..."

#case: gro file
if args.xtcfilename == "no":	
	detect_pep_conf(0)

#case: xtc file
else:
	for f_index in range(0,nb_frames_to_process):
		ts = U.trajectory[frames_to_process[f_index]]
		if ts.time/float(1000) > args.t_end:
			break
		progress = '\r -processing frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
				
		#frame properties
		f_time = ts.time/float(1000)
		f_nb = ts.frame
		frames_nb[f_index] = f_nb
		frames_time[f_index] = f_time
		box_dim = U.trajectory.ts.dimensions
		
		#process ff lipids
		detect_pep_conf(f_time)
		
	print ''
		
#=========================================================================================
# produce outputs
#=========================================================================================
print "\nWriting outputs..."
write_xvg()
					
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
