from typing import List, Dict
import os, h5py, json, time, math, re
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import tqdm

from dataset.utility import (
	load_disorder_dbs,
	find_disorder_regions )


# Define common paths to be used downstream.
version = "v_19"
db_path = "../database"
v19_dir_path = os.path.join( db_path, "v_19/" )
merged_complexes_dir = os.path.join(
	v19_dir_path, "merged_binary_complexes" )
train_file_path = os.path.join(
	v19_dir_path, "prot_1-2_train_v_19.csv" )
ood_file_path = os.path.join(
	v19_dir_path, "prot_1-2_test_v_19.csv" )
disprot_path = os.path.join( db_path, "input_files/DisProt.csv" )
ideal_path = os.path.join( db_path, "input_files/IDEAL.csv" )
mobidb_path = os.path.join( db_path, "input_files/MobiDB.csv" )


def read_file( file_name: str ):
	with open( file_name, "r" ) as f:
		data = f.readlines()[0].split( "," )
	if data[-1] == "":
		return data[:-1]
	else:
		return data


def write_to_file( plist: List, file_name: str ):
	with open( file_name, "w" ) as w:
		w.writelines( ",".join( plist ) )


def split_entry_id( entry_id: str, return_pos: bool = False ):
	uni_id1, uni_id2 = entry_id.split( "--" )
	uni_id1, start1, end1 = uni_id1.split( ":" )
	uni_id2, copy_num = uni_id2.split( "_" )
	uni_id2, start2, end2 = uni_id2.split( ":" )
	if return_pos:
		return uni_id1, int( start1 ), int( end1 ), uni_id2, int( start2 ), int( end2 ), copy_num
	else:
		return uni_id1, uni_id2, copy_num

################################################################################
def get_dor_ddr_complexes():
	"""
	Given the OOD dataset, obtain segregate
		the DOR and DDR complexes.
	For this we use the summed contact map for the complex.
	"""
	logs = {}
	for name, file in zip( ["train", "ood"], [train_file_path, ood_file_path] ):
		print( f"DOR/DDR complexes from {name} set..." )
		complexes_v19 = read_file( file_name = file )
		dor_complexes, ddr_complexes = [], []
		if name == "ood":
			dor_complexes_file = os.path.join( f"./ooddor_{version}.csv" )
			ddr_complexes_file = os.path.join( f"./oodddr_{version}.csv" )

		for i, entry_id in enumerate( complexes_v19 ):
			uni_id1, uni_id2, copy_num = split_entry_id( entry_id = entry_id )

			# These Uniprot pairs are sequence redundant with PDB70 at 20% seq identity.
			# 	Ignoring these from evaluation.
			if f"{uni_id1}--{uni_id2}_{copy_num}" in [
									"P0DTC9--P0DTD1_2", "Q96PU5--Q96PU5_0", "P0AG11--P0AG11_4", 
									"Q9IK92--Q9IK91_0", "Q16236--O15525_0", "P12023--P12023_0",
									"O85041--O85043_0", "P25024--P10145_0"]:
				continue
			# Ignoring this entry, as AF2-multimer crashed for this.
			if entry_id == "P0DTD1:1743:1808--P0DTD1:1565:1641_1":
				continue

			hf = h5py.File( os.path.join(
				merged_complexes_dir,
				f"{uni_id1}--{uni_id2}_{copy_num}.h5" ), "r" )

			summed_cmap = np.array( hf["summed_cmap"] )
			total_conformers = total_conformers = int( np.array( hf["conformers"] ) )

			contacts_idx = np.where( summed_cmap > 0 )
			unique_contacts = summed_cmap[contacts_idx]

			if np.all( unique_contacts == total_conformers ):
				dor_complexes.append( entry_id )
			else:
				ddr_complexes.append( entry_id )

		print( f"DOR: {len( dor_complexes )} \t DDR: {len( ddr_complexes )}" )
		logs[name] = {
			"dor": len( dor_complexes ),
			"ddr": len( ddr_complexes )
		}
		if name == "ood":
			write_to_file( dor_complexes, dor_complexes_file )
			write_to_file( ddr_complexes, ddr_complexes_file )
	return logs


################################################################################
disprot, ideal, mobidb = load_disorder_dbs(
	disprot_path = disprot_path,
	ideal_path = ideal_path,
	mobidb_path = mobidb_path )

def get_disordered_regions( uni_id: str ):
	disordered_residues = find_disorder_regions( disprot = disprot,
												ideal = ideal,
												mobidb = mobidb,
												uni_ids = [uni_id],
												min_len = 1, return_ids = False )
	return disordered_residues


def get_frac_disordered( uni_id: str, start: int, end: int ):
	disordered_regions = get_disordered_regions( uni_id )

	uni_res = np.arange( start, end + 1, 1 )
	disorder_uni_res = []
	for reg in disordered_regions:
		overlap = set( reg ).intersection( set( uni_res ) )
		disorder_uni_res += list( overlap )

	total_uni = end - start + 1

	frac_disorder = len( disorder_uni_res )/total_uni

	return frac_disorder

################################################################################
def get_full_idr_complexes():
	"""
	Given the OOD dataset, obtain complexes that contain 100%
		disordered residues.
	"""
	logs = {}
	for name, file in zip( ["train", "ood"], [train_file_path, ood_file_path] ):
		print( f"100% IDR complexes in {name} set..." )
		complexes_v19 = read_file( file_name = file )

		full_idr_complexes = []
		if name == "ood":
			full_idr_complexes_file = os.path.join( f"./oodidr_{version}.csv" )

		for i, entry_id in enumerate( complexes_v19 ):
			uni_id1, s1, e1, uni_id2, s2, e2, copy_num = split_entry_id(
				entry_id = entry_id,
				return_pos = True )

			# These Uniprot pairs are sequence redundant with PDB70 at 20% seq identity.
			# 	Ignoring these from evaluation.
			if f"{uni_id1}--{uni_id2}_{copy_num}" in [
									"P0DTC9--P0DTD1_2", "Q96PU5--Q96PU5_0", "P0AG11--P0AG11_4", 
									"Q9IK92--Q9IK91_0", "Q16236--O15525_0", "P12023--P12023_0",
									"O85041--O85043_0", "P25024--P10145_0"]:
				continue
			# Ignoring this entry, as AF2-multimer crashed for this.
			if entry_id == "P0DTD1:1743:1808--P0DTD1:1565:1641_1":
				continue

			frac_disorder = get_frac_disordered(
				uni_id = uni_id1,
				start = s1,
				end = e1 )
			if frac_disorder == 1:
				full_idr_complexes.append( entry_id )

		print( f"Full IDR complexes = {len( full_idr_complexes )}" )

		logs[name] = {"total": len( full_idr_complexes )}
		if name == "ood":
			write_to_file( full_idr_complexes, full_idr_complexes_file )
	return logs


################################################################################
################################################################################
logs1 = get_dor_ddr_complexes()

logs2 = get_full_idr_complexes()

with open( "Logs_side_analysis.txt", "w" ) as w:
	w.writelines( "DOR vs DDR complexes\n" )
	w.writelines( f"\tTrain: DOR ({logs1['train']['dor']}) |  DDR ({logs1['train']['ddr']})\n" )
	w.writelines( f"\tOOD: DOR ({logs1['ood']['dor']}) |  DDR ({logs1['ood']['ddr']})\n" )
	w.writelines( "\n" )
	w.writelines( "Complexes with 100% disordered residues in prot1\n" )
	w.writelines( f"\tTrain: ({logs2['train']['total']}) |  OOD ({logs2['ood']['total']})\n" )
