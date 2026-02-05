"""
######## ------>"May the Force serve u well..." <------#########
Script to obtain Disobind or Disobind+AF2 predictions for a user.
Inputs:
	Disobind needs the UniProt IDs for both the proteins and the protein 1 (1st in the pair) must be an IDR.
	The input must be provided in a csv file formatas shown below:
	For Disobind only
		Uni_ID1,start_res1,end_res1,Uni_ID2,start_res2
	For Disobind+AF2
		Uni_ID1,start_res1,end_res1,Uni_ID2,start_res2,end_res2,af2_struct_file,af2_json_file,chain1,chain2,offset1,offset2

Outputs:
	Disobind, AF2, AF2+Disobind predictions for all tasks and all CG resolutions.
"""

############# One above all #############
##-------------------------------------##
import math, os, subprocess, argparse, json, time, pickle as pkl
import warnings
import numpy as np
import pandas as pd
from omegaconf import OmegaConf
from multiprocessing import Pool
import tqdm
import Bio
from Bio.PDB import PDBParser, MMCIFParser


import torch
from torch import nn

from typing import List, Dict

import dataset
from dataset.from_APIs_with_love import get_uniprot_seq
from dataset.utility import get_contact_map
from dataset.create_input_embeddings import Embeddings
from src.models.get_model import get_model
from src.utils import prepare_input
from analysis.params import parameter_files


class Disobind():
	def __init__( self, args ):
		"""
		Get Disobind or Disobind+AF2 predictions for an input pairs of proteins,
			at least one of which is an IDR.
		"""
		# Input file containing the prot1/2 entry_ids.
		self.input_file = args.f
		# No. of CPU cores for parallelism.
		self.cores = args.c
		# Predict cmaps.
		self.predict_cmap = args.cm
		# Coarse grained resolutions.
		self.required_cg = args.cg
		# Name for output directory.
		self.output_dir = args.o

		# Embedding type to be used for prediction.
		self.embedding_type = "T5"
		# Use global/local embeddings.
		self.scope = "global"
		# Device to be used for running Disobind.
		self.device = args.d
		# Contact probability threshold.
		self.threshold = 0.5
		# Distance cutoff for defining a contact.
		self.dist_threshold = 8
		# pLDDT cutoff.
		self.plddt_threshold = 70
		# PAE cutoff.
		self.pae_threshold = 5
		# Batch size for obtaining embeddings.
		self.batch_size = 200
		# Objective settings to be used for prediction.
		self.objective = ["", "", "", ""]
		# Load a dict storing paths for each model.
		self.parameters = parameter_files( 19 )
		# Dict to store predictions for all tasks.
		self.predictions = {}

		# Dict to store Uniprot sequences.
		self.uniprot_seq = {}

		# Dir containing the Disobind models.
		self.model_dir = "./models/"
		# Filename to store predictions.
		self.output_filename = "Predictions"

		if os.path.exists( self.output_dir ):
			reply = input( "Output directory already exists. Abort? (Y/n)\t" )
			if reply == "Y":
				exit()
		else:
			os.makedirs( self.output_dir, exist_ok = True )
		
		# Absolute path to the analysis dir.
		self.abs_path = os.getcwd()
		# Path for the FASTA file for prot1/2 sequences.
		self.fasta_file = f"{self.abs_path}/{self.output_dir}/p1_p2_test.fasta"
		# Path fo rthe h5 file for prot1/2 embeddings.
		self.emb_file = f"{self.abs_path}/{self.output_dir}/p1_p2_test.h5"
		# Uniprot seq file.
		self.uni_seq_file = f"{self.abs_path}/{self.output_dir}/UniProt_seq.json"


	def forward( self ):
		"""
		Get the Uniprot sequences for all unique Uniprot IDs amongst prot1/2.
			For test mode,
				Use the downloaded Uniprot sequences.
			For predict mode,
				Download the sequences.
		"""
		entry_ids, af_dict = self.read_csv_input()

		prot_pairs = self.process_input_pairs( entry_ids )

		self.get_predictions( prot_pairs, af_dict )

		np.save( f"./{self.output_dir}/{self.output_filename}.npy", self.predictions )



	def get_required_tasks( self ):
		"""
		Get all tasks required by the user.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		tasks --> list of user required tasks.		
		"""
		tasks = []

		if self.required_cg not in [0, 1, 5, 10]:
			raise ValueError( "Incorrect coarse-grained resolution specified. \n" + 
								"Choose from [0, 1, 5, 10]..." )

		if self.required_cg in [5, 10]:
			warnings.warn( f"\n\033[1m Using coarse-grained model with kernel size = {self.required_cg}.\033[0m\n"+
							f"\033[1m C-terminal residues would be lost if protein length\033[0m"+
							"\033[1m is not a multiple of kernel size...\n\033[0m" )

		for obj in ["interaction", "interface"]:
			if obj == "interaction" and not self.predict_cmap:
				continue
			
			for cg in [1, 5, 10]:
				if self.required_cg == 0:
					tasks.append( f"{obj}_{cg}" )

				elif cg == self.required_cg:
					tasks.append( f"{obj}_{cg}" )

		print( f"\033[1mRunning Disobind for the followinng tasks: \n{tasks}\033[0m" )
		return tasks


	def get_predictions( self, prot_pairs, af_dict ):
		"""

		Input:
		----------
		prot_pairs --> dict containig list of fragment pairs for all protein pairs.

		Returns:
		----------
		None
		"""
		total_pairs = len( prot_pairs )

		required_tasks = self.get_required_tasks()


		for start in np.arange( 0, total_pairs, self.batch_size ):
			t_start = time.time()
			if start + self.batch_size >= total_pairs:
				end = total_pairs
			else:
				end = start + self.batch_size
			
			print( f"\n{start}:{end}-----------------------------------------" )
			batch = prot_pairs[start:end]
			
			print( "Creating global embeddings for the input sequences..." )
			self.create_embeddings( batch )

			os.chdir( self.abs_path )

			print( "Running predictions..." )
			batch_preds = self.predict( required_tasks, af_dict )
			self.predictions.update( batch_preds )

			t_end = time.time()
			print( f"Time taken for batch {start}-{end} = {( t_end - t_start )/60} minutes\n" )

			subprocess.call( ["rm", f"{self.emb_file}", f"{self.fasta_file}"] )


###################################################################################
##-------------------------------------------------------------------------------##
###################################################################################
	def read_csv_input( self ):
		"""
		Read from a csv file containing prot1 and prot2 info as:
			Uni_ID1,start1,end1,Uni_ID2,start2,end2
		To combine AF2 and Disobind predictions, provide the path to the AF2 model and .pkl file as:
			Uni_ID1,start1,end1,Uni_ID2,start2,end2,model_file_path,pkl_file_path,chain1,chain2,offset1,offset2

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		entry_ids --> list of entry_ids for al binary complexes.
		"""
		entry_ids = []
		af_dict = {}
		with open( self.input_file, "r" ) as f:
			input_pairs = f.readlines()

		# Validate if the input is in the correct format:
		for pair in input_pairs:
			pair = pair.strip()
			# Ignore empty string if present in the input.
			if len( pair.split( "," ) ) == 6:
				uni_id1, start1, end1, uni_id2, start2, end2 = pair.split( "," )
				af_struct_file, af_json_file = None, None
				chain1, chain2, offset1, offset2 = None, None, None, None

			elif len( pair.split( "," ) ) == 12:
				( uni_id1, start1, end1, uni_id2, start2, end2,
					af_struct_file, af_json_file, chain1, chain2, offset1, offset2 ) = pair.split( "," )

			else:
				raise ValueError( f"Incorrect input format..." )

			entry_id = f"{uni_id1}:{start1}:{end1}--{uni_id2}:{start2}:{end2}_0"
			entry_ids.append( entry_id )
			af_dict[entry_id] = {
						"struct_file": af_struct_file,
						"json_file": af_json_file,
						"required_chains": {
										"chains": [chain1, chain2],
										"offsets": [offset1, offset2]
										 }
									}

		return entry_ids, af_dict


###################################################################################
##-------------------------------------------------------------------------------##
###################################################################################
	def parallelize_uni_seq_download( self, uni_id ):
		"""
		Obtain all unique UniProt IDs in the provided input.
		Download all unique UniProt sequences.

		Input:
		----------
		uni_id --> (str) UniProt ID protein.

		Returns:
		----------
		uni_id --> (str) UniProt ID protein.
		seq --> (str) UniProt seq for protein.
		"""
		seq = get_uniprot_seq( uni_id, max_trials = 10, wait_time = 20, return_id = False )

		return uni_id, seq


	def get_unique_uni_ids( self, entry_ids ):
		"""
		Given a list of all protein pairs, get all the unique Uniprot IDs.
		"""
		unique_uni_ids = []
		for head in entry_ids:
			head1, head2 = head.split( "--" )
			uni_id1 = head1.split( ":" )[0]
			uni_id2 = head2.split( ":" )[0]

			if uni_id1 not in unique_uni_ids:
				unique_uni_ids.append( uni_id1 )
			if uni_id2 not in unique_uni_ids:
				unique_uni_ids.append( uni_id2 )

		return unique_uni_ids



	def download_uniprot_seq( self, entry_ids ):
		"""
		Obtain all unique UniProt IDs in the provided input.
		Download all unique UniProt sequences.

		Input:
		----------
		entry_ids --> list of entry_id for a binary complexes.

		Returns:
		----------
		None
		"""
		unique_uni_ids = self.get_unique_uni_ids( entry_ids )

		with Pool( self.cores ) as p:
			for result in tqdm.tqdm( p.imap_unordered( 
													self.parallelize_uni_seq_download, 
													unique_uni_ids ), 
										total = len( unique_uni_ids ) ):
				uni_id, seq = result

				if len( seq ) != 0:
					self.uniprot_seq[uni_id] = seq
				else:
					raise RuntimeError( f"Unable to download seq for Uniprot ID: {uni_id}. Please retry..." )

		print( "Unique Uniprot IDs = ", len( unique_uni_ids ) )
		print( "Total Uniprot sequences obtained = ", len( self.uniprot_seq ) )


###################################################################################
##-------------------------------------------------------------------------------##
###################################################################################
	def process_input_pairs( self, entry_ids ):
		"""
		For all input pairs, convert them into a uniform format:
			UniID1:start1:end1--UniID2:start2:end2_num

		Input:
		----------
		entry_ids --> list of entry_id for a binary complexes.

		Returns:
		----------
		processed_pairs --> list of protein 1/2 fragment pairs to obtain Disobind predictions for.
		"""
		print( "\nDownloading UniProt sequences..." )
		if not os.path.exists( self.uni_seq_file ):
			self.download_uniprot_seq( entry_ids )

			with open( self.uni_seq_file, "w" ) as w:
				json.dump( self.uniprot_seq, w )
		else:
			with open( self.uni_seq_file, "r" ) as f:
				self.uniprot_seq = json.load( f )

		prot_pairs = []
		for head in entry_ids:
			head1, head2 = head.split( "--" )
			uni_id1, start1, end1 = head1.split( ":" )
			uni_id2, start2, end2 = head2.split( ":" )

			if uni_id1 not in self.uniprot_seq.keys() or uni_id2 not in self.uniprot_seq.keys():
				print( f"Skipping the input pair -- {head}" )		

			else:
				prot_pairs.append( head )
		
		return prot_pairs



###################################################################################
##-------------------------------------------------------------------------------##
###################################################################################
	def create_embeddings( self, entry_ids ):
		"""
		Use the Shredder() class to:
			Create fasta files and get embeddings.

		Input:
		----------
		entry_ids --> list of entry_id for a binary complexes.

		Returns:
		----------
		None
		"""
		self.prot1_emb, self.prot2_emb = Embeddings( scope = self.scope, embedding_type = self.embedding_type, 
													uniprot_seq = self.uniprot_seq, base_path = self.output_dir, 
													fasta_file = self.fasta_file, 
													emb_file = self.emb_file, 
													headers = entry_ids, 
													load_cmap = False,
													eval_ = True
													 ).initialize( return_emb = True )



	def load_model( self, model_ver, params_file ):
		"""
		Load pre-trained model in evaluation mode for making predictions.

		Input:
		----------
		model_ver --> version of the model to be loaded.
		params_file --> model parameters file.

		Returns:
		----------
		model --> trained model in eval() mode.
		"""
		model_config = f"./params/Model_config_{model_ver}.yml"
		model_config = OmegaConf.load( model_config )
		model_config = model_config.conf.model_params

		model_ver = model_ver.split( "_" )
		# Model name.
		mod = "_".join( model_ver[:-1] )
		# model version.
		ver = model_ver[-1]
		# Load Disobind model.
		model = get_model( model_config ).to( self.device )
		model.load_state_dict( 
							torch.load( f"{self.model_dir}{mod}/Version_{ver}/{params_file}.pth", 
										map_location = self.device )
							 )
		model.eval()

		return model



	def apply_settings( self, obj, cg ):
		"""
		Load models and modify the self.objective depending on the model type.
			self.objective includes 4 hparams to be set:
			0. objective: interaction or interface
			1. bin_size: level of coarse graining.
			2. bin_input: if bin_size > 1, average the input embeddings.
			3. single_output: if True, use single output prediction task.

		Input:
		----------
		key1 --> specifies the model to be used.
		key2 --> specifies the coarse graining to be used.

		Returns:
		----------
		model --> model corresponding to the specified settings. 
		"""
		print( f"\nLoading model for {obj}: CG = {cg}..." )
		model = self.load_model( *self.parameters[obj][cg] )

		if "interaction" in obj:
			self.objective[0] = "interaction"
			if "mo" in obj:
				self.objective[3] = False
			elif "so" in obj:
				self.objective[3] = True
		
		elif "interface" in obj:
			self.objective[0] = "interface"
			self.objective[3] = False

		_, cg = cg.split( "_" )
		self.objective[1] = int( cg )
		if int( cg ) != 1:
			self.objective[2] = True
		else:
			self.objective[2] = False
		return model



	def get_input_tensors( self, prot1, prot2 ):
		"""
		Prepare prot1, prot2, target tensors according to the objective.
			If self.pad is True:
				Pad the prot1/2 and cmap arrys.
				Create a binary target_mask.
		Convert the input arrays to torch.tensor and load on the specified device.
		**Note: If not using padding, terminal residues might be lost upon coarse graining 
				if protein length not divisible by bin_size.

		Input:
		----------
		prot1 --> (np.array) prot1 embedding of dimension [L1, C].
		prot1 --> (np.array) prot2 embedding of dimension [L2, C].

		Returns:
		----------
		prot1 --> (torch.tensor) prot1 embedding of dimension [N, L1, C].
		prot2 --> (torch.tensor) prot2 embedding of dimension [N, L2, C].
		target --> (torch.tensor) contact map of dimension [N, L1, L2].
			If padding, L1 == L2; N = 1.
		target_mask --> (torch.tensor) binary mask of dimension [N, L1, L2].
		interaction_mask --> (torch.tensor) binary mask to
							ignore pads in interface_block [N, L1, L2].
		eff_len --> effective length of prot1/2 (eff_len) post coarse graining.
		"""
		num_res1 = prot1.shape[0]
		num_res2 = prot2.shape[0]

		if self.objective[1] == 1:
			eff_len = [num_res1, num_res2]
		else:
			cg = self.objective[1]
			eff_len = [num_res1//cg, num_res2//cg]

		target_mask = np.ones( ( 1, num_res1, num_res2 ) )

		target_mask = torch.from_numpy( target_mask )

		prot1 = torch.from_numpy( prot1 ).unsqueeze( 0 ).float()
		prot2 = torch.from_numpy( prot2 ).unsqueeze( 0 ).float()

		( prot1, prot2, target,
			target_mask, interaction_mask ) = prepare_input( prot1, prot2, None, 
															[True, target_mask], 
															objective = self.objective[0], 
															bin_size = self.objective[1], 
															bin_input = self.objective[2], 
															single_output = self.objective[3] )
		prot1 = prot1.to( self.device ).float()
		prot2 = prot2.to( self.device ).float()
		target_mask = target_mask.to( self.device ).float()
		interaction_mask = interaction_mask.to( self.device ).float()

		return prot1, prot2, target, target_mask, interaction_mask, eff_len


	def predict( self, required_tasks, af_dict ):
		"""
		Predict cmap for the input protein pair from all models.
		Store all predictions in a nested  dict:
			self.predictions{
				pair_id: {
						entry_id: {
							{obj}_{cg}: {
									"Disobind_preds",
									"AF2_preds",
									"Diso+AF2_preds",
									"Final_Disobind preds",
									"Final_AF2 preds",
									"Final_Diso+AF2 preds"
							}
						}
				}
			}
			header --> identifier for input pair as specified in input file.
			cg --> coarse grainng level.
		In test mode, also plot the distribution of contact counts, 
			prot1/2 lengths for the OOD set.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		None
		"""
		predictions = {}
		
		# For all entries in OOD set.
		for idx, entry_id in enumerate( self.prot1_emb.keys() ):
			head1, head2 = entry_id.split( "--" )
			head2, num = head2.split( "_" )
			uni_id1, start1, end1 = head1.split( ":" )
			uni_id2, start2, end2 = head2.split( ":" )

			# Pair does not have the _0 at the end.
			pair_id = f"{uni_id1}::--{uni_id2}::"
			if pair_id not in predictions.keys():
				predictions[pair_id] = {}

			predictions[pair_id][entry_id] = {}
			
			prot1_emb = np.array( self.prot1_emb[entry_id] )
			prot2_emb = np.array( self.prot2_emb[entry_id] )

			# For all objectives (interaction, interface).
			for obj in self.parameters.keys():
				
				# For all coarse grainings.
				for cg in self.parameters[obj].keys():
					if f"{obj}_{cg.split( '_' )[1]}" not in required_tasks:
						continue
					
					model = self.apply_settings( obj, cg )

					( prot1, prot2, target,
						target_mask, interaction_mask,
						eff_len ) = self.get_input_tensors( prot1 = prot1_emb, 
															prot2 = prot2_emb )
					
					# Get model predictions.
					with torch.no_grad():
						output = model( prot1, prot2, interaction_mask )
						
						# output = output*target_mask
					
					output = output.detach().cpu().numpy()
					# target_mask = target_mask.detach().cpu().numpy()

					# Shape must be [m, n] --> m and n are padded lengths of prot1/2.
					if "interaction" in self.objective[0]:
						output = output.reshape( eff_len )
					# Shape must be [m+n, 1].
					elif "interface" in self.objective[0]:
						output = output.reshape( sum( eff_len ), 1 )

					_, cg = cg.split( "_" )

					# Get AF2 pred for the entry.
					model_file = af_dict[entry_id]["struct_file"]
					pae_file = af_dict[entry_id]["json_file"]
					required_chains = af_dict[entry_id]["required_chains"]
					# If AF2 input is provided.
					if model_file is not None:
						af_obj = AfPrediction( struct_file_path = model_file, data_file_path = pae_file,
												required_chains = required_chains )

						af2_pred = af_obj.get_confident_interactions( [int( start1 ), int( end1 )],
																	[int( start2 ), int( end2 )] )

						af2_pred = self.process_af2_pred( af2_pred )

						# Get Disobind+AF2 output.
						m, n = output.shape
						# Here Disobind output is not binary. But doesn't matter as we are just taking a max.
						diso_af2 = np.stack( [output.reshape( -1 ), af2_pred.reshape( -1 )], axis = 1 )
						diso_af2 = np.max( diso_af2, axis = 1 ).reshape( m, n )

						af2_pred, df_af2 = self.extract_model_output( entry_id,
																		af2_pred,
																		eff_len,
																		"af2" )
						diso_af2, df_diso_af2 = self.extract_model_output( entry_id,
																			diso_af2,
																			eff_len,
																			"diso_af2" )

					else:
						diso_af2 = np.array( [] )
						af2_pred = np.array( [] )
						df_af2 = pd.DataFrame( {} )
						df_diso_af2 = pd.DataFrame( {} )

					output, df_diso = self.extract_model_output( entry_id, output, eff_len, "diso" )

					predictions[pair_id][entry_id][f"{obj}_{cg}"] = {
																	"Disobind": np.float32( output ),
																	"AF2": np.float32( af2_pred ),
																	"Diso+AF2": np.float32( diso_af2 ),
																	"Final_diso_preds": df_diso,
																	"Final_af2_preds": df_af2,
																	"Final_af2_diso_preds": df_diso_af2
																		}
					print( f"{idx} ------------------------------------------------------------\n" )
		return predictions


	def process_af2_pred( self, af2_pred ):
		"""
		Get contact map or interface residues from AF2 
			predicted contact map for the required CG.

		Input:
		----------
		af2_pred --> (np.array) AF2 predicted contact map.

		Returns:
		----------
		"""
		if self.objective[1] > 1:
			m = nn.MaxPool2d( kernel_size = self.objective[1], stride = self.objective[1] )
			af2_pred = torch.from_numpy( af2_pred ).float()
			af2_pred = m( af2_pred.unsqueeze( 0 ).unsqueeze( 0 ) )
			af2_pred = af2_pred.squeeze( [0, 1] ).numpy()

		if "interface" in self.objective[0]:
			# af2_pred is laready padded to [100, 100].
			p1 = np.zeros( ( af2_pred.shape[0], 1 ) )
			p2 = np.zeros( ( af2_pred.shape[1], 1 ) )

			idx = np.where( af2_pred == 1 )

			p1[idx[0]] = 1
			p2[idx[1]] = 1
			af2_pred = np.concatenate( [p1, p2], axis = 0 )

		return af2_pred


	def get_beads( self, cg, prot1, prot2 ):
		"""
		Get the prot1/2 beads corresponding to a residue (cg 1) 
			or a set of residues (cg 5/10).

		Input:
		----------
		cg --> (int) resolution of coarse-graining.
		prot1 --> (str) identifier comprising "uni_id1:start1:end1"
				start1, end1 are the UniProt positions.
		prot2 --> (str) identifier comprising "uni_id2:start2:end2"
				start2, end2 are the UniProt positions.

		Returns:
		----------
		beads1 --> list of baeds at the apprpriate resolution for prot1.
		beads2 --> list of baeds at the apprpriate resolution for prot2.
		"""
		uni_id1, start1, end1 = prot1.split( ":" )
		uni_id2, start2, end2 = prot2.split( ":" )
		start1, end1 = int( start1 ), int( end1 )
		start2, end2 = int( start2 ), int( end2 )

		len_p1 = end1 - start1 + 1
		len_p2 = end2 - start2 + 1
		
		if cg == 1:
			beads1 = np.arange( start1, end1 + 1, 1 )
			beads1 = np.array( list( map( str, beads1 ) ) )
			beads2 = np.arange( start2, end2 + 1, 1 ) 
			beads2 = np.array( list( map( str, beads2 ) ) )

			p1_cg_len = len_p1
			p2_cg_len = len_p2

		else:
			beads1, beads2 = [], []
			for s in np.arange( start1, end1 + 1, cg ):
				e = s + cg - 1 if s + cg - 1 < end1 else end1
				beads1.append( f"{s}-{e}" )
			beads1 = np.array( beads1 )
			
			for s in np.arange( start2, end2 + 1, cg ):
				e = s + cg - 1 if s + cg - 1 < end2 else end2
				beads2.append( f"{s}-{e}" )
			beads2 = np.array( beads2 )

			# H_out = ( ( H_in + 2*padding - dilation( kernel_size - 1 ) )/stride )
			# For us padding = 0, dilation = 1, and kernel_size = stride = cg.
			# if len_p1 == 100:
			# 	p1_cg_len = ( ( len_p1 - ( cg - 1 ) )//cg )
			# else:
			p1_cg_len = math.ceil( ( ( len_p1 - ( cg - 1 ) )/cg ) )
			
			# if len_p2 == 100:
			# 	p2_cg_len = ( ( len_p2 - ( cg - 1 ) )//cg )
			# else:
			p2_cg_len = math.ceil( ( ( len_p2 - ( cg - 1 ) )/cg ) )

		return p1_cg_len, beads1, p2_cg_len, beads2


	def extract_model_output( self, entry_id: str, output: np.array, eff_len: List, name: str ):
		"""
		Reshape the model output and target into the required shape 
			i.e. [L1, L2] for interaction.
				[L1+L2, 1] for interface.
		eff_len specifies the effective size post coarse graining.

		Input:
		----------
		output --> (np.array) model predictions.
		target --> (np.array) target cmap.
		eff_len -->  list containing effective lengths for prot1/2.
		name --> identifier for Disobind/AF2/AF2-Disobind output.

		Returns:
		----------
		output --> (np.array) model predictions.
		target --> (np.array) target cmap.
		"""
		obj = self.objective[0]
		cg = self.objective[1]
		prot1, prot2 = entry_id.split( "--" )
		prot2, _ = prot2.split( "_" )

		p1_cg_len, beads1, p2_cg_len, beads2 = self.get_beads( cg, prot1, prot2 )
		
		if "interaction" in obj:
			output = output[:len_p1, :len_p2]
			idx = np.where( output >= self.threshold )
			df = pd.DataFrame()
			if len( idx[0] ) != 0:
				df["Protein1"] = [prot1]*len( idx[0] )
				df["Residue1"] = beads1[idx[0]]
			else:
				df["Protein1"] = []
				df["Residue1"] = []
			if len( idx[1] ) != 0:
				df["Protein2"] = [prot2]*len( idx[1] )
				df["Residue2"] = beads2[idx[1]]
			else:
				df["Protein2"] = []
				df["Residue2"] = []

		elif "interface" in  obj:
			# prot1 residues are first.
			interface1 = output[:p1_cg_len]
			interface2 = output[p1_cg_len:]

			idx1 = np.where( interface1 >= self.threshold )[0]
			idx2 = np.where( interface2 >= self.threshold )[0]

			interface1_beads = beads1[idx1] if len( idx1 ) != 0 else []
			interface2_beads = beads2[idx2] if len( idx2 ) != 0 else []

			if len( interface1_beads ) > len( interface2_beads ):
				for _ in range( len( interface1_beads ) - len( interface2_beads ) ):
					interface2_beads = np.append( interface2_beads, "" )
			else:
				for _ in range( len( interface2_beads ) - len( interface1_beads ) ):
					interface1_beads = np.append( interface1_beads, "" )

			df = pd.DataFrame()
			df["Protein1"] = interface1_beads
			df["Protein2"] = interface2_beads
			output = np.concatenate( ( interface1, interface2 ), axis = 0 )
		df.to_csv( f"{self.abs_path}/{self.output_dir}/{name}_{entry_id}_{obj}_cg{cg}.csv" )

		return output, df


###################################################################################
##-------------------------------------------------------------------------------##
###################################################################################
class AfPrediction():
	def __init__( self, struct_file_path: str, data_file_path: str, required_chains: Dict ):
		# AF2/3 structure file path.
		self.struct_file_path = struct_file_path
		# AF2/3 structure data file path.
		self.data_file_path = data_file_path
		self.chains = required_chains["chains"]
		self.offsets = list( map( int, required_chains["offsets"] ) )

		# Distance threshold in (Angstorm) to define a contact between residue pairs.
		self.dist_threshold = 8
		# pLDDt cutoff to consider a confident prediction.
		self.plddt_cutoff = 70
		# PAE cutoff to consider a confident prediction.
		self.pae_cutoff = 5

		# Biopython Structure object.
		self.structure = self.get_structure( 
									self.get_parser()
									 )
		# Average PAE matrix.
		self.get_pae()



	def get_parser( self ):
		"""
		Get the required parser (PDB/CIF) for the input file.
		"""
		ext = os.path.splitext( self.struct_file_path )[1]

		if "pdb" in ext:
			parser = PDBParser()
		elif "cif" in ext:
			parser = MMCIFParser()
		else:
			raise ValueError( "Incorrect file format.. Suported .pdb/.cif only." )

		return parser


	def get_structure( self, parser: Bio.PDB.PDBParser ):
		"""
		Return the Biopython Structure object for the input file.
		"""
		basename = os.path.basename( self.struct_file_path )
		structure = parser.get_structure( basename, self.struct_file_path )

		return structure


	def get_chains( self ):
		"""
		A generator that yields all Chain objects from the structure.
		"""
		for chain in self.structure[0]:
			yield chain


	def get_residues( self, chain: str ):
		"""
		Get all residues in the specified chains from the structure.
		"""
		for residue in self.structure[0][chain]:
			yield residue


	def extract_perresidue_quantity( self, residue, quantity: str ): 
		"""
		Given the Biopython residue object, return the specified quantity:
			1. residue position
			2. Ca-coordinate
			3. Ca-pLDDT
		"""
		# Use Ca-atom for all other amino acids.
		rep_atom = "CA"

		if quantity == "res_pos":
			result = residue.id[1]

		elif quantity == "coords":
			coords = residue[rep_atom].coord
			result = coords
		
		elif quantity == "plddt":
			plddt = residue[rep_atom].bfactor
			result = plddt
		
		else:
			raise ValueError( f"Specified quantity: {quantity} does not exist..." )

		return result


	def get_data_dict( self ):
		"""
		Parse the AF2/3 data file.
			AF2 data file is saved as a .pkl file 
				whereas for AF3 it's stored as .json.
		"""
		ext = os.path.splitext( self.data_file_path )[1]

		if "json" in ext:
			with open( self.data_file_path, "rb" ) as f:
				data = json.load( f )
		else:
			raise ValueError( "Incorrect file format.. Suported .json only." )

		return data[0]


	def get_pae( self ):
		"""
		Return the PAE matrix from the data dict.
			AF2/3 PAE matrix is asymmetric.
			Hence, we consider the average PAE: ( PAE + PAE.T )/2.
		"""
		data = self.get_data_dict()
		# For AF2.
		if "predicted_aligned_error" in data.keys():
			pae = np.array( data["predicted_aligned_error"] )
		else:
			raise ValueError( "PAE matrix not found..." )

		self.pae = ( pae + pae.T )/2


	def add_offset( self, curr_res: int, offset: int ):
		"""
		Add the specified offset to the residue position.
		"""
		return curr_res + offset


	def is_start_residue( self, curr_res: int, prot_res: List ):
		"""
		Check if a residue is the fragment start residue.
		Assuming that the offset is already added.
		"""
		return curr_res == prot_res[0]


	def is_fragment_residue( self, curr_res: int, prot_res: List ):
		"""
		A residue is valid if it belongs to the required protein fragment.
		Assuming that the offset is already added.
		"""
		start, end = prot_res
		if curr_res >= start and curr_res <= end:
			valid = True
		else:
			valid = False
		return valid


	def get_required_residues( self, prot1_res: List, prot2_res: List ):
		"""
		The AF2 prediction may contain:
			Structure for a non-binary complex containing protein1 and protein2.
			Structure for the full length protein1 and protein2.
			Structure for a fragment containing residues specified by prot1_res and prot2_res.
			Structure for the fragment specified by prot1_res and prot2_res.
		A generator that yields Residue objects for the required fragment residues.
		"""
		res_dict = {}
		# For both the protein fragments.
		for i, prot_res in enumerate( [prot1_res, prot2_res] ):
			chain = self.chains[i]
			offset = self.offsets[i]
			for residue in self.get_residues( chain ):
				res_id = self.extract_perresidue_quantity( residue, "res_pos" )
				offset_res_id = self.add_offset( res_id, offset )
				# Select only residues part of the fragment.
				if self.is_fragment_residue( offset_res_id, prot_res ):
					if chain not in res_dict.keys():
						res_dict[chain] = [residue]
					else:
						res_dict[chain].append( residue )
		return res_dict


	def get_indices_for_pae( self, prot1_res: List, prot2_res: List ):
		"""
		Get the indices for the fragment start and end residues to 
			select required elements in PAE matrix.
		Essentially count all residues upto the start and end residue of the fragment.
		"""
		frag_index_dict = {k:[] for k in ["prot1", "prot2"]}
		# Count the no. of residue till the start residue of a fragment.
		total_length = 0
		# For all chains.
		for chain in self.get_chains():
			chain_id = chain.id
			
			for residue in self.get_residues( chain_id ):
				res_id = self.extract_perresidue_quantity( residue, "res_pos" )
				
				# For both the protein fragments.
				for i, prot_res in enumerate( [prot1_res, prot2_res] ):
					key = f"prot{i+1}"
					# Do not check residues in a chain again if the indices have been accounted.
					if frag_index_dict[key] == []:
						if self.chains[i] == chain_id:
							offset = self.offsets[i]
							offset_res_id = self.add_offset( res_id, offset )
							if self.is_start_residue( offset_res_id, prot_res ):
								frag_len = prot_res[1] - prot_res[0] + 1
								# Start index is the no. of residues till the start residue.
								start_idx = total_length
								# End index is the start index plus the fragment length.
								end_idx = start_idx + frag_len
								frag_index_dict[key] = [start_idx, end_idx]
								break
				total_length += 1

		return frag_index_dict


	def get_required_coords( self, res_dict: Dict ):
		"""
		Get the coordinates for all Ca atoms of all residues.
		"""
		coords_dict = {}
		for chain in res_dict:
			for residue in res_dict[chain]:
				coords = self.extract_perresidue_quantity( residue, "coords" )
				if chain not in coords_dict.keys():
					coords_dict[chain] = np.array( coords )
				else:
					coords_dict[chain] = np.append( coords_dict[chain], coords )

		coords_dict = {k: v.reshape( -1, 3 ) for k, v in coords_dict.items()}
		return coords_dict



	def get_required_plddt( self, res_dict: Dict ):
		"""
		Get the pLDDT score for all Ca atoms of all residues.
		"""
		plddt_dict = {}
		for chain in res_dict:
			for residue in res_dict[chain]:
				plddt = self.extract_perresidue_quantity( residue, "plddt" )
				if chain not in plddt_dict.keys():
					plddt_dict[chain] = np.array( [plddt] )
				else:
					plddt_dict[chain] = np.append( plddt_dict[chain], plddt )

		plddt_dict = {k: v.reshape( -1, 1 ) for k, v in plddt_dict.items()}
		return plddt_dict


	def get_required_pae( self, prot1_res: List, prot2_res: List ):
		"""
		Get the PAE matrix for the interacting region.
			For this we need the cumulative residue index 
				uptil the required residue position.
		"""
		frag_index_dict = self.get_indices_for_pae( prot1_res, prot2_res )
		start1_idx, end1_idx = frag_index_dict["prot1"]
		start2_idx, end2_idx = frag_index_dict["prot2"]

		required_pae = self.pae[start1_idx:end1_idx, start2_idx:end2_idx]

		return required_pae


	def get_interaction_data( self, res_dict: Dict, prot1_res: List, prot2_res: List ):
		"""
		Get the interaction map, pLDDT, and PAE for the interacting region.
		"""
		chain1, chain2 = self.chains
		coords_dict = self.get_required_coords( res_dict )
		coords1 = coords_dict[chain1]
		coords2 = coords_dict[chain2]
		contact_map = get_contact_map( coords1, coords2, self.dist_threshold )
		
		plddt_dict = self.get_required_plddt( res_dict )
		plddt1 = plddt_dict[chain1]
		plddt2 = plddt_dict[chain2]

		pae = self.get_required_pae( prot1_res, prot2_res )

		return contact_map, plddt1, plddt2, pae


	def apply_confidence_cutoffs( self, plddt1, plddt2, pae ):
		"""
		mask low-confidence interactions.
		"""
		plddt1 = np.where( plddt1 >= self.plddt_cutoff, 1, 0 )
		plddt2 = np.where( plddt2 >= self.plddt_cutoff, 1, 0 )
		plddt_matrix = plddt1 * plddt2.T

		pae = np.where( pae <= self.pae_cutoff, 1, 0 )

		return plddt_matrix, pae


	def validate_input( self, res_dict: Dict, prot1_res: List, prot2_res: List ):
		"""
		Check if the no. of residues selected is the same as the length of the fragment.
		"""
		for i, prot_res in enumerate( [prot1_res, prot2_res] ):
			chain = self.chains[i]
			start, end = prot_res
			frag_len = end - start + 1

			res_selected = len( res_dict[chain] )
			if res_selected > frag_len:
				raise ValueError( "It's raining residues.\n" + 
							f"No. of residues selected - {res_selected} exceeds the fragment length - {frag_len}..." )
			elif res_selected < frag_len:
				raise ValueError( "Missing residues in predicted structures.\n" + 
							f"No. of residues selected - {res_selected} is less than the fragment length - {frag_len}..." )



	def get_confident_interactions( self, prot1_res: List, prot2_res: List ):
		"""
		For the specified regions in the predicted structure, 
			obtain all confident interacting residue pairs.
		Assuming that chain 1 and 2 will correspond to protein 1 and 2 respectively.
		"""
		res_dict = self.get_required_residues( prot1_res, prot2_res )
		
		self.validate_input( res_dict, prot1_res, prot2_res )
		
		chain1, chain2 = self.chains
		interacting_region = {}
		interacting_region[chain1] = prot1_res
		interacting_region[chain2] = prot2_res
		contact_map, plddt1, plddt2, pae = self.get_interaction_data( res_dict, prot1_res, prot2_res )
		plddt_matrix, pae = self.apply_confidence_cutoffs( plddt1, plddt2, pae )
		confident_interactions = contact_map * plddt_matrix * pae

		# m, n = confident_interactions.shape
		# pad = np.zeros( ( self.max_len, self.max_len ) )
		# pad[:m, :n] = confident_interactions
		# confident_interactions = pad

		return confident_interactions


#################################################################
if __name__ == "__main__":
	tic = time.time()
	parser = argparse.ArgumentParser( description = "Script to obtain Disobind predictions for the specified protein pairs.")

	parser.add_argument( 
						"--input", "-f", dest = "f", 
						help = "Provide protein pairs as -- Uni_ID1,start_res1,end_res1,Uni_ID2,start_res2,end_res2 -- in a csv file format.", 
						required = True )
	
	parser.add_argument( 
						"--max_cores", "-c", dest = "c", 
						help = "No. of cores to be used.", 
						type = int, required = False, default = 2 )

	parser.add_argument( 
						"--output_dir", "-o", dest = "o", 
						help = "Name of the output directory.", 
						type = str, required = False, default = "output" )

	parser.add_argument( 
						"--device", "-d", dest = "d", 
						help = "Device to be used (cpu/cuda.", 
						type = str, required = False, default = "cpu" )

	parser.add_argument( 
						"--cmaps", "-cm", dest = "cm", 
						help = "If this flag is provided, will predict contact maps.", 
						action = "store_true", required = False, default = False )

	parser.add_argument( 
						"--coarse", "-cg", dest = "cg", 
						help = "Provid the resolution of coarse graining. If 0, will run for all CG.", 
						type = int, required = False, default = 1 )


	args = parser.parse_args()

	Disobind( args ).forward()
	toc = time.time()
	print( f"Time taken = {toc - tic} seconds" )
	print( "May the Force be with you..." )
#################################################################


