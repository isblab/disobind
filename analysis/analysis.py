"""
Evaluate Disobind, AF2, AF3 on the OOD test set and Misc dataset.
Also assess the performance of AIUPred, MORFchibi, and DeepDISOBind.
Assess Disobind, AF2, and Disobind+AF2 predictions on:
	Disordered vs ordered regions.
	Disorderpromoting, aromatic, hydrophobuc, and polar residues.
	LIPs
"""
import os, random, json, h5py
from typing import List, Tuple, Dict, Iterator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
from torch import nn

from dataset.utility import ranges

from src.metrics import torch_metrics

MAX_LEN_DICT = {"ood": {19: 100}, "ooddor": {19: 100}, "oodddr": {19: 100}, "oodidr": {19: 100}, "misc": {19: 200}}


class JudgementDay():
	"""
	Pipeline all analysis to be performed.
	"""
	def __init__( self ):
		# Dataset version.
		self.data_version = 19
		self.model_version = 19
		# path for the dataset dir.
		self.base_path = f"../database/"
		# Seed for the PRNGs.
		self.global_seed = 11
		# Cutoff to dfine contact for disobind predictions.
		self.contact_threshold = 0.5
		# ipTM cutoff for confident predictions.
		self.iptm_cutoff = 0.0
		self.pad = True
		self.mode = "ood" # "ood" or "misc"
		# Max prot1/2 lengths.
		self.max_len = MAX_LEN_DICT[self.mode][self.model_version]
		self.prec = 2
		self.device = "cuda"

		if "ood" in self.mode:
			self.base_af_dir = f"./AF_preds_v{self.data_version}/"
			self.diso_base_dir = f"./Predictions_{self.mode}_v_{self.data_version}/"
			# OOD set target contact maps file.
			self.target_cmap = f"{self.base_path}v_{self.data_version}/Target_bcmap_test_v_{self.data_version}.h5"
			# Disobind predictions file.
			self.disobind_preds = f"{self.diso_base_dir}Disobind_Predictions.npy"

		elif self.mode == "misc":
			self.base_af_dir = f"./AF_misc_preds_v{self.data_version}/"
			self.diso_base_dir = f"./Predictions_{self.mode}_v_{self.data_version}/"
			# OOD set target contact maps file.
			self.target_cmap = f"{self.base_path}Misc/misc_test_target.h5"
			# Disobind predictions file.
			self.disobind_preds = f"{self.diso_base_dir}Disobind_Predictions_misc.npy"
			self.summary_file = "../database/Misc/Summary.json"

		else:
			raise ValueError( "Incorrect mode specified (ood/misc supported)..." )

		# AF2 predictions file.
		self.af2m_preds = f"{self.base_af_dir}Predictions_af2m_results_{self.iptm_cutoff}.npy"
		# AF3 predictions file.
		self.af3_preds = f"{self.base_af_dir}Predictions_af3_results_{self.iptm_cutoff}.npy"

		# Predictions from AIUPred and DeepDisoBind
		self.other_methods = "./other_methods/other_methods.npy"

		# OOD set target contact maps file.
		self.fraction_positives = f"{self.base_path}v_{self.data_version}/T5/global-None/fraction_positives.json"
		# File containing info about the merged binary complexes in the dataset.
		self.merged_binary_complexes_dir = f"{self.base_path}v_{self.data_version}/merged_binary_complexes/"

		self.output_dir = f"./Analysis_{self.mode.upper()}_{self.data_version}_{self.iptm_cutoff}/"

		# File to store all results.
		self.full_results_file = f"{self.output_dir}Results_OOD_set_{self.data_version}.csv"
		self.subset_results_file = f"{self.output_dir}Results_OOD_set_subset_{self.data_version}.csv"
		self.other_methods_result_file = f"{self.output_dir}Results_other_methods_{self.data_version}.csv"
		
		# Files for the plots and raw data.
		self.af_conf_pred_counts_file = f"{self.output_dir}Confident_AF_preds_{self.data_version}.txt"
		self.af_confidence_file = f"{self.output_dir}AF_confidence_plot_{self.data_version}"
		self.af_confidence_scores = f"{self.output_dir}AF_confidence_scores_{self.data_version}.csv"
		self.sparsity_file = f"{self.output_dir}Sparsity_F1_plot_{self.data_version}"
		self.case_specific_analysis_file = f"{self.output_dir}Case_sp_analysis_{self.data_version}"
		self.top_preds_diso_af2_file = f"{self.output_dir}Top_preds_Diso_AF2_{self.data_version}"
		self.top_preds_file = f"{self.output_dir}Top_preds_{self.data_version}"
		self.misc_dict_file = f"{self.output_dir}Misc_dict_{self.data_version}.json"


	def seed_worker( self ):
		# Set the seeds for PRNG.
		torch.manual_seed( self.global_seed )
		torch.cuda.manual_seed_all( self.global_seed )
		np.random.seed( self.global_seed )
		random.seed( self.global_seed )


	def forward( self ):
		if not os.path.exists( self.output_dir ):
			os.makedirs( self.output_dir )
		self.seed_worker()

		# Load all predictions.
		self.af2m_preds = np.load( self.af2m_preds, allow_pickle = True ).item()
		self.af3_preds = np.load( self.af3_preds, allow_pickle = True ).item()
		self.disobind_preds = np.load( self.disobind_preds, allow_pickle = True ).item()
		self.target_cmap = h5py.File( self.target_cmap, "r" )
		if self.mode == "ood":
			self.other_methods = np.load( self.other_methods, allow_pickle = True ).item()
		self.fraction_positives = json.load( open( self.fraction_positives, "r" ) )

		self.count_confident_AF_predictions()
		self.eval_performance()


	def get_tasks( self ):
		tasks = []
		for obj in ["interaction", "interface"]:
			for cg in [1, 5, 10]:
				tasks.append( f"{obj}_{cg}" )

		return tasks


	def prepare_target( self, target, task ):
		if self.pad:
			max_len = [self.max_len, self.max_len]
		else:
			max_len = [target.shape]
		
		mask = np.zeros( max_len )
		h, w = target.shape
		mask[:h,:w] = target
		target = mask

		target = torch.from_numpy( target )

		_, cg = task.split( "_" )
		cg = int( cg )
		if cg > 1:
			# MaxPool2d requires input of shape ( N, C, L1, L2 ).
			# 	target: ( N, L1, L2 ) --> ( N, 1, L1, L2 )
			eff_len = [max_len[0]//cg, max_len[1]//cg]
			
			with torch.no_grad():
				target = target.unsqueeze( 0 )
				m = nn.MaxPool2d( kernel_size = cg, stride = cg )
				# print( cg, "\t", target.size() )
				target = m( target )
				# Remove the extra C dimension.
				target = target.squeeze( 0 )
		else:
			eff_len = max_len

		if "interface" in task:
			# Find indexes of conatct elements.
			idx = torch.where( target == 1 )
			# Get interacting residues for prot1.
			p1_target = torch.zeros( ( target.shape[0], 1 ) )
			p1_target[idx[0], :] = 1

			# Get interacting residues for prot2.
			p2_target = torch.zeros( ( target.shape[1], 1 ) )
			p2_target[idx[1], :] = 1
			
			# target: ( N, L1, L2 ) --> ( N, L1+L2 )
			target = torch.cat( ( p1_target, p2_target ), axis = 0 )

		return target


	def get_random_baseline( self, target, task ):
		# Fraction of contacts in train set.
		obj, cg = task.split( "_" )
		p = self.fraction_positives[obj][cg]
		n = 1 - p

		random_preds = np.random.choice( [0, 1], ( target.shape ), [n, p] )

		return random_preds


	def create_ood_set_tensors( self ) -> Iterator[Tuple[str, Dict[str, np.array], List]]:
		"""
		Generator function that loads all predictions, masks, target cmap, disorder matrices 
			from Disobind, AF2, AF3 OOD prediction files, create a batch for all OOD entries 
			and yields an OOD dict for each task.
		This function returns the following:
			task --> mme of the task (e.g. interaction_1).
			ood_dict: dict containing all required tensors for the dataset
						(ood/misc) for a task.
			entry_ids --> a list of all the dataset entry_id.
		"""
		# For all tasks.
		for task in self.get_tasks():
			counts = [0, 0, 0]
			ood_keys = ["AF2_pLDDT_PAE", "AF3_pLDDT_PAE",
					"Disobind", "Random_baseline", "target_mask", 
					"IDR-IDR", "order",
					"lips", "disorder_promoting_aa", "aromatic_aa",
					"hydrophobic_aa", "polar_aa", "targets"]
			if self.mode  == "ood":
				ood_keys.extend( ["Aiupred", "Deepdisobind", "Morfchibi"] )

			# All fields to be included for the results.
			ood_dict = {key:[] for key in ood_keys}
			
			# For all entries.
			entry_ids = []
			for idx, entry_id in enumerate( self.disobind_preds.keys() ):
				target = np.array( self.target_cmap[entry_id] )
				target = self.prepare_target( target, task )

				# Ignoring this entry, as AF2-multimer crashed for this.
				if entry_id == "P0DTD1:1743:1808--P0DTD1:1565:1641_1":
					continue

				u1, u2 = entry_id.split( "--" )
				u1 = u1.split( ":" )[0]
				u2, c = u2.split( "_" )
				u2 = u2.split( ":" )[0]

				entry_ids.append( entry_id )

					# These Uniprot pairs are sequence redundant with PDB70 at 20% seq identity.
					# 	Ignoring these from evaluation.
				if f"{u1}--{u2}_{c}" in ["P0DTC9--P0DTD1_2", "Q96PU5--Q96PU5_0", "P0AG11--P0AG11_4", 
										"Q9IK92--Q9IK91_0", "Q16236--O15525_0", "P12023--P12023_0",
										"O85041--O85043_0", "P25024--P10145_0"]:
					continue

				# For all fields in the prediction dict.
				for ood_key in ood_dict.keys():
					if "AF2" in ood_key:
						ood_dict[ood_key].append( self.af2m_preds[entry_id][task][ood_key] )
					
					elif "AF3" in ood_key:
						ood_dict[ood_key].append( self.af3_preds[entry_id][task][ood_key] )

					elif "Disobind" in ood_key:
						ood_dict[ood_key].append( self.disobind_preds[entry_id][task][ood_key] )

					elif ood_key in ["Aiupred", "Morfchibi", "Deepdisobind"]:
						if "interface" in task and self.mode == "ood":
							method = self.other_methods[ood_key.lower()]
							ood_dict[ood_key].append( method[f"{u1}--{u2}_{c}"][task] )
						# Empty arrays for all other tasks.
						else:
							dummy = self.disobind_preds[entry_id][task]["Disobind"]
							ood_dict[ood_key].append( np.zeros( dummy.shape ) )

					elif ood_key in ["target_mask", "IDR-IDR", "order"]:
						ood_dict[ood_key].append( self.disobind_preds[entry_id][task][ood_key] )

					elif ood_key == "targets":
						# Apply the padding amsk.
						ood_dict[ood_key].append( target*self.disobind_preds[entry_id][task]["target_mask"] )

					elif "Random" in ood_key:
						random_preds = self.get_random_baseline( target, task )
						random_preds = random_preds*self.disobind_preds[entry_id][task]["target_mask"]

						ood_dict[ood_key].append( random_preds )

					elif "lips" in ood_key:
						if self.mode == "ood":
							lips1 = self.disobind_preds[entry_id][task]["prot1_lips_mask"]
							lips2 = self.disobind_preds[entry_id][task]["prot2_lips_mask"]

							if "interaction" in task:
								# Consider LIP-any interactions.
								lips2 = np.ones( lips2.shape[0] )
								lips_mat = lips1*lips2.T
							elif task == "interface_1":
								lips_mat = np.concatenate( ( lips1, lips2 ), axis = 0 )
							else:
								continue
								# raise ValueError( "Only interface task supported..." )
							ood_dict[ood_key].append( lips_mat )
						# For non ood dataset, just create empty tensors.
						else:
							dummy = self.disobind_preds[entry_id][task]["Disobind"]
							ood_dict[ood_key].append( np.zeros( dummy.shape ) )

					elif ood_key in ["disorder_promoting_aa", "aromatic_aa", "hydrophobic_aa", "polar_aa"]:
						if self.mode == "ood":
							aa1_mask = self.disobind_preds[entry_id][task]["prot1_aa_mask"][ood_key]
							aa2_mask = self.disobind_preds[entry_id][task]["prot2_aa_mask"][ood_key]
							if "interaction" in task:
								aa_mask = aa1_mask*aa2_mask.T
							elif task == "interface_1":
								aa_mask = np.concatenate( ( aa1_mask, aa2_mask ), axis = 0 )
							else:
								continue
							ood_dict[ood_key].append( aa_mask )
						# For non ood dataset, just create empty tensors.
						else:
							dummy = self.disobind_preds[entry_id][task]["Disobind"]
							ood_dict[ood_key].append( np.zeros( dummy.shape ) )

			for key in ood_dict.keys():
				if len( ood_dict[key] ) == 0:
					continue
				else:
					ood_dict[key] = np.stack( ood_dict[key] )

			yield task, ood_dict, entry_ids


	def get_base_model_preds( self, ood_dict: Dict[str, np.array] ):
		"""
		Get predictions and the respective target for Disobind, AF2, AF3, and Random baseline.
		"""
		base_models_dict = {}
		target = ood_dict["targets"]
		for model_key in ["AF2_pLDDT_PAE", "AF3_pLDDT_PAE", "Disobind", "Random_baseline"]:
			preds = ood_dict[model_key]

			base_models_dict[model_key] = {
					"pred": preds, "target": target
				}
		return base_models_dict


	def combine_diso_af_preds( self, ood_dict: Dict[str, np.array],
								af_model: str ) -> Dict[str, np.array]:
		"""
		Combine Disobind and AF2/3 predictions by taking the max over either.
		Return a dict containing the combine dmodel prediction and the respective target.
		"""
		combined_model_name = f"Disobind_{af_model}"
		af_key = f"{af_model}_pLDDT_PAE"
		diso_key = "Disobind"

		diso = ood_dict[diso_key]
		af = ood_dict[af_key]
		b, m, n = diso.shape
		af_diso = np.stack( [af.reshape( b, m*n ), diso.reshape( b, m*n )], axis = 1 )
		af_diso = np.max( af_diso, axis = 1 ).reshape( b, m, n )

		target = ood_dict["targets"]
		return {combined_model_name: {"pred":af_diso, "target": target}}


	def get_preds_for_disorder_order_residues( self, ood_dict: Dict[str, np.array],
												preds_dict: Dict[str, np.array]
												) -> Dict[str, Dict[str, np.array]]:
		"""
		Obtain predictions focusing specifically on disordered or ordered residues.
		Get predictions from Disobind, AF2, AF3, Diso+AF2, Diso+AF3.
		"""
		interactions_dict = {}
		target = ood_dict["targets"]
		for interaction_name in ["IDR-IDR", "order"]:
			for model_key in ["Disobind", "AF2_pLDDT_PAE", "AF3_pLDDT_PAE", "Disobind_AF2"]:
				if model_key in ood_dict:
					preds = ood_dict[model_key]
				elif model_key in preds_dict:
					preds = preds_dict[model_key]["pred"]
				else:
					raise ValueError( f"{model_key} not present in ood_dict and preds_dict..." )
			
				interaction_mat = ood_dict[interaction_name]

				interaction_pred = preds*interaction_mat
				interaction_target = target*interaction_mat

				inetraction_key = f"{model_key}_{interaction_name}"
				interactions_dict[inetraction_key] = {
						"pred": interaction_pred,
						"target": interaction_target
					}
		return interactions_dict


	def get_preds_for_interaction_types( self, ood_dict: Dict[str, np.array],
											preds_dict: Dict[str, np.array]
											) -> Dict[str, Dict[str, np.array]]:
		target = ood_dict["targets"]
		interaction_types_dict = {}
		for mask_name in ["disorder_promoting_aa", "aromatic_aa", "hydrophobic_aa", "polar_aa", "lips"]:
			for model_name in ["Disobind", "AF2_pLDDT_PAE", "AF3_pLDDT_PAE", "Disobind_AF2"]:

				if model_name in ood_dict:
					preds = ood_dict[model_name]
				elif model_name in preds_dict:
					preds = preds_dict[model_name]["pred"]
				else:
					print( mask_name, "  ", model_name )
					raise ValueError( f"{mask_name}/{model_name} not " +
										"present in ood_dict and preds_dict...")

				
				mask_mat = ood_dict[mask_name]
				interaction_pred = preds*mask_mat
				interaction_target = target*mask_mat

				interaction_key = f"{model_name}_{mask_name}"
				interaction_types_dict[interaction_key] = {
						"pred": interaction_pred,
						"target": interaction_target
					}
		return interaction_types_dict



	def get_other_method_preds( self, ood_dict: Dict[str, np.array]
								) -> Dict[str, Dict[str, np.array]]:

		"""
		Get the predictions and the respective targets for
			Aiupred, Deepdisobind, Morfchibi.
		Considering only protein prediction.
		"""
		other_methods_dict = {}
		target = ood_dict["targets"]
		p1_target = target[:self.max_len]
		for model_name in ["Aiupred", "Morfchibi", "Deepdisobind"]:
			preds = ood_dict[model_name]

			other_methods_dict[model_name] = {
						"pred": preds,
						"target": p1_target
					}
		return other_methods_dict



	def eval_performance( self ):
		"""
		Calculate metrics for OOD set across all tasks and models 
			and save the results on disk.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		None
		"""
		results_dict = {key:[] for key in ["Objective", "CG", "Model", "Recall", "Precision", "F1-score"]}

		for task, ood_dict, entry_ids in self.create_ood_set_tensors():
			print( f"Task {task}..." )
			preds_dict = {}

			base_models_dict = self.get_base_model_preds( ood_dict = ood_dict )
			preds_dict.update( base_models_dict )

			af2_diso = self.combine_diso_af_preds( ood_dict = ood_dict, af_model = "AF2" )
			preds_dict.update( af2_diso )

			if task in ["interaction_1", "interface_1"]:
				interactions_dict = self.get_preds_for_disorder_order_residues( ood_dict = ood_dict,
																				preds_dict = preds_dict )
				preds_dict.update( interactions_dict )

				interaction_types_dict = self.get_preds_for_interaction_types( ood_dict = ood_dict,
																				preds_dict = preds_dict )
				preds_dict.update( interaction_types_dict )


			if "interface" in task:
				if self.mode == "ood":
					other_methods_dict = self.get_other_method_preds( ood_dict = ood_dict )
					preds_dict.update( other_methods_dict )

				if self.mode == "misc":
					# Case specific analysis for interface 1 prediction.
					self.case_specific_analysis( entry_ids, preds_dict )

			# now calculate the metrics for all predictions.
			for model_key in preds_dict.keys():
				obj, cg = task.split( "_" )
				# print( f"{model_key}" )
				preds = torch.from_numpy( preds_dict[model_key]["pred"] ).float()
				target = torch.from_numpy( preds_dict[model_key]["target"] ).float()
				metrics = self.calculate_metrics( pred = preds.to( self.device ),
													target = target.to( self.device ),
													multidim_avg = "global",
													contact_threshold = self.contact_threshold )
				results_dict["Objective"].append( obj.title() )
				results_dict["CG"].append( cg )
				results_dict["Model"].append( model_key )
				results_dict["Recall"].append( metrics[0] )
				results_dict["Precision"].append( metrics[1] )
				results_dict["F1-score"].append( metrics[2] )


			# Just adding a empty row for separating different tasks.
			for key in results_dict.keys():
				results_dict[key].append( "" )

		self.create_sparsity_f1_plots( results_dict )

		# Dump all calculated metrics on disk.
		df = pd.DataFrame( results_dict )
		df.to_csv( self.full_results_file )


	# def get_samplewise_metrics( self, entries, diso_preds, af2_preds, diso_af2_preds, target ):
	def case_specific_analysis( self, entries, preds_dict: torch.tensor ):
		"""
		Save the samplewise F1-score for interface_1 on the Misc dataset.
		"""
		with open( self.summary_file, "r" ) as f:
			summary = json.load( f )
		pdb_ids = []
		for entry_id in entries:
			pdb_ids.append( summary[entry_id]["pdb_id"] )

		entries = np.array( entries )
		metrics_dict = {}
		for model_key in ["Disobind", "AF2_pLDDT_PAE", "Disobind_AF2"]:
			metrics_dict[model_key] = {}

			preds = torch.from_numpy( preds_dict[model_key]["pred"] ).float()
			target = torch.from_numpy( preds_dict[model_key]["target"] ).float()

			metrics = torch_metrics( preds = preds.to( self.device ), 
									target = target.to( self.device ), 
									threshold = self.contact_threshold,
									multidim_avg = "samplewise_none",
									device = self.device )
			met = metrics[2].detach().cpu().numpy()
			metrics_dict[model_key]["F1"] = np.round( met, self.prec )

		# Save samplewsie metrics on OOD set.
		df = pd.DataFrame()
		df["Entry ID"] = entries
		df["PDB ID"] = pdb_ids
		for model_name in metrics_dict:
			for met_name in metrics_dict[model_name]:
				df[f"{model_name}-{met_name}"] = metrics_dict[model_name][met_name]

		df.to_csv( f"{self.case_specific_analysis_file}.csv", index = False )

		misc_dict = {}
		af2_diso = preds_dict["Disobind_AF2"]["pred"].copy()
		target = preds_dict["Disobind_AF2"]["target"].copy()
		for i, entry_id in enumerate( entries ):
			misc_dict[entry_id] = {}
			misc_dict[entry_id]["AF2+Diso"] = list( af2_diso[i].reshape( -1 ) )
			misc_dict[entry_id]["Target"] = list( target[i].astype( float ).reshape( -1 ) )
		with open( self.misc_dict_file, "w" ) as w:
			json.dump( misc_dict, w )


	def eval_other_methods( self, diso_preds: torch.tensor, aiupred_preds: torch.tensor,
							deepdisobind_preds: torch.tensor, morfchibi_preds: torch.tensor,
							target: torch.tensor ):
		"""
		Compare Disobind performance with AIUPred, DeepDisoBind, and MORfchibi.
		The later 3 are partner-independent interface predictors for IDRs and
			so we evaluate performance only for the prot1 in all OOD entries
			which is an IDR.
		"""
		# Need to slice out the protein1 interface only from diso_preds and target.
		diso_p1 = diso_preds[:,:self.max_len].to( self.device )
		aiupred_p1 = aiupred_preds[:,:self.max_len].to( self.device )
		deepdisobind_p1 = deepdisobind_preds[:,:self.max_len].to( self.device )
		morfchibi_p1 = morfchibi_preds[:,:self.max_len].to( self.device )
		target_p1 = target[:,:self.max_len].to( self.device )

		diso_metrics = self.calculate_metrics( pred = diso_p1,
													target = target_p1,
													multidim_avg = "global",
													contact_threshold = self.contact_threshold )
		aiupred_metrics = self.calculate_metrics( pred = aiupred_p1,
													target = target_p1,
													multidim_avg = "global",
													contact_threshold = self.contact_threshold )
		deepdisobind_metrics = self.calculate_metrics( pred = deepdisobind_p1,
														target = target_p1,
														multidim_avg = "global",
														contact_threshold = self.contact_threshold )
		morfchibi_metrics = self.calculate_metrics( pred = morfchibi_p1,
													target = target_p1,
													multidim_avg = "global",
													contact_threshold = 0.775 )


		df = pd.DataFrame()
		df["Model"] = ["Disobind_uncal", "Aiupred", "Deepdisobind", "Morfchibi"]
		df["Recall"] = [diso_metrics[0], aiupred_metrics[0], deepdisobind_metrics[0], morfchibi_metrics[0]]
		df["Precision"] = [diso_metrics[1], aiupred_metrics[1], deepdisobind_metrics[1], morfchibi_metrics[1]]
		df["F1-score"] = [diso_metrics[2], aiupred_metrics[2], deepdisobind_metrics[2], morfchibi_metrics[2]]
		df.to_csv( self.other_methods_result_file, index = False )


	def count_confident_AF_predictions( self ):
		"""
		Count the no. of confident AF2 and AF3 predictions based on:
			iptm score for AF2 and AF3.
		Plot AF2-ipTM vs AF3-ipTM.
		Save the raw data on disk.

		Input:
		----------
		Does not take any arguments.

		Returns:
		----------
		None
		"""
		counts1, counts2, counts3, counts4, counts5, counts6 = 0, 0, 0, 0, 0, 0
		af2_score, af3_score = [], []

		selected_entries = []
		for key in self.af2m_preds.keys():
			if key == "P0DTD1:1743:1808--P0DTD1:1565:1641_1":
				continue

			u1, u2 = key.split( "--" )
			u1 = u1.split( ":" )[0]
			u2, c = u2.split( "_" )
			u2 = u2.split( ":" )[0]

			selected_entries.append( key )
			af2_score.append( self.af2m_preds[key]["scores"][0] )
			af3_score.append( self.af3_preds[key]["scores"][0] )

			if af2_score[-1] <= self.iptm_cutoff:
				counts1 += 1
			if af3_score[-1] <= self.iptm_cutoff:
				counts2 += 1
			if af2_score[-1] >= 0.8:
				counts3 += 1
			if af3_score[-1] >= 0.8:
				counts4 += 1
			if af2_score[-1] >= self.iptm_cutoff and af3_score[-1] <= self.iptm_cutoff:
				counts5 += 1
			if af2_score[-1] >= af3_score[-1]:
				counts6 += 1

		plt.plot( [0, 1], [0, 1], color = "gray" )
		plt.scatter( af2_score, af3_score )
		plt.axvline( self.iptm_cutoff, color = "red" )
		plt.axhline( self.iptm_cutoff, color = "red" )
		plt.xlabel( "AF2 iptm" )
		plt.ylabel( "AF3 iptm" )
		plt.savefig( f"{self.af_confidence_file}.png" , dpi = 300 )
		plt.close()

		df = pd.DataFrame()
		df["Entry_id"] = selected_entries
		df["AF2_iptm"] = af2_score
		df["AF3_iptm"] = af3_score
		df.to_csv( f"{self.af_confidence_file}.csv" )
		
		with open( self.af_conf_pred_counts_file, "w" ) as w:
			w.writelines( f"AF2 preds with ipTM <= {self.iptm_cutoff} = {counts1}\n" )
			w.writelines( f"AF3 preds with ipTM <= {self.iptm_cutoff} = {counts2}\n" )
			w.writelines( f"AF2 preds with ipTM >= 0.8 = {counts3}\n" )
			w.writelines( f"AF3 preds with ipTM >= 0.8 = {counts4}\n" )
			w.writelines( f"AF2  ipTM >= {self.iptm_cutoff} and AF3  ipTM <= {self.iptm_cutoff} = {counts5}\n" )
			w.writelines( f"AF2  ipTM >= AF3  ipTM = {counts6}\n" )


	def create_sparsity_f1_plots( self, results_dict ):
		"""
		Create the plot for sparsity in dataset vs the model performance.
		Save the raw data on disk.

		Input:
		----------
		results_dict --> dict containing metrics for all models on all the tasks.

		Returns:
		----------
		None
		"""
		print( "\nCreating aparsity vs F1 score plot... plot..." )

		sparsity = []
		print( self.fraction_positives )
		for entry_id in self.fraction_positives.keys():
			for key2 in self.fraction_positives[entry_id].keys():
				sparsity.append( round( 1 - self.fraction_positives[entry_id][key2], 4 )*100 )

		f1 = []
		for i in range( len( results_dict["Objective"] ) ):
			if results_dict["Model"][i] == "Disobind":
				f1.append( results_dict["F1-score"][i] )

		plt.plot( sparsity, f1, marker = "o", color = "blue" )
		plt.xlabel( "Sparsity" )
		plt.ylabel( "F1-score" )
		plt.savefig( f"{self.sparsity_file}.png", dpi = 300 )
		plt.close()

		df = pd.DataFrame()
		df["Sparsity"] = sparsity
		df["F1"] = f1
		df.to_csv( f"{self.sparsity_file}.csv" )



	def calculate_metrics( self, pred, target, multidim_avg, contact_threshold ):
		"""
		Calculate the following metrics:
			Recall, Precision, F1score, AvgPrecision, MCC, AUROC, Accuracy.

		Input:
		----------
		pred --> (torch.tensor) model output.
		target --> (torch.tensor) binary output labels.
		multidim_avg --> averaging scheme - global/samplewise/samplewise-none.

		Returns:
		----------
		metric_array --> np.array containing the calculated metric values in order:
			Recall, Precision, F1score, AvgPrecision, MCC, AUROC, Accuracy.	
		"""
		metrics = torch_metrics( pred, target, contact_threshold, multidim_avg, self.device )
		metric_array = np.array( [
					metrics[0].item(),
					metrics[1].item(),
					metrics[2].item(),
					metrics[3].item(),
					metrics[4].item(),
					metrics[5].item(),
					metrics[6].item()
					] ) #.reshape( 1, 7 )

		return np.round( metric_array, 2 )


if __name__ == "__main__":
	JudgementDay().forward()
	print( "May the Force be with you..." )
