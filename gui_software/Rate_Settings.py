# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 08:04:31 2020

@author: Naylor

for the moment we are going to create two settings menus.  this is for the
sake of quick creation and separating the values
Will not use all of the settings, can swap in or out as necessary
Also since .yaml files can be opened by notepad and are not encrypted
we will not give the user the option to change the defaults in the gui.
we'll start with the rate settings as I am more familiar with
"""

import os

from PyQt5 import uic, QtWidgets

import deuterater.settings as settings
from utils.useful_classes import setting_numerical_info, setting_string_info


#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
paired_settings = [["minimum_allowed_sequence_rate","maximum_allowed_sequence_rate"], 
                   ["lowest_allowed_norm_isotope", "highest_allowed_norm_isotope"]
                   ]

settings_file_ui_location = os.path.join(location, "ui_files", "Settings_Menu.ui")
loaded_ui = uic.loadUiType(settings_file_ui_location)[0]

class Rate_Setting_Menu(QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent = None, current_setting_file = None):
        super(Rate_Setting_Menu, self).__init__(parent)
        settings.load(current_setting_file)
        self.current_setting_file =current_setting_file
        #$this is needed to slim things down a bit
        self.setWindowTitle("Rate Settings Menu")
        self.setupUi(self)
        self.all_settings=[
            setting_string_info(self.recognize_available_cores, "recognize_available_cores",
                                 settings.recognize_available_cores, True),
            setting_numerical_info(self.default_cores, "n_processors",
                                   settings.n_processors, True),
            setting_numerical_info(self.min_allowed_timepoints_enrichment, "min_allowed_timepoints_enrichment",
                                   settings.min_allowed_timepoints_enrichment, True),
            setting_numerical_info(self.starting_enrichment_table_timepoints, "starting_enrichment_table_timepoints",
                                   settings.starting_enrichment_table_timepoints, True),
            setting_string_info(self.rt_unit, "id_file_rt_unit",
                                settings.id_file_rt_unit, False),
            setting_numerical_info(self.time_window, "time_window",
                                    settings.time_window, False),
            setting_numerical_info(self.ppm_error, "ppm_window",
                                    settings.ppm_window, True),
            setting_string_info(self.use_chromatography_division,
                                "use_chromatography_division",
                                settings.use_chromatography_division,
                                False),
            setting_numerical_info(self.mz_prox_filter,
                                   "mz_proximity_tolerance",
                                   settings.mz_proximity_tolerance,
                                   True),
            setting_numerical_info(self.rt_prox_filter,
                                   "rt_proximity_tolerance",
                                   settings.rt_proximity_tolerance,
                                   False),
            setting_string_info(self.label_key, "label_key",
                                settings.label_key, False),
            setting_numerical_info(self.min_AA_length, "min_aa_sequence_length",
                                    settings.min_aa_sequence_length, True),
            setting_numerical_info(self.min_allowed_n_value, "min_allowed_n_values",
                                    settings.min_allowed_n_values, True),
            setting_numerical_info(self.minimum_nonzero_points_rate, "min_non_zero_timepoints_rate",
                                   settings.min_non_zero_timepoints_rate, True),
            setting_numerical_info(self.min_allowed_rate, "minimum_allowed_sequence_rate",
                                   settings.minimum_allowed_sequence_rate, False),
            setting_numerical_info(self.max_allowed_rate, "maximum_allowed_sequence_rate",
                                   settings.maximum_allowed_sequence_rate, False),
            setting_string_info(self.error_graph_option, "error_estimation",
                                settings.error_estimation, False),
            setting_numerical_info(self.minimum_sequences_to_combine_for_protein_rate,
                                   "minimum_sequences_to_combine_for_protein_rate",
                                   settings.minimum_sequences_to_combine_for_protein_rate,
                                   True),
            setting_numerical_info(self.lowest_allowed_norm_isotope, "lowest_allowed_norm_isotope",
                                   settings.lowest_allowed_norm_isotope, False),
            setting_numerical_info(self.highest_allowed_norm_isotope, "highest_allowed_norm_isotope",
                                   settings.highest_allowed_norm_isotope, False),
            setting_numerical_info(self.m0_decreasing_allowed_noise, "m0_decreasing_allowed_noise",
                                   settings.m0_decreasing_allowed_noise, False),
            setting_numerical_info(self.median_absolute_residuals_cutoff_value_single_point,
                                   "median_absolute_residuals_cutoff_single_point",
                                   settings.median_absolute_residuals_cutoff_single_point,
                                   False),
            setting_numerical_info(self.median_absolute_residuals_cutoff_value_two_points,
                                   "median_absolute_residuals_cutoff_two_points",
                                   settings.median_absolute_residuals_cutoff_two_points,
                                   False),
            setting_numerical_info(self.median_absolute_residuals_cutoff_value_general,
                                   "median_absolute_residuals_cutoff_general",
                                   settings.median_absolute_residuals_cutoff_general,
                                   False),
            setting_string_info(self.graph_file_type, 
                                "graph_output_format",
                                settings.graph_output_format,
                                False)

            
            ]
        for setting_object in self.all_settings:
            setting_object.set_object_value()
        self.SaveButton.clicked.connect(self.save_settings)
        self.ExitButton.clicked.connect(self.close)
        
    def save_settings(self):
        #$we need to provide the values that are not altred for the dump
        save_value_dict = Rate_Setting_Menu._get_filters()
        for setting_object in self.all_settings:
            name, value = setting_object.save_value()
            save_value_dict[name] = value
            
        for pair_of_settings in paired_settings:
            if save_value_dict[pair_of_settings[0]] >= save_value_dict[pair_of_settings[1]]:
                error_msg = f"{pair_of_settings[0]} must be less than {pair_of_settings[1]}. Please correct to save."
                QtWidgets.QMessageBox.information(self, "Error", error_msg)
                return False
            
        settings.freeze(self.current_setting_file, save_value_dict)
        return True
        
    def check_for_changes(self):
        for setting in self.all_settings:
            if not setting.compare_value():
                return False
        return True
            
    
    #$should overwrite the close of the exit button and the red x in the corner  
    def closeEvent(self, event):
        if self.check_for_changes():
            event.accept()
        else:
            reply = QtWidgets.QMessageBox.question(self, "Unsaved Changes", 
                ("There are unsaved changes.  Would you like to save before"
                " exiting?"), 
                QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No|QtWidgets.QMessageBox.Cancel)
            if reply == QtWidgets.QMessageBox.No:
                event.accept()        
            elif reply == QtWidgets.QMessageBox.Yes:
                allowed_to_save = self.save_settings()
                if allowed_to_save:
                    event.accept()
                else:
                    event.ignore()  
            #$hit red x or cancel just don't exit
            else: 
                event.ignore()  
      
    #$the point of this is to get the values from the settings .yaml 
    #$that the user is not altering. don't need special classes since 
    #$ we have to officially declare and all we need to do is get the value
    @staticmethod
    def _get_filters():
        unalterable_settings = {
            "debug_level" : settings.debug_level,
            "trim_ids_to_mzml_bounds" : settings.trim_ids_to_mzml_bounds,
            "chunk_size" : settings.chunk_size,
            "chunking_method_threshold" : settings.chunking_method_threshold,
            "peak_ratio_denominator" : settings.peak_ratio_denominator,
            "aa_labeling_sites_path" : settings.aa_labeling_sites_path,
            "peak_lookback" : settings.peak_lookback,
            "peak_lookahead" : settings.peak_lookahead,
            "baseline_lookback" : settings.baseline_lookback,
            "min_envelopes_to_combine" : settings.min_envelopes_to_combine,
            "zscore_cutoff": settings.zscore_cutoff,
            "desired_points_for_optimization_graph": settings.desired_points_for_optimization_graph,
            "max_valid_angle": settings.max_valid_angle,
            "ms_level": settings.ms_level,
            "intensity_filter": settings.intensity_filter,
            "rel_height": settings.rel_height,
            "sampling_rate": settings.sampling_rate,
            "smoothing_width": settings.smoothing_width,
            "smoothing_order": settings.smoothing_order,
            "allowed_peak_variance_min": settings.allowed_peak_variance_min,
            "allowed_neutromer_peak_variance": settings.allowed_neutromer_peak_variance,
            "adduct_weight": settings.adduct_weight,
            "variance_weight": settings.variance_weight,
            "ID_weight": settings.ID_weight,
            "intensity_weight": settings.intensity_weight,
            "how_divided": settings.how_divided,
            "max_enrichment_allowed": settings.max_enrichment_allowed
            }
        return unalterable_settings
        