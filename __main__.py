# -*- coding: utf-8 -*-
"""
Copyright (c) 2021 Bradley Naylor, Christian Andersen, Michael Porter, Kyle Cutler, Chad Quilling, J.C. Price, and Brigham Young University
All rights reserved.
Redistribution and use in source and binary forms,
with or without modification, are permitted provided
that the following conditions are met:
    * Redistributions of source code must retain the
      above copyright notice, this list of conditions
      and the following disclaimer.
    * Redistributions in binary form must reproduce
      the above copyright notice, this list of conditions
      and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors
      may be used to endorse or promote products derived
      from this software without specific prior written
      permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


#$we will of course need to expand things later, but we'll sort that out later
import os
import sys
import multiprocessing as mp
import csv
import pandas as pd
from tqdm import tqdm

from PyQt5 import uic, QtWidgets, QtCore, QtGui
from shutil import copyfile, rmtree

from deuteconvert.peaks85 import Peaks85
from deuteconvert.peaksXplus import PeaksXplus
from deuteconvert.peaksXpro import PeaksXpro
from deuterater.extractor import Extractor
from gui_software.Time_Table import TimeWindow
from gui_software.Enrichment_Table import EnrichmentWindow
from utils.enrichment_class_fitter import PerformEnrichmentClass
from deuterater.combine_extracted_files import CombineExtractedFiles
from deuterater.initial_intensity_calculator import theoretical_enrichment_calculator
from deuterater.rate_calculator import RateCalculator
from deuterater.protein_rate_combiner import Peptides_to_Proteins
from utils.chromatography_division import ChromatographyDivider
from utils.useful_classes import deuterater_step, deuteconvert_peaks_required_headers
import deuterater.settings as settings
from gui_software.Rate_Settings import Rate_Setting_Menu

#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.abspath(__file__))

#$set locations for the various settings files
rate_settings_file = os.path.join(location, "resources","temp_settings.yaml")
default_rate_settings = os.path.join(location, "resources","settings.yaml")
guide_settings_file = os.path.join(location, "resources", 
                                   "temp_guide_settings.yaml")
default_guide_settings = os.path.join(location, "resources", 
                                      "guide_settings.yaml")

#$each analysis step gets a deuterater_step object.  first argument is the output name
#$so we can check for overwriting (except extract which has dynamic output names)
#$second argument is columns required of the input so we don't have to check
#$when we're actually doing the calculations.
Extract_object = deuterater_step("", ['Precursor Retention Time (sec)',
                    'Precursor m/z','Identification Charge', 'Sequence',
                    'Protein ID', "cf"])
Time_Enrich_object = deuterater_step("time_enrichment_data.tsv", [
    "Precursor Retention Time (sec)", "Protein ID", "Protein Name", "Precursor m/z",
    "Identification Charge", "Homologous Proteins", "n_isos", 
    "Sequence", "cf", "abundances", "mzs"])
Combine_object = deuterater_step("combined_extracted_files_output.tsv", [
    "Filename", "Time","Subject ID", "Spacing Column 1",
    "Subject ID Enrichment", "Time Enrichment", "Enrichment"])
delta_by_enrichment = deuterater_step("delta_by_enrichment.tsv", [
    "Precursor Retention Time (sec)", "Protein ID", "Protein Name", "Precursor m/z",
    "Identification Charge", "Homologous Proteins", "n_isos", "time","literature_n",
    "Sequence", "cf", "abundances", "mzs", "sample_id", "Time Enrichment", "Enrichment Values"])
sequence_rate_calculation = deuterater_step("rate_by_sequence.csv", ["Protein ID", "Protein Name",
    "Sequence","n_isos", "time", "sample_id", "Time Enrichment", "Enrichment Values",
    "abundances", "Theoretical Unlabeled Normalized Abundances", "n_isos", "literature_n"
    ])

protiein_rate_combination = deuterater_step("Final_Protein_Rates.csv", [
    "Subject ID", "Protein ID", "Protein Name", "Sequence", "Abundance rate", 
    "Unique Timepoints", "Number of Measurements", "Approximate Variance",
    "mean of all absolute residuals", "num times",
    "n_isos", "num measurements", "time values",
    "dropped points", "M0 constantly decreasing",
    "Error column"
])
#$create a dictionary so we can call the deuterater_step objects
step_object_dict = {
    "Extract":Extract_object,
    "Provide Time and Enrichment":Time_Enrich_object,
    "Combine Extracted Files": Combine_object,
    "Calculate Delta by Enrichment":delta_by_enrichment, 
    "Rate Calculation": sequence_rate_calculation,
    "Combine Sequence Rates": protiein_rate_combination
    }
#$ converter only does something other than make a file with a header for 
#$id files made with Peaks.  each version is sligthly different so needs a differnt
#$analysis
convert_options = {
    "Peaks 8.5": Peaks85,
    "Peaks X+": PeaksXplus,
    "Peaks XPro": PeaksXpro,
    "Template": ""
    }

#$need 3 different output files from peaks, make sure the needed headers are there 
convert_needed_headers = {
    "Peaks 8.5": deuteconvert_peaks_required_headers(
        ['Accession', '#Peptides', '#Unique', 'Description'],
        ['Peptide', 'Quality', 'Avg. ppm', 'Start', 'End', 'PTM'],
        ['Peptide', 'RT range', 'RT mean', 'm/z', 'z', 'Accession', 'PTM']
        ),
    "Peaks X+": deuteconvert_peaks_required_headers(
        ['Accession', '#Peptides', '#Unique', 'Description'],
        ['Peptide', 'ppm', 'Start', 'End', 'PTM'],
        ['DB Peptide', 'Denovo Peptide', 'RT Begin', 'RT End',
            'RT', 'm/z', 'z', 'Accession']
        ),
    "Peaks XPro": deuteconvert_peaks_required_headers(
        ['Accession', '#Peptides', '#Unique', 'Description'],
        ['Peptide', 'ppm', 'Start', 'End', 'PTM'],
        ['DB Peptide', 'Denovo Peptide', 'RT Begin', 'RT End',
            'RT', 'm/z', 'z', 'Accession']
        )
    }
#$default converter option
default_converter = "Template"
#$ get the header we want if using Peaks
converter_header = PeaksXplus.correct_header_order


#$prepare the gui
main_file_ui_location = os.path.join(location, "ui_files", "Main_Menu.ui")
loaded_ui = uic.loadUiType(main_file_ui_location)[0]
class MainGuiObject(QtWidgets.QMainWindow, loaded_ui):
    def __init__(self, parent = None):
        QtWidgets.QMainWindow.__init__(self,parent)
        #$allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint | 
                            QtCore.Qt.WindowMaximizeButtonHint)
        
        self.setupUi(self)
        make_temp_file(default_rate_settings, rate_settings_file)
        make_temp_file(default_guide_settings, guide_settings_file)
        self.file_loc = location
        
        #$set up the guide file options
        self.guide_file_options.addItems(convert_options.keys())
        #$set the default value for the converter
        index = self.guide_file_options.findText(default_converter)
        self.guide_file_options.setCurrentIndex(index)
        
        self.GuideFileButton.clicked.connect(self.create_guide_file)
        self.RateCalculationButton.clicked.connect(self._calc_rates)
        self.actionSettings.triggered.connect(self.change_settings)
        
        #$make the logo show up
        #$use of command from http://stackoverflow.com/questions/8687723/pyqthow-do-i-display-a-image-properly 
        #$first answer accesed 5/27/2016
        myPixmap = QtGui.QPixmap(os.path.join(location, "resources", "Logo.PNG")) 
        self.Logo.setPixmap(myPixmap)
        self.Logo.setScaledContents(True)
        self.setWindowTitle("DeuteRater")
    
    #$protmpt the user for necessary Peaks outputs
    def Peaks_File_Collection(self, header_checker_object):
         #$ get the files we need
        #TODO switch to reading from a folder instead of individual files?$
        QtWidgets.QMessageBox.information(self, "Info", ("Please select the "
            "files you would like to turn into a guide file. The order is \""
            "proteins.csv\", \"protein-peptides.csv\", \"feature.csv\""))
            
        protein_file, file_type = QtWidgets.QFileDialog.getOpenFileName(self, 
                "Choose protein file to Load", self.file_loc, "CSV (*.csv)")
        if protein_file == "": 
            return ""
        else:
            has_needed_columns = header_checker_object.protein_file_check(protein_file)
            if not has_needed_columns:
                QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                    "missing needed columns. Please correct and try again".format(
                    protein_file)))
                return ""
            self.file_loc = os.path.dirname(protein_file)
            
        protein_peptide_file, file_type = QtWidgets.QFileDialog.getOpenFileName(
            self, "Choose protein_peptide file to Load", self.file_loc, 
            "CSV (*.csv)")
        if protein_peptide_file == "":  
            return ""
        else:
            has_needed_columns = header_checker_object.protein_peptide_check(protein_peptide_file)
            if not has_needed_columns:
                QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                    "missing needed columns. Please correct and try again".format(
                    protein_peptide_file)))
                return ""
            
            self.file_loc = os.path.dirname(protein_peptide_file)
            
        feature_file, file_type = QtWidgets.QFileDialog.getOpenFileName(self, 
                "Choose features file to Load", self.file_loc, "CSV (*.csv)")
        if feature_file == "": return ""
        else:
            has_needed_columns = header_checker_object.features_check(feature_file)
            if not has_needed_columns:
                QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                    "missing needed columns. Please correct and try again".format(
                    feature_file)))
                return ""
            
            self.file_loc = os.path.dirname(feature_file)
        return [protein_file, protein_peptide_file, feature_file]
    
    #$this is to govern the different guide file functions
    def create_guide_file(self):
        guide_file_type = str(self.guide_file_options.currentText())
        #$collect any guide files needed
        #$template doesn't need one since it just needs one output
        if guide_file_type in ["Peaks 8.5", "Peaks X+", "Peaks XPro"]:
             input_files = self.Peaks_File_Collection(convert_needed_headers[guide_file_type])
        
        #$guide_file_type has to be first or input_files may not be defined
        if guide_file_type != "Template" and input_files =="":
            return
        if guide_file_type != "Template":
             #$do the actual calculations
            converter = convert_options[guide_file_type](input_files,
                guide_settings_file)
            converter.convert()                               
        
        #$get output file
        QtWidgets.QMessageBox.information(self, "Info", ("Your guide file was "
            "created. Please select the output file location"))
        while (True):
            save_file, filetype = QtWidgets.QFileDialog.getSaveFileName(self, 
                "Provide Save File", self.file_loc, "CSV (*.csv)")
            if save_file == "": return
            try:
                if guide_file_type != "Template":
                    converter.write(save_file)
                else:
                    df = pd.DataFrame(columns =converter_header )
                    df.to_csv(save_file, sep=',', index=False)
                break
            except IOError:
                QtWidgets.QMessageBox.information(self, "Error", 
                    ("File {} is open in another program. Please close it and "
                    "try again or select a different file".format(save_file)))
        self.file_loc = os.path.dirname(save_file)
        QtWidgets.QMessageBox.information(self, "Success", 
                "Guide file successfully saved")
    
    #$adjusts the text fo the "Rate Calculation" button
    def _calc_rates(self):
        try:
            
            self.RateCalculationButton.setText("Currently Processing... Please Wait.")
            self.run_rate_workflow()
        finally:
            self.RateCalculationButton.setText("Rate Calculation")
     
    #$calls the various functions to calculate turnover rates
    def run_rate_workflow(self):
        #$load the settings so we can use them
        settings.load(rate_settings_file)
        
        #$first we need to check which steps are checked 
        worklist = self.check_table_checklist()
        #$ only proceed if we have steps to perform
        if worklist == []:
            QtWidgets.QMessageBox.information(self, "Error", ("No options were "
                "checked. Please check steps to perform and try again"))
            return
        #$worklist is a string in the case of an error
        elif type(worklist) == str:
            QtWidgets.QMessageBox.information(self, "Error", worklist)
            return
        #$second we need to an output folder and check it for the ability to write
        #$ and if we will be overwriting anything.
        QtWidgets.QMessageBox.information(self, "Info", ("Please select folder "
                "for output"))
        output_folder = QtWidgets.QFileDialog.getExistingDirectory(
            self, 
            "Select an Output Folder", 
            self.file_loc, 
            QtWidgets.QFileDialog.ShowDirsOnly)
        if output_folder == "": return
        #$change location we start asking for things at
        #$don't change after this point since all output is going in here
        self.file_loc = output_folder
        
        #don't care if overwrite rate_settings.yaml but should check if want to use the settings already in the folder
        if os.path.exists(os.path.join(output_folder, "rate_settings.yaml")):
            comp_result = settings.compare(rate_settings_file, os.path.join(output_folder, "rate_settings.yaml"))
            if comp_result != "MATCH":
                if comp_result == "Mismatched Keys":
                    qBox = QtWidgets.QMessageBox(self)
                    qBox.setWindowTitle("Question")
                    question = "A settings file already exists in this output folder. Would you like to use those settings,or overwrite them?"
                    qBox.setText(question)
                    qBox.setIcon(QtWidgets.QMessageBox.Question)
                    qBox.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel)

                    yButton = qBox.button(QtWidgets.QMessageBox.Yes)
                    yButton.setText("Use Settings")
                    nButton = qBox.button(QtWidgets.QMessageBox.No)
                    nButton.setText("Overwrite")
                    response = qBox.exec_()
                    if response == QtWidgets.QMessageBox.Yes:
                        settings.load(os.path.join(output_folder, "rate_settings.yaml"))
                        settings.freeze(rate_settings_file)
                    elif response == QtWidgets.QMessageBox.No:
                        if self.check_file_removal([os.path.join(output_folder, "rate_settings.yaml")], ask_permission = False):
                            settings.freeze(os.path.join(output_folder, "rate_settings.yaml"))
                        else:
                            return
                    else:
                        return
                else:
                    #$no point asking if we can delete the file if it is the same anyway.  still want to overwrite as part of checking permissions.
                    #$may not overwrite later
                    if self.check_file_removal([os.path.join(output_folder, "rate_settings.yaml")], ask_permission = False):
                        settings.freeze(os.path.join(output_folder, "rate_settings.yaml"))
                    else:
                        return
        else:
            settings.freeze(os.path.join(output_folder, "rate_settings.yaml"))
        
        
        
        #$then need to check if the files exist. if so warn the user. 
        #$extract should be analyzed after the user chooses the mzmls. also we're leaving graphs out of this
        no_extract_list = [w for w in worklist if w !="Extract"]
        outputs_to_check = []
        for worklist_step in no_extract_list:
            step_object_dict[worklist_step].complete_filename(self.file_loc)
            outputs_to_check.append(step_object_dict[worklist_step].full_filename)

        proceed = self.check_file_removal(outputs_to_check)
        if not proceed: 
            return
        
        #$now we need to get input and do the work. each step can only occur 
        #$once and they occur in order. so we will write them in order
        previous_output_file = "" 
        extracted_files =[]
        #$if the user is doing a full worklist it makes since to run the 2nd step
        #$(requesting their input) first so we don't bother stall partway through
        #$to ask for input.  this boolean governs that
        make_table_in_order = True
        for analysis_step in worklist:
            if analysis_step == "Extract":
                #$no if for this one, if extract is here it is the start
                id_file = self.collect_single_file("ID", "Extract", "CSV (*.csv)")
                if id_file == "": return
                #$always check if is good since it is first
                infile_is_good = self.check_input(step_object_dict[analysis_step],
                                                                   id_file)
                if not infile_is_good:  return
                
                mzml_files = self.collect_multiple_files("Centroided Data",
                                analysis_step, "mzML (*.mzML)")
                if mzml_files ==[]: return
                
                mzml_filenames = [os.path.basename(filename) for filename in
                                 mzml_files]
                extracted_files = [filename.replace(".mzML", ".tsv") for
                                        filename in mzml_filenames]
                extracted_files = [os.path.join(output_folder, filename) for
                                   filename in extracted_files]
                if settings.use_chromatography_division != "No":
                    extracted_intermediate_files = [filename.replace(".tsv", "_no_division.tsv") for
                                                    filename in extracted_files]
                else:
                    extracted_intermediate_files = extracted_files
                needed_files = list(set(extracted_files + extracted_intermediate_files))
                #$won't know what files to check until this point so need to check again here
                #$and ask permission to delete them of course
                proceed = self.check_file_removal(needed_files)
                if not proceed:  return
                #$need to run the table if necessary. taken from the 
                #$"Provide Time and Enrichment" elif lower down
                if "Provide Time and Enrichment" in worklist:
                    previous_output_file = step_object_dict[
                        "Provide Time and Enrichment"].full_filename
                    self.get_data_table = TimeWindow(self, 
                            extracted_files, previous_output_file)
                    self.get_data_table.exec_()
                    
                    if os.path.isfile(previous_output_file):
                        self.get_data_table2 = EnrichmentWindow(self, 
                                settings.min_allowed_timepoints_enrichment, 
                                settings.starting_enrichment_table_timepoints,
                                previous_output_file, settings.max_enrichment_allowed)
                        self.get_data_table2.exec_()
                    else: return
                    if os.path.isfile(previous_output_file):
                        enrichment_graph_folder = os.path.dirname(previous_output_file)
                        enrichment_graph_folder = os.path.join(enrichment_graph_folder,
                                                               "Enrichment_Graphs")
                        self.make_folder(enrichment_graph_folder)
                        enrichment_class = PerformEnrichmentClass(previous_output_file, enrichment_graph_folder, settings.graph_output_format)
                        enrichment_class.perform_calculations()
                        spline_error = enrichment_class.report_error()
                        if spline_error != "":
                            QtWidgets.QMessageBox.information(self, "Error", spline_error)
                            return 
                    else: return
                    
                    #$don't make the table twice
                    make_table_in_order = False
                    #$now that the table is done we need to confirm the user
                    #$hit the proceed button on the table (same check as in
                    #$elif analysis_step == "Theory Generation" )
                    if not os.path.exists(previous_output_file): return
                #$actually run the extraction
                for m in tqdm(range(len(mzml_files)), total=len(mzml_files), desc="Extracting mzml files: "):
                    extractor = Extractor(
                        id_path=os.path.join(self.file_loc, id_file),
                        mzml_path=mzml_files[m],
                        out_path=extracted_intermediate_files[m],
                        settings_path=rate_settings_file,
                    )
                    extractor.load()
                    extractor.run()
                    extractor.write()
                    #$need to delete classes when they're done or they may linger in RAM
                    del extractor
                #$this takes awhile so only bother with it if necessary
                if settings.use_chromatography_division != "No":
                    divider = ChromatographyDivider(settings_path=rate_settings_file,
                                                    input_paths=extracted_intermediate_files,
                                                    out_paths=extracted_files,
                                                    )
                    divider.divide()
                    del divider
            #$if we've already done this (make_table_in_order is false) there is no need
            #$ to 
            elif analysis_step == "Provide Time and Enrichment" and make_table_in_order:
                #$if coming right after a list
                if extracted_files == []:
                    extracted_files = self.collect_multiple_files(
                        "Extracted Data", 
                        analysis_step,
                        "TSV (*.tsv)"
                        )
                    if extracted_files == []: return
                    #$ensure the input files are good. only need to deal with
                    #$if the user just selected instead of extracting fresh
                    for e_file in extracted_files:
                        infile_is_good = self.check_input(
                            step_object_dict[analysis_step], e_file)
                        if not infile_is_good: return
                #$there are two tables that must be handled in order
                
                #$ now that we have the extracted files we can make a table
                #$the talbe will handle the output
                previous_output_file = step_object_dict[
                    analysis_step].full_filename
                self.get_data_table = TimeWindow(self, 
                        extracted_files, previous_output_file)
                self.get_data_table.exec_()
                #$here we can add the call to the next table. since
                #$ we are just using the same outfile can just keep checking
                #$ if it exists
                if os.path.isfile(previous_output_file):
                    self.get_data_table2 = EnrichmentWindow(self, 
                            settings.min_allowed_timepoints_enrichment, 
                            settings.starting_enrichment_table_timepoints,
                            previous_output_file, settings.max_enrichment_allowed)
                    self.get_data_table2.exec_()
                else: return
                
                #$now that we have all the data we can graph the results in case the user wishes to view the data
                if os.path.isfile(previous_output_file):
                        enrichment_graph_folder = os.path.dirname(previous_output_file)
                        enrichment_graph_folder = os.path.join(enrichment_graph_folder,
                                                               "Enrichment_Graphs")
                        self.make_folder(enrichment_graph_folder)
                        enrichment_class = PerformEnrichmentClass(previous_output_file, enrichment_graph_folder, settings.graph_output_format)
                        enrichment_class.perform_calculations()
                        spline_error = enrichment_class.report_error()
                        if spline_error != "":
                            QtWidgets.QMessageBox.information(self, "Error", spline_error)
                            return 
                            
                else: return
            
            #$now we need to merge the various extracted files into one file
            #$for future analysis
            elif analysis_step == "Combine Extracted Files":
                #$since the files are in the table can just read that in
                if previous_output_file == "":
                    previous_output_file = self.collect_single_file(
                        "time and enrichment", 
                        analysis_step, 
                        "spreadsheet (*.csv *.tsv)"
                    )
                    if previous_output_file == "": return
                    infile_is_good = self.check_input(
                        step_object_dict[analysis_step], previous_output_file)
                    if not infile_is_good: return
                #$else is to deal with a failed write from the previous table
                #$ don't need an error message just return
                elif not os.path.exists(previous_output_file): 
                    return
                
                #$final check to see if all of the files in the input table 
                #$still exist.  don't want to error out in the middle of 
                #$multiprocessing
                final_proceed =  self.check_files_from_files(
                    previous_output_file, 0)
                if not final_proceed: return
                combiner = CombineExtractedFiles(
                    enrichment_path=previous_output_file,
                    out_path= step_object_dict[analysis_step].full_filename,
                    settings_path=rate_settings_file,
                    needed_columns = step_object_dict["Calculate Delta by Enrichment"].required_columns
                )
                combiner.prepare()
                combiner.write()
                del combiner
                previous_output_file = step_object_dict[analysis_step].full_filename
                
            #$this step is to get all the numberical values necessary to do the fitting in the following steps
            elif analysis_step == "Calculate Delta by Enrichment":
                if previous_output_file == "":
                    previous_output_file = self.collect_single_file(
                        "combine output", 
                        analysis_step, 
                        "spreadsheet (*.csv *.tsv)"
                    )
                    if previous_output_file == "": return
                    infile_is_good = self.check_input(
                        step_object_dict[analysis_step], 
                        previous_output_file)
                    if not infile_is_good: return
                    #$not sure why this would happen but we'll put it here
                #$to avoid future error
                elif not os.path.exists(previous_output_file):
                    return
                
                enrich_delta = theoretical_enrichment_calculator(
                    prepared_data_path = previous_output_file, 
                    out_path = step_object_dict[analysis_step].full_filename,
                    settings_path = rate_settings_file
                    )
                enrich_delta.prepare()
                enrich_delta.write()
                del enrich_delta
                previous_output_file = step_object_dict[analysis_step].full_filename
                
            #$this step calls the fitter and runs it
            elif analysis_step == "Rate Calculation":
                
                if previous_output_file == "":
                    previous_output_file = self.collect_single_file(
                        "calculate delta by enrichment", 
                        analysis_step, 
                        "spreadsheet (*.csv *.tsv)"
                    )
                    if previous_output_file == "": return

                    
                    infile_is_good = self.check_input(
                        step_object_dict[analysis_step], 
                        previous_output_file)
                    if not infile_is_good: return
                
                #$need to get a graph folders and ensure they exist
                #$don't worry about overwriting files. not deleting all just causes confusion.
                GraphFolder_isotopes = os.path.join(self.file_loc, "Graph_Folder_Isotopes")
                self.make_folder(GraphFolder_isotopes)
                if settings.error_estimation == "none":
                    GraphFolder_optimization = ""
                else:
                    GraphFolder_optimization = os.path.join(self.file_loc, "Graph_Folder_Optimization")
                    self.make_folder(GraphFolder_optimization)
                
                ratecalc = RateCalculator(
                    model_path = previous_output_file,
                    out_path = step_object_dict[analysis_step].full_filename,
                    graph_folder_isotopes = GraphFolder_isotopes,
                    graph_folder_optimization = GraphFolder_optimization,
                    settings_path = rate_settings_file
                )
                ratecalc.calculate()
                ratecalc.write()
                del ratecalc
                previous_output_file = step_object_dict[analysis_step].full_filename
                
            #$last thing to do is to roll up the proteins
            elif analysis_step == "Combine Sequence Rates":
                if previous_output_file == "":
                    previous_output_file = self.collect_single_file(
                        "Rate Calculation", 
                        analysis_step, 
                        "spreadsheet (*.csv *.tsv)"
                    )
                    if previous_output_file == "": return

                    
                    infile_is_good = self.check_input(
                        step_object_dict[analysis_step], 
                        previous_output_file)
                    if not infile_is_good: return
                GraphFolder_averages = os.path.join(self.file_loc, "Graph_Folder_Protein_Averages")
                self.make_folder(GraphFolder_averages)
                
                combine_rate_calc = Peptides_to_Proteins(
                    model_path = previous_output_file,
                    out_path = step_object_dict[analysis_step].full_filename,
                    settings_path = rate_settings_file,
                    graph_folder_path = GraphFolder_averages
                    )
                combine_rate_calc.calculate()
                combine_rate_calc.write()
                del combine_rate_calc
                
        QtWidgets.QMessageBox.information(self, "Success", 
                "Analysis completed successfully")
        
        
       
    #$this checks which parts of the worklist table are checked and that there are no blanks    
    def check_table_checklist(self):
        current_worklist =[]
        error_check = []
        for i in range(self.RequestedStepsTable.rowCount()):
            if self.RequestedStepsTable.item(i,0).checkState() == QtCore.Qt.Checked:
                current_worklist.append(str(
                    self.RequestedStepsTable.item(i,0).text()))
                error_check.append(i)
        #$the point of the error check is to ensure ther are no gaps in the 
        #$checklist since then we will be missing critical info for the next
        #$step
        if error_check == []: #$max does not work on empty lists
            return current_worklist
        proper_length = max(error_check) - min(error_check) + 1
        if len (current_worklist) < proper_length:
            return ("There are gaps in your worklist.  Please check all boxes "
                    "between the first and last checked box")
        return current_worklist

    #$ask the user for detailed information.  mostly for detailing which
    #$files will be overwritten
    def large_text_question_for_use(self, title, infoText, detailedText):
        question = QtWidgets.QMessageBox(self)
        question.setWindowTitle("Question")
        question.setText(title)
        question.setIcon(QtWidgets.QMessageBox.Question)
        
        question.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        question.setDefaultButton(QtWidgets.QMessageBox.No)
        
        question.setDetailedText(detailedText)
        question.setInformativeText(infoText)
        
        question.setStyleSheet("QMessageBox{min-width:650 px;}")
        response = question.exec_()
        
        if response == QtWidgets.QMessageBox.Yes:
            return True
        else:
            return False
    
    #$give information to the user
    def large_text_information(self, title, text, detailedText):
        info = QtWidgets.QMessageBox(self)
        info.setWindowTitle(title)
        info.setText(text)
        info.setIcon(QtWidgets.QMessageBox.Information)
        
        info.setStandardButtons(QtWidgets.QMessageBox.Ok)
        info.setDefaultButton(QtWidgets.QMessageBox.Ok)
        
        info.setDetailedText(detailedText)
        
        info.setStyleSheet("QMessageBox{min-width:650 px;}")
        info.exec_()

    #$ask user for one file
    def collect_single_file(self, to_load, step_name, load_type):
        QtWidgets.QMessageBox.information(self, "Info", ("Please choose the {} "
                "file to load for {} step".format(to_load, step_name)))
        filename, file_type = QtWidgets.QFileDialog.getOpenFileName(self, 
                "Choose {} file to load".format(to_load), self.file_loc, 
                load_type)
        return filename
    
    #$ask user for potentially more than one file
    def collect_multiple_files(self, to_load, step_name, load_type):
        QtWidgets.QMessageBox.information(self, "Info", ("Please choose the {} "
                "files to load for {} step".format(to_load, step_name)))
        filenames, file_type = QtWidgets.QFileDialog.getOpenFileNames(self, 
                "Choose {} file to load".format(to_load), self.file_loc, 
                load_type)
        return filenames
    
    #$pull up the settings menu.
    def change_settings(self):
        self.set_menu = Rate_Setting_Menu(self, rate_settings_file)
        self.set_menu.show()
     
    #$checks if a file exists, and if we can actually do so.  if we can't
    #$an error will be triggered later so we'll sort that out here.  
    #$has the advantage of deleting the files so the user can't open a file between now and when we write it
    #$ask the user about this for everything but settings (that is dealt with better in the run_rate_worklist function)
    def check_file_removal(self, list_of_filenames, ask_permission = True):
        files_to_remove = []
        open_files = []
        #$find files that exist
        for filename in list_of_filenames:
            if os.path.exists(filename):
                files_to_remove.append(filename)
        # $let the user decide if we should continue.
        if files_to_remove != []:
            if ask_permission:
                proceed = self.large_text_question_for_use("Some files already exist and will be overwritten.",
                                                           "Do you still wish to proceed?",
                                                           "Files to be overwritten:\n" + ",\n".join(files_to_remove))
                if not proceed:
                    return False
            for filename in files_to_remove:
                try:
                    os.remove(filename)
                except PermissionError:
                    open_files.append(filename)
            if open_files != []:
                self.large_text_information("Error", "Some files cannot be overwritten.\n\n "
                                                     "They are likely open in another program. Please close "
                                                     "and try again.",
                                            "Files unable to be opened:\n" + ",\n".join(open_files))

                return False
        # $will return true if no files already exist or the user wants to
        # $overwrite and they can be removed so we have permission
        return True
      
    #$check a given input file has the headers we need
    def check_input (self, relevant_object, filename):
        has_needed_columns = relevant_object.check_input_file(filename)
        if not has_needed_columns:
            QtWidgets.QMessageBox.information(self, "Error", ("File {} is "
                "missing needed columns. Please correct and try again".format(
                filename)))
        return has_needed_columns
    
    #$there are cases (specifically the enrichment table referencing extracted files) where there
    #$is a chance that the files referenced may not exist
    #$this will check for that
    def check_files_from_files(self, input_file, filename_column):
        with open (input_file, 'r') as infile:
            if input_file[-4:] == ".tsv":
                reader = csv.reader(infile, delimiter = "\t")
            elif input_file[-4:] == ".csv":
                reader = csv.reader(infile)
            next(reader)#$skip header
            #$may wish to add a check that rows are of appropriate length
            for row in reader:
                #$ the != "" is needed for cases where the water data is longer than the measured timepoints
                if row[filename_column] != "" and not os.path.exists(row[filename_column]):
                    QtWidgets.QMessageBox.information(self, "Error", (f"File {row[filename_column]} "
                        "could not be found. Please correct input file and try "
                        "again"
                        ))
                    return False
        return True
    
    #$ensures a folder exists by making it if it does not.
    def make_folder(self, folder, non_graph = False):
        if non_graph: #$ don't care aobut overwriting graph folder, that is necessary.   the main output folder might have necessary things in it the user wants to keep
            if not os.path.isdir(folder):
                os.makedirs(folder)    
        else:
            if os.path.isdir(folder):
                rmtree(folder)
            os.makedirs(folder)
            
    
#$since we have to load settings in each file, and need a way to adjust
#$settings, we'll 
def make_temp_file(filename, new_filename):
    copyfile(filename, new_filename)

#$the actual main to run the Gui
if __name__ == '__main__':
    #$needed for windows multiprocessing which will happen at some point
    import sys
    mp.freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationDisplayName("DeuteRater-H")
    app.setApplicationName("DeuteRater-H")
    app.setWindowIcon(QtGui.QIcon(os.path.join(location, "resources", "Logo_64_clean.PNG")))
    gui_object = MainGuiObject(None)
    gui_object.show()
    app.exec_()
