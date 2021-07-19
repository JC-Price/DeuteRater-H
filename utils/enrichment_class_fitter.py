# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 11:17:40 2020

deal with splines and fitting the equations
will also do any error checking if necessary 

@author: JCPrice
"""

import os
import numpy as np
import pandas as pd
import scipy.interpolate as si

from utils.graphing_tools import enrichment_graph

#$the spline equation will be separate for easier calling in other modules
def spline_interp(t, d):
    std_d = np.std(d)
    w_y = [3/std_d for i in range(len(t))]
    spl = si.splrep(t, d, w=w_y)
    return lambda x: si.splev(x, spl)

#$this just holds some data for ease
class subject_data(object):
    def __init__(self, subject_name, x, y):
        self.subject_name = subject_name
        self.x = x
        self.y = y
    def calculate_spline(self):
        self.spline =  spline_interp(self.x, self.y)
    
    def calculate_theory_data(self):
        self.theory_x = np.arange(0, max(self.x), .1)
        self.theory_y = self.spline(self.theory_x)

old_header =["Subject ID Enrichment", "Time Enrichment", "Enrichment"]
#$we'll use a class for consistency with other parts
class PerformEnrichmentClass(object):
    def __init__(self, input_filename, graph_folder_name):
        temp_df = pd.read_csv(input_filename, sep = "\t")
        self.subject_ids = list(temp_df[old_header[0]].unique())
        self.load_previous_data(input_filename)
        self.graph_folder_name = graph_folder_name
        
    def perform_calculations(self):
        for subject_key in self.subject_dictionary.keys():
            temp_subject_class = self.subject_dictionary[subject_key]
            temp_subject_class.calculate_spline()
            temp_subject_class.calculate_theory_data()
            temp_graph_name = os.path.join(self.graph_folder_name, temp_subject_class.subject_name +".pdf")
            enrichment_graph(temp_subject_class.x, temp_subject_class.y, 
                             temp_subject_class.theory_x, temp_subject_class.theory_y,
                             temp_subject_class.subject_name, temp_graph_name)
            
        
    #$need to load the file and prepare for fitting
    def load_previous_data(self, input_filename):
        self.subject_dictionary = {}
        df = pd.read_csv(input_filename, delimiter=("\t"))
        df = df[old_header]
        for subject in self.subject_ids:   
            temp_df = df[df[old_header[0]] == subject]
            #$need to pass numpy array values for the fits
            self.subject_dictionary[subject] = subject_data(
                subject,
                temp_df[old_header[1]].to_numpy(),
                temp_df[old_header[2]].to_numpy()
                )
        
if __name__ == '__main__':
    test = PerformEnrichmentClass("", "")
    