# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:10:22 2020

@author: JCPrice

holds  basic functions for graphing. should be called as needed, no need to 
run as the main.
we'll call directly from the calculator, which will simplify the code here
by a significant amount
"""

import matplotlib.pyplot as plt
import numpy as np


error_line_symbol = 'k--'
data_points_symbol = 'ro'

multi_colors=['r','g','b','c','m', 'y']

#$need to graph the results of the binomial fit
#$adapted from code by Christian in the Transtrum lab
def graph_rate_results(n_isos, save_file_name, time_data, 
                       normed_isotope_data, mval, theory_times, theory_zero,
                       subject_sequence):
    plt.title(f"Isotope levels {subject_sequence}")
    plt.xlabel("Time")
    plt.ylabel("Fraction")
    for j in range(n_isos):
        plt.plot(theory_times, mval[:,j], multi_colors[j]+'-', label = "M"+str(j))
    for j in range(n_isos):
         plt.plot(time_data,normed_isotope_data[:,j], multi_colors[j]+"o")
         plt.plot(0, theory_zero[j], multi_colors[j] + "d")
    plt.legend(loc = "upper right")
    plt.savefig(save_file_name)
    plt.clf()

#$need to graph the optimization of error for troubleshooting purposes
#$adapted from code by Christian in the Transtrum lab
def graph_optimization_of_error(k_value, theoretical_k, cost, save_file_name, subject_sequence, legend_name):
    plt.title(f"Finding optimal k {subject_sequence}")
    plt.ylabel("Sum-square Error")
    plt.xlabel("k")
    plt.plot(theoretical_k,cost,"bo-", label = legend_name)
    plt.plot(np.repeat(k_value,2), [min(cost),max(cost)])
    plt.legend()
    plt.savefig(save_file_name)
    plt.clf()
    
def enrichment_graph(x_values, y_values, predicted_x, 
                     predicted_y, subject_name, save_file):
    plt.title(subject_name)
    plt.xlabel("Time")
    plt.ylabel("Enrichment")
    plt.plot(x_values, y_values, data_points_symbol)
    plt.plot(predicted_x, predicted_y, error_line_symbol)
    
    plt.savefig(save_file)
    plt.clf()
    
def graph_average_boxplots(values, save_file_name, subject_protein_name):
    #$give the values some jitter so they are easily visible from each other
    #$jitter is entirely random.  if we decide to remove, just make all  x_values 1
    x_values = np.random.normal(1, .035, size = len(values))
    average = np.mean(values)

    plt.plot(x_values, values, 'r*', markersize=11, label = "Peptide Rates")
    #$since the linewidth is an argument, is needs a value, but if we try and declare it as a dict,
    #$it will not be defined.  declaring its valuea and forcing to dict gets around that
    plt.boxplot(values,
            whiskerprops = dict(linewidth=1.5), 
            boxprops= dict(linewidth=3.0), 
            capprops = dict(linewidth = 2.0), 
            medianprops = dict(linewidth = 2.0)
            )
    #$the boxplot does have an option to show the average, but getting it into the legend
    #$is irritating
    plt.plot([1.00], [average], 
             markerfacecolor = "deepskyblue", 
             markeredgecolor = "deepskyblue", 
             markersize=10, 
             marker = "d", 
             label = "Mean", 
             linestyle="None")
    plt.legend()
    plt.title(subject_protein_name)
    plt.xticks([]) #$remove x-axis labels to avoid confusion
    plt.ylabel('Sequence Rate', fontsize=12)
    plt.savefig(save_file_name)
    plt.clf()