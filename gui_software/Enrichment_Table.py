# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 12:07:49 2020

@author: JCPrice
"""

import os
import pandas as pd

from PyQt5 import uic, QtWidgets, QtCore, QtGui

import gui_software.general_table_functions as gtf
#import general_table_functions as gtf

#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ui_file = os.path.join(location, "ui_files", "Enrichment_Table.ui")

loaded_ui = uic.loadUiType(ui_file)[0]

#$ it is tricky to actually get the header out of the qtablewidget and they
#$ need different checks anyway so we'll just declare it here
current_columns = ["Subject ID", "Time", "Enrichment"]

class EnrichmentWindow (QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent = None, min_allowed_times = None, 
                 starting_num_timepoints = None, outfile = None):
        super(EnrichmentWindow, self).__init__(parent)
        #$allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint | 
                            QtCore.Qt.WindowMaximizeButtonHint)
        self.early_exit = True
        self.setupUi(self)
        
        self.outfile = outfile
        self.min_allowed_times = min_allowed_times
        temp_df = pd.read_csv(outfile, sep = "\t")
        self.subject_ids = list(temp_df[current_columns[0]].unique())
        #self.set_up_table(subject_ids, starting_num_timepoints)
        gtf.set_up_table(self.EnrichmentTable, self.subject_ids, 
                         starting_num_timepoints)
        
        #$ just set the minimum from a setting to avoid needing multiple
        #$if statements
        self.enrichments_for_subject.setMinimum(self.min_allowed_times)
        self.enrichments_for_subject.setValue(starting_num_timepoints)
        self.current_timepoints = starting_num_timepoints
        
        self.CancelButton.clicked.connect(self.close)
        self.ProceedButton.clicked.connect(self.check_and_save)
        self.UpdateTableButton.clicked.connect(self.update_table_size)
        
        #$make some easy shortcuts. setting up undo and redo are the hardest
        #$conly do backspace, del, copy and paste for right now
        QtWidgets.QShortcut(QtGui.QKeySequence('Backspace'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.EnrichmentTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Del'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.EnrichmentTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+c'),
                            self).activated.connect(lambda: gtf.Copy_to_Clipboard(
                                self.EnrichmentTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+v'),
                            self).activated.connect(lambda:gtf.Paste(
                                self.EnrichmentTable))

    #$functional
    def update_table_size(self):
        desired_value = int(self.enrichments_for_subject.value())
        #$no point in doing anything if not adding or deleting
        if desired_value == self.current_timepoints:
            return
        elif desired_value < self.current_timepoints:
            row_to_delete = []
            for s in range(len(self.subject_ids)):
                #$ offset finds the indicies of the first instance of the 
                #$subject id
                offset = self.current_timepoints * s
                row_to_delete.extend([offset + x for x in range(desired_value,
                                              self.current_timepoints)])
            #$go backwards to ensure that we don't move critical rows into 
            #$delete positions
            for index_value in reversed(row_to_delete):
                self.EnrichmentTable.removeRow(index_value)
        else:
            for s in range(len(self.subject_ids)):
                offset = desired_value * s
                for index_value in range(self.current_timepoints, desired_value):
                    true_index = index_value + offset
                    self.EnrichmentTable.insertRow(true_index)
                    name_item = QtWidgets.QTableWidgetItem(self.subject_ids[s])
                    name_item.setFlags(QtCore.Qt.ItemIsSelectable | 
                                   QtCore.Qt.ItemIsEnabled)
                    self.EnrichmentTable.setItem(true_index, 0, name_item)
                    for i in range(1, self.EnrichmentTable.columnCount()):
                        #$makes an empty string in all other cells
                        self.EnrichmentTable.setItem(true_index, i, 
                                        QtWidgets.QTableWidgetItem(""))
        self.current_timepoints = desired_value
    
   
    def check_and_save (self):
        #$don't have to check the ids, those are out of the user's control
        class_dict = {}
        for r in range(self.EnrichmentTable.rowCount()):
            #$first column is not editable so can be ignored
            current_id = str(self.EnrichmentTable.item(r,0).text())
            #$ can't just loop as with Time_Table.py since we need to analyze
            #$ together and need to store in a class
            time_string = self.EnrichmentTable.item(r,1).text()
            enrichment_string = self.EnrichmentTable.item(r,2).text()
            #$blanks are okay as long as all are blank
            if  time_string == "" and enrichment_string == "":
                continue
            time_value = gtf.basic_number_error_check(
                    time_string, current_columns[1], r)
            enrichment_value =  gtf.basic_number_error_check(
                    enrichment_string, current_columns[2], r)
            for value in [time_value, enrichment_value]: 
                if type(value) == str:
                    QtWidgets.QMessageBox.information(self, "Error", value)
                    return 
            if current_id in class_dict.keys():
                class_dict[current_id].add_data(time_value, enrichment_value)
            else:
                class_dict[current_id] = subject(current_id, time_value,
                                                    enrichment_value)
        #$if all are blank (or all for one subject are blank) it will hit
        #$continue and never be added.  this is somehting to complain about
        for id_value in self.subject_ids:
            if id_value not in class_dict.keys():
                QtWidgets.QMessageBox.information(self, "Error", ("Subject ID {}"
                            " has no values. Correct to proceed".format(id_value)))
                return 
        
        #$now that we have the data in the subject object and let it do the 
        #$error checking
        for sample_id in class_dict.keys():
            return_value = class_dict[sample_id].error_checking(
                self.min_allowed_times)
            if type(return_value) == str:
                QtWidgets.QMessageBox.information(self, "Error", return_value)
                return
        self.make_save_file(class_dict)
            
    def make_save_file(self, class_dict):
        #$csv would be faster but since the files may be of different lengths
        #$it would likely be more complex than necessary
        full_names = []
        full_times = []
        full_enrichment = []
        blank_column = [""]
        for key in class_dict.keys():
            full_names.extend([key for i in range(len(class_dict[key].times))])
            full_times.extend(class_dict[key].times)
            full_enrichment.extend(class_dict[key].enrichments)
        blank_column = blank_column * len(full_names)
        df = pd.read_csv(self.outfile, delimiter ="\t")
        df2 = pd.DataFrame({"Spacing Column 1": blank_column, 
                            "Subject ID Enrichment":full_names,
                            "Time Enrichment": full_times,
                            "Enrichment": full_enrichment
                            })
        df = pd.concat([df, df2], axis =1)
        df.to_csv(self.outfile, sep ="\t", index = False)
        self.early_exit = False
        self.close()    
    
    #$need to get rid of the output if exiting early to avoid problems 
    #$with the file exiting but not be completed for other tables and 
    #$analysis step
    def closeEvent(self, event):
        if self.early_exit:
            os.remove(self.outfile)
        event.accept()
        
   
#$we can shift this to useful classes or general_table_functions later, but
#$for now put it here
class subject(object):
    def __init__(self, name, time, enrichment):
        self.name = name
        self.times = [time]
        self.enrichments = [enrichment]
    def add_data(self, time, enrichment):
        self.times.append(time)
        self.enrichments.append(enrichment)
        
    #$ need to do a small number of basic tests
    def error_checking(self, minimum_time_points):
        if len(set(self.times))< minimum_time_points:
            return "Subject {} has only {} unique time points. {} are requied.".format(
                self.name, len(set(self.times)), minimum_time_points)
        #$for now disallow all enrichments for a sample being 0
        if len(set(self.enrichments)) == 1  and self.enrichments[0] == 0.0:
            return "Subject {} has no enrichment.".format(self.name)
        return 1


if __name__ == '__main__':
    #$needed for windows multiprocessing which will happen at some point
    import sys
    app = QtWidgets.QApplication(sys.argv)
    output_file = "C:\\Software\\Testing\\Deuterater-H_test\\test_table_out.tsv"
    gui_object = EnrichmentWindow(None, ["A", "B", "C"], 3, 5, output_file)
    gui_object.show()
    app.exec_()
