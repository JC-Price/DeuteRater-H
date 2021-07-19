# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 08:21:34 2020

@author: JCPrice

The point of this is to hold the time and enrichment table
any adjustments to the table should be made here or in the qt designer
"""
import os
import csv
import pandas as pd

from pathlib import Path

from PyQt5 import uic, QtWidgets, QtCore, QtGui

import gui_software.general_table_functions as gtf

#location = os.path.dirname(os.path.abspath(sys.executable))
location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ui_file = os.path.join(location, "ui_files", "Time_Table.ui")

loaded_ui = uic.loadUiType(ui_file)[0]

#$ it is tricky to actually get the header out of the qtablewidget and they
#$ need different checks anyway so we'll just declare it here
current_columns = ["Filename", "Time", "Subject ID"]

class TimeWindow(QtWidgets.QDialog, loaded_ui):
    def __init__(self, parent = None, filenames = [], outfile = None):
        super(TimeWindow, self).__init__(parent)
        #$allows a maximize button
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint | 
                            QtCore.Qt.WindowMaximizeButtonHint)
        
        self.setupUi(self)
        
        self.outfile = outfile
        self.filepaths = [Path(f) for f in filenames]
        gtf.set_up_table(self.TimeTable, self.filepaths, rows_are_paths = True)
        
        self.CancelButton.clicked.connect(self.close)
        self.ProceedButton.clicked.connect(self.check_and_save)
        
        #$make some easy shortcuts. setting up undo and redo are the hardest
        #$conly do backspace, del, copy and paste for right now
        QtWidgets.QShortcut(QtGui.QKeySequence('Backspace'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.TimeTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Del'),
                            self).activated.connect(lambda: gtf.Clear_Contents(
                                self.TimeTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+c'),
                            self).activated.connect(lambda: gtf.Copy_to_Clipboard(
                                self.TimeTable))
        QtWidgets.QShortcut(QtGui.QKeySequence('Ctrl+v'),
                            self).activated.connect(lambda:gtf.Paste(
                                self.TimeTable))
        
    #$we can just write out a .csv or .tsv with the csv module
    #$first we need to do some checks
    def check_and_save(self):
        results =[current_columns]
        for r in range(self.TimeTable.rowCount()):
            #$filename is not editable so can be ignored
            #current_row = [str(self.TimeTable.item(r,0).text())]
            current_row = [self.filepaths[r].resolve()]
            for i in range(1, len(current_columns)-1):
                test_value = gtf.basic_number_error_check(
                    self.TimeTable.item(r,i).text(),
                    current_columns[i], r)
                if type(test_value) == str:
                    QtWidgets.QMessageBox.information(self, "Error", test_value)
                    return 
                current_row.append(test_value)
            #get the sample group name (no need for it to be )
            test_value, error_code = TimeWindow._basic_string_check(
                    self.TimeTable.item(r,len(current_columns)-1).text(),
                    current_columns[-1], r)
            if error_code:
                QtWidgets.QMessageBox.information(self, "Error", test_value)
                return
            current_row.append(test_value)
            results.append(current_row)
        #$ if we've gotten here, we're good to write out
        with open(self.outfile, "w", newline ='') as temp_out:
            writer = csv.writer(temp_out, delimiter ='\t')
            writer.writerows(results)
            #for row in results:
            #    writer.writerow(row)
        self.close()
     
    @staticmethod    
    def _basic_string_check(text_value, column_name, row_number):
        append_to_error = " at \"{}\" column, row {}. Correct to proceed.".format(
            column_name, row_number)
        if text_value == "":
            return "Blank present" + append_to_error,True
        else: 
            return text_value, False
    

if __name__ == '__main__':
    #$needed for windows multiprocessing which will happen at some point
    import sys
    app = QtWidgets.QApplication(sys.argv)
    output_file = "C:\\Software\\Testing\\DeuteRater_Initial\\test_table_out.csv"
    gui_object = TimeWindow(None, ["A", "B", "C"], output_file)
    gui_object.show()
    app.exec_()
