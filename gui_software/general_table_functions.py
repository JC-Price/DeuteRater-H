# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:46:01 2020

@author: JCPrice

this will have functions needed for qdialog windows and the set up table
this could be done with inheritance but trying to fiddle with inheriting
from different uic files, making a general class and so on is much more
work and time than just a list of functions we can just call

"""

from PyQt5 import QtWidgets, QtCore

def set_up_table(table_object, row_names, repeats = 1, rows_are_paths = False):
    table_object.setRowCount(len(row_names)*repeats)
    current_row = 0
    for name in row_names:
        for i in range(repeats):
            if rows_are_paths:
                name_item  = QtWidgets.QTableWidgetItem(str(name.stem))
            else:
                name_item = QtWidgets.QTableWidgetItem(str(name))
            #$ this makes the item uneditable so the user can't mess it up
            #$from https://stackoverflow.com/questions/17104413/pyqt4-how-to-
            #$select-table-rows-and-disable-editing-cells anser 2 accessed 9/25/20
            name_item.setFlags(QtCore.Qt.ItemIsSelectable | 
                               QtCore.Qt.ItemIsEnabled)
            table_object.setItem(current_row, 0, name_item)
            for i in range(1, table_object.columnCount()):
                #$makes an empty string in all other cells
                table_object.setItem(current_row, i, 
                                QtWidgets.QTableWidgetItem(""))
            current_row += 1
            
#$ gets the top left selected cell or all selected rows or columns
def HighlightedCells(table_object, all_data = False, no_names = True):
    rows = []
    columns = []
    for selected_cell in table_object.selectedIndexes():
        rows.append(int(selected_cell.row()))
        columns.append(int(selected_cell.column()))
    #$ don't let the user paste over or delete the names (copy is fine)
    if no_names:
        columns = [c for c in columns if c!= 0]
    if all_data:
        #$only get unique rows and colums and sort them (set is unsorted)
        rows = sorted(list(set(rows)))
        columns = sorted(list(set(columns)))
        return rows, columns
    else:
        return rows[0], columns[0]

def Clear_Contents(table_object):
    selected_rows, selected_columns = HighlightedCells(table_object, True)
    for r in selected_rows:
        for c in selected_columns:
            table_object.item(r,c).setText("")
            
def Copy_to_Clipboard(table_object):
        rows, columns = HighlightedCells(table_object, 
                                         all_data = True, no_names = False)
        string_for_clipboard = ""
        for r in rows:
            for c in columns:
                string_for_clipboard += str(table_object.item(r,c).text())
                string_for_clipboard += "\t"
            #$ swap out last tab for a newline
            string_for_clipboard = string_for_clipboard[:-1] + "\n"
        clipboard = QtWidgets.QApplication.clipboard()
        #$ trim off the final new line and get it on the clipboard
        clipboard.setText(string_for_clipboard[:-1]) 
    
def Paste(table_object):
    start_row, start_column = HighlightedCells(table_object, no_names = False)
    #$ do not paste over the first column
    if start_column == 0:
        return
    text = str(QtWidgets.QApplication.clipboard().text())
    text_lines = text.split("\n")
    text_grid = [line.split("\t") for line in text_lines]
    #excel sometimes adds whitespace to end of clipboard.  must remove 
    if text_grid[-1] == ['']: del text_grid[-1]
    line_lengths = [len(text_list) for text_list in text_grid]
    needed_columns = max(line_lengths)
    needed_rows =len(text_grid)
    #$ check that we don't need more rows or columns that we have
    #$ separate lines for easy readability
    if start_row + needed_rows > table_object.rowCount():
        return
    if start_column + needed_columns > table_object.columnCount():
        return
    for c in range(needed_columns):
        for r in range(needed_rows):
            table_object.item(start_row +r,
                start_column + c).setText(text_grid[r][c])
            
#$does some basic error checking for numerical data           
def basic_number_error_check(text_value, column_name, row_number):
        append_to_error = " at \"{}\" column, row {}. Correct to proceed.".format(
            column_name, row_number+1)
        if text_value == "":
            return "Blank present" + append_to_error
        try:
            num_value = float(text_value)
        except ValueError:
            return "Non-numerical value" + append_to_error
        if num_value < 0:
            return "Negative value" + append_to_error
        #$we could also check if it is under some maximum, but for now
        #$trust the user
        return num_value