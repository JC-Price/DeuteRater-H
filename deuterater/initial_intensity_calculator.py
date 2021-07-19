# -*- coding: utf-8 -*-
"""
delta_enrichment_calculator

this calls emass and fits the change in abundance and spacing as  a function 
of enrichment
"""



from tqdm import tqdm
import pandas as pd
import multiprocessing as mp
import numpy as np

from pathlib import Path
from functools import partial

from utils.emass import emass
import deuterater.settings as settings


max_isos = 5 #$constant based on the n_isos based on the mass (done in the extractor)
p0_guess = 1 #$seems to work for most fits. if it causes problems we can adjust

class theoretical_enrichment_calculator(object):
    def __init__(self, prepared_data_path, out_path, settings_path):
        settings.load(settings_path)
        self.settings_path = settings_path
        
        self.prepared_data_path = Path(prepared_data_path)
        self.out_path = out_path
        
        if self.prepared_data_path.suffix == '.tsv':
            self.data_df = pd.read_csv(
                filepath_or_buffer=str(self.prepared_data_path),
                sep='\t'
            )
        elif self.prepared_data_path.suffix == '.csv':
            self.data_df = pd.read_csv(
                filepath_or_buffer=str(self.prepared_data_path),
                sep=','
            )
        if settings.recognize_available_cores is True:
            self._n_partitions = mp.cpu_count()
        else:
            self._n_partitions = settings.n_processors
            
    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )    
        
    def prepare(self):
        unique_sequnces_df = self.data_df.drop_duplicates(subset = ["Sequence"])
    
        new_columns = theoretical_enrichment_calculator._make_new_columns()
        func = partial(theoretical_enrichment_calculator._individual_process, 
                       new_columns = new_columns,
                       minimum_n_value = settings.min_allowed_n_values,
                       minimum_sequence_length = settings.min_aa_sequence_length)
        df_split = np.array_split(unique_sequnces_df, self._n_partitions)
                     
        mp_pools = mp.Pool(self._n_partitions)
        
        final_df = pd.concat(tqdm(mp_pools.imap(func, df_split), total = self._n_partitions),axis =1)
        mp_pools.close()
        mp_pools.join()
        final_df = final_df.T
        final_df = final_df.set_index("Sequence")
        self.model = pd.merge(self.data_df, final_df, left_on= "Sequence", right_index = True)
        
    @staticmethod
    def _individual_process(df, new_columns, 
                            minimum_n_value,minimum_sequence_length):
         variable_list = []
         for row in df.itertuples():
            output_series = pd.Series(index = new_columns, dtype = "object")
            output_series["Sequence"] = row.Sequence
            #$emass takes longer now as does the graphing. 
            #$if there is a reason to drop let's do it now
            if len(row.Sequence) < minimum_sequence_length:
                variable_list.append(theoretical_enrichment_calculator._error_message_results(
                    f"Sequence is less than {minimum_sequence_length} amino acids",
                    output_series))
                continue
            if row.literature_n < minimum_n_value:
                variable_list.append(theoretical_enrichment_calculator._error_message_results(
                    f"less than {minimum_n_value} labeling sites",
                    output_series))
                continue
            intensity_values = \
                theoretical_enrichment_calculator._fit_emass(row.cf,
                      row.n_isos
                )

            output_series["Theoretical Unlabeled Normalized Abundances"] = ", ".join(intensity_values)
            variable_list.append(output_series)
         return(pd.concat(variable_list,axis =1))
     

    #$if an error happens it is most efficient to have a easy function
    def _error_message_results(error_message, output_series):
        #$don't need to know which names are which or how many columns there are, 
        #$just need python to fill all non-Sequence columns
        #$position 0 is the sequence name which we don't wish to overwrite
        for index_name in output_series.index[1:]:
            output_series[index_name] = error_message
        return output_series

                
    

    #$calculate unlabeled intensity, if we need to return  m/z values or
    #$adjust for different n_values, do it here or in emass itself.
    def _fit_emass(sequence, n_isos):
        intensity_values = emass(
                    sequence,
                    n_isos
                )
        return [str(i) for i in intensity_values]

    #$this creates the header for the variables. is a function in case we need 
    #$to add various columns  (if we want to graph emass output or something)
    @staticmethod
    def _make_new_columns():
        new_columns = ["Sequence"]
        new_columns.extend(["Theoretical Unlabeled Normalized Abundances"])
        return new_columns
    