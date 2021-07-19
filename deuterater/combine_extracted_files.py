'''Theoretical Value Preparation

This purpose of this module is to calculate the expected theoretical values
of each identified analyte.

'''

import pandas as pd
import numpy as np  # noqa: 401

import traceback  # noqa: 401
from pathlib import Path
from functools import partial
import multiprocessing as mp

from tqdm import tqdm  # noqa: 401

import deuterater.settings as settings

import deuteconvert.peptide_utils as peputils

literature_n_name = "literature_n"

class CombineExtractedFiles():
    def __init__(self, enrichment_path, out_path, settings_path, needed_columns):
        settings.load(settings_path)
        self.settings_path = settings_path
        self.enrichment_path = Path(enrichment_path)
        self.needed_columns = needed_columns
        aa_label_df = pd.read_csv(settings.aa_labeling_sites_path, sep='\t')
        aa_label_df.set_index('study_type', inplace=True)
        self.aa_labeling_dict = aa_label_df.loc[settings.label_key, ].to_dict()
        
        
        #$pull in two sub tables from the output table
        if self.enrichment_path.suffix == '.tsv':
            self._file_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep='\t',
                usecols = ["Filename", "Time", "Subject ID"]
            )
            self._enrichment_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep='\t',
                usecols = ["Subject ID Enrichment", "Time Enrichment", "Enrichment"]
            )
        elif self.enrichment_path.suffix == '.csv':
            self._file_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep=',',
                usecols = ["Filename", "Time", "Subject ID"]
            )
            self._enrichment_data = pd.read_csv(
                filepath_or_buffer=str(self.enrichment_path),
                sep=',',
                usecols = ["Subject ID Enrichment", "Time Enrichment", "Enrichment"]
            )
            
        #$since the multiple sub tables can have different length, get rid
        #$of the rows that are empty
        self._file_data.dropna(inplace = True, how = "all")
        self._enrichment_data.dropna(inplace = True, how = "all")
        self._data_dict = self.collect_enrichment_data()
        
        if settings.recognize_available_cores is True:
            self._n_partitions = mp.cpu_count()
        else:
            self._n_partitions = settings.n_processors

        self._mp_pool = mp.Pool(self._n_partitions)
        self.out_path = out_path
        self.model = None

    def write(self):
        self.model.to_csv(
            path_or_buf=self.out_path,
            sep='\t',
            index=False
        )
        
    def collect_enrichment_data(self):
        data_dict = {}
        for subject, subject_df in self._enrichment_data.groupby("Subject ID Enrichment"):
            x_values = ", ".join([str(x) for x in subject_df["Time Enrichment"]])
            y_values = ", ".join([str(y) for y in subject_df["Enrichment"]])
            data_dict[subject] = [x_values, y_values]
        return data_dict
            

    def prepare(self):
        if settings.debug_level == 0:
            args_list = self._file_data.to_records(index=False).tolist()
            func = partial(CombineExtractedFiles._mp_prepare, self.settings_path, self._data_dict, self.aa_labeling_dict)
            results = list(
                tqdm(
                    self._mp_pool.imap_unordered(func, args_list),
                    total=len(self._file_data)
                )
            
            )
            
            
        elif settings.debug_level >= 1:
            print('Beginning single-processor theory preparation.')
            results = []
            for row in tqdm(self._file_data.itertuples(),
                            total=len(self._enrichment_df)):
                # TODO: how to handle functions. Default I would think
                df = pd.read_csv(filepath_or_buffer=row.file, sep='\t')
                df = CombineExtractedFiles._apply_filters(df)
                if literature_n_name not in df.columns:
                    if self.aa_labeling_dict != "":
                        df = df.apply(CombineExtractedFiles._calculate_literature_n, axis =1 , args = (self.aa_labeling_dict,))
                
                df['time'] = row.time
                df["sample_id"]  = row.sample_id
                df["Time Enrichment"] = self._data_dict[row.sample_id][0]
                df["Enrichment Values"] = self._data_dict[row.sample_id][1]
                results.append(df)

        self.model = pd.concat(results)
        self.model = self.model.drop(columns=['drop'])
        #$the above drop is rows that have a problem.  now we need to filter columns
        #$otherwise the carry forward increases file size quite a bit
        #$by doing here it should not affect anything.
        self.model = self.model[self.needed_columns]
        

        self._mp_pool.close()
        self._mp_pool.join()

    @staticmethod
    def _mp_prepare(settings_path, data_dict,  aa_labeling_dict, args):
        settings.load(settings_path)
        #file_path, time, enrichment = args
        file_path, time, sample_id = args
        df = pd.read_csv(filepath_or_buffer=file_path, sep='\t')
        df = CombineExtractedFiles._apply_filters(df)
        
        if literature_n_name not in df.columns:
            if aa_labeling_dict != "":
                df = df.apply(CombineExtractedFiles._calculate_literature_n, axis =1 , args = (aa_labeling_dict,))
        df['time'] = time
        df["sample_id"]  = sample_id
        df["Time Enrichment"] = data_dict[sample_id][0]
        df["Enrichment Values"] = data_dict[sample_id][1]
            
        return df
    
    
    @staticmethod
    def _calculate_literature_n(row, aa_labeling_dict):
        aa_counts = {}
        for aa in row["Sequence"]:
            if aa not in aa_counts.keys():
                aa_counts[aa] = 0
            aa_counts[aa] += 1
        literature_n = peputils.calc_add_n(aa_counts, aa_labeling_dict)
        row[literature_n_name] = literature_n
        return row

    # def _load(self):
    #     '''Pulls in the relevant data from the model.

    #     Parameters
    #     ----------
    #     model : :obj:`pandas.Dataframe`
    #         The unmodified data model

    #     Returns
    #     -------
    #     :obj:`pandas.Dataframe`
    #         A dataframe containing a copy of the relevant columns

    #     '''
    #     # TODO: Basic checks like whether the data looks right need done
    #     if not isinstance(self.model, pd.DataFrame):
    #         try:
    #             # TODO: csv/tsv flexibility?
    #             self.model = pd.read_csv(self.model, sep='\t')
    #         except Exception as e:
    #             # TODO: better exception logging
    #             print(e)
    #             traceback.print_tb(e.__traceback__)
    #             raise
    #     try:
    #         # TODO: csv/tsv flexibility?
    #         self._enrichment_df = pd.read_csv(self.enrichment_path, sep='\t')
    #     except Exception as e:
    #         # TODO: better exception logging
    #         print(e)
    #         traceback.print_tb(e.__traceback__)
    #         raise

    @staticmethod
    def _apply_filters(df):
        '''filters the internal dataframe

        This function does not modify the dataframe in place.

        Parameters
        ----------
        df : :obj:`pandas.Dataframe`
            The internal dataframe of the theory_value_prep function

        Returns
        -------
        :obj:`pandas.Dataframe`
            The filtered dataframe. Does not modify in place.
        '''
        # This is 'clean_up_data' in the old deuterater
        # This is a
        data = df.dropna(
            axis='index',
            subset=['mzs', 'abundances']
        ).copy()
        data['drop'] = False
        for row in data.itertuples():
            mask = ((data['mz'] - row.mz).abs() <
                    settings.mz_proximity_tolerance)
            data.loc[mask, 'drop'] = True
        data.to_csv("C:\\Data\\test_deuterater_h_output\\final_algorithm_initial_test\\deuterater_.6_results\\test_errors\\check_filters.csv")
        data = data[~data['drop']]

        # TODO: Check to see if no data went through
        return data
    

def main():
    print('please use the main program interface')


if __name__ == '__main__':
    main()
