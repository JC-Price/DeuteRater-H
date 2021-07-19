import yaml
import traceback
import os

from pathlib import Path


from utils.exc import InvalidSettingsWarning
from utils.exc import InvalidSettingsError  # noqa: 401 

# TODO: How would I dynamically load a different settings file?
# TODO: ^^^This really should be able to be passed in
# TODO: should different steps have different settings files?
# TODO: add reasonable constraints on settings
# TODO: verbose exception output
# TODO: add documentation on what this file is for
# TODO: Shorten variable names where possible
# TODO: add error checking where applicable
# TODO: set type annotations and default values
# TODO: determine good defaults
# TODO: type annotations for path variable
# TODO: Discuss this style vs settings class style
# TODO: Figure out how to reduce redundancy. Like a better singleton

# NOTE: we can use LibYAML c bindings if we need more speed

location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
resource_location = os.path.join(location, "resources")

debug_level: int
recognize_available_cores: bool
n_processors: int
id_file_rt_unit: str
trim_ids_to_mzml_bounds: bool
chunk_size: int
chunking_method_threshold: int
time_window: float
ppm_window: int
label_key: str
aa_labeling_sites_path: str
peak_lookback: int
peak_lookahead: int
baseline_lookback: int
peak_ratio_denominator: int
zscore_cutoff: int
mz_proximity_tolerance: float
min_aa_sequence_length: int
min_allowed_n_values: int
starting_enrichment_table_timepoints: int
error_estimation: str
min_non_zero_timepoints_rate: int
min_allowed_timepoints_enrichment: int
max_rate_check_step_size: float
minimum_allowed_sequence_rate: float
maximum_allowed_sequence_rate: float
minimum_sequences_to_combine_for_protein_rate: int
lowest_allowed_norm_isotope: float
highest_allowed_norm_isotope: float
m0_decreasing_allowed_noise: float
median_absolute_residuals_cutoff_single_point: float
median_absolute_residuals_cutoff_two_points: float
median_absolute_residuals_cutoff_general: float
desired_points_for_optimization_graph: int

# TODO: add quick explanation of how this works, inc. 'global' doc link
def load(settings_path):
    # NOTE: Look at the python documentation for the 'global' statement if you
    #       Want to understand how this module works
    try:
        settings_path = Path(settings_path)
        with settings_path.open('r') as f:
            s = yaml.load(f, Loader=yaml.FullLoader)
        global debug_level
        debug_level = s['debug_level']
        if debug_level not in [0, 1, 2]:
            raise InvalidSettingsWarning(
                'Invalid debug level value given'
            )
            print('Running with debug level 0')
            debug_level = 0

        global recognize_available_cores
        recognize_available_cores = s['recognize_available_cores']

        global n_processors
        n_processors = s['n_processors']

        global id_file_rt_unit
        id_file_rt_unit = s['id_file_rt_unit']

        global trim_ids_to_mzml_bounds
        trim_ids_to_mzml_bounds = s['trim_ids_to_mzml_bounds']

        global chunk_size
        chunk_size = s['chunk_size']

        global chunking_method_threshold
        chunking_method_threshold = s['chunking_method_threshold']

        global time_window
        time_window = s['time_window']

        global ppm_window
        ppm_window = s['ppm_window']

        global label_key
        label_key = s["label_key"]
        
        global aa_labeling_sites_path
        aa_labeling_sites_path = os.path.join(resource_location, s["aa_labeling_sites_path"])

        global peak_lookback
        peak_lookback = s['peak_lookback']

        global peak_lookahead
        peak_lookahead = s['peak_lookahead']

        global baseline_lookback
        baseline_lookback = s['baseline_lookback']

        global min_envelopes_to_combine
        min_envelopes_to_combine = s['min_envelopes_to_combine']

        global peak_ratio_denominator
        peak_ratio_denominator = s['peak_ratio_denominator']

        global zscore_cutoff
        zscore_cutoff = s['zscore_cutoff']
        
        global mz_proximity_tolerance
        mz_proximity_tolerance = s["mz_proximity_tolerance"]
        
        global min_aa_sequence_length
        min_aa_sequence_length = s["min_aa_sequence_length"]
        
        global min_allowed_n_values
        min_allowed_n_values = s["min_allowed_n_values"]
        
        global starting_enrichment_table_timepoints
        starting_enrichment_table_timepoints = s["starting_enrichment_table_timepoints"]
        
        global error_estimation
        error_estimation = s["error_estimation"]
        
        global min_non_zero_timepoints_rate
        min_non_zero_timepoints_rate = s["min_non_zero_timepoints_rate"]
        
        global min_allowed_timepoints_enrichment
        min_allowed_timepoints_enrichment = s["min_allowed_timepoints_enrichment"]
        
        global max_rate_check_step_size
        max_rate_check_step_size = s["max_rate_check_step_size"]
        
        global minimum_allowed_sequence_rate
        minimum_allowed_sequence_rate = s["minimum_allowed_sequence_rate"]
        
        global maximum_allowed_sequence_rate
        maximum_allowed_sequence_rate = s["maximum_allowed_sequence_rate"]
        
        global minimum_sequences_to_combine_for_protein_rate
        minimum_sequences_to_combine_for_protein_rate = s["minimum_sequences_to_combine_for_protein_rate"]
        
        global lowest_allowed_norm_isotope
        lowest_allowed_norm_isotope = s["lowest_allowed_norm_isotope"]
        
        global highest_allowed_norm_isotope
        highest_allowed_norm_isotope = s["highest_allowed_norm_isotope"]
        
        global m0_decreasing_allowed_noise
        m0_decreasing_allowed_noise = s["m0_decreasing_allowed_noise"]
        
        global median_absolute_residuals_cutoff_single_point
        median_absolute_residuals_cutoff_single_point = s["median_absolute_residuals_cutoff_single_point"]
        
        global median_absolute_residuals_cutoff_two_points
        median_absolute_residuals_cutoff_two_points = s["median_absolute_residuals_cutoff_two_points"]
        
        global median_absolute_residuals_cutoff_general
        median_absolute_residuals_cutoff_general = s["median_absolute_residuals_cutoff_general"]
        
        global desired_points_for_optimization_graph
        desired_points_for_optimization_graph = s["desired_points_for_optimization_graph"]

    except Exception as e:
        print(e)
        traceback.print_tb(e.__traceback__)


def freeze(path=None, settings_dict = None):
    if not settings_dict:
        settings_dict = {
            'debug_level': debug_level,
            'recognize_available_cores': recognize_available_cores,
            'n_processors': n_processors,
            'id_file_rt_unit': id_file_rt_unit,
            'trim_ids_to_mzml_bounds': trim_ids_to_mzml_bounds,
            'chunk_size': chunk_size,
            'chunking_method_threshold': chunking_method_threshold,
            'time_window': time_window,
            'ppm_window': ppm_window,
            "label_key": label_key,
            "aa_labeling_sites_path": aa_labeling_sites_path,
            'peak_lookback': peak_lookback,
            'peak_lookahead': peak_lookahead,
            'baseline_lookback': baseline_lookback,
            'min_envelopes_to_combine': min_envelopes_to_combine,
            'peak_ratio_denominator': peak_ratio_denominator,
            'zscore_cutoff': zscore_cutoff,
            "min_aa_sequence_length": min_aa_sequence_length,
            "mz_proximity_tolerance":mz_proximity_tolerance,
            "min_allowed_n_values": min_allowed_n_values,
            "starting_enrichment_table_timepoints": starting_enrichment_table_timepoints,
            "error_estimation": error_estimation,
            "min_non_zero_timepoints_rate": min_non_zero_timepoints_rate,
            "min_allowed_timepoints_enrichment": min_allowed_timepoints_enrichment,
            "max_rate_check_step_size": max_rate_check_step_size,
            "minimum_allowed_sequence_rate": minimum_allowed_sequence_rate,
            "maximum_allowed_sequence_rate": maximum_allowed_sequence_rate,
            "minimum_sequences_to_combine_for_protein_rate": minimum_sequences_to_combine_for_protein_rate,
            "lowest_allowed_norm_isotope": lowest_allowed_norm_isotope,
            "highest_allowed_norm_isotope": highest_allowed_norm_isotope,
            "m0_decreasing_allowed_noise": m0_decreasing_allowed_noise,
            "median_absolute_residuals_cutoff_single_point": median_absolute_residuals_cutoff_single_point,
            "median_absolute_residuals_cutoff_two_points": median_absolute_residuals_cutoff_two_points,
            "median_absolute_residuals_cutoff_general": median_absolute_residuals_cutoff_general,
            "desired_points_for_optimization_graph": desired_points_for_optimization_graph
            
        }
    if path:
        with open(path, 'w') as frozen_settings_file:
            yaml.dump(
                data=settings_dict,
                stream=frozen_settings_file,
                canonical=False
            )
    else:
        print(yaml.dump(data=settings_dict, canonical=False))

