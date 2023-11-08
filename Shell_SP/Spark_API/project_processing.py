import json
import os
import time
from datetime import datetime, timedelta
from typing import List, Generator, Tuple
from warnings import warn

import numpy as np
import openpyxl as pyxl
from pymcr.mcr import McrAR
import scipy
from scipy import interpolate
from scipy import signal
from scipy.optimize import curve_fit
from scipy.stats import norm as stat_norm
from scipy.stats import truncnorm as stat_truncnorm

import data_repository_interface as dri
import spectrometer_project as sp
from constants import WELL_IDS, WELLPLATE_COLUMNS, WELLPLATE_ROWS

_WELLPLATE_COLUMNS = [f":{num}" for num in WELLPLATE_COLUMNS]
_WELLPLATE_ROWS = [f"{letter}:" for letter in WELLPLATE_ROWS]

SPARK_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
MANIFEST_TIME_FORMAT = sp.TIME_FORMAT  # "%Y-%m-%d--%H:%M:%S"
REF_TIME = datetime(2020, 1, 1)


class SkipIteration(Exception):
    """ Used for skipping an outer loop iteration from a nested loop """
    pass


class DynamicInterpolator:
    """ Saves results from prior interpolations to save time """
    def __init__(self):
        """ Creates a dynamic programming version of interp1d & bisplrep """
        self.library = dict()

    def interp1d(self, x, y, **kwargs):
        """
        Initialize a 1-D linear interpolation class

        Ref: :class:`interpolate.interp1d`

        :param x: Any
        :param y: Any
        :keyword kind: str = 'linear'
        :keyword axis: int = -1
        :keyword copy: bool = True
        :keyword bounds_error: Any = None
        :keyword fill_value: Any = np.nan
        :keyword assume_sorted: bool = False
        :return: interpolate.interp1d(x, y, **kwargs)
        """
        key = ('interp1d', tuple(x), tuple(y))
        if key in self.library:
            return self.library[key]
        new_eval = interpolate.interp1d(x, y, **kwargs)
        self.library[key] = new_eval
        return new_eval

    def bisplrep(self, x1, x2, y, **kwargs):
        """
        Find a bivariate B-spline representation of a surface.

        Ref: :class:`scipy.interpolate.bisplrep`

        Given a set of data points (x[i], y[i], z[i]) representing a surface z=f(x,y), compute a B-spline
        representation of the surface. Based on the routine SURFIT from FITPACK

        :param x1: ndarray
        :param x2: ndarray
        :param y: ndarray
        :keyword w: ndarray = None
        :keyword xb: float = None
        :keyword xe: float = None
        :keyword yb: float = None
        :keyword ye: float = None
        :keyword kx: int = 3
        :keyword ky: int = 3
        :keyword task: int = 0
        :keyword s: float = None
        :keyword eps: float = 1e-16
        :keyword tx: ndarray = None
        :keyword ty: ndarray = None
        :keyword full_output: int = 0
        :keyword nxest: int = None
        :keyword nyest: int = None
        :keyword quiet: int = 1
        :return: scipy.interpolate.bisplrep(x1, x2, y, **kwargs)
        """
        key = ('bisplrep', tuple(x1), tuple(x2), tuple(y))
        if key in self.library:
            return self.library[key]
        new_eval = scipy.interpolate.bisplrep(x1, x2, y, **kwargs)
        self.library[key] = new_eval
        return new_eval


dynamic_scipy = DynamicInterpolator()


def cast_well_data_to_table(project: sp.SparkProjectFileHandler, times: List[timedelta],
                            well_id, use_signal='y-axis (processed)') -> np.ndarray:
    """ Creates a table of data

    :param project: The project being read
    :param times: A list to populate the time axis of the table
    :param well_id: Which well (lookup in project)
    :param use_signal: Specifies which signal to use
    :return: A 2D numpy array (table[0, 1:] - the X axis, table[1:, 0] - the Time axis, table[1:, 1:] -  the Y data)
    """
    well_data, method_specs = project.get_data_by_well(well_id)
    n_time_points = len(well_data)
    lambda_values: list = well_data[0][3]
    n_wavelength_points = len(lambda_values)

    table = np.empty(shape=[n_time_points + 1, n_wavelength_points + 1], dtype=np.float64) * np.nan
    table[0, 1:] = lambda_values
    table[1:, 0] = np.cumsum([t.total_seconds() for t in times], dtype=np.float64)

    for i, d, m, x in well_data:
        table[int(i) + 1, 1:] = d[use_signal]

    return table


def gauss(x, *p):
    """
    Calculates a gaussian function based on mean, std.dev, and amplitude

    :param x: x
    :param p: [mean, std.dev, amp]
    :return: gauss(x; mean, std.dev, amp)
    """
    mu, sigma, amp = p
    return amp * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def norm(x, mu, sigma, bounds: tuple = None):
    """
    Returns a normal distribution.  The distribution can be truncated.

    :param x: The domain over which the distribution is cast
    :param mu: The mean
    :param sigma: The standard deviation
    :param bounds: (unspecified/None) - unbound, normal gauss_norm,
      (if specified: (lower, upper)) - binds the domain to a value or the min/max of x if bound is None.
    :return: The pdf of a (truncated) normal distribution over x
    """
    if bounds is None:
        return stat_norm(x, mu, sigma)
    left_bound, right_bound = bounds
    if left_bound is None:
        left_bound = min(x)
    if right_bound is None:
        right_bound = max(x)
    a, b = (left_bound - mu) / sigma, (right_bound - mu) / sigma
    return stat_truncnorm.pdf(x, a, b, loc=mu, scale=sigma)


def moving_average(x, y, n=3, func=lambda _y: _y, _binding=(None, None)):
    """
    Calculates the (truncated gaussian-convoluted) moving average of a dataset y over the domain x.

    :param x: The domain to which y belongs
    :param y: The function being averaged
    :param n: The standard deviation of the gaussian used to average the signal
    :param func: Applied to each element of y before averaging
      (e.g. `abs` for L-MAE, `lambda _y: np.power(_y, 2)` for L-MSE)
    :param _binding: (See 'norm' for controlling the truncation of the gaussian)
    :return: The gaussian-convoluted mean value of y at x
    """
    y = func(y)
    return np.array([np.sum(norm(x, x[i], n, _binding) * y) for i, _ in enumerate(y)])


def average_temperature(metadata_iteration_doc):
    """ Calculates the time-average temperature across scans inside the plate reader

    :param metadata_iteration_doc: Metadata from a Spark datafile
    :return: Average temperature as indicated by a metadata iteration document
    """
    try:
        temperature = metadata_iteration_doc['Temperature [°C]']
        start_times = metadata_iteration_doc['Start Time']  # SPARK_TIME_FORMAT
        end_times = metadata_iteration_doc['End Time']
    except KeyError:
        return None
    try:
        temperature = [float(t) for t in temperature]
        start_times = [datetime.strptime(_s, SPARK_TIME_FORMAT) for _s in start_times]
        end_times = [datetime.strptime(_e, SPARK_TIME_FORMAT) for _e in end_times]
    except (TypeError, ValueError):
        return None
    if len(start_times) != len(end_times):
        return None
    interval_times = [(_e - _s).total_seconds() for _s, _e in zip(start_times, end_times)]
    temperature = np.array(temperature)
    interval_times = np.array(interval_times)
    mean_temp = np.nansum(np.multiply(temperature, interval_times)) / np.nansum(interval_times)
    if np.isnan(mean_temp):
        return None
    return mean_temp


def trapz_error(f: np.ndarray, x: np.ndarray, dx=3.0):
    f = f.flatten()
    x = x.flatten()
    if f.size != x.size:
        raise ValueError("Integration requires f and x to be of the same size")
    indicies = range(0, f.size)
    integral = 0
    integral_error = 0
    for i in indicies:
        try:
            f[i + 1]
        except IndexError:
            break
        y_term = 0.5 * (f[i + 1] + f[i])
        x_term = x[i + 1] - x[i]
        y_err = np.sqrt(0.5 * (max(0.005, 0.005 * f[i + 1]) ** 2 + max(0.005, 0.005 * f[i]) ** 2))
        x_err = 0.3 * np.sqrt(2)
        e_factor = x_term / dx
        z = y_term * x_term
        z_error = np.sqrt(e_factor) * np.sqrt((x_term * y_err) ** 2 + (y_term * x_err) ** 2)
        # z distributed inside radical to avoid dividing my zero when the signal is 0
        if not np.isfinite(z_error):
            warn(f"Non-finite error term encountered at: "
                 f"(i; x_i, f_i, x_i+1, f_i+1) = "
                 f"({i}; {x[i]}, {f[i]}, {x[i + 1]}, {f[i + 1]})")
            z_error = 0
        integral += z
        integral_error += z_error ** 2
    return integral, np.sqrt(integral_error)


def power_mean(numbers):
    """ Returns the geometric mean of a set of numbers

    :param numbers: an iterable of numbers
    :return: (Product[n])^(1/N)
    """
    return np.exp(np.sum([np.log(n) for n in numbers]) / len(numbers))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class DataTamer:
    """ For processing raw signals into usable data """
    BLANK = "Blank"
    WELLPLATE = "Wellplate"
    WAVELENGTH = "Wavel."
    SIGNAL = "signal"
    UNCERTAINTY = "sigma"
    ABS = "abs"
    PL = "pl"
    LUM = "lum"
    MEASUREMENT_TYPES = [ABS, PL, LUM]

    def __init__(self, db_path=None, db_prebuilt=None):
        """ Generates a set of databases to hold wellplate and solvent information

        :param db_path: Build database from data stored at this location
        :param db_prebuilt: Use existing database
        """
        self._exp_database = dict()
        """
        'C1COCCO1': {
            ABS: {
                '200': [0.1674, 0.1091, 0.1106, ...],\n
                '100': [0.1374, 0.0741, 0.1042, ...],\n
                ...\n
                WAVELENGTHS: [200, 202, 204, 206, ...]
            },\n
            PL: {...}
        }, ...\n
        'Std96well_PP': {
            ABS: {
                P_EAV: [3.3633, 3.3795, 3.3709, 3.3120, ...],\n
                P_ESD: [0.0243, 0.0288, 0.0202, 0.0204, ...],\n
                WAVELENGTHS: [200, 202, 204, 206, ...]
            },\n
            PL: {...}
        }, ...\n
        BLANK: {
            ABS: {
                P_EAV: [0.0001, 0.0008, 0.0001, -0.0001, ...],\n
                P_ESD: [0.0002, 0.0010, 0.0006, 0.0005, ...],\n
                WAVELENGTHS: [200, 201, 202, 203, ...]
            },\n
            PL: {}
        }
        """
        self._lit_database = dict()
        """
        'O': [ [ [wavelengths], [extinctions-coeff] ] ],
        'CCCCCCCCO': [ [ [wavelengths], [extinctions-coeff] ] ],
        'CCO': [[[],[]]],
        'CS(=O)C': [[[],[]]]
        """

        self._loaded_wellplates = list()
        self._analysis_range = None

        if db_prebuilt:
            self._exp_database.update(db_prebuilt[0])
            self._lit_database.update(db_prebuilt[1])
        else:
            reload = self.load_database_recovery(db_path)
            load_lit = self.build_literature_reference_database(db_path)
            load_exp = self.build_experimental_reference_database(db_path)
            load_success = reload or (load_lit and load_exp)
            if load_success:
                with open(os.path.join(db_path, 'spark_database.json'), 'w+') as fh:
                    try:
                        json.dump(self._exp_database, fh)
                    except (TypeError, ValueError) as je:
                        print(repr(je))
                        print("Failed to save spark database")
            else:
                raise RuntimeError("Spark Reference values for data processing not found")

    def load_database_recovery(self, db_path: str):
        for dir_path, _, files in os.walk(db_path):
            for filename in files:
                file_io_str = os.path.join(dir_path, filename)
                if filename != 'spark_database.json':
                    continue
                try:
                    with open(file_io_str, 'r') as fh:
                        try:
                            self._exp_database.update(json.load(fh))
                            return True
                        except (ValueError, json.decoder.JSONDecodeError) as dle:
                            print(repr(dle))
                            return False
                except (FileNotFoundError, PermissionError) as fae:
                    print(repr(fae))
                    return False
        return False

    def build_literature_reference_database(self, db_path: str):
        for dir_path, _, files in os.walk(db_path):
            for filename in files:
                file_io_str = os.path.join(dir_path, filename)
                if filename != 'solvents.json':
                    continue
                try:
                    with open(file_io_str) as fh:
                        try:
                            self._lit_database = json.load(fh)
                            return True
                        except (ValueError, json.decoder.JSONDecodeError) as dle:
                            print(repr(dle))
                            return False
                except (FileNotFoundError, PermissionError) as fae:
                    print(repr(fae))
                    return False
        return False

    def build_experimental_reference_database(self, db_path: str):
        is_read = None
        for dir_path, _, files in os.walk(db_path):
            for filename in files:
                file_io_str = os.path.join(dir_path, filename)
                if '.json' in filename:
                    continue
                if '.xlsx' not in filename:
                    continue
                if '~$' in filename:  # temp files made by excel
                    continue
                # Load the workbook
                try:
                    workbook = pyxl.load_workbook(file_io_str)
                except Exception as e:
                    print(repr(e))
                    continue

                # Cycle through each sheet (a measurement type: abs, pl, lum)
                for sheet_id in self.MEASUREMENT_TYPES:
                    worksheet = workbook[sheet_id]
                    wavelengths = list()
                    classification = None
                    # Read row by row
                    for row in worksheet.iter_rows():
                        active_row = [cell.value for cell in row if cell.value is not None]
                        active_length = len(active_row)
                        # Skip emtpy rows
                        if active_length == 0:
                            continue
                        # Single-element rows are headers for the measurement class (blank, solvent, wellplate)
                        elif active_length == 1:
                            if active_row[0] == self.WELLPLATE:
                                classification = self.WELLPLATE
                            else:
                                classification = None
                            continue
                        # Multi-element rows are data
                        else:
                            d_header, *d_values = active_row
                            d_values = [float(v) if not isinstance(v, str) else None for v in d_values]
                            if d_header == self.WAVELENGTH:
                                wavelengths = d_values
                            else:
                                if not wavelengths:
                                    print(f"Warning in file '{filename}': "
                                          f"File nor properly formatted, expected wavelengths entry "
                                          f"before data signal entry")
                                    is_read = False
                                    continue
                                if not (len(wavelengths) == len(d_values)):
                                    print(f"Warning in file '{filename}': "
                                          f"Length of data signal ({len(d_values)}) does not match "
                                          f"length of associated wavelengths ({len(wavelengths)}) ")
                                    is_read = False
                                    continue
                                d_entry, d_type = d_header.rsplit("_", 1)
                                self._exp_database.setdefault(d_entry, dict())
                                self._exp_database[d_entry].setdefault(sheet_id, dict())
                                self._exp_database[d_entry][sheet_id][d_type] = d_values
                                self._exp_database[d_entry][sheet_id][self.WAVELENGTH] = wavelengths
                                is_read = True if is_read is None else is_read
                                if classification == self.WELLPLATE:
                                    if d_entry not in self._loaded_wellplates:
                                        self._loaded_wellplates.append(d_entry)
        return is_read

    # # # #

    @staticmethod
    def apply_floor(data_array, floor):
        """ Replaces all data less than floor with floor

        compare: floor_signal() which is the high-level command that actually commits the update to the project

        :param data_array: An iterable of data
        :param floor: The floor to be used
        :return: [floored data]
        """
        def compare_to_floor(x, _floor):
            if x is None:
                return True
            elif np.isnan(x):
                return True
            else:
                return x > _floor

        return [v if compare_to_floor(v, floor) else floor for v in data_array]

    def convert_metadata(self, _project: sp.SparkProjectFileHandler, i_iter):
        """ Translates metadata formats for a plate key and scan mode into the form used by the DataTamer

        :param _project: A project file handler
        :param i_iter: The iteration number
        :return: the wellplate key, the scan mode (either "abs" or "pl")
        """
        _well_plate_key = None
        _scan_mode = None
        _plate_id, _mode = _project.get_spark_metadata(i_iter, key=('Plate', 'Mode'))
        for k in self._loaded_wellplates:
            if k in _plate_id:
                _well_plate_key = k
                break
        _scan_mode = {
            'absorbance': 'abs',
            'fluorescence': 'pl',
        }.get(_mode.lower(), None)
        return _well_plate_key, _scan_mode

    @staticmethod
    def background_correction(waves, ints, solvent_waves, solvent_ints):
        """ Removes the background signal of some solvent

        :param waves: data x-axis
        :param ints: data y-axis
        :param solvent_waves: reference x-axis
        :param solvent_ints: reference y-axis
        :return: the signal of the reference as present in the data
        """
        interpl_values = (dynamic_scipy.interp1d(solvent_waves, solvent_ints))(waves)
        peaks, _ = scipy.signal.find_peaks(interpl_values)
        prominence = scipy.signal.peak_prominences(interpl_values, peaks)
        peak_prominence_index = peaks[np.where(prominence[0] == np.amax(prominence[0]))[0]]
        solvent_multiplier = ints[peak_prominence_index] / interpl_values[peak_prominence_index]
        solvent_background = solvent_multiplier * interpl_values
        return solvent_background, solvent_multiplier

    @staticmethod
    def is_number(_str):
        """ Helper method, checks if a string is numerical

        :param _str: A string to test
        :return: True - can be cast to a float, False - otherwise
        """
        try:
            float(_str)
            return True
        except (ValueError, TypeError):
            return False

    @staticmethod
    def consolidate_solvent_data(solvent_data) -> List[Tuple[float, str]]:
        """ Takes solvent well data and consolidates duplicate entries

        It is possible for a well context to say [["DMSO", 90], ["DMSO", 90]] instead of [["DMSO", 180]]

        :param solvent_data: Well context data for solvents
        :return: A consolidated version
        """
        memory = dict()
        for solvent_datum in solvent_data:
            _, volume, smiles = solvent_datum
            memory.setdefault(smiles, 0)
            memory[smiles] += volume
        return [[v, k] for k, v in memory.items()]  # noqa (want to hint list structure but list can't)

    # # # #

    def project_iterator(self, project: sp.SparkProjectFileHandler, method, library,
                         on_wells=None, on_iterations=None, floor=None, options=None) -> Generator:
        """
        Helper method: many processes are looped over all wells and iterations, so this allows that to be automated

        **method arguments**
          * self - A copy of the :class:`DataTaner` object making the call
          * project - The :class:`sp.SparkProjectFileHandler` object in use
          * library - Either the 'self._exp_database' or 'self._lit_database' properties <dict>
          * iteration_data - The iteration data element from project.get_data-by_well
          * x - The x-axis data
          * measure_mode - The 'mode' detail from project.spark_method_details()
          * well_plate_id - The wellplate type being processed
          * well_id - The ID of the well processed
          * i - The iteration #
          * options - A dictionary of additional kwargs that a method may need

        :param project: A copy of the project (a SparkProjectFileHandler)
        :param method: A function handle expecting arguments in the order and type described above
        :param library: Either the 'self._exp_database' or 'self._lit_database' properties
        :param on_wells: A list of wells over which the iterator should iterate
        :param on_iterations: A list of iterations <int> over which the iterator should iterate
        :param floor: Replaces all values less than floor with floor (None - do not floor)
        :param options: A dictionary of any additional arguments required
        :return: Generator (well_ID, iteration #, datum) for updating the project
        """
        if on_wells is None:
            on_wells = []
        if on_iterations is None:
            on_iterations = []

        measure_mode = project.get_spark_method_details().get('mode', None)
        if measure_mode is None:
            raise ValueError("Measurement modes is not specified")

        for well_id in project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue
            data, _ = project.get_data_by_well(well_id)
            for i, d, _, x in data:
                if on_iterations and (int(i) not in on_iterations):
                    continue
                well_plate_type, check_mode = self.convert_metadata(project, i)

                # Quickly check that things are compatible and that the analysis can be performed
                if check_mode != measure_mode:
                    raise ValueError(f"Measurement modes ({measure_mode}) in metadata "
                                     f"and method_specifications ({check_mode}) do not agree")

                try:
                    method(self, project, library, d, x, measure_mode, well_plate_type, well_id, i, options)
                except SkipIteration:
                    continue

                yield well_id, i, d

                if floor is not None:
                    d['y-axis (processed)'] = self.apply_floor(d['y-axis (processed)'], floor)

    def subtract_wellplate_background(self, project: sp.SparkProjectFileHandler, *,
                                      on_wells=None, on_iterations=None, floor=None):
        """
        Subtracts the contribution of the wellplate from the signal.

        References: :class:`sp.SparkProjectFileHandler`

        :param project: A SparkProjectFileHandler object
        :param on_wells: A list of wells which should be processed
        :param on_iterations: A list of iterations (int) which should be processed
        :param floor: If not None: Replaces all values less than floor with floor
        :return: None
        """

        #                     self,  project,  library,  d,  x,  measure_mode,  well_plate_type, well_id, i, options
        def _internal_method(_self, _project, _library, _d, _x, _measure_mode, _well_plate_type, *_, **__):
            if _well_plate_type is None:
                raise ValueError("Wellplate type not identified")

            ref_val = np.array(_library[_well_plate_type][_measure_mode][_self.SIGNAL])
            ref_waves = np.array(_library[_well_plate_type][_measure_mode][_self.WAVELENGTH])
            exp_val = np.array(_d.get('y-axis (processed)', _d['y-axis']), dtype=np.float)
            exp_waves = np.array(_x)
            bkg_val = (dynamic_scipy.interp1d(ref_waves, ref_val))(exp_waves)
            # print("EXP")
            # pprint(exp_val.flatten())
            # print("BKG")
            # pprint(bkg_val.flatten())
            exp_val = exp_val.flatten() - bkg_val.flatten()
            _d['y-axis (processed)'] = exp_val.tolist()

        project.commit_data_by_generator(
            self.project_iterator(project, _internal_method, self._exp_database, on_wells, on_iterations, floor),
            log=f"subtract_wellplate_background(on_wells={on_wells}, on_iterations={on_iterations}, floor={floor})"
        )

    def subtract_signal_baseline(self, project: sp.SparkProjectFileHandler, *,
                                 on_wells=None, on_iterations=None, floor=None):
        """
        Subtracts the contribution of the Spark's bias from the signal.

        References: :class:`sp.SparkProjectFileHandler`

        :param project: A SparkProjectFileHandler object
        :param on_wells: A list of wells which should be processed
        :param on_iterations: A list of iterations (int) which should be processed
        :param floor: If not None: Replaces all values less than floor with floor
        :return: None
        """

        #                     self,  project,  library,  d,  x,  measure_mode, well_plate_type, well_id, i, options
        def _internal_method(_self, _project, _library, _d, _x, _measure_mode, *_, **__):
            ref_val = np.array(_library[self.BLANK][_measure_mode][_self.SIGNAL])
            ref_waves = np.array(_library[self.BLANK][_measure_mode][_self.WAVELENGTH])
            ref_val = moving_average(ref_waves, ref_val, n=10)
            exp_val = np.array(_d.get('y-axis (processed)', _d['y-axis']))
            exp_waves = np.array(_x)
            bkg_val = (dynamic_scipy.interp1d(ref_waves, ref_val))(exp_waves)
            exp_val = exp_val.flatten() - bkg_val.flatten()
            _d['y-axis (processed)'] = exp_val.tolist()

        project.commit_data_by_generator(
            self.project_iterator(project, _internal_method, self._exp_database, on_wells, on_iterations, floor),
            log=f"subtract_signal_baseline(on_wells={on_wells}, on_iterations={on_iterations}, floor={floor})"
        )

    def subtract_background_by_solvent_volume(self, project: sp.SparkProjectFileHandler, *,
                                              on_wells=None, on_iterations=None, floor=None):
        """
        Subtracts the contribution of the solvent from the signal by using the volume well context details

        References: :class:`sp.SparkProjectFileHandler`

        :param project: A SparkProjectFileHandler object
        :param on_wells: A list of wells which should be processed
        :param on_iterations: A list of iterations (int) which should be processed
        :param floor: If not None: Replaces all values less than floor with floor
        :return: None
        """

        #                     self,  project,  library,  d,  x, measure_mode, well_plate_type, well_id, i, options
        def _internal_method(_self, _project, _library, _d, _x, _measure_mode, *_, **__):
            solvent_data = _self.consolidate_solvent_data(_d['well_context'].get('solvents', []))
            for solvent_datum in solvent_data:
                volume_ul, sol_smiles = solvent_datum
                sol_background_data = _library.get(sol_smiles, {}).get(_measure_mode, None)
                if sol_background_data is None:
                    print(f"Missing solvent data for '{sol_smiles}'")
                    continue
                volumes = [float(vol) for vol in sol_background_data.keys() if _self.is_number(vol)]
                volumes += [0, ]
                volumes = np.array(volumes)

                wavelengths = sol_background_data[_self.WAVELENGTH]
                wavelengths = np.array(wavelengths)
                vx, wx = np.meshgrid(volumes, wavelengths)
                ax = [v for k, v in sol_background_data.items() if _self.is_number(k)]
                ax += [[0, ] * len(ax[0])]
                ax = np.array(ax)
                tck, _, ier, msg = dynamic_scipy.bisplrep(vx.flatten(), wx.flatten(), ax.flatten(),
                                                          full_output=True)
                if ier > 0:
                    print("scipy.interpolate.bisplrep error: ", msg)
                    continue
                exp_values = np.array(_d.get('y-axis (processed)', _d['y-axis']))
                exp_waves = np.array(_x)
                volume_ul = np.array([volume_ul, ])
                predicted_solvent_spectrum = scipy.interpolate.bisplev(volume_ul, exp_waves, tck)
                exp_values = exp_values.flatten() - predicted_solvent_spectrum.flatten()
                _d['y-axis (processed)'] = exp_values.tolist()

        project.commit_data_by_generator(
            self.project_iterator(project, _internal_method, self._exp_database, on_wells, on_iterations, floor),
            log=f"subtract_background_by_solvent_volume("
                f"on_wells={on_wells}, on_iterations={on_iterations}, floor={floor})"
        )

    def subtract_background_by_solvent_peaks(self, project: sp.SparkProjectFileHandler, *,
                                             on_wells=None, on_iterations=None, floor=None):
        """
        Subtracts the contribution of the solvent from the signal by using solvent peak fitting

        References: :class:`sp.SparkProjectFileHandler`

        :param project: A SparkProjectFileHandler object
        :param on_wells: A list of wells which should be processed
        :param on_iterations: A list of iterations (int) which should be processed
        :param floor: If not None: Replaces all values less than floor with floor
        :return: None
        """

        #                     self,  project,  library,  d,  x, measure_mode, well_plate_type, well_id, i, options
        def _internal_method(_self, _project, _library, _d, _x, *_, **__):
            solvent_data = _self.consolidate_solvent_data(_d['well_context'].get('solvents', []))
            for solvent_datum in solvent_data:
                volume_ul, sol_smiles = solvent_datum
                sol_background_data = _library.get(sol_smiles, None)
                if sol_background_data is None:
                    print(f"Missing solvent data for '{sol_smiles}'")
                    continue
                waves, epsilon = sol_background_data[0][0]
                waves = np.array(waves)
                epsilon = np.array(epsilon)
                exp_values = np.array(_d.get('y-axis (processed)', _d['y-axis']))
                exp_waves = np.array(_x)
                sol_bkg, _ = _self.background_correction(exp_waves, exp_values, waves, epsilon)
                exp_values = exp_values.flatten() - sol_bkg.flatten()
                _d['y-axis (processed)'] = exp_values.tolist()

        project.commit_data_by_generator(
            self.project_iterator(project, _internal_method, self._lit_database, on_wells, on_iterations, floor),
            log=f"subtract_background_by_solvent_peaks("
                f"on_wells={on_wells}, on_iterations={on_iterations}, floor={floor})"
        )

    def subtract_background_by_reference_well(self, project: sp.SparkProjectFileHandler, *,
                                              on_wells=None, on_iterations=None, floor=None,
                                              ref_opts: dict = None):
        """
        Subtracts a reference signal from the processed wells.

        References: :class:`sp.SparkProjectFileHandler`

        **ref_opts**
          * ref_well - The well ID used as the reference (":" can be used to match the corresponding field of the
            well being processed, "::" will reference all wells against themselves)
          * iter_mode - 'static' if the reference is always iteration 0 or 'dynamic' if the iteration should match
            that of the iteration being processed
          * src_mode - 'y-axis' to use the raw signal or 'y-axis (processed)' to use the processed signal as reference
            (note 'y-axis (processed)' can now be abbreviated 'processed')

        :param project: A SparkProjectFileHandler object
        :param on_wells: A list of wells which should be processed
        :param on_iterations: A list of iterations (int) which should be processed
        :param floor: If not None: Replaces all values less than floor with floor
        :param ref_opts: A dictionary for the options of how the reference well is to be treated
        :return: None
        """
        if ref_opts is None:
            ref_opts = dict()

        #                     self,  project,  library,  d,  x,  measure_mode, well_plate_type,  well_id,  i,  options
        def _internal_method(_self, _project, _library, _d, _x, _measure_mode, _, _well_id, _i, _options):
            reference_well_id = _options.get('ref_well', 'not specified')
            iteration_mode = _options.get('iter_mode', 'static')
            source_mode = _options.get('src_mode', 'y-axis')

            if reference_well_id == "::":
                reference_well_id = _well_id
            elif reference_well_id in _WELLPLATE_COLUMNS:
                sel_col = reference_well_id[1:]
                sel_row = _well_id[:1]
                reference_well_id = f"{sel_row}{sel_col}"
            elif reference_well_id in _WELLPLATE_ROWS:
                sel_row = reference_well_id[:1]
                sel_col = _well_id[1:]
                reference_well_id = f"{sel_row}{sel_col}"
            elif reference_well_id in WELL_IDS:
                pass
            else:
                raise ValueError(f"kwarg 'ref_well' must be a well ID not '{reference_well_id}'")

            if iteration_mode == 'static':
                iter_index = 0
            elif iteration_mode == 'dynamic':
                iter_index = _i
            else:
                raise ValueError(f"kwarg 'iter_mode' must be 'static' or 'dynamic' not '{iteration_mode}'")

            ref_d, _, _, x_data = project.get_datum_by_well(reference_well_id, iter_index)

            if source_mode == 'y-axis':
                ref_values = np.array(ref_d['y-axis'])
            elif 'processed' in source_mode:
                ref_values = np.array(ref_d.get('y-axis (processed)', ref_d['y-axis']))
                if (reference_well_id == _well_id) and (iter_index == _i):
                    # if using processed signal, we need to prevent a well from zeroing itself (as this will affect
                    # all subsequent attempts to reference this well).
                    raise SkipIteration
            else:
                raise ValueError(f"kwarg 'src_mode' must be 'y-axis' or 'y-axis (processed)' not '{iteration_mode}'")

            ref_waves = np.array(x_data)

            exp_values = np.array(_d.get('y-axis (processed)', _d['y-axis']))
            exp_waves = np.array(_x)
            if ref_waves == exp_waves:
                exp_values = exp_values - ref_values
            else:
                sol_bkg, _ = _self.background_correction(exp_waves, exp_values, ref_waves, ref_values)
                exp_values = exp_values.flatten() - sol_bkg.flatten()
            _d['y-axis (processed)'] = exp_values.tolist()

        project.commit_data_by_generator(
            self.project_iterator(project, _internal_method, self._lit_database,
                                  on_wells, on_iterations, floor, ref_opts),
            log=f"subtract_background_by_reference_well("
                f"on_wells={on_wells}, on_iterations={on_iterations}, floor={floor}, ref_opts={ref_opts})"
        )

    def estimate_signal_uncertainty(self, measure_mode, well_iteration_data):
        exp_values = np.array(well_iteration_data.get('y-axis (processed)', well_iteration_data['y-axis']))
        exp_waves = np.array(well_iteration_data['x-axis'])  # This will be a pointer not the true x-axis

        ref_val = np.array(self._exp_database[self.BLANK][measure_mode][self.UNCERTAINTY])
        ref_waves = np.array(self._exp_database[self.BLANK][measure_mode][self.WAVELENGTH])
        ref_val = moving_average(ref_waves, ref_val, n=10)
        ref_val = (dynamic_scipy.interp1d(ref_waves, ref_val))(exp_waves)

        return np.sqrt(np.power(0.005 * exp_values, 2) + np.power(ref_val, 2) + np.power(0.005, 2))

    def floor_signal(self, project: sp.SparkProjectFileHandler, *,
                     on_wells=None, on_iterations=None, floor=None):
        """
        Subtracts the contribution of the wellplate from the signal.

        References: :class:`sp.SparkProjectFileHandler`

        compare: apply_floor() which actually does the math

        :param project: A SparkProjectFileHandler object
        :param on_wells: A list of wells which should be processed
        :param on_iterations: A list of iterations (int) which should be processed
        :param floor: If not None: Replaces all values less than floor with floor
        :return: None
        """

        def _internal_method(*_, **__):
            pass

        project.commit_data_by_generator(
            self.project_iterator(project, _internal_method, self._exp_database, on_wells, on_iterations, floor),
            log=f"floor_signal(on_wells={on_wells}, on_iterations={on_iterations}, floor={floor})"
        )


class Peaks:
    """ A Data processor for static Peak-related experiments """
    def __init__(self, project: sp.SparkProjectFileHandler, tamer: DataTamer):
        """ Create a processor for Peak data processing

        :param project: A project file
        :param tamer: A data tamer for processing signals
        """
        self.project = project
        self.tamer = tamer

    @staticmethod
    def _process_peaks(wavelength_axis, signal_axis, sort_by: str, **find_peaks_args) -> list:
        """ Calculates the peaks and their prominences for a given data set

        :param wavelength_axis: The x-axis (nm)
        :param signal_axis: the y-data
        :param sort_by: How the results should be sorted ('prominence' - sort most prominent to least prominent,
          'wavelength' - sort red to blue, 'product' - sort by prominence * wavelength, 'none' - use default ordering)
        :param find_peaks_args: Args passed to scipy.signal.find_peaks() method
        :return: A list of [wavelength, prominence] entries
        """
        _prominence = 'prominence'
        _wavelength = 'wavelength'
        _product = 'product'
        _none = 'none'
        sort_by = sort_by.lower()
        if sort_by not in [_prominence, _wavelength, _product, _none]:
            print(f"Warning: sort_by '{sort_by}' not recognized")
            sort_by = _none

        if isinstance(signal_axis, list):
            signal_axis = np.array(signal_axis)
        if isinstance(wavelength_axis, list):
            wavelength_axis = np.array(wavelength_axis)
        mask = ~np.isnan(signal_axis.astype(float))
        peak_indices, *_ = scipy.signal.find_peaks(signal_axis[mask], **find_peaks_args)
        peak_lambda: np.ndarray = wavelength_axis[mask][peak_indices]
        prominence, *_ = scipy.signal.peak_prominences(signal_axis[mask], peak_indices)
        result = np.array([peak_lambda, prominence])
        if sort_by == _prominence:
            result = result[:, (-result[1]).argsort()]
        elif sort_by == _wavelength:
            result = result[:, (-result[0]).argsort()]
        elif sort_by == _product:
            result = result[:, (-result[0] * result[1]).argsort()]
        else:
            pass
        ret_val = result.tolist()
        if isinstance(ret_val, list):
            return ret_val
        else:
            return [ret_val, ]

    def process_peaks(self, on_wells, iteration: int, sort_by: str = 'none', **find_peaks_args):
        """ Determines the peak data for a project on the given wells for the specified iteration

        Wraps private version which is called on each well.

        :param on_wells: Which wells to process (Falsey - process All wells)
        :param iteration: Which iteration to use
        :param sort_by: How peaks should be sorted ('prominence' - sort most prominent to least prominent,
          'wavelength' - sort red to blue, 'product' - sort by prominence * wavelength, 'none' - use default ordering)
        :param find_peaks_args: Additional arguments for the find_peaks() method
        :return: None
        """
        properties = {
            well_id: {
                'wavelengths': None,
                'prominences': None,
            } for well_id in self.project.well_ids()
        }

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue
            data, meta, method_spec, x_data = self.project.get_datum_by_well(well_id, iteration)

            result_lambda_prom = self._process_peaks(x_data, data.get('y-axis (processed)', data['y-axis']),
                                                     sort_by, **find_peaks_args)
            properties[well_id]['wavelengths'] = result_lambda_prom[0]
            properties[well_id]['prominences'] = result_lambda_prom[1]

        self.project.commit_property("peaks", [iteration, ], properties,
                                     log=f"process_peaks("
                                         f"on_wells={on_wells}, iteration={iteration}, sort_by={sort_by}, "
                                         f"find_peaks_args={find_peaks_args})")

    def commit_to_database(self, name='SP', n_retry=4, wait=1, iteration=0, on_wells=None, *,
                           peaks=False, spectra=True, campaign=None):
        """ Uploads processed data to the database

        :param name: Logon name (default: 'SP')
        :param n_retry: Number of times a commit will be retried
        :param wait: Time between attempts (in seconds)
        :param iteration: The iteration of the data being committed
        :param on_wells: The wells being committed to the database (falsey - all wells)
        :param peaks: Include peak prominence data
        :param spectra: Include spectral data
        :param campaign: The name of the associated campaign
        :return: success_log, failure_log
        """
        if on_wells is None:
            on_wells = []
        commit_log = list()
        iteration = str(int(iteration))

        project = self.project.pull_project_copy()
        measurement_mode = project['method_specifications']['mode']
        data_origin = ('measured', 'amd_platform')

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue

            # Preload general/contextual information
            iteration_data = project.get(well_id, {}).get(iteration, {})
            well_context = iteration_data.get('well_context', {})
            try:
                smiles = ".".join(
                    well_context.get(
                        'target_molecules',  # The smiles for the target analytes are in this list
                        [well_context['reagents'][0][2], ]  # Backup: Assume first reagent is the target
                    )
                )
            except (KeyError, IndexError):
                commit_log.append((well_id, None, "target_analytes field failed lookup"))
                continue
            try:
                solvent = well_context['solvents']
            except KeyError:
                commit_log.append((well_id, None, "solvents field failed lookup"))
                solvent = "Solvent not reported"
            meta_data = {'well_context': well_context,
                         'history': project.get('history', ""),
                         'source': self.project.directory}
            try:
                temperature = [average_temperature(project['abscissa']['metadata'][iteration]), "C"]
            except (KeyError, TypeError):
                temperature = ["Temperature not reported", "C"]

            # If peak data requested
            if peaks:
                data_name = "_uv-vis_spectral_peaks"
                type_tag = f"exp_{measurement_mode}_stick"

                try:
                    data = project['calculated_properties']['peaks'][iteration][well_id]
                    value = {
                        'xdatalabel': ['wavelength', 'nm'],
                        'xdata': data['wavelengths'],
                        'ydatalabel': ['prominence', 'a.u.'],
                        'ydata': [data['prominences'], ],
                        'temperature': temperature,
                    }
                except KeyError as ke:
                    commit_log.append((well_id, type_tag,
                                       f"Failed to pull data: key error on '{next(iter(ke.args), '<unknown>')}'"))
                    continue

                for _ in range(n_retry):
                    code, _, _ = dri.add_data(name, smiles=smiles, data_name=data_name, campaign_name=campaign,
                                              type_tag=type_tag, data_origin=data_origin,
                                              value=value, solvent=solvent, meta_data=meta_data,
                                              get_details=True)
                    if code == "Success":
                        commit_log.append((well_id, type_tag, False))
                        break
                    elif code == "Connection Failed":
                        commit_log.append((well_id, type_tag, "Database not online"))
                    else:
                        time.sleep(wait)
                        continue
                else:
                    commit_log.append((well_id, type_tag, f"Max retries reached for database push"))

            # If spectral data requested
            if spectra:
                data_name = "_uv-vis_spectrum"
                type_tag = f"exp_{measurement_mode}_spec"

                x_axis_pointer = iteration_data.get('x-axis', None)
                try:
                    y_axis = iteration_data['y-axis (processed)']
                    meta_data['signal'] = "y-axis (processed)"
                except KeyError:
                    try:
                        y_axis = iteration_data['y-axis']
                        meta_data['signal'] = "y-axis"
                    except KeyError:
                        commit_log.append((well_id, type_tag, f"No spectral information found"))
                        continue
                try:
                    x_axis = project['abscissa']['x-axis'][x_axis_pointer]
                except KeyError:
                    commit_log.append((well_id, type_tag, f"Failed to find x-axis data"))
                    continue
                if (not y_axis) or (not x_axis):
                    commit_log.append((well_id, type_tag, f"Axis data is bad"))
                    continue

                value = {
                    'xdatalabel': ['wavelength', 'nm'],
                    'xdata': x_axis,
                    'ydatalabel': ['absorbance', 'a.u.'],
                    'ydata': [y_axis, ],
                    'temperature': temperature,
                }

                for _ in range(n_retry):
                    code, _, _ = dri.add_data(name, smiles=smiles, data_name=data_name, campaign_name=campaign,
                                              type_tag=type_tag, data_origin=data_origin,
                                              value=value, solvent=solvent, meta_data=meta_data)
                    if code == "Success":
                        commit_log.append((well_id, type_tag, False))
                        break
                    elif code == "Connection Failed":
                        commit_log.append((well_id, type_tag, "Database not online"))
                    else:
                        time.sleep(wait)
                        continue
                else:
                    commit_log.append((well_id, type_tag, f"Max retries reached for database push"))

        success_log = [item for item in commit_log if not item[2]]
        failure_log = [item for item in commit_log if item[2]]
        return success_log, failure_log


class SignalCheck(Peaks):
    """ Extension of Peaks to add signal-checking """
    def __init__(self, project: sp.SparkProjectFileHandler, tamer: DataTamer, **method_kwargs):
        """ Create a processor for Peak data processing

        :param project: A project file
        :param tamer: A data tamer for processing signals
        :keyword lb: The lower bound (default: 0.5)
        :keyword up: The upper bound (default: 2)
        :keyword wavelengths: Defines range which is checked (default: [200, 1000] nm)
        :keyword target: The desired signal value (default: geometric mean of the upper and lower bounds)
        """
        super(SignalCheck, self).__init__(project, tamer)
        self.lower_bound = method_kwargs.get('lb', 0.5)
        self.upper_bound = method_kwargs.get('ub', 2)
        self.wavelength_range = method_kwargs.get('wavelengths', [200, 1000])
        self.target = method_kwargs.get('target', power_mean([self.lower_bound, self.upper_bound]))

    def in_range(self, v):
        """ Helper method, checks if the signal is within the

        :param v: A number to be checked against self.wavelength_range --> Lower bound & Upper bound
        :return: Lower bound < v < Upper bound
        """
        return self.wavelength_range[0] < v < self.wavelength_range[1]

    def check_bounds(self, on_wells):
        """ Wrapper to check peaks in multiple wells for signals being "in range"

        :param on_wells: A list of wells to restrict processing to (falsey - use all wells)
        :return: The wells which were out of range
        """
        bad_wells = dict()
        self.process_peaks(on_wells, 0, sort_by='prominence', prominence=self.lower_bound / 100)
        data = self.project.get_property('peaks')
        data = data.get("0", dict())
        for well_id, well_data in data:
            if on_wells and (well_id not in on_wells):
                continue
            prominences = well_data.get('prominences', [])
            prom_waves = well_data.get('wavelengths', [])
            # Reduce list of prominences to only those in the inspected range
            prominences = [prominences[i] for i, v in enumerate(prom_waves) if self.in_range(v)]
            # If none remain, assume signal was very weak
            if not prominences:
                bad_wells[well_id] = self.lower_bound / 100
                continue
            # If there are signals, check them against the qualifying range
            #  The max signal is used as the representative
            elected_prominence = max(prominences)
            if not (self.lower_bound < elected_prominence < self.upper_bound):
                bad_wells[well_id] = elected_prominence
        return bad_wells


class Kinetics:
    """ A Data processor for dynamic experiments """
    def __init__(self, project: sp.SparkProjectFileHandler, tamer: DataTamer):
        """ Create a processor for Kinetic data processing

        :param project: A project file
        :param tamer: A data tamer for processing signals
        """
        self.project = project
        self.tamer = tamer
        self.times: List[timedelta] = list()

    @property
    def times_log_entry(self):
        """ Helper method for internal logging purposes """
        return {'Interval times (min) set to': [t.total_seconds() / 60.0 for t in self.times]}

    def use_given_interval(self, interval_minutes=0):
        """
        Used to set intervals for the time between scans.

        When the data is cattywampus and you just need to set the interval time manually

        :return: A list of the intervals in seconds
        """
        default_interval = timedelta(minutes=interval_minutes)
        project_metadata = self.project.get_spark_metadata(_all=True)
        iterations = sorted([int(ii) for ii in project_metadata.keys()])

        intervals = list()
        for _ in iterations:
            intervals.append(default_interval)
        intervals[0] = timedelta(seconds=0)

        self.times = intervals

    def use_interval_times(self):
        """ Used to set intervals for the time between scans.

        The default value returned will be the interval specified by the method, but if start and end times for two
        adjacent runs are present, the difference between their average will be given.  When retry attempts are
        involved, the average is taken from the first start time and last end time.

        :return: A list of the intervals in seconds

        :raises ValueError: If there are fewer than 2 iterations or there are missing iterations
        """
        project_metadata = self.project.get_spark_metadata(_all=True)
        method_data = self.project.get_spark_method_details()

        default_interval = method_data.get('interval', 0)
        if (default_interval <= 0) or (default_interval is None):
            default_interval = timedelta(minutes=0)
        else:
            default_interval = timedelta(minutes=default_interval)

        iterations = sorted([int(ii) for ii in project_metadata.keys()])
        if len(iterations) < 2:
            raise ValueError(f"Interval calculation requires at least two entries")

        check_list = list(range(min(iterations), max(iterations) + 1))
        if iterations != check_list:
            raise ValueError(f"Missing interval key(s), check: {check_list} vs {min(iterations)}:{max(iterations)}")

        intervals = list()
        prev_time = REF_TIME
        for ii in iterations:
            try:
                start_time = project_metadata[str(ii)]['Start Time'][0]
                end_time = project_metadata[str(ii)]['End Time'][-1]
                _s = datetime.strptime(start_time, SPARK_TIME_FORMAT) - REF_TIME
                _e = datetime.strptime(end_time, SPARK_TIME_FORMAT) - REF_TIME
                mean_time = REF_TIME + (_s + _e) / 2.0
            except (TypeError, KeyError, IndexError, ValueError):
                current_time = None
            else:
                current_time = mean_time
            if prev_time is None or current_time is None:
                intervals.append(default_interval)
            else:
                intervals.append(current_time - prev_time)
            prev_time = current_time
        intervals[0] = timedelta(seconds=0)

        self.times = intervals

    def use_exposure_times(self, q_name, assay):
        """ Attempts to determine the actual time under exposure from the log files

        :param q_name: The associated queue
        :param assay: The associated assay
        :return: None, but sets self.times to a List[timedelta]
        """
        project_metadata = self.project.get_spark_metadata(_all=True)
        method_data = self.project.get_spark_method_details()

        def get_exposure_times(_q_name, _assay):
            manifest_document = list()
            # Collate all log files (a project may span multiple steps)
            for dir_path, _, files in os.walk(self.project.directory):
                for filename in files:
                    if not (("Spark_job_" in filename)
                            and (".log" in filename)
                            and (_q_name in filename)
                            and (_assay in filename)):
                        continue
                    spark_log_file = os.path.join(dir_path, filename)
                    try:
                        with open(spark_log_file, 'r') as fh:
                            # Extract the lines about Photoreactor stuff
                            for line in fh:
                                try:
                                    _manifest_time, _manifest_msg = line.strip().split(' - ')
                                    _manifest_time = datetime.strptime(_manifest_time, MANIFEST_TIME_FORMAT)
                                except (ValueError, TypeError):
                                    continue
                                if "Photoreactor" in _manifest_msg:
                                    manifest_document.append((_manifest_time, _manifest_msg))
                    except (FileNotFoundError, PermissionError):
                        continue
            # Sort in chronological order
            manifest_document.sort(key=lambda x: x[0])

            # An exposure times is ("Unload event" minus "most recent Load event")
            in_rxtr = False
            durations = list()
            memory = None
            for manifest_time, manifest_msg in manifest_document:
                if not in_rxtr:
                    if 'Photoreactor Load' in manifest_msg:
                        memory = manifest_time
                else:
                    if 'Photoreactor Unload' in manifest_msg:
                        time_stamp = manifest_time
                        if memory is not None:
                            durations.append(time_stamp - memory)
                            memory = None
            return durations

        default_exposure = method_data.get('interval', 0)
        if (default_exposure <= 0) or (default_exposure is None):
            default_exposure = timedelta(minutes=0)
        else:
            default_exposure = timedelta(minutes=default_exposure)

        iterations = sorted([int(ii) for ii in project_metadata.keys()])
        if len(iterations) < 2:
            raise ValueError(f"Interval calculation requires at least two entries")

        check_list = list(range(min(iterations), max(iterations) + 1))
        if iterations != check_list:
            raise ValueError(f"Missing interval key(s), check: {check_list} vs {min(iterations)}:{max(iterations)}")

        default_times = [default_exposure, ] * len(iterations)
        pulled_times = get_exposure_times(q_name, assay)

        if (not pulled_times) or (len(pulled_times) != len(default_times)):
            self.times = default_times
        else:
            self.times = [timedelta(seconds=0), ] + pulled_times

    @staticmethod
    def _simple_regression(well_data_table: np.ndarray, start_from):
        t_axis: np.ndarray = well_data_table[1:, 0]  # [1, t]
        x_axis: np.ndarray = well_data_table[0, 1:]  # [1, l]
        y_axis: np.ndarray = well_data_table[1:, 1:]  # [t, l]

        # start_from = 350  # lower lambda in nm
        # convert to index
        for index, wavelength in enumerate(x_axis):
            if wavelength >= start_from:
                start_from = index
                break
        else:
            start_from = 0

        try:
            initial_signal = np.nanmean(y_axis[0, :], axis=0)
            signal_max_lambda = np.nanargmax(initial_signal[start_from:]) + start_from
            # max lambda at t=0, [1, ]  ## or at least t_min for full set provided it isn't too big
        except ValueError as ve:
            print("")
            raise ve

        if signal_max_lambda == y_axis.shape[1]:
            print("Warning: Signal maximum is on the high-wavelength edge of the scan")
            raise ValueError("Signal argmax exceeds max wavelength")

        signal_max_absorbance = initial_signal[signal_max_lambda]  # max abs at t=0, [1, ]
        y_axis = np.log(y_axis[:, signal_max_lambda] / signal_max_absorbance)  # [t, 1]

        mean_signal = np.nanmean(y_axis)  # [1, ]
        mean_time = np.nanmean(np.multiply(t_axis, np.divide(y_axis, y_axis)))

        delta_s = y_axis - mean_signal
        delta_t = t_axis - mean_time

        beta = np.nansum(np.multiply(delta_t, delta_s)) / np.nansum(np.multiply(delta_t, delta_t))
        alpha = mean_signal - beta * mean_time

        epsilon = y_axis - (beta * t_axis + alpha)
        epsilon = np.multiply(epsilon, epsilon)
        if len(t_axis) > 2:
            var_time = np.nansum(np.power(delta_t, 2))
            sigma_b = np.sqrt(
                np.divide(
                    np.multiply(
                        1 / (len(t_axis) - 2.0),
                        np.nansum(epsilon)
                    ),
                    var_time
                )
            )
            sigma_a = sigma_b * np.sqrt(np.nanmean(t_axis))
        else:
            sigma_a = None
            sigma_b = None

        t_residual = np.sqrt(np.nanmean(epsilon))

        u_residual = well_data_table[1:, signal_max_lambda] - np.exp(alpha) * np.exp(
            t_axis * beta) * signal_max_absorbance
        u_residual = np.multiply(u_residual, u_residual)
        u_residual = np.sqrt(np.nanmean(u_residual))

        return {"lambda max": x_axis[signal_max_lambda],
                "abs max": signal_max_absorbance,
                "intercept": alpha,
                "slope": beta,
                "intercept uncertainty": sigma_a,
                "slope uncertainty": sigma_b,
                "rmse_transformed": t_residual,
                "rmse_untransformed": u_residual}

    def simple_regression(self, assay_name, on_wells, start_from=350):
        """
        Performs regression of the integrated area of the signal

        :param assay_name: Name of the assay (determined how the data is saved in the project)
        :param on_wells: List of wells on which to perform the regression
        :param start_from: Defines a lower cutoff for wavelength (to avoid integrating the wellplate)
        :return: None
        """
        properties = {well_id: dict() for well_id in self.project.well_ids()}

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue
            data_table = cast_well_data_to_table(self.project, self.times, well_id)
            properties[well_id] = self._simple_regression(data_table, start_from)

        log_entry = f"simple_regression(assay_name={assay_name}, on_wells={on_wells})"
        self.project.commit_property(assay_name, None, properties, log=[log_entry, self.times_log_entry])

    @staticmethod
    def _advanced_regression(well_data_table: np.ndarray, options=None):
        if options is None:
            options = dict()

        lower_lambda_cutoff = options.get("lambda_min", 300)
        upper_lambda_cutoff = options.get("lambda_max", 650)
        max_time_index = options.get("max_time_points", None)
        max_time_total = options.get("max_time", None)
        time_resolve = options.get("time_mode", 'min')
        min_points_to_fit = options.get("min_points", 5)
        max_iterations = options.get('max_iter', 100)
        regressor = options.get('st_regr', 'NNLS')
        verbose = options.get('verbose', False)
        initializer = options.get('components', lambda x, y: [x, y[-1, :] - x, ])
        killjoys = McrAR(max_iter=max_iterations, st_regr=regressor)

        t_axis: np.ndarray = well_data_table[1:, 0].reshape(-1, 1)  # [t, 1]
        l_axis: np.ndarray = well_data_table[0, 1:].reshape(1, -1)  # [1, l]
        signal: np.ndarray = well_data_table[1:, 1:]  # [t, l]  # noqa

        if max_time_total is not None:
            _t_max_i = np.argmax(t_axis > max_time_total)
            _t_max_i = t_axis.size if _t_max_i == 0 else _t_max_i
        else:
            _t_max_i = t_axis.size
        if (max_time_index is None) and (max_time_total is None):
            max_time = 9 + 2
        elif max_time_total is None:
            max_time = max_time_index + 2
        elif max_time_index is None:
            max_time = _t_max_i
        else:
            if time_resolve == "min":
                max_time = min(max_time_index + 2, _t_max_i)
            elif time_resolve == "max":
                max_time = max(max_time_index + 2, _t_max_i)
            else:
                max_time = max_time_index + 2
        max_time = max(max_time, 1 + min_points_to_fit)

        selected_wavelengths = np.where((l_axis[0, :] >= lower_lambda_cutoff) & (l_axis[0, :] <= upper_lambda_cutoff),
                                        True,
                                        False)
        signal = signal[:, selected_wavelengths]  # noqa
        l_axis = l_axis[0, selected_wavelengths].flatten()
        bad_signal_mask = np.where(np.any(np.isnan(signal), axis=0), False, True)
        signal = signal[:, bad_signal_mask] # noqa
        l_axis = l_axis[bad_signal_mask]

        empty_signal_t_mask = np.where(np.all(np.isnan(signal), axis=1), False, True)
        t_axis = t_axis[empty_signal_t_mask, 0].reshape(-1, 1)  # [t, 1]
        signal = signal[empty_signal_t_mask, :]  # [t, l]  # noqa

        t_axis = t_axis[0:max_time, 0]  # [t*, 1]
        signal = signal[0:max_time, :]  # [t*, l]  # noqa
        signal[np.where(signal < 0)] = 0
        init_signal = signal[0, :]

        init_signal = np.array(initializer(init_signal, signal))
        init_signal[init_signal < 0] = 0
        killjoys.fit(signal, ST=init_signal, verbose=verbose)
        extracted_signal = killjoys.ST_opt_[0, :].flatten()
        extracted_signal = extracted_signal  # /np.nanmax(extracted_signal)
        alt_signal = killjoys.ST_opt_[1, :].flatten()
        alt_signal = alt_signal  # /np.nanmax(alt_signal)
        integrated_signal = killjoys.C_opt_[:, 0]
        fitting_iterations = killjoys.n_iter

        # t_axis
        # integrated_signal
        _x_axis = t_axis.flatten()  # [t, ]
        _y_axis = integrated_signal  # [t, 1]
        _y_axis = _y_axis.flatten()  # [t, ]

        negative_mask = np.where(_y_axis < 0, False, True)
        _x_axis = _x_axis[negative_mask]  # [t, ]
        _y_axis = _y_axis[negative_mask]  # [t, ]

        ref_x = _x_axis.tolist()
        ref_y = _y_axis.tolist()

        result = {"fit_limit": 0,
                  "rmse_transformed": 1e12,
                  "intercept": 0,
                  "slope": "Insufficient points for fitting",
                  "x_data": [],
                  "x_data_label": "s",
                  "y_data": [],
                  "y_data_label": "a.u.",
                  "slope uncertainty": 1e12,
                  "extracted_reference": [l_axis.tolist(), extracted_signal.tolist(), alt_signal.tolist()],
                  "MCR_iter": f"{fitting_iterations}/{max_iterations}",
                  }

        for fitting_limit in range(1 + min_points_to_fit, len(_x_axis)):
            # Extract subset of data to fit
            _x_data_raw = _x_axis[1:fitting_limit]  # [t, ]
            _y_data_raw = _y_axis[1:fitting_limit]  # [t, ]
            # Normalize and take ln
            _lny_data_raw = np.log(_y_data_raw / _y_axis[0])
            # Expunge NaN values from x and y
            nan_mask = np.where(np.isnan(_lny_data_raw), False, True)
            _x_data_clipped = _x_data_raw[nan_mask]
            _y_data_clipped = _lny_data_raw[nan_mask]

            fitting_x_data = _x_data_clipped
            fitting_y_data = _y_data_clipped

            if fitting_x_data.size == 0:
                continue

            fit_info = np.polyfit(fitting_x_data, fitting_y_data, deg=1, full=True)
            beta = fit_info[0][0]
            alpha = fit_info[0][1]
            _n = len(fitting_x_data)
            rmse = np.sqrt(fit_info[1][0] / _n)
            r_var_time = np.sqrt(np.nanvar(fitting_x_data))

            if len(fitting_x_data) > 2:
                sigma_m = ((rmse / r_var_time) / np.sqrt(_n - 2))
                # sigma = sqrt( 1/(n-2) * Sum(e^2) // Sum((x - <x>)^2) )
                #       = sqrt( n/(n-2) * (1/n) * Sum(e^2) // Sum((x - <x>)^2) )
                #       = [RMSE] * sqrt( n/(n-2) // Sum((x - <x>)^2) )
                #       = [RMSE] * sqrt( 1/(n-2) // (Sum((x - <x>)^2)/n) )
                #       = [RMSE] * sqrt( 1/(n-2) // Var(x))
                #       = ([RMSE]/sqrt(Var(x))) * sqrt(1/(n-2))
            else:
                sigma_m = 1e13  # must be larger than 1e12 (default value) or the "dry" runs will overwrite defaults

            if sigma_m <= result['slope uncertainty']:
                result.update({"fit_limit": fitting_limit,
                               "rmse_transformed": rmse,
                               "intercept": alpha,
                               "slope": beta,
                               "x_data": ref_x,
                               "y_data": ref_y,
                               "slope uncertainty": sigma_m})
        return result

    def advanced_regression(self, assay_name, on_wells, options=None):
        """
        Performs spectral deconvolution on the data in order to perform 1st order degradation kinetic analysis on the
        pristine dye spectrum.

        **Options**
          * lambda_min (Default: 350) - Ignore signal from wavelengths below lambda_min
          * lambda_max (Default: 650) - Ignore signal from wavelengths above lambda_max
          * max_time_points (Default: 9) - Sets maximum number of time points over which to fit by limiting the index
            (1st time point is used for normalization and is not included in fitting)
          * max_time_total (Default: 10000) - Sets maximum number of time points over which to fit by limiting the value
            of the time axis (1st time point is used for normalization and is not included in fitting)
          * time_mode (Default: 'min') - Used to resolve having both max_time_points and max_time_total set
          * min_points (Default: 5) - Sets the minimum number of time points ver which to fit
          * max_iter (Default: 100) - [McrAR kwarg] Maximum number of fitting iterations to attempt
          * st_regr (Default: NNLS') - [McrAR kwarg] Regression to use on signal (default is Non-negative least squares)
          * verbose (Default: False) - [McrAR.fit kwarg] Boolean flag for iteration details to be printed

        References: :class:`McrAR`

        :param assay_name: Name of the assay (determined how the data is saved in the project)
        :param on_wells: List of wells on which to perform the regression
        :param options: A dictionary of options for the regression
        :return: None
        """
        properties = {well_id: dict() for well_id in self.project.well_ids()}

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue
            data_table = cast_well_data_to_table(self.project, self.times, well_id)
            properties[well_id] = self._advanced_regression(data_table, options)

        log_entry = f"advanced_regression(assay_name={assay_name}, on_wells={on_wells}, options={options})"
        self.project.commit_property(assay_name, None, properties, log=[log_entry, self.times_log_entry])

    def commit_to_database(self, name='AH', assay_name=None, n_retry=4, wait=1,
                           on_wells=None, _iteration=0, campaign=None):
        """ Uploads processed data to the database

        :param name: Logon name (default: 'AH')
        :param assay_name: The name of the assay (for finding all the separate scans)
        :param n_retry: Number of times a commit will be retried
        :param wait: Time between attempts (in seconds)
        :param _iteration: The iteration of the data being committed
        :param on_wells: The wells being committed to the database (falsey - all wells)
        :param campaign: The name of the associated campaign
        :return: success_log, failure_log
        """
        if assay_name is None:
            raise ValueError("Missing required argument 'assay_name'")
        if on_wells is None:
            on_wells = []
        commit_log = list()
        iteration = str(int(_iteration))

        project = self.project.pull_project_copy()
        measurement_mode = project['method_specifications']['mode']
        data_origin = ('measured', 'amd_platform')

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue

            # Preload general/contextual information
            iteration_data = project.get(well_id, {}).get(iteration, {})
            well_context = iteration_data.get('well_context', {})
            try:
                smiles = ".".join(
                    well_context.get(
                        'target_molecules',  # The smiles for the target analytes are in this list
                        [well_context['reagents'][0][2], ]  # Backup: Assume first reagent is the target
                    )
                )
            except (KeyError, IndexError):
                commit_log.append((well_id, None, "target_molecules field failed lookup"))
                continue
            try:
                solvent = well_context['solvents']
            except KeyError:
                commit_log.append((well_id, None, "solvents field failed lookup"))
                solvent = "Solvent not reported"
            meta_data = {'well_context': well_context,
                         'history': project.get('history', ""),
                         'source': self.project.directory}
            try:
                if assay_name == 'photodeg':
                    temperature = [project['method_specifications']['photo_temp'], "C"]
                else:
                    temperature = [average_temperature(project['abscissa']['metadata'][iteration]), "C"]
            except (KeyError, TypeError):
                temperature = ["Temperature not reported", "C"]

            # If peak data requested
            data_name = f"_{assay_name}_kinetic_parameters"
            type_tag = f"exp_kinetic_{measurement_mode}"

            try:
                data = project['calculated_properties'][assay_name][well_id]
                value = {
                    'fit_information': data,
                    'temperature': temperature,
                    'fit_type': "log-linearized_LS",
                    'kinetic_rate': data['slope'] if isinstance(data['slope'], str) else -data['slope'],
                    'kinetic_rate_label': "1/s",
                    'dataset_type': assay_name
                }
            except KeyError as ke:
                commit_log.append((well_id, type_tag,
                                   f"Failed to pull data: key error on '{next(iter(ke.args), '<unknown>')}'"))
                value = {
                    'processed_data': "Error: Data not found",
                    'temperature': temperature,
                }

            for _ in range(n_retry):
                code, _, _ = dri.add_data(name, smiles=smiles, data_name=data_name, campaign_name=campaign,
                                          type_tag=type_tag, data_origin=data_origin,
                                          value=value, solvent=solvent, meta_data=meta_data,
                                          get_details=True)
                if code == "Success":
                    commit_log.append((well_id, type_tag, False))
                    break
                elif code == "Connection Failed":
                    commit_log.append((well_id, type_tag, "Database not online"))
                else:
                    time.sleep(wait)
                    continue
            else:
                commit_log.append((well_id, type_tag, f"Max retries reached for database push"))

        success_log = [item[:2] for item in commit_log if not item[2]]
        failure_log = [item for item in commit_log if item[2]]
        return success_log, failure_log


class LogD:
    def __init__(self, project: sp.SparkProjectFileHandler, tamer: DataTamer):
        self.project = project
        self.tamer = tamer

    @staticmethod
    def consolidate_solvent_data(solvent_data: List[Tuple[str, float, str]]) -> List[Tuple[str, float, str]]:
        memory = dict()
        for solvent_datum in solvent_data:
            _, volume, smiles = solvent_datum
            memory.setdefault(smiles, 0)
            memory[smiles] += volume
        return [[k, v, k] for k, v in memory.items()]  # noqa (want to hint the structure of the list but list can't)

    @staticmethod
    def _process_log_d_by_peaks(org_data, org_x_data, aqu_data, aqu_x_data, **find_peaks_args):
        lower_bound = max(min(org_x_data), min(min(aqu_x_data)))
        upper_bound = min(max(org_x_data), max(min(aqu_x_data)))
        wavelength_axis = np.arange(lower_bound, upper_bound + 1)

        org_x_data = np.array(org_x_data)
        org_y_data = np.array(org_data['y-axis (processed)'])
        aqu_x_data = np.array(aqu_x_data)
        aqu_y_data = np.array(aqu_data['y-axis (processed)'])

        mask = ~np.isnan(org_y_data.astype(float))
        org_signal = (dynamic_scipy.interp1d(org_x_data[mask], org_y_data[mask]))(wavelength_axis)
        mask = ~np.isnan(aqu_y_data.astype(float))
        aqu_signal = (dynamic_scipy.interp1d(aqu_x_data[mask], aqu_y_data[mask]))(wavelength_axis)
        total_signal = org_signal + aqu_signal

        mask = ~np.isnan(total_signal.astype(float))
        peak_indices, *_ = scipy.signal.find_peaks(total_signal[mask], **find_peaks_args)
        peak_lambda: np.ndarray = wavelength_axis[mask][peak_indices]
        prominence, *_ = scipy.signal.peak_prominences(total_signal[mask], peak_indices)

        temp = np.array([peak_lambda, prominence])
        temp = temp[:, (-temp[1]).argsort()]
        main_peak_wavelength = temp[0][0]
        main_peak_prominence = temp[1][0]

        p0 = [10.0, main_peak_prominence]

        mask = ~np.isnan(org_signal.astype(float))
        coeff_org, var_org = curve_fit(lambda x, *p: gauss(x, main_peak_wavelength, *p),
                                       wavelength_axis[mask], org_signal[mask], p0=p0)
        var_org = np.diag(var_org)

        mask = ~np.isnan(aqu_signal.astype(float))
        coeff_aqu, var_aqu = curve_fit(lambda x, *p: gauss(x, main_peak_wavelength, *p),
                                       wavelength_axis[mask], aqu_signal[mask], p0=p0)
        var_aqu = np.diag(var_aqu)

        distribution_coefficient = (coeff_org[0] * coeff_org[1]) / (coeff_aqu[0] * coeff_aqu[1])
        rel_sigma_d = np.sqrt(
            var_org[0] / np.power(coeff_org[0], 2)
            + var_org[1] / np.power(coeff_org[1], 2)
            + var_aqu[0] / np.power(coeff_aqu[0], 2)
            + var_aqu[1] / np.power(coeff_aqu[1], 2)
        )
        log_d = np.log10(distribution_coefficient)
        sigma_log_d = rel_sigma_d * (1.0 / np.log(10))

        return log_d, sigma_log_d

    def process_log_d_by_peaks(self, on_wells, **find_peaks_args):
        properties = {
            well_id: {
                'logD': None,
                'sigma': None,
            } for well_id in self.project.well_ids()
        }

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue
            data_org, _, _, x_data_aqu = self.project.get_datum_by_well(well_id, 0)
            data_aqu, _, _, x_data_org = self.project.get_datum_by_well(well_id, 1)

            result = self._process_log_d_by_peaks(data_org, x_data_org, data_aqu, x_data_aqu, **find_peaks_args)
            properties[well_id]['logD'] = result[0]
            properties[well_id]['sigma'] = result[1]

        self.project.commit_property("LogD", None, properties,
                                     log=f"process_log_d_by_peaks("
                                         f"on_wells={on_wells}, find_peaks_args={find_peaks_args}")

    @staticmethod
    def _process_log_d_by_integral(org_data, org_x_data, aqu_data, aqu_x_data):
        lower_bound = max(min(org_x_data), min(min(aqu_x_data)))
        upper_bound = min(max(org_x_data), max(min(aqu_x_data)))
        wavelength_axis = np.arange(lower_bound, upper_bound + 1)

        org_x_axis = np.array(org_x_data)
        aqu_x_axis = np.array(aqu_x_data)
        org_y_axis = np.array(org_data['y-axis (processed)'])
        aqu_y_axis = np.array(aqu_data['y-axis (processed)'])

        mask = ~np.isnan(org_y_axis.astype(float))
        org_signal = (dynamic_scipy.interp1d(org_x_axis[mask], org_y_axis[mask]))(wavelength_axis)
        mask = ~np.isnan(aqu_y_axis.astype(float))
        aqu_signal = (dynamic_scipy.interp1d(aqu_x_axis[mask], aqu_y_axis[mask]))(wavelength_axis)

        delta_x = np.diff(wavelength_axis)

        org_integral, org_error = trapz_error(org_signal, wavelength_axis, min(delta_x))
        aqu_integral, aqu_error = trapz_error(aqu_signal, wavelength_axis, min(delta_x))

        distribution_coefficient = (org_integral / aqu_integral)
        rel_sigma_d = np.sqrt(np.power(aqu_error / aqu_integral, 2) + np.power(org_error / org_integral, 2))
        log_d = np.log10(distribution_coefficient)
        sigma_log_d = rel_sigma_d * (1.0 / np.log(10))

        return log_d, sigma_log_d

    def process_log_d_by_integral(self, on_wells):
        properties = {
            well_id: {
                'logD': None,
                'sigma': None,
            } for well_id in self.project.well_ids()
        }

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue
            data_org, _, _, x_data_org = self.project.get_datum_by_well(well_id, 0)
            data_aqu, _, _, x_data_aqu = self.project.get_datum_by_well(well_id, 1)

            result = self._process_log_d_by_integral(data_org, x_data_org, data_aqu, x_data_aqu)
            properties[well_id]['logD'] = result[0]
            properties[well_id]['sigma'] = result[1]

        self.project.commit_property("LogD", None, properties,
                                     log=f"process_log_d_by_integral(on_wells={on_wells})")

    def commit_to_database(self, name='AH', n_retry=4, wait=1, on_wells=None, campaign=None):
        if on_wells is None:
            on_wells = []
        commit_log = list()

        project = self.project.pull_project_copy()
        data_origin = ('measured', 'amd_platform')

        for well_id in self.project.well_ids():
            if on_wells and (well_id not in on_wells):
                continue

            # Preload general/contextual information
            iteration_data_0 = project.get(well_id, {}).get("0", {})
            well_context_0 = iteration_data_0.get('well_context', {})
            iteration_data_1 = project.get(well_id, {}).get("1", {})
            well_context_1 = iteration_data_1.get('well_context', {})
            try:
                smiles_0 = ".".join(
                    sorted(
                        well_context_0.get(
                            'target_molecules',  # The smiles for the target analytes are in this list
                            [well_context_0['reagents'][0][2], ]  # Backup: Assume first reagent is the target
                        )
                    )
                )
                smiles_1 = ".".join(
                    sorted(
                        well_context_1.get(
                            'target_molecules',  # The smiles for the target analytes are in this list
                            [well_context_1['reagents'][0][2], ]  # Backup: Assume first reagent is the target
                        )
                    )
                )
            except (KeyError, IndexError):
                commit_log.append((well_id, None, "target_molecules field failed lookup"))
                continue
            if smiles_0 != smiles_1:
                commit_log.append((well_id, None, "target molecule(s) mismatch between phases"))
            smiles = smiles_0
            try:
                solvent_0: list = well_context_0['solvents']
                solvent_1: list = well_context_1['solvents']
            except KeyError:
                commit_log.append((well_id, None, "solvents field failed lookup"))
                solvent = "Solvent not reported"
            else:
                solvent = self.consolidate_solvent_data(solvent_0 + solvent_1)
            meta_data = {'well_context (phase 0)': well_context_0,
                         'well_context (phase 1)': well_context_1,
                         'history': project.get('history', ""),
                         'source': self.project.directory}
            try:
                temperature_0 = [average_temperature(project['abscissa']['metadata']["0"]), "C"]
            except (KeyError, TypeError):
                temperature_0 = ["Temperature not reported", "C"]
            try:
                temperature_1 = [average_temperature(project['abscissa']['metadata']["1"]), "C"]
            except (KeyError, TypeError):
                temperature_1 = ["Temperature not reported", "C"]
            meta_data["temperatures"] = {"org": temperature_0,
                                         "aqu": temperature_1}

            # If peak data requested
            data_name = f"_shakeflask_logD"
            type_tag = f"exp_logp_abs"

            try:
                data = project['calculated_properties']['LogD'][well_id]
                log_d = data['logD']
                sigma = data['sigma']
                pH = 7.4  # noqa (that is pH)
                value = {
                    'logp_value': [pH, log_d, sigma],
                    'temperature': [None, "C"],
                }
            except KeyError as ke:
                commit_log.append((well_id, type_tag,
                                   f"Failed to pull data: key error on '{next(iter(ke.args), '<unknown>')}'"))
                value = {
                    'logp_value': "Error: Data not found",
                    'temperature': [None, "C"],
                }

            for _ in range(n_retry):
                code, _, _ = dri.add_data(name, smiles=smiles, data_name=data_name, campaign_name=campaign,
                                          type_tag=type_tag, data_origin=data_origin,
                                          value=value, solvent=solvent, meta_data=meta_data,
                                          get_details=True)
                if code == "Success":
                    commit_log.append((well_id, type_tag, False))
                    break
                elif code == "Connection Failed":
                    commit_log.append((well_id, type_tag, "Database not online"))
                else:
                    time.sleep(wait)
                    continue
            else:
                commit_log.append((well_id, type_tag, f"Max retries reached for database push"))

        success_log = [item[:2] for item in commit_log if not item[2]]
        failure_log = [item for item in commit_log if item[2]]
        return success_log, failure_log


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def sigma_f(f, var, sig, dx=0.01):
    accumulator = 0
    for i, _ in enumerate(var):
        v_p2 = list(var)
        v_p1 = list(var)
        v_m1 = list(var)
        v_m2 = list(var)
        v_p2[i] = v_p2[i] + 2 * dx  # noqa
        v_p1[i] = v_p1[i] + dx
        v_m1[i] = v_m1[i] - dx
        v_m2[i] = v_m2[i] - 2 * dx  # noqa
        df_dvar = (- 1.0 * f(v_p2)  # noqa
                   + 8.0 * f(v_p1)
                   - 8.0 * f(v_m1)
                   + 1.0 * f(v_m2)
                   ) / (dx * 12)
        accumulator += np.power(df_dvar * sig[i], 2)
    return np.sqrt(accumulator)


def subtract_background(self, background_keys: list, bkg_sub_factor=1.00, use_zero_floor=False):
    background_experiments = {k: v for k, v in self.data.items() if k[0] in background_keys}
    mean_background_signals = dict()  # [(phase, sample_volume), __]

    for k in background_experiments.keys():
        index = (k[1], k[2])
        mean_background_signals[index] = mean_background_signals.get(index, dict())
        mean_background_signals[index]['sum'] = mean_background_signals[index].get('sum', 0)
        mean_background_signals[index]['sum'] += background_experiments[k]
        mean_background_signals[index]['count'] = mean_background_signals[index].get('count', 0)
        mean_background_signals[index]['count'] += 1
    for k in mean_background_signals.keys():
        mean_background_signals[k]['mean'] = mean_background_signals[k]['sum'] / mean_background_signals[k]['count']

    for k, v in self.data.items():
        sel_bkg = mean_background_signals[(k[1], k[2])]['mean']
        self.data[k][:, 1] = self.data[k][:, 1] - bkg_sub_factor*sel_bkg[np.isin(sel_bkg[:, 0], self.data[k][:, 0]), 1]
        if use_zero_floor:
            self.data[k][np.where(self.data[k][:, 1] <= 0), 1] = 0


def background_subtraction(self, signals: dict, background_keys: list, bkg_sub_factor=1.00, use_zero_floor=False):
    if signals is None:
        signals = self.time_series_data
    t_axis = signals[background_keys[0]][1:, 0]
    l_axis = signals[background_keys[0]][0, 1:]
    height = len(background_keys)

    background_signals = np.zeros([height, len(t_axis), len(l_axis)])

    for i, k in enumerate(background_keys):
        background_signals[i][:, :] = signals[k][1:, 1:]

    mean_background = np.nanmean(background_signals, axis=0)

    t_bar = np.nanmean(np.multiply(t_axis.reshape(-1, 1), np.divide(mean_background, mean_background)), axis=0)
    t_minus_t_bar = t_axis.reshape(-1, 1) - t_bar.reshape(1, -1)

    s_bar = np.nanmean(mean_background, axis=0)
    s_minus_s_bar = mean_background - s_bar

    beta = np.nansum(np.multiply(t_minus_t_bar, s_minus_s_bar), axis=0) / np.nansum(
        np.multiply(t_minus_t_bar, t_minus_t_bar), axis=0)  # noqa
    alpha = s_bar - beta * t_bar

    for k in signals.keys():
        signals[k][1:, 1:] -= (alpha.reshape(1, -1) + beta * t_axis.reshape(-1, 1)) * bkg_sub_factor
        if use_zero_floor:
            signals[k][1:, 1:] = np.where(signals[k][1:, 1:] < 0, 0, signals[k][1:, 1:])
