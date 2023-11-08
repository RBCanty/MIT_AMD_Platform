#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Based on code by MattMcDonald (Created on Wed Jan 12 11:32:05 2022)
Uses interface specification by Brent Koscher
Adapted for MCN by Ben Canty

Generalized Data DB interface code similar to database_interface.py but for the data database.
"""

import datetime
import json
import requests
import custom_exceptions as cexc
from pprint import pformat
from mcn_status import MCN_CFG, S_DATABASES
from typing import List, Union, Tuple

Number = Union[int, float]

DATABASE_URL, DATABASE_PORT = MCN_CFG[S_DATABASES]['Data']


def update_database(url, port):
    """ Updates the URL and PORT of the database for the interface

    :param url: The new url for the database (e.g. 'http://01.23.45.678:12345')
    :param port: The new port for the database (e.g. '67890')
    :return: None
    """
    if (not url) or (not port):
        return
    global DATABASE_URL, DATABASE_PORT
    try:
        DATABASE_URL = str(url)
        DATABASE_PORT = str(int(port))
    except (TypeError, ValueError):
        DATABASE_URL, DATABASE_PORT = MCN_CFG[S_DATABASES]['Data']


def database_request(request: dict, user: str, domain=None, port=None, get_details=False):
    """ Convenience function for making database requests

    References: :class:`cexc.DatabaseRequestError`

    :param request: A dictionary of the proper request
    :param user: Used for authentication
    :param domain: (e.g. 'http://01.23.45.678:12345')
    :param port: (e.g. 67890)
    :param get_details: If false, raise DatabaseRequestError; if True, return requests post responses
    :return: Document, Details, Response Code
    """
    if domain is None:
        domain = DATABASE_URL
    if port is None:
        port = DATABASE_PORT
    port = str(port)
    user = str(user).lower()
    username = f"{user}_computer"
    password = f"amd_{user}_user"

    incoming_request = {'auth': {'port': port, 'user': username, 'password': password},
                        'request': request}

    try:
        resp = requests.post(domain, json=incoming_request, headers={'content-type': 'application/authjson'})
    except requests.exceptions.ConnectionError:
        if get_details:
            return "Connection Failed", "", ""
        else:
            raise cexc.DatabaseRequestError(f"Failed to connect to database, check connection")

    req_doc, err_det = json.loads(resp.content.decode('utf-8'))
    resp_code = resp.status_code

    if (resp_code != 200) or (req_doc == "Error"):
        if get_details:
            return req_doc, err_det, resp.status_code
        else:
            raise cexc.DatabaseRequestError(f"Request {pformat(request, compact=True)} by {user} failed:\n"
                                            f'----Response Object: Start----\n'
                                            f'{pformat(req_doc, compact=True)}\n\n'
                                            f'{pformat(err_det, compact=True)}\n'
                                            f'----Response Object: End----')

    return req_doc, err_det, resp.status_code


def complete_queue(campaign_name: str, model_details: dict):
    """ Signals the data database that a queue is complete for proper model retraining

    :param campaign_name: Name of the campaign the data belongs to
    :param model_details: Details for updated the database (type: dict) with keys for each model (type: string).
      Basically the only field currently is status (type: string) which the HTTPServer uses to keep track of the job.
    :return:
    """
    dr_request = {
        'request_type': 'complete_data_campaign',  # Name of the HTTPServer function (type: string)
        'campaign_name': campaign_name,
        'model_details': model_details,
    }
    return database_request(dr_request, "MC")


def add_data(name, *, smiles: str, data_name: str,
             type_tag: str, data_origin: Tuple[str, str],
             value: dict, solvent: List[Tuple[str, Number, str]],
             molecular_uncertainty: str = None,
             meta_data: dict = None, campaign_name: str = "",
             get_details=False):
    """
    Used to add data to the data repository

    Type Tags:
      * **Spectral:**
      * calc_abs_spec
      * calc_abs_stick
      * calc_pl_spec
      * calc_pl_stick
      * calc_ir_spec
      * calc_ir_stick
      * exp_abs_spec
      * exp_pl_spec
      * exp_ir_spec
      * exp_ir_stick
      * **LogP:**
      * exp_logp_abs
      * exp_logp_hplc
      * exp_logp_lit
      * calc_logp

    Data Origin:
      - 'literature', 'measured', or 'calculated'
      -  location specifier (ex. 'DOI:XX', 'amd_platform', etc.)

    Value requirements:
      * **Spectral:**
      * xdatalabel - [label name, units]
      * xdata - List
      * ydatalabel - [label name, units]
      * ydata - List of Lists
      * ir_support_id - (Optional, string: for use with IR)
      * **LogP:**
      * logp_value - List of Lists[pH, value, uncertainty]
      * temperature - (Optional, string)
      * **Kinetic**
      * fit_information - A dictionary of the regression details (should at least have slope, intercept, and rmse)
      * temperature - [temperature, units]
      * fit_type - How the data was fit (e.g. "log-linearized_LS")
      * kinetic_rate - A number (or text explaining why no number)
      * kinetic_rate_label: (e.g. "1/s")

    References: :class:`cexc.DatabaseRequestError`

    :param name: Name of system making request (for log-on)
    :param smiles: SMILES string for the compound being attributed data
    :param data_name: Name for the data (internally prefixed with SMILES string)
    :param campaign_name: Name of the campaign with which the data is associated
    :param type_tag: Identifies the data (pl spectra vs logD vs hyperpolarizability etc.)
    :param data_origin: Tag for the source of the data
    :param value: A dictionary of the required data fields for a given measurement
    :param solvent: A list of solvent entries [name, volume, smiles]
    :param molecular_uncertainty: (Optional) Allows to specify purchased, synthesized, mass-hit, etc.
    :param meta_data: (Optional) Custom field
    :param get_details: (Optional boolean) True - returns failure details; False - raises DatabaseRequestError
      on bad request
    :return: (Data document, details, response code)

    :raises cexc.DatabaseRequestError: If connection fails or if bad request (from database_request() method)
    """
    if meta_data is None:
        meta_data = {}
    my_request = {'request_type': 'add_data',
                  'molecule': [smiles, 'smiles'],
                  'add_data': True,
                  'data_to_add': {'data_name': f"{smiles}{data_name}",
                                  'campaign_name': campaign_name,
                                  'type_tag': type_tag,
                                  'data_origin': data_origin,
                                  'origin_date': datetime.datetime.now().strftime('%m/%d/%Y'),
                                  'solvent': solvent,
                                  'metadata': meta_data}}
    my_request['data_to_add'].update(value)
    if molecular_uncertainty:
        my_request['data_to_add']['molecular_uncertainty'] = molecular_uncertainty
    return database_request(my_request, name, get_details=get_details)


def test_connection():
    """ Tests connectivity to database

    :return: "Connected" or "No connection"
    """
    test_request = {'request_type': 'add_reaction_bad',
                    'reaction_type': 'synthesis',
                    'reaction_smiles': 'garbage',
                    'reaction_details': 'garbage'}
    try:
        result, _, _ = database_request(test_request, "MC", get_details=True)
    except Exception as e:  # noqa
        print("Results Database Test Exception", repr(e))
        return "No connection"
    if result == "Connection Failed":
        return "No connection"
    return "Connected"


if __name__ == '__main__':
    print(test_connection())
