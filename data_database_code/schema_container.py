# -*- coding: utf-8 -*-
"""
Schema templates for the AMD Database and query capabilities
Created on Wed Dec 18 09:41:37 2019

@author: Brent
"""


def get_schema(schema_type):
    if schema_type == "initial_entry":
        schema_to_return = {
                "type": "object",
                "properties": {
                        "compound_inchi_key" : {"type" : "string"},
                        "date_created" : {"type" : "string"},
                        "last_modified" : {"type" : "string"},
                        "basic_properties": {
                                "type": "object",
                                "properties": {
                                            "compound_inchi": {"type" : "string"},
                                            "compound_inchi_fixedh": {"type" : "string"},
                                            "compound_smiles": {"type": "string"},
                                            "compound_mformula": {"type": "string"},
                                            "compound_mweight": {"type": "number"},
                                            "compound_cas": {"type": "string"},
                                            "other_identifiers": {"type": "array"},
                                            "compound_iupac_name": {"type": "string"}
                                            },
                                "required": ["compound_inchi","compound_smiles","compound_mformula","compound_mweight"],
                                "additionalProperties": False
                                },
                        "phys_chem_props": {
                                "type": "object",
                                "properties":{
                                        "boiling_point": {"type": "array"},
                                        "density": {"type": "array"},
                                        "melting_point": {"type": "array"},
                                        "refractive_index": {"type": "array"},
                                        "vapor_pressure": {"type": "array"}
                                        },
                                "additionalProperties": False
                                },
                        "associated_reactions": {"type": "object"}
                },
                "required": ["compound_inchi_key","date_created","last_modified"]
                }
    
    elif schema_type == 'reaction_entry':
        schema_to_return = {
                "type": "object",
                "properties": {
                        }
                }
    elif schema_type == 'add_reaction_details':
        schema_to_return = {
                "type": "object",
                "properties": {
                        "reaction_reactants": {"type": "array"},
                        "reaction_reagents": {"type": "array"},
                        "reaction_solvents": {"type": "array"},
                        "reaction_catalysts": {"type": "array"},
                        "reaction_products": {"type": "array"},
                        "reaction_temperature": {"type": "number"},
                        "reaction_duration": {"type": "string"},
                        "reaction_template": {"type": ["string", "array"]},
                        "reaction_location": {"type": "array"},
                        "reaction_date": {"type": "string"},
                        "reaction_metadata": {"type": ["string", "array", "object"]}
                        },
                "required": ["reaction_date", "reaction_duration", "reaction_location",
                             "reaction_reagents", "reaction_solvents", "reaction_temperature",
                             "reaction_reactants", "reaction_products"],
                "additionalProperties": False
                }
                            
    elif schema_type == 'spectral_data':
        schema_to_return = {
                "type": "object",
                "properties": {
                        "data_name": {"type" : "string"},
                        "campaign_name": {"type": "string"},
                        "date_added": {"type" : "string"},
                        "type_tag": {"type" : "string"},
                        "molecule_certainty": {"type" : "string"},
                        "data_origin": {"type" : "array"},
                        "origin_date": {"type" : "string"},
                        "temperature": {"type" : ["array", "number"]},
                        "xdatalabel": {"type" : "array"},
                        "xdata": {"type" : "array"},
                        "ydatalabel": {"type" : "array"},
                        "ydata": {"type" : "array"},
                        "metadata": {"type" : ["string", "array", "object"]},
                        "solvent": {"type" : ["array", "string"]},
                        "ir_support_id": {"type" : "string"}
                        },
                "required": ["data_name","type_tag","data_origin","origin_date","xdatalabel","xdata","ydatalabel","ydata","solvent"],
                "additionalProperties": False
                }
    elif schema_type == 'kinetic_data':
        schema_to_return = {
                "type": "object",
                "properties": {
                        "data_origin": {"type" : "array"},
                        "origin_date": {"type" : "string"},
                        "campaign_name": {"type": "string"},
                        "source": {"type" : "string"},
                        "data_name": {"type" : "string"},
                        "date_added": {"type" : "string"},
                        "type_tag": {"type" : "string"},
                        "solvent": {"type" : ["array", "string"]},
                        "fit_information": {"type": "object"},
                        "temperature": {"type" : ["array", "number"]},
                        "fit_type": {"type" : "string"},
                        "kinetic_rate": {"type" : ["array", "number", "string"]},
                        "kinetic_rate_label": {"type" :["array", "string"]},
                        "dataset_type": {"type" : "string"},
                        "molecule_certainty": {"type" : "string"},
                        "metadata": {"type" : ["string", "array", "object"]},
                        },
                "required": ["data_name", "type_tag", "data_origin", "origin_date",
                            "solvent", "fit_information", "fit_type", "temperature",
                            "kinetic_rate", "kinetic_rate_label", "dataset_type"],
                "additionalProperties": False
                }
    elif schema_type == 'logp_data':
        schema_to_return = {
                "type": "object",
                "properties": {
                        "data_name": {"type" : "string"},
                        "date_added": {"type" : "string"},
                        "campaign_name": {"type": "string"},
                        "type_tag": {"type" : "string"},
                        "molecule_certainty": {"type" : ["string", "array"]},
                        "data_origin": {"type" : "array"},
                        "origin_date": {"type" : "string"},
                        "logp_value": {"type" : "array"},
                        "temperature": {"type" : ["array", "number"]},
                        "metadata": {"type" : ["string", "array", "object"]},
                        "solvent": {"type" : ["string", "array"]}
                        },
                "required": ["data_name","type_tag","data_origin","origin_date","logp_value","solvent"],
                "additionalProperties": False
                }
    
    elif schema_type == 'loge_data':
        schema_to_return = {
                "type": "object",
                "properties": {
                        "data_name": {"type" : "string"},
                        "date_added": {"type" : "string"},
                        "campaign_name": {"type": "string"},
                        "type_tag": {"type" : "string"},
                        "data_origin": {"type" : "array"},
                        "origin_date": {"type" : "string"},
                        "loge_value": {"type": "object"},
                        "metadata": {"type" : ["string", "array", "object"]}
                        },
                "required": ["data_name","type_tag","origin_date","loge_value"],
                "additionalProperties": False
                }            
    
    elif schema_type == 'query_data_database':
        schema_to_return = {
                "type": "object",
                "properties":{
                        "request_type": {"type": "string"},
                        "molecule": {"type": "array"},
                        "return_data": {"type": "boolean"}
                        },
                "required": ["request_type","molecule","return_data"]
                }
    return schema_to_return



