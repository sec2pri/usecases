#!/usr/bin/env python
"""
EuropePMC Abstract Bibliometric Search Script for HGNC symbols

This script performs HGNC gene symbol searches in article abstracts
using the EuropePMC web services and stores the literature hit 
metadata as results in JSON format.

Author: Javier Millan Acosta
ORCID: 0000-0002-4166-7093
Affiliation: Maastricht University
License: #TODO

Usage:
    python bibliometric_search.py [database_choice]

    [database_choice]: Comma-separated list of database choices defined in the configuration file.

"""

import requests
import pandas as pd
import os
import time
import json
import sys
import yaml
import argparse

CONFIG = 'config.yml'
USER_AGENT = "Mozilla/5.0"
page_size = 1000
cursor_mark = '*'
format_type = 'json'
base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ABSTRACT:\" {} \"&resultType=core&cursorMark={}&pageSize={}&format={}"



def check_config_file(choices, config_file_path=CONFIG):
    """
    Check if the specified sources exist in the config file.

    Args:
        choices (list): List of source choices to check.
        config_file_path (str): Path to the config file.

    Raises:
        SystemExit: Exits with an error message if a source is not found in the config.
    """
    # Check if the config file exists
    if not os.path.isfile(config_file_path):
        print(f"Config file '{config_file_path}' does not exist.")
        exit(1)
    else:
        config_file = open(config_file_path).read()
        config_data = yaml.safe_load(config_file)
        for choice in choices:
            if choice not in config_data.keys():
                print(f'Invalid source: {choice}')
                exit(1)
            
def read_config_file(yaml_file_path, choice):
    """
    Read and validate the configuration data for a specific choice.

    Args:
        yaml_file_path (str): Path to the YAML config file.
        choice (str): The selected choice.

    Returns:
        dict: The validated configuration data.

    Raises:
        SystemExit: Exits with an error message if the configuration is invalid.
    """
    try:
        with open(yaml_file_path, 'r') as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        # Check if required keys exist in the YAML data
        required_keys = ['path', 'firstColumn', 'firstColumn_alt', 'secondColumn', 'secondColumn_alt', 'keep_columns', 'result']
        for key in required_keys:
            if key not in yaml_data[choice]:
                print(f"Required key '{key}' not found in the YAML data.")
                exit(1)
        return yaml_data

    except Exception as e:
        print(f"An error occurred while processing the config file: {str(e)}")
        exit(1)


def ensure_directory_exists(directory_path, config=CONFIG):
    # Check if the directory already exists
    if os.path.exists(directory_path):
        return True  # Directory already exists
    else:
        print(f"The directory '{directory_path}' does not exist.")
        try:
            os.makedirs(directory_path)
            print(f"Directory '{directory_path}' created successfully.")
        except Exception as e:
            print(f"Error creating directory: {e}")
            exit(1)


def search_identifier(symbol, keep_columns, base_url=base_url, cursor_mark=cursor_mark, page_size=page_size, format_type=format_type, ):
    """
    Searches EuropePMC for the provided identifier and returns a DataFrame with matching results.

    Args:
        symbol (str): Identifier to search for.
        keep_columns (list): List of columns to keep in the result DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing search results.
    """
    # Initialize the list to store retracted articles
    matches = []
    # Make the initial request to get the total number of results
    url = base_url.format(symbol, cursor_mark, page_size, format_type)
    try:
        response = requests.get(url).json()
        total = response['hitCount']

        if 'resultList' not in response:
            print("No results found.")
            return pd.DataFrame()
        # Calculate the number of requests needed to retrieve all results
        matches.extend(response['resultList']['result'])
        num_requests = (total + page_size - 1) // page_size
        i = 0
        # Iterate through each page and append the results to the list
        while response['nextCursorMark'] is not None: 
            i+=1
            cursor_mark = response['nextCursorMark']
            url = base_url.format(symbol, cursor_mark, page_size, format_type)
            try:
                response = requests.get(url).json()
                matches.extend([i for i in response['resultList']['result'] if 'Human' in str(response['resultList']['result'][i]['meshHeadingList'])])
            except Exception as e:
                #print("An error occurred: " + str(e))
                break   
    except Exception as e:
        #print("An error occurred: " + str(e), {symbol})
        pass
        # Get the intersection of columns between the DataFrame and keep_columns
        df = pd.DataFrame(matches)
        res = df.fillna("NA")

        # Define a function to retrieve nested column values
        def get_nested_value(row, column):
            keys = column.split('.')
            nested_value = row
            for key in keys:
                if key in nested_value:
                    nested_value = nested_value[key]
                else:
                    return "NA"  
            return nested_value

        # Loop through the keep_columns and apply the get_nested_value function
        for col in keep_columns:
            res[col] = res.apply(lambda row: get_nested_value(row, col), axis=1)

        # Drop the original columns that were replaced by nested values
        keep_columns = [i if '.' not in i else i.split('.')[-1] for i in keep_columns]
        for col in res.columns:
            # Check if the column is not in keep_columns and is not a nested column
            if col not in keep_columns and '.' not in col:
                res = res.drop(col, axis=1)
        print(f'Columns retrieved: {", ".join(res.columns)}\n_____')
        return res

def multiple_search(df, results_path, config):
    """
    Searches for all symbols in the provided dataset and builds a DataFrame with unique hits.

    Args:
        df (pd.DataFrame): Input dataset.
        results_path (str): Path to save the results.
        config (dict): Configuration data.
    """
    start_time = time.time()
    exclude = ['NA', 'Entry Withdrawn']
    firstColumn_alt = config['firstColumn_alt']
    secondColumn_alt = config['secondColumn_alt']
    firstColumn = config['firstColumn']
    secondColumn = config['secondColumn']
    keep_columns = config['keep_columns']
    try:
        if secondColumn_alt != "":
            all_sym = list(set(df[firstColumn_alt])) + list(set(df[secondColumn_alt]))
        else:
            all_sym = list(set(df[firstColumn_alt]))
    except KeyError as e:
        print(f'Key Error {e}: {df.keys()}')
        exit(0)
    results = {
        'version': 0,
        'count': 0,
        'result': {'firstColumn_alt': [], 'secondColumn_alt': []}}
    
    total = len(all_sym)
    i = 0

    if not os.path.exists('results'):
        os.mkdir('results')

    firstColumn_alt_dir = 'results/firstColumn_alt'
    if not os.path.exists(firstColumn_alt_dir):
        os.mkdir(firstColumn_alt_dir)

    secondColumn_alt_dir = 'results/secondColumn_alt'
    if not os.path.exists(secondColumn_alt_dir):
        os.mkdir(secondColumn_alt_dir)

    print('Searching...')

    for symbol in all_sym:
        iteration_time = time.time()
        elapsed_time = iteration_time - start_time
        sys.stdout.flush()
        print(f'Symbol {i} of {total}\n{symbol}\n', end='\n', flush=True)
        i += 1

        if symbol not in exclude:
            search = search_identifier(symbol, keep_columns=keep_columns,)
            
            if symbol in df[firstColumn_alt].values:
                stype = 'firstColumn_alt'
            else:
                stype = 'secondColumn_alt'
            try:
                # Process the search results here
                results['result'][stype].append({symbol: search.transpose().to_dict()})
                results['count'] += 1
            except Exception as e:
                print(f'No results ({e})')

    with open(results_path, 'w') as f:
        json.dump(results, f)

    print(f'\nElapsed time: {round(elapsed_time, 0)}s                               ', flush = True)

def main():
    print('bibliometric_sec2prim.py\n_______________\n')
    parser = argparse.ArgumentParser(description="Check if a database key exists in the config file.")
    parser.add_argument("database", type=str, help="Database key to check")
    sources = parser.parse_args().database.split(',')
    # Check directories, files
    check_config_file(sources)
    for choice in sources:
        config_choice = read_config_file(choice=choice, yaml_file_path=CONFIG)[choice]
        data = config_choice['path']
        results = config_choice['result']
        # Load the dataset
        df = pd.read_csv(data, sep='\t',)
        print("_______________\nBegin search")
        # Searching for publications matching all identifiers in the table
        multiple_search(df=df, results_path=results, config=config_choice,)  # Pass keep_columns as an argument
        print("Done searching\n", flush=True)
        print(f'Results saved to {results}')

if __name__ == '__main__':
    main()
