import requests
import pandas as pd
import os
import time
import json
import sys


USER_AGENT = "Mozilla/5.0"
page_size = 1000
cursor_mark = '*'
format_type = 'json'
base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ABSTRACT:\" {} \"&resultType=core&cursorMark={}&pageSize={}&format={}"
relevant_columns = ['pmid', 'pmcid', 'doi', 'title', 'pubYear']


def search_identifier(symbol, base_url=base_url, cursor_mark=cursor_mark, page_size=page_size, format_type=format_type, relevant_columns=relevant_columns, test=False):
    """Searches EuropePMC for the provided identifier, returns a list with all data for matches"""
    # Initialize the list to store retracted articles
    matches = []
    # Make the initial request to get the total number of results
    url = base_url.format(symbol, cursor_mark, page_size, format_type)
    try:
        response = requests.get(url).json()
        total = response['hitCount']
        # Check if there are any results
        if test==True:
            res_test = len(response['resultList']['result'])
            if res_test >0:
                return f'Getting OK response for test: {symbol}'
            else:
                return
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
    # Get the intersection of columns between the DataFrame and relevant_columns
    df = pd.DataFrame(matches)
    common_columns = list(set(df.columns) & set(relevant_columns))
    # Create a new DataFrame with only the common columns
    try:
        res = df.loc[:, common_columns]
    except KeyError as e:
        #print(f"Error: {e} column not found.")
        res = df
    res = res.fillna("NA")
    return res




def multiple_search(hgnc_df):
    """
    Searches for all HGNCs in the provided dataset and builds a DataFrame with unique hits.
    """
    start_time = time.time()
    exclude = ['NA', 'Entry Withdrawn']
    all_sym = list(set(hgnc_df['primary_symbol'])) + list(set(hgnc_df['secondary_symbol']))
    results = {
        'version': 0,
        'count': 0,
        'result': {'primary_symbol': [], 'secondary_symbol': []}}
    
    total = len(all_sym)
    i = 0

    if not os.path.exists('results'):
        os.mkdir('results')

    primary_symbol_dir = 'results/primary_symbol'
    if not os.path.exists(primary_symbol_dir):
        os.mkdir(primary_symbol_dir)

    secondary_symbol_dir = 'results/secondary_symbol'
    if not os.path.exists(secondary_symbol_dir):
        os.mkdir(secondary_symbol_dir)

    print('Searching...')

    for symbol in all_sym:
        iteration_time = time.time()
        elapsed_time = iteration_time - start_time
        sys.stdout.flush()
        print(f'Symbol {i} of {total}  {symbol}                                         ', end='\r', flush=True)
        i += 1

        if symbol not in exclude:
            search = search_identifier(symbol)
            
            if symbol in hgnc_df['primary_symbol'].values:
                stype = 'primary_symbol'
            else:
                stype = 'secondary_symbol'

            # Process the search results here
            results['result'][stype].append({symbol: search.transpose().to_dict()})
            results['count'] += 1

    with open('results/hgnc_in_abstracts.json', 'w') as f:
        json.dump(results, f)

    print(f'\nElapsed time: {round(elapsed_time, 0)}s                               ', flush = True)

def main():
    # Load the dataset
    hgnc_df = pd.read_csv('data/hgnc.tsv', sep='\t')
    print('TEST: MTOR')
    print(search_identifier('MTOR', test = True))
    print("--------------------\nBegin search")
    # Searching for publications matching all identifiers in the table
    multiple_search(hgnc_df=hgnc_df)
    print("Done searching\n", flush = True)
    print('Results saved to results/hgnc_in_abstracts.json')

if __name__ == '__main__':
    main()
