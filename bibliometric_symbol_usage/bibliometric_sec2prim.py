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
    python bibliometric_search.py [database_choice] [-c/--continue]

    [database_choice]: Comma-separated list of database choices defined in the configuration file.
    -c/--continue: Continue from where the last search left off

"""

import requests
import pandas as pd
import os
import time
import json
import sys
import yaml
import argparse
import logging
import numpy as np
from datetime import datetime
from urllib.parse import quote_plus

CONFIG = "config.yml"
page_size = 1000
cursor_mark = "*"
format_type = "json"
base_url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ABSTRACT:" {} "&resultType=core&cursorMark={}&pageSize={}&format={}'

# HGNC API configuration
url_genenames = 'https://rest.genenames.org/fetch/symbol/{}'
date_format = '%Y-%m-%d'


# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def check_config_file(choices, config_file_path=CONFIG):
    """
    Check if the specified sources exist in the config file.

    Args:
        choices (list): List of source choices to check.
        config_file_path (str): Path to the config file.

    Raises:
        SystemExit: Exits with an error message if a source is not found in the config.
    """
    if not os.path.isfile(config_file_path):
        print(f"Config file '{config_file_path}' does not exist.")
        exit(1)

    with open(config_file_path, "r") as f:
        config_data = yaml.safe_load(f)

    for choice in choices:
        if choice not in config_data.keys():
            print(f"Invalid source: {choice}")
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
        with open(yaml_file_path, "r") as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)

        # Check if required keys exist in the YAML data
        required_keys = [
            "path",
            "firstColumn",
            "firstColumn_alt",
            "secondColumn",
            "secondColumn_alt",
            "keep_columns",
            "result",
        ]
        for key in required_keys:
            if key not in yaml_data[choice]:
                print(f"Required key '{key}' not found in the YAML data.")
                exit(1)
        return yaml_data

    except Exception as e:
        print(f"An error occurred while processing the config file: {str(e)}")
        exit(1)


def load_existing_results(results_path):
    """Load existing results if they exist."""
    if os.path.exists(results_path):
        try:
            with open(results_path, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Could not load existing results: {e}")

    return {
        "version": 0,
        "count": 0,
        "result": {"firstColumn_alt": [], "secondColumn_alt": []},
        "processed_symbols": [],
    }


def save_results_safely(results, results_path):
    """Save results with backup to prevent data loss."""
    # Ensure parent directory exists
    os.makedirs(os.path.dirname(results_path), exist_ok=True)

    backup_path = results_path + ".backup"
    try:
        # Save to backup first
        with open(backup_path, "w") as f:
            json.dump(results, f, indent=2)

        # If backup successful, save to main file
        with open(results_path, "w") as f:
            json.dump(results, f, indent=2)

        # Remove backup if main save successful
        if os.path.exists(backup_path):
            os.remove(backup_path)

    except Exception as e:
        logger.error(f"Error saving results: {e}")
        if os.path.exists(backup_path):
            logger.info(f"Backup saved to {backup_path}")


def search_identifier(
    symbol,
    keep_columns,
    base_url=base_url,
    cursor_mark=cursor_mark,
    page_size=page_size,
    format_type=format_type,
):
    """
    Searches EuropePMC for the provided identifier and returns a DataFrame with matching results.

    Args:
        symbol (str): Identifier to search for.
        keep_columns (list): List of columns to keep in the result DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing search results.
    """
    matches = []
    url = base_url.format(symbol, cursor_mark, page_size, format_type)
    try:
        response = requests.get(url).json()
        total = response["hitCount"]

        if "resultList" not in response:
            return pd.DataFrame()
        matches.extend(response["resultList"]["result"])
        num_requests = (total + page_size - 1) // page_size
        i = 0
        while response["nextCursorMark"] is not None:
            i += 1
            cursor_mark = response["nextCursorMark"]
            url = base_url.format(symbol, cursor_mark, page_size, format_type)
            try:
                response = requests.get(url).json()
                matches.extend(
                    [
                        i
                        for i in response["resultList"]["result"]
                        if "Human"
                        in str(response["resultList"]["result"][i]["meshHeadingList"])
                    ]
                )
            except Exception:
                break
    except Exception:
        pass

    df = pd.DataFrame(matches).fillna("NA")

    def get_nested_value(row, column):
        keys = column.split(".")
        nested_value = row
        for key in keys:
            if key in nested_value:
                nested_value = nested_value[key]
            else:
                return "NA"
        return nested_value

    for col in keep_columns:
        df[col] = df.apply(lambda row: get_nested_value(row, col), axis=1)

    keep_columns = [i if "." not in i else i.split(".")[-1] for i in keep_columns]
    for col in df.columns:
        if col not in keep_columns and "." not in col:
            df = df.drop(col, axis=1)
    return df


def print_progress_bar(iteration, total, length=50):
    """Print a simple progress bar."""
    percent = f"{100 * (iteration / float(total)):.1f}"
    filled_length = int(length * iteration // total)
    bar = "█" * filled_length + "-" * (length - filled_length)
    print(f"\rProgress |{bar}| {percent}% ({iteration}/{total})", end="", flush=True)


def sanitize_filename(filename):
    """Sanitize filename to remove invalid characters."""
    import re

    # Replace invalid characters with underscore
    return re.sub(r'[<>:"/\\|?*]', "_", str(filename))


def multiple_search(df, results_path, config, continue_search=False):
    """
    Searches for all symbols in the provided dataset and builds a DataFrame with unique hits.

    Args:
        df (pd.DataFrame): Input dataset.
        results_path (str): Path to save the results.
        config (dict): Configuration data.
        continue_search (bool): Whether to continue from previous search.
    """
    start_time = time.time()
    last_progress_time = start_time
    exclude = ["NA", "Entry Withdrawn"]
    firstColumn_alt = config["firstColumn_alt"]
    secondColumn_alt = config["secondColumn_alt"]
    # Use keep_columns as-is from config
    keep_columns = config["keep_columns"]
    logger.info(f"Using keep_columns: {keep_columns}")

    try:
        all_sym = list(set(df[firstColumn_alt]))
        if secondColumn_alt and secondColumn_alt in df.columns:
            all_sym.extend(list(set(df[secondColumn_alt])))
    except KeyError as e:
        logger.error(f"Key Error {e}: {df.keys()}")
        exit(0)

    # Remove duplicates, excluded values, and gene names of length 1
    all_sym = [
        str(sym) for sym in set(all_sym) if sym not in exclude and len(str(sym)) > 1
    ]

    # Load existing results or create new
    results = load_existing_results(results_path)

    # Get processed symbols to skip if continuing
    processed_symbols = set(results.get("processed_symbols", []))

    # Create directories using actual column names
    primary_dir = f"results/{firstColumn_alt}"
    secondary_dir = (
        f"results/{secondColumn_alt}" if secondColumn_alt else "results/secondary"
    )

    os.makedirs(primary_dir, exist_ok=True)
    os.makedirs(secondary_dir, exist_ok=True)

    if continue_search:
        # Check existing individual files
        existing_primary_files = set()
        existing_secondary_files = set()

        if os.path.exists(primary_dir):
            existing_primary_files = {
                f.replace(".json", "")
                for f in os.listdir(primary_dir)
                if f.endswith(".json")
            }

        if os.path.exists(secondary_dir):
            existing_secondary_files = {
                f.replace(".json", "")
                for f in os.listdir(secondary_dir)
                if f.endswith(".json")
            }

        # Combine all existing files (need to handle sanitized filenames)
        all_existing_files = existing_primary_files | existing_secondary_files

        # Also check for sanitized versions of symbols
        existing_symbols = set()
        for symbol in all_sym:
            safe_name = sanitize_filename(symbol)
            if safe_name in all_existing_files or symbol in all_existing_files:
                existing_symbols.add(symbol)

        total_existing_files = len(existing_primary_files) + len(
            existing_secondary_files
        )

        logger.info(f"Continue mode enabled:")
        logger.info(
            f"  - Found {len(processed_symbols)} processed symbols in main results file"
        )
        logger.info(
            f"  - Found {len(existing_primary_files)} existing files in {primary_dir}"
        )
        logger.info(
            f"  - Found {len(existing_secondary_files)} existing files in {secondary_dir}"
        )
        logger.info(f"  - Total existing individual files: {total_existing_files}")
        logger.info(
            f"  - Matching symbols from existing files: {len(existing_symbols)}"
        )

        # Combine both sources of already processed symbols
        all_processed = processed_symbols | existing_symbols

        # Debug: Show overlap analysis
        overlap = processed_symbols & existing_symbols
        only_in_main = processed_symbols - existing_symbols
        only_in_files = existing_symbols - processed_symbols

        logger.info(f"Overlap analysis:")
        logger.info(
            f"  - Symbols in both main file and individual files: {len(overlap)}"
        )
        logger.info(
            f"  - Symbols only in main file (missing individual files): {len(only_in_main)}"
        )
        logger.info(
            f"  - Symbols only as individual files (missing from main): {len(only_in_files)}"
        )

        if only_in_main:
            logger.warning(
                f"Found {len(only_in_main)} symbols in main file without individual files - will reprocess these!"
            )
        if only_in_files:
            logger.info(
                f"Found {len(only_in_files)} individual files not recorded in main file - will trust individual files"
            )

        # Only skip symbols that have individual files (the proof of complete processing)
        # Reprocess symbols that are only in main file but missing individual files
        truly_processed = (
            existing_symbols  # Only trust individual files as proof of completion
        )

        if truly_processed:
            remaining_symbols = [s for s in all_sym if s not in truly_processed]
            logger.info(
                f"Will reprocess {len(only_in_main)} symbols missing individual files"
            )
            print(
                f"Continuing search: {len(truly_processed)} truly processed (with files), {len(remaining_symbols)} remaining (including {len(only_in_main)} to reprocess)"
            )
            all_sym = remaining_symbols
        elif total_existing_files > 0:
            logger.info(
                f"No processed symbols in main file, but found {total_existing_files} individual symbol files"
            )

        # Clear processed_symbols list in results since we're reprocessing some
        if only_in_main:
            logger.info(
                f"Resetting processed_symbols list to only include symbols with individual files"
            )
            results["processed_symbols"] = list(existing_symbols)

    total = len(all_sym)
    if total == 0:
        print("No symbols to process.")
        return

    print(f"Searching {total} symbols...")
    save_frequency = 50  # Save every 50 symbols

    for i, symbol in enumerate(all_sym):
        current_time = time.time()

        # Check if we've been stuck on this symbol for too long (no progress in 10 minutes)
        if current_time - last_progress_time > 600:  # 10 minutes
            logger.warning(
                f"No progress detected for 10 minutes on symbol {symbol}. Waiting 5 minutes before continuing..."
            )
            time.sleep(300)  # 5 minutes
            last_progress_time = current_time

        try:
            logger.info(f"Processing symbol {i+1}/{total}: {symbol}")
            print_progress_bar(i + 1, total)

            search = search_identifier(symbol, keep_columns=keep_columns)

            # Update progress time after successful search
            last_progress_time = time.time()
            logger.info(f"Completed search for symbol {symbol}")

            if not search.empty:
                stype = (
                    firstColumn_alt
                    if symbol in df[firstColumn_alt].values
                    else secondColumn_alt
                )
                symbol_dir = (
                    primary_dir
                    if symbol in df[firstColumn_alt].values
                    else secondary_dir
                )

                try:
                    # Save individual symbol result as JSON file with sanitized filename
                    safe_symbol_name = sanitize_filename(symbol)
                    symbol_file = os.path.join(symbol_dir, f"{safe_symbol_name}.json")
                    symbol_data = {symbol: search.transpose().to_dict()}

                    with open(symbol_file, "w") as f:
                        json.dump(symbol_data, f, indent=2)

                    # Also add to main results
                    if stype not in results["result"]:
                        results["result"][stype] = []
                    results["result"][stype].append(symbol_data)
                    results["count"] += 1

                except Exception as e:
                    logger.error(f"Error processing results for {symbol}: {e}")

            # Add to processed symbols
            if "processed_symbols" not in results:
                results["processed_symbols"] = []
            results["processed_symbols"].append(symbol)

            # Save periodically
            if (i + 1) % save_frequency == 0:
                save_results_safely(results, results_path)
                logger.info(f"Saved progress: {i + 1}/{total} symbols processed")
                last_progress_time = time.time()  # Update progress time after save

        except Exception as e:
            logger.error(f"Unexpected error processing {symbol}: {e}")
            # Wait a bit before continuing to avoid cascading errors
            time.sleep(10)
            continue

    # Final save
    save_results_safely(results, results_path)

    elapsed_time = time.time() - start_time
    print(
        f'\nCompleted in {round(elapsed_time, 0)}s. Found {results["count"]} results.'
    )
    print(f"Individual symbol files saved to: {primary_dir} and {secondary_dir}")


def generate_overall_data(df, json_file):
    try:
        # Load JSON file
        with open(json_file, "r") as json_data:
            json_content = json.load(json_data)

        # Create a DataFrame for primary and secondary symbols
        symbol_columns = ['primarySymbol', 'secondarySymbol', 'mapping_cardinality_sec2pri']
        primary_df = df[symbol_columns].drop_duplicates()
        secondary_df = df[symbol_columns].drop_duplicates()

        # Initialize tracking variables
        overall_data = []
        seen_primaries_years = {}
        hgncs = df  # Reference to the input dataframe

        # Process primarySymbol and secondarySymbol data
        length = len(json_content.get('primarySymbol', [])) + len(json_content.get('secondarySymbol', []))
        i = -1
        for alt_data, symbol_type in [(json_content.get('primarySymbol', []), 'primary'), (json_content.get('secondarySymbol', []), 'secondary')]:
            
            for entry in alt_data:
                symbol, article_data = entry.popitem()
                symbol_df = primary_df if symbol_type == 'primary' else secondary_df
                symbol_row = symbol_df[symbol_df[symbol_type + 'Symbol'] == symbol]
                if not symbol_row.empty:
                    symbol_row = symbol_row.iloc[0]
                    for article in article_data.values():
                        accession = ""
                        primary_symbol = symbol_row['primarySymbol'] if symbol_type == 'secondary' else symbol
                        secondary_symbol = symbol if symbol_type == 'secondary' else symbol_row['secondarySymbol']
                        journal_title = article.get('journalInfo.journal.title', 'NA')
                        pubYear = article.get('pubYear', 'NA')
                        cardinality = symbol_row['mapping_cardinality_sec2pri']
                        pmid = article.get('pmid', 'NA')
                        
                        if primary_symbol not in seen_primaries_years.keys():
                            # Look for symbol change date using the HGNC API
                            url_symbol = url_genenames.format(symbol)
                            r = requests.get(url_symbol, headers={'Accept':'application/json'})
                            try:
                                accession = r.json()['response']['docs'][0]['hgnc_id']
                                # Check that this is the same accession
                                if 'primaryID' in hgncs.columns and secondary_symbol == hgncs[hgncs['primaryID'] == accession]['secondarySymbol'].values[0]: 
                                    try:
                                        date = r.json()['response']['docs'][0]['date_symbol_changed']
                                        date_object = datetime.strptime(date, date_format)
                                        year_change = date_object.year
                                    except Exception:
                                        year_change = np.nan
                            except:
                                year_change = np.nan
                            seen_primaries_years[primary_symbol] = year_change
                        else:
                            year_change = seen_primaries_years[primary_symbol]
                        
                        i += 1
                        # Check for 'NA' values and skip the loop iteration if any is found
                        if 'NA' in [primary_symbol, secondary_symbol, journal_title, pubYear, cardinality, pmid]:
                            continue
                        
                        overall_data.append({
                            'primarySymbol': primary_symbol,
                            'secondarySymbol': secondary_symbol,
                            'journal_title': journal_title,
                            'pubYear': pubYear,
                            'mapping_cardinality_sec2pri': cardinality,
                            'symbol_type': symbol_type,
                            'pmid': pmid,
                            'year_change': year_change
                        })
                        
                        if not np.isnan(year_change):
                            print(f'pmid # {i}: {pmid}, primary: {primary_symbol} --({year_change})--> secondary: {secondary_symbol}, accession: {accession}', end='\r', flush=True)
                        else:
                            print(f'pmid # {i}: {pmid}:, {primary_symbol}, no year_change available, {secondary_symbol}, accession: {accession}', end='\r')

        # Create a DataFrame from the overall_data list
        overall_df = pd.DataFrame(overall_data)
        if not overall_df.empty:
            overall_df = overall_df[overall_df['pubYear'] != np.nan]
            overall_df = overall_df[overall_df['pubYear'] != 'NA']
            overall_df['pubYear'] = pd.to_numeric(overall_df['pubYear'])
            overall_df['symbol_type'] = overall_df['symbol_type'].astype('category')
            overall_df['symbol_code'] = ['0' if i=='primary' else '1' for i in overall_df['symbol_type']]
            overall_df['symbol'] = np.where(overall_df['symbol_type'] == 'primary', overall_df['primarySymbol'], overall_df['secondarySymbol'])
            overall_df = overall_df[overall_df['primarySymbol'] != 'Entry Withdrawn']
            
            # Save to pickle and CSV
            overall_df.to_pickle('overall_df.pkl')
            overall_df.to_csv('overall_df.csv', index=False)
            logger.info("Saved overall_df to pickle and CSV files")
        
        return overall_df

        # Preserve original fallback logic for compatibility
        if "primarySymbol" not in json_content:
            raise KeyError("Missing 'primarySymbol' key in JSON file.")

        primary_symbol_data = json_content["primarySymbol"]
        original_overall_data = []

        # Process each primary symbol entry (original logic preserved)
        for primary_symbol, articles in primary_symbol_data.items():
            symbol_row = df[df["primarySymbol"] == primary_symbol]

            if symbol_row.empty:
                continue

            symbol_row = symbol_row.iloc[0]
            secondary_symbol = symbol_row["secondarySymbol"]
            cardinality = symbol_row["mapping_cardinality_sec2pri"]

            for article_id, article_data in articles.items():
                journal_title = article_data.get("journalInfo.journal.title", "NA")
                pub_year = article_data.get("pubYear", "NA")
                pmid = article_data.get("pmid", "NA")

                # Skip invalid entries
                if "NA" in [
                    primary_symbol,
                    secondary_symbol,
                    journal_title,
                    pub_year,
                    cardinality,
                    pmid,
                ]:
                    continue

                original_overall_data.append(
                    {
                        "primarySymbol": primary_symbol,
                        "secondarySymbol": secondary_symbol,
                        "journal_title": journal_title,
                        "pubYear": pub_year,
                        "mapping_cardinality_sec2pri": cardinality,
                        "symbol_type": "primary",
                        "pmid": pmid,
                    }
                )

        # Create DataFrame from collected data (original structure preserved)
        if not overall_data:  # If new processing didn't work, fall back to original
            overall_df = pd.DataFrame(
                original_overall_data,
                columns=[
                    "primarySymbol",
                    "secondarySymbol",
                    "journal_title",
                    "pubYear",
                    "mapping_cardinality_sec2pri",
                    "symbol_type",
                    "pmid",
                ],
            )
        return overall_df

    except FileNotFoundError:
        print(f"File not found: {json_file}")
        return pd.DataFrame(
            columns=[
                "primarySymbol",
                "secondarySymbol",
                "journal_title",
                "pubYear",
                "mapping_cardinality_sec2pri",
                "symbol_type",
                "pmid",
            ]
        )
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON file '{json_file}': {e}")
        return pd.DataFrame(
            columns=[
                "primarySymbol",
                "secondarySymbol",
                "journal_title",
                "pubYear",
                "mapping_cardinality_sec2pri",
                "symbol_type",
                "pmid",
            ]
        )
    except KeyError as e:
        print(f"Key error: {e}")
        return pd.DataFrame(
            columns=[
                "primarySymbol",
                "secondarySymbol",
                "journal_title",
                "pubYear",
                "mapping_cardinality_sec2pri",
                "symbol_type",
                "pmid",
            ]
        )
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return pd.DataFrame(
            columns=[
                "primarySymbol",
                "secondarySymbol",
                "journal_title",
                "pubYear",
                "mapping_cardinality_sec2pri",
                "symbol_type",
                "pmid",
            ]
        )


def main():
    print("bibliometric_sec2prim.py\n_______________\n")
    parser = argparse.ArgumentParser(
        description="Check if a database key exists in the config file."
    )
    parser.add_argument("database", type=str, help="Database key to check")
    parser.add_argument(
        "-c",
        "--continue",
        action="store_true",
        dest="continue_search",
        help="Continue from where the last search left off",
    )
    args = parser.parse_args()

    sources = args.database.split(",")

    check_config_file(sources)

    for choice in sources:
        config_choice = read_config_file(choice=choice, yaml_file_path=CONFIG)[choice]
        data = config_choice["path"]
        results = config_choice["result"]

        try:
            df = pd.read_csv(data, sep="\t")
            print(f"Processing dataset: {choice}")
            multiple_search(
                df=df,
                results_path=results,
                config=config_choice,
                continue_search=args.continue_search,
            )
            print(f"Results saved to {results}\n")
        except Exception as e:
            logger.error(f"Error processing dataset {choice}: {e}")
            continue


if __name__ == "__main__":
    main()
