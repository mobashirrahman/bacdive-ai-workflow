import pickle
import os
import argparse
import numpy as np
import csv
import logging
import sys
import datetime
import traceback
import json
from pathlib import Path

def setup_logger(log_level=logging.INFO, log_file=None):
    """
    Set up and configure logger with console and file handlers
    
    Parameters:
    -----------
    log_level : int
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    log_file : str, optional
        Path to log file. If None, logs will only be written to console
        
    Returns:
    --------
    logger : logging.Logger
        Configured logger instance
    """
    # Create logger
    logger = logging.getLogger("BacDive-AI")
    logger.setLevel(log_level)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Create console handler
    console_handler = logging.StreamHandler(sys.stderr)  # Changed to stderr so stdout can be used for JSON output
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if log_file is provided
    if log_file:
        # Ensure log directory exists
        log_dir = os.path.dirname(log_file)
        if log_dir:
            Path(log_dir).mkdir(parents=True, exist_ok=True)
            
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def main():
    # Set up logging
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = os.path.join(os.path.dirname(__file__), "logs", "predict")
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    log_file = os.path.join(log_dir, f"prediction_{timestamp}.log")
    
    logger = setup_logger(log_level=logging.INFO, log_file=log_file)
    logger.info("Starting BacDive-AI prediction")
    
    try:
        # set up possible arguments
        logger.debug("Setting up argument parser")
        parser = argparse.ArgumentParser(description="BacDive-AI package",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("trait", help="Trait to predict, all to see all models")
        parser.add_argument("file", help="Interproscan file or file with list of Pfams")
        parser.add_argument("--debug", action="store_true", help="Enable debug logging")
        parser.add_argument("--log-file", help="Custom log file path")
        parser.add_argument("--text-output", action="store_true", help="Output human-readable text instead of JSON")

        # load configuration from command line arguments
        logger.debug("Parsing command line arguments")
        args = parser.parse_args()
        config = vars(args)
        
        # Update logging level if debug flag is set
        if config.get('debug'):
            logger.setLevel(logging.DEBUG)
            logger.debug("Debug logging enabled")
        
        # Use custom log file if provided
        if config.get('log_file'):
            for handler in logger.handlers[:]:
                if isinstance(handler, logging.FileHandler):
                    logger.removeHandler(handler)
                    
            custom_log_path = config.get('log_file')
            file_handler = logging.FileHandler(custom_log_path)
            formatter = logging.Formatter(
                '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            logger.info(f"Using custom log file: {custom_log_path}")

        # define constants
        logger.debug("Setting up constants")
        EVALUE = 1e-20
        MODELPATH = os.path.dirname(__file__)+"/models/"
        logger.info(f"Model directory: {MODELPATH}")

        # set up all available models
        logger.debug("Configuring available trait models")
        traits = {
            "acidophile": 'Acidophilic',
            "gram-positive": 'Gram-positive',
            "spore-forming": 'Spore-forming',
            "aerobic": 'Aerobic',
            "anaerobic": 'Anaerobic',
            "thermophile": 'Thermophilic',
            "psychrophile": 'Psychrophilic',
            "motile2+": 'Flagellated motility',
        }
        logger.debug(f"Available traits: {', '.join(traits.keys())}")

        # check if selected trait is supported
        trait = config.get('trait')
        logger.info(f"Selected trait: {trait}")
        
        if trait != 'all' and trait not in traits:
            logger.error(f"Trait '{trait}' is not supported")
            logger.info(f"Supported traits are: {', '.join(list(traits.keys()))}")
            print(trait, 'is not supported')
            print('Supported traits are:', ', '.join(list(traits.keys())))
            return 1

        # load pfam data set
        filename = config.get('file')
        logger.info(f"Loading Pfams from file: {filename}")
        
        try:
            pfams = []
            with open(filename, 'r', encoding='utf-8') as f:
                logger.debug(f"Parsing file {filename} with E-value threshold {EVALUE}")
                csv_file = csv.reader(f, delimiter='\t')
                line_count = 0
                included_count = 0
                
                for line in csv_file:
                    line_count += 1
                    try:
                        evalue = float(line[8])
                        if evalue > EVALUE:
                            continue
                        pfams.append(line[4])
                        included_count += 1
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing line {line_count}: {e}")
                        logger.debug(f"Problematic line content: {line}")
                
                logger.info(f"Processed {line_count} lines, included {included_count} Pfams below E-value threshold")
            
            pfams = set(pfams)
            logger.info(f"Found {len(pfams)} unique Pfams for analysis")
            
            # if specific trait has been selected, reduce to this
            if trait != 'all': 
                logger.debug(f"Using only {trait} model")
                traits = {trait: traits[trait]}
            else:
                logger.info("Running prediction for all available traits")

            # go through all remaining traits
            results = {}
            for trait in traits:
                label = traits[trait]
                logger.info(f"Processing trait: {label} ({trait})")

                # load model
                model_path = MODELPATH + trait + "_data.p"
                logger.debug(f"Loading model from {model_path}")
                
                try:
                    dump = pickle.load(open(model_path, "rb"))
                    clf = dump.get('model')
                    categories = dump.get("categories")
                    strains = dump.get('strains')
                    logger.debug(f"Model loaded with {len(categories)} categories and {len(strains) if strains else 0} strains")
                    
                    # transform sample
                    logger.debug("Transforming input data")
                    lst = dict(zip(*np.unique(list(pfams), return_counts=True)))
                    X = [lst.get(k) if k in lst else 0 for k in categories]
                    logger.debug(f"Feature vector created with {len(X)} dimensions")
                    
                    # predict with probability
                    logger.debug("Making prediction")
                    proba = clf.predict_proba([X])
                    y = proba.argmax(axis=1)
                    y = y[0]
                    proba = proba[0]
                    true_index = list(clf.classes_).index(y)
                    confidence = proba[true_index] * 100
                    
                    # store result
                    text_result = f"{label}: {bool(y)} ({round(confidence, 2)}%)"
                    results[label] = {"value": bool(y), "confidence": round(confidence, 2)}
                    logger.info(f"Prediction result: {text_result}")
                    
                    # Print to stderr for human-readable output during processing
                    if config.get('text_output'):
                        print(text_result)
                    
                except FileNotFoundError:
                    logger.error(f"Model file not found: {model_path}")
                    sys.stderr.write(f"Error: Model file not found for trait '{trait}'\n")
                except Exception as e:
                    logger.error(f"Error processing trait {trait}: {str(e)}")
                    logger.debug(f"Exception details: {traceback.format_exc()}")
                    sys.stderr.write(f"Error processing trait '{trait}': {str(e)}\n")
            
            logger.info("Prediction completed successfully")
            
            # Output the final results as JSON to stdout
            if not config.get('text_output'):
                # Format results for JSON output
                output = {
                    "genome_file": os.path.basename(filename),
                    "predictions": {}
                }
                
                # Add all prediction results
                for trait_name, result in results.items():
                    output["predictions"][trait_name] = {
                        "prediction": result["value"],
                        "confidence": result["confidence"]
                    }
                
                # Output JSON to stdout
                print(json.dumps(output, indent=2))
            else:
                # Text output was already printed trait-by-trait above
                pass
                
            return 0
            
        except FileNotFoundError:
            logger.error(f"Input file not found: {filename}")
            sys.stderr.write(f"Error: File not found: {filename}\n")
            return 1
        except Exception as e:
            logger.error(f"Error loading Pfams: {str(e)}")
            logger.debug(f"Exception details: {traceback.format_exc()}")
            sys.stderr.write(f"Error loading Pfams: {str(e)}\n")
            return 1
            
    except Exception as e:
        logger.critical(f"Unhandled exception: {str(e)}")
        logger.debug(f"Exception details: {traceback.format_exc()}")
        sys.stderr.write(f"Critical error: {str(e)}\n")
        return 1
    finally:
        logger.info("BacDive-AI prediction finished")

if __name__ == "__main__":
    sys.exit(main())
