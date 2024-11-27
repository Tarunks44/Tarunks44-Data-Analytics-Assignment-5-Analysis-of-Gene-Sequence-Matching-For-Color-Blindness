import numpy as np
import time
import math
from typing import List, Tuple, Dict
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor
from itertools import islice

class GeneAnalysisError(Exception):
    """Custom exception for gene analysis errors"""
    pass

def validate_input_files(required_files: List[str]) -> None:
    """Validate presence and readability of input files"""
    for file_path in required_files:
        if not Path(file_path).is_file():
            raise GeneAnalysisError(f"Required file not found: {file_path}")
        try:
            with open(file_path, 'r') as f:
                f.readline()
        except Exception as e:
            raise GeneAnalysisError(f"Error reading file {file_path}: {str(e)}")

def calculate_rank_array(LastCol: str, counts: dict) -> np.ndarray:
    """Optimized rank array calculation using numpy"""
    try:
        # Convert to uint8 for memory efficiency
        last_col_array = np.frombuffer(LastCol.encode(), dtype=np.uint8)
        rank = np.zeros(len(LastCol), dtype=np.int32)
        
        # Use lookup table for faster conversion
        lookup = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3, ord('$'): 4}
        
        # Calculate cumulative counts
        cumsum = np.zeros(5, dtype=np.int32)
        cumsum[1:] = np.cumsum([counts['A'], counts['C'], counts['G'], counts['T']])
        
        # Vectorized rank calculation
        for i, (char, start) in enumerate(lookup.items()):
            mask = last_col_array == char
            count = np.sum(mask)
            if count > 0:
                rank[mask] = np.arange(count, dtype=np.int32) + cumsum[i]
        
        return rank
    except Exception as e:
        raise GeneAnalysisError(f"Error in rank array calculation: {str(e)}")

def get_exact_match_idcs(read: str, FirstCol: str, LastCol: str, Rank: np.ndarray, 
                        Map: np.ndarray, lookup_cache: Dict = {}) -> List[int]:
    """Optimized exact match finding with caching"""
    try:
        read = read.replace("N", "A")
        
        # Check cache
        cache_key = read[-1]
        if cache_key in lookup_cache:
            batch_start, batch_end = lookup_cache[cache_key]
        else:
            batch_start, batch_end = FirstCol.index(read[-1]), FirstCol.rindex(read[-1])
            lookup_cache[cache_key] = (batch_start, batch_end)
        
        for i in range(2, len(read)+1):
            last = LastCol[batch_start:batch_end+1]
            try:
                sub_batch_start = Rank[batch_start + last.index(read[-i])]
                sub_batch_end = Rank[batch_start + last.rindex(read[-i])]
            except ValueError:
                return []
            batch_start, batch_end = sub_batch_start, sub_batch_end
        
        return Map[batch_start:batch_end+1].tolist()
    except Exception as e:
        raise GeneAnalysisError(f"Error in exact match finding: {str(e)}")

def match_read_to_loc_optimized(read: str, FirstCol: str, LastCol: str, Rank: np.ndarray, 
                              Map: np.ndarray, RefSeq: str, lookup_cache: Dict) -> List[int]:
    """Optimized read matching with vectorized operations"""
    try:
        l = len(read)
        third = l // 3
        
        # Get matches for each third using numpy operations
        parts = [
            (read[:third], 0),
            (read[third:2*third], third),
            (read[2*third:], 2*third)
        ]
        
        all_positions = []
        for part_read, offset in parts:
            positions = np.array(get_exact_match_idcs(part_read, FirstCol, LastCol, Rank, Map, lookup_cache))
            if len(positions) > 0:
                all_positions.append(positions - offset)
        
        if not all_positions:
            return []
        
        # Efficient position combination and filtering
        positions = np.unique(np.concatenate(all_positions))
        positions = positions[(positions >= 0) & (positions + l <= len(RefSeq))]
        
        # Vectorized mismatch calculation
        read_array = np.frombuffer(read.encode(), dtype=np.uint8)
        valid_positions = []
        
        for pos in positions:
            ref_slice = np.frombuffer(RefSeq[pos:pos+l].encode(), dtype=np.uint8)
            mismatches = np.sum(read_array != ref_slice)
            if mismatches <= 2:
                valid_positions.append(int(pos))
        
        return valid_positions
    except Exception as e:
        raise GeneAnalysisError(f"Error in read matching: {str(e)}")

def count_exons_optimized(positions: List[int], red_pos: np.ndarray, 
                         green_pos: np.ndarray) -> np.ndarray:
    """Optimized exon counting using vectorized operations"""
    try:
        counts = np.zeros(12)
        if not positions:
            return counts
        
        positions = np.array(positions)[:, np.newaxis]
        
        # Vectorized range checks
        red_matches = (positions >= red_pos[:, 0]) & (positions <= red_pos[:, 1])
        green_matches = (positions >= green_pos[:, 0]) & (positions <= green_pos[:, 1])
        
        # Calculate counts
        for i in range(6):
            red_match = np.any(red_matches[:, i])
            green_match = np.any(green_matches[:, i])
            
            if red_match and green_match:
                counts[i] = counts[i+6] = 0.5
            elif red_match:
                counts[i] = 1.0
            elif green_match:
                counts[i+6] = 1.0
        
        return counts
    except Exception as e:
        raise GeneAnalysisError(f"Error in exon counting: {str(e)}")

def process_read_batch(batch: List[str], FirstCol: str, LastCol: str, Rank: np.ndarray,
                      Map: np.ndarray, RefSeq: str, red_pos: np.ndarray, 
                      green_pos: np.ndarray, lookup_cache: Dict) -> Tuple[np.ndarray, np.ndarray, int]:
    """Process a batch of reads in parallel"""
    exact_counts = np.zeros(12)
    mismatch_counts = np.zeros(12)
    total = 0
    
    for read in batch:
        try:
            positions = match_read_to_loc_optimized(read, FirstCol, LastCol, Rank, Map, RefSeq, lookup_cache)
            total += len(positions)
            
            if positions:
                counts = count_exons_optimized(positions, red_pos, green_pos)
                mismatches = np.sum(np.frombuffer(RefSeq[positions[0]:positions[0]+len(read)].encode(), dtype=np.uint8) 
                                  != np.frombuffer(read.encode(), dtype=np.uint8))
                
                if mismatches == 0:
                    exact_counts += counts
                mismatch_counts += counts
        except Exception as e:
            print(f"Warning: Error processing read: {str(e)}", file=sys.stderr)
            continue
    
    return exact_counts, mismatch_counts, total

def main():
    try:
        # Validate input files
        required_files = ["chrX_last_col.txt", "chrX.fa", "reads", "chrX_map.txt"]
        validate_input_files(required_files)
        
        # Gene positions
        RedExonPos = np.array([
            [149249757, 149249868],
            [149256127, 149256423],
            [149258412, 149258580],
            [149260048, 149260213],
            [149261768, 149262007],
            [149264290, 149264400]
        ])
        
        GreenExonPos = np.array([
            [149288166, 149288277],
            [149293258, 149293554],
            [149295542, 149295710],
            [149297178, 149297343],
            [149298898, 149299137],
            [149301420, 149301530]
        ])
        
        # Load data with error handling
        start_time = time.time()
        print("Loading data...")
        
        try:
            LastCol = ''.join(np.loadtxt("chrX_last_col.txt", dtype=str))
            RefSeq = ''.join(np.loadtxt("chrX.fa", dtype=str)[1:])
            Reads = np.loadtxt("reads", dtype=str)
            Map = np.loadtxt("chrX_map.txt", dtype=int)
        except Exception as e:
            raise GeneAnalysisError(f"Error loading input files: {str(e)}")
        
        # Initialize lookup cache
        lookup_cache = {}
        
        # Calculate counts and create FirstCol
        counts = {'A': LastCol.count('A'), 'C': LastCol.count('C'),
                 'G': LastCol.count('G'), 'T': LastCol.count('T'), '$': LastCol.count('$')}
        FirstCol = 'A'*counts['A'] + 'C'*counts['C'] + 'G'*counts['G'] + 'T'*counts['T'] + '$'
        
        # Calculate Rank array
        Rank = calculate_rank_array(LastCol, counts)
        
        # Print initial stats
        print(f"Length of the map data: {len(LastCol)}")
        print(counts)
        print(f"Found $ at position: {LastCol.index('$')}")
        
        # Process test string
        test_string = "GAGGACAGCACCCAGTCCAGCATCTTCACCTACACCAACAGCAACTCCACCAGAGGTGAGCCAGCAGGCCCGTGGAGGCTGGGTGGCTGCACTGGGGGCCA"
        test_positions = match_read_to_loc_optimized(test_string, FirstCol, LastCol, Rank, Map, RefSeq, lookup_cache)
        
        if test_positions:
            print(f"Map value for test string: {test_positions[0]}")
            matches = get_exact_match_idcs(test_string, FirstCol, LastCol, Rank, Map, lookup_cache)
            if matches:
                print(f"Band first match index: {matches[0]}")
                print(f"Band last match index: {matches[-1]}")
        
        # Process reads in parallel batches
        print("\nProcessing reads...")
        batch_size = 1000
        exact_red = np.zeros(6)
        exact_green = np.zeros(6)
        mismatch_red = np.zeros(6)
        mismatch_green = np.zeros(6)
        total_matches = 0
        
        reads_to_process = Reads[2936000:2948000]
        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = []
            for i in range(0, len(reads_to_process), batch_size):
                batch = reads_to_process[i:i+batch_size]
                future = executor.submit(process_read_batch, batch, FirstCol, LastCol, 
                                      Rank, Map, RefSeq, RedExonPos, GreenExonPos, lookup_cache)
                futures.append(future)
            
            for future in futures:
                try:
                    exact_counts, mismatch_counts, batch_total = future.result()
                    exact_red += exact_counts[:6]
                    exact_green += exact_counts[6:]
                    mismatch_red += mismatch_counts[:6]
                    mismatch_green += mismatch_counts[6:]
                    total_matches += batch_total
                except Exception as e:
                    print(f"Warning: Error processing batch: {str(e)}", file=sys.stderr)
        
        # Print results
        print(f"Counts for an exact match of red gene= {exact_red.tolist()}")
        print(f"Counts for an exact match of green gene = {exact_green.tolist()}")
        print(f"Count of red exons considering up to two mismatches = {mismatch_red.tolist()}")
        print(f"Count of green exons considering up to two mismatches = {mismatch_green.tolist()}")
        
        # Calculate and print log probabilities
        configs = [
            ([0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5]),
            ([1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0]),
            ([0.33, 0.33, 1.0, 1.0], [0.67, 0.67, 0.0, 0.0]),
            ([0.33, 0.33, 0.33, 1.0], [0.67, 0.67, 0.67, 0.0])
        ]
        
        print("\nLog probabilities for configurations:")
        config_names = [
            "(50%, 50%, 50%, 50%)",
            "(100%, 100%, 0%, 0%)",
            "(33%, 33%, 100%, 100%)",
            "(33%, 33%, 33%, 100%)"
        ]
        
        for (red_prob, green_prob), config_name in zip(configs, config_names):
            try:
                log_prob = 0
                for i in range(4):
                    if mismatch_red[i+1] > 0 and red_prob[i] > 0:
                        log_prob += mismatch_red[i+1] * np.log(red_prob[i])
                    if mismatch_green[i+1] > 0 and green_prob[i] > 0:
                        log_prob += mismatch_green[i+1] * np.log(green_prob[i])
                
                print(f"Config {config_name}: Log Probability {log_prob:.1f} & Probability: {np.exp(log_prob):.6e}")
            except Exception as e:
                print(f"Warning: Error calculating probability for {config_name}: {str(e)}", 
                      file=sys.stderr)
        
        print(f"\nTime taken :: {time.time() - start_time:.6f}")
        print(f"Total matches :: {total_matches}")
        
    except GeneAnalysisError as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()