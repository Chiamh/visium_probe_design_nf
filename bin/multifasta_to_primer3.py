#!/usr/bin/env python
import primer3
from Bio import SeqIO
import argparse
import sys
import os

#The delta G balue is calculated at 37 degrees.
def calculate_hairpin_and_homodimer_tm(
    sequence,
    mv_conc,
    dv_conc,
    dntp_conc,
    dna_conc,
    dmso_conc=0.0,
    dmso_fact=0.6,
    formamide_conc=0.0,
    annealing_temp_c=50.0,
    temp_c=37.0
):
    """
    Calculates the self-hairpin Tm and dG for a single sequence.
    """

    # Create the ThermoAnalysis object with user-specified concentrations
    tm_analyser = primer3.thermoanalysis.ThermoAnalysis(
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        dmso_conc=dmso_conc,
        dmso_fact=dmso_fact,
        formamide_conc=formamide_conc,
        annealing_temp_c=annealing_temp_c,
        temp_c=temp_c
    )

    
    
    # Create the ThermoAnalysis object for hairpins and homodimers with user-specified concentrations
    secondary_analyser = primer3.thermoanalysis.ThermoAnalysis(
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        temp_c=temp_c
    )

    try:
        oligo_tm = tm_analyser.calc_tm(str(sequence)) #a float
        hairpin_analysis = secondary_analyser.calc_hairpin(str(sequence))
        homodimer_analysis = secondary_analyser.calc_homodimer(str(sequence))
        #returns a five element tuple
        return oligo_tm, hairpin_analysis.tm, hairpin_analysis.dg, homodimer_analysis.tm, homodimer_analysis.dg
    except Exception as e:
        # Handle calculation errors (e.g., sequence too short, invalid characters)
        return "ERROR", f"Primer3 Error: {e}"
		
		
def process_fasta(fasta_file, args):
    """
    Reads a multi-fasta file and generates a list of results.
    """
    results = []
    
    # Define the TSV header
    header = ["Sequence_ID", "Sequence_Length", "Sequence", "MeltTemp", "Hairpin_Tm_C", "Hairpin_dG_cal_per_mol", "Homodimer_Tm_C", "Homodimer_dG_cal_per_mol"]
    results.append("\t".join(header))

    # Print the conditions to stderr for context/logging
    sys.stderr.write(f"\n--- Calculating Tm and deltaG under conditions, at 37 degrees Celsius ---\n")
    sys.stderr.write(f"Monovalent Cations (M+): {args.mv} mM\n")
    sys.stderr.write(f"Divalent Cations (Mg2+): {args.dv} mM\n")
    sys.stderr.write(f"Oligo Concentration (DNA): {args.dna} nM\n")
    sys.stderr.write(f"----------------------------------------\n")

    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id
            sequence = str(record.seq).upper()
            seq_len = len(sequence)
            
            # Pass all arguments to the hairpin calculation function
            oligo_tm, hairpin_tm, hairpin_dg, homodimer_tm, homodimer_dg = calculate_hairpin_and_homodimer_tm(
                sequence,
                mv_conc=args.mv,
                dv_conc=args.dv,
                dntp_conc=args.dntp,
                dna_conc=args.dna
            )
			
            # Format results for the TSV output
            row = [
                seq_id,
                str(seq_len),
                sequence,
                #The value is formatted to two decimal places if it's a number. If it was an "ERROR" string, the string is kept.
                f"{oligo_tm:.2f}" if isinstance(oligo_tm, float) else oligo_tm,
                f"{hairpin_tm:.2f}" if isinstance(hairpin_tm, float) else hairpin_tm,
                f"{hairpin_dg:.2f}" if isinstance(hairpin_dg, float) else hairpin_dg,
                f"{homodimer_tm:.2f}" if isinstance(homodimer_tm, float) else homodimer_tm,
                f"{homodimer_dg:.2f}" if isinstance(homodimer_dg, float) else homodimer_dg
            ]
            results.append("\t".join(row))

    except FileNotFoundError:
        sys.stderr.write(f"Error: Input file '{fasta_file}' not found.\n")
        return None
    except Exception as e:
        sys.stderr.write(f"An unexpected error occurred: {e}\n")
        return None
        
    return "\n".join(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate self-hairpin and homodimer Tm for sequences in a multi-FASTA file using primer3-py."
    )
    
    # Required argument: Input FASTA file
    parser.add_argument("fasta_file", help="Path to the input multi-FASTA file.")
    
    # Optional arguments for hybridization and concentration settings (with typical defaults)
    parser.add_argument(
        "--mv", type=float, default=165.0, 
        help="Monovalent cation (e.g., Na+, K+) concentration in mM (Default: 165.0 mM for SSC buffer)."
    )
    parser.add_argument(
        "--dv", type=float, default=0, 
        help="Divalent cation (e.g., Mg2+) concentration in mM (Default: 0 for SSC buffer)."
    )
    parser.add_argument(
        "--dntp", type=float, default=0, 
        help="dNTP concentration in mM (affects free Mg2+ calculation) (Default: 0)."
    )
    parser.add_argument(
        "--dna", type=float, default=2.4, 
        help="Oligonucleotide concentration in nM (Default: 2.4 nM per probe for 10X Visium probe mix)."
    )

    args = parser.parse_args()
    
    tsv_output = process_fasta(args.fasta_file, args)

    if tsv_output:
        print(tsv_output)
		
		
#python multifasta_to_primer3.py oligos.fasta --mv 75.0 --dv 3.0 > custom_results.tsv
