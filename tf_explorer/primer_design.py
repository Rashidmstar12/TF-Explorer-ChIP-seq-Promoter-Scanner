import primer3

def design_primers(sequence, product_size_min=70, product_size_max=300, primer_opt_size=20, primer_min_tm=57.0, primer_max_tm=63.0, target_junction=None, included_region=None):
    """
    Generates primers using primer3-py.
    """
    seq_args = {
        'SEQUENCE_ID': 'example',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': included_region if included_region else [0, len(sequence)]
    }
    
    # If a target junction is provided (position in cDNA), we force primers to overlap it.
    if target_junction:
        # Target 20bp around the junction to ensure product spans it
        # Ensure we don't go out of bounds
        start = max(1, target_junction - 10)
        length = 20
        if start + length > len(sequence):
            length = len(sequence) - start
            
        seq_args['SEQUENCE_TARGET'] = [start, length] 

    global_args = {
        'PRIMER_OPT_SIZE': primer_opt_size,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': primer_min_tm,
        'PRIMER_MAX_TM': primer_max_tm,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[product_size_min, product_size_max]],
        'PRIMER_EXPLAIN_FLAG': 1,
    }

    try:
        results = primer3.bindings.designPrimers(seq_args, global_args)
        return results
    except Exception as e:
        return {"error": str(e)}
