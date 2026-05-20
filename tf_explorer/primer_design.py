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

def check_primer_safety(fp_seq, rp_seq):
    """
    Calculates secondary structures (hairpins, homodimers, heterodimers) for a primer pair.
    Returns safety status and thermodynamics.
    """
    try:
        # Hairpin calculations
        hp_f = primer3.bindings.calc_hairpin(fp_seq)
        hp_r = primer3.bindings.calc_hairpin(rp_seq)
        
        # Homodimer calculations
        hd_f = primer3.bindings.calc_homodimer(fp_seq)
        hd_r = primer3.bindings.calc_homodimer(rp_seq)
        
        # Heterodimer calculations
        het = primer3.bindings.calc_heterodimer(fp_seq, rp_seq)
        
        # Convert dg to kcal/mol
        dg_hp_f = hp_f.dg / 1000.0 if hp_f.structure_found else 0.0
        dg_hp_r = hp_r.dg / 1000.0 if hp_r.structure_found else 0.0
        dg_hd_f = hd_f.dg / 1000.0 if hd_f.structure_found else 0.0
        dg_hd_r = hd_r.dg / 1000.0 if hd_r.structure_found else 0.0
        dg_het = het.dg / 1000.0 if het.structure_found else 0.0
        
        # Check risks
        # Hairpin risk: Tm >= 40.0 C
        has_hp_f = hp_f.structure_found and hp_f.tm >= 40.0
        has_hp_r = hp_r.structure_found and hp_r.tm >= 40.0
        
        # Dimer risk: dg < -8.0 kcal/mol
        has_hd_f = hd_f.structure_found and dg_hd_f < -8.0
        has_hd_r = hd_r.structure_found and dg_hd_r < -8.0
        has_het = het.structure_found and dg_het < -8.0
        
        if has_hp_f or has_hp_r:
            safety = "🔴 High Hairpin Risk"
        elif has_hd_f or has_hd_r or has_het:
            safety = "⚠️ Dimer Risk"
        else:
            safety = "🟢 Safe"
            
        return {
            "fp_hairpin_tm": hp_f.tm if hp_f.structure_found else 0.0,
            "fp_hairpin_dg": dg_hp_f,
            "rp_hairpin_tm": hp_r.tm if hp_r.structure_found else 0.0,
            "rp_hairpin_dg": dg_hp_r,
            "fp_homodimer_tm": hd_f.tm if hd_f.structure_found else 0.0,
            "fp_homodimer_dg": dg_hd_f,
            "rp_homodimer_tm": hd_r.tm if hd_r.structure_found else 0.0,
            "rp_homodimer_dg": dg_hd_r,
            "heterodimer_tm": het.tm if het.structure_found else 0.0,
            "heterodimer_dg": dg_het,
            "safety_status": safety,
            "has_hp_f": has_hp_f,
            "has_hp_r": has_hp_r,
            "has_hd_f": has_hd_f,
            "has_hd_r": has_hd_r,
            "has_het": has_het
        }
    except Exception as e:
        return {
            "safety_status": "⚠️ Error evaluating safety",
            "error": str(e)
        }
