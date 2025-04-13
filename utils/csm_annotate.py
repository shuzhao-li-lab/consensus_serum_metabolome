import os
from mining import *
from csm_align import *
from khipu_custom import *


KCD_formula_coordinate = build_KCD_from_formula_coordinate(formula_coordinate)

def annotate_asari_table_by_csm(
                        asari_table, 
                        KCD_formula_coordinate, 
                        isotope_search_patterns,
                        adduct_search_patterns,
                        extended_adducts,
                        bmDict,                 # from cpd DB
                        reflib,                 # e.g. lib_rppos, 
                        neuDict,                # from neutral mass records
                        primary_ions_ordered,   # ion mode dependent
                        mode='pos', 
                        method='rppos',         # not used
                        mz_tolerance_ppm=5, 
                        rt_tolerance=2,
                        quality_filter=False,
                        ):
    '''
    Takes an asari feature table and produce CSM annotation. 
    Output anno dict, header, list of primary_fids, stats. 
    A separate step writes recommended anno table and full anno table. 
    
    The relationships among features need to be established before annotation.
    Step 1. khipu pre-annotation.
    Step 2. mapping unique epds to CSM.
    Step 3. Complement annotation by L4 in asari (can be replaced later).
    
    Future add: if method not in methods - find closest match
    '''
    list_khipus, list_singletons, list_features = run_khipu_insitu(asari_table,
                                KCD_formula_coordinate, 
                                isotope_search_patterns,
                                adduct_search_patterns,
                                extended_adducts,
                                mode = mode,
                                mz_tolerance_ppm=mz_tolerance_ppm, 
                                rt_tolerance=rt_tolerance,
                                separate_singletons=True,
                                )
    # list_singletons from run_khipu_insitu is list of IDs
    list_singletons = [f for f in list_features if f['id'] in list_singletons]
    
    return annotate_lists_by_csm(
                        list_khipus, 
                        list_singletons,
                        list_features,
                        bmDict, 
                        reflib, 
                        neuDict, 
                        primary_ions_ordered,
                        mz_tolerance_ppm,
                        quality_filter=quality_filter
                        )

def annotate_epdJSON_by_csm(
                        json_epd_file, 
                        bmDict,                 # from cpd DB
                        reflib,                 # e.g. lib_rppos, 
                        neuDict,                # from neutral mass records
                        primary_ions_ordered,   # ion mode dependent
                        mz_tolerance_ppm=5, 
                        quality_filter=False,
                        ):
    '''
    Similar to annotate_asari_table_by_csm but takes khipu result JSON as input. 
    Output CSM annotation dict, header, list of primary_fids, stats. 
    '''
    list_khipus, list_singletons = get_all_features_from_epdlist(
        json.load(open(json_epd_file))
    )
    list_features = list_khipus + list_singletons
    return annotate_lists_by_csm(
                        list_khipus, 
                        list_singletons,
                        list_features,
                        bmDict, 
                        reflib, 
                        neuDict, 
                        primary_ions_ordered,
                        mz_tolerance_ppm,
                        quality_filter=quality_filter
                        )

def annotate_lists_by_csm(
                        list_khipus, 
                        list_singletons, 
                        list_features, 
                        bmDict, 
                        reflib,
                        neuDict, 
                        primary_ions_ordered,
                        mz_tolerance_ppm,
                        quality_filter=False
                        ):
    '''
    returns anno dict, header, list of primary_fids, stats
    
    '''
    clean_dict_match, featureDict, primary_fids = csmf_annotate_from_separate_lists(
        list_khipus, list_singletons, list_features, reflib, 
        primary_ions_ordered=primary_ions_ordered, 
        mz_tolerance_ppm=mz_tolerance_ppm,
    )
    anno, header = export_csmf_annotation(
        clean_dict_match, featureDict, bmDict=bmDict,
                                        reflib=reflib, 
                                        neuDict=neuDict,
                                        )
    stats = {
        'Number_features': len(featureDict),
        'Number_khipus': len(list_khipus),
        'Number_singletons': len(list_singletons),
        'Compounds_in_DB': len(bmDict),
        'Size_CSM_ref_lib': len(reflib),
        'Number_annotated_khipus': len([k for k in anno if k in primary_fids]),
        'Number_annotated_features': len(anno),
    }
    # Filter. Do not use if output recommended and full anno files separately
    if quality_filter:
            new = {}
            for k,v in anno.items():
                    if v['ion_csm'] in primary_ions:
                            new[k] = v
            anno = new
            
    return anno, header, primary_fids, stats

def write_csm_anno_2files(anno, header, extra_anno_dict, primary_fids,
                          outdir='', outstr='test_csmf_anno.tsv', sep='\t'):
    '''
    Write 2 tables of annotation, following columns in header, for full and recommended_.
    anno as dict of fid: dict
    featureDict : all features
    extra_anno_dict : extracted from asari output, default HMDB level 4 annotation.
    
    Features can match other DBs but not CSM. That's okay. 
    '''
    header = header + ['L4_guess']
    for k,v in anno.items():
            v['L4_guess'] = ''
            if k in extra_anno_dict:
                v['L4_guess'] = extra_anno_dict[k]
                
    # Full list of features
    s = sep.join(header) + '\n'
    s += '\n'.join(
        [sep.join([str(x[ii]) for ii in header]) for x in anno.values()]
    )
    with open(os.path.join(outdir, outstr.replace('.tsv', '_full.tsv')), 'w') as O:
        O.write(s)

    # recommended_ list of features
    s = sep.join(header) + '\n'
    s += '\n'.join(
        [sep.join([str(x[ii]) for ii in header]) for k,x in anno.items() if k in primary_fids]
    )
    with open(os.path.join(outdir, outstr.replace('.tsv', '_recommended.tsv')), 'w') as O:
        O.write(s)

def get_asari_default_anno(infile, ii=9):
    '''From asari Feature_annotation.tsv'''
    d = {}
    for line in open(infile).readlines()[1:]:
        a = line.split('\t')
        if len(a) > ii:
            d[a[0]] = a[ii]
    return d

def load_csm_full(db_cpds_serum, 
                  neu_csm, 
                  lib_hilicpos,
                  lib_rppos,
                  lib_hilicneg,
                  lib_rpneg,
                  ):
    '''
    Load from JSON files full CSM components: DB_CPDS_SERUM, neuCSM, method specific libraries.
    returns bmDict, neuCSM, csmlib_dict
    '''
    BMDB = json.load(open(db_cpds_serum))
    bmDict = {}
    for x in BMDB: bmDict[x['id']] = x

    neuCSM = json.load(open(neu_csm))
    lib_hilicpos = json.load(open(lib_hilicpos))
    lib_rppos = json.load(open(lib_rppos))
    lib_hilicneg = json.load(open(lib_hilicneg))
    lib_rpneg = json.load(open(lib_rpneg))
    csmlib_dict = {
        'hilicpos': lib_hilicpos,
        'rppos': lib_rppos,
        'hilicneg': lib_hilicneg,
        'rpneg': lib_rpneg
    }

    return bmDict, neuCSM, csmlib_dict

        
def test_main():
    '''
    Test/template for using annotation functions.
    Hard coded path. 
    '''
    db_cpds_serum = '/Users/lish/li.proj/consensus_serum_metabolome/v7/DB_CPDS_SERUM/blood_metabolites_20250212.json'
    CSM_PATH = '/Users/lish/li.proj/consensus_serum_metabolome/v7/DB_REF_CSM/r1/'
    neu_csm = CSM_PATH + 'r1_neu_mass_registries_annotated_1012.json'

    mDict, neuCSM, csmlib_dict = load_csm_full(db_cpds_serum, 
                    neu_csm, 
                    lib_hilicpos = CSM_PATH + 'r1_ref_hilic_pos_csmfs_20250304.json',
                    lib_rppos = CSM_PATH + 'r1_ref_rp_pos_csmfs_20250304.json',
                    lib_hilicneg = CSM_PATH + 'r1_ref_hilic_neg_csmfs_20250304.json',
                    lib_rpneg = CSM_PATH + 'r1_ref_rp_neg_csmfs_20250304.json',
                    )
    
    indir = '/Users/lish/li.collaborations/Patel_HHEAR/2_03082025/'
    list_infiles = [
        ('ST000909/study/HILICpos/asari_ST000909_HILICpos__ppm5_37222026/export/full_Feature_table.tsv', 'hilicpos',),
        ('ST000909/study/RPneg/asari_ST000909_RPneg__ppm5_37231831/export/full_Feature_table.tsv', 'rpneg',),
        ('ST002118/study/RPpos/asari_ST002118_RPpos__ppm5_3819562/export/full_Feature_table.tsv', 'rppos'),
    ]

    all_stats = {}
    for f, method in list_infiles:
        _fname_handle_ = os.path.split(f.replace('/export/full_Feature_table.tsv', ''))[1]
        outfile = indir + 'csmf_anno_hhear_' + _fname_handle_ + '.tsv'
        feature_table_outfile = indir + _fname_handle_ + 'full_Feature_table.tsv'
        
        extra_anno_dict = get_asari_default_anno(indir + f.replace("export/full_Feature_table", "Feature_annotation"))
        
        mode='pos'
        adduct_search_patterns=adduct_search_patterns_pos
        primary_ions_ordered = primary_ions_pos_ordered
        if 'neg' in method:
            mode = 'neg'
            adduct_search_patterns=adduct_search_patterns_neg
            primary_ions_ordered = primary_ions_neg_ordered
                
        anno, header, primary_fids, stats = annotate_asari_table_by_csm(
                                indir + f, 
                                KCD_formula_coordinate, 
                                isotope_search_patterns,
                                adduct_search_patterns,
                                extended_adducts,
                                mDict,
                                csmlib_dict[method],
                                neuCSM,
                                primary_ions_ordered = primary_ions_ordered,
                                mode=mode, 
                                method=method,
                                mz_tolerance_ppm=5, rt_tolerance=2, 
                                quality_filter=False,
                )
        
        all_stats[_fname_handle_] = stats
        # This outputs two annotation files.
        write_csm_anno_2files(
            anno, header, extra_anno_dict, primary_fids, "", outfile)
        
        # make new copy of feature table
        # with open(feature_table_outfile, 'w') as O:
        #    O.write(open(indir + f).read())

    # Write stats
    with open(os.path.join(indir, '__stats.json'), 'w', encoding='utf-8') as O:
        json.dump(all_stats, O,  ensure_ascii=False, indent=2) 
        