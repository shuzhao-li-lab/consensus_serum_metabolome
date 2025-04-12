import os
from mining import *
from csm_align import *
from khipu_custom import *


KCD_formula_coordinate = build_KCD_from_formula_coordinate(formula_coordinate)

def get_asari_default_anno(infile, ii=9):
    '''From asari Feature_annotation.tsv'''
    d = {}
    for line in open(infile).readlines()[1:]:
        a = line.split('\t')
        if len(a) > ii:
            d[a[0]] = a[ii]
    return d

def annotate_asari_table_by_csm(
                        asari_table, 
                        KCD_formula_coordinate, 
                        isotope_search_patterns,
                        adduct_search_patterns,
                        extended_adducts,
                        bmDict, #=bmDict,
                        reflib, #=lib_rppos, 
                        neuDict, #=neuCSM,
                        primary_ions_ordered, 
                        mode = 'pos', method='rppos', 
                        mz_tolerance_ppm=5, rt_tolerance=2,
                        quality_filter=False,
                        ):
    '''
    The relationships among features need to be established before annotation.
    Step 1. khipu pre-annotation.
    Step 2. mapping unique epds to CSM.
    Step 3. Complement annotation by L4 in asari (can be replaced later).
    
    Output stats, recommended anno table and full anno table. 
    
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
                        bmDict, #=bmDict,
                        reflib, #=lib_rppos, 
                        neuDict, #=neuCSM,
                        primary_ions_ordered,
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
                        # mode = 'pos', method='rppos',
                        quality_filter=False
                        ):
    '''
    returns anno dict, header, list of primary_fids, stats
    
    '''
    clean_dict_match, featureDict, primary_fids = csmf_annotate_from_separate_lists(
        list_khipus, list_singletons, list_features, reflib, 
        primary_ions_ordered=primary_ions_ordered,
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
    with open(os.path.join(outdir, 'full_'+outstr), 'w') as O:
        O.write(s)

    # recommended_ list of features
    s = sep.join(header) + '\n'
    s += '\n'.join(
        [sep.join([str(x[ii]) for ii in header]) for k,x in anno.items() if k in primary_fids]
    )
    with open(os.path.join(outdir, 'recommended_'+outstr), 'w') as O:
        O.write(s)


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
    BMDB = json.load(open('blood_metabolites_20250212.json'))
    neuCSM = json.load(open('DB/r1_neu_mass_registries_annotated_1012.json'))

    lib_hilicpos = json.load(open('DB/r1_ref_hilic_pos_csmfs_20250304.json'))
    lib_rppos = json.load(open('DB/r1_ref_rp_pos_csmfs_20250304.json'))

    lib_hilicneg = json.load(open('DB/r1_ref_hilic_neg_csmfs_20250304.json'))
    lib_rpneg = json.load(open('DB/r1_ref_rp_neg_csmfs_20250304.json'))
    
    
    '''


    indir = '/Users/lish/li.collaborations/Patel_HHEAR/2_03082025/'
    list_infiles = [
        ('ST000909/study/HILICpos/asari_ST000909_HILICpos__ppm5_37222026/export/full_Feature_table.tsv', 'hilicpos',),
        ('ST000909/study/RPneg/asari_ST000909_RPneg__ppm5_37231831/export/full_Feature_table.tsv', 'rpneg',),
        ('ST001004/study/HILICpos/asari_ST001004_HILICpos__ppm5_37235537/export/full_Feature_table.tsv', 'hilicpos',),
        ('ST001004/study/RPneg/asari_ST001004_RPneg__ppm5_380114/export/full_Feature_table.tsv', 'rpneg',),
        ('ST001209/study/HILICpos/asari_ST001209_HILICpos__ppm5_3802438/export/full_Feature_table.tsv', 'hilicpos'),
        ('ST001209/study/RPneg/asari_ST001209_RPneg__ppm5_3812124/export/full_Feature_table.tsv', 'rpneg',),
        ('ST001428/study/HILICpos/asari_ST001428_HILICpos__ppm5_382420/export/full_Feature_table.tsv', 'hilicpos'),
        ('ST001428/study/RPneg/asari_ST001428_RPneg__ppm5_384620/export/full_Feature_table.tsv', 'rpneg',),
        ('ST001932/study/HILICpos/asari_ST001932_HILICpos__ppm5_3852110/export/full_Feature_table.tsv', 'hilicpos'),
        ('ST001932/study/RPneg/asari_ST001932_RPneg__ppm5_38141935/export/full_Feature_table.tsv', 'rpneg'),
        ('ST002118/study/RPpos/asari_ST002118_RPpos__ppm5_3819562/export/full_Feature_table.tsv', 'rppos'),
    ]


    for f, method in list_infiles:
        _fname_handle_ = os.path.split(f.replace('/export/full_Feature_table.tsv', ''))[1]
        outfile = indir + 'csmf_anno_hhear_' + _fname_handle_ + '.tsv'
        feature_table_outfile = indir + _fname_handle_ + 'full_Feature_table.tsv'
        
        extra_anno_dict = get_asari_default_anno(indir + f.replace("export/full_Feature_table", "Feature_annotation"))
        
        mode='pos'
        adduct_search_patterns=adduct_search_patterns_pos
        # dict_member_size = pos_dict_member_size
        if 'neg' in method:
            mode = 'neg'
            adduct_search_patterns=adduct_search_patterns_neg
            # dict_member_size = neg_dict_member_size
        
        
        
        
        anno, header = annotate_asari_table_by_csm(
                            indir + f, 
                            KCD_formula_coordinate, 
                            isotope_search_patterns,
                            adduct_search_patterns,
                            extended_adducts,
                            # dict_member_size=dict_member_size,
                            mode = mode, 
                            method=method,
                            mz_tolerance_ppm=5, rt_tolerance=2,
                            bmDict=bmDict,
                            reflib= csmlib_dict[method], 
                            neuDict=neuCSM,
                            outfile=outfile,
                            sep='\t'
        )
        write_anno_csmf_format(anno, header, extra_anno_dict, outfile)
        
        # move feature table
        with open(feature_table_outfile, 'w') as O:
            O.write(open(indir + f).read())

