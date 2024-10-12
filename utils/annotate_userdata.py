import json

from jms.dbStructures import ExperimentalEcpdDatabase
from asari.default_parameters import adduct_search_patterns, \
    adduct_search_patterns_neg, isotope_search_patterns, extended_adducts
    
from .mining import *

def CSM_annoate_FeatureTable(ftable, 
                             neuMR,
                             mode = 'pos_HILIC',
                             mz_tolerance_ppm=5, 
                             rt_tolerance=2
                             ):
    '''
    neuMR : annotation dict of neuMR, neutral mass registry in CSM.
    return Fullfeature2anno dictionary, which can be ported to JSON and better .tsv for nonprogrammers.
    
    Test:
    f = '../annotation_datasets/MT02/FeatureTables/full_Feature_table.tsv'
    num_samples, mt02 = read_features_from_asari_table(open(f).read())
    mt02 = get_empCpds_from_FeatureList(mt02)
    
    neuMR = json.load(open('DB/r1_neu_mass_registries_annotated.json'))
    
    '''
    num_samples, plist = read_features_from_asari_table(open(ftable).read())
    epds = get_empCpds_from_FeatureList(plist, mode=mode, 
                                        mz_tolerance_ppm=mz_tolerance_ppm, 
                                        rt_tolerance=rt_tolerance)
    
    Fullfeature2anno, Fullepd2anno = annotate_epds_on_csm(list(epds.values()),
                                              neuMR,
                                              mode,
                                              5,
                                              )

    return Fullfeature2anno


def get_empCpds_from_FeatureList(FeatureList,
                                    mode='pos', 
                                    mz_tolerance_ppm=5, 
                                    rt_tolerance=2
                                ):
    '''
    
    There's an early version in minging.features2cpds
    
    calculate_neutral_mass_pos, calculate_neutral_mass_neg from mining
    '''
    EED = ExperimentalEcpdDatabase(mode=mode, 
                                    mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=rt_tolerance)
    if mode == 'pos':
        EED.adduct_patterns = adduct_search_patterns
        calculate_neutral_mass_ = calculate_neutral_mass_pos
        default_ion = 'M0,M+H+'
    else:
        EED.adduct_patterns = adduct_search_patterns_neg
        calculate_neutral_mass_ = calculate_neutral_mass_neg
        default_ion = 'M0,M-H-'
    EED.isotope_search_patterns = isotope_search_patterns
    EED.extended_adducts = extended_adducts

    EED.build_from_list_peaks(FeatureList,
                              charges=[1, 2, 3],
                              has_parent_masstrack=True, # False if FeatureList not from asari
                              )
    # this got EED.dict_empCpds. Next do singletons
    singletons = {}
    for k, peak in EED.dict_peaks.items():
        if k not in EED.peak_to_empCpd:
            peak['ion_relation'] = default_ion
            neutral_formula_mass = calculate_neutral_mass_(peak['mz'], default_ion)
            singletons[k] = {
                'interim_id': '_singleton_' + k,
                'neutral_formula_mass': neutral_formula_mass,
                'MS1_pseudo_Spectra': [peak],
            }
    EED.dict_empCpds.update(singletons)
    return EED.dict_empCpds

def make_dicts_features2cpds(list_empCpds):
    '''
    Make two dicts of {feature ID: (epd ID, ion)}, {epd ID: (feature ID, ion)}
    '''
    f2c, c2f = {}, {}
    for epd in list_empCpds:
        c2f[epd['interim_id']] = []
        for f in epd['MS1_pseudo_Spectra']:
            f2c[f['id']] = (epd['interim_id'], f['ion_relation'])
            c2f[epd['interim_id']].append((f['id'], f['ion_relation']))
            
    return f2c, c2f

def annotate_epds_on_csm(list_empCpds, 
                         csm, 
                         method='pos_HILIC', 
                         mz_ppm=5, 
                         rt_tolerance=1e20):
    '''
    list_empCpds : including both khipus and singletons, with neutral mass
    csm : in the format of neuMR with all annotations
    
    # from asari.tools import match_features as mf
    '''
    list_csm, list_epd, _result = [], [], []
    for x in csm.values():
        list_csm.append({'mz': x['mass'], 'rtime': 0, 'id': x['id']})
    for x in list_empCpds:
        if x['neutral_formula_mass']:
            list_epd.append({'mz': x['neutral_formula_mass'], 'rtime': 0, 'id': x['interim_id']})
            
    matched = mf.list_match_lcms_features(
        list_epd, list_csm, mz_ppm=mz_ppm, rt_tolerance=rt_tolerance)
        # [('kp267_218.1152', ['r1_neu_218.115403'])]
    for k,v in matched.items():
        for x in v:
            a, n = fetch_neumr_anno(csm[x], method=method)
            _result.append((k, a, n))
            
    _d, feature2anno, epd2anno = {}, {}, {}
    for x in _result:
        k, Anno, N = x
        _d[k] = recommend_cpds(Anno, N)
    
    f2c, c2f = make_dicts_features2cpds(list_empCpds)
    
    for k in _d:
        if _d[k][0]: # has recommendation
            epd2anno[k] = _d[k]
            for f in c2f[k]:
                feature2anno[f[0]] = (f[1], _d[k])
                
    return feature2anno, epd2anno
    
    
def rank_MDB_records(MDB_records):
    '''
    Rank by BloodPaperCount, blood_conc
    '''
    L1, L2, L3 = [], [], []
    for x in MDB_records:
        if 'BloodPaperCount' in x and x['BloodPaperCount']:
            L1.append(x)
        elif 'blood_conc' in x and x['blood_conc']:
            L2.append(x)
        else:
            L3.append(x)
            
    return sorted(L1, key=lambda x: x['BloodPaperCount'], reverse=True
                  ) + sorted(L2, key=lambda x: x['blood_conc'], reverse=True) + L3


def recommend_cpds(annodict, N=1, max_recommend=3):
    '''
    Prototyping annotation recommendation.
    
    Will implement using feature number on mass track later.
    
        result['User_annotations'] = UA
    result['MDB_records'] = MR.get('MDB_records', [])
    if not UA and not result['MDB_records']:
        result['formula_matches'] = MR.get('formula_matches', [])
    
    '''
    def clean_user_anno(LL):
        '''
                    if isinstance(LL[0], dict):
                for x in LL:
                    new.append(x)
            else:
        
        '''
        new = []
        for L in LL:
            if isinstance(L[1], dict):
                new.append(L[1])
            elif isinstance(L[1], list):
                new += L[1]
        return new
    
    top, recomm3 = None, []
    if not annodict['User_annotations'] and not annodict['MDB_records']:
        return top, annodict.get('formula_matches', [])
    
    elif annodict['User_annotations'] and annodict['MDB_records']:
        return clean_user_anno(annodict['User_annotations'][:N]), \
                    rank_MDB_records(annodict['MDB_records'])[:max_recommend]
    
    elif annodict['User_annotations']:
        rr = clean_user_anno(annodict['User_annotations'])
        return rr[:N], rr[N: max_recommend]
    
    elif annodict['MDB_records']:
        rr = rank_MDB_records(annodict['MDB_records'])
        return rr[:N], rr[N: max_recommend]
    

def fetch_neumr_anno(MR, method='pos_HILIC'):
    '''
    Fetch annotation from MR, preferring the method method.
    
    Priority: 
    User_annotations
    MDB_records
    formula_matches
    
    Protyping function
    
    '''
    result = {}
    ionmode, column = method.split('_')
    try:
        number_epds = len(MR['Methods'][method]['epds'])
    except KeyError:
        number_epds = 0
    
    UA = []
    if MR['User_annotations']:
        for k,v in MR['User_annotations'].items():
            if ionmode in k and column in k:
                # perfect method match
                UA += v
        if not UA:
            if ionmode in k: 
                UA += v
    
    result['User_annotations'] = UA
    result['MDB_records'] = MR.get('MDB_records', [])
    if not UA and not result['MDB_records']:
        result['formula_matches'] = MR.get('formula_matches', [])
    
    return result, number_epds









