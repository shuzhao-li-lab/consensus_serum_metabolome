'''
CSM includes reference metabolomes per method, i.e. HILIC+, RP+, HILIC-, RP- .

Select CSM features with RTI, blood conc/papers, and likely unambiguous

Per dataset:
Primary features present in CSM. Regression to get RT dict.

DB: blood_metabolites (DB_CPD_SERUM) should have inchikey and major DB identifiers.
'''
import json
import numpy as np
from scipy.stats import linregress

import matplotlib.pyplot as plt
import seaborn as sn

from asari.tools import match_features as mf
from mining import get_M0, get_M1

#
# I/O and conversions
#

def get_primary_feature_from_a_khipu(khipu, primary_ions_ordered):
    '''
    Find first M0 feature in primary_ions_ordered.
    Different from the mining.get_M0 function, which looks for highest peak in M0 features. 
    
    E.g. primary_ions_pos_ordered = ['M0,M+H+', 'M0,Na/H', 
                            'M0,M+H+, 2x charged', 'M0,Na/H, 2x charged', 
                            'M0,M+H+, 3x charged', 'M0,Na/H, 3x charged']
    '''
    _d = {}
    for f in khipu['MS1_pseudo_Spectra']:   # f is not necessarily ordered
        if f['ion_relation'] in primary_ions_ordered:
            _d[f['ion_relation']] = f
            
    new = [_d[ii] for ii in primary_ions_ordered if ii in _d]
    if new:
        return new[0]
    else:    # backup if no common pri ion is present
        return get_M0(khipu['MS1_pseudo_Spectra'])

def get_primary_features_from_epds(infile_json, 
                                   natural_ratio_limit=0.5, 
                                   primary_ion='M0,M+H+'):
    full = json.load(open(infile_json))
    landmark_features, primary_features, singletons = get_primary_features_from_epdlist(
        full, natural_ratio_limit, primary_ion
    )
    return landmark_features, primary_features, singletons

def get_all_features_from_epdlist(epdlist):
    '''
    This function is used when input data are read from JSON annotation, 
    and some types need fix.
    Direct results from in situ khipu avoid type problems. 
    '''
    khipu_featuress, singletons = [], []
    for epd in epdlist:
        features = epd["MS1_pseudo_Spectra"]
        for f in features:
            f['id'] = f['id_number']
            f['snr'] = float(f['snr'])
            f['goodness_fitting'] = float(f['goodness_fitting'])
            f['representative_intensity'] = f['peak_area'] = float(f['peak_area'])
        if len(features) == 1:
            singletons.append(features[0])
        else:
            khipu_featuress += features
    return khipu_featuress, singletons

def extract_12c_landmark_features(list_khipus, ion='M0,M+H+'):
    '''
    Simple function returns landmark_features by selecting ion type. 
    Function get_primary_features_from_epdlist applies intensity ratios.
    '''
    landmark_features = []
    for k in list_khipus:
        m0 = [f for f in k['MS1_pseudo_Spectra'] if f['ion_relation']==ion]
        if m0:
            landmark_features += m0
    return landmark_features

def get_primary_features_from_epdlist(epdlist, 
                                   natural_ratio_limit=0.5, 
                                   primary_ion='M0,M+H+'):
    '''
    Use khipu JSON result as input.
    Use primary_ion for pos or neg ionization data.
    Returns landmark_features, primary_features, singletons
    
    This is used for selecting landmarks for alignment etc. 
    The selection of M0, M1 here is not complete (ions prioritized by max intensity), 
    not suitable for complete feature annotation.
    Use get_all_features_from_epdlist for complete annotation.
    '''
    # fix typing
    for epd in epdlist:
        for f in epd["MS1_pseudo_Spectra"]:
            f['id'] = f['id_number']
            f['snr'] = float(f['snr'])
            f['goodness_fitting'] = float(f['goodness_fitting'])
            f['representative_intensity'] = f['peak_area'] = float(f['peak_area'])
            
    print("Total feature number", len(epdlist))
    landmark_features, primary_features, singletons = [], [], []
    for epd in epdlist:
        features = epd["MS1_pseudo_Spectra"]
        if len(features) == 1:
            singletons.append(features[0])
        else:
            M0, M1 = get_M0(features), get_M1(features)
            primary_features.append(M0)
            if M0 and M1 and M0['modification']==M1['modification'] and M0['ion_relation']==primary_ion:
                if float(M1['representative_intensity'])/(1 + float(
                    M0['representative_intensity'])) < natural_ratio_limit:
                    landmark_features.append( M0 )

    print("landmark_features, primary, singletons:", 
          len(landmark_features), len(primary_features), len(singletons))
    
    return landmark_features, primary_features, singletons

def convert_track2features(features):
    '''returns 
    tracks, and dict of massTrack to [features]'''
    tracks = []
    d = {}
    for f in features:
        t = f['parent_masstrack_id']
        if t in d:
            d[t].append(f['id'])
        else:
            d[t] = [f['id']]
            tracks.append({
                'id': t,
                'mz': f['mz'],
                'rtime': 0,
                'members': []
            })
                
    for track in tracks:
        track['members'] = d[track['id']]
    return tracks, d

def convert_csmfs2tracks(csmfs):
    '''
    Example ref_csmfs :     {'id': 'r1_pos_1274.858866_HILIC_0',
                                'ion': '',
                                'neuMR': 'r1_neu_1273.851586',
                                'mz': 1274.858866,
                                'rti': None,
                                'number_isomers': 1,
                                'isomer_elution_order': 0,
                                'annotation_reported': 'MDB0094924r1'}
    '''
    tracks = []
    d = {}
    for f in csmfs:
        t = '_'.join(f['id'].split('_')[:-2])
        if t in d:
            d[t].append(f['id'])
        else:
            d[t] = [f['id']]
            tracks.append({
                'id': t,
                'mz': f['mz'],
                'rtime': 0,
                'members': []
            })
    for track in tracks:
        track['members'] = d[track['id']]
    return tracks, d

#
# CSMF reference based alignment and annotation
#

def align_user_features_to_csm(primary_features, singletons, 
                               list_csmf, mz_ppm=5
                               ):
    '''
    return fully matched list of CSMFs for annotation of user data. 
    Primary features, excluding other features in khipus, 
        are used because CSMFs are mostly primary features, and 
        annotation of a primary ion covers the full khipu. 
        Other features in Khipu get transferred annotation outside this function.
    
    Ideas considered but likely poor ones: 
        slope and intercept can be used to calculate feature RT to match to RTIs of CSMFs.
        exact_isomer_numbers : force identical number of isomers. 
        Not enforced here as top N features will be selected for matching.
    
    The annotation is two-steps:  massTracks then features.
    massTracks matches should be unique. 
    Remaining features (annotate_without_aligned) can get matches somewhere, e.g. L4 in asari default HMDB match.

    returns mapping dict.
    '''
    fTracks, dict_fTracks = convert_track2features(primary_features + singletons)
    csmfTracks, dict_csmfTracks = convert_csmfs2tracks(list_csmf)

    dict_features, dict_csmfs = {}, {}
    for f in primary_features + singletons:
        dict_features[f['id']] = f
    for c in list_csmf:
        dict_csmfs[c['id']] = c
    
    matched = mf.best_mz_match_lcms_features(
                    csmfTracks, fTracks, mz_ppm=mz_ppm, rt_tolerance=1e20
                    )  # returns mass track level matches as dict
    
    clean_dict_match = get_clean_dict_match(
        matched, dict_csmfTracks, dict_fTracks, dict_csmfs, dict_features
        )
    return clean_dict_match

def annotate_without_aligned(primary_features, singletons, clean_dict_match, list_neumr):
    '''
    Placeholder.
    The user features that failed in alignment to CSMFs can be annotated by m/z matching. 
    Lower confidence, but final step to round up the annotation.
    remaining = [f for f in primary_features + singletons if f['id'] not in clean_dict_match]
    # find in neuMR list
    '''
    pass

def csmf_annotate_from_epd_json(json_epd_file, 
                                reflib,
                                method,
                                primary_ion='M0,M+H+'
                                ):
    '''
    Obsolete. Use csm_annotate.annotate_epdJSON_by_csm. 
    '''
    epdlist = json.load(open(json_epd_file))
    clean_dict_match, featureDict = csmf_annotate_from_epd_list(epdlist, 
                                reflib, # method, primary_ion
                                )
    return clean_dict_match, featureDict

def csmf_annotate_from_epd_list(epdlist, 
                                reflib,
                                ):
    '''
    Obsolete. Use csm_annotate.annotate_epdJSON_by_csm. 
    Returns mappd dict btw user features and CSM features.
    
    Some functions like this one are for conversion only. 
    Because list of empCpds is passed with both khipus and singletons, repacking is needed. 
    We should avoid these conversions if we can.
    
    Parameters removed: method, primary_ion='M0,M+H+'
    '''
    this_khipu_features, this_singletons = get_all_features_from_epdlist(epdlist)
    featureDict = {}
    for f in this_khipu_features + this_singletons:
        featureDict[f['id']] = f
    clean_dict_match =  align_user_features_to_csm(
                    this_khipu_features, this_singletons, reflib.values(), 
                    ) 
    return clean_dict_match, featureDict

def csmf_annotate_from_separate_lists(list_khipus, list_singletons, list_features,
                                csm_reflib,
                                primary_ions_ordered,
                                mz_tolerance_ppm
                                ):
    '''
    Returns mappd dict btw user features and CSM features, and a feature dict. 
    list_khipus, list_singletons, list_features = run_khipu_insitu(asari_table, ...)
    '''
    featureDict = {}
    for f in list_features:
        featureDict[f['id']] = f
        
    # separate primary features, redundant khipu features, singletons
    # feature['parent_epd_id'] points back to khipu
    list_pri_features = [
        get_primary_feature_from_a_khipu(kk, primary_ions_ordered) for kk in list_khipus
    ]
    clean_dict_match =  align_user_features_to_csm(
                    list_pri_features, list_singletons, csm_reflib.values(), 
                    mz_ppm=mz_tolerance_ppm
                    )   # feature ID to CSM_feature ID
    
    # expand anno to redundant khipu features
    k2anno = {}
    for k,v in clean_dict_match.items():
        if 'parent_epd_id' in featureDict[k]:
            k2anno[featureDict[k]['parent_epd_id']] = v
    for kk in list_khipus:
        for f in kk['MS1_pseudo_Spectra']:
            if kk['interim_id'] in k2anno:                              # not all khipus have a match
                clean_dict_match[f['id']] = k2anno[kk['interim_id']]    # okay to overwrite pri features
    
    return clean_dict_match, featureDict, [f['id'] for f in list_pri_features]


def get_formula_from_neumr(n, neuDict):
    if 'formula_matches' in neuDict[n] and neuDict[n]['formula_matches']:
        return neuDict[n]['formula_matches'][0]['formula']
    else:
        return ''

def export_csmf_annotation(clean_dict_match, 
                           featureDict,
                           bmDict,
                           reflib,
                           neuDict,
                           ):
    '''
    returns annotation dict and header.
    Delivery format via dict: 
        Feature ID, mz, rtime, ion_relation, CSMF ID, 'ion_csm', 
        neuMR ID, formula, top recommendation (MDB ID, name), top3 [(MDB ID, name), ...], 
        score dict. 
    '''
    bmDict['unknown'] = {'name': 'unknown'}
    anno = {}
    header = ['feature_ID', 'mz', 'rtime', 'ion_relation', 
              'CSMF_ID', 'ion_csm', 'neuMR_ID', 'formula', 
              'top_recommendation', 'top_recommendation_id', 'top_recommendation_name', 'top_recommendation_score', 
              'top3', 'full_score_dict']
    
    for fid, csmf_id in clean_dict_match.items():
        _neu = reflib[csmf_id]['neuMR']
        _feature = featureDict[fid]
        top1, top3, scores = '', '', ''
        top1_id, top1_name, top1_score = '', '', ''
            
        if 'score_dict' in reflib[csmf_id]:
            scores = sorted([(v,k) for k,v in reflib[csmf_id]['score_dict'].items()], 
                               reverse=True)
            scores = [x for x in scores if 'unknown' != x[1]]
            top3 = [
                (x[0], x[1], bmDict[x[1]]['name']) for x in scores[:3]
            ]
            if top3:
                top1 = top3[0]
                (top1_score, top1_id, top1_name) = top1
            else:
                top3, scores = '', ''
        
        anno[fid] = dict(zip(header, (fid, _feature['mz'], _feature['rtime'], 
                _feature.get('ion_relation', ''), # may not exist for singletons
                csmf_id, reflib[csmf_id]['ion'], _neu, get_formula_from_neumr(_neu, neuDict),
                top1, top1_id, top1_name, top1_score, top3, scores 
                )))

    return anno, header
    
def write_anno_csmf_format_old(anno, header, outfile='test_csmf_anno.tsv', sep='\t'):
    # anno as dict of fid: dict
    s = sep.join(header) + '\n'
    s += '\n'.join(
        [sep.join([str(x[ii]) for ii in header]) for x in anno.values()]
    )
    with open(outfile, 'w') as O:
        O.write(s)

#
# other utils
#

def get_clean_dict_match(matched, dict_csmfTracks, dict_fTracks, 
                         dict_csmfs, dict_features):
    '''
    To resolve matches to feature level. 
    The input matched as dict (m/z matches), at mass track level. 
    '''
    new = {}
    for k,v in matched.items():
        # at track level, massTracks matches should be unique
        csmfs = dict_csmfTracks[k]
        features = dict_fTracks[v]
        resolved_pairs = resolve_multiple_match(
            csmfs, features, dict_csmfs, dict_features) # at feature level
        for x,y in resolved_pairs:
            new[y] = x   # dict feature -> csmf
    return new

def match_csmf_features_equal_n(csmfs, features, dict_csmfs, dict_features):
    '''
    match csmfs and features by elution order.
    Exact RT and RTI are not used here, but will be used elsewhere and in future improvements.
    Returns list of pairs.
    '''
    # 
    ordered_csmfs = sorted(csmfs, key=lambda x: dict_csmfs[x]['isomer_elution_order'])
    ordered_features = sorted(features, key=lambda x: dict_features[x]['rtime'])
    return list(zip(ordered_csmfs, ordered_features))

def resolve_multiple_match(csmfs, features, dict_csmfs, dict_features):
    '''
    csmfs : csm feature IDs
    features : feature IDs to be annotated
    
    resolve multiple matches btw features and CSMFs on the same track.
    This is based on 
    - matching elution order
    - prioritization of features, annotations and MDB records by signal levels.
    '''
    N_csmf, N_features = len(csmfs), len(features)
    if N_csmf == 1 and N_features == 1:
        # we can add RT check here in the future, to reduce false matches
        return [(csmfs[0], features[0])]
    
    else: # get equal numbers of CSMFs and user features to use
        if N_csmf > N_features:
            sel_features = features
            csmfs = select_topN_csmfs(csmfs, N_features, dict_csmfs)
        elif N_csmf < N_features:
            # N_csmf < N_features, sort features by abundance; take first matched N
            sel_features = sorted(features, 
                        key=lambda x: dict_features[x]['peak_area'], reverse=True
                        )[:N_csmf]
        else:    # N_csmf == N_features
            sel_features = features
        
        return match_csmf_features_equal_n(
            csmfs, sel_features, dict_csmfs, dict_features
            )

def select_topN_csmfs(csmfs, N, dict_csmfs):
    '''
    Select top N CSMFs on the same m/z, ranked by number of input member features.
    No longer using dict_member_szie. number_contributing_members in CSMF after 20250304.
    
    : how many input members contributed to this consensus CSMF.
    '''
    return sorted(csmfs, key=lambda x: dict_csmfs[x]['number_contributing_members'], reverse=True)[:N]


def select_topN_csmfs_retired(csmfs, N, dict_csmfs, dict_member_size):
    '''
    Select top N CSMFs on the same m/z, ranked by number of input member features.
    dict_member_szie : how many input members contributed to this consensus CSMF.
    '''
    return sorted(csmfs, key=lambda x: dict_member_size[dict_csmfs[x]['id']], reverse=True)[:N]

def get_MDB_records(NMR):
    if 'MDB_records' in NMR:
        return [x['id'] for x in NMR['MDB_records']]
    else:
        return []

def get_hmdbinchikey(hid, hmdb2inchikey):
    r = hmdb2inchikey.get(hid, '')
    if not r:
        r = hmdb2inchikey.get(hid.replace('HMDB', 'HMDB00'), '') # try fixing ID version
    return r

def get_inchikey_match(m, hmdb2inchikey):
    '''
    get_inchikey from userAnnotation, 
    which could be Metabolon with old HMDB ID or inchikey from ADAP
    '''
    if 'inchikey' in m and m['inchikey']:
        return m['inchikey']
    elif 'hmdb' in m and m['hmdb']:
        return get_hmdbinchikey(m['hmdb'], hmdb2inchikey)
    else:
        return ''

def select_mdb(LL, 
               bmDict,
               method_key='RTI_HILIC_RETIP', 
               priorities = ['BloodPaperCount', 'blood_conc', 'xlogP'],
               ):
    '''
    return best guess of bmDict members that match to a CSMF .
    
    input e.g. 'MDB': ['MDB0000196r1', 'MDB0000550r1', 'MDB0006678r1', ...]
    priorities not fully implemented here. 
    This is for preparing data to get RTI regression. 
    A detailed annotation method will be later elsewhere.
    '''
    LL = [bmDict[x] for x in LL]
    prioritized = [x for x in LL if method_key in x]
    if not prioritized:
        prioritized = [x for x in LL if priorities[0] in x]
    else:
        prioritized = [x for x in prioritized if priorities[0] in x]
        
    prioritized.sort(key=lambda x: x[priorities[0]], reverse=True)
    if prioritized:
        return prioritized[0]
    else:
        return None
    
def test_match_feature_number(study_number_ontrack, list_matched_cmsfs):
    return [
        x for x in list_matched_cmsfs if x['number_isomers'] == study_number_ontrack
    ]

def select_user_features_match_csm(landmark_features, pos_csmfsTree):
    '''Get uniquely matched landmark_features, used for regression of RTI.
    '''
    dict_mtrack_density = {}
    _mtracks = [x['parent_masstrack_id'] for x in landmark_features]
    for x in landmark_features:
        dict_mtrack_density[x['id']] = _mtracks.count(x['parent_masstrack_id'])
    matched = []
    for x in landmark_features:
        if dict_mtrack_density[x['id']] == 1:
            _m = find_all_matches_centurion_indexed_list(x['mz'], 
                                pos_csmfsTree, 5)
            if _m and len(_m)==1:
                matched.append((x, _m))
    return matched
    
def model_regress_rti_1step(list_matched, bmDict, rti='RTI_HILIC_RETIP', do_plot=False, plot_outfile='rg.pdf'):
    '''regression step to compute RTI, not using tracks with multiple isomers
    Some blood metabolite entries in bmDict has RTI values, e.g. RTI_HILIC_RETIP.
    Returns slope, intercept
    
    Note:
    Function model_regress_rti_2step is preferred to filter outliers.
    These functions use known RTI in blood metabolite collection to do regression.
    '''
    xx2_rtime, yy2_rti = [], []
    for m, LL in list_matched:
        csmf = LL[0]
        if 'MDB' in csmf:
            sel = select_mdb(csmf['MDB'], bmDict, method_key=rti)
            if sel and rti in sel:
                yy2_rti.append(sel[rti])
                xx2_rtime.append(m['rtime'])
                    
    RG = linregress(xx2_rtime, yy2_rti)
    if do_plot:
        plt.figure(figsize=(3,3))
        sn.regplot(x=xx2_rtime, y=yy2_rti)
        plt.savefig(plot_outfile)
        
    return RG.slope, RG.intercept

def get_res_stdev(xx, yy, slope, intercept):
    errors = []
    for ii in range(len(yy)):
        errors.append( yy[ii] - slope*xx[ii] - intercept)
    return np.std(errors)

def model_regress_rti_2step(list_matched, bmDict, rti='RTI_HILIC_RETIP', do_plot=False, plot_outfile='rg.pdf'):
    '''regression step to compute RTI, not using tracks with multiple isomers
    Some blood metabolite entries in bmDict has RTI values, e.g. RTI_HILIC_RETIP.

    To do 2-step regression
    rm data outliers after 1st round of regression
    
    Returns slope, intercept
    '''
    xx2_rtime, yy2_rti = [], []
    for m, LL in list_matched:
        csmf = LL[0]
        if 'MDB' in csmf:
            sel = select_mdb(csmf['MDB'], bmDict, method_key=rti)
            if sel and rti in sel:
                yy2_rti.append(sel[rti])
                xx2_rtime.append(m['rtime'])
                    
    RG = linregress(xx2_rtime, yy2_rti)
    err = get_res_stdev(xx2_rtime, yy2_rti, RG.slope, RG.intercept)
    # round 2
    newxx, newyy = [], []
    for ii in range(len(yy2_rti)):
        if abs(yy2_rti[ii] - RG.slope*xx2_rtime[ii] - RG.intercept) < 2 * err:
            newxx.append(xx2_rtime[ii])
            newyy.append(yy2_rti[ii])
    RG = linregress(newxx, newyy)
    
    if do_plot:
        plt.figure(figsize=(3,3))
        sn.regplot(x=newxx, y=newyy)
        plt.savefig(plot_outfile)
        
    return RG.slope, RG.intercept
