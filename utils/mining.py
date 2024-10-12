import os
import copy
import json
import csv
from collections import namedtuple

import numpy as np
from scipy.signal import find_peaks 


import statsmodels.api as sm

from mass2chem.lib.formula_coordinate import formula_coordinate

from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from asari.tools import match_features as mf

from asari.default_parameters import adduct_search_patterns, \
    adduct_search_patterns_neg, isotope_search_patterns, extended_adducts


# Primary ions used to select primary vTracks
primary_ions_pos_ordered = ['M0,M+H+', 'M0,Na/H', 
                            'M0,M+H+, 2x charged', 'M0,Na/H, 2x charged', 
                            'M0,M+H+, 3x charged', 'M0,Na/H, 3x charged']
primary_ions_neg_ordered = ["M0,M-H-", 
                            "M0,M-H-, 2x charged", 
                            "M0,M-H-, 3x charged"]
primary_ions_pos = set(primary_ions_pos_ordered)
primary_ions_neg = set(primary_ions_neg_ordered)


#
# io related
#
def read_features_from_asari_table(text_table, 
                        id_col=0, mz_col=1, rtime_col=2, 
                        left_base=3,
                        right_base=4,
                        parent_masstrack_id=5, 
                        peak_area=6, cSelectivity=7, goodness_fitting=8, snr=9, detection_counts=10,
                        delimiter="\t"):
    '''
    Read a feature table from asari result.
    Returns
    List of peaks.
    '''
    featureLines = text_table.splitlines()
    header = featureLines[0].split(delimiter)
    num_features = len(featureLines)-1
    num_samples = len(header) - 11    # samples start at col 11 in asari output
    # sanity check
    print("table header looks like: \n  ", header[:20])
    print("Read %d feature lines" %num_features)
    L = []
    for ii in range(1, num_features+1):
        if featureLines[ii].strip():
            a = featureLines[ii].split(delimiter)
            L.append({
                'id_number': a[id_col], 
                'id': a[id_col],            # deliberate 
                'mz': float(a[mz_col]), 
                'rtime': float(a[rtime_col]),
                'apex': float(a[rtime_col]),
                'left_base': float(a[left_base]),
                'right_base': float(a[right_base]),
                'parent_masstrack_id': a[parent_masstrack_id],
                'peak_area': a[peak_area],
                'cSelectivity': a[cSelectivity],
                'goodness_fitting': a[goodness_fitting],
                'snr': a[snr],
                'detection_counts': a[detection_counts],
            })
    return num_samples, L

def export_json_to_table(j, outfile, sep='\t'):
    '''
    Convert btw JSON list and SQL styles when compatible, as tables are more readable and smaller in file size.
    The JSON fields are used as table columns.
    '''
    fields = list(j[0].keys())
    if 'id' in fields:
        fields.remove('id')
        fields = ['id'] + fields
        
    s = sep.join(fields) + '\n'
    for line in j:
        s += sep.join([ str(line[ii]) for ii in fields ]) + '\n'
        
    with open(outfile, 'w') as O:
        O.write(s)

def read_table_to_json(file, sep='\t'):
    '''
    Read btw SQL style table into JSON list, assuming 1st line as header.
    The table columns are used as JSON fields.
    Input is treated as str as types are not detected.
    '''
    LL = []
    lines = open(file).readlines()
    header = lines[0].rstrip().split(sep)
    num_fields = len(header)
    for line in lines[1:]:
        a = line.rstrip().split(sep)
        LL.append(dict(zip(header, a)))
        
    return LL

def parse_msp_to_listdict(file, field_separator=': ', return_peaks=False):
    '''
    parse a msp file within reasonable size.
    return list of dictionaries, [{key, value pairs}, ...]
    '''
    features = []
    w = open(file).read().rstrip().split('\n\n')
    for block in w:
        d = {}
        lines = block.splitlines()
        data = []
        for line in lines:
            if field_separator in line:
                x = line.split(field_separator)
                d[x[0]] = x[1]
            elif line.strip():
                data.append(line.split())
        if return_peaks:
            try:
                d['peaks'] = [(float(x[0]), float(x[1])) for x in data]
            except(TypeError):
                print("Failed to convert peaks: ", block)
        features.append(d)
        
    return features

def archive_expt_features(short_name, num_samples, mass_accuracy_ratio, list_features, outdir):
    '''
    Record 'id_number', 'mz', 'rtime', 'SNR', 'detection_counts' etc in text.
    '''
    use_cols = [ 'id_number', 'mz', 'rtime', 'parent_masstrack_id', 
                    'peak_area', 'cSelectivity', 'goodness_fitting', 'snr', 'detection_counts' ] 
    s = '#num_samples: ' + str(num_samples) + '\n'
    s += '#mass_accuracy_ratio: ' + str(mass_accuracy_ratio) + '\n'
    s += '\t'.join(use_cols) + '\n'
    
    for f in list_features:
        s += '\t'.join([str(f[ii]) for ii in use_cols
            ]) + '\n'
        
    with open(os.path.join(outdir, short_name+"_feature_archive.tsv"), 'w') as O:
        O.write(s)
        
def check_good_peak(peak, snr=5, shape=0.9, area=1e6):
    '''
    Decide if peak is good quality based on SNR and peak shape.
    Not using peak area now.
    '''
    if peak['snr'] > snr and peak['goodness_fitting'] > shape:
        # Not using peak['peak_area'] > area:
        return True
    else:
        return False

# extract list_features from empCpd JSON file
def epd2featurelist(list_epds):
    list_features = []
    for epd in list_epds:
        list_features += epd["MS1_pseudo_Spectra"]

    return list_features

def get_feature2epd_dict(list_epds, mode):
    _d = {}
    for epd in list_epds:
        primay_feature = get_primary_epd_feature(epd, mode)
        for f in epd["MS1_pseudo_Spectra"]:
            _d[f["id"]] = {
                'parent_epd_id': epd['interim_id'],
                'neutral_formula_mass': epd['neutral_formula_mass'],
                'ion_relation': f['ion_relation'],
                # establish primary ion feature here, i.e. link btw this and primary feature
                'primary_ion_feature': primay_feature,
                'mz': f['mz'],
                'rtime': f['rtime'],
                'snr': f['snr'],
                'is_good_peak': f['is_good_peak'],
                'detection_counts': f['detection_counts']
            }
    return _d


def get_primary_epd_feature(epd, mode='pos'):
    '''
    Get feature ID of the primary ion in this epd.
    '''
    if mode == 'pos':
        primary_ions = primary_ions_pos_ordered
    elif mode == 'neg':
        primary_ions = primary_ions_neg_ordered
    else:
        raise ValueError("ion mode error.")
    
    _primary = [f for f in epd["MS1_pseudo_Spectra"] if f['ion_relation'] in primary_ions]
    if _primary:
        _ordered = [
            (primary_ions.index(f['ion_relation']), f['id']) for f in _primary
        ]
        return sorted(_ordered)[0][1]   # first f['id']
    else:
        return ''
    

def get_epd_stats(list_epds, natural_ratio_limit=0.5):
    '''
    Get numbers of khipus and singletons, and isopairs and good khipus.
    Good khipu = isopair and M0 being a good feature.
    Has to run after epd2featurelist_from_file.
    Returns dict.
    '''
    khipus_isopairs, num_isopair_mtracks, good_khipus = get_isopairs_good_khipus(
                                list_epds, natural_ratio_limit)

    # M0 must be always the first item in MS1_pseudo_Spectra
    return {
        'num_empcpds': len(list_epds),
        'num_khipus_isopairs': len(khipus_isopairs),
        'num_isopair_mtracks': num_isopair_mtracks,
        'num_good_khipus': len(good_khipus),
        'num_singletons': count_singletons(list_epds),
        'good_khipus': good_khipus
    }

def epd2featurelist_from_file(file, mode='pos', snr=5, shape=0.9):
    '''
    returns list_features, feature2epd_dict, epd_summary
    '''
    list_epds = json.load(open(file))
    for epd in list_epds:
        for f in epd["MS1_pseudo_Spectra"]:
            # fix typing
            f['snr'] = float(f['snr'])
            f['goodness_fitting'] = float(f['goodness_fitting'])
            f['peak_area'] = float(f['peak_area'])
            f['is_good_peak'] = check_good_peak(f, snr, shape)
            
    return epd2featurelist(list_epds), get_feature2epd_dict(list_epds, mode), get_epd_stats(list_epds)


def read_master_datasets_records(infile='r1_datasets_stats.tsv', sep='\t'):
    '''
    returns
    dict of namedTuples
    '''
    Dataset = namedtuple('Dataset', ['dataset_id_int',
                                     'feature_table_id',
                                    'mode',
                                    'chromatography',
                                    'num_samples',
                                    'num_features',
                                    'num_good_features',
                                    'num_empcpds',
                                    'num_khipus_isopairs',
                                    'num_isopair_mtracks',
                                    'num_good_khipus',
                                    'num_singletons',
                                    'mz_calibration_ratio',
                                    'num_csm_matched'])
    d = {}
    for dts in map(Dataset._make, csv.reader(open(infile, "r"), delimiter=sep)):
        d[dts.feature_table_id] = dts
    _ = d.pop('feature_table_id') # remove header item
    return d


def short_vtrack_summary(d):
    # ['id', 'numb_good_features', 'median_numb_hilic_features_perdataset', 'median_num_rp_features_perdataset', 'num_all_features',
    # 'num_datasets', 'num_hilic_datasets', 'num_rp_datasets',  'num_studies', 'max_good_features_per_method',
    # 'is_primary_track', d['primary_ion_list'
    # ]
    return (
        d['id'], d['numb_good_features'], d['median_features_HILIC'], d['median_features_RP'], 
        d['num_all_features'], d['num_datasets'], d['number_datasets_HILIC'], d['number_datasets_RP'], d['num_studies'], 
        d['max_good_features_per_method'],
        d['is_primary_track'], d['primary_ion_list'],
    )

def export_table_by_fields(summary_tally, fields, outfile='short_vtrack_summary_vTracks_consensus_orbi_pos.tsv'):
    # export summary
    s = '\t'.join(fields) + '\n'

    for x in summary_tally:
        s += '\t'.join([
            str(ii) for ii in short_vtrack_summary(x)
        ]) + '\n'
        
    with open(outfile, 'w') as O:
        O.write(s)



def export_json_consensus_vtracks(summary_tally, dataset_dict, outfile='vTracks_consensus_summary_pos_r1.json'):
    '''
    summary_tally : e.g. [{
                        "id": "r1_pos_70.073284",
                        "is_primary_track": true,
                        "num_all_features": 469,
                        "numb_good_features": 119,
                        "num_datasets": 68,
                        "num_studies": 28,
                            "number_datasets_HILIC": 25,
                        "median_features_HILIC": 1.0,
                        "stdev_features_HILIC": 0.574108003776293,
                        "number_datasets_RP": 43,
                        "median_features_RP": 1.0,
                        "stdev_features_RP": 1.7415483931952611,
                        "max_good_features_per_method": 9,
                        "ions_stats": [
                        [
                            74,
                            ""
                        ],
                        [
                            35,
                            "13C/12C,M+H+"
                        ],
                        [
                            7,
                            "M0,M+H+"
                        ],
                        [
                            2,
                            "M0,M+H+, 3x charged"
                        ],
                        [
                            1,
                            "M0,ACN"
                        ]
                        ],
                        "features_per_dataset": [
                        [
                            "ST002016_RPpos_B1_ppm5_351391",
                            2
                        ], ...}, ...]

    outfile : name of JSON file to export the list of cmTracks
    dataset_dict : map dataset id to full name of each source dataset.
    
    '''
    def _reformat(max_dataset):
        new = {}
        fields = ['id', 'rtime', 'parent_epd_id', 'snr']
        if max_dataset:
            new['dataset_id'] = max_dataset[0]['dataset_id']
            new['variables'] = fields
            new['features'] = [[x[ii] for ii in fields] for x in max_dataset]
        return new
        
    new = []
    # deepcopy not to modify input data
    for cmt in summary_tally:
        tmp = copy.deepcopy(cmt) 
        tmp['features_per_dataset'] = [(dataset_dict[x[0]], x[1]) for x in cmt['features_per_dataset']] 
        tmp['max_dataset'] = _reformat(cmt['max_dataset'])
        new.append(tmp)
            
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(new, f,  ensure_ascii=False, indent=2)  # cls=NpEncoder,



#
# Annotation related
#
# khipu needs update on better calculation of multi-charged ions. But this is a start:
#

def build_consensus_masslist(List1, List2, ppm=1):
    '''
    To build a list of consensus mass values from two input lists.
    Redundance is expected within and in between lists. 
    
    example input:
        [('r1_pos_191.04317_RP_1', 168.05387000000002),
        ('r1_pos_191.044982_HILIC_0', 168.055682),
        ('r1_pos_191.044982_RP_0', 168.055682),
        ('r1_pos_191.048392_HILIC_0', 190.041112),
        ('r1_pos_191.048392_RP_0', 190.041112),
        ('r1_pos_191.053244_RP_0', 168.063944)]
    returns
        Unique list of mass values, and corresponding input feature IDs. 
        
    ref: asari.mass_functions.bin_by_median
    '''
    def tol(x, ppm=ppm):
        return x * 0.000001 * ppm
    
    List_of_tuples = sorted([(x[1], x[0]) for x in List1 + List2])
    
    new = [[List_of_tuples[0], ], ]
    for X in List_of_tuples[1:]:
        if X[0]-np.median([ii[0] for ii in new[-1]]) < tol(X[0]):       
            # median moving with list change
            new[-1].append(X)
        else:
            new.append([X])
    
    return new


def calculate_neutral_mass_pos(mz, ion):
    '''
    ion in ['M0,Na/H', 'M0,M+H+', 'M0,M+H+, 2x charged', 'M0,Na/H, 2x charged', 
                        'M0,M+H+, 3x charged', 'M0,Na/H, 3x charged']
    '''
    if ion == 'M0,M+H+':
        return mz - 1.00728
    elif ion == 'M0,Na/H':
        return mz - 22.9893
    elif ion == 'M0,M+H+, 2x charged':
        return (mz - 1.00728) * 2
    elif ion == 'M0,M+H+, 3x charged':
        return (mz - 1.00728) * 3
    elif ion == 'M0,Na/H, 2x charged':
        return (mz - 22.9893) * 2
    elif ion == 'M0,Na/H, 3x charged':
        return (mz - 22.9893) * 3
    else:
        return None
    
    
def calculate_neutral_mass_neg(mz, ion):
    '''
    ion in [["M0,Na/H", "M0,M-H-", "M0,Na/H, 2x charged", "M0,M-H-, 2x charged", 
                        "M0,Na/H, 3x charged", "M0,M-H-, 3x charged"]]
    Some not so frequent - maybe adjust later.
    '''
    if ion == 'M0,M-H-':
        return mz + 1.00728
    elif ion == 'M0,Na/H':
        return mz - 21.9820
    elif ion == 'M0,M-H-, 2x charged':
        return (mz + 1.00728) * 2
    elif ion == 'M0,M-H-, 3x charged':
        return (mz + 1.00728) * 3
    elif ion == 'M0,Na/H, 2x charged':
        return (mz - 21.9820) * 2
    elif ion == 'M0,Na/H, 3x charged':
        return (mz - 21.9820) * 3
    else:
        return None
    

# mass calibration
def custom_mz_calibrate(features, 
                        kcd_instance, 
                        mode='pos', 
                        mz_tolerance_ppm=30, 
                        required_calibrate_threshold=0.000002):
    '''
    mass calibration, better done after khipu because this may shift all m/z values in features.
    Modified from asari.experiment. 
    mz_tolerance_ppm : leave a bigger window only to prefilter data points into the regression.
    '''
    mass_accuracy_ratio = kcd_instance.evaluate_mass_accuracy_ratio(
        [f['mz'] for f in features], 
        mode=mode, mz_tolerance_ppm=mz_tolerance_ppm)
    
    if mass_accuracy_ratio:
        if abs(mass_accuracy_ratio) > required_calibrate_threshold:
            print("Mass shift is greater than %2.1f ppm. Correction applied." 
                    %(required_calibrate_threshold*1000000))
            _correction = mass_accuracy_ratio + 1
            for F in features:
                F['mz'] = F['mz'] / _correction
                F['mz_corrected_by_division'] = _correction
        else:
            print("No mz correction is needed.")
    else:
        print("Mass accuracy check is skipped, mass_accuracy_ratio computing unsuccessful.")
    
    return mass_accuracy_ratio, features
    

def features2cpds(features, 
                    mode='pos', mz_tolerance_ppm=5,
                    rt_tolerance=2,  # seconds
                    outfile="exported_empCpds.json"):
    '''
    Export to JSON empCpds plus singletons
    '''
    EED = ExperimentalEcpdDatabase(mode=mode, 
                                mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=rt_tolerance)
    # passing patterns from .default_parameters
    if mode == 'pos':
        EED.adduct_patterns = adduct_search_patterns
    else:
        EED.adduct_patterns = adduct_search_patterns_neg
    EED.isotope_search_patterns = isotope_search_patterns
    EED.extended_adducts = extended_adducts
    
    # First khipu organized empCpds
    EED.build_from_list_peaks(features)
    
    new_id_start = len(EED.dict_empCpds)
    singletons = [p for p in EED.dict_peaks if p not in EED.peak_to_empCpd.keys()]
    for p in singletons:
        new_id_start += 1
        interim_id = '_singleton_' + str(new_id_start)
        feature = EED.dict_peaks[p]
        feature.update({ "isotope": "M0",
            "modification": "",
            "ion_relation": "",}
            )
        EED.dict_empCpds[interim_id] = {'interim_id': interim_id,
                'neutral_formula_mass': None, 'neutral_formula': None,
                "Database_referred": [],
                "identity": [],
                'MS1_pseudo_Spectra': [ feature ],    
        }
    
    # export empCpds here to JSON, including singletons.
    EED.export_empCpds(outfile)
    
def custom_kcd_annotate(full_list_empCpds, kcd_instance, 
                        mode='pos', mz_tolerance_ppm=10,
                        rt_tolerance=2,  # seconds
                        ):
    '''
    Get annotation on full_list_empCpds(list_empCpds plus singletons) using kcd_instance.
    Summary on matches to kcd_instance by C12/13 isopairs and by all (empCpds + singleton).
    list_empCpds are prebuilt so that the result can be compared on different KCD instances.
    
    EED.extend_empCpd_annotation() can handle singletons.
    '''
    EED = ExperimentalEcpdDatabase(mode=mode, 
                                mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=rt_tolerance)
    EED.build_from_list_empCpds(full_list_empCpds)
    # This matches empCpds to kcd_instance
    EED.extend_empCpd_annotation(kcd_instance)

    return EED.dict_empCpds



def line2cpd_dict(header, line):
    a = line.split('\t')
    d = dict(zip(header, a))
    if 'accession' in d:
        d['primary_id'] = d['accession']
    elif 'id' in d:
        d['primary_id'] = d['id']
    else:
        d['primary_id'] = line[:20]     # hack 
        
    if d['monisotopic_molecular_weight']:
        d['neutral_formula_mass'] = float(d['monisotopic_molecular_weight'])
        d['neutral_formula'] = d['chemical_formula']
    return d

def build_KCD_from_formula_coordinate(formula_coordinate):
    list_compounds_formula_coordinate = []
    for x in formula_coordinate:
        list_compounds_formula_coordinate.append(
            {
                'primary_id': x[1],
                'primary_db': x[2],
                "neutral_formula": x[1],
                "neutral_formula_mass": x[0], 
            }
        )
        
    KCD_formula_coordinate = knownCompoundDatabase()
    KCD_formula_coordinate.mass_index_list_compounds(list_compounds_formula_coordinate)
    KCD_formula_coordinate.build_emp_cpds_index()
    
    return KCD_formula_coordinate

def build_KCD_from_hmdb_tsv(tsv_file):
    # Build KCD HMDB5 serum subset
    # serum = '/Users/lish/li.proj/Resources/HMDB-5/parsed_serum_metabolites.tsv'
    table = open(tsv_file).readlines()                                                                                    
    header = table[0].split('\t')   

    cpds = []
    for line in table[1:]:
        c = line2cpd_dict(header, line)
        if 'neutral_formula_mass' in c and 'monisotopic_molecular_weight' in c:
            cpds.append(c)
            
    print(len(cpds), cpds[55])
    print(header)

    KCD_x = knownCompoundDatabase()
    KCD_x.mass_index_list_compounds(cpds)
    KCD_x.build_emp_cpds_index()
    
    return KCD_x

def build_KCD_from_GEM(gem_json):
    '''
    Build KCD humanGEM
    Returns
    KCD instance
    
    [print(x) for x in [
        M['meta_data'],
        len(M['list_of_reactions']),
        M['list_of_reactions'][222],
        len(M['list_of_compounds']),
        M['list_of_compounds'][300], ]
    ]

    '''
    M = json.load(open(gem_json))
    new = []
    for x in M['list_of_compounds']:
        if x['neutral_mono_mass']:
            x['neutral_formula_mass'] = x['neutral_mono_mass']
            x['primary_id'] = x['id']
            new.append(x)

    KCD_gem = knownCompoundDatabase()
    KCD_gem.mass_index_list_compounds(new)
    KCD_gem.build_emp_cpds_index()
    
    return KCD_gem


def custom_export_peak_annotation(dict_empCpds, kcd_instance, export_file_name):
    '''
    Export feature annotation to tab delimited tsv file, where interim_id is empCpd id; and JSON.
    modified from asari.experiment 
    
    Parameters
    ----------
    dict_empCpds : dict
        dictionary of empirical compounds, using interim_id as key, as seen in JMS.
    kcd_instance : KnownCompoundDatabase instance
        the known compound database that was used in annotating the empirical compounds.
    export_file_name : str
        to used in output file name.
        
    Usage
    -----
    full_list_empCpds  = json.load(open(json_empcpd))
    dict_empCpds = custom_kcd_annotate(full_list_empCpds, KCD_hmdb, mode='pos', mz_tolerance_ppm=10,
                            rt_tolerance=2,)  # seconds)
    custom_export_peak_annotation(dict_empCpds, KCD_hmdb, 
                            export_file_name='rppos407_KCD_hmdb_kcd_annotate_result.tsv')
    
    '''
    s = "[peak]id_number\tmz\trtime\tapex(scan number)\t[EmpCpd]interim_id\
        \t[EmpCpd]ion_relation\tneutral_formula\tneutral_formula_mass\
        \tname_1st_guess\tmatched_DB_shorts\tmatched_DB_records\n"
        
    match_dict = {}
    
    for _, V in dict_empCpds.items():
        name_1st_guess, matched_DB_shorts, matched_DB_records = '', '', ''
        if 'list_matches' in V:
            list_matches = V['list_matches']
            if list_matches:
                match_dict[V['interim_id']] = [x['primary_id'] for x in 
                                               kcd_instance.mass_indexed_compounds[list_matches[0][0]]['compounds']]
                
                name_1st_guess = kcd_instance.mass_indexed_compounds[list_matches[0][0]]['compounds'][0]['name']
                matched_DB_shorts = ", ".join([ "(" + kcd_instance.short_report_emp_cpd(xx[0]) + ")"  for xx in list_matches])
                matched_DB_records = ", ".join([str(xx) for xx  in list_matches])

        for peak in V['MS1_pseudo_Spectra']:
            s += '\t'.join([str(x) for x in [
                peak['id_number'], peak['mz'], peak['rtime'], peak['apex'], 
                V['interim_id'], peak.get('ion_relation', ''),
                V['neutral_formula'], V['neutral_formula_mass'],
                name_1st_guess, matched_DB_shorts, matched_DB_records]]) + "\n"

    with open(export_file_name, encoding='utf-8', mode='w') as O:
        O.write(s)
        
    # export JSON 
    outfile = export_file_name.replace('.tsv', '') + '.json'
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(match_dict, f,  ensure_ascii=False, indent=2)  # cls=NpEncoder,
        
    print("\nAnnotation of %d Empirical compounds was written to %s.\n\n" %(len(dict_empCpds), export_file_name))


def extract_userAnnotation(indir, infile, 
                           dict_datasets_int_id,
                           dict_f2csmf,
                           ):
    '''
    infile : JSON file annotation to a raw feature.
    e.g.     "F4730": {
                    "name": "3-phenylpropionate (hydrocinnamate)",
                    "id": "3-phenylpropionate (hydrocinnamate)",
                    "platform": "LC/MS Neg",
                    "mz": 149.0608,
                    "ri": 2860.0,
                    "rtime": 2860.0,
                    "superclass": "Xenobiotics",
                    "subclass": "Benzoate Metabolism",
                    "cas": null,
                    "pubchem": "107",
                    "hmdb": "HMDB00764",
                    "kegg": "C05629",
                    "monoisotopic_mass": 150.068079564
                },
    
    returns {csmf: [(study_id, cpd_dict), (), ...]}
    
    '''
    new = {}
    _int_id_ = dict_datasets_int_id[infile.split('.')[0]]
    dict_anno = json.load(open(os.path.join(indir, infile)))
    for k,v in dict_anno.items():
        _k_ = str(_int_id_) + '_' + k
        matched_csmf = dict_f2csmf.get(_k_)
        if matched_csmf:
            if matched_csmf in new:
                new[matched_csmf].append(
                    (_int_id_, v)
                )
            else:
                new[matched_csmf] = [(_int_id_, v)]
    print(
        "Number of annotated raw features: ", len(dict_anno), 
        "\n",
        "Number of annotated CSM features: ", len(new)
    )
    return new



#
# Get empCpds by filtering 13C/12C pattern
#

def get_feature_of_max_intensity(featureList):
    '''
    To get feature of max intensity, and avoid errors by sorting.
    e.g. sorted(M0, reverse=True)[0][1] leas to `TypeError: '<' not supported between instances of 'dict' and 'dict'`
    Use np.argmax here, which is okay with ties.
    '''
    ints = [f['representative_intensity'] for f in featureList]
    idx = np.argmax(ints)
    return featureList[idx]

def get_M0(MS1_pseudo_Spectra):
    '''returns M0 feature with highest representative_intensity.
    Without verifying which ion form.'''
    M0 = [f for f in MS1_pseudo_Spectra if f['isotope']=='M0']
    if M0:
        return get_feature_of_max_intensity(M0)
    else:
        return []
    
def get_M1(MS1_pseudo_Spectra):
    '''returns M+1 feature with highest representative_intensity.
    Without verifying which ion form.'''
    M = [f for f in 
          MS1_pseudo_Spectra if f['isotope']=='13C/12C']
    if M:
        return get_feature_of_max_intensity(M)
    else:
        return []
     
def get_highest_13C(MS1_pseudo_Spectra):
    '''returns 13C labeled feature with highest representative_intensity.
    Without verifying which ion form. Because the label goes with sepecific atoms depending on pathway.
    '''
    M = [f for f in 
          MS1_pseudo_Spectra if '13C/12C' in f['isotope']]
    if M:
        return get_feature_of_max_intensity(M)
    else:
        return []

    
def filter_khipus(list_empCpds, natural_ratio_limit=0.5):
    '''
    returns 
    isopair_empCpds = with good natural 13C ratio, based on M1/M0, not checking adduct form.
 
    Usage
    -----
    full_list_empCpds  = json.load(open(json_empcpd))
    isopair_empCpds = filter_khipus(full_list_empCpds)

    '''
    isopair_empCpds = []
    for epd in list_empCpds:
        # interim_id = v['interim_id']
        M0, M1 = get_M0(epd['MS1_pseudo_Spectra']), get_M1(epd['MS1_pseudo_Spectra'])
        if M0 and M1:
            if float(M1['representative_intensity'])/(1 + float(M0['representative_intensity'])) < natural_ratio_limit:
                isopair_empCpds.append( epd['interim_id'] )

    return isopair_empCpds
    
    
def get_isopairs_good_khipus(list_empCpds, natural_ratio_limit=0.5):
    '''
    returns 
    Two lists of khipus, isopair_empCpds_ids (IDs only), good_khipus (full dict), and number of isopair_mtracks.
    isopair_empCpds = with good natural 13C ratio, based on M1/M0, not checking adduct form.
    good_khipus = isopair_empCpds and M0 being a good feature.
    
    Some inline MS/MS expts cause split MS1 peaks. Thus isopair_mtracks are more indicative of the data coverage.
    
    Usage
    -----
    full_list_empCpds  = json.load(open(json_empcpd))
    isopair_empCpds, num_isopair_mtracks, good_khipus = get_isopairs_good_khipus(full_list_empCpds)

    Use filter_khipus if not considering is_good_peak.
    '''
    isopair_empCpds_ids, isopair_mtracks, good_khipus = [], [], []
    for epd in list_empCpds:
        # interim_id = v['interim_id']
        M0, M1 = get_M0(epd['MS1_pseudo_Spectra']), get_M1(epd['MS1_pseudo_Spectra'])
        if M0 and M1:
            if float(M1['representative_intensity'])/(1 + float(M0['representative_intensity'])) < natural_ratio_limit:
                isopair_empCpds_ids.append( epd['interim_id'] )
                # if epd["MS1_pseudo_Spectra"][0]['is_good_peak']: 
                if M0['is_good_peak']: # not assuming first feature
                    good_khipus.append( epd )
                    isopair_mtracks.append( M0['parent_masstrack_id'] )

    return isopair_empCpds_ids, len(set(isopair_mtracks)), good_khipus
    
    
def count_singletons(list_empCpds):
    return len([epd for epd in list_empCpds if len(epd['MS1_pseudo_Spectra'])==1])
    
    
def get_isopair_features(full_list_empCpds, isopair_empCpds):
    '''
    Not clean, including more features than desired..
    '''
    isopair_features = []
    for epd in full_list_empCpds:
        if epd['interim_id'] in isopair_empCpds:
            isopair_features += epd['MS1_pseudo_Spectra']

    return isopair_features
    


def compare(list1, list2, mz_ppm=5, rt_tolerance=1e9):
    '''compare matches and print unmatched in list1. 
    rt_tolerance can be relaxed to avoid RT comparison when different methods are not compatible.
    
    To extract matched from list1, use mf.list_match_lcms_features.
    '''
    print("\n  Best match comparisons:")
    valid_matches, dict1, dict2 = mf.bidirectional_best_match(list1, list2, mz_ppm, rt_tolerance)

    print("\n  List based inclusive comparisons:")
    dict1, dict2 = mf.bidirectional_match(list1, list2, mz_ppm, rt_tolerance)

    unmatched = [p for p in list1 if p['id'] not in dict1]
    print("\n\nUnmatched features ****** ", len(unmatched), "*******\n")
    # [p for p in list1 if p['id'] not in [x[0] for x in valid_matches]]
    # print(unmatched)
    
def compare2(list1, list2, mz_ppm=5, rt_tolerance=1e9):
    '''
    No use other than wrapping mf.bidirectional_match
    '''
    print("\n  List based inclusive comparisons:")
    dict1, dict2 = mf.bidirectional_match(list1, list2, mz_ppm, rt_tolerance)
    return dict1, dict2 

#
# KDE functions
#

def get_kde_peaks(x_kde_support, y_kde_density, 
                  height=0.5,
                  distance=10,
                  prominence=0.25,
                  width=10,
                  wlen=200,
                  ):
    '''
    Find peaks in KDE density.
    Returns list [(m/z peak, kde density), ...]
    
    We fit KDE with bandwith of 1 ppm.
    
    distance : 10, # about 1 ppm on this dataset
    prominence : half of min height
    wlen : about 20 ppm
    '''
    peaks, properties = find_peaks(y_kde_density, 
                                    height=height, 
                                    distance=distance,
                                    prominence=prominence,
                                    width=width, 
                                    wlen=wlen,
                                    ) 
    real_apexes = [x_kde_support[ii] for ii in peaks]
    return list(zip(real_apexes, properties['peak_heights']))


def extract_segment_kde_peaks(sorted_mz_values,
                                  distance=10, width=10, wlen=200,
                                  relative_height=0.1,
                              ):
    '''
    Run KDE and KDE peak detection for a large list (e.g. 18912000 data points in posList) sorted_mz_values.
    Returns list [(m/z peak, kde density), ...]
    
    This only establishes the consensus peaks as features.
    The density value is relative. 
    We should remap each dataset to the consensus peaks to get reliable stats.
    Preannotation is done within each dataset too.
    '''
    mz_partitions = list(range(50, 400, 50)) + list(range(400, 2100, 100))
    result = []
    
    for ii in range(len(mz_partitions)-1):
        mzstart, mzend = mz_partitions[ii], mz_partitions[ii+1]
        mzdata = [x for x in sorted_mz_values if mzstart<x<mzend]
        bandwidth = mzstart/1e6

        kde = sm.nonparametric.KDEUnivariate(mzdata)
        kde.fit(bw=bandwidth) 

        threshold = relative_height * kde.density.max()
        # len(mzdata)/ kde.support.shape[0]
        prominence = 0.5 * threshold

        peaks_density = get_kde_peaks(kde.support, kde.density, height=threshold, 
                                            distance=distance,
                                            prominence=prominence,
                                            width=width, 
                                            wlen=wlen)
        
        print(mzstart, len(mzdata), kde.support.shape[0], threshold, len(peaks_density))
        result += peaks_density
        
    return result

# --------------------------------------------------------------------------------
# Similar functions below now in extDataModels
#

def summarize_ions(cfeature):
    '''
    returns dict of ions in a tally_consensus_feature
    '''
    ion_relations = [f['ion_relation'] for f in cfeature]
    _d = []
    for ion in set(ion_relations):
        _d.append( (ion_relations.count(ion), ion) )
    _d.sort(reverse=True)
    
    return _d


def summarize_vtrack(id_vtrack, vtrack, mode='pos', filtering=True):
    '''
    Calculate per ref virtual track: #features, #datasets, #studies, #epds, #max_features_per_method
    vtrack : virtual track, defined by KDE peaks in cumulative data. Renamed to consensus mass track later.
    filtering : use good peaks only for stats.
    vtrack example:   tally_consensus_features_pos['pos_100.028708'] = [
                        {'id': 'F1705',
                        'dataset_id': 'ST001930_RPpos_B16_ppm5_24125850',
                        'parent_epd_id': 'kp1939_198.5473',
                        'neutral_formula_mass': 198.54731080027364,
                        'ion_relation': 'M0,M+H+, 2x charged'}, 
                        ...]
    '''
    num_features, numb_good_features = len(vtrack), None
    if filtering:
        vtrack = [f for f in vtrack if f['is_good_peak']]
        numb_good_features = len(vtrack)
        
    _dataset2features = {}
    for f in vtrack:
        if f['dataset_id'] in _dataset2features:
            _dataset2features[f['dataset_id']].append(f)
        else:
            _dataset2features[f['dataset_id']] = [f]
            
    # feature per dataset is an estimate of isomers
    features_per_dataset = [(k,len(v)) for k,v in _dataset2features.items()]
    
    ions = summarize_ions(vtrack)
    
    HILIC, RP = [], []
    max_good_features_per_method, max_dataset = 0, None
    for k,v in _dataset2features.items():
        if 'RP' in k:
            RP.append(len(v))
        elif 'HILIC' in k:
            HILIC.append(len(v))
            
        good_features = [f for f in v if f['is_good_peak']] # redundant but independent from filtering
        if len(good_features) > max_good_features_per_method:
            max_good_features_per_method = len(good_features)
            max_dataset = k
    
    is_primary_track = False
    primary_ion_list, alternative_ion_list = get_epds_determine_primary_track(
        vtrack, max_good_features_per_method, ions, mode=mode)
    if primary_ion_list:
        is_primary_track = True
    
    studies = set([f['dataset_id'].split('_')[0] for f in vtrack])
    if max_dataset:
        max_dataset = _dataset2features[max_dataset]
    else:
        max_dataset = []
    
    return {
        'id': id_vtrack,
        'is_primary_track': is_primary_track,
        'num_all_features': num_features,
        'numb_good_features': numb_good_features,
        'num_datasets': len(_dataset2features),
        'num_studies': len(studies),
        'features_per_dataset': features_per_dataset,
        'number_datasets_HILIC': len(HILIC),
        'median_features_HILIC': np.median(HILIC),
        'stdev_features_HILIC': np.std(HILIC),
        'number_datasets_RP': len(RP),
        'median_features_RP': np.median(RP),
        'stdev_features_RP': np.std(RP),
        'max_good_features_per_method': max_good_features_per_method,
        # 'ions_stats': ions,
        'max_dataset': max_dataset,
        'primary_ion_list': primary_ion_list, 
        'alternative_ion_list': alternative_ion_list,
    }
   

def get_epds_determine_primary_track(vtrack, max_good_features_per_method, ions, mode):
    '''Determine if this is a primary track based on if the top ions contain any primary ion forms.
    Get matched best feature (highest SNR) to get neutral_formula_mass in its linked empCpd, each per ion form.
    Returns primary_list (with best_features), alternatives (remaining top ions, including '')
    primary_list will be [] if not a primary track
    
    ions : sorted list of counts of ion forms, e.g. 
                        [(438, ''),
                        (351, 'M0,Na/H'),
                        (158, 'M0,M+H+'),
                        (42, 'M0,HCl'),
                        (28, '13C/12C, 3x charged,M+H+, 3x charged'),
                        (17, 'M0,K/H'),
                        (12, 'M0,M+H+, 2x charged'),
                        (4, 'M0,Na/H, 2x charged'),
                        (3, 'M0,HCl, 2x charged'),
                        (2, 'M0,ACN'),
                        (1, 'M0,M+H+, 3x charged'),
                        (1, '13C/12C,K/H'),
                        (1, '13C/12C,ACN'), ...]
    '''
    if mode == 'pos':
        primary = primary_ions_pos
    elif mode == 'neg':
        primary = primary_ions_neg
        
    primary_list, alternatives = [], []
    for x in ions[:max_good_features_per_method]:
        if x[1] in primary:
            _match = [(f['snr'], str(f), f) for f in vtrack if f['is_good_peak'] and x[1] == f['ion_relation']]  # str f as tiebreaker
            primary_list.append(
                [x[1], x[0], sorted(_match, reverse=True)[0][2]]
            )
        else:
            alternatives.append(
                [x[1], x[0]]
            )
    
    return primary_list, alternatives
    
    

