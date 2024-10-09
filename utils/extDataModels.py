'''
pre-annotation algorithms implemented in class cmTrack.
How to get ion annotations from many source datasets by popular votes.
'''
import pickle
from itertools import permutations
import numpy as np

# from mining import primary_ions_pos_ordered, primary_ions_neg_ordered

def no_tie_int_median(a):
    '''
    Parameter a : input as a list, not array.
    Used to get consensus number of features on a mass registry. We don't want tied median here.
    Returns the value above median if len(a) is even; force int.
    
    Note :
    1. np.percentile changed argument, not to use to avoid version confusion.
    2. isinstance(m, int) doesn't work. This requires checking odd or even number as len(a), 
    since even number can return an integer too.
    '''
    if len(a) % 2 == 0: # even number
        return int(np.median( [-1] + a ))
    else :
        return int(np.median(a))


class extended_Dataset:
    '''
    # Placeholder, a conceptual class for a dataset.
    See mining.read_master_datasets_records to get a dict of namedTuples for this purpose.
        
    desired:
        'feature_table_id',
        'mode',
        'chromatography',
        'num_samples',
        'num_features',
        'num_good_features',
        'mz_calibration_ratio',
        'num_empcpds',
        'num_khipus_isopairs',
        'num_isopair_mtracks',
        'num_good_khipus',
        'num_singletons',
        'good_khipus',
        'num_features_matched_csm'
    '''
    def __init__(self):
        id = ''
        feature_table_id = ''
        mode = 'pos'   # or 'neg'
        chromatography = 'HILIC'   # or 'RP' etc
        
        num_samples = 0
        num_features = 0
        num_good_features = 0   # defined by e.g. SNR > 5, goodness_fitting > 0.9
        mz_calibration_ratio = 0
        
        num_khipus = 0
        num_khipus_isopairs = 0
        num_singletons = 0
        
    def export_json(self):
        pass
    
    
class cmRegistry:
    '''
    A consensus mass registry from large-scale data mining.
    Unknown number of compounds may exist on a cmTrack.
    The initial cmRegistry is established by KDE peaks, either pos or neg ion mode. 
    We keep list of good khipus here to facilitate preannotation.
    This implements popular votes to get represnetative empCpds on each chromatography and ionization mode.
    
    '''
    def __init__(self, 
                 id, 
                 mode, 
                 matched_features, 
                 dict_datasets_int_id,
                 chromatography_types=['HILIC', 'RP']
                 ):
        '''
        After a cmTrack is determined (via KDE), all matched features are input here and analyzed.
        matched_features are from tally_consensus_features_neg/pos.
        We also need an external reference of good khipus from all datasets.
        chromatography_types are 'HILIC' and 'RP' for now.
        KDE value is not included here, but can be looked up elsewhere.
        '''
        self.id = id        # e.g., 'r1_pos_302.232463'; m/z is coded in id.
        self.mode = mode    # 'pos' or 'neg'. 'neu' will be a separate class
        
        self.matched_features = matched_features
        self.good_features = [f for f in matched_features if f['is_good_peak']]

        self.chromatography_types = chromatography_types
        self.dict_datasets_int_id = dict_datasets_int_id
        self.is_primary_track = None
        self.num_all_features = len(self.matched_features)
        self.numb_good_features = len(self.good_features)
        self.dict_datasets = {} # separate data by methods
        self.dict_summary = {}
        self.dict_csmf_members = {}
        
    def summarize_features(self, datasets_meta_dict, concise=False):
        '''
        We get list of matched features to a cmTrack, e.g. [   
            {'id': 'F36583',
            'rtime': 2259.24,
            'is_good_peak': False,
            'snr': 2.0,
            'dataset_id': 'ST002845_RPpos__ppm5_3505225',
            'parent_epd_id': '_singleton_48607',
            'neutral_formula_mass': None,
            'ion_relation': ''}, ...,
            {'id': 'F1705',
            'dataset_id': 'ST001930_RPpos_B16_ppm5_24125850',
            'parent_epd_id': 'kp1939_198.5473',
            'neutral_formula_mass': 198.54731080027364,
            'ion_relation': 'M0,M+H+, 2x charged'}, 
            ...]
        
        This gets stats on ions and empCpds on each chromatography_group.
        '''
        self.ions = self.summarize_ions(self.matched_features)
        self.dict_datasets = self.group_data_by_methods(datasets_meta_dict)
        self.summarize_chromatography_groups(concise)

    def group_data_by_methods(self, datasets_meta_dict):
        '''
        Use good_features only.
        
        datasets_meta_dict['MTBLS136_RPpos_B19_ppm5_3581551']:
            Dataset(feature_table_id='MTBLS136_RPpos_B19_ppm5_3581551', mode='pos', 
            chromatography='RP', num_samples='144', num_features='29442', num_good_features='13338', 
            num_empcpds='22235', num_khipus_isopairs='1533', num_isopair_mtracks='906', num_good_khipus='1058', num_singletons='17291', 
            mz_calibration_ratio='-2.4537297844409523e-07', num_features_matched_csm='5910')
        '''
        dict_datasets = {}
        for c in self.chromatography_types:
            dict_datasets[c] = {}
        
        for f in self.good_features:
            iid = f['dataset_id']
            c = datasets_meta_dict[iid].chromatography
            if iid in dict_datasets[c]:
                dict_datasets[c][iid].append(f)
            else:
                dict_datasets[c][iid] = [f]
        
        return dict_datasets
        
        
    def summarize_chromatography_groups(self, concise=False):
        '''
        update dict_summary with number_ dataset stats and list of empCpds.
        Use good features only.
        '''
        for c in self.chromatography_types:
            group_dict, dict_csmf_members = self.analyze_group_features(c, self.dict_datasets[c])
            self.dict_csmf_members.update(dict_csmf_members)
            self.dict_summary[c] = self.make_concise_mr( group_dict, concise )
        
        
    def make_concise_mr(self, d, concise):
        '''
        Get a concise version of list of epds. Not in use now.
        Input is like:
        {'chromatography': 'HILIC',
                    'mode': 'pos',
                    'list_empCpds': [{'representative_feature': {'id': 'F1151',
                                    'rtime': 73.55,
                                    'is_good_peak': True,
                                    'snr': 3407.0,
                                    'dataset_id': 'ST001004_HILICpos_HILICpos_batch2_ppm5_361020',
                                    'parent_epd_id': 'kp11_85.0526',
                                    'neutral_formula_mass': 85.05264603323,
                                    'ion_relation': '13C/12C,M+H+'},
                                    'ion_relation': '13C/12C,M+H+',
                                    'ion_count': 41,
                                    'id': 'r1_pos_87.063389_epd_0'},
                                {'representative_feature': {'id': 'F1145',
                                    'rtime': 3.12,
                                    'is_good_peak': True,
                                    'snr': 25.0,
                                    'dataset_id': 'ST001004_HILICpos_HILICpos_batch2_ppm5_361020',
                                    'parent_epd_id': '_singleton_5461',
                                    'neutral_formula_mass': None,
                                    'ion_relation': ''},
                                    'ion_relation': '13C/12C,M+H+',
                                    'ion_count': 42,
                                    'id': 'r1_pos_87.063389_epd_1'}],
                    'representative': 'ST001004_HILICpos_HILICpos_batch2_ppm5_361020',
                    'median_feature_number': 2}
        '''
        def shrink_epd(epd):
            return {
                'id': epd['id'],
                'voted_ion_relation': epd['voted_ion_relation'],
                'representative_feature': {
                    'id': epd['representative_feature']['id'],
                    'rtime': epd['representative_feature']['rtime'],
                    'snr': epd['representative_feature']['snr'],
                }
            }
        if concise:
            return {
                    'representative': d['representative'],
                    'median_feature_number': d['median_feature_number'],
                    'number_datasets': d['number_datasets'],
                    'list_csm_features': [shrink_epd(e) for e in d['list_csm_features']]
                }
        else:
            return d
        
    def analyze_group_features(self, chromatography, datasets):
        '''
        datasets : dict of dataset_id to features, in a single method
        
            # list_good_features can be 100s to 1000s
            #for fl in datasets.values():
            #    new['list_good_features'] += [(x['mz'], x['rtime'], x['snr'], x['dataset_id']) for x in fl]
            # datasets[representative]
            
        dict_csmf_members : Plan to store members outside pos/neg_mass_registries.json, too bulky. 
        
        '''
        new = {'chromatography': chromatography,
               'mode': self.mode,
               'median_feature_number': None,
               'number_studies': None,
               'number_datasets': None,
               'representative': None,
               # 'list_good_features': [], 
               'list_csm_features': [],  # containing supporting features now
               }
        dict_csmf_members = {}
        
        if datasets:
            dataset_ids = list(datasets.keys())
            study_ids = set([x.split('_')[0] for x in dataset_ids])     # assuming id format 'MTBLS136_RPpos_B9_ppm5_359148'
            new['number_studies'] = len(study_ids)
            new['number_datasets'] = len(dataset_ids)
 
            median_feature_number = no_tie_int_median([len(v) for v in datasets.values()])
            new['median_feature_number'] = median_feature_number
            conformed_datasets_ids = [k for k,v in datasets.items() if len(v)==median_feature_number]
            
            # determine consensus features by popular votes; 
            list_csm_features, dict_csmf_members = self.vote_consensus_features(
                median_feature_number, [datasets[ii] for ii in conformed_datasets_ids], chromatography
            )
            new['list_csm_features'] = list_csm_features
            
            new['representative'] = self.get_representative_dataset(conformed_datasets_ids, datasets)
            new['representative_mass_track'] = []
            
        else:  
            print(self.id, datasets)
            
        return new, dict_csmf_members
        
    def get_representative_dataset(self, conformed_datasets_ids, datasets):
        '''
        Get representative dataset with best ave SNR. Only from datasets with consensus number of features.
        
        # Not showing representative_feature anymore; relevant info in representative_mass_track
        for feature in representative_dataset:
            list_csm_features.append(
                {'representative_feature' : feature}
            )
        '''
        ave_snrs = [
            np.mean([x['snr'] for x in datasets[ii]]) for ii in conformed_datasets_ids
        ]
        idx = np.argmax(ave_snrs)
        return conformed_datasets_ids[idx]

    def vote_consensus_features(self, median_feature_number, conformed_datasets, chromatography):
        '''
        return list_csm_features, dict_csmf_members 
        
        Including member raw features and inferred neutral mass.
        
        Ordered by RT now; 
        will add ranking by intensity orders and other methods later.
        Since the votes are limited to conformed mass tracks, not a major concern.
        
        '''
        def get_top_N_features(dataset, N):
            # return top N features in rtime order
            LL = sorted(dataset, reverse=True, key=lambda x: x['snr'])
            LL = sorted(LL[:N], key=lambda x: x['rtime'])
            return LL
            
        def get_voted_ion(list_features):
            # return top ion that is not '', e.g. (555, 'M0,M+H+'), if '' not the only one
            # input sensitive to rtime order
            r = self.summarize_ions(list_features)
            if r[0][1]:
                return r[0]
            elif len(r) > 1:
                return r[1]
            else:
                return r[0]
            
        list_csm_features = []
        dict_csmf_members = {}
        N = median_feature_number       #len(representative_dataset)
        # use up to N best (snr) features per dataset
        for ii in range(N):
            list_csm_features.append({})
            
        selected_features = []
        # conformed_datasets have same N
        for dataset in conformed_datasets:
            selected_features.append(get_top_N_features(dataset, N))
                
        for ii in range(N):
            _fid_ = self.id + '_' + chromatography + '_' + str(ii)
            members = [x[ii] for x in selected_features]
            count, ion = get_voted_ion(members)
            list_csm_features[ii]['id'] = _fid_
            list_csm_features[ii]['voted_ion_relation'] = ion
            list_csm_features[ii]['votes'] = count
            list_csm_features[ii]['neutral_mass'] = self.summarize_neutral_mass(members, voted_ion=ion)
            list_csm_features[ii]['num_members'] = len(conformed_datasets)
            
            # should save separately as members are many
            dict_csmf_members[_fid_] = [self.get_short_feature_tuple(x) for x in members]
            
        return list_csm_features, dict_csmf_members

    def summarize_ions(self, features):
        '''
        returns ions and counts in features, e.g.
            [(1126, ''),
            (555, 'M0,M+H+'),
            (61, 'M0,M+H+, 2x charged'),
            (49, 'M0,ACN'),
            (36, '13C/12C*2,M+H+'),
            (21, 'M0,HCl'),
            (14, 'M0,M+H+, 3x charged'),
            (10, 'M0,Na/H, 3x charged'),
            (10, 'M0,Na/H, 2x charged'),
            (7, '13C/12C,HCl'),
            (3, 'M0,K/H, 2x charged'),
            (3, 'M0,HCl, 2x charged'),
            (2, 'M0,K/H'),
            (2, '13C/12C,M+H+'),
            (2, '13C/12C, 2x charged,M+H+, 2x charged'),
            (2, '13C/12C*2,ACN'),
            (1, 'M0,K/H, 3x charged'),
            (1, 'M0,ACN, 3x charged'),
            (1, '13C/12C, 2x charged,HCl, 2x charged')]
        '''
        ion_relations = [f['ion_relation'] for f in features]
        _d = []
        for ion in set(ion_relations):
            _d.append( (ion_relations.count(ion), ion) )
        _d.sort(reverse=True)
        
        return _d

    def summarize_neutral_mass(self, features, voted_ion):
        '''
        Enforceing voted_ion since there can be discrepency btw median and voted ions.
        '''
        all = [x['neutral_formula_mass'] 
                          for x in features if x['neutral_formula_mass'] 
                          and x['ion_relation']==voted_ion
                          ]
        if all:
            return np.median(all)
        else:
            return None
        
    def get_short_feature_tuple(self, x):
        '''
        returns ('dataset_id', 'int_id', 'mz', 'rtime', 'snr')
        Use self.dict_datasets_int_id to convert int_id.
        '''
        return (
            self.dict_datasets_int_id[x['dataset_id']], 
            x['id'], 
            x['mz'], x['rtime'], x['snr']
        )

    def export_json(self):
        return self.dict_summary



class neutralMassRegistry:
    '''
    A unit for a specific neutral mass to include
    a) experimental observations by methods, e.g. [neg_RP, neg_HILIC, pos_RP, pos_HILIC, GCMS, MSMS]
    Each LC method may have multiple csm_features.
    b) organization of csm_features to csm_emp_cpds.
    c) list of theoretical compounds of this molecular weight (annotation by previous database or literature).
    d) user supplied annotation from metabolomics datasets.
    e) algorithm to assign annotations.
    
    '''
    def __init__(self, 
                 id, 
                 mass, 
                 primary_ions = ['M0,M+H+',
                                'M0,Na/H',
                                'M0,M+H+, 2x charged',
                                'M0,Na/H, 2x charged',
                                'M0,M+H+, 3x charged',
                                'M0,Na/H, 3x charged',
                                'M0,M-H-',
                                'M0,M-H-, 2x charged',
                                'M0,M-H-, 3x charged']
                 ):
        self.id = id
        self.neutral_mass = mass
        self.primary_ions = primary_ions
        self.data_by_method = {}            # assembled data here
        
        # self.csm_features = {}              # details of each csm_feature as dict
        self.list_theoretical_cpds = []     # from database records
        self.empcpd_annotation = {}            # dict for each empCpd
        
    def assemble_lcms_data(self, 
                      list_csm_features, 
                      dict_csm_fatures, 
                      dict_feature_relations,
                      expt_methods=['neg_RP', 'neg_HILIC', 'pos_RP', 'pos_HILIC'],
                      ):
        '''
        Updates self.data_by_method.
        
                'csm_features': [['r1_neg_415.099438_HILIC_1', 'M0,M-H-'],
                ['r1_neg_415.099438_RP_0', 'M0,M-H-'],
                ['r1_neg_415.099438_RP_1', 'M0,M-H-'],
                ['r1_neg_415.099438_RP_2', 'M0,M-H-'],
                ['r1_neg_416.102796_RP_3', '13C/12C,M-H-'],
                ['r1_neg_416.107064_HILIC_1', '13C/12C,M-H-']]
        '''
        data_by_method = self.organize_by_method(list_csm_features, expt_methods)
        for k, v in data_by_method.items():
            if v: # leave out empty lists
                self.data_by_method[k] = self.organize_csm_epds(
                                            v, dict_csm_fatures, dict_feature_relations)
            
    def organize_csm_epds(self, LL_csm_features, dict_csm_fatures, dict_feature_relations):
        '''
        Organize LL_csm_features in a single method into csm_empCpds.
        A csm_empCpd is similar to empCpd but based on csm_features (meta features summarized from many datasets).
        
        1. Get relationship pairs btw all csm_features.
        2. Determine which csm_features are primary, thus number of empCpds.
        3. Return empCpds and any unassigned csm_features (not primary by definition).
        
        method: 
        {
            epds: [{'primary': ['r1_neg_415.099438_RP_0', 'M0,M-H-'],
                    'others': [...],
                }, ...],
            unassigned: [...],
        }
        
        '''
        relationship_pairs  = []
        for pair in permutations(LL_csm_features, 2):
            if self.test_csm_features_belong(pair[0][0], pair[1][0], dict_csm_fatures, dict_feature_relations):
                relationship_pairs.append(pair)
        
        # can be used for diagnosis
        self.relationship_pairs = relationship_pairs
        # separate primary from secondary
        children = [x[1][0] for x in relationship_pairs]       # possible x -> y -> z, y is garanteed as child not paren
        primary = [x for x in LL_csm_features if x[1] in self.primary_ions and x[0] not in children]
                    
        epds, assigned = [], []
        for f in primary:
            matched_children = [pair[1][0] for pair in relationship_pairs if pair[0][0]==f[0]]
            matched_children_level2 = [pair[1][0] for pair in relationship_pairs if pair[0][0] in matched_children]
            mm = matched_children + matched_children_level2
            assigned += mm + [f[0]]
            epds.append(
                {
                    'primary': f,
                    'others': [x for x in LL_csm_features if x[0] in mm],
                }
            )
        unassigned = [x for x in LL_csm_features if x[0] not in assigned]
        return {
            'epds': epds,
            'unassigned': unassigned,
        }
        
    def summarize(self):
        '''
        Summarize info to export.
        GC-MS and MS/MS data are not handled here, but outside in current process.
        '''
        num_epds, num_epds_withkhipu = self.count_epd_khipus()
        
        return {
            'id': self.id,
            'mass': self.neutral_mass,
            'num_epds': num_epds,
            'num_epds_withkhipu': num_epds_withkhipu,
            'Methods': self.data_by_method,   
        }
        
    def count_epd_khipus(self):
        '''
        Report how many csm_epds and how many with khipus 
        '''
        num_epds, num_epds_withkhipu = 0, 0
        for k,v in self.data_by_method.items():
            for epd in v['epds']:
                num_epds += 1
                if len(epd['others']) > 0:  # valid khipu
                    num_epds_withkhipu += 1
                    
        return num_epds, num_epds_withkhipu
        
        
    def test_csm_features_belong(self, csmf, csmf2, dict_csm_fatures, dict_feature_relations):
        '''
        Test if csmf2 belong to csmf's compound.
        csmf1, csmf2 : 'r1_neg_415.099438_RP_0', 'r1_neg_416.102796_RP_3'
        If a primary feature in one csmf points to another, get a relationship.
        '''
        members = [(x[0], x[1]) for x in dict_csm_fatures[csmf]['members']]
        return set(members).intersection(self.get_parents_of_memberfeatures(
                dict_csm_fatures[csmf2]['members'], dict_feature_relations)
        )      
        
        
    def get_parents_of_memberfeatures(self, members, dict_feature_relations):
        '''
         'members': [[961, 'F46351', 415.1015, 431.23, 59.0],
            [850, 'F109715', 415.0989, 301.31, 30.0],
            [801, 'F56009', 415.0992, 318.63, 13.0]],
            
        dict_feature_relations : {(464, 'F535'): (464, 'F88'),
            (464, 'F536'): (464, 'F89'),
            (2, 'F693'): (2, 'F17')}
        '''
        r = [dict_feature_relations.get((m[0], m[1]), None) for m in members]
        return set([x for x in r if x])

        
    def organize_by_method(self, list_csm_features, expt_methods):
        '''
        Hardcoded in csm_features ID for now.
        expt_methods=['neg_RP', 'neg_HILIC', 'pos_RP', 'pos_HILIC']
        '''
        def assign(feature_id_str, expt_methods):
            segments = feature_id_str.split('_')
            matched = None
            for _m in expt_methods:
                method_keywords = _m.split('_')
                # checking two keywords only here
                if method_keywords[0] in segments and method_keywords[1] in segments:
                    matched = _m
            if not matched:
                print("No method match: ", feature_id_str)
            return matched
        
        d = {}
        for md in expt_methods:
            d[md] = []
        for f in list_csm_features:
            _m = assign(f[0], expt_methods)
            if _m:
                d[_m].append(f)
        return d
        
        
        
    def get_match_user_libraries(self, user_libraries):
        for LL in user_libraries:
            self.from_user_libraries.append(
                self.get_match_from_a_library(LL)
            )
    
    def get_match_from_a_library(self, library):
        
        pass
        
    def map_annotation(self):
        '''
        Mapping experimental data to annotation.
        Rationale:
        1. if only one empCpd from user expt and only one theoretical_cpd, treat as a match.
           If ionization is mismatch, lower confidence.
        2. if multiple empCpds and/or multiple theoretical_cpds, match rank order, by
           a) blood_conc or paper counts vs peak_area or snr,
           b) xlogP vs rtime
           c) to be added/improved methods
           
        returns recommendation per empCpd
        '''
        result = {}
        if len(self.list_theoretical_cpds) == 0:
            for x in self.list_empCpds:
                    result[x['interim_id']] = {
                    'best': None,
                    'others': []
                    }
        else:
            if len(self.list_empCpds) == 1:
                
                if len(self.list_theoretical_cpds) == 1:
                    result[self.list_empCpds[0]['interim_id']] = {
                        'best': self.list_theoretical_cpds[0],
                        'others': []
                    }
                else:
                    # pick best 
                    _best = self.top_cpd_by_blood_conc(self.list_theoretical_cpds)
                    if not _best:
                        _best = self.top_cpd_by_BloodPaperCount(self.list_theoretical_cpds)
                        if not _best:     # 1st one, as good as random 
                            _best = self.list_theoretical_cpds[0]
                    remaining = [x for x in self.list_theoretical_cpds if x['inchikey'] != _best['inchikey']]
                    result[self.list_empCpds[0]['interim_id']] = {
                        'best': _best,
                        'others': remaining
                    }
            else:       # Multiple input user empCpds
                if len(self.list_theoretical_cpds) == 1:
                    for x in self.list_empCpds:
                        result[x['interim_id']] = {
                        'best': self.list_theoretical_cpds[0],
                        'others': []
                        }
                else:
                    # will need develop this; should be method dependent
                    ranked_cpds = self.rank_cpds_(self.list_theoretical_cpds)
                    ranked_epds = self.rank_epds_by_intensity(self.list_empCpds)
                    NC, NE = len(ranked_cpds), len(ranked_epds)
                    if NC > NE:
                        # each epd at least gets one match
                        for ii, x in enumerate(ranked_epds):
                            result[x['interim_id']] = {
                                'best': ranked_cpds[ii],
                                'others': [x for x in self.list_theoretical_cpds if x['inchikey'] != ranked_cpds[ii]['inchikey']]
                            }
                    else:
                        for ii, y in enumerate(ranked_cpds):
                            result[ranked_epds[ii]['interim_id']] = {
                                'best': y,
                                'others': [x for x in self.list_theoretical_cpds if x['inchikey'] != y['inchikey']]
                            }
                        for ii in range(NE, NC):
                            result[ranked_epds[ii]['interim_id']] = {
                                'best': None,
                                'others': self.list_theoretical_cpds
                            }
                        
        self.annotation = result
    
    def map_annotation_singletons(self):
        # with list of singleton features
        # will rewrite
        #
        result = {}
        for x in self.list_empCpds:
            result[x['id']] = {
            'best': None,
            'others': self.list_theoretical_cpds
            }
        self.annotation = result
    
    
    def top_cpd_by_blood_conc(self, LL):
        new = [x for x in LL if x['blood_conc']]
        if new:
            return sorted(new, key=lambda x: x['blood_conc'])[-1]
        else:
            return None
        
    def top_cpd_by_BloodPaperCount(self, LL):
        new = [x for x in LL if 'BloodPaperCount' in x and x['BloodPaperCount']]
        if new:
            return sorted(new, key=lambda x: x['BloodPaperCount'])[-1]
        else:
            return None
        
    def rank_cpds_(self, LL):
        new = [x for x in LL if x['blood_conc']]
        others = [x for x in LL if not x['blood_conc']]
        new2, others2 = [], []
        for x in others: 
            if 'BloodPaperCount' in x and x['BloodPaperCount']:
                new2.append(x)
            else:
                others2.append(x)
        
        return sorted(new, reverse=True, key=lambda x: x['blood_conc']) + sorted(
            new2, reverse=True, key=lambda x: x['BloodPaperCount']
            ) + others2
        
    def rank_epds_by_intensity(self, LL, sort_key='peak_area'):
        # key can be ppeak_area, snr etc
        return sorted(LL, reverse=True, key=lambda x: sum([y[sort_key] for y in x['MS1_pseudo_Spectra']]))
        
        
        
    def export_json(self):
        # self.dict_summary
        return {
            'list_theoretical_cpds': self.list_theoretical_cpds,
            'list_CSM_cpds': self.list_CSM_cpds,
            'list_empCpds': self.list_empCpds,
            'annotation': self.annotation,
        }




#
# ----------------------------------------------
#

if __name__ == "__main__":
    
    import sys
    sys.path.append("/Users/lish/li.github/consensus_serum_metabolome/utils")
    from mining import read_master_datasets_records
    
    #with open('r1_dict_good_khipus.pickle', 'rb') as f:
    #    datasets_good_khipus = pickle.load(f)
        
    with open('r1_dict_cmt_matched_features_pos.pickle', 'rb') as f:
        dict_cmt = pickle.load(f)
        
    datasets_meta_dict = read_master_datasets_records(infile='r1_datasets_stats.tsv', sep='\t')
        
    test_cmt = 'r1_pos_302.232463'
    CMT = cmTrack(test_cmt,
               mode='pos',
               matched_features=dict_cmt[test_cmt],
               chromatography_types=['HILIC', 'RP']
               )
    
    print("\n\n~~~~~~~~~ start ~~~~~~~~~~\n\n")
    print(test_cmt)
    
    CMT.summarize_features( datasets_meta_dict )
    
    # print to test
    print("\n\n~~~~~~~~~ summary ~~~~~~~~~~\n\n")
    for k,v in CMT.dict_summary.items():
        print(k, v, '\n')
   