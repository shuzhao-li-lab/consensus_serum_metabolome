'''
pre-annotation algorithms implemented in class cmTrack.
How to get ion annotations from many source datasets by popular votes.
'''
import pickle
import numpy as np

# Primary ions used to select primary vTracks
primary_ions_pos = set(['M0,Na/H', 'M0,M+H+', 
                        'M0,M+H+, 2x charged', 'M0,Na/H, 2x charged', 
                        'M0,M+H+, 3x charged', 'M0,Na/H, 3x charged'])
primary_ions_neg = set(["M0,Na/H", "M0,M-H-", 
                        "M0,Na/H, 2x charged", "M0,M-H-, 2x charged", 
                        "M0,Na/H, 3x charged", "M0,M-H-, 3x charged"])



class extDataset:
    '''
    # Placeholder, 
    See read_master_datasets_records to get a dict of namedTuples for thsi purpose.
        
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
    
    
class cmTrack:
    '''
    A consensus mass track from large-scale data mining.
    Unknown number of compounds may exist on a cmTrack.
    The initial cmTrack is established by KDE peaks. 
    We keep list of good khipus here to facilitate preannotation.
    This implements popular votes to get represnetative empCpds on each chromatography and ionization mode.
    
    '''
    def __init__(self, 
                 id, mode, matched_features, 
                 chromatography_types=['HILIC', 'RP']
                 ):
        '''
        After a cmTrack is determined (via KDE), all matched features are input here and analyzed.
        matched_features are from tally_consensus_features_neg/pos.
        We also need an external reference of good khipus from all datasets.
        chromatography_types are 'HILIC' and 'RP' for now.
        '''
        self.id = id        # e.g., 'r1_pos_302.232463'; m/z is coded in id.
        self.mode = mode    # 'pos'   # or 'neg'
        self.matched_features = matched_features
        self.good_features = [f for f in matched_features if f['is_good_peak']]
        self.chromatography_types = chromatography_types
        self.is_primary_track = None
        self.num_all_features = len(self.matched_features)
        self.numb_good_features = len(self.good_features)
        self.dict_datasets = {} # separate data by methods
        self.dict_summary = {}
        
    def summarize_features(self, datasets_meta_dict):
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
        self.summarize_chromatography_groups()

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
        
        
    def summarize_chromatography_groups(self):
        '''
        update dict_summary with number_ dataset stats and list of empCpds.
        Use good features only.
        '''
        for c in self.chromatography_types:
            self.dict_summary[c] = self.make_concise_epds(
                                        self.analyze_group_features(c, self.dict_datasets[c])
                                        )
        
        
    def make_concise_epds(self, d):
        '''
        Get a concise version of list of epds. Input is like:
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
                'ion_count': epd['ion_count'],
                'ion_relation': epd['ion_relation'],
                'representative_feature': {
                    'id': epd['representative_feature']['id'],
                    'rtime': epd['representative_feature']['rtime'],
                    'snr': epd['representative_feature']['snr'],
                }
            }
        return {
                    'representative': d['representative'],
                    'median_feature_number': d['median_feature_number'],
                    'list_empCpds': [shrink_epd(e) for e in d['list_empCpds']]
                }
        
        
    def analyze_group_features(self, chromatography, datasets):
        '''
        datasets : dict of dataset_id to features, in a single method
        '''
        new = {'chromatography': chromatography,
               'mode': self.mode,
               'list_empCpds': [],
               'representative': None,
               'median_feature_number': None
               }
        if datasets:
            median_feature_number = np.median([len(v) for v in datasets.values()])
            median_feature_number = round(median_feature_number + 0.01) # covers the case of tied median, e.g. 2.5
            
            conformed_datasets_ids = [k for k,v in datasets.items() if len(v)==median_feature_number]
            if conformed_datasets_ids:
                representative = self.get_representative_dataset(conformed_datasets_ids, datasets)
            else:
                # representative is not matched to any; use the closest
                conformed_datasets_ids = sorted([(abs(len(v)-median_feature_number), k) for k,v in datasets.items()])
                representative = conformed_datasets_ids[0][1]
                conformed_datasets_ids = [representative]
                median_feature_number = len(datasets[representative])
                
            new['representative'] = representative
            new['median_feature_number'] = median_feature_number
            # determine empCpds by popular votes
            list_empCpds = self.vote_consensus_epds(
                datasets[representative], [datasets[ii] for ii in conformed_datasets_ids], chromatography
            )
            new['list_empCpds'] = list_empCpds
        else:  
            print(self.id, datasets)
            
        return new
        
        
    def get_representative_dataset(self, conformed_datasets_ids, datasets):
        '''
        Get representative dataset with best ave SNR.
        '''
        ave_snrs = [
            np.mean([x['snr'] for x in datasets[ii]]) for ii in conformed_datasets_ids
        ]
        idx = np.argmax(ave_snrs)
        return conformed_datasets_ids[idx]

    
    
    def vote_consensus_epds(self, representative_dataset, list_datasets, chromatography):
        '''
        return list_empCpds
        
        Ranked by RT now; 
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
            r = self.summarize_ions(list_features)
            if r[0][1]:
                return r[0]
            elif len(r) > 1:
                return r[1]
            else:
                return r[0]
            
        list_empCpds = []
        N = len(representative_dataset)
        # use up to N best (snr) features per dataset
        for feature in representative_dataset:
            list_empCpds.append(
                {'representative_feature' : feature}
            )
            
        selected_features = []
        # conformed_datasets have same N
        for dataset in list_datasets:
            selected_features.append(get_top_N_features(dataset, N))
                
        for ii in range(N):
            count, ion = get_voted_ion([x[ii] for x in selected_features])
            list_empCpds[ii]['ion_relation'] = ion
            list_empCpds[ii]['ion_count'] = count
            list_empCpds[ii]['id'] = self.id + '_' + chromatography + '_' + str(ii)

        return list_empCpds


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


    def determine_primary_track(self):
        pass

    def export_json(self):
        return self.dict_summary

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
   