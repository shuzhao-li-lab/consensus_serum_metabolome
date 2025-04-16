from mass2chem.lib.formula_coordinate import formula_coordinate
from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from mining import read_features_from_asari_table, custom_mz_calibrate

# ion patterns, refer to https://www.biorxiv.org/content/10.1101/2025.02.04.636472v2
isotope_search_patterns_pos = [ (1.003355, '13C/12C', (0, 0.8)),
                            (2.00671, '13C/12C*2', (0, 0.8)),
                            (3.010065, '13C/12C*3', (0, 0.8)),
                            # (3.9948, '44Ca/40Ca', (0, 0.1)), # 2%
                            (1.9970, '37Cl/35Cl', (0.1, 0.8)), # 24.24%
                            ]

isotope_search_patterns_neg = [ (1.003355, '13C/12C', (0, 0.8)),
                            (2.00671, '13C/12C*2', (0, 0.8)),
                            (3.010065, '13C/12C*3', (0, 0.8)),
                            (1.9970, '37Cl/35Cl', (0.1, 0.8)), # 24.24%
                            (1.9958, '32S/34S', (0, 0.1)), # 4%
                            ]

adduct_search_patterns_pos = [  # initial patterns are relative to M+H+
                            (21.98194, 'Na/H'),
                            (41.026549, 'ACN'),     # Acetonitrile
                            (67.987424, 'NaCOOH'),
                            (37.955882, 'K/H'),
                            (32.026215, 'CH3OH')
                            ]
adduct_search_patterns_neg = [  
                            (21.98194, 'Na/H'), 
                            (67.987424, 'NaCOOH'),
                            (82.0030, 'NaCH2COOH'),
                            # (1.99566, 'F <-> OH'), 
                            (41.026549, 'ACN'),
                            (37.955882, 'K/H'),
                            ]
extended_adducts = [  # can also exclude neutral loss here as a step after khipu
                            (1.0078, 'H'),
                            (17.02655, 'NH3'),
                            (18.0106, 'H2O'),      # easy to confuse with bio reactions
                            (18.033823, 'NH4'),
                            (27.01089904, 'HCN'),
                            (27.99492, 'CO'),
                            (32.026215, 'CH3OH'),
                            (35.9767, 'HCl'),
                            (37.94694, 'Ca/H2'),
                            (43.96389, 'Na2/H2'),
                            (46.00548, 'CO2H2'),
                            (67.987424, 'NaCOOH'),
                            (83.961361, 'KCOOH'),
                            (97.96737927, 'H2SO4'),
                            (97.97689507, 'H3PO4'),
                            # fragments
                            (-14.0155,  '\'14.01565\', "± CH2, alkane chains, waxes, fatty acids, methylation; or \'-'),
                            (-18.0104, '\'-18.010565\', \'H2O\', "{\'H\': -2, \'O\': -1}"'),
                            (-2.0156, '\'2.01565\', \'± 2H, opening or forming of double bond\', "{\'H\': 2}"'),
                            (-28.0312, '\'28.0313\', \'± C2H4, natural alkane chains such as fatty acids\', "{\'C\': 2, \'H\': 4}"'),
                            (-15.9948, '\'15.99492\', \'± O, e.g. oxidation/reduction\', "{\'O\': 1}"'),
                            (-17.0264, '\'-17.026549\', \'NH3\', "{\'N\': -1, \'H\': -3}"'),
                            (-26.0155, "' C2H2'"),
                            (-27.9948, '\'27.99492\', \'± CO\', "{\'C\': 1, \'O\': 1}"'),
                            (-32.026, '\'32.026215\', \'MeOH\', "{\'C\': 1, \'H\': 4, \'O\': 1}"'),
                            (-42.0104, '\'42.01057\', \'± COCH2\', "{\'C\': 2, \'O\': 1, \'H\': 2}"'),
                            (-67.9872, '\'67.987424\', \'NaCOOH\', "{\'C\': 1, \'O\': 2, \'Na\': 1, \'H\': 1}"'),
                            (-13.9791, '\'13.97927\', \'O <-> 2H, e.g. Oxidation follwed by H2O elimination\', "{\'H\': -2, \'O\': 1}"'),
                            (-42.0468, '\'42.04695\', \'± C3H6, propylation\', "{\'C\': 3, \'H\': 6}"'),
                            (-46.0053, '\'-46.005479\', \'H2O+CO\', "{\'C\': -1, \'H\': -2, \'O\': -2}"')
]


def run_khipu_insitu(asari_table, 
              KCD_formula_coordinate, 
              isotope_search_patterns,
              adduct_search_patterns,
              extended_adducts,
              mode = 'pos',
              mz_tolerance_ppm=5, rt_tolerance=2,
              out_json='epds.json',   # use None to switch off
              separate_singletons=False,
              ):
    '''
    returns list of khipus, [singletons,] features.
    Features are updated with parent_epd_id. 
    
    mz_tolerance_ppm : m/z tolerance in matching, user supplied parameter.
    rt_tolerance : retention time tolerance in seconds for matching, user supplied parameter.
    returns list of empCpds (singletons included) and list of features
    export json  if out_json
    '''
    _n, list_features_ = read_features_from_asari_table(open(asari_table).read())
    for f in list_features_:
        f['representative_intensity'] = f['peak_area']
  
    mass_accuracy_ratio, features = custom_mz_calibrate(list_features_, KCD_formula_coordinate, 
                                                    mode=mode, 
                                                    mz_tolerance_ppm=30, # this is for initial filter only
                                                    required_calibrate_threshold=0.000002)

    EED = ExperimentalEcpdDatabase(mode=mode, mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=rt_tolerance)
    EED.adduct_patterns = adduct_search_patterns
    EED.isotope_search_patterns = isotope_search_patterns
    EED.extended_adducts = extended_adducts

    # First khipu organized empCpds
    EED.build_from_list_peaks(features)
    singletons = [p for p in EED.dict_peaks if p not in EED.peak_to_empCpd.keys()]
    list_khipus = list(EED.dict_empCpds.values())
    # EED.peak_to_empCpd[P['id_number']] = interim_id
    for f in list_features_:
        f['parent_epd_id'] = EED.peak_to_empCpd.get(f['id'], '')
    
    EED = eed_combined_epdList(EED, singletons)
    if out_json:     # export empCpds here to JSON, including singletons.
        EED.export_empCpds(out_json)
    
    if separate_singletons:
        return list_khipus, singletons, list_features_
    else:
        return list(EED.dict_empCpds.values()), list_features_


def eed_combined_epdList(EED, singletons):
    new_id_start = len(EED.dict_empCpds)
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
    return EED

