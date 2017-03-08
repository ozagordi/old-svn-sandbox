#PI wild type
PIWT = {
    10: 'L',
    11: 'V',
    16: 'G',
    20: 'K',
    24: 'L',
    30: 'D',
    32: 'V',
    33: 'L',
    34: 'E',
    36: 'M',
    46: 'M',
    47: 'I',
    48: 'G',
    50: 'I',
    53: 'F',
    54: 'I',
    60: 'D',
    62: 'I',
    63: 'L',
    64: 'I',
    71: 'A',
    73: 'G',
    76: 'L',
    77: 'V',
    82: 'V',
    84: 'I',
    85: 'I',
    88: 'N',
    89: 'L',
    90: 'L',
    93: 'I'
    }

# PIs drug resistant mutation
"""
PIRM = {
    'ataza': [10, 16, 20, 24, 32, 33, 34, 36, 46, 48,
              50, 53, 54, 60, 62, 64, 71, 73, 82, 84,
              85, 88, 90, 93],
    
    'fosam': [10, 32, 46, 47, 50, 54, 73, 76, 82, 84,
              90],

    'darun': [11, 32, 33, 47, 50, 54, 73, 76, 84, 89],

    'indin': [10, 20, 24, 32, 36, 46, 54, 71, 73, 76,
              77, 82, 84, 90],

    'lopin': [10, 20, 24, 32, 33, 46, 47, 50, 53, 54,
              63, 71, 73, 76, 82, 84, 90],

    'nelfi': [10, 30, 36, 46, 71, 77, 82, 84, 88, 90],

    'saqui': [10, 24, 48, 54, 62, 71, 73, 77, 82, 84,
              90],

    'tipra': [10, 13, 20, 33, 35, 36, 43, 46, 47, 54,
              58, 69, 74, 82, 83, 84, 90]
}
"""
PIRM = {
    'ataza': { 10: 'IFVC',
               16: 'E',
               20: 'RMITV',
               24: 'I',
               32: 'I',
               33: 'IFV',
               34: 'Q',
               36: 'ILV',
               46: 'IL',
               48: 'V',
               50: 'L',
               53: 'LY',
               54: 'LVMTA',
               60: 'E',
               62: 'V',
               64: 'LMV',
               71: 'VITL',
               73: 'CSTA',
               82: 'ATFI',
               84: 'V',
               85: 'V',
               88: 'S',
               90: 'M',
               93: 'LM' },
    
    'fosam': { 10: 'FIRV',
               32: 'I',
               46: 'IL',
               47: 'V',
               50: 'V',
               54: 'LVM',
               73: 'S',
               76: 'V',
               82: 'AFST',
               84: 'V',
               90: 'M' },

    'darun': { 11: 'I',
               32: 'I',
               33: 'F',
               47: 'V',
               50: 'V',
               54: 'ML',
               73: 'S',
               76: 'V',
               84: 'V',
               89: 'V' },

    'indin': { 10: 'IRV',
               20: 'MR',
               24: 'I',
               32: 'I',
               36: 'I',
               46: 'IL',
               54: 'V',
               71: 'VT',
               73: 'SA',
               76: 'V',
               77: 'I',
               82: 'AFT',
               84: 'V',
               90: 'M' },

    'lopin': { 10: 'FIRV',
               20: 'MR',
               24: 'I',
               32: 'I',
               33: 'F',
               46: 'IL',
               47: 'VA',
               50: 'V',
               53: 'L',
               54: 'VLAMTS',
               63: 'P',
               71: 'VT',
               73: 'S',
               76: 'V',
               82: 'AFTS',
               84: 'V',
               90: 'M' },

    'nelfi': { 10: 'FI',
               30: 'N',
               36: 'I',
               46: 'IL',
               71: 'VT',
               77: 'I',
               82: 'AFTS',
               84: 'V',
               88: 'DS',
               90: 'M' },

    'saqui': { 10: 'IRV',
               24: 'I',
               48: 'V',
               54: 'VL',
               62: 'V',
               71: 'VT',
               73: 'S',
               77: 'I',
               82: 'AFTS',
               84: 'V',
               90: 'M' },

    'tipra': { 10: 'V',
               13: 'V',
               20: 'MR',
               33: 'F',
               35: 'G',
               36: 'I',
               43: 'T',
               46: 'L',
               47: 'V',
               54: 'AMV',              
               58: 'E',
               69: 'K',
               74: 'P',
               82: 'LT',
               83: 'D',
               84: 'V',
               90: 'M' }
    }

drug_names = {
    'ataza': 'Atazanavir',
    'fosam': 'Fosamprenavir',
    'darun': 'Darunavir',
    'indin': 'Indinavir',
    'lopin': 'Lopinavir',
    'nelfi': 'Nelfinavir',
    'saqui': 'Saquinavir',
    'tipra': 'Tipranavir'
    }