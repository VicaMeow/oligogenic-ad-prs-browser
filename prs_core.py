import pandas as pd
import numpy as np

# 基于PGS000334的完整SNP数据 - 22个阿尔茨海默病相关SNP
# 包含从ad_snp_database_final.py提取的MAF数据
SNP_DATA = {
    'rs6656401': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': 0.14,
        'chromosome': '1',
        'position': 207692049,
        'locus_name': 'CR1',
        'ref_allele': 'A',
        'alt_allele': 'G', 
        'eur_freq_alt_allele': 0.8035
    },
    'rs744373': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': -0.14,
        'chromosome': '2',
        'position': 127894615,
        'locus_name': 'BIN1',
        'ref_allele': 'A',
        'alt_allele': 'G',
        'eur_freq_alt_allele': 0.228
    },
    'rs7419666': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.06,
        'chromosome': '2',
        'position': 234003359,
        'locus_name': 'INPP5D',
        'ref_allele': 'T',
        'alt_allele': 'C',
        'eur_freq_alt_allele': 0.4762
    },
    'rs9381563': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.08,
        'chromosome': '6',
        'position': 47432637,
        'locus_name': '',
        'ref_allele': 'C',
        'alt_allele': 'T',
        'eur_freq_alt_allele': 0.645
    },
    'rs1476679': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': 0.09,
        'chromosome': '7',
        'position': 100004446,
        'locus_name': 'ZCWPW1',
        'ref_allele': 'C',
        'alt_allele': 'T',
        'eur_freq_alt_allele': 0.7035
    },
    'rs7791765': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': 0.09,
        'chromosome': '7',
        'position': 143099107,
        'locus_name': 'EPHA1',
        'ref_allele': 'T',
        'alt_allele': 'G',
        'eur_freq_alt_allele': 0.201
    },
    'rs17057043': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': 0.08,
        'chromosome': '8',
        'position': 27220310,
        'locus_name': 'PTK2B',
        'ref_allele': 'G',
        'alt_allele': 'A',
        'eur_freq_alt_allele': 0.3363
    },
    'rs11136000': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.11,
        'chromosome': '8',
        'position': 27464519,
        'locus_name': 'CLU',
        'ref_allele': 'T',
        'alt_allele': 'C',
        'eur_freq_alt_allele': 0.6085
    },
    'rs7920721': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': -0.07,
        'chromosome': '10',
        'position': 11720308,
        'locus_name': 'AL512631.1',
        'ref_allele': 'A',
        'alt_allele': 'G',
        'eur_freq_alt_allele': 0.3895
    },
    'rs12292911': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': 0.06,
        'chromosome': '11',
        'position': 47449072,
        'locus_name': 'PSMC3',
        'ref_allele': 'G',
        'alt_allele': 'A',
        'eur_freq_alt_allele': 0.383
    },
    'rs7935829': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': 0.09,
        'chromosome': '11',
        'position': 59942815,
        'locus_name': 'MS4A6A',
        'ref_allele': 'A',
        'alt_allele': 'G',
        'eur_freq_alt_allele': 0.3836
    },
    'rs3851179': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.12,
        'chromosome': '11',
        'position': 85868640,
        'locus_name': 'RNU6-560P',
        'ref_allele': 'T',
        'alt_allele': 'C',
        'eur_freq_alt_allele': 0.6458
    },
    'rs11218343': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': 0.21,
        'chromosome': '11',
        'position': 121435587,
        'locus_name': 'SORL1',
        'ref_allele': 'T',
        'alt_allele': 'C',
        'eur_freq_alt_allele': 0.03893
    },
    'rs17125944': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.11,
        'chromosome': '14',
        'position': 53400629,
        'locus_name': 'FERMT2',
        'ref_allele': 'T',
        'alt_allele': 'C',
        'eur_freq_alt_allele': 0.09511
    },
    'rs941648': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': -0.07,
        'chromosome': '14',
        'position': 92931737,
        'locus_name': 'SLC24A4',
        'ref_allele': 'G',
        'alt_allele': 'A',
        'eur_freq_alt_allele': 0.773
    },
    'rs593742': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': 0.07,
        'chromosome': '15',
        'position': 59045774,
        'locus_name': 'ADAM10',
        'ref_allele': 'A',
        'alt_allele': 'G',
        'eur_freq_alt_allele': 0.2955
    },
    'rs7225151': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': 0.10,
        'chromosome': '17',
        'position': 5137047,
        'locus_name': 'SCIMP',
        'ref_allele': 'G',
        'alt_allele': 'A',
        'eur_freq_alt_allele': 0.1194
    },
    'rs9896864': {
        'effect_allele': 'A',
        'other_allele': 'G',
        'effect_weight': -0.21,
        'chromosome': '17',
        'position': 61536308,
        'locus_name': 'AC005828.5',
        'ref_allele': 'A',
        'alt_allele': 'G',
        'eur_freq_alt_allele': 0.02011
    },
    'rs3795065': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.08,
        'chromosome': '19',
        'position': 1039444,
        'locus_name': 'CNN2',
        'ref_allele': 'C',
        'alt_allele': 'T',
        'eur_freq_alt_allele': 0.6502
    },
    'rs7412': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.44,
        'chromosome': '19',
        'position': 45412079,
        'locus_name': 'APOE',
        'ref_allele': 'C',
        'alt_allele': 'T',
        'eur_freq_alt_allele': 0.07671
    },
    'rs429358': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -1.13,
        'chromosome': '19',
        'position': 45411941,
        'locus_name': 'APOE',
        'ref_allele': 'T',
        'alt_allele': 'C',
        'eur_freq_alt_allele': 0.1486
    },
    'rs6064392': {
        'effect_allele': 'T',
        'other_allele': 'C',
        'effect_weight': -0.11,
        'chromosome': '20',
        'position': 54984768,
        'locus_name': 'CASS4',
        'ref_allele': 'G',
        'alt_allele': 'T',
        'eur_freq_alt_allele': 0.08648
    }
}

def get_effect_allele_frequency(rsid, snp_info):
    """获取GWAS effect_allele在欧洲人群中的频率"""
    effect_allele = snp_info['effect_allele']
    alt_allele = snp_info['alt_allele']
    alt_freq = snp_info['eur_freq_alt_allele']
    
    if effect_allele == alt_allele:
        return alt_freq
    else:
        # effect_allele是ref_allele
        return 1 - alt_freq

def generate_realistic_genotype(rsid, snp_info):
    """根据MAF和Hardy-Weinberg平衡生成现实的基因型"""
    effect_freq = get_effect_allele_frequency(rsid, snp_info)
    other_freq = 1 - effect_freq
    
    effect_allele = snp_info['effect_allele']
    other_allele = snp_info['other_allele']
    
    # Hardy-Weinberg期望频率
    prob_effect_homo = effect_freq ** 2        # AA
    prob_hetero = 2 * effect_freq * other_freq  # Aa
    prob_other_homo = other_freq ** 2          # aa
    
    # 根据概率随机选择
    rand = np.random.random()
    if rand < prob_effect_homo:
        return effect_allele + effect_allele
    elif rand < prob_effect_homo + prob_hetero:
        return effect_allele + other_allele
    else:
        return other_allele + other_allele

def get_genotype_options(effect_allele, other_allele):
    """获取基因型选项"""
    options = [
        f"{other_allele}{other_allele}",  # 纯合子（其他等位基因）
        f"{effect_allele}{other_allele}",  # 杂合子
        f"{effect_allele}{effect_allele}"   # 纯合子（效应等位基因）
    ]
    return options

def calculate_genotype_score(genotype, effect_allele, other_allele, weight):
    """计算单个基因型的得分"""
    if genotype == f"{other_allele}{other_allele}":
        return 0  # 参考基因型
    elif genotype == f"{effect_allele}{other_allele}" or genotype == f"{other_allele}{effect_allele}":
        return weight  # 杂合子
    elif genotype == f"{effect_allele}{effect_allele}":
        return 2 * weight  # 纯合子
    else:
        return 0  # 默认值

def calculate_prs(genotypes):
    """计算多基因风险评分（PRS）"""
    total_score = 0
    
    for rsid, genotype in genotypes.items():
        if rsid in SNP_DATA:
            snp_info = SNP_DATA[rsid]
            score = calculate_genotype_score(
                genotype,
                snp_info['effect_allele'],
                snp_info['other_allele'],
                snp_info['effect_weight']
            )
            total_score += score
    
    return total_score

def initialize_default_genotypes():
    """初始化默认基因型 - 使用基于MAF的现实化随机生成"""
    genotypes = {}
    
    for rsid, snp_info in SNP_DATA.items():
        genotypes[rsid] = generate_realistic_genotype(rsid, snp_info)
    
    return genotypes

def generate_realistic_genotypes():
    """生成基于Hardy-Weinberg平衡的现实化基因型集合"""
    genotypes = {}
    
    for rsid, snp_info in SNP_DATA.items():
        genotypes[rsid] = generate_realistic_genotype(rsid, snp_info)
    
    return genotypes

def get_risk_interpretation(prs_score):
    """解释PRS分数的风险含义"""
    # 基于PGS000334的分数范围调整风险分层
    if prs_score > 0.5:
        return {
            'level': 'High Risk',
            'color': 'red',
            'description': 'Genetic risk significantly above average'
        }
    elif prs_score > 0.0:
        return {
            'level': 'Moderate Risk',
            'color': 'orange', 
            'description': 'Genetic risk slightly above average'
        }
    elif prs_score > -0.5:
        return {
            'level': 'Average Risk',
            'color': 'blue',
            'description': 'Genetic risk near population average'
        }
    else:
        return {
            'level': 'Low Risk',
            'color': 'green',
            'description': 'Genetic risk below average'
        }

def get_snp_summary_stats():
    """获取SNP汇总统计"""
    total_snps = len(SNP_DATA)
    positive_weights = sum(1 for snp in SNP_DATA.values() if snp['effect_weight'] > 0)
    negative_weights = sum(1 for snp in SNP_DATA.values() if snp['effect_weight'] < 0)
    
    return {
        'total_snps': total_snps,
        'risk_increasing': positive_weights,
        'protective': negative_weights
    }

def get_frequency_stats():
    """获取频率统计信息，用于调试和验证"""
    stats = {}
    
    for rsid, snp_info in SNP_DATA.items():
        effect_freq = get_effect_allele_frequency(rsid, snp_info)
        stats[rsid] = {
            'effect_allele': snp_info['effect_allele'],
            'effect_freq': effect_freq,
            'other_freq': 1 - effect_freq,
            'weight': snp_info['effect_weight']
        }
    
    return stats