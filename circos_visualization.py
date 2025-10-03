import numpy as np
import matplotlib.pyplot as plt
from prs_core import SNP_DATA, calculate_prs
import streamlit as st

# 染色体长度信息
CHROMOSOME_LENGTHS = {
    '1': 249250621, '2': 242193529, '3': 198295559, '4': 191154276,
    '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022,
    '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
    '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753,
    '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
    '21': 48129895, '22': 51304566
}

# 现代化配色方案
CHROMOSOME_COLORS = {
    '1': '#B5C3D7', '2': '#B6C9C0', '3': '#F5D8B7', '4': '#F1D0C6',
    '5': '#DDAEAB', '6': '#D0E0EF', '7': '#E9CDDF', '8': '#C8BFD9',
    '9': '#E4E3BF', '10': '#6E8FB2', '11': '#7DA494', '12': '#EAB67A',
    '13': '#E5A79A', '14': '#C16E71', '15': '#ABC8E5', '16': '#D8A0C1',
    '17': '#9F8DB8', '18': '#D0D63E', '19': '#B5C3D7', '20': '#B6C9C0',
    '21': '#F5D8B7', '22': '#F1D0C6'
}

# 主题色彩
THEME_COLORS = {
    'primary': '#6E8FB2',
    'secondary': '#7DA494', 
    'accent': '#EAB67A',
    'danger': '#C16E71',
    'warning': '#E5A79A',
    'success': '#7DA494',
    'info': '#ABC8E5',
    'light': '#F8F9FA',
    'dark': '#495057',
    'muted': '#6C757D'
}

def create_circos_plot(genotypes, selected_snp=None, figsize=(6, 6)):
    """
    创建优化的Circos图
    """
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection='polar'))
    
    # 计算总长度
    total_length = sum(CHROMOSOME_LENGTHS.values())
    
    # 从12点钟方向开始，顺时针排列
    current_angle = np.pi / 2
    
    # 绘制所有22条染色体，形成完整圆形
    for chrom in sorted(CHROMOSOME_LENGTHS.keys(), key=int):
        chrom_length = CHROMOSOME_LENGTHS[chrom]
        angle_span = (chrom_length / total_length) * 2 * np.pi
        end_angle = current_angle - angle_span
        
        # 绘制染色体弧
        theta = np.linspace(current_angle, end_angle, 30)
        r_inner, r_outer = 0.75, 0.9
        
        # 检查是否有SNP决定颜色
        has_snp = any(snp['chromosome'] == chrom for snp in SNP_DATA.values())
        if has_snp:
            color = CHROMOSOME_COLORS.get(chrom, '#CCCCCC')
            alpha = 0.8
        else:
            color = '#F0F0F0'
            alpha = 0.6
        
        ax.fill_between(theta, r_inner, r_outer,
                       color=color, alpha=alpha,
                       edgecolor='white', linewidth=0.8)
        
        # 为所有染色体添加标签
        mid_angle = (current_angle + end_angle) / 2
        label_r = 0.95
        ax.text(mid_angle, label_r, f"Chr{chrom}", 
               ha='center', va='center', fontsize=6, 
               weight='bold' if has_snp else 'normal',
               color='#333333' if has_snp else '#666666')
        
        current_angle = end_angle
    
    # 绘制SNP点（现在带有更丰富的注释信息）
    current_angle = np.pi / 2
    
    for chrom in sorted(CHROMOSOME_LENGTHS.keys(), key=int):
        chrom_length = CHROMOSOME_LENGTHS[chrom]
        angle_span = (chrom_length / total_length) * 2 * np.pi
        end_angle = current_angle - angle_span
        
        for rsid, snp_info in SNP_DATA.items():
            if snp_info['chromosome'] == chrom:
                position = snp_info['position']
                effect_weight = snp_info['effect_weight']
                relative_pos = position / CHROMOSOME_LENGTHS[chrom]
                snp_angle = current_angle - (relative_pos * angle_span)
                
                current_genotype = genotypes.get(rsid, 'Unknown')
                
                if effect_weight > 0:
                    color = THEME_COLORS['danger']
                    effect_type = 'Risk'
                else:
                    color = THEME_COLORS['info']
                    effect_type = 'Protective'
                
                size = max(40, min(150, abs(effect_weight) * 300))
                
                if rsid == selected_snp:
                    color = THEME_COLORS['accent']
                    size *= 1.3
                    edge_color = THEME_COLORS['primary']
                    edge_width = 3
                    alpha = 1.0
                    zorder = 100
                else:
                    edge_color = 'white'
                    edge_width = 1.5
                    zorder = 10
                    
                    effect_allele_count = current_genotype.count(snp_info['effect_allele'])
                    alpha = 0.5 + 0.25 * effect_allele_count
                
                # 绘制SNP点
                scatter = ax.scatter(snp_angle, 0.7, s=size, c=color, 
                          alpha=alpha, zorder=zorder,
                          edgecolors=edge_color, linewidths=edge_width)
                
                # 为选中的SNP添加详细标注
                if rsid == selected_snp:
                    annotation_text = (f"{rsid}\n"
                                     f"Gene: {snp_info.get('locus_name', 'Unknown')}\n"
                                     f"Genotype: {current_genotype}\n"
                                     f"Effect: {effect_weight:+.3f} ({effect_type})\n"
                                     f"Chr{chrom}:{position:,}")
                    
                    ax.annotate(annotation_text,
                               xy=(snp_angle, 0.7), xycoords='data',
                               xytext=(0.4, 0.4), textcoords='data',
                               fontsize=8, ha='center', va='center',
                               bbox=dict(boxstyle="round,pad=0.3", 
                                        facecolor='white', 
                                        edgecolor=THEME_COLORS['primary'],
                                        alpha=0.9),
                               arrowprops=dict(arrowstyle='->', 
                                             connectionstyle='arc3,rad=0.2',
                                             color=THEME_COLORS['primary']))
        
        current_angle = end_angle
    
    # 设置图形范围
    ax.set_ylim(0, 1.0)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)
    
    # 添加中心信息
    current_prs = calculate_prs(genotypes)
    effect_snps = sum(1 for rsid, genotype in genotypes.items() 
                     if SNP_DATA[rsid]['effect_allele'] in genotype)
    
    ax.text(0, 0, f"PRS\n{current_prs:.3f}", ha='center', va='center',
           fontsize=12, fontweight='bold', color=THEME_COLORS['primary'],
           bbox=dict(boxstyle="round,pad=0.25", facecolor='white', 
                    edgecolor=THEME_COLORS['primary'], linewidth=2))
    
    ax.text(0, -0.12, f"Effect SNPs: {effect_snps}/22", ha='center', va='center',
           fontsize=9, color=THEME_COLORS['muted'],
           bbox=dict(boxstyle="round,pad=0.15", facecolor='white', 
                    edgecolor=THEME_COLORS['muted']))
    
    # 添加图例
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=THEME_COLORS['danger'], 
                   markersize=7, alpha=0.8, label='Risk SNPs'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=THEME_COLORS['info'], 
                   markersize=7, alpha=0.8, label='Protective SNPs')
    ]
    if selected_snp:
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=THEME_COLORS['accent'], 
                      markersize=8, label='Selected SNP')
        )
    
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(-0.1, 1.0), fontsize=8)
    
    plt.title('AD PRS Genome Browser - Circos View', pad=10, fontsize=12, fontweight='bold')
    
    return fig

def display_circos_in_streamlit(genotypes, selected_snp=None):
    """
    在Streamlit中显示Circos图 - 使用容器控制大小
    """
    container = st.container()
    with container:
        col1, col_circos, col2 = st.columns([0.1, 1, 0.1])
        with col_circos:
            try:
                fig = create_circos_plot(genotypes, selected_snp, figsize=(6, 6))
                st.pyplot(fig, use_container_width=True)
                plt.close(fig)
            except Exception as e:
                st.error(f"Circos visualization error: {str(e)}")
                st.info("Please check your Python environment and dependencies.")

def check_pycircos_availability():
    """
    检查pyCircos是否可用 - 简化版本
    """
    return True