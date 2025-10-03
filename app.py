import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math
from prs_core import (
    calculate_prs, 
    generate_realistic_genotypes,
    generate_realistic_genotype,
    get_genotype_options, 
    SNP_DATA, 
    get_risk_interpretation
)

# 导入新的circos可视化模块
from circos_visualization import display_circos_in_streamlit, check_pycircos_availability

# 现代化主题色彩
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

st.set_page_config(
    page_title="AD PRS Genome Browser",
    page_icon="🧠",
    layout="wide"
)

# 自定义CSS样式
st.markdown("""
<style>
    /* 主容器样式 */
    .main > div {
        padding-top: 1rem;
    }
    
    /* 卡片样式 */
    .custom-card {
        background: white;
        padding: 1.5rem;
        border-radius: 12px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        border: 1px solid #E9ECEF;
        margin-bottom: 1rem;
    }
    
    /* 标题样式 */
    .section-title {
        color: #495057;
        font-size: 1.1rem;
        font-weight: 600;
        margin-bottom: 1rem;
        display: flex;
        align-items: center;
        gap: 0.5rem;
    }
    
    /* 按钮组样式 */
    .button-group {
        display: flex;
        gap: 0.5rem;
        flex-wrap: wrap;
    }
    
    /* 统计卡片 */
    .stat-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
        padding: 1rem;
        border-radius: 8px;
        text-align: center;
        border: 1px solid #dee2e6;
    }
    
    .stat-value {
        font-size: 1.5rem;
        font-weight: bold;
        color: #495057;
    }
    
    .stat-label {
        font-size: 0.9rem;
        color: #6c757d;
        margin-top: 0.25rem;
    }
    
    /* Streamlit组件优化 */
    .stSelectbox > div > div > select {
        border-radius: 8px;
        border: 2px solid #E9ECEF;
    }
    
    .stButton > button {
        border-radius: 8px;
        font-weight: 500;
        transition: all 0.2s ease;
    }
    
    .stButton > button:hover {
        transform: translateY(-1px);
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }
    
    /* 隐藏Streamlit默认元素 */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

def load_disclaimer():
    """加载disclaimer内容"""
    try:
        with open('disclaimer.md', 'r', encoding='utf-8') as file:
            return file.read()
    except FileNotFoundError:
        return """
# Educational Disclaimer

**⚠️ IMPORTANT: FOR EDUCATIONAL AND RESEARCH PURPOSES ONLY**

This application is developed as a coursework project and is intended solely for educational demonstration and learning purposes.

## Medical Disclaimer
• NOT medical advice, diagnosis, or treatment recommendations
• NOT predictive of actual Alzheimer's disease risk or health outcomes
• Consult qualified healthcare providers for any health concerns

## Technical Limitations
• Uses simulated genotypes for educational demonstration only
• Population percentiles are simulated and do NOT represent actual population distributions
• Software provided "AS IS" without warranty of any kind

**By using this application, you acknowledge this is an educational tool and agree not to use it for medical decision-making.**
        """

def show_disclaimer_page():
    """显示disclaimer页面"""
    st.markdown("""
    <div style="text-align: center; padding: 2rem 0;">
        <h1 style="color: #495057;">🧠 AD PRS Genome Browser</h1>
        <p style="color: #6C757D; font-size: 1.1rem;">Educational Polygenic Risk Score Tool</p>
    </div>
    """, unsafe_allow_html=True)
    
    # 显示disclaimer内容
    disclaimer_content = load_disclaimer()
    st.markdown(disclaimer_content)
    
    # 同意按钮
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 1, 1])
    with col2:
        if st.button("✅ I have read and agree to the terms", 
                    use_container_width=True, 
                    type="primary"):
            st.session_state.disclaimer_accepted = True
            st.rerun()

def render_snp_dropdown():
    """渲染SNP下拉选择器"""
    st.markdown('<div class="section-title">🎯 SNP Selection</div>', unsafe_allow_html=True)
    
    # 创建SNP选项列表
    snp_options = ["Select SNP..."]
    for chrom in sorted(set(info['chromosome'] for info in SNP_DATA.values()), key=int):
        for rsid, info in SNP_DATA.items():
            if info['chromosome'] == chrom:
                effect_indicator = "🔴" if info['effect_weight'] > 0 else "🔵"
                snp_options.append(f"{effect_indicator} Chr{chrom} - {rsid}")
    
    # 下拉选择器
    selected_option = st.selectbox(
        "Select SNP to edit:",
        snp_options,
        key="snp_dropdown",
        help="Red=Risk type, Blue=Protective type"
    )
    
    # 处理选择
    if selected_option != "Select SNP...":
        # 从选项中提取rsid
        rsid = selected_option.split(" - ")[-1]
        if rsid != st.session_state.get('selected_snp'):
            st.session_state.selected_snp = rsid
            st.rerun()

def render_compact_editor():
    """渲染紧凑的编辑器"""
    selected_snp = st.session_state.get('selected_snp', None)
    
    if selected_snp and selected_snp in SNP_DATA:
        snp_info = SNP_DATA[selected_snp]
        current_genotype = st.session_state.genotypes[selected_snp]
        
        # 简化的配色方案
        CHROMOSOME_COLORS = {
            '1': '#B5C3D7', '2': '#B6C9C0', '3': '#F5D8B7', '4': '#F1D0C6',
            '5': '#DDAEAB', '6': '#D0E0EF', '7': '#E9CDDF', '8': '#C8BFD9',
            '9': '#E4E3BF', '10': '#6E8FB2', '11': '#7DA494', '12': '#EAB67A',
            '13': '#E5A79A', '14': '#C16E71', '15': '#ABC8E5', '16': '#D8A0C1',
            '17': '#9F8DB8', '18': '#D0D63E', '19': '#B5C3D7', '20': '#B6C9C0',
            '21': '#F5D8B7', '22': '#F1D0C6'
        }
        
        chrom_color = CHROMOSOME_COLORS.get(snp_info['chromosome'], '#CCCCCC')
        
        # 编辑器标题
        st.markdown(f"""
        <div style="background: {chrom_color}; 
                    color: white; padding: 12px; border-radius: 8px; margin-bottom: 15px; 
                    text-align: center; box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
            <h4 style="margin: 0; font-size: 1.1rem;">✏️ {selected_snp}</h4>
            <small style="opacity: 0.9;">Chr{snp_info['chromosome']} | {snp_info.get('locus_name', 'Unknown')}</small>
        </div>
        """, unsafe_allow_html=True)
        
        # 基因型编辑
        st.markdown("**Genotype Selection**")
        options = get_genotype_options(snp_info['effect_allele'], snp_info['other_allele'])
        
        new_genotype = st.selectbox(
            "Select genotype:",
            options,
            index=options.index(current_genotype) if current_genotype in options else 0,
            key=f"edit_{selected_snp}",
            help=f"Effect allele: {snp_info['effect_allele']}"
        )
        
        # 效应预测
        st.markdown("**Effect Prediction**")
        for genotype in options:
            effect_count = genotype.count(snp_info['effect_allele'])
            effect_score = effect_count * snp_info['effect_weight']
            is_current = (genotype == new_genotype)
            
            # 状态指示器
            if effect_count == 2:
                status = "🔴 High"
                status_color = THEME_COLORS['danger']
            elif effect_count == 1:
                status = "🟡 Med"
                status_color = THEME_COLORS['warning']
            else:
                status = "⚪ None"
                status_color = THEME_COLORS['muted']
            
            # 背景色
            if is_current:
                bg_color = "rgba(234, 182, 122, 0.2)"
                border = f"2px solid {THEME_COLORS['accent']}"
            else:
                bg_color = THEME_COLORS['light']
                border = "1px solid #DEE2E6"
            
            st.markdown(f"""
            <div style="background: {bg_color}; padding: 8px 12px; border-radius: 6px; 
                        margin-bottom: 6px; border: {border}; font-size: 14px;">
                <strong>{genotype}</strong> {status} → <strong style="color: {status_color};">{effect_score:+.3f}</strong>
                {'👈' if is_current else ''}
            </div>
            """, unsafe_allow_html=True)
        
        # 应用更改按钮
        if new_genotype != current_genotype:
            if st.button("✅ Apply Changes", use_container_width=True, type="primary"):
                st.session_state.genotypes[selected_snp] = new_genotype
                st.session_state.selected_snp = None
                st.success(f"✅ {selected_snp} Updated")
                st.rerun()
        
        # 快捷操作
        col1, col2 = st.columns(2)
        with col1:
            if st.button("🎲", use_container_width=True, help="MAF-based random genotype"):
                st.session_state.genotypes[selected_snp] = generate_realistic_genotype(selected_snp, snp_info)
                st.session_state.selected_snp = None
                st.rerun()
        with col2:
            if st.button("✖", use_container_width=True, help="Close editor"):
                st.session_state.selected_snp = None
                st.rerun()
    
    else:
        # 未选中状态
        st.markdown(f"""
        <div class="custom-card" style="text-align: center; color: {THEME_COLORS['muted']};">
            <div class="section-title" style="justify-content: center;">✏️ Genotype Editor</div>
            <p style="font-size: 14px; margin: 0;">Select SNP from the left dropdown menu</p>
        </div>
        """, unsafe_allow_html=True)

def create_percentile_chart():
    """创建基于真实欧洲人群PRS分布的百分位图表"""
    current_prs = calculate_prs(st.session_state.genotypes)
    
    # 真实的欧洲人群PRS分布参数
    POPULATION_MEAN = -2.4101
    POPULATION_STD = 0.6444
    
    # 理论极值
    THEORETICAL_MIN = -5.46
    THEORETICAL_MAX = 1.86
    
    # 计算当前PRS在人群中的百分位数
    from scipy.stats import norm
    percentile = norm.cdf(current_prs, loc=POPULATION_MEAN, scale=POPULATION_STD) * 100
    # 限制百分位数范围在合理区间内
    percentile = max(0.1, min(99.9, percentile))
    
    st.markdown('<div class="section-title">📊 Population Percentile</div>', unsafe_allow_html=True)
    
    # 创建正态分布曲线图
    fig = go.Figure()
    
    # 生成基于真实参数的正态分布数据，使用完整理论范围
    x_min, x_max = THEORETICAL_MIN, THEORETICAL_MAX
    x_range = np.linspace(x_min, x_max, 500)
    
    # 计算正态分布概率密度
    y_normal = norm.pdf(x_range, loc=POPULATION_MEAN, scale=POPULATION_STD)
    
    # 绘制正态分布曲线
    fig.add_trace(go.Scatter(
        x=x_range,
        y=y_normal,
        mode='lines',
        line=dict(color=THEME_COLORS['info'], width=2),
        fill='tozeroy',
        fillcolor='rgba(171, 200, 229, 0.3)',
        name='European Population Distribution',
        showlegend=False,
        hovertemplate='PRS: %{x:.3f}<br>Density: %{y:.4f}<extra></extra>'
    ))
    
    # 如果当前PRS在显示范围内，添加标记
    if x_min <= current_prs <= x_max:
        # 计算当前PRS对应的概率密度
        user_y = norm.pdf(current_prs, loc=POPULATION_MEAN, scale=POPULATION_STD)
        
        # 添加当前位置的垂直线
        fig.add_trace(go.Scatter(
            x=[current_prs, current_prs],
            y=[0, user_y],
            mode='lines',
            line=dict(color=THEME_COLORS['primary'], width=3),
            showlegend=False,
            hoverinfo='skip'
        ))
        
        # 添加当前位置点
        fig.add_trace(go.Scatter(
            x=[current_prs],
            y=[user_y],
            mode='markers',
            marker=dict(
                size=10,
                color=THEME_COLORS['primary'],
                symbol='circle',
                line=dict(color='white', width=2)
            ),
            showlegend=False,
            hovertemplate=f'Your PRS: {current_prs:.3f}<br>Percentile: {percentile:.1f}%<extra></extra>'
        ))
        
        # 添加百分位数标注
        fig.add_annotation(
            x=current_prs,
            y=user_y,
            text=f"{percentile:.1f}%",
            showarrow=True,
            arrowhead=2,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor=THEME_COLORS['primary'],
            ax=0,
            ay=-25,
            font=dict(size=11, color=THEME_COLORS['primary'], family="Arial Bold"),
            bgcolor="white",
            bordercolor=THEME_COLORS['primary'],
            borderwidth=1
        )
    
    # 添加理论极值标记
    # 理论最小值标记
    fig.add_trace(go.Scatter(
        x=[THEORETICAL_MIN],
        y=[0],
        mode='markers',
        marker=dict(
            size=8,
            color=THEME_COLORS['success'],
            symbol='triangle-up',
            line=dict(color='white', width=1)
        ),
        showlegend=False,
        hovertemplate=f'Theoretical Minimum: {THEORETICAL_MIN}<extra></extra>'
    ))
    
    # 理论最大值标记
    fig.add_trace(go.Scatter(
        x=[THEORETICAL_MAX],
        y=[0],
        mode='markers',
        marker=dict(
            size=8,
            color=THEME_COLORS['danger'],
            symbol='triangle-up',
            line=dict(color='white', width=1)
        ),
        showlegend=False,
        hovertemplate=f'Theoretical Maximum: {THEORETICAL_MAX}<extra></extra>'
    ))
    
    # 添加极值文本标注
    fig.add_annotation(
        x=THEORETICAL_MIN,
        y=max(y_normal) * 0.1,
        text="Theoretical<br>Minimum",
        showarrow=False,
        font=dict(size=8, color=THEME_COLORS['success']),
        xanchor='center',
        bgcolor='rgba(255,255,255,0.8)',
        bordercolor=THEME_COLORS['success'],
        borderwidth=1,
        borderpad=2
    )
    
    fig.add_annotation(
        x=THEORETICAL_MAX,
        y=max(y_normal) * 0.1,
        text="Theoretical<br>Maximum",
        showarrow=False,
        font=dict(size=8, color=THEME_COLORS['danger']),
        xanchor='center',
        bgcolor='rgba(255,255,255,0.8)',
        bordercolor=THEME_COLORS['danger'],
        borderwidth=1,
        borderpad=2
    )
    
    # 添加分布统计信息标注
    fig.add_annotation(
        x=x_min + 0.3,
        y=max(y_normal) * 0.8,
        text=f"European Population<br>Mean: {POPULATION_MEAN:.3f}<br>SD: {POPULATION_STD:.3f}",
        showarrow=False,
        font=dict(size=9, color=THEME_COLORS['muted']),
        xanchor='left',
        bgcolor='rgba(255,255,255,0.9)',
        bordercolor=THEME_COLORS['muted'],
        borderwidth=1,
        borderpad=4
    )
    
    fig.update_layout(
        xaxis=dict(
            title="PRS Score",
            showgrid=True,
            gridcolor='#F0F0F0',
            range=[x_min, x_max],
            dtick=1.0  # 每1.0显示一个刻度
        ),
        yaxis=dict(
            title="Probability Density",
            showgrid=False,
            showticklabels=True,
            tickformat='.3f'
        ),
        height=250,
        margin=dict(t=15, b=40, l=45, r=15),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    st.plotly_chart(fig, use_container_width=True, key="percentile_chart")
    
    # 添加百分位数解释
    if percentile < 5:
        percentile_desc = "Very Low (Bottom 5%)"
        percentile_color = THEME_COLORS['success']
    elif percentile < 25:
        percentile_desc = "Below Average (Bottom 25%)"
        percentile_color = THEME_COLORS['info']
    elif percentile < 75:
        percentile_desc = "Average Range (25-75%)"
        percentile_color = THEME_COLORS['muted']
    elif percentile < 95:
        percentile_desc = "Above Average (Top 25%)"
        percentile_color = THEME_COLORS['warning']
    else:
        percentile_desc = "Very High (Top 5%)"
        percentile_color = THEME_COLORS['danger']
    
    # 显示解释信息
    st.markdown(f"""
    <div style="text-align: center; padding: 8px; background: rgba(171, 200, 229, 0.1); 
                border-radius: 6px; border-left: 3px solid {percentile_color};">
        <small style="color: {THEME_COLORS['dark']};">
            <strong>Your Position:</strong><br>
            {percentile:.1f}th percentile<br>
            <span style="color: {percentile_color};">{percentile_desc}</span>
        </small>
    </div>
    """, unsafe_allow_html=True)

def render_control_panel():
    """渲染控制面板"""
    st.markdown('<div class="section-title">🎛️ Global Controls</div>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("🎲 Randomize All", use_container_width=True, help="MAF-based realistic randomization"):
            st.session_state.genotypes = generate_realistic_genotypes()
            st.session_state.selected_snp = None
            st.rerun()
    
    with col2:
        if st.button("🛡️ Low Risk", use_container_width=True, help="Maximum protective configuration"):
            for rsid, snp_info in SNP_DATA.items():
                effect_allele = snp_info['effect_allele']
                other_allele = snp_info['other_allele']
                if snp_info['effect_weight'] < 0:
                    st.session_state.genotypes[rsid] = effect_allele + effect_allele
                else:
                    st.session_state.genotypes[rsid] = other_allele + other_allele
            st.session_state.selected_snp = None
            st.rerun()

    with col3:
        if st.button("⚡ High Risk", use_container_width=True, help="Maximum risk configuration"):
            for rsid, snp_info in SNP_DATA.items():
                effect_allele = snp_info['effect_allele']
                other_allele = snp_info['other_allele']
                if snp_info['effect_weight'] > 0:
                    st.session_state.genotypes[rsid] = effect_allele + effect_allele
                else:
                    st.session_state.genotypes[rsid] = other_allele + other_allele
            st.session_state.selected_snp = None
            st.rerun()

def render_summary_stats():
    """渲染摘要统计信息"""
    current_prs = calculate_prs(st.session_state.genotypes)
    
    # 计算统计信息
    effect_snps = sum(1 for rsid, genotype in st.session_state.genotypes.items() 
                     if SNP_DATA[rsid]['effect_allele'] in genotype)
    
    protective_snps = sum(1 for rsid, genotype in st.session_state.genotypes.items() 
                         if SNP_DATA[rsid]['effect_weight'] < 0 and 
                         SNP_DATA[rsid]['effect_allele'] in genotype)
    
    risk_snps = sum(1 for rsid, genotype in st.session_state.genotypes.items() 
                   if SNP_DATA[rsid]['effect_weight'] > 0 and 
                   SNP_DATA[rsid]['effect_allele'] in genotype)
    
    # 计算百分位数
    from scipy.stats import norm
    POPULATION_MEAN = -2.4101
    POPULATION_STD = 0.6444
    percentile = norm.cdf(current_prs, loc=POPULATION_MEAN, scale=POPULATION_STD) * 100
    percentile = max(0.1, min(99.9, percentile))
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-value" style="color: {THEME_COLORS['primary']}">{current_prs:.3f}</div>
            <div class="stat-label">PRS Score</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-value" style="color: {THEME_COLORS['accent']}">{percentile:.1f}%</div>
            <div class="stat-label">Percentile</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-value" style="color: {THEME_COLORS['danger']}">{risk_snps}</div>
            <div class="stat-label">Risk SNPs</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div class="stat-card">
            <div class="stat-value" style="color: {THEME_COLORS['success']}">{protective_snps}</div>
            <div class="stat-label">Protective SNPs</div>
        </div>
        """, unsafe_allow_html=True)

def show_app_content():
    """显示主应用内容"""
    # 顶部简化警告
    st.markdown("""
    <div style="background-color: #FFF3CD; padding: 0.5rem 1rem; border-radius: 0.25rem; 
                border-left: 3px solid #FFC107; margin-bottom: 1rem;">
        <small style="color: #856404;">
            ⚠️ Educational Tool Only - Not for Medical Use
        </small>
    </div>
    """, unsafe_allow_html=True)
    
    # 检查pyCircos可用性并显示状态
    if not check_pycircos_availability():
        st.warning("""
        📦 **pyCircos not detected** - Using fallback visualization mode.  
        For enhanced Circos plots, install: `pip install python-circos`
        """)
    
    # 原始应用内容
    render_summary_stats()
    
    st.markdown("---")
    
    # 三列布局：SNP选择 | Circos图 | 编辑器+百分位
    col_select, col_circos, col_right = st.columns([1, 2.5, 1])
    
    with col_select:
        render_snp_dropdown()
        st.markdown("---")
        create_percentile_chart()
    
    with col_circos:
        # 使用新的circos可视化模块
        st.markdown('<div class="section-title" style="text-align: center;">🧠 Genome Browser - Circos View</div>', unsafe_allow_html=True)
        display_circos_in_streamlit(
            st.session_state.genotypes, 
            st.session_state.get('selected_snp', None)
        )
    
    with col_right:
        render_compact_editor()
    
    st.markdown("---")
    
    # 控制面板
    render_control_panel()
    
    # 页脚信息
    st.markdown(f"""
    <div style="text-align: center; color: {THEME_COLORS['muted']}; font-size: 12px; margin-top: 2rem; padding: 1rem;">
        <p>🧬 Alzheimer's Disease Polygenic Risk Score based on PGS000334 | 22 Associated SNPs</p>
        <p style="font-size: 10px;">⚠️ For research and educational purposes only, not medical advice</p>
        <p style="font-size: 10px;">Educational tool developed for PSYC301 coursework</p>
    </div>
    """, unsafe_allow_html=True)

def main():
    # 初始化disclaimer状态
    if 'disclaimer_accepted' not in st.session_state:
        st.session_state.disclaimer_accepted = False
    
    # 根据disclaimer状态显示不同页面
    if not st.session_state.disclaimer_accepted:
        show_disclaimer_page()
    else:
        # 初始化会话状态
        if 'genotypes' not in st.session_state:
            st.session_state.genotypes = generate_realistic_genotypes()
        if 'selected_snp' not in st.session_state:
            st.session_state.selected_snp = None
        
        show_app_content()

if __name__ == "__main__":
    main()