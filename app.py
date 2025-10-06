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

from circos_visualization import display_circos_in_streamlit, check_pycircos_availability

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
    'muted': '#6C757D',
    'purple': '#9F8DB8',
    'yellow': '#D0D63E'
}

st.set_page_config(
    page_title="AD PRS Genome Browser",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
<style>
    .main {
        background: linear-gradient(135deg, #E9ECEF 0%, #D4DBE8 50%, #C8D5E3 100%);
        background-attachment: fixed;
    }
    
    .main > div {
        padding-top: 0 !important;
        padding-bottom: 0;
    }
    
    .block-container {
        padding-top: 1rem !important;
        padding-bottom: 1rem;
    }
    
    .glass-card {
        background: rgba(255, 255, 255, 0.75);
        backdrop-filter: blur(12px);
        -webkit-backdrop-filter: blur(12px);
        border-radius: 16px;
        border: 1px solid rgba(255, 255, 255, 0.4);
        box-shadow: 0 8px 32px rgba(110, 143, 178, 0.15);
        padding: 1.5rem;
        margin-bottom: 1rem;
        transition: all 0.3s ease;
    }
    
    .glass-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 12px 40px rgba(110, 143, 178, 0.2);
    }
    
    .section-header {
        color: #495057;
        font-size: 0.95rem;
        font-weight: 700;
        margin-bottom: 1rem;
        text-transform: uppercase;
        letter-spacing: 1px;
        padding-bottom: 8px;
        border-bottom: 2px solid rgba(110, 143, 178, 0.3);
    }
    
    .stat-card-glass {
        background: linear-gradient(135deg, 
            rgba(255, 255, 255, 0.8) 0%, 
            rgba(248, 249, 250, 0.7) 100%);
        backdrop-filter: blur(10px);
        -webkit-backdrop-filter: blur(10px);
        padding: 0.8rem 1rem;
        border-radius: 12px;
        text-align: center;
        border: 1px solid rgba(255, 255, 255, 0.5);
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        transition: all 0.3s ease;
    }
    
    .stat-card-glass:hover {
        transform: translateY(-2px) scale(1.01);
        box-shadow: 0 6px 24px rgba(110, 143, 178, 0.2);
    }
    
    .stat-value {
        font-size: 1.5rem;
        font-weight: bold;
        margin-bottom: 0.2rem;
    }
    
    .stat-label {
        font-size: 0.75rem;
        color: #6c757d;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        font-weight: 600;
    }
    
    .snp-header-glass {
        backdrop-filter: blur(8px);
        -webkit-backdrop-filter: blur(8px);
        color: white;
        padding: 14px;
        border-radius: 10px;
        margin-bottom: 1rem;
        text-align: center;
        box-shadow: 0 4px 16px rgba(0, 0, 0, 0.15);
        border: 1px solid rgba(255, 255, 255, 0.3);
    }
    
    .genotype-option {
        background: rgba(255, 255, 255, 0.6);
        backdrop-filter: blur(8px);
        -webkit-backdrop-filter: blur(8px);
        padding: 10px 14px;
        border-radius: 8px;
        margin-bottom: 8px;
        font-size: 14px;
        border: 1px solid rgba(222, 226, 230, 0.6);
        transition: all 0.2s ease;
    }
    
    .genotype-option:hover {
        background: rgba(255, 255, 255, 0.8);
        transform: translateX(3px);
    }
    
    .warning-glass {
        background: rgba(255, 243, 205, 0.85);
        backdrop-filter: blur(10px);
        -webkit-backdrop-filter: blur(10px);
        padding: 0.8rem 1.2rem;
        border-radius: 8px;
        border-left: 4px solid #EAB67A;
        margin-bottom: 1.2rem;
        box-shadow: 0 2px 12px rgba(234, 182, 122, 0.2);
    }
    
    .snp-badge {
        display: inline-block;
        padding: 3px 10px;
        border-radius: 12px;
        font-size: 10px;
        font-weight: 700;
        letter-spacing: 0.5px;
        margin-left: 6px;
    }
    
    .badge-risk {
        background: rgba(193, 110, 113, 0.2);
        color: #C16E71;
        border: 1px solid #C16E71;
    }
    
    .badge-protective {
        background: rgba(125, 164, 148, 0.2);
        color: #7DA494;
        border: 1px solid #7DA494;
    }
    
    .circos-container {
        height: 1px;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    .stSelectbox > div > div {
        background: rgba(255, 255, 255, 0.7);
        backdrop-filter: blur(8px);
        -webkit-backdrop-filter: blur(8px);
        border-radius: 8px;
        border: 2px solid rgba(110, 143, 178, 0.3);
    }
    
    .stButton > button {
        background: rgba(255, 255, 255, 0.8);
        backdrop-filter: blur(10px);
        -webkit-backdrop-filter: blur(10px);
        border-radius: 8px;
        border: 2px solid rgba(110, 143, 178, 0.3);
        font-weight: 600;
        transition: all 0.2s ease;
        letter-spacing: 0.3px;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(110, 143, 178, 0.25);
        border-color: #6E8FB2;
        background: rgba(255, 255, 255, 0.95);
    }
    
    .stButton > button[kind="primary"] {
        background: linear-gradient(135deg, #6E8FB2 0%, #7DA494 100%);
        color: white;
        border: none;
    }
    
    .stButton > button[kind="primary"]:hover {
        background: linear-gradient(135deg, #5d7a98 0%, #6a8d7d 100%);
    }
    
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

def load_disclaimer():
    try:
        with open('disclaimer.md', 'r', encoding='utf-8') as file:
            return file.read()
    except FileNotFoundError:
        return """
# Educational Disclaimer

**IMPORTANT: FOR EDUCATIONAL AND RESEARCH PURPOSES ONLY**

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
    st.markdown("""
    <div style="text-align: center; padding: 2rem 0;">
        <h1 style="color: #495057;">AD PRS Genome Browser</h1>
        <p style="color: #6C757D; font-size: 1.1rem;">Educational Polygenic Risk Score Tool</p>
    </div>
    """, unsafe_allow_html=True)
    
    disclaimer_content = load_disclaimer()
    st.markdown(disclaimer_content)
    
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 1, 1])
    with col2:
        if st.button("I have read and agree to the terms", 
                    use_container_width=True, 
                    type="primary"):
            st.session_state.disclaimer_accepted = True
            st.rerun()

def render_snp_dropdown():
    st.markdown('<div class="section-header">SNP Selection</div>', unsafe_allow_html=True)
    
    snp_options = ["— Select SNP —"]
    for chrom in sorted(set(info['chromosome'] for info in SNP_DATA.values()), key=int):
        for rsid, info in SNP_DATA.items():
            if info['chromosome'] == chrom:
                badge = "RISK" if info['effect_weight'] > 0 else "PROT"
                snp_options.append(f"Chr{chrom} | {rsid} | {badge}")
    
    selected_option = st.selectbox(
        "Choose variant:",
        snp_options,
        key="snp_dropdown",
        label_visibility="collapsed"
    )
    
    if selected_option != "— Select SNP —":
        rsid = selected_option.split(" | ")[1]
        if rsid != st.session_state.get('selected_snp'):
            st.session_state.selected_snp = rsid
            st.rerun()

def render_compact_editor():
    selected_snp = st.session_state.get('selected_snp', None)
    
    if selected_snp and selected_snp in SNP_DATA:
        snp_info = SNP_DATA[selected_snp]
        current_genotype = st.session_state.genotypes[selected_snp]
        
        CHROMOSOME_COLORS = {
            '1': '#B5C3D7', '2': '#B6C9C0', '3': '#F5D8B7', '4': '#F1D0C6',
            '5': '#DDAEAB', '6': '#D0E0EF', '7': '#E9CDDF', '8': '#C8BFD9',
            '9': '#E4E3BF', '10': '#6E8FB2', '11': '#7DA494', '12': '#EAB67A',
            '13': '#E5A79A', '14': '#C16E71', '15': '#ABC8E5', '16': '#D8A0C1',
            '17': '#9F8DB8', '18': '#D0D63E', '19': '#B5C3D7', '20': '#B6C9C0',
            '21': '#F5D8B7', '22': '#F1D0C6'
        }
        
        chrom_color = CHROMOSOME_COLORS.get(snp_info['chromosome'], '#CCCCCC')
        
        st.markdown(f"""
        <div class="snp-header-glass" style="background: linear-gradient(135deg, {chrom_color} 0%, {chrom_color}dd 100%);">
            <h4 style="margin: 0; font-size: 1.1rem;">{selected_snp}</h4>
            <small style="opacity: 0.95;">Chr{snp_info['chromosome']} | {snp_info.get('locus_name', 'Unknown')}</small>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("**Genotype Selection**")
        options = get_genotype_options(snp_info['effect_allele'], snp_info['other_allele'])
        
        new_genotype = st.selectbox(
            "Select genotype:",
            options,
            index=options.index(current_genotype) if current_genotype in options else 0,
            key=f"edit_{selected_snp}"
        )
        
        st.markdown("**Effect Prediction**")
        for genotype in options:
            effect_count = genotype.count(snp_info['effect_allele'])
            effect_score = effect_count * snp_info['effect_weight']
            is_current = (genotype == new_genotype)
            
            if effect_count == 2:
                status = "HIGH"
                status_color = THEME_COLORS['danger']
            elif effect_count == 1:
                status = "MED"
                status_color = THEME_COLORS['warning']
            else:
                status = "NONE"
                status_color = THEME_COLORS['muted']
            
            bg_style = ""
            if is_current:
                bg_style = f"background: rgba(234, 182, 122, 0.25); border: 2px solid {THEME_COLORS['accent']};"
            else:
                bg_style = "background: rgba(255, 255, 255, 0.5); border: 1px solid #DEE2E6;"
            
            st.markdown(f"""
            <div class="genotype-option" style="{bg_style}">
                <strong>{genotype}</strong> 
                <span class="snp-badge badge-{'risk' if effect_count > 0 else 'protective'}">{status}</span>
                → <strong style="color: {status_color};">{effect_score:+.3f}</strong>
                {'← Current' if is_current else ''}
            </div>
            """, unsafe_allow_html=True)
        
        if new_genotype != current_genotype:
            if st.button("Apply Changes", use_container_width=True, type="primary"):
                st.session_state.genotypes[selected_snp] = new_genotype
                st.session_state.selected_snp = None
                st.success(f"{selected_snp} Updated")
                st.rerun()
        
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Randomize", use_container_width=True, help="MAF-based random genotype"):
                st.session_state.genotypes[selected_snp] = generate_realistic_genotype(selected_snp, snp_info)
                st.session_state.selected_snp = None
                st.rerun()
        with col2:
            if st.button("Close", use_container_width=True, help="Close editor"):
                st.session_state.selected_snp = None
                st.rerun()
    
    else:
        st.markdown(f"""
        <div class="glass-card" style="text-align: center; color: {THEME_COLORS['muted']};">
            <div class="section-header" style="border: none; text-align: center;">Genotype Editor</div>
            <p style="font-size: 14px; margin: 0;">Select SNP from dropdown</p>
        </div>
        """, unsafe_allow_html=True)

def create_percentile_chart():
    current_prs = calculate_prs(st.session_state.genotypes)
    
    POPULATION_MEAN = -2.4101
    POPULATION_STD = 0.6444
    THEORETICAL_MIN = -5.46
    THEORETICAL_MAX = 1.86
    
    from scipy.stats import norm
    percentile = norm.cdf(current_prs, loc=POPULATION_MEAN, scale=POPULATION_STD) * 100
    percentile = max(0.1, min(99.9, percentile))
    
    st.markdown('<div class="section-header">Population Percentile</div>', unsafe_allow_html=True)
    
    fig = go.Figure()
    
    x_min, x_max = THEORETICAL_MIN, THEORETICAL_MAX
    x_range = np.linspace(x_min, x_max, 500)
    y_normal = norm.pdf(x_range, loc=POPULATION_MEAN, scale=POPULATION_STD)
    
    fig.add_trace(go.Scatter(
        x=x_range,
        y=y_normal,
        mode='lines',
        line=dict(color=THEME_COLORS['info'], width=2),
        fill='tozeroy',
        fillcolor='rgba(171, 200, 229, 0.3)',
        name='European Population',
        showlegend=False
    ))
    
    if x_min <= current_prs <= x_max:
        user_y = norm.pdf(current_prs, loc=POPULATION_MEAN, scale=POPULATION_STD)
        
        fig.add_trace(go.Scatter(
            x=[current_prs, current_prs],
            y=[0, user_y],
            mode='lines',
            line=dict(color=THEME_COLORS['primary'], width=3),
            showlegend=False
        ))
        
        fig.add_trace(go.Scatter(
            x=[current_prs],
            y=[user_y],
            mode='markers',
            marker=dict(size=10, color=THEME_COLORS['primary'], symbol='circle'),
            showlegend=False
        ))
    
    fig.update_layout(
        xaxis=dict(title="PRS Score", showgrid=True, gridcolor='rgba(0,0,0,0.05)'),
        yaxis=dict(title="Density", showgrid=False),
        height=220,
        margin=dict(t=10, b=40, l=45, r=15),
        plot_bgcolor='rgba(255,255,255,0.5)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    
    st.plotly_chart(fig, use_container_width=True, key="percentile_chart")
    
    if percentile < 5:
        percentile_desc = "Very Low Risk"
        percentile_color = THEME_COLORS['success']
    elif percentile < 25:
        percentile_desc = "Below Average"
        percentile_color = THEME_COLORS['info']
    elif percentile < 75:
        percentile_desc = "Average Range"
        percentile_color = THEME_COLORS['muted']
    elif percentile < 95:
        percentile_desc = "Above Average"
        percentile_color = THEME_COLORS['warning']
    else:
        percentile_desc = "Very High Risk"
        percentile_color = THEME_COLORS['danger']
    
    st.markdown(f"""
    <div style="text-align: center; padding: 12px; margin-top: 8px;">
        <div style="display: inline-block; padding: 16px 24px; 
                    background: linear-gradient(135deg, rgba(255,255,255,0.9) 0%, rgba(248,249,250,0.8) 100%);
                    backdrop-filter: blur(10px);
                    border-radius: 12px;
                    border: 2px solid {percentile_color};
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08);">
            <div style="font-size: 13px; color: {THEME_COLORS['dark']}; margin-bottom: 4px;">
                <strong>Your Position</strong>
            </div>
            <div style="font-size: 26px; font-weight: bold; color: {percentile_color}; margin: 8px 0;">
                {percentile:.1f}<span style="font-size: 16px;">th</span>
            </div>
            <div style="font-size: 12px; color: {THEME_COLORS['muted']}; margin-bottom: 4px;">
                percentile
            </div>
            <div style="font-size: 13px; color: {percentile_color}; font-weight: 600; 
                        padding: 6px 12px; background: rgba(255,255,255,0.5); 
                        border-radius: 6px; margin-top: 8px;">
                {percentile_desc}
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

def render_control_panel():
    st.markdown('<div class="section-header">Global Controls</div>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("⟳ Randomize All", use_container_width=True):
            st.session_state.genotypes = generate_realistic_genotypes()
            st.session_state.selected_snp = None
            st.rerun()
    
    with col2:
        if st.button("↓ Minimize Risk", use_container_width=True):
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
        if st.button("↑ Maximize Risk", use_container_width=True):
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
    current_prs = calculate_prs(st.session_state.genotypes)
    
    effect_snps = sum(1 for rsid, genotype in st.session_state.genotypes.items() 
                     if SNP_DATA[rsid]['effect_allele'] in genotype)
    
    protective_snps = sum(1 for rsid, genotype in st.session_state.genotypes.items() 
                         if SNP_DATA[rsid]['effect_weight'] < 0 and 
                         SNP_DATA[rsid]['effect_allele'] in genotype)
    
    risk_snps = sum(1 for rsid, genotype in st.session_state.genotypes.items() 
                   if SNP_DATA[rsid]['effect_weight'] > 0 and 
                   SNP_DATA[rsid]['effect_allele'] in genotype)
    
    from scipy.stats import norm
    POPULATION_MEAN = -2.4101
    POPULATION_STD = 0.6444
    percentile = norm.cdf(current_prs, loc=POPULATION_MEAN, scale=POPULATION_STD) * 100
    percentile = max(0.1, min(99.9, percentile))
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.markdown(f"""
        <div class="stat-card-glass">
            <div class="stat-value" style="color: {THEME_COLORS['primary']}; font-size: 1.5rem; margin-bottom: 0.2rem;">{current_prs:.3f}</div>
            <div class="stat-label" style="font-size: 0.75rem;">PRS SCORE</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div class="stat-card-glass">
            <div class="stat-value" style="color: {THEME_COLORS['accent']}; font-size: 1.5rem; margin-bottom: 0.2rem;">{percentile:.1f}%</div>
            <div class="stat-label" style="font-size: 0.75rem;">PERCENTILE</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown(f"""
        <div class="stat-card-glass">
            <div class="stat-value" style="color: {THEME_COLORS['danger']}; font-size: 1.5rem; margin-bottom: 0.2rem;">{risk_snps}</div>
            <div class="stat-label" style="font-size: 0.75rem;">RISK SNPS</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div class="stat-card-glass">
            <div class="stat-value" style="color: {THEME_COLORS['success']}; font-size: 1.5rem; margin-bottom: 0.2rem;">{protective_snps}</div>
            <div class="stat-label" style="font-size: 0.75rem;">PROTECTIVE SNPS</div>
        </div>
        """, unsafe_allow_html=True)

def show_app_content():
    st.markdown("""
    <div class="warning-glass" style="margin-top: 0; margin-bottom: 0.8rem;">
        <small style="color: #856404;">
            <strong>Educational Tool Only</strong> – PRS of Late-Onset Alzheimer's Disease Visualization Tool: Not for medical use or clinical decision-making
        </small>
    </div>
    """, unsafe_allow_html=True)
    
    render_summary_stats()
    
    st.markdown("---")
    
    col_select, col_circos, col_right = st.columns([1, 2.5, 1])
    
    with col_select:
        render_snp_dropdown()
        st.markdown("---")
        create_percentile_chart()
    
    with col_circos:
        display_circos_in_streamlit(
            st.session_state.genotypes, 
            st.session_state.get('selected_snp', None)
        )
    
    with col_right:
        render_compact_editor()
        
        # 将 Global Controls 移到这里
        st.markdown("---")
        render_control_panel()
    
    st.markdown(f"""
        <div style="text-align: center; color: {THEME_COLORS['muted']}; 
             font-size: 12px; margin-top: 2rem; padding: 1rem;
             background: rgba(255, 255, 255, 0.5);
             backdrop-filter: blur(10px);
             border-radius: 8px;">
            <p><strong>Bocheng Shi</strong> • Student #81442386</p>
            <p style="font-size: 10px;">PSYC 301 Coursework Project • UBC 2025W1</p>
            <p style="font-size: 10px;">LOAD PRS Visualization | PGS000334 | 22 SNPs</p>
        </div>
    """, unsafe_allow_html=True)
      
def main():
    if 'disclaimer_accepted' not in st.session_state:
        st.session_state.disclaimer_accepted = False
    
    if not st.session_state.disclaimer_accepted:
        show_disclaimer_page()
    else:
        if 'genotypes' not in st.session_state:
            st.session_state.genotypes = generate_realistic_genotypes()
        if 'selected_snp' not in st.session_state:
            st.session_state.selected_snp = None
        
        show_app_content()

if __name__ == "__main__":
    main()