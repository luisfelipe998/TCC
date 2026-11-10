"""Independent arithmetic checks of printed source values; not velocity validation."""
import hashlib
import json
from pathlib import Path

import CoolProp
from CoolProp.CoolProp import PropsSI

ROOT = Path(__file__).resolve().parent
T = 300.15
RL = PropsSI('D', 'T', T, 'Q', 0, 'R600a')
RG = PropsSI('D', 'T', T, 'Q', 1, 'R600a')
SIGMA = PropsSI('I', 'T', T, 'Q', 0, 'R600a')
HLV = PropsSI('H', 'T', T, 'Q', 1, 'R600a') - PropsSI('H', 'T', T, 'Q', 0, 'R600a')

def ra(x, G):
    return (x/RG) / ((1+0.12*(1-x))*(x/RG+(1-x)/RL)
        + 1.18*(1-x)*(9.81*SIGMA*(RL-RG))**0.25/(G*RL**0.5))

table = [('bubbles small', .012, .36, .93, .7),
         ('bubbles larger', .025, .50, .91, 1.8),
         ('plug', .056, .68, .88, 2.4),
         ('slug', .076, .73, .86, 3.5),
         ('churn', .22, .85, .76, 10.5),
         ('wavy annular', .30, .88, .73, 11.7)]
data = {
    'interpretation': 'Equation/table consistency diagnostics only. Neither experimental accuracy nor a calibration of our records.',
    'CoolProp_version': CoolProp.__version__,
    'fluid': 'R600a', 'T_K': T, 'rho_l_kg_m3': RL, 'rho_g_kg_m3': RG,
    'sigma_N_m': SIGMA, 'h_lv_J_kg': HLV,
    'deoliveira2021_churn_range_check': {
        'source': 'Section 4.4, PDF p.17; Eq.17, PDF p.14',
        'printed_x_range': [.35, .62], 'printed_alpha_range': [.61, .85],
        'printed_S_range': [1.40, 2.34],
        'alpha_implied_by_x_and_S_min': 1/(1+2.34*((1-.35)/.35)*(RG/RL)),
        'alpha_implied_by_x_and_S_max': 1/(1+1.4*((1-.62)/.62)*(RG/RL)),
        'minimum_S_from_printed_x_and_alpha': .35*RL*(1-.85)/(.85*RG*(1-.35)),
    },
    'dias2025knn_table2_check_at_stated_G516_T27C': [
        {'pattern': name, 'printed_x': x, 'printed_alpha': a, 'RA_alpha': ra(x,516),
         'printed_jl_m_s': jl, 'calculated_jl_m_s': 516*(1-x)/RL,
         'printed_jg_m_s': jg, 'calculated_jg_m_s': 516*x/RG}
        for name,x,a,jl,jg in table],
    'Xu_Fang_secondary_source_Froude_check': [
        {'nominal_G_kg_m2_s': G, 'Fr_lo': G*G/(9.81*.001*RL*RL),
         'exceeds_secondary_source_Fr145_max': G*G/(9.81*.001*RL*RL)>145}
        for G in [380,510,640,760,890]],
    'rho_g_over_rho_l': RG/RL,
    'source_hashes_sha256': {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted((ROOT/'extracts').glob('*')) if p.is_file()},
}
(ROOT/'source_numeric_checks.json').write_text(json.dumps(data, indent=2)+'\n')
print(json.dumps(data, indent=2))
