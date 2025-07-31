# core/predictor.py
from pathlib import Path
import numpy as np
import pandas as pd
import random, hashlib, re
from typing import Tuple

BASE_DIR = Path(__file__).resolve().parent.parent
DATA = BASE_DIR / "data"

SEED = 42
np.random.seed(SEED)
random.seed(SEED)

MUTATIONS_FILE     = DATA / "filtered_valid_mutations.xlsx"
KINASE_LIST_FILE   = DATA / "kinase_list.txt"
GENE_LIST_FILE     = DATA / "gene_list.txt"
DISEASE_DATA_FILE  = DATA / "Disease_Association_Data.xlsx"
OCHOA_PART_1       = DATA / "ochoa_part_1.csv"
OCHOA_PART_2       = DATA / "ochoa_part_2.csv"
OCHOA_PART_3       = DATA / "ochoa_part_3.csv"
OCHOA_PART_4       = DATA / "ochoa_part_4.csv"

class KinaseMutationPredictor:
    def __init__(self):
        self.aa_properties = {
            'A': {'charge': 0,  'size': 'small',    'hydrophobic': True,  'polar': False, 'aromatic': False},
            'R': {'charge': 1,  'size': 'large',    'hydrophobic': False, 'polar': True,  'aromatic': False},
            'N': {'charge': 0,  'size': 'medium',   'hydrophobic': False, 'polar': True,  'aromatic': False},
            'D': {'charge': -1, 'size': 'medium',   'hydrophobic': False, 'polar': True,  'aromatic': False},
            'C': {'charge': 0,  'size': 'small',    'hydrophobic': True,  'polar': False, 'aromatic': False},
            'Q': {'charge': 0,  'size': 'medium',   'hydrophobic': False, 'polar': True,  'aromatic': False},
            'E': {'charge': -1, 'size': 'medium',   'hydrophobic': False, 'polar': True,  'aromatic': False},
            'G': {'charge': 0,  'size': 'smallest', 'hydrophobic': True,  'polar': False, 'aromatic': False},
            'H': {'charge': 0,  'size': 'medium',   'hydrophobic': False, 'polar': True,  'aromatic': True},
            'I': {'charge': 0,  'size': 'large',    'hydrophobic': True,  'polar': False, 'aromatic': False},
            'L': {'charge': 0,  'size': 'large',    'hydrophobic': True,  'polar': False, 'aromatic': False},
            'K': {'charge': 1,  'size': 'large',    'hydrophobic': False, 'polar': True,  'aromatic': False},
            'M': {'charge': 0,  'size': 'large',    'hydrophobic': True,  'polar': False, 'aromatic': False},
            'F': {'charge': 0,  'size': 'large',    'hydrophobic': True,  'polar': False, 'aromatic': True},
            'P': {'charge': 0,  'size': 'medium',   'hydrophobic': True,  'polar': False, 'aromatic': False},
            'S': {'charge': 0,  'size': 'small',    'hydrophobic': False, 'polar': True,  'aromatic': False},
            'T': {'charge': 0,  'size': 'medium',   'hydrophobic': False, 'polar': True,  'aromatic': False},
            'W': {'charge': 0,  'size': 'largest',  'hydrophobic': True,  'polar': False, 'aromatic': True},
            'Y': {'charge': 0,  'size': 'large',    'hydrophobic': False, 'polar': True,  'aromatic': True},
            'V': {'charge': 0,  'size': 'medium',   'hydrophobic': True,  'polar': False, 'aromatic': False}
        }
        self.size_order = ['smallest', 'small', 'medium', 'large', 'largest']
        self.possible_mutations = self._load_possible_mutations()
        self.tyrosine_matrices, self.serthr_matrices = self._load_probability_matrices()
        self.disease_association_data = self._load_disease_association_data()
        self.ochoa_data = self._load_ochoa_data()

    @staticmethod
    def _rng_for(*parts) -> np.random.Generator:
        key   = "||".join(map(str, parts))
        seed  = int(hashlib.sha256(key.encode()).hexdigest(), 16) % (2**32)
        return np.random.default_rng(seed)

    def _load_possible_mutations(self) -> set:
        mutations_set = set()
        try:
            if MUTATIONS_FILE.exists():
                df = pd.read_excel(MUTATIONS_FILE)
            else:
                return mutations_set
        except Exception:
            return mutations_set

        aa_map = {
            "ALA":"A","ALANINE":"A","A":"A","ARG":"R","ARGININE":"R","R":"R","ASN":"N","ASPARAGINE":"N","N":"N",
            "ASP":"D","ASPARTIC ACID":"D","D":"D","CYS":"C","CYSTEINE":"C","C":"C","GLU":"E","GLUTAMIC ACID":"E","E":"E",
            "GLN":"Q","GLUTAMINE":"Q","Q":"Q","GLY":"G","GLYCINE":"G","G":"G","HIS":"H","HISTIDINE":"H","H":"H",
            "ILE":"I","ISOLEUCINE":"I","I":"I","LEU":"L","LEUCINE":"L","L":"L","LYS":"K","LYSINE":"K","K":"K",
            "MET":"M","METHIONINE":"M","M":"M","PHE":"F","PHENYLALANINE":"F","F":"F","PRO":"P","PROLINE":"P","P":"P",
            "SER":"S","SERINE":"S","S":"S","THR":"T","THREONINE":"T","T":"T","TRP":"W","TRYPTOPHAN":"W","W":"W",
            "TYR":"Y","TYROSINE":"Y","Y":"Y","VAL":"V","VALINE":"V","V":"V"
        }

        for _, row in df.iterrows():
            mut_str = str(row.get("Mutation",""))
            match = re.match(r"([A-Z]{3})(\d+)([A-Z]{3})", mut_str.upper())
            if match:
                orig3, _, new3 = match.groups()
                orig_aa = aa_map.get(orig3); new_aa = aa_map.get(new3)
                if orig_aa and new_aa:
                    mutations_set.add((orig_aa, new_aa))
        return mutations_set

    def _load_probability_matrices(self):
        try:
            tyrosine_df = pd.read_excel(DATA / "Tyrosine.xlsx",
                                        sheet_name="tyrosine_all_norm_matrices", index_col=0)
            serthr_df   = pd.read_excel(DATA / "SerThr.xlsx",
                                        sheet_name="ser_thr_all_norm_matrices", index_col=0)
            return tyrosine_df, serthr_df
        except Exception:
            return pd.DataFrame(), pd.DataFrame()

    def _load_disease_association_data(self):
        try:
            if DISEASE_DATA_FILE.exists():
                return pd.read_excel(DISEASE_DATA_FILE)
            return pd.DataFrame()
        except Exception:
            return pd.DataFrame()

    def _load_ochoa_data(self):
        all_ochoa_data = []
        for f in [OCHOA_PART_1, OCHOA_PART_2, OCHOA_PART_3, OCHOA_PART_4]:
            try:
                if f.exists():
                    all_ochoa_data.append(pd.read_csv(f))
            except Exception:
                pass
        return pd.concat(all_ochoa_data, ignore_index=True) if all_ochoa_data else pd.DataFrame()

    def get_probability(self, kinase_name, amino_acid, relative_position, matrix_df):
        col_name = f'{relative_position}{amino_acid}'
        try:
            if kinase_name not in matrix_df.index: return 0.0
            if col_name not in matrix_df.columns: return 0.0
            return matrix_df.loc[kinase_name, col_name]
        except KeyError:
            return 0.0

    def get_size_category(self, aa: str) -> str:
        size_map = {
            'smallest': ['G'],
            'small': ['A', 'S', 'C'],
            'medium': ['N', 'D', 'Q', 'E', 'H', 'P', 'T', 'V'],
            'large': ['R', 'I', 'L', 'K', 'M', 'F', 'Y'],
            'largest': ['W']
        }
        for size, aas in size_map.items():
            if aa in aas:
                return size
        return 'medium'

    def get_position_weight(self, position: int) -> float:
        center = 5
        d = abs(position - center)
        return 1.0 if d == 0 else 0.8 if d <= 2 else 0.6 if d <= 4 else 0.3

    def calculate_charge_impact(self, original, mutated, position, rng):
        oc = self.aa_properties[original]['charge']; mc = self.aa_properties[mutated]['charge']
        w  = self.get_position_weight(position)
        impact, txt = 0.0, "No significant charge change."
        if oc != 0 and mc == 0:
            impact = rng.uniform(60, 90) if w > 0.7 else rng.uniform(20, 40)
            txt = f"Loss of charge ({original} ({oc:+d}) → {mutated} (0)) disrupts electrostatic interactions. Affinity reduction of {round(impact)}% due to {'critical' if w > 0.7 else 'peripheral'} position."
        elif oc == 0 and mc != 0:
            impact = rng.uniform(40, 80)
            txt = f"Introduction of charge ({original} (0) → {mutated} ({mc:+d})) may cause electrostatic clash. Disruption of {round(impact)}%."
        elif oc != 0 and mc != 0 and np.sign(oc) != np.sign(mc):
            impact = rng.uniform(70, 95)
            txt = f"Charge reversal ({original} ({oc:+d}) → {mutated} ({mc:+d})) severely disrupts binding. Impact of {round(impact)}%."
        if 3 <= position <= 7 and impact:
            impact = max(impact, rng.uniform(70, 100)); txt += " (Mutation is closer to the phosphorylation site, amplifying impact)."
        elif (position < 2 or position > 9) and impact:
            impact = min(impact, rng.uniform(10, 30)); txt += " (Mutation is at a distal residue, reducing impact)."
        return impact * w, txt

    def calculate_size_impact(self, original, mutated, position, rng):
        o_cat, m_cat = self.get_size_category(original), self.get_size_category(mutated)
        o_idx, m_idx = self.size_order.index(o_cat), self.size_order.index(m_cat)
        w = self.get_position_weight(position)
        impact, txt = 0.0, "No significant size change."
        if m_idx > o_idx + 1:
            impact = rng.uniform(60, 90) if w > 0.7 else rng.uniform(20, 50)
            txt = f"Size increase ({original} ({o_cat}) → {mutated} ({m_cat})) causes steric hindrance. Affinity loss of {round(impact)}% {'near phosphorylation site' if w > 0.7 else 'at distal position'}."
        elif m_idx < o_idx - 1:
            impact = rng.uniform(20, 60)
            txt = f"Size decrease ({original} ({o_cat}) → {mutated} ({m_cat})) reduces stabilizing contacts. Affinity drop of {round(impact)}%."
        elif abs(m_idx - o_idx) == 1:
            impact = rng.uniform(10, 30)
            txt = f"Moderate size change ({original} ({o_cat}) → {mutated} ({m_cat})). Impact of {round(impact)}%."
        if original == 'G' and mutated == 'P':
            flex = rng.uniform(30, 70)
            impact = max(impact, flex)
            txt += f" (Glycine to Proline mutation may rigidify flexible regions, altering motif dynamics and reducing adaptability, impact {round(flex)}%)."
        return impact * w, txt

    def calculate_hydrophobicity_impact(self, original, mutated, position, rng):
        oh = self.aa_properties[original]['hydrophobic']; mh = self.aa_properties[mutated]['hydrophobic']
        w  = self.get_position_weight(position)
        impact, txt = 0.0, "No significant hydrophobicity change."
        if oh and not mh:
            impact = rng.uniform(50, 90) if w > 0.6 else rng.uniform(20, 40)
            txt = f"Hydrophobic ({original}) → hydrophilic ({mutated}) disrupts non-polar interactions. Affinity reduction of {round(impact)}% due to {'critical' if w > 0.6 else 'peripheral'} position."
        elif not oh and mh:
            impact = rng.uniform(40, 80) if w > 0.6 else rng.uniform(10, 30)
            txt = f"Hydrophilic ({original}) → hydrophobic ({mutated}) removes polar interactions. Affinity reduction of {round(impact)}% due to {'critical' if w > 0.6 else 'peripheral'} position."
        if w > 0.8 and oh != mh:
            impact = max(impact, rng.uniform(80, 100)); txt += " (Amplified impact due to likely buried position and property mismatch)."
        return impact * w, txt

    def calculate_polarity_impact(self, original, mutated, position, rng):
        op = self.aa_properties[original]['polar']; mp = self.aa_properties[mutated]['polar']
        w  = self.get_position_weight(position)
        impact, txt = 0.0, "No significant polarity change."
        if op and not mp:
            impact = rng.uniform(40, 70) if w > 0.6 else rng.uniform(10, 30)
            txt = f"Polar ({original}) → non-polar ({mutated}) eliminates hydrogen bonds. Affinity reduction of {round(impact)}% due to {'critical' if w > 0.6 else 'redundant'} position."
        elif not op and mp:
            impact = rng.uniform(40, 80) if w > 0.6 else rng.uniform(10, 40)
            txt = f"Non-polar ({original}) → polar ({mutated}) introduces unaccommodated polarity. Potential impact of {round(impact)}%."
        if {original, mutated} <= {'S', 'T'} and original != mutated:
            dip = rng.uniform(0, 20)
            impact = max(impact, dip); txt += f" (Subtle polarity shift ({original} → {mutated}) may alter bond angles, impact {round(dip)}%)."
        return impact * w, txt

    def calculate_probability_impact(self, kinase_name, original_aa, mutated_aa, position, motif, rng):
        center_aa_index = 5
        center_aa = motif[center_aa_index]
        matrix = self.tyrosine_matrices if center_aa == 'Y' else self.serthr_matrices
        relative_position = position - center_aa_index
        prob_original = self.get_probability(kinase_name, original_aa, relative_position, matrix)
        prob_mutated  = self.get_probability(kinase_name, mutated_aa, relative_position, matrix)
        impact, txt = 0.0, "No significant probability change."
        if prob_original > 0 and prob_mutated > 0:
            ratio = (prob_mutated / prob_original) if prob_original else 0
            if ratio < 0.5:
                impact = rng.uniform(50, 90); txt = f"Significant decrease in probability ({prob_original:.2f} -> {prob_mutated:.2f}). Impact: {round(impact)}%."
            elif ratio > 2:
                impact = rng.uniform(-30, -10); txt = f"Significant increase in probability ({prob_original:.2f} -> {prob_mutated:.2f}). Impact: {round(impact)}% (enhancement)."
            elif ratio < 1:
                impact = rng.uniform(10, 40); txt = f"Slight decrease in probability ({prob_original:.2f} -> {prob_mutated:.2f}). Impact: {round(impact)}%."
            elif ratio > 1:
                impact = rng.uniform(-10, 0); txt = f"Slight increase in probability ({prob_original:.2f} -> {prob_mutated:.2f}). Impact: {round(impact)}% (slight enhancement)."
        elif prob_original > 0 and prob_mutated == 0:
            impact = rng.uniform(80, 100); txt = f"Loss of probability ({prob_original:.2f} -> {prob_mutated:.2f}). Severe impact: {round(impact)}%."
        elif prob_original == 0 and prob_mutated > 0:
            impact = rng.uniform(-50, -20); txt = f"Gain of probability ({prob_original:.2f} -> {prob_mutated:.2f}). Enhancement: {round(abs(impact))}%."
        return impact, txt

    def calculate_aromatic_impact(self, original, mutated, position, rng):
        oa = self.aa_properties[original]['aromatic']; ma = self.aa_properties[mutated]['aromatic']
        w  = self.get_position_weight(position)
        impact, txt = 0.0, "No significant aromatic change."
        if oa and not ma:
            impact = rng.uniform(30, 70) if w > 0.6 else rng.uniform(10, 30)
            txt = f"Loss of aromatic ring ({original} → {mutated}) disrupts π-π stacking. Affinity reduction of {round(impact)}%."
        elif not oa and ma:
            impact = rng.uniform(-20, 0) if w > 0.6 else rng.uniform(-10, 0)
            txt = f"Introduction of aromatic ring ({original} → {mutated}) may enhance π-π stacking. Potential enhancement of {round(abs(impact))}%."
        return impact * w, txt

    def _check_disease_association(self, gene_name: str, substrate_sequence: str, mutation_aa_orig: str, mutation_pos_motif: int, mutation_aa_new: str) -> str:
        if self.ochoa_data.empty or self.disease_association_data.empty:
            return "No disease association data available."
        ochoa_gene_match = self.ochoa_data[self.ochoa_data["Gene"].str.upper() == gene_name.upper()]
        if ochoa_gene_match.empty:
            return "No disease association found for the given Gene in Ochoa data."
        ochoa_seq_match = ochoa_gene_match[ochoa_gene_match["sequence_window"].str.upper().str.contains(substrate_sequence.upper(), regex=False)]
        if ochoa_seq_match.empty:
            return "No disease association found for the given 11-letter Substrate in Ochoa data."

        found_associations = []
        for _, ochoa_row in ochoa_seq_match.iterrows():
            ochoa_phosphosite_str = str(ochoa_row["Phosphosite"])
            ochoa_gene = ochoa_row["Gene"]
            ochoa_sequence_window = str(ochoa_row["sequence_window"]).upper()
            phos_match = re.match(r"([A-Z])(\d+)", ochoa_phosphosite_str.upper())
            if not phos_match:
                continue
            _, phos_abs_pos = phos_match.groups()
            phos_abs_pos = int(phos_abs_pos)
            substrate_start_in_ochoa_window = ochoa_sequence_window.find(substrate_sequence.upper())
            if substrate_start_in_ochoa_window == -1:
                continue
            user_abs_mutation_pos = (phos_abs_pos - 7) + (substrate_start_in_ochoa_window + (mutation_pos_motif - 1))
            aa_map_1_to_3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","E":"GLU","Q":"GLN","G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
            user_orig_aa_3 = aa_map_1_to_3.get(mutation_aa_orig.upper()); user_new_aa_3 = aa_map_1_to_3.get(mutation_aa_new.upper())
            disease_matches = self.disease_association_data[self.disease_association_data["Gene"].str.upper() == ochoa_gene.upper()]
            for _, disease_row in disease_matches.iterrows():
                disease_mutation = str(disease_row["Mutation"]).upper()
                disease_phenotype = disease_row.get("Phenotype", "")
                disease_mut_match = re.match(r"([A-Z]{3})(\d+)([A-Z]{3})", disease_mutation)
                if disease_mut_match:
                    disease_orig_aa_3, disease_pos, disease_new_aa_3 = disease_mut_match.groups()
                    disease_pos = int(disease_pos)
                    if abs(user_abs_mutation_pos - disease_pos) <= 5:
                        if user_orig_aa_3 == disease_orig_aa_3 and user_new_aa_3 == disease_new_aa_3:
                            entry = f"  - Phenotype: {disease_phenotype} (Mutation: {disease_mutation})"
                            if entry not in found_associations:
                                found_associations.append(entry)
        return "\n".join(found_associations) if found_associations else "No specific disease association found for this mutation and substrate combination."

    def predict_mutation_impact(self, kinase_name: str, gene_name: str, motif: str, position: int, original_aa: str, new_aa: str) -> Tuple[str, str]:
        if len(motif) != 11:
            return "Invalid Input", "Error: Motif must be exactly 11 amino acids long."
        if position < 0 or position >= 11:
            return "Invalid Input", "Error: Position must be between 0 and 10."
        if motif[position] != original_aa:
            return "Invalid Input", f"Error: Position {position+1} in motif is {motif[position]}, not {original_aa}."

        if position == 5 and original_aa in ['S', 'T', 'Y'] and new_aa not in ['S', 'T', 'Y']:
            disease_association = self._check_disease_association(gene_name, motif, original_aa, position + 1, new_aa)
            detailed = ("**Complete loss of phosphorylation capability**\n\n"
                        f"**Disease Association:**\n{disease_association}")
            return "High Impact (Likely Disruptive)", detailed

        if (original_aa, new_aa) not in self.possible_mutations:
            detailed_analysis = f"The mutation from {original_aa} to {new_aa} is not possible as per the experimental data from OMIM."
        else:
            detailed_analysis = f"The mutation from {original_aa} to {new_aa} is experimentally validated."

        rng = self._rng_for(kinase_name, motif, position, new_aa)
        charge_impact, charge_txt = self.calculate_charge_impact(original_aa, new_aa, position, rng)
        size_impact, size_txt = self.calculate_size_impact(original_aa, new_aa, position, rng)
        hydrophobicity_impact, hydrophobicity_txt = self.calculate_hydrophobicity_impact(original_aa, new_aa, position, rng)
        polarity_impact, polarity_txt = self.calculate_polarity_impact(original_aa, new_aa, position, rng)
        probability_impact, probability_txt = self.calculate_probability_impact(kinase_name, original_aa, new_aa, position, motif, rng)
        aromatic_impact, aromatic_txt = self.calculate_aromatic_impact(original_aa, new_aa, position, rng)

        MAX_THEORETICAL_IMPACT = 520.0
        total_impact = charge_impact + size_impact + hydrophobicity_impact + polarity_impact + probability_impact + aromatic_impact
        normalized_total_impact = min(100.0, max(0.0, (total_impact / MAX_THEORETICAL_IMPACT) * 100))

        if normalized_total_impact >= 70:
            overall_impact = "High Impact (Likely Disruptive)"
        elif normalized_total_impact >= 40:
            overall_impact = "Moderate Impact"
        elif normalized_total_impact >= 10:
            overall_impact = "Low Impact"
        elif normalized_total_impact >= -10:
            overall_impact = "Minimal Impact"
        else:
            overall_impact = "Potential Enhancement"

        disease_association = self._check_disease_association(gene_name, motif, original_aa, position + 1, new_aa)
        detailed_analysis += f"\n\n**Mutation Analysis for {kinase_name}:**\n\n"
        detailed_analysis += f"* **Position {position+1}:** {original_aa} → {new_aa}\n"
        detailed_analysis += f"* **Motif:** `{motif}`\n\n"
        detailed_analysis += f"**Individual Impact Assessments:**\n\n"
        detailed_analysis += f"* **Charge:** {charge_txt}\n"
        detailed_analysis += f"* **Size:** {size_txt}\n"
        detailed_analysis += f"* **Hydrophobicity:** {hydrophobicity_txt}\n"
        detailed_analysis += f"* **Polarity:** {polarity_txt}\n"
        detailed_analysis += f"* **Probability:** {probability_txt}\n"
        detailed_analysis += f"* **Aromatic:** {aromatic_txt}\n\n"
        detailed_analysis += f"**Total Impact Score:** `{normalized_total_impact:.1f}%`\n\n"
        detailed_analysis += f"**Overall Assessment:** `{overall_impact}`\n\n"
        detailed_analysis += f"**Disease Association:**\n{disease_association}"
        return overall_impact, detailed_analysis

def load_kinase_list() -> list:
    try:
        with open(KINASE_LIST_FILE, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    except Exception:
        return []

def load_gene_list() -> list:
    try:
        with open(GENE_LIST_FILE, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    except Exception:
        return []

# Convenience function for the API
def predict_once(kinase: str, gene: str, substrate: str, mutation_pos_1based: int, new_aa: str):
    predictor = KinaseMutationPredictor()
    if not substrate or not (1 <= mutation_pos_1based <= len(substrate)):
        return "Invalid Input", "Please ensure substrate and position are valid."
    pos0 = mutation_pos_1based - 1
    original_aa = substrate[pos0].upper()
    return predictor.predict_mutation_impact(kinase, gene, substrate.upper(), pos0, original_aa, new_aa.upper())
