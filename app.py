from collections import Counter
from io import StringIO

import numpy as np
import streamlit as st
import pandas as pd
from filterframes import from_dta_select_filter
from peptacular.sequence import strip_mods, convert_ip2_sequence
import matplotlib.pyplot as plt
from constants import COMPARISON_AA_FREQUENCIES


def get_aa_counts(peptides: list[str]) -> dict[str, int]:
    aa_counter = Counter()
    for peptide in peptides:
        aa_counter.update(peptide)  # Update with each peptide string
    return aa_counter

@st.cache_data()
def get_dta_select_counts(files: list[st.file_uploader]) -> dict[str, int]:
    peptides = []
    for file in files:
        file_io = StringIO(file.getvalue().decode("utf-8"))
        _, peptide_df, _, _ = from_dta_select_filter(file_io)
        peptides.extend([strip_mods(convert_ip2_sequence(peptide)) for peptide in peptide_df['Sequence'].tolist()])

    return get_aa_counts(peptides)


def calculate_log2fold_change(observed_freqs: dict[str, float], baseline_freqs: dict[str, float]):
    all_aa = set(observed_freqs.keys()).union(baseline_freqs.keys())
    log2fold_changes = {}
    for aa in all_aa:
        observed_freq = observed_freqs.get(aa, 0)
        baseline_freq = baseline_freqs.get(aa, 0)
        if baseline_freq == 0:
            log2fold_changes[aa] = np.inf
        else:
            log2fold_changes[aa] = np.log2(observed_freq / baseline_freq)

    return log2fold_changes


with st.sidebar:
    st.title("AA Freq Explorer")

    filter_files = st.file_uploader("Upload DtaSelect-filter.txt files", accept_multiple_files=True, type='.txt')
    baseline_aa_freq = COMPARISON_AA_FREQUENCIES[st.selectbox('Select baseline AA frequency', COMPARISON_AA_FREQUENCIES.keys())]

    run = st.button('Run', use_container_width=True, type='primary')

if run:

    st.subheader('Amino Acid Frequencies Comparison')

    if not filter_files:
        st.write('No files uploaded')
        st.stop()

    aa_counts = get_dta_select_counts(filter_files)

    # add all aa's from baseline to aa_freqs
    for aa in baseline_aa_freq:
        if aa not in aa_counts:
            aa_counts[aa] = 0

    aa_to_remove = set()
    for aa in aa_counts:
        if aa not in baseline_aa_freq:
            aa_to_remove.add(aa)

    for aa in aa_to_remove:
        del aa_counts[aa]

    total_aa = sum(aa_counts.values())
    aa_freqs = {aa: count / total_aa for aa, count in aa_counts.items()}

    log2fold_changes = calculate_log2fold_change(aa_freqs, baseline_aa_freq)

    df = pd.DataFrame.from_dict(log2fold_changes, orient='index', columns=['Log2Fold Change'])
    df['Observed Frequency'] = df.index.map(aa_freqs.get)
    df['Baseline Frequency'] = df.index.map(baseline_aa_freq.get)
    df['Observed Count'] = df.index.map(aa_counts.get)
    df['Expected Count'] = df['Baseline Frequency'] * total_aa
    df['Expected Count'] = df['Expected Count'].astype(int)


    # fill inf
    df.replace([np.inf, -np.inf], np.nan, inplace=True)


    # Apply a gradient color map to the Log2Fold Change column
    styled_df = df.style.background_gradient(subset=['Log2Fold Change'], cmap='coolwarm')

    # Display the styled dataframe in Streamlit
    st.dataframe(styled_df, column_config={
        'Log2Fold Change': st.column_config.NumberColumn(format="%.2f"),
    }, use_container_width=True)
    st.subheader('AA Frequency Log2Fold Change')
    st.bar_chart(df['Log2Fold Change'], x_label='Amino Acid', y_label='Log2Fold Change')

