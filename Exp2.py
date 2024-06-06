import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import numpy as np
from scipy.signal import find_peaks

# Load the Excel file and specify the sheet name
file_path = r"C:\Users\17244\OneDrive\Grad School\NUCE 597\Exp2\Data.xlsx"
df = pd.read_excel(file_path, sheet_name='Sheet1')

count_time = 300  # Count time in seconds
target_energy = 1332.50  # Target energy in keV

# Simplified reference list of isotopes and their gamma energies (in keV)
isotope_reference = {
    "Co-57": [14.41, 122.06, 136.47],
    "Sn-119m": [23.87],
    "Th-233m": [29.78],
    "Pb-199m": [35.49],
    "I-125": [35.49],
    "Te-127m": [41.89],
    "Tb-161": [48.92],
    "Sn-113": [49.53],
    "I-125m": [49.53],
    "Pb-203": [51.64],
    "Ba-133": [53.16, 79.61, 276.40, 302.85, 356.01, 383.85],
    "Rh-104m": [77.53, 306.62],
    "W-187": [68.74, 122.06, 136.47],
    "I-131": [80.19, 364.49, 637.00, 723.30, 772.58],
    "Pt-195m": [130.89],
    "Hg-203": [279.20],
    "Ra-226": [186.21],
    "Bi-214": [609.31, 1120.29, 1238.11, 1377.67, 1764.49, 2204.21],
    "Tl-208": [583.19, 2614.53],
    "Cd-115m": [335.20],
    "Y-91": [687.67, 888.02],
    "U-238": [49.55, 92.38, 1001.03],
    "Th-234": [63.29],
    "Pa-234m": [1001.03],
    "Pb-214": [295.22, 351.93],
    "Cs-137": [661.66],
    "Co-60": [1173.23, 1332.50],
    "Am-241": [59.54],
    "Na-22": [511.00, 1274.53],
    "Eu-152": [121.78, 244.70, 344.28, 411.12, 444.01, 778.90, 964.08, 1085.87, 1112.07, 1408.01],
    "Mn-54": [834.85],
    "Zn-65": [1115.54],
    "K-40": [1460.82],
    "Fe-59": [1099.25],
    "Co-58": [810.77],
    "Ir-192": [316.51, 468.07, 604.41],
    "In-111": [171.28, 245.35, 508.05],
    "Ga-67": [93.31, 184.58, 296.97, 393.53],
    "Tc-99m": [140.51],
    "I-123": [159.00],
    "Xe-133": [80.99, 161.63, 223.21, 528.27],
    "Ce-141": [145.44],
    "Cs-134": [569.31, 604.70, 795.85],
    "Zn-69m": [438.63, 907.60],
    "W-188": [48.2],
    "Pd-109": [88.03],
    "Ho-166": [80.6],
    "Re-186": [137.15],
    "Tm-171": [66.0],
    "Gd-153": [41.0, 102.0],
    "Sn-113": [28.0],
    "Sc-46": [889.3],
    "Na-22": [511.0],
    "F-18": [511.0],
    "Ga-68": [1077.3],
    "Ba-137m": [661.7],
    "Cr-51": [320.08],
    "Bi-210m": [205.0],
    "Y-88": [898.0],
    "Au-198": [411.8],
    "Xe-133m": [233.2],
    "Mo-99": [140.51],
    "Tc-99m": [140.51],
    "Hg-203": [279.2],
    "Ra-224": [240.99],
    "Ac-228": [911.20, 968.97],
    "Th-228": [84.39, 129.07, 231.0, 338.3, 911.2, 964.08, 1588.19, 2614.53],
    "Ra-224": [240.99, 241.0],
    "Rn-220": [510.6],
    "Po-216": [146.3],
    "Pb-212": [238.6, 300.1, 510.8],
    "Bi-212": [727.3, 785.4, 860.6, 1620.5],
    "Tl-208": [583.19, 860.6, 2614.53]
}

def calculate_statistics(counts, count_time):
    count_rate = counts / count_time
    count_uncertainty = np.sqrt(counts.values)
    count_rate_uncertainty = count_uncertainty / count_time
    
    stats_df = pd.DataFrame({
        'Counts': counts.values,
        'Count Rate (counts/sec)': count_rate.values,
        'Count Uncertainty': count_uncertainty,
        'Count Rate Uncertainty (counts/sec)': count_rate_uncertainty
    })
    return stats_df

def find_closest_energy_row(df, energy_col, target_energy):
    energies = pd.to_numeric(df[energy_col][2:], errors='coerce')
    closest_index = (energies - target_energy).abs().idxmin()
    return closest_index

def compile_statistics_for_target_energy(df, target_energy, count_time):
    detector_starts = {
        'CZT Kromtek': 0,
        'Falcon 5000': 8,
        'Identifinder': 16,
        'Inspector 1000': 24,
        'Microdetective': 32,
        'RIIDE': 40
    }
    
    compiled_stats = []

    for detector_name, start_col in detector_starts.items():
        energy_col = df.columns[start_col + 1]
        closest_index = find_closest_energy_row(df, energy_col, target_energy)
        
        sample_col = df.columns[start_col + 3]  # Co-60 is the second sample column
        counts = pd.to_numeric(df.at[closest_index, sample_col], errors='coerce')
        stats = calculate_statistics(pd.Series([counts]), count_time).loc[0]
        stats['Detector'] = detector_name
        stats['Sample'] = 'Co-60'
        stats['Energy (keV)'] = df.at[closest_index, energy_col]
        compiled_stats.append(stats)
    
    compiled_df = pd.DataFrame(compiled_stats)
    return compiled_df

def calculate_energy_resolution(df, energy_col, sample_col, target_energy):
    energy = pd.to_numeric(df[energy_col][2:], errors='coerce').reset_index(drop=True)
    counts = pd.to_numeric(df[sample_col][2:], errors='coerce').reset_index(drop=True)
    
    closest_index = (energy - target_energy).abs().idxmin()
    peak_energy = energy.at[closest_index]
    peak_count = counts.at[closest_index]
    half_max = peak_count / 2
    
    try:
        left_half_max_indices = np.where(counts[:closest_index] <= half_max)[0]
        right_half_max_indices = np.where(counts[closest_index:] <= half_max)[0]
        
        if len(left_half_max_indices) > 0 and len(right_half_max_indices) > 0:
            left_half_max = left_half_max_indices[-1]
            right_half_max = right_half_max_indices[0] + closest_index
            fwhm = energy.at[right_half_max] - energy.at[left_half_max]
        else:
            fwhm = np.nan
    except Exception as e:
        fwhm = np.nan  # Return NaN if FWHM cannot be calculated
    
    return peak_energy, fwhm

def compare_detectors(df, target_energy, count_time):
    detector_starts = {
        'CZT Kromtek': 0,
        'Falcon 5000': 8,
        'Identifinder': 16,
        'Inspector 1000': 24,
        'Microdetective': 32,
        'RIIDE': 40
    }
    
    results = []

    for detector_name, start_col in detector_starts.items():
        energy_col = df.columns[start_col + 1]
        sample_col = df.columns[start_col + 3]  # Co-60 is the second sample column
        
        counts = pd.to_numeric(df[sample_col][2:], errors='coerce')
        closest_index = find_closest_energy_row(df, energy_col, target_energy)
        stats = calculate_statistics(counts, count_time).loc[closest_index - 2]
        peak_energy, fwhm = calculate_energy_resolution(df, energy_col, sample_col, target_energy)
        
        result = {
            'Detector': detector_name,
            'Peak Energy (keV)': peak_energy,
            'FWHM (keV)': fwhm,
            'Count Rate (counts/sec)': stats['Count Rate (counts/sec)'],
            'Count Rate Uncertainty (counts/sec)': stats['Count Rate Uncertainty (counts/sec)']
        }
        results.append(result)
    
    results_df = pd.DataFrame(results)
    return results_df

def find_prominent_peaks(counts, energy, num_peaks=5):
    peaks, properties = find_peaks(counts, height=counts.max()/10)
    prominent_peaks_indices = peaks[np.argsort(properties['peak_heights'])][-num_peaks:]
    prominent_peaks_energies = energy[prominent_peaks_indices]
    return prominent_peaks_energies, prominent_peaks_indices

def match_peaks_to_isotopes(peaks_energies, isotope_reference, tolerance=5.0):
    matched_isotopes = {}
    for peak in peaks_energies:
        for isotope, energies in isotope_reference.items():
            for energy in energies:
                if abs(peak - energy) <= tolerance:
                    if peak not in matched_isotopes:
                        matched_isotopes[peak] = []
                    matched_isotopes[peak].append(isotope)
    return matched_isotopes

isotope_results = []
def parse_and_plot_all(df, normalize=False):
    detector_starts = {
        'CZT Kromtek': 0,
        'Falcon 5000': 8,
        'Identifinder': 16,
        'Inspector 1000': 24,
        'Microdetective': 32,
        'RIIDE': 40
    }
    
    sample_types = ['Background', 'Co-60', 'Monazite', 'Pitchblende', 'Uslug']
    
    for sample_type in sample_types:
        traces = []
        for detector_name, start_col in detector_starts.items():
            energy_col = df.columns[start_col + 1]
            sample_col = df.columns[start_col + 2 + sample_types.index(sample_type)]
            
            energy = pd.to_numeric(df[energy_col].iloc[2:], errors='coerce')
            counts = pd.to_numeric(df[sample_col].iloc[2:], errors='coerce')
            
            if normalize:
                counts = counts / counts.max()
            
            trace = go.Scatter(x=energy, y=counts, mode='lines', name=f'{detector_name} {sample_type} Counts')
            traces.append(trace)
            
            # Find prominent peaks
            peaks_energies, peaks_indices = find_prominent_peaks(counts, energy)
            matched_isotopes = match_peaks_to_isotopes(peaks_energies, isotope_reference)
            
            for peak_energy, isotopes in matched_isotopes.items():
                print(f"Detector: {detector_name}, Sample: {sample_type}, Peak Energy: {peak_energy:.2f} keV, Likely Isotopes: {', '.join(isotopes)}")
                isotope_results.append({
                    'Detector': detector_name,
                    'Sample': sample_type,
                    'Peak Energy (keV)': peak_energy,
                    'Likely Isotopes': ', '.join(isotopes)
                })
            
            stats_df = calculate_statistics(counts, count_time)
       
        layout = go.Layout(
            title=f'{sample_type} Counts vs. Energy for All Detectors' + (' (Normalized) Without Background' if normalize else ''),
            xaxis=dict(title='Energy'),
            yaxis=dict(title='Counts'),
            showlegend=True,
            legend=dict(orientation='h', yanchor='bottom', y=-0.1, xanchor='center', x=0.5)
        )
        fig = go.Figure(data=traces, layout=layout)
        
        pio.show(fig)

compiled_statistics_df = compile_statistics_for_target_energy(df, target_energy, count_time)
comparison_df = compare_detectors(df, target_energy, count_time)
print(comparison_df)
parse_and_plot_all(df, normalize=True)

# Save identified isotopes to a CSV file
isotope_df = pd.DataFrame(isotope_results)
isotope_df.to_csv("identified_isotopes.csv", index=False)