import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import numpy as np

# Load the Excel file and specify the sheet name
file_path = r"C:\Users\17244\OneDrive\Grad School\NUCE 597\Exp2\Data.xlsx"
df = pd.read_excel(file_path, sheet_name='Sheet1')

count_time = 300  # Count time in seconds
target_energy = 1332.50  # Target energy in keV

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