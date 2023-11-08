import joblib

def get_wavelength(predicted_array,vacuum_tddft_peak_scaler):
    peakwavs_max_tddft_pred = 1240/vacuum_tddft_peak_scaler.inverse_transform(predicted_array)
    return peakwavs_max_tddft_pred
