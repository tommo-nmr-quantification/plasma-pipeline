"""Reference-guided peak-height extraction for N,N-dimethylglycine."""
import bisect
import csv
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from . import core

def roi_columns(axis, roi):
    """Input TSV columns with ID at zero; stop is exclusive."""
    high, low = roi
    if not axis[0] >= high > low >= axis[-1]:
        raise ValueError("TCN42 ROI must be covered by the descending axis")
    start = len(axis) - bisect.bisect_right(axis[::-1], high)
    stop = len(axis) - bisect.bisect_left(axis[::-1], low) + 2
    if start < 1 or stop > len(axis) + 1:
        raise ValueError("TCN42 ROI endpoints fall outside spectral columns")
    return start, stop

def four_peaks(indices):
    for comb in combinations(indices, 4):
        if all(abs(d - target) <= tol for d, target, tol in
               zip(np.diff(comb), (27, 76, 27), (2, 7, 2))):
            return list(comb)  # Use the first matching four-peak pattern.
    candidates = []
    for distances, tolerances, missing in [((27, 103), (2, 9), 3),
                                            ((27, 76), (2, 7), 4),
                                            ((76, 27), (7, 2), 1),
                                            ((103, 27), (9, 2), 2)]:
        for a, b, c in combinations(indices, 3):
            if all(abs(d - target) <= tol for d, target, tol in
                   zip((b-a, c-b), distances, tolerances)):
                candidates.append({3: [a, b, c-27, c], 4: [a, b, c, c+27],
                                   1: [a-27, a, b, c], 2: [a, a+27, b, c]}[missing])
    return min(candidates, key=lambda x: sum((x[i+1]-x[i]-d)**2
               for i, d in enumerate((27, 76, 27)))) if candidates else None

def gaussian(x, baseline, amplitude, center, width):
    return baseline + amplitude * np.exp(-((x-center)**2) / (2*width**2))

def gaussian_fit(section, guess):
    x = np.array(section.index)
    y = section.where((section.index >= guess[2]-7) &
                      (section.index <= guess[2]+7), guess[0])
    fitted, _ = curve_fit(gaussian, x, y, guess)
    return [guess[0], *fitted[1:4]]

def double_gaussian(x, baseline, height1, center1, width1, height2, center2, width2):
    # Each Gaussian term includes the baseline parameter.
    return (gaussian(x, baseline, height1, center1, width1) +
            gaussian(x, baseline, height2, center2, width2))

def double_gaussian_fit(section, guess):
    split = int(guess[2] - section.index[0])
    left, right = section.iloc[:split], section.iloc[split:]
    li = np.where(left < guess[0])[0][-1] if np.any(left < guess[0]) else np.argmin(left)
    ri = np.where(right < guess[0])[0][0] if np.any(right < guess[0]) else np.argmin(right)
    y = section.where((section.index >= left.index[li]) &
                      (section.index <= right.index[ri]), guess[0])
    fitted, _ = curve_fit(double_gaussian, np.array(section.index), y, guess)
    return [guess[0], *fitted[1:4]], [guess[0], *fitted[4:]]

def validate_reference(reference, start, stop):
    if isinstance(reference, bool) or not isinstance(reference, (int, np.integer)):
        raise ValueError("reference_pos must be an integer input-TSV column (ID=0)")
    if reference-21 < start or reference+21 >= stop:
        raise ValueError("reference_pos +/- 21 must lie within the extracted ROI")

def extract_height(section, reference):
    """One float64 spectrum ROI, indexed by consecutive input TSV columns (ID=0)."""
    start, stop = int(section.index[0]), int(section.index[-1])+1
    validate_reference(reference, start, stop)
    if list(section.index) != list(range(start, stop)) or not np.isfinite(section.values).all():
        raise ValueError("Expected finite values on consecutive input TSV columns")
    peaks, _ = find_peaks(section, height=section.quantile(0.75), distance=10)
    near = section.loc[reference-21:reference+21]
    sub, _ = find_peaks(near, height=near.max()*0.35, distance=2)
    indices = np.union1d(start+peaks, reference-21+sub)
    quartet = four_peaks(indices)
    result = dict(height=None, mode="unresolved", status="no_quartet",
                  target_column=None, quartet=[], baseline=None, fallback_reason="")
    if quartet is None:
        return result
    if any(k < start or k >= stop for k in quartet):
        raise ValueError("Reconstructed quartet falls outside the ROI")
    target = indices[abs(indices-reference).argmin()]
    baseline = section.quantile(0.3)
    result.update(target_column=int(target), quartet=[int(k) for k in quartet],
                  baseline=float(baseline))
    overlap = 2 if abs(quartet[2]-target) <= 10 else 3 if abs(quartet[3]-target) <= 10 else None
    if overlap is None:
        guess = [baseline, section.loc[target], target, 2.7]
        height = gaussian_fit(section.loc[target-18:target+18], guess)[1]
        mode = "isolated"
    else:
        neighbor = quartet[overlap]
        mode = "overlap%d_subtraction" % (overlap+1)
        if abs(neighbor-target) > 5:
            try:
                guess = [baseline, section.loc[neighbor], neighbor, 2.7,
                         section.loc[target], target, 2.7]
                a, b = double_gaussian_fit(section.loc[min(neighbor,target)-18:
                                                       max(neighbor,target)+18], guess)
                dmg, other = (b, a) if abs(neighbor-a[2]) < abs(neighbor-b[2]) else (a, b)
                if (abs(1-other[1]/(section.loc[neighbor]-baseline)) > 0.5 or
                        min(other[1], dmg[1]) < baseline):
                    raise ValueError("double Gaussian quality rule failed")
                result.update(height=float(dmg[1]), mode="overlap%d_partial" % (overlap+1), status="ok")
                return result
            except (RuntimeError, ValueError) as exc:
                result["fallback_reason"] = type(exc).__name__ + ": " + str(exc)
        else:
            result["fallback_reason"] = "peak separation <= 5 columns"
        other_column = quartet[3 if overlap == 2 else 2]
        guess = [baseline, section.loc[other_column], other_column, 2.7]
        fitted = gaussian_fit(section.loc[other_column-18:other_column+18], guess)
        fitted[2] += -27 if overlap == 2 else 27
        height = section.loc[target] - gaussian(target, *fitted)
    result.update(height=float(height), mode=mode, status="ok")
    return result

def spectra_roi(path, axis, fmt, start, stop):
    """Check the complete row width, then use the pandas C float parser.

    A Python float()/round-trip parser can change the last bits of the feature
    and hence the nonlinear fit. Keep pandas' default C-parser conversion.
    Only the ROI is held in memory; input validation does not convert BINs.
    """
    skip = {'raw8':8,'header':1,'none':0}[fmt]
    count = 0
    with Path(path).open(newline='') as stream:
        reader = csv.reader(stream,delimiter='\t')
        for i in range(skip):
            row = next(reader,None)
            if row is None: raise ValueError('BIN is shorter than the specified preamble')
            if fmt=='header':
                if core.axis_hash(tuple(float(x) for x in row[1:])) != core.axis_hash(axis):
                    raise ValueError('BIN header and model axis differ')
        for row in reader:
            if len(row)!=len(axis)+1 or not row[0]:
                raise ValueError('BIN width mismatch or empty sample ID')
            count+=1
    if not count: raise ValueError('No spectra')
    for chunk in pd.read_csv(path,sep='\t',header=None,skiprows=skip,
                             usecols=[0]+list(range(start,stop)),dtype={0:str},
                             chunksize=64,engine='c',float_precision=None,keep_default_na=False):
        for _,row in chunk.iterrows():
            yield str(row.iloc[0]),row.iloc[1:].astype('float64')
