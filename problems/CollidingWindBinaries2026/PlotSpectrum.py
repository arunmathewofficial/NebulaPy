import os
import numpy as np
import matplotlib.pyplot as plt

SpectrumFile = "/Users/tony/Desktop/Bowshock-Xray/Post-Processing/XraySpectrumTest/Ostar_mhd-nemo-dep_d2n0128l3_Optical_1.613175e+02.txt"
CatalogFile  = "/Users/tony/Desktop/Bowshock-Xray/Post-Processing/XraySpectrumTest/Ostar_mhd-nemo-dep_d2n0128l3_LineCatalog_Optical.txt"

OutDir = "/Users/tony/Desktop/Bowshock-Xray/Post-Processing/XraySpectrumTest"
Filename = "OpticalSpectrum_with_lines"

match_tolerance = 2.0
max_labels = 60
min_peak_separation = 2
peak_prominence_factor = 1.01


def read_line_catalog(catalog_file):
    catalog = []

    with open(catalog_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#") or line.startswith("Spectral line") or line.startswith("-"):
                continue

            parts = line.split()

            try:
                catalog.append({
                    "ion": f"{parts[0]} {parts[1]}",
                    "wavelength": float(parts[2]),
                    "energy": float(parts[3]),
                    "avalue": float(parts[4]),
                    "text": line
                })
            except Exception:
                continue

    return catalog


def find_spectrum_peaks(wavelength, spectrum, min_separation=2, prominence_factor=1.01):
    peaks = []

    for i in range(1, len(spectrum) - 1):

        if spectrum[i] > spectrum[i - 1] and spectrum[i] > spectrum[i + 1]:

            i1 = max(0, i - min_separation)
            i2 = min(len(spectrum), i + min_separation + 1)

            local_region = np.concatenate(
                (spectrum[i1:i], spectrum[i + 1:i2])
            )

            local_continuum = np.median(local_region)

            if local_continuum > 0.0 and spectrum[i] / local_continuum >= prominence_factor:
                peaks.append(i)

    return np.array(peaks, dtype=int)


def match_peaks_to_catalog(wavelength, spectrum, peak_indices, catalog, tolerance):
    matches = []

    catalog_wvl = np.array([line["wavelength"] for line in catalog])

    for idx in peak_indices:
        peak_wvl = wavelength[idx]
        peak_lum = spectrum[idx]

        nearest_idx = np.argmin(np.abs(catalog_wvl - peak_wvl))
        best_line = catalog[nearest_idx]
        best_sep = abs(best_line["wavelength"] - peak_wvl)

        if best_sep <= tolerance:
            matches.append({
                "peak_wavelength": peak_wvl,
                "peak_luminosity": peak_lum,
                "catalog_ion": best_line["ion"],
                "catalog_wavelength": best_line["wavelength"],
                "catalog_energy": best_line["energy"],
                "catalog_avalue": best_line["avalue"],
                "separation": best_sep,
                "catalog_entry": best_line["text"]
            })

    return matches


wavelength, spectrum = np.loadtxt(SpectrumFile, unpack=True)
catalog = read_line_catalog(CatalogFile)

peak_indices = find_spectrum_peaks(
    wavelength,
    spectrum,
    min_separation=min_peak_separation,
    prominence_factor=peak_prominence_factor
)

matches = match_peaks_to_catalog(
    wavelength,
    spectrum,
    peak_indices,
    catalog,
    tolerance=match_tolerance
)

matches = sorted(matches, key=lambda x: x["peak_luminosity"], reverse=True)

matched_file = os.path.join(OutDir, f"{Filename}_matched_lines.txt")

with open(matched_file, "w") as f:
    f.write(
        "# Peak_wavelength[A] "
        "Peak_Lambda[erg s^-1 A^-1] "
        "Ion "
        "Catalog_wavelength[A] "
        "Energy[keV] "
        "A-value "
        "Separation[A] "
        "Catalog_entry\n"
    )

    for m in matches:
        f.write(
            f"{m['peak_wavelength']:.8e} "
            f"{m['peak_luminosity']:.8e} "
            f"{m['catalog_ion']} "
            f"{m['catalog_wavelength']:.8e} "
            f"{m['catalog_energy']:.8e} "
            f"{m['catalog_avalue']:.8e} "
            f"{m['separation']:.8e} "
            f"{m['catalog_entry']}\n"
        )


fig, ax = plt.subplots(figsize=(10, 5))

ax.set_xlabel(r"Wavelength [$\AA$]", fontsize=12)
ax.set_ylabel(r"$L_\lambda$ [erg s$^{-1}$ $\AA^{-1}$]", fontsize=12)

ax.plot(
    wavelength,
    spectrum,
    color="green",
    linewidth=1.4,
    label="Spectrum"
)

ax.set_yscale("log")

positive_spectrum = spectrum[spectrum > 0.0]
ymin = np.min(positive_spectrum)
ymax = np.max(positive_spectrum)

ax.set_ylim(ymin * 0.8, ymax * 20.0)

for i, m in enumerate(matches[:max_labels]):

    peak_x = m["peak_wavelength"]
    peak_y = m["peak_luminosity"]

    line_label = f"{m['catalog_ion']} {m['catalog_wavelength']:.2f}"

    yoffset = 1.15 + 0.25 * (i % 4)

    ax.text(
        peak_x,
        peak_y * yoffset,
        line_label,
        rotation=90,
        fontsize=6.5,
        va="bottom",
        ha="center",
        color="black",
        clip_on=False
    )

ax.minorticks_on()

ax.tick_params(
    axis='both',
    which='major',
    direction='in',
    top=True,
    right=True,
    length=6,
    width=1.2,
    labelsize=11
)

ax.tick_params(
    axis='both',
    which='minor',
    direction='in',
    top=True,
    right=True,
    length=3,
    width=1.0
)

for spine in ax.spines.values():
    spine.set_linewidth(1.2)

ax.legend(loc='best', frameon=False, fontsize=10)

fig.tight_layout()

outfile = os.path.join(OutDir, f"{Filename}.png")

fig.savefig(outfile, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f" Read spectrum from {SpectrumFile}")
print(f" Read catalog from {CatalogFile}")
print(f" Number of catalog lines: {len(catalog)}")
print(f" Found {len(peak_indices)} spectral peaks")
print(f" Matched {len(matches)} peaks with catalog lines")
print(f" Saved matched line list to {matched_file}")
print(f" Saved spectrum plot to {outfile}")