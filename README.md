# ALMA Obs Calc

`ALMA_obs_calc.py` is a small Python utility for estimating line fluxes and observing properties for ALMA targets. It converts galaxy properties such as star formation rate, redshift, line width, and source classification into predicted quantities for CO, [CII] 158 um, and [OIII] 88 um.

The code is based on a script from Gareth Jones and was modified by Jan Scholtz.

## What it does

The main `ALMA_calc` class can:

- Estimate molecular gas mass from star formation rate for main-sequence or starburst galaxies.
- Convert molecular gas mass into CO line luminosity and flux density.
- Estimate [CII] 158 um and [OIII] 88 um flux densities from star formation rate.
- Compute a simple Gaussian line profile and sample it into spectral channels.
- Report representative frequency, channel width, peak amplitude, and spectral window placement.

## Requirements

- Python 3
- `numpy`
- `astropy`

Install the dependencies with:

```bash
pip install numpy astropy
```

## Input data

The calculator expects a galaxy dictionary with fields such as:

- `NAME` # Name of the object 
- `SFR` # SFR - used for calculation of [CII], [OIII] or CO (for CO only if MH2 is not provided)
- `ZRED` # Redshift
- `MSorSB`
- `alpha`
- `FWHM`
- `numchan`
- `flipper`
- `ZP`
- `Cosmo`
- `PSF_n`
- `mu`
- `CO_ladder`

Optional fields include:

- `MH2`

## Example usage

```python
from astropy.cosmology import FlatLambdaCDM
from ALMA_obs_calc import ALMA_calc

gal = {
    'NAME': 'Example Galaxy',
    'SFR': 50,
    'ZRED': 2.5,
    'MSorSB': 'MS',
    'alpha': 4.3,
    'FWHM': 300,
    'numchan': 9,
    'flipper': [False, False],
    'ZP': 1.0,
    'Cosmo': FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725),
    'PSF_n': 1,
    'mu': 1,
    'CO_ladder': 'SF',
}

calc = ALMA_calc(gal)

calc.CO_calc([1, 2, 3])
calc.CII_calc()
calc.OIII_calc()
```

## Notes

- The script prints results directly to stdout and does not currently write files or expose a command-line interface.
- Several calculations depend on empirical relations from the literature, so the outputs are intended as estimates rather than exact predictions.
- The CO ladder can be switched between star-forming and AGN-like assumptions with the `CO_ladder` field.

## Output

Each calculation prints values such as:

- Line luminosity or integrated flux
- Representative observing frequency
- Gaussian peak amplitude
- Per-channel flux estimates
- Suggested spectral window centers

## Reference

The script includes relations drawn from the following sources mentioned in the code comments:

- Sargent et al. 2014
- Narayanan et al. 2012
- Carilli & Walter 2013
- De Looze et al. 2014
- Herrera-Camus et al. 2021
