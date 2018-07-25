# pbpolcoadd

Copied from [s4cmb](https://github.com/JulienPeloton/s4cmb) by Julien Peloton, based on AnalysisBackend binned map-maker.

## Install

Run
```
pip install . --user
```
in main directory.

## Usage

```
from pbpolcoadd import polcoadd

# T+P binned map-making
polcoadd.tod2map_pol(d0, d4r, d4i, w0, w4, nhit, waferi1d, waferpa, waferts, weight4, weight0, nt, wafermask_pixel, nch, w0.shape[0])

# Only T map-making
polcoadd.tod2map_temp(d0, w0, nhit, waferi1d, waferts, weight0, nt, wafermask_pixel, nch, w0.shape[0])

# Simple map-making (e.g. for ground template)
polcoadd.tod2map_simple(signal_map=signal_map, hits_map=hits_map, pointing=pointing, array=array.flatten(), npix=nch, nt=nt, mask=mask.flatten(), nskypix=npix)
```
