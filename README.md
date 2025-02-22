- [Astrometry](#astrometry)
- [Get started](#get-started)
- [Examples](#examples)
  - [Use the solver without `with`](#use-the-solver-without-with)
  - [Provide size and position hints](#provide-size-and-position-hints)
  - [Print progress information (download and solve)](#print-progress-information-download-and-solve)
  - [Print field stars metadata](#print-field-stars-metadata)
  - [Calculate field stars pixel positions with astropy](#calculate-field-stars-pixel-positions-with-astropy)
  - [Print series description and size (without downloading them)](#print-series-description-and-size-without-downloading-them)
  - [Disable tune-up and distortion](#disable-tune-up-and-distortion)
  - [Stop the solver early using the log-odds callback](#stop-the-solver-early-using-the-log-odds-callback)
    - [Return after the first match](#return-after-the-first-match)
    - [Return early if the best log-odds are larger than 100.0](#return-early-if-the-best-log-odds-are-larger-than-1000)
    - [Return early if there are at least ten matches](#return-early-if-there-are-at-least-ten-matches)
    - [Return early if the three best matches are similar](#return-early-if-the-three-best-matches-are-similar)
- [Choosing series](#choosing-series)
- [Documentation](#documentation)
  - [Solver](#solver)
  - [SizeHint](#sizehint)
  - [PositionHint](#positionhint)
  - [Action](#action)
  - [Parity](#parity)
  - [SolutionParameters](#solutionparameters)
  - [Solution](#solution)
  - [Match](#match)
  - [Star](#star)
  - [Series](#series)
  - [batches\_generator](#batches_generator)
  - [SupportsFloatMapping](#supportsfloatmapping)
- [Contribute](#contribute)
  - [Clone the repository](#clone-the-repository)
  - [Format and lint](#format-and-lint)
  - [Build a local version](#build-a-local-version)
  - [Publish](#publish)
  - [MSVC compatibility (work in progress)](#msvc-compatibility-work-in-progress)

# Astrometry

Astrometry turns a list of star positions into a pixel-to-sky transformation (WCS) by calling C functions from the Astrometry.net library (https://astrometry.net).

Astrometry.net star index files ("series") are automatically downloaded when required.

This package is useful for solving plates from a Python script, comparing star extraction methods, or hosting a simple local version of Astrometry.net with minimal dependencies. See https://github.com/dam90/astrometry for a more complete self-hosting solution.

Unlike Astrometry.net, Astrometry does not include FITS parsing or image pre-processing algorithms. Stars must be provided as a list of pixel positions.

This library works on Linux and macOS but it requires WSL on Windows.

We are not the authors of the Astrometry.net library. You should cite works from https://astrometry.net/biblio.html if you use the Astrometry.net algorithm via this package.

# Get started

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install astrometry
```

```py
import astrometry

with astrometry.Solver(
    astrometry.series_5200.index_files(
        cache_directory="astrometry_cache",
        scales={6},
    )
) as solver:

    stars = [
        [388.9140568247906, 656.5003281719216],
        [732.9210858972549, 473.66395545775106],
        [401.03459504299843, 253.788113189415],
        [312.6591868096163, 624.7527729425295],
        [694.6844564647456, 606.8371776658344],
        [741.7233477959561, 344.41284826261443],
        [867.3574610200455, 672.014835980283],
        [1063.546651153479, 593.7844603550848],
        [286.69070190952704, 422.170016812049],
        [401.12779619355155, 16.13543616977013],
        [205.12103484692776, 698.1847350789413],
        [202.88444768690894, 111.24830187635557],
        [339.1627757703069, 86.60739435924549],
    ]

    solution = solver.solve(
        stars=stars,
        size_hint=None,
        position_hint=None,
        solution_parameters=astrometry.SolutionParameters(),
    )

    if solution.has_match():
        print(f"{solution.best_match().center_ra_deg=}")
        print(f"{solution.best_match().center_dec_deg=}")
        print(f"{solution.best_match().scale_arcsec_per_pixel=}")
```

`solve` is thread-safe. It can be called any number of times from the same `Solver` object.

# Examples

## Use the solver without `with`

```py
import astrometry

solver = astrometry.Solver(...)
solver.solve(...) # solve can be called multiple times
solver.close() # close the index files used by the solver
```

## Provide size and position hints

```py
import astrometry

with astrometry.Solver(...) as solver:
    solution = solver.solve(
        stars=...,
        size_hint=astrometry.SizeHint(
            lower_arcsec_per_pixel=1.0,
            upper_arcsec_per_pixel=2.0,
        ),
        position_hint=astrometry.PositionHint(
            ra_deg=65.7,
            dec_deg=36.2,
            radius_deg=1.0,
        ),
        solution_parameters=...
    )
```

## Print progress information (download and solve)

```py
import astrometry
import logging

logging.getLogger().setLevel(logging.INFO)

with astrometry.Solver(...) as solver:
    solution = solver.solve(...)
```

## Print field stars metadata

Astrometry extracts metadata from the star index ("series"). See [Choosing series](#choosing-series) for a description of the available data.

```py
import astrometry

with astrometry.Solver(...) as solver:
    solution = solver.solve(...)

    if solution.has_match():
        for star in solution.best_match().stars:
            print(f"{star.ra_deg}º, {star.dec_deg}º:", star.metadata)
```

## Calculate field stars pixel positions with astropy

```py
import astrometry

with astrometry.Solver(...) as solver:
    solution = solver.solve(...)

    if solution.has_match():
        wcs = solution.best_match().astropy_wcs()
        pixels = wcs.all_world2pix(
            [[star.ra_deg, star.dec_deg] for star in solution.best_match().stars],
            0,
        )
        # pixels is a len(solution.best_match().stars) x 2 numpy array of float values
```

`astropy.wcs.WCS` provides many more functions to probe the transformation properties and convert from and to pixel coordinates. See https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html for details. Astropy (https://pypi.org/project/astropy/) must be installed to use this method.

## Print series description and size (without downloading them)

```py
import astrometry

print(astrometry.series_5200_heavy.description)
print(astrometry.series_5200_heavy.size_as_string({2, 3, 4}))
```

See [Choosing Series](#choosing-series) for a list of available series.

## Disable tune-up and distortion

```py
import astrometry

with astrometry.Solver(...) as solver:
    solution = solver.solve(
        stars=...,
        size_hint=...,
        position_hint=...,
        solution_parameters=astrometry.SolutionParameters(
            sip_order=0,
            tune_up_logodds_threshold=None,
        ),
    )
```

## Stop the solver early using the log-odds callback

### Return after the first match

```py
import astrometry

with astrometry.Solver(...) as solver:
    solution = solver.solve(
        stars=...,
        size_hint=...,
        position_hint=...,
        solution_parameters=astrometry.SolutionParameters(
            logodds_callback=lambda logodds_list: astrometry.Action.STOP,
        ),
    )
```

### Return early if the best log-odds are larger than 100.0

```py
import astrometry

with astrometry.Solver(...) as solver:
    solution = solver.solve(
        stars=...,
        size_hint=...,
        position_hint=...,
        solution_parameters=astrometry.SolutionParameters(
            logodds_callback=lambda logodds_list: (
                astrometry.Action.STOP
                if logodds_list[0] > 100.0
                else astrometry.Action.CONTINUE
            ),
        ),
    )
```

### Return early if there are at least ten matches

```py
import astrometry

with astrometry.Solver(...) as solver:
    solution = solver.solve(
        stars=...,
        size_hint=...,
        position_hint=...,
        solution_parameters=astrometry.SolutionParameters(
            logodds_callback=lambda logodds_list: (
                astrometry.Action.STOP
                if len(logodds_list) >= 10.0
                else astrometry.Action.CONTINUE
            ),
        ),
    )
```

### Return early if the three best matches are similar

```py
import astrometry

def logodds_callback(logodds_list: list[float]) -> astrometry.Action:
    if len(logodds_list) < 3:
        return astrometry.Action.CONTINUE
    if logodds_list[1] > logodds_list[0] - 10 and logodds_list[2] > logodds_list[0] - 10:
        return astrometry.Action.STOP
    return astrometry.Action.CONTINUE


with astrometry.Solver(...) as solver:
    solution = solver.solve(
        stars=...,
        size_hint=...,
        position_hint=...,
        solution_parameters=astrometry.SolutionParameters(
            logodds_callback=logodds_callback,
        ),
    )
```

# Choosing series

This library downloads series from http://data.astrometry.net. A solver can be instantiated with multiple series and scales as follows:

```py
import astrometry

with astrometry.Solver(
    astrometry.series_5200.index_files(
        cache_directory="astrometry_cache",
        scales={4, 5, 6},
    )
    + astrometry.series_4200.index_files(
        cache_directory="astrometry_cache",
        scales={6, 7, 12},
    )
) as solver:
    ...
```

Astrometry.net gives the following recommendations to choose a scale:

> Each index file contains a large number of “skymarks” (landmarks for the sky) that allow our solver to identify your images. The skymarks contained in each index file have sizes (diameters) within a narrow range. You probably want to download index files whose quads are, say, 10% to 100% of the sizes of the images you want to solve.
>
> For example, let’s say you have some 1-degree square images. You should grab index files that contain skymarks of size 0.1 to 1 degree, or 6 to 60 arcminutes. Referring to the table below, you should [try index files with scales 3 to 9]. You might find that the same number of fields solve, and faster, using just one or two of the index files in the middle of that range - in our example you might try [5, 6 and 7].
>
> -- _http://astrometry.net/doc/readme.html_

| Scale                     | 0               | 1               | 2               | 3               | 4            | 5             | 6             | 7             | 8             | 9             | 10            | 11             | 12              | 13              | 14              | 15              | 16              | 17               | 18                | 19                |
| ------------------------- | --------------- | --------------- | --------------- | --------------- | ------------ | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------- | --------------- | --------------- | --------------- | --------------- | --------------- | ---------------- | ----------------- | ----------------- |
| Skymark diameter (arcmin) | [2.0,&nbsp;2.8] | [2.8,&nbsp;4.0] | [4.0,&nbsp;5.6] | [5.6,&nbsp;8.0] | [8,&nbsp;11] | [11,&nbsp;16] | [16,&nbsp;22] | [22,&nbsp;30] | [30,&nbsp;42] | [42,&nbsp;60] | [60,&nbsp;85] | [85,&nbsp;120] | [120,&nbsp;170] | [170,&nbsp;240] | [240,&nbsp;340] | [340,&nbsp;480] | [480,&nbsp;680] | [680,&nbsp;1000] | [1000,&nbsp;1400] | [1400,&nbsp;2000] |

The table below lists series sizes and properties (we copied the descriptions from http://data.astrometry.net). You can access a series' object with `astrometry.series_{name}` (for example `astrometry.series_4200`).

| Name       | Total&nbsp;size | Scales       | Description                                                                                                                                                                                                                                          | Metadata                                                                                                                                                                                                                                                                                                                                                                                  |
| ---------- | --------------- | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 4100       | 0.36&nbsp;GB    | [7,&nbsp;19] | built from the Tycho-2 catalog, good for images wider than 1 degree, recommended                                                                                                                                                                     | MAG_BT:&nbsp;float<br />MAG_VT:&nbsp;float<br />MAG_HP:&nbsp;float<br />MAG:&nbsp;float                                                                                                                                                                                                                                                                                                   |
| 4200       | 33.78&nbsp;GB   | [0,&nbsp;19] | built from the near-infared 2MASS survey, runs out of stars at the low end, most users will probably prefer 4100 or 5200                                                                                                                             | j_mag:&nbsp;float                                                                                                                                                                                                                                                                                                                                                                         |
| 5000       | 76.24&nbsp;GB   | [0,&nbsp;7]  | an older version from Gaia-DR2 but without Tycho-2 stars merged in, our belief is that series_5200 will work better than this one                                                                                                                    | source_id:&nbsp;int<br />phot_g_mean_mag:&nbsp;float<br />phot_bp_mean_mag:&nbsp;float<br />phot_rp_mean_mag:&nbsp;float<br />parallax:&nbsp;float<br />parallax_error:&nbsp;float<br />pmra:&nbsp;float<br />pmra_error:&nbsp;float<br />pmdec:&nbsp;float<br />pmdec_error:&nbsp;float<br />ra:&nbsp;float<br />dec:&nbsp;float<br />ref_epoch:&nbsp;float                              |
| 5200       | 36.14&nbsp;GB   | [0,&nbsp;6]  | LIGHT version built from Tycho-2 + Gaia-DR2, good for images narrower than 1 degree, combine with 4100-series for broader scale coverage, the LIGHT version contains smaller files with no additional Gaia-DR2 information tagged along, recommended | -                                                                                                                                                                                                                                                                                                                                                                                         |
| 5200_heavy | 79.67&nbsp;GB   | [0,&nbsp;6]  | HEAVY version same as 5200, but with additional Gaia-DR2 information (magnitude in G, BP, RP, proper motions and parallaxes), handy if you want that extra Gaia information for matched stars                                                        | ra:&nbsp;float<br />dec:&nbsp;float<br />mag:&nbsp;float<br />ref_cat:&nbsp;str<br />ref_id:&nbsp;int<br />pmra:&nbsp;float<br />pmdec:&nbsp;float<br />parallax:&nbsp;float<br />ra_ivar:&nbsp;float<br />dec_ivar:&nbsp;float<br />pmra_ivar:&nbsp;float<br />pmdec_ivar:&nbsp;float<br />parallax_ivar:&nbsp;float<br />phot_bp_mean_mag:&nbsp;float<br />phot_rp_mean_mag:&nbsp;float |
| 6000       | 1.20&nbsp;GB    | [4,&nbsp;6]  | very specialized, uses GALEX Near-UV measurements, and only a narrow range of scales                                                                                                                                                                 | fuv_mag:&nbsp;float<br />nuv_mag:&nbsp;float                                                                                                                                                                                                                                                                                                                                              |
| 6100       | 1.58&nbsp;GB    | [4,&nbsp;6]  | very specialized, uses GALEX Far-UV measurements, and only a narrow range of scales                                                                                                                                                                  | fuv_mag:&nbsp;float<br />nuv_mag:&nbsp;float                                                                                                                                                                                                                                                                                                                                              |

The table below indicates the total file size for each scale (most series have multiple index files per scale).

| Name       | 0             | 1             | 2             | 3            | 4              | 5              | 6              | 7              | 8             | 9             | 10            | 11            | 12           | 13           | 14           | 15             | 16             | 17             | 18             | 19             |
| ---------- | ------------- | ------------- | ------------- | ------------ | -------------- | -------------- | -------------- | -------------- | ------------- | ------------- | ------------- | ------------- | ------------ | ------------ | ------------ | -------------- | -------------- | -------------- | -------------- | -------------- |
| 4100       | -             | -             | -             | -            | -              | -              | -              | 165.00&nbsp;MB | 94.55&nbsp;MB | 49.77&nbsp;MB | 24.87&nbsp;MB | 10.21&nbsp;MB | 5.30&nbsp;MB | 2.73&nbsp;MB | 1.38&nbsp;MB | 740.16&nbsp;kB | 408.96&nbsp;kB | 247.68&nbsp;kB | 187.20&nbsp;kB | 144.00&nbsp;kB |
| 4200       | 14.22&nbsp;GB | 9.25&nbsp;GB  | 5.06&nbsp;GB  | 2.63&nbsp;GB | 1.31&nbsp;GB   | 659.09&nbsp;MB | 328.25&nbsp;MB | 165.44&nbsp;MB | 81.84&nbsp;MB | 41.18&nbsp;MB | 20.52&nbsp;MB | 8.02&nbsp;MB  | 4.17&nbsp;MB | 2.16&nbsp;MB | 1.10&nbsp;MB | 596.16&nbsp;kB | 339.84&nbsp;kB | 213.12&nbsp;kB | 164.16&nbsp;kB | 132.48&nbsp;kB |
| 5000       | 34.79&nbsp;GB | 20.19&nbsp;GB | 10.74&nbsp;GB | 5.44&nbsp;GB | 2.71&nbsp;GB   | 1.36&nbsp;GB   | 676.79&nbsp;MB | 340.73&nbsp;MB | -             | -             | -             | -             | -            | -            | -            | -              | -              | -              | -              | -              |
| 5200       | 17.20&nbsp;GB | 9.49&nbsp;GB  | 4.86&nbsp;GB  | 2.45&nbsp;GB | 1.22&nbsp;GB   | 614.89&nbsp;MB | 307.72&nbsp;MB | -              | -             | -             | -             | -             | -            | -            | -            | -              | -              | -              | -              | -              |
| 5200_heavy | 36.46&nbsp;GB | 21.20&nbsp;GB | 11.29&nbsp;GB | 5.72&nbsp;GB | 2.85&nbsp;GB   | 1.43&nbsp;GB   | 714.56&nbsp;MB | -              | -             | -             | -             | -             | -            | -            | -            | -              | -              | -              | -              | -              |
| 6000       | -             | -             | -             | -            | 892.55&nbsp;MB | 457.66&nbsp;MB | 233.23&nbsp;MB | -              | -             | -             | -             | -             | -            | -            | -            | -              | -              | -              | -              | -              |
| 6100       | -             | -             | -             | -            | 599.33&nbsp;MB | 384.09&nbsp;MB | 214.79&nbsp;MB | -              | -             | -             | -             | -             | -            | -            | -            | -              | -              | -              | -              | -              |

# Documentation

## Solver

```py
class Solver:
    def __init__(self, index_files: list[pathlib.Path]): ...

    def solve(
        self,
        stars: typing.Iterable[SupportsFloatMapping],
        size_hint: typing.Optional[SizeHint],
        position_hint: typing.Optional[PositionHint],
        solution_parameters: SolutionParameters,
    ) -> Solution: ...

    def __enter__(self) -> Solver: ...

    def __exit__(
        self,
        exception_type: typing.Optional[typing.Type[BaseException]],
        value: typing.Optional[BaseException],
        traceback: typing.Optional[types.TracebackType],
    ) -> bool: ...
```

`solve` is thread-safe and can be called any number of times.

-   `index_files`: List of index files to use for solving. The list need not come from a `Series` object. Series subsets and combinations are possible as well.
-   `star`: iterator over a list of pixel coordinates for the input stars.
-   `size_hint`: Optional angular pixel size range ([SizeHint](#sizehint)). Significantly speeds up `solve` when provided. If `size_hint` is `None`, the range `[0.1, 1000.0]` is used. This default range can be changed by setting `astrometry.DEFAULT_LOWER_ARCSEC_PER_PIXEL` and `astrometry.DEFAULT_UPPER_ARCSEC_PER_PIXEL` to other values.
-   `position_hint`: Optional field center Ra/Dec coordinates and error radius ([PositionHint](#positionhint)). Significantly speeds up `solve` when provided. If `position_hint` is None, the entire sky is used (`radius_deg = 180.0`).
-   `solution_parameters`: Advanced solver parameters ([SolutionParameters](#solutionparameters))

## SizeHint

```py
@dataclasses.dataclass
class SizeHint:
    lower_arcsec_per_pixel: float
    upper_arcsec_per_pixel: float
```

`lower_arcsec_per_pixel` and `upper_arcsec_per_pixel` must be larger than `0` and `upper_arcsec_per_pixel` must be smaller than or equal to `upper_arcsec_per_pixel`.

## PositionHint

```py
@dataclasses.dataclass
class PositionHint:
    ra_deg: float
    dec_deg: float
    radius_deg: float
```

-   `ra_deg` must be in the range `[0.0, 360.0[`.
-   `dec_deg` must be in the range `[-90.0, 90.0]`.
-   `radius` must be larger than or equal to zero.

All values are in degrees and must use the same frame of reference as the index files. Astrometry.net index files use J2000 FK5 (https://docs.astropy.org/en/stable/api/astropy.coordinates.FK5.html). ICRS and FK5 differ by less than 0.1 arcsec (https://www.iers.org/IERS/EN/Science/ICRS/ICRS.html).

## Action

```py
class Action(enum.Enum):
    STOP = 0
    CONTINUE = 1
```

## Parity

```py
class Parity(enum.IntEnum):
    NORMAL = 0
    FLIP = 1
    BOTH = 2
```

## SolutionParameters

```py
@dataclasses.dataclass
class SolutionParameters:
    solve_id: typing.Optional[str] = None
    uniformize_index: bool = True
    deduplicate: bool = True
    sip_order: int = 3
    sip_inverse_order: int = 0
    distance_from_quad_bonus: bool = True
    positional_noise_pixels: float = 1.0
    distractor_ratio: float = 0.25
    code_tolerance_l2_distance: float = 0.01
    minimum_quad_size_pixels: typing.Optional[float] = None
    minimum_quad_size_fraction: float = 0.1
    maximum_quad_size_pixels: float = 0.0
    maximum_quads: int = 0
    maximum_matches: int = 0
    parity: Parity = Parity.BOTH
    tune_up_logodds_threshold: typing.Optional[float] = 14.0
    output_logodds_threshold: float = 21.0
    slices_generator: typing.Callable[[int], typing.Iterable[tuple[int, int]]] = astrometry.batches_generator(25)
    logodds_callback: typing.Callable[[list[float]], Action] = lambda _: Action.CONTINUE
```

-   `solve_id`: Optional plate identifier used in logging messages. If `solve_id` is `None`, it is automatically assigned a unique integer. The value can be retrieved from the Solution object (`solution.solve_id`).
-   `uniformize_index`: Uniformize field stars at the matched index scale before verifying a match.
-   `deduplicate`: De-duplicate field stars before verifying a match.
-   `sip_order`: Polynomial order of the Simple Imaging Polynomial distortion (see https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf). `0` disables SIP distortion. `tune_up_logodds_threshold` must be `None` if `sip_order` is `0`.
-   `sip_inverse_order`: Polynomial order of the inversee Simple Polynomial distortion. Usually equal to `sip_order`. `0` means "equal to `sip_order`".
-   `distance_from_quad_bonus`: Assume that stars far from the matched quad will have larger positional variance.
-   `positional_noise_pixels`: Expected error on the positions of stars.
-   `distractor_ratio`: Fraction of distractors in the range `]0, 1]`.
-   `code_tolerance_l2_distance`: Code tolerance in 4D codespace L2 distance.
-   `minimum_quad_size_pixels`: Minimum size of field quads to try, `None` calculates the size automatically as `minimum_quad_size_fraction * min(Δx, Δy)`, where `Δx` (resp. `Δy`) is the maximum `x` distance (resp. `y` distance) between stars in the field.
-   `minimum_quad_size_fraction`: Only used if `minimum_quad_size_pixels` is `None` (see above).
-   `maximum_quad_size_pixels`: Maximum size of field quads to try, `0.0` means no limit.
-   `maximum_quads`: Number of field quads to try, `0` means no limit.
-   `maximum_matches`: Number of quad matches to try, `0` means no limit.
-   `parity`: Parity.NORMAL does not flip the axes, Parity.FLIP does, and Parity.BOTH tries flipped and non-flipped axes (at the cost of doubling computations).
-   `tune_up_logodds_threshold`: Matches whose log-odds are larger than this value are tuned-up (SIP distortion estimation) and accepted if their post-tune-up log-odds are larger than `output_logodds_threshold`. `None` disables tune-up and distortion estimation (SIP). The default Astrometry.net value is `math.log(1e6)`.
-   `output_logodds_threshold`: Matches whose log-odds are larger than this value are immediately accepted (added to the solution matches). The default Astrometry.net value is `math.log(1e9)`.
-   `slices_generator`: User-provided function that takes a number of stars as parameter and returns an iterable (such as list) of two-elements tuples representing ranges. The first tuple item (start) is included while the second tuple item (end) is not. The returned ranges can have a variable size and/or overlap. The algorithm compares each range with the star catalogue sequentially. Small ranges significantly speed up the algorithm but increase the odds of missing matches. `astrometry.batches_generator(n)` generates non-overlapping batches of `n` stars.
-   `logodds_callback`: User-provided function that takes a list of matches log-odds as parameter and returns an `astrometry.Action` object. `astrometry.Action.CONTINUE` tells the solver to keep searching for matches whereas `astrometry.Action.STOP` tells the solver to return the current matches immediately. The log-odds list is sorted from highest to lowest value and should not be modified by the callback function.

Accepted matches are always tuned up, even if they hit `tune_up_logodds_threshold` and were already tuned-up. Since log-odds are compared with the thresholds before the tune-up, the final log-odds are often significantly larger than `output_logodds_threshold`. Set `tune_up_logodds_threshold` to a value larger than or equal to `output_logodds_threshold` to disable the first tune-up, and `None` to disable tune-up altogether. Tune-up logic is equivalent to the following Python snippet:

```py
# This (pseudo-code) snippet assumes the following definitions:
# match: candidate match object
# log_odds: current match log-odds
# add_to_solution: appends the match to the solution list
# tune_up: tunes up a match object and returns the new match and the new log-odds
if tune_up_logodds_threshold is None:
    if log_odds >= output_logodds_threshold:
        add_to_solution(match)
else:
    if log_odds >= output_logodds_threshold:
        tuned_up_match, tuned_up_loggods = tune_up(match)
        add_to_solution(tuned_up_match)
    elif log_odds >= tune_up_logodds_threshold:
        tuned_up_match, tuned_up_loggods = tune_up(match)
        if tuned_up_loggods >= output_logodds_threshold:
            tuned_up_twice_match, tuned_up_twice_loggods = tune_up(tuned_up_match)
            add_to_solution(tuned_up_twice_match)
```

Astrometry.net gives the following description of the tune-up algorithm. See `tweak2` in _astrometry.net/solver/tweak2.c_ for the implementation.

> Given an initial WCS solution, compute SIP polynomial distortions using an annealing-like strategy. That is, it finds matches between image and reference catalog by searching within a radius, and that radius is small near a region of confidence, and grows as you move away. That makes it possible to pick up more distant matches, but they are downweighted in the fit. The annealing process reduces the slope of the growth of the matching radius with respect to the distance from the region of confidence.
>
> -- _astrometry.net/include/astrometry/tweak2.h_

## Solution

```py
@dataclasses.dataclass
class Solution:
    solve_id: str
    matches: list[Match]

    def has_match(self) -> bool: ...

    def best_match(self) -> Match: ...

    def to_json(self) -> str: ...

    @classmethod
    def from_json(cls, solution_as_json: str) -> Solution: ...

```

`matches` are sorted in descending log-odds order. `best_match` returns the first match in the list. `to_json` and `from_json` may be used to save and load solutions.

## Match

```py
@dataclasses.dataclass
class Match:
    logodds: float
    center_ra_deg: float
    center_dec_deg: float
    scale_arcsec_per_pixel: float
    index_path: pathlib.Path
    stars: tuple[Star, ...]
    quad_stars: tuple[Star, ...]
    wcs_fields: dict[str, tuple[typing.Any, str]]

    def astropy_wcs(self) -> astropy.wcs.WCS: ...
```

-   `logodds`: Log-odds (https://en.wikipedia.org/wiki/Logit) of the match.
-   `center_ra_deg`: Right ascension of the stars bounding box's center, in degrees and in the frame of reference of the index (J200 FK5 for Astrometry.net series).
-   `center_dec_deg`: Declination of the stars bounding box's center in degrees and in the frame of reference of the index (J200 FK5 for Astrometry.net series).
-   `scale_arcsec_per_pixel`: Pixel scale in arcsec per pixel.
-   `index_path`: File system path of the index file used for this match.
-   `stars`: List of visible index stars. This list is almost certainly going to differ from the input stars list.
-   `quad_stars`: The index stars subset (usually 4 but can be 3 or 5) used in the hash code search step (see https://arxiv.org/pdf/0910.2233.pdf, 2. Methods).
-   `wcs_fields`: WCS fields describing the transformation between pixel coordinates and world coordinates. This dictionary can be passed directly to `astropy.wcs.WCS`.

`astropy_wcs` generates an Astropy WCS object. Astropy (https://pypi.org/project/astropy/) must be installed to use this method. See [Calculate field stars pixel positions with astropy](#calculate-field-stars-pixel-positions-with-astropy) for details.

## Star

```py
@dataclasses.dataclass
class Star:
    ra_deg: float
    dec_deg: float
    metadata: dict[str, typing.Any]
```

`ra_deg` and `dec_deg` are in degrees and use the same frame of reference as the index files. Astrometry.net index files use J2000 FK5 (https://docs.astropy.org/en/stable/api/astropy.coordinates.FK5.html). ICRS and FK5 differ by less than 0.1 arcsec (https://www.iers.org/IERS/EN/Science/ICRS/ICRS.html).

The contents of `metadata` depend on the data available in index files. See [Series](#series) for details.

## Series

```py
@dataclasses.dataclass
class Series:
    name: str
    description: str
    scale_to_sizes: dict[int, tuple[int, ...]]
    url_pattern: str

    def size(self, scales: typing.Optional[set[int]] = None): ...

    def size_as_string(self, scales: typing.Optional[set[int]] = None): ...

    def index_files(
        self,
        cache_directory: typing.Union[bytes, str, os.PathLike],
        scales: typing.Optional[set[int]] = None,
    ) -> list[pathlib.Path]: ...
```

-   `name` defines the cache subdirectory name.
-   `description` is a copy of the text description in http://data.astrometry.net.
-   `scale_to_sizes` maps each available HEALPix resolution to index files sizes in bytes. The smaller the scale, the larger the number of HEALPix subdivisions.
-   `url_pattern` is the base pattern used to generate file download links.
-   `size` returns the cumulative file sizes for the given scales in bytes. If `scales` is `None`, all the scales available for the series are used.
-   `size_as_string` returns a human-readable string representation of `size`.
-   `index_files` returns index files paths for the given scales (or all available scales if `scales` is `None`). This function downloads files that are not already in the cache directory. `cache_directory` is created if it does not exist. Download automatically resumes for partially downloaded files.

Change the constants `astrometry.CHUNK_SIZE`, `astrometry.DOWNLOAD_SUFFIX` and `astrometry.TIMEOUT` to configure the downloader parameters.

## batches_generator

```py
def batches_generator(
    batch_size: int,
) -> typing.Callable[[int], typing.Iterable[tuple[int, int]]]:
    ...
```

-   `batch_size` sets the size of the generated batches.

`batches_generator` returns a slices generator compatible with `SolutionParameters.slices_generator`. The slices are non-overlapping and non-full slices are ignored. For instance, a batch size of `25` over `83` stars would generate the slices `(0, 25)`, `(25, 50)`, and `(50, 75)`.

## SupportsFloatMapping

```py
class SupportsFloatMapping(typing.Protocol):
    def __getitem__(self, index: typing.SupportsIndex, /) -> typing.SupportsFloat:
        ...
```

# Contribute

## Clone the repository

Clone this repository and pull its submodule:

```sh
git clone --recursive https://github.com/neuromorphicsystems/astrometry.git
cd astrometry
cd astrometry.net
git reset --hard 04e97d5365525ac429dca69b344f9f4d9070b739
git apply ../astrometry.net.patch
cd ..
```

or

```sh
git clone https://github.com/neuromorphicsystems/astrometry.git
cd astrometry
git submodule update --recursive
cd astrometry.net
git reset --hard 04e97d5365525ac429dca69b344f9f4d9070b739
git apply ../astrometry.net.patch
cd ..
```

After building, reset the patch changes to keep the submodule in a clean state

```sh
git reset --hard 04e97d5365525ac429dca69b344f9f4d9070b739
```

## Format and lint

```sh
clang-format -i astrometry_extension/astrometry_extension.c astrometry_extension/astrometry_extension_utilities.h
pip install '.[lint]'
isort .; black .; pyright .
```

## Build a local version

```sh
pip install -e . # use 'CC="ccache clang" pip install -e .' to speed up incremental builds
```

## Publish

1. Bump the version number in _setup.py_.

2. Create a new release on GitHub.

## MSVC compatibility (work in progress)

-   _fitsbin.c_, _kdtree_internal.c_, _kdtree_internal_fits.c_, _solver.c_: replace Variable Length Arrays (VAL) with \_alloca (`type name[size]` -> `type* name = _alloca(size)`)
-   _fitsbin.c_, _fitsioutils.c_, _fitstable.c_: cast `void*` to `char*` to enable pointer arithmetic
-   _anqfits.c_, _verify.c_: `#define debug(args...)` -> `#define debug(...)`
-   _qfits_time.c_: remove `#include <pwd.h>`
-   _log.h_, _keywords.h_: remove `__attribute__` directives
-   _bl.c_: remove redundant bl.inc and bl_nl.c includes
-   replace all POSIX calls (file IO, network, select...). This requires significant effort when targeting the entire Astrometry.net library. It might be less complicated with Astrometry, which uses only a subset of Astrometry.net.
