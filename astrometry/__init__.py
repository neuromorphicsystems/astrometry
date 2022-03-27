from __future__ import annotations
import astrometry_extension
import dataclasses
import enum
import operator
import pathlib
import threading
import typing
from .series import CHUNK_SIZE as CHUNK_SIZE
from .series import DOWNLOAD_SUFFIX as DOWNLOAD_SUFFIX
from .series import TIMEOUT as TIMEOUT
from .series import size_to_string as size_to_string
from .series import Series as Series
from .series import series_4100 as series_4100
from .series import series_4200 as series_4200
from .series import series_5000 as series_5000
from .series import series_5200 as series_5200
from .series import series_5200_heavy as series_5200_heavy
from .series import series_6000 as series_6000
from .series import series_6100 as series_6100


DEFAULT_LOWER_ARCSEC_PER_PIXEL = 0.1
DEFAULT_UPPER_ARCSEC_PER_PIXEL = 1000.0


@dataclasses.dataclass
class SizeHint:
    lower_arcsec_per_pixel: float
    upper_arcsec_per_pixel: float

    def __post_init__(self):
        assert self.lower_arcsec_per_pixel > 0
        assert self.upper_arcsec_per_pixel > 0
        assert self.lower_arcsec_per_pixel <= self.upper_arcsec_per_pixel


@dataclasses.dataclass
class PositionHint:
    ra_deg: float
    dec_deg: float
    radius_deg: float

    def __post_init__(self):
        assert self.ra_deg >= 0.0 and self.ra_deg < 360.0
        assert self.dec_deg >= -90.0 and self.dec_deg <= 90.0
        assert self.radius_deg >= 0


@dataclasses.dataclass
class Star:
    ra_deg: float
    dec_deg: float
    metadata: dict[str, typing.Any]


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


@dataclasses.dataclass
class Solution:
    solve_id: str
    matches: list[Match]

    def __post_init__(self):
        if len(self.matches) > 0:
            self.matches.sort(key=operator.attrgetter("logodds"), reverse=True)

    def has_match(self) -> bool:
        return len(self.matches) > 0

    def best_match(self) -> Match:
        return self.matches[0]


class Action(enum.Enum):
    STOP = 0
    CONTINUE = 1

    def __bool__(self):
        return self == Action.CONTINUE


class Solver(astrometry_extension.Solver):
    def __init__(self, index_files: list[pathlib.Path]):
        super().__init__([str(path.resolve()) for path in index_files])
        self.solve_id_lock = threading.Lock()
        self.solve_id = 0

    def solve(
        self,
        stars_xs: typing.Iterable[float],
        stars_ys: typing.Iterable[float],
        size_hint: typing.Optional[SizeHint],
        position_hint: typing.Optional[PositionHint],
        solve_id: typing.Optional[str],
        tune_up_logodds_threshold: typing.Optional[float],
        output_logodds_threshold: float,
        logodds_callback: typing.Callable[[list[float]], Action],
    ) -> Solution:
        with self.solve_id_lock:
            self.solve_id += 1
        if solve_id is None:
            solve_id = str(self.solve_id)
        if size_hint is None:
            size_hint = SizeHint(
                lower_arcsec_per_pixel=DEFAULT_LOWER_ARCSEC_PER_PIXEL,
                upper_arcsec_per_pixel=DEFAULT_UPPER_ARCSEC_PER_PIXEL,
            )
        raw_solution = super().solve(
            stars_xs if isinstance(stars_xs, list) else list(stars_xs),
            stars_ys if isinstance(stars_ys, list) else list(stars_ys),
            size_hint.lower_arcsec_per_pixel,
            size_hint.upper_arcsec_per_pixel,
            None
            if position_hint is None
            else (
                position_hint.ra_deg,
                position_hint.dec_deg,
                position_hint.radius_deg,
            ),
            solve_id,
            tune_up_logodds_threshold,
            output_logodds_threshold,
            logodds_callback,
        )
        if raw_solution is None:
            return Solution(solve_id=solve_id, matches=[])
        stars = {
            key: Star(ra_deg=ra_deg, dec_deg=dec_deg, metadata=metadata)
            for key, (ra_deg, dec_deg, metadata) in raw_solution[0].items()
        }
        solution = Solution(
            solve_id=solve_id,
            matches=[
                Match(
                    logodds=logodds,
                    center_ra_deg=center_ra_deg,
                    center_dec_deg=center_dec_deg,
                    scale_arcsec_per_pixel=scale_arcsec_per_pixel,
                    index_path=pathlib.Path(index_path),
                    stars=tuple(stars[key] for key in stars_keys),
                    quad_stars=tuple(stars[key] for key in quad_stars_keys),
                    wcs_fields=wcs_fields,
                )
                for (
                    logodds,
                    center_ra_deg,
                    center_dec_deg,
                    scale_arcsec_per_pixel,
                    index_path,
                    stars_keys,
                    quad_stars_keys,
                    wcs_fields,
                ) in raw_solution[1]
            ],
        )
        return solution
