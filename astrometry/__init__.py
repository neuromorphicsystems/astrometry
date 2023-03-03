from __future__ import annotations
import astrometry_extension
import dataclasses
import enum
import json
import math
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


class SupportsFloatMapping(typing.Protocol):
    def __getitem__(self, index: typing.SupportsIndex, /) -> typing.SupportsFloat:
        ...


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


class Action(enum.IntEnum):
    STOP = 0
    CONTINUE = 1

    def __bool__(self):
        return self == Action.CONTINUE


class Parity(enum.IntEnum):
    NORMAL = 0
    FLIP = 1
    BOTH = 2


def batches_generator(
    batch_size: int,
) -> typing.Callable[[int], typing.Iterable[tuple[int, int]]]:
    def slices_generator(stars_count: int) -> typing.Iterable[tuple[int, int]]:
        if stars_count < batch_size:
            yield (0, stars_count)
        for start in range(0, stars_count - batch_size + 1, batch_size):
            yield (start, start + batch_size)

    return slices_generator


@dataclasses.dataclass
class SolutionParameters:
    solve_id: typing.Optional[str] = None
    uniformize_index: bool = True  # uniformize field stars at the matched index scale before verifying a match
    deduplicate: bool = True  # de-duplicate field stars before verifying a match
    sip_order: int = 3  # 0 means "no SIP distortion"
    sip_inverse_order: int = 0  # 0 means "same as sip_order"
    distance_from_quad_bonus: bool = True  # Assume that stars far from the matched quad will have larger positional variance?
    positional_noise_pixels: float = 1.0
    distractor_ratio: float = 0.25
    code_tolerance_l2_distance: float = 0.01
    minimum_quad_size_pixels: typing.Optional[
        float
    ] = None  # None means "use minimum_quad_size_fraction"
    minimum_quad_size_fraction: float = 0.1
    maximum_quad_size_pixels: float = 0.0  # 0.0 means "infinity"
    maximum_quads: int = 0  # number of field quads to try, 0 means no limit
    maximum_matches: int = 0  # number of quad matches to try, 0 means no limit
    parity: Parity = Parity.BOTH
    tune_up_logodds_threshold: typing.Optional[float] = 14.0  # None means "no tune-up"
    output_logodds_threshold: float = 21.0
    slices_generator: typing.Callable[
        [int], typing.Iterable[tuple[int, int]]
    ] = batches_generator(25)
    logodds_callback: typing.Callable[[list[float]], Action] = lambda _: Action.CONTINUE

    def __post_init__(self):
        assert self.sip_order >= 0
        assert self.sip_inverse_order >= 0
        assert self.positional_noise_pixels >= 0.0
        assert self.distractor_ratio > 0.0 and self.distractor_ratio <= 1.0
        assert self.code_tolerance_l2_distance >= 0.0
        assert (
            self.minimum_quad_size_pixels is None
            or self.minimum_quad_size_pixels >= 0.0
        )
        assert self.minimum_quad_size_fraction >= 0.0
        assert self.maximum_quad_size_pixels >= 0.0
        assert self.sip_order > 0 or self.tune_up_logodds_threshold is None


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

    def to_json(self):
        solution_as_dict = dataclasses.asdict(self)
        for match in solution_as_dict["matches"]:
            match["index_path"] = str(match["index_path"])
        return json.dumps(solution_as_dict)


class Solver(astrometry_extension.Solver):
    def __init__(self, index_files: list[pathlib.Path]):
        super().__init__([str(path.resolve()) for path in index_files])
        self.solve_id_lock = threading.Lock()
        self.solve_id = 0

    def solve(
        self,
        stars: typing.Iterable[SupportsFloatMapping],
        size_hint: typing.Optional[SizeHint],
        position_hint: typing.Optional[PositionHint],
        solution_parameters: SolutionParameters,
    ) -> Solution:
        with self.solve_id_lock:
            self.solve_id += 1

        stars_xs: list[float] = []
        stars_ys: list[float] = []
        for star in stars:
            stars_xs.append(float(star[0]))
            stars_ys.append(float(star[1]))
        slices_starts: list[int] = []
        slices_ends: list[int] = []
        for star_slice in solution_parameters.slices_generator(len(stars_xs)):
            assert star_slice[0] >= 0
            assert star_slice[1] <= len(stars_xs)
            assert star_slice[0] < star_slice[1]
            slices_starts.append(star_slice[0])
            slices_ends.append(star_slice[1])
        assert len(slices_starts) > 0
        solve_id = (
            str(self.solve_id)
            if solution_parameters.solve_id is None
            else solution_parameters.solve_id
        )
        if size_hint is None:
            size_hint = SizeHint(
                lower_arcsec_per_pixel=DEFAULT_LOWER_ARCSEC_PER_PIXEL,
                upper_arcsec_per_pixel=DEFAULT_UPPER_ARCSEC_PER_PIXEL,
            )
        if solution_parameters.minimum_quad_size_pixels is None:
            star_x_minimum = math.inf
            star_x_maximum = 0.0
            for star_x in stars_xs:
                star_x_minimum = min(star_x_minimum, star_x)
                star_x_maximum = max(star_x_maximum, star_x)
            star_y_minimum = math.inf
            star_y_maximum = 0.0
            for star_y in stars_ys:
                star_y_minimum = min(star_y_minimum, star_y)
                star_y_maximum = max(star_y_maximum, star_y)
            extent = min(
                star_x_maximum - star_x_minimum, star_y_maximum - star_y_minimum
            )
            if extent < 0.0:
                minimum_quad_size_pixels = 0.0
            else:
                minimum_quad_size_pixels = (
                    solution_parameters.minimum_quad_size_fraction * extent
                )
        else:
            minimum_quad_size_pixels = solution_parameters.minimum_quad_size_pixels
        assert (
            solution_parameters.maximum_quad_size_pixels == 0.0
            or minimum_quad_size_pixels < solution_parameters.maximum_quad_size_pixels
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
            solution_parameters.uniformize_index,
            solution_parameters.deduplicate,
            solution_parameters.sip_order,
            solution_parameters.sip_inverse_order,
            solution_parameters.distance_from_quad_bonus,
            solution_parameters.positional_noise_pixels,
            solution_parameters.distractor_ratio,
            solution_parameters.code_tolerance_l2_distance,
            minimum_quad_size_pixels,
            solution_parameters.maximum_quad_size_pixels,
            solution_parameters.maximum_quads,
            solution_parameters.maximum_matches,
            int(solution_parameters.parity),
            solution_parameters.tune_up_logodds_threshold,
            solution_parameters.output_logodds_threshold,
            slices_starts,
            slices_ends,
            solution_parameters.logodds_callback,
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
