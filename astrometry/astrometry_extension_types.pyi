import typing

class Solver:
    def __init__(self, paths: list[str]): ...
    def close(self): ...
    def solve(
        self,
        stars_xs: list[float],
        stars_ys: list[float],
        scale_lower: float,
        scale_upper: float,
        position_hint: typing.Optional[tuple[float, float, float]],
        solve_id: str,
        uniformize_index: int,
        deduplicate: bool,
        sip_order: int,
        sip_inverse_order: int,
        distance_from_quad_bonus: int,
        positional_noise_pixels: float,
        distractor_ratio: float,
        code_tolerance_l2_distance: float,
        minimum_quad_size_pixels: float,
        maximum_quad_size_pixels: float,
        maximum_quads: int,
        maximum_matches: int,
        parity: int,
        tune_up_logodds_threshold: typing.Optional[float],
        output_logodds_threshold: float,
        slices_starts: list[int],
        slices_ends: list[int],
        logodds_callback: typing.Callable[[list[float]], typing.Any],
    ) -> typing.Optional[
        tuple[
            # stars
            dict[
                tuple[int, int],  # (index_id, star_id)
                tuple[
                    float,  # ra_deg
                    float,  # dec_deg
                    dict[str, typing.Any],  # metadata
                ],
            ],
            # matches
            tuple[
                tuple[
                    float,  # logodds
                    float,  # center_ra_deg
                    float,  # center_dec_deg
                    float,  # scale
                    str,  # index_path
                    tuple[
                        tuple[int, int], ...
                    ],  # stars_keys ((index_id, star_id), ...)
                    tuple[
                        tuple[int, int], ...
                    ],  # quad_stars_keys ((index_id, star_id), ...)
                    dict[str, tuple[typing.Any, str]],  # wcs_fields
                ],
                ...,
            ],
        ]
    ]: ...
