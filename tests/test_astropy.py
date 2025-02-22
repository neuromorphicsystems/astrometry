import pathlib

import astrometry

dirname = pathlib.Path(__file__).resolve().parent

with astrometry.Solver(
    astrometry.series_5200.index_files(
        cache_directory=dirname / "astrometry_cache",
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

    print(solution.best_match().astropy_wcs())
