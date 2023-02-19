import astrometry
import logging
import pathlib

logging.getLogger().setLevel(logging.INFO)

solver = astrometry.Solver(
    astrometry.series_5200.index_files(
        cache_directory=pathlib.Path(__file__).resolve().parent / "astrometry_cache",
        scales={6},
    )
)

stars = [
    [283.7181611390329, 250.9437436782978],
    [388.9140568247906, 656.5003281719216],
    [732.9210858972549, 473.66395545775106],
    [401.03459504299843, 253.788113189415],
    [312.6591868096163, 624.7527729425295],
    [694.6844564647456, 606.8371776658344],
    [741.7233477959561, 344.41284826261443],
    [867.3574610200455, 672.014835980283],
    [1178.0701325872994, 120.66335820053426],
    [1063.546651153479, 593.7844603550848],
    [1266.479933601124, 478.6594707983611],
    [286.69070190952704, 422.170016812049],
    [401.12779619355155, 16.13543616977013],
    [393.1902113796317, 485.8601927863897],
    [865.3547639559572, 614.3599340062373],
    [205.12103484692776, 698.1847350789413],
    [504.2664247977979, 214.23557044789035],
    [-67.78648235582016, 646.7413890177716],
    [202.88444768690894, 111.24830187635557],
    [747.2580778813443, 116.51880571011176],
    [339.1627757703069, 86.60739435924549],
    [592.1438288540525, 508.6376406353861],
]

for sip_order in range(0, 7):
    print(f"{sip_order=}")
    solution = solver.solve(
        stars=stars,
        size_hint=astrometry.SizeHint(
            lower_arcsec_per_pixel=1.0,
            upper_arcsec_per_pixel=2.0,
        ),
        position_hint=astrometry.PositionHint(
            ra_deg=65.7, dec_deg=36.2, radius_deg=1.0
        ),
        solution_parameters=astrometry.SolutionParameters(
            sip_order=sip_order,
            tune_up_logodds_threshold=None if sip_order == 0 else 14.0,
        ),
    )
    if solution.has_match():
        print(solution.best_match().wcs_fields)
    else:
        print(
            "Astrometry could not find a match.",
            "Reduce output_logodds_threshold in solution_parameters,",
            "or use another set of index files",
        )
