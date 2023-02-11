import astrometry
import json
import logging
import pathlib

logging.getLogger().setLevel(logging.INFO)

solver = astrometry.Solver(
    astrometry.series_4100.index_files(
        cache_directory=pathlib.Path(__file__).resolve().parent / "astrometry_cache",
        scales={7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
    )
)

# x, y, radius (pixels)
with open("apod4_stars.json") as stars_file:
    stars = sorted(json.load(stars_file), key=lambda star: -star[2])

stars = stars[:50]

print("Solve with default parameters")
solution = solver.solve(
    stars_xs=[star[0] for star in stars],
    stars_ys=[star[1] for star in stars],
    size_hint=astrometry.SizeHint(
        lower_arcsec_per_pixel=163.0,
        upper_arcsec_per_pixel=172.0,
    ),
    position_hint=astrometry.PositionHint(ra_deg=185.8, dec_deg=55.5, radius_deg=10.0),
    solution_parameters=astrometry.SolutionParameters(),
)
if solution.has_match():
    for star in solution.best_match().stars:
        print(star)
else:
    print(
        "Astrometry could not find a match.",
        "Reduce output_logodds_threshold in solution_parameters,",
        "or use another set of index files",
    )
