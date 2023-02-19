import astrometry
import json
import logging
import pathlib

dirname = pathlib.Path(__file__).resolve().parent

logging.getLogger().setLevel(logging.INFO)

solver = astrometry.Solver(
    astrometry.series_5000.index_files(
        cache_directory=dirname / "astrometry_cache",
        scales={5},
    )
)

with open(dirname / "astrosat_stars.json") as stars_file:
    stars = json.load(stars_file)

print("Solve with default parameters")
solution = solver.solve(
    stars=stars,
    size_hint=astrometry.SizeHint(
        lower_arcsec_per_pixel=0.40,
        upper_arcsec_per_pixel=0.43,
    ),
    position_hint=astrometry.PositionHint(
        ra_deg=53.0523, dec_deg=-27.8858, radius_deg=0.33
    ),
    solution_parameters=astrometry.SolutionParameters(
        sip_order=0,
        tune_up_logodds_threshold=None,
        output_logodds_threshold=20.0,
        logodds_callback=lambda logodds_list: astrometry.Action.STOP,
        slices_generator=astrometry.batches_generator(50),
    ),
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
