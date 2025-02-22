import pathlib
import random

import numpy
import psutil

import astrometry

dirname = pathlib.Path(__file__).resolve().parent

SOLVES: int = 100

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

process = psutil.Process()


def size_to_string(size: int) -> str:
    if size < 0:
        size = -size
        prefix = "-"
    else:
        prefix = ""
    if size < 1000:
        return f"{prefix}{size} B"
    if size < 1000000:
        return f"{prefix}{size / 1e3:.2f} kB"
    if size < 1000000000:
        return f"{prefix}{size / 1e6:.2f} MB"
    if size < 1000000000000:
        return f"{prefix}{size / 1e9:.2f} GB"
    return f"{prefix}{size / 1e12:.2f} TB"


memory_samples = numpy.zeros(SOLVES, dtype=numpy.int64)

for index in range(0, SOLVES):
    with astrometry.Solver(
        astrometry.series_5200.index_files(
            cache_directory=dirname / "astrometry_cache",
            scales={6},
        )
    ) as solver:
        solution = solver.solve(
            stars=stars,
            size_hint=None,
            position_hint=None,
            solution_parameters=astrometry.SolutionParameters(),
        )
        assert solution.has_match()
        del solution
    memory_sample = process.memory_full_info().uss
    memory_samples[index] = memory_sample
    print(
        f"(repeated solver allocation with context manager) solve {index + 1} / {SOLVES}, memory usage: {size_to_string(memory_sample)}"
    )
deltas = numpy.diff(memory_samples)
repeated_solver_allocation_with_context_manager_leak = (
    size_to_string(int(numpy.min(deltas[deltas >= 0]))),
    size_to_string(int(round((memory_samples[-1] - memory_samples[0]) / SOLVES))),
    size_to_string(int(numpy.max(deltas[deltas >= 0]))),
)
print(f"{repeated_solver_allocation_with_context_manager_leak=}")

print("")
for index in range(0, SOLVES):
    solver = astrometry.Solver(
        astrometry.series_5200.index_files(
            cache_directory=dirname / "astrometry_cache",
            scales={6},
        )
    )
    solution = solver.solve(
        stars=stars,
        size_hint=None,
        position_hint=None,
        solution_parameters=astrometry.SolutionParameters(),
    )
    assert solution.has_match()
    solver.close()
    del solution
    del solver
    memory_sample = process.memory_full_info().uss
    memory_samples[index] = memory_sample
    print(
        f"(repeated solver allocation) solve {index + 1} / {SOLVES}, memory usage: {size_to_string(memory_sample)}"
    )
deltas = numpy.diff(memory_samples)
repeated_solver_allocation_leak = (
    size_to_string(int(numpy.min(deltas[deltas >= 0]))),
    size_to_string(int(round((memory_samples[-1] - memory_samples[0]) / SOLVES))),
    size_to_string(int(numpy.max(deltas[deltas >= 0]))),
)
print(f"{repeated_solver_allocation_leak=}")

print("")
with astrometry.Solver(
    astrometry.series_5200.index_files(
        cache_directory=dirname / "astrometry_cache",
        scales={6},
    )
) as solver:
    for index in range(0, SOLVES):
        solution = solver.solve(
            stars=stars,
            size_hint=None,
            position_hint=None,
            solution_parameters=astrometry.SolutionParameters(
                output_logodds_threshold=3.0,
            ),
        )
        assert solution.has_match()
        del solution
        memory_sample = process.memory_full_info().uss
        memory_samples[index] = memory_sample
        print(
            f"(repeated solves) solve {index + 1} / {SOLVES}, memory usage: {size_to_string(memory_sample)}"
        )
deltas = numpy.diff(memory_samples)
repeated_solves_leak = (
    size_to_string(int(numpy.min(deltas[deltas >= 0]))),
    size_to_string(int(round((memory_samples[-1] - memory_samples[0]) / SOLVES))),
    size_to_string(int(numpy.max(deltas[deltas >= 0]))),
)
print(f"{repeated_solves_leak=}")

print("")
with astrometry.Solver(
    astrometry.series_5200.index_files(
        cache_directory=dirname / "astrometry_cache",
        scales={6},
    )
) as solver:
    for index in range(0, SOLVES):
        stars = []
        for _ in range(0, 12):
            stars.append([random.random() * 1000.0, random.random() * 1000.0])
        solution = solver.solve(
            stars=stars,
            size_hint=None,
            position_hint=None,
            solution_parameters=astrometry.SolutionParameters(
                output_logodds_threshold=3.0,
            ),
        )
        del solution
        del stars
        memory_sample = process.memory_full_info().uss
        memory_samples[index] = memory_sample
        print(
            f"(repeated solves, random stars) solve {index + 1} / {SOLVES}, memory usage: {size_to_string(memory_sample)}"
        )
deltas = numpy.diff(memory_samples)
repeated_solves_random_stars_leak = (
    size_to_string(int(numpy.min(deltas[deltas >= 0]))),
    size_to_string(int(round((memory_samples[-1] - memory_samples[0]) / SOLVES))),
    size_to_string(int(numpy.max(deltas[deltas >= 0]))),
)
print(f"{repeated_solves_random_stars_leak=}")

print("")
print(f"{repeated_solver_allocation_with_context_manager_leak=}")
print(f"{repeated_solver_allocation_leak=}")
print(f"{repeated_solves_leak=}")
print(f"{repeated_solves_random_stars_leak=}")
