import argparse
import astropy.wcs
import json
import numpy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", help="path to a JSON file contanining WCS fields")
parser.add_argument("output", help="path to the output SVG file")
parser.add_argument("left", type=float, help="left coordinate in pixels")
parser.add_argument("bottom", type=float, help="bottom coordinate in pixels")
parser.add_argument("top", type=float, help="top coordinate in pixels")
parser.add_argument("right", type=float, help="right coordinate in pixels")
parser.add_argument(
    "--grid-spacing-horizontal",
    type=float,
    default=50,
    help="distortion grid x spacing in pixels",
)
parser.add_argument(
    "--grid-spacing-vertical",
    type=float,
    default=50,
    help="distortion grid y spacing in pixels",
)
parser.add_argument(
    "--subdivisions-horizontal",
    type=int,
    default=50,
    help="number of x subdivisions between grid vertices",
)
parser.add_argument(
    "--subdivisions-vertical",
    type=int,
    default=50,
    help="number of y subdivisions between grid vertices",
)
parser.add_argument(
    "--foc2pix",
    action="store_true",
    help="plot the foc2pix (part of world2pix) distortion instead of the pix2foc (part of pix2world) distortion",
)
parser.add_argument(
    "--padding",
    type=float,
    default=20.0,
    help="padding in pixels",
)
parser.add_argument(
    "--stroke-width",
    type=float,
    default=3.0,
    help="stroke width in pixels",
)
parser.add_argument(
    "--distorted-stroke-width",
    type=float,
    default=3.0,
    help="stroke width in pixels",
)
parser.add_argument(
    "--exaggerate",
    type=float,
    default=1.0,
    help="multiply the distortion by this factor",
)
args = parser.parse_args()

assert args.left < args.right
assert args.bottom < args.top

with open(args.input) as wcs_fields_file:
    wcs = astropy.wcs.WCS(
        {key: tuple(value) for key, value in json.load(wcs_fields_file).items()}
    )

if args.foc2pix:
    transformation = wcs.sip.foc2pix
else:
    transformation = wcs.sip.pix2foc

lines: list[numpy.ndarray] = []
distorted_lines: list[numpy.ndarray] = []
x_range = [numpy.Infinity, -numpy.Infinity]
y_range = [numpy.Infinity, -numpy.Infinity]
offset = transformation(
    numpy.array([[(args.left + args.right) / 2, (args.bottom + args.top) / 2]]), 0
)
for y in numpy.arange(
    args.bottom, args.top + args.grid_spacing_vertical, args.grid_spacing_vertical
):
    x = numpy.arange(
        args.left,
        args.right,
        args.grid_spacing_horizontal / args.subdivisions_horizontal,
    )
    lines.append(numpy.stack((x, numpy.repeat(y, len(x)))).transpose())
    distorted_lines.append(transformation(lines[-1], 0) - offset)
    if args.exaggerate != 1.0:
        distorted_lines[-1] = (
            lines[-1] + numpy.subtract(distorted_lines[-1], lines[-1]) * args.exaggerate
        )
    x_range[0] = min(x_range[0], lines[-1][:, 0].min())
    x_range[1] = max(x_range[1], lines[-1][:, 0].max())
    y_range[0] = min(y_range[0], lines[-1][:, 1].min())
    y_range[1] = max(y_range[1], lines[-1][:, 1].max())
    x_range[0] = min(x_range[0], distorted_lines[-1][:, 0].min())
    x_range[1] = max(x_range[1], distorted_lines[-1][:, 0].max())
    y_range[0] = min(y_range[0], distorted_lines[-1][:, 1].min())
    y_range[1] = max(y_range[1], distorted_lines[-1][:, 1].max())
for x in numpy.arange(
    args.left, args.right + args.grid_spacing_horizontal, args.grid_spacing_horizontal
):
    y = numpy.arange(
        args.left,
        args.right,
        args.grid_spacing_vertical / args.subdivisions_vertical,
    )
    lines.append(numpy.stack((numpy.repeat(x, len(y)), y)).transpose())
    distorted_lines.append(transformation(lines[-1], 0) - offset)
    if args.exaggerate != 1.0:
        distorted_lines[-1] = (
            lines[-1] + numpy.subtract(distorted_lines[-1], lines[-1]) * args.exaggerate
        )
    x_range[0] = min(x_range[0], lines[-1][:, 0].min())
    x_range[1] = max(x_range[1], lines[-1][:, 0].max())
    y_range[0] = min(y_range[0], lines[-1][:, 1].min())
    y_range[1] = max(y_range[1], lines[-1][:, 1].max())
    x_range[0] = min(x_range[0], distorted_lines[-1][:, 0].min())
    x_range[1] = max(x_range[1], distorted_lines[-1][:, 0].max())
    y_range[0] = min(y_range[0], distorted_lines[-1][:, 1].min())
    y_range[1] = max(y_range[1], distorted_lines[-1][:, 1].max())
if (y_range[1] - y_range[0] + args.padding * 2) * 1280 < (
    x_range[1] - x_range[0] + args.padding * 2
) * 720:
    width = 1280
    height = (
        (y_range[1] - y_range[0] + args.padding * 2)
        / (x_range[1] - x_range[0] + args.padding * 2)
        * 1280
    )
else:
    width = (
        (x_range[1] - x_range[0] + args.padding * 2)
        / (y_range[1] - y_range[0] + args.padding * 2)
        * 720
    )
    height = 720


with open(args.output, "w") as output_file:
    output_file.write(
        f'<svg version="1.1" viewBox="{x_range[0] - args.padding} {y_range[0] - args.padding} {x_range[1] - x_range[0] + 2 * args.padding} {y_range[1] - y_range[0] + 2 * args.padding}" style="width: {width}px; height: {height}px" xmlns="http://www.w3.org/2000/svg">\n'
    )
    for line in lines:
        points = " L".join(f"{point[0]},{point[1]}" for point in line)
        output_file.write(
            f'<path d="M{points}" fill="transparent" stroke="#888888" stroke-width="{args.stroke_width}" />'
        )
    for line in distorted_lines:
        points = " L".join(f"{point[0]},{point[1]}" for point in line)
        output_file.write(
            f'<path d="M{points}" fill="transparent" stroke="#EA4335" stroke-width="{args.distorted_stroke_width}" />'
        )
    output_file.write("</svg>")
