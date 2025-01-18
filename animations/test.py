import numpy as np
from manim import *
from interval import *
from number_line import *

from manim.utils.tex_templates import _new_ams_template

oswald = _new_ams_template()
oswald.description = "Oswald"
oswald.add_to_preamble(
    r"""
\usepackage[T1]{fontenc}
\usepackage{Oswald}
\renewcommand{\familydefault}{\sfdefault}
\usepackage[frenchmath]{mathastext}
""",
)

# light
background_color = "#ffffff"
default_color = "#000001"
config.background_color = background_color
# dark
# background_color = "#000001"
# default_color = "#ffffff"
# config.background_color = background_color

center = [4.5, 3, 0]

stroke_width = 4

intervals_color = ["#517b77", "#caa000", "#f66000"]
# intervals_color = ["#f98a00", "#009a8e", "#750013"]
# intervals_color = ["#A1C181", "#FCCA46", "#FE7F2D"]
# intervals_color = ["#095473", "#5B95AA", "#349B90", "#B9D6BC", "#F2D2A2"]

offset_color = "#9b001a"

intervals = {
    0: [[-2, 1], [7, 9]],
    1: [[2, 5], [10, 14]],
    2: [[10, 20]],
}


def with_ghosts(intervals):
    intervals_with_stencil_ghost = {}
    for k, v in intervals.items():
        intervals_with_stencil_ghost[k] = []
        for l in v:
            intervals_with_stencil_ghost[k].append([l[0] - 1, l[0]])
            intervals_with_stencil_ghost[k].append([l[1], l[1] + 1])
    return intervals_with_stencil_ghost


def projection(interval):
    intervals_projection = {}
    for k, v in intervals.items():
        if k > 0:
            intervals_projection[k - 1] = []
            for l in v:
                intervals_projection[k - 1].append([(l[0] >> 1) - 1, (l[1] >> 1) + 1])
    return intervals_projection


with register_font("../public/theme/fonts/Oswald-Regular.ttf"):
    Text.set_default(font="Oswald")

MathTex.set_default(tex_template=oswald, font_size=30, color=default_color)

Interval.set_default(stroke_width=stroke_width, color=default_color, tick_size=0.05)


def init_mesh(interval=intervals, stroke_width=stroke_width):
    mesh = []
    colors = []
    levels = []

    for k, v in interval.items():
        mesh.append(VGroup())
        colors.append(intervals_color[k])
        levels.append(k)
        for i in v:
            mesh[-1].add(Interval(k, range=i, stroke_width=stroke_width))

    return VGroup(*mesh), colors, levels


class mesh(Scene):
    def construct(self):
        self.camera.frame_width = 14.2
        self.camera.resize_frame_shape()
        self.camera.frame_center = [4, 2, 0]
        mesh, colors, levels = init_mesh(stroke_width=8)
        # [m.set_color(c) for m, c in zip(mesh, colors)]
        [m.shift(2 * l * UP) for m, l in zip(mesh, levels)]

        self.next_section()
        self.add(mesh)

        meshg, colorsg, levelsg = init_mesh(with_ghosts(intervals))
        [m.set_color(intervals_color[1]) for m in meshg]
        [m.shift(2 * l * UP) for m, l in zip(meshg, levelsg)]
        self.add(meshg)

        self.next_section()

        mesh, colors, levels = init_mesh(projection(intervals))
        [m.set_color(intervals_color[0]) for m in mesh]
        [m.shift(2 * l * UP) for m, l in zip(mesh, levels)]
        self.add(mesh)

        self.next_section()
