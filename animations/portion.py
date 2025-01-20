from tkinter import RIGHT
from turtle import color
import numpy as np
from manim import *
from interval import *
from number_line import *

from manim.utils.tex_templates import _new_ams_template

from PIL import Image

color_palette = [
    "#1565C0",
    "#D32F2F",
    "#388E3C",
    "#FBC02D",
    "#7B1FA2",
]


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
default_color = "#9E9E9E"
config.background_color = background_color
# dark
# background_color = "#000001"
# default_color = "#ffffff"
# config.background_color = background_color

center = [2, 3, 0]

stroke_width = 2

# intervals_color = ["#517b77", "#caa000", "#f66000"]
# intervals_color = ["#f98a00", "#009a8e", "#750013"]
# intervals_color = ["#A1C181", "#FCCA46", "#FE7F2D"]
intervals_color = ["#095473", "#5B95AA", "#349B90", "#B9D6BC", "#F2D2A2"]

offset_color = "#9b001a"

domain_i = {
    0: [[-2, 3]],
    1: [[-4, 6]],
    2: [[-8, 12]],
    3: [[-16, 24]],
    4: [[-32, 48]],
}

intervals = {
    0: [[-2, 3]],
    1: [[-2, 4]],
    2: [[-2, 6]],
    3: [[-2, 10]],
    4: [[-2, 18]],
}

first = {
    0: [[-2, 2]],
    1: [[-2, 2]],
    2: [[-2, 2]],
    3: [[-1, 2]],
    4: [[0, 1]],
}

second = {
    0: [[-1, 3]],
    1: [[0, 4]],
    2: [[2, 6]],
    3: [[6, 9]],
    4: [[15, 16]],
}

with register_font("fonts/Oswald-Regular.ttf"):
    Text.set_default(font="Oswald")

MathTex.set_default(tex_template=oswald, font_size=30, color=default_color)

Interval.set_default(stroke_width=stroke_width, color=default_color, tick_size=0.02)


def init_mesh(color, intervals, stroke_width=2, opacity=1, move=0.5):
    mesh = []
    levels = []

    for k, v in intervals.items():
        mesh.append(VGroup())
        levels.append(k)
        for i in v:
            mesh[-1].add(
                Interval(k, range=i).set_stroke(width=stroke_width, opacity=opacity)
            )

    [m.set_color(color) for m in mesh]
    [m.shift(move * l * UP) for m, l in zip(mesh, levels)]
    return VGroup(*mesh)


class portion(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 8
        self.camera.resize_frame_shape()

        limits = VGroup(
            Line([0, -1, 0], [0, 5, 0], color=color_palette[2]).set_stroke(
                width=1.5, opacity=0.6
            ),
            Line([1, -1, 0], [1, 5, 0], color=color_palette[2]).set_stroke(
                width=1.5, opacity=0.6
            ),
        )
        # self.add(line_1, line_2)

        domain = init_mesh(default_color, domain_i, stroke_width=1.0, move=1)
        # self.add(domain)

        mesh = init_mesh(color_palette[0], intervals, stroke_width=1.5, move=1)
        mesh.z_index = 10
        mesh[0].z_index = 10

        center = mesh[0].get_center()
        self.camera.frame_center = center
        self.play(Create(mesh[0], run_time=4))
        self.wait()

        self.next_section()
        center = mesh.get_center()

        self.play(
            AnimationGroup(
                Create(domain), self.camera.frame.animate.move_to(center), run_time=2
            )
        )
        self.wait()

        self.play(
            Create(limits),
        )
        self.wait()

        ex_1 = init_mesh(color_palette[4], first, stroke_width=2.0, move=1)
        for e in ex_1:
            e.z_index = 20

        g = VGroup(ex_1)

        self.play(
            Create(ex_1[-1]),
        )
        self.wait()

        self.next_section()

        lines_left = [
            ([0, 4, 0], [-1 / (1 << 3), 3, 0]),
            ([-1 / (1 << 3), 3, 0], [-2 / (1 << 2), 2, 0]),
            ([-2 / (1 << 2), 2, 0], [-2 / (1 << 1), 1, 0]),
            ([-2 / (1 << 1), 1, 0], [-2 / (1 << 0), 0, 0]),
        ]
        lines_right = [
            ([1 / (1 << 4), 4, 0], [2 / (1 << 3), 3, 0]),
            ([2 / (1 << 3), 3, 0], [2 / (1 << 2), 2, 0]),
            ([2 / (1 << 2), 2, 0], [2 / (1 << 1), 1, 0]),
            ([2 / (1 << 1), 1, 0], [2 / (1 << 0), 0, 0]),
        ]

        lines = VGroup()
        for l_l, l_r in zip(lines_left, lines_right):
            lines.add(Line(*l_l, color=color_palette[4], stroke_width=1.3))
            lines.add(Line(*l_r, color=color_palette[4], stroke_width=1.3))

        i = 0
        for j in range(0, len(lines), 2):
            self.play(
                Create(lines[j]),
                Create(lines[j + 1]),
            )
            self.play(Create(ex_1[-2 - i]))
            i += 1
        self.wait()

        self.next_section()

        ex_2 = init_mesh(color_palette[3], second, stroke_width=2.0, move=1)
        for e in ex_2:
            e.z_index = 20

        lines_left = [
            ([15 / (1 << 4), 4, 0], [6 / (1 << 3), 3, 0]),
            ([6 / (1 << 3), 3, 0], [2 / (1 << 2), 2, 0]),
            ([2 / (1 << 2), 2, 0], [0 / (1 << 1), 1, 0]),
            ([0 / (1 << 1), 1, 0], [-1 / (1 << 0), 0, 0]),
        ]
        lines_right = [
            ([16 / (1 << 4), 4, 0], [9 / (1 << 3), 3, 0]),
            ([9 / (1 << 3), 3, 0], [6 / (1 << 2), 2, 0]),
            ([6 / (1 << 2), 2, 0], [4 / (1 << 1), 1, 0]),
            ([4 / (1 << 1), 1, 0], [3 / (1 << 0), 0, 0]),
        ]

        lines2 = VGroup()
        for l_l, l_r in zip(lines_left, lines_right):
            lines2.add(Line(*l_l, color=color_palette[3], stroke_width=1.3))
            lines2.add(Line(*l_r, color=color_palette[3], stroke_width=1.3))

        self.play(
            Create(ex_2[-1]),
        )
        self.wait()
        i = 0
        for j in range(0, len(lines2), 2):
            self.play(
                Create(lines2[j]),
                Create(lines2[j + 1]),
                #     Line(*l_r, color=default_color, stroke_width=1.3),
                # )
            )
            self.play(Create(ex_2[-2 - i]))
            i += 1
        self.wait()
