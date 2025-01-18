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

intervals = {
    0: [[0, 1]],
    1: [[2, 3]],
    2: [[6, 7]],
    3: [[14, 16]],
    4: [[32, 48]],
}

int2coarse = {4: [[40, 42]]}
int2new = {3: [[20, 21]]}

new_intervals = {
    0: [[0, 1]],
    1: [[2, 3]],
    2: [[6, 7]],
    3: [[14, 16], [20, 21]],
    4: [[32, 40], [42, 48]],
}

proj_i = {
    1: [[5, 6]],
}

ghosts_i = {}
for k, v in intervals.items():
    ghosts_i[k] = []
    for i in v:
        ghosts_i[k].append([i[0] - 1, i[0]])
        ghosts_i[k].append([i[1], i[1] + 1])

pred_i_1 = {}
for k, v in intervals.items():
    if k > 0:
        pred_i_1[k - 1] = []
        for i in v:
            pred_i_1[k - 1].append([(i[0] >> 1) - 1, (i[1] - 1 >> 1) + 2])

pred_i_2 = {}
for k, v in intervals.items():
    if k > 0:
        pred_i_2[k - 1] = []
        for i in v:
            pred_i_2[k - 1].append([(i[0] >> 1) - 1, (i[1] - 1 >> 1) + 2])
    if k > 1:
        for i in v:
            pred_i_2[k - 2].append([(i[0] >> 2) - 1, (i[1] - 1 >> 2) + 2])

unions = {4: []}
for k in range(4, 0, -1):
    unions[k - 1] = [[i[0] >> 1, (i[1] - 1 >> 1) + 1] for i in intervals[k]]
    unions[k - 1].extend([[i[0] >> 1, (i[1] - 1 >> 1) + 1] for i in unions[k]])

all_domain = {0: [[0, 3]], 1: [[0, 6]], 2: [[0, 12]], 3: [[0, 24]], 4: [[0, 48]]}

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


class mesh(MovingCameraScene):
    def construct(self):
        self.camera.frame_width = 5
        self.camera.resize_frame_shape()

        domain = init_mesh(default_color, all_domain, stroke_width=1, opacity=0.2)

        mesh = init_mesh(color_palette[0], intervals, stroke_width=1.5, move=0)
        center = mesh.get_center()
        self.camera.frame_center = center
        self.add(mesh)

        frame = self.renderer.get_frame()
        img = Image.fromarray(frame)
        img.save("intial_mesh.png")

        self.next_section()

        self.play(
            FadeIn(domain, run_time=2),
            *[
                m.animate(run_time=3).shift(0.5 * l * UP)
                for m, l in zip(mesh, range(5))
            ],
            self.camera.frame.animate.move_to(domain)
        )
        self.wait()

        self.next_section()

        ghosts = init_mesh(color_palette[4], ghosts_i, stroke_width=1)
        self.play(
            FadeIn(ghosts),
        )
        self.wait()

        self.next_section()

        union = init_mesh(color_palette[3], unions, stroke_width=1)

        self.play(FadeOut(ghosts), FadeIn(union))
        self.wait()

        self.next_section()

        pred_1 = init_mesh(color_palette[4], pred_i_1, stroke_width=1)

        self.play(FadeOut(union), FadeIn(pred_1))

        self.wait()

        self.next_section()

        coarse = init_mesh(color_palette[2], int2coarse, stroke_width=2)
        tmp = init_mesh(color_palette[0], int2new, stroke_width=1)
        new_int = init_mesh(color_palette[0], new_intervals, stroke_width=1.5)

        g = VGroup(coarse, tmp)

        self.play(FadeIn(coarse), self.camera.frame.animate.scale(0.6).move_to(g))
        self.wait()

        self.next_section()

        arrow = Arrow(
            start=coarse.get_bottom() - np.array([0, 0.025, 0]),
            end=tmp.get_top(),
            color=default_color,
            max_tip_length_to_length_ratio=0.1,
            stroke_width=1,
        )
        self.play(Create(arrow))
        self.play(FadeIn(new_int), FadeOut(arrow, pred_1, mesh, coarse))
        self.wait()

        self.next_section()

        self.play(
            FadeOut(new_int),
            FadeIn(mesh),
            self.camera.frame.animate.scale(2).move_to(domain),
        )
        self.wait()

        self.next_section()

        pred_2 = init_mesh(color_palette[4], pred_i_2, stroke_width=1)

        self.play(FadeIn(pred_2))
        self.wait()

        self.next_section()

        proj = init_mesh(color_palette[2], proj_i, stroke_width=1)

        self.play(FadeIn(proj))
        self.wait()
