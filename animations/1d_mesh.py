from tkinter import RIGHT
from turtle import color
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

center = [2, 3, 0]

stroke_width = 6

intervals_color = ["#517b77", "#caa000", "#f66000"]
# intervals_color = ["#f98a00", "#009a8e", "#750013"]
# intervals_color = ["#A1C181", "#FCCA46", "#FE7F2D"]
# intervals_color = ["#095473", "#5B95AA", "#349B90", "#B9D6BC", "#F2D2A2"]

offset_color = "#9b001a"

intervals = {
    0: [[0, 2], [5, 6]],
    1: [[4, 5], [6, 8], [9, 10]],
    2: [[10, 12], [16, 18]],
}

intervals_ghost = {
    0: [[-1, 0], [2, 3], [4, 5], [6, 7]],
    1: [[3, 4], [5, 6], [8, 9], [10, 11]],
    2: [[9, 10], [12, 13], [15, 16], [18, 19]],
}

intervals_all = {
    0: [[-1, 3], [4, 7]],
    1: [[3, 11]],
    2: [[9, 13], [15, 19]],
}

with register_font("fonts/Oswald-Regular.ttf"):
    Text.set_default(font="Oswald")

MathTex.set_default(tex_template=oswald, font_size=30, color=default_color)

Interval.set_default(stroke_width=stroke_width, color=default_color, tick_size=0.1)

def init_mesh(interval=intervals):
    mesh = []
    colors = []
    levels = []

    for k, v in interval.items():
        mesh.append(VGroup())
        colors.append(intervals_color[k])
        levels.append(k)
        for i in v:
            mesh[-1].add(Interval(k, range=i))

    return VGroup(*mesh), colors, levels

def init_axe():
    level_axe = MyNumberLine(x_range=[0, 2.5], unit_size=2, include_tip=True, tip_length=0.1, rotation=90*DEGREES, label_direction=LEFT, color=default_color)
    level_axe.add_numbers()
    level_axe.shift(LEFT + 2.5*UP)

    level_label = MathTex("Level")
    level_label.next_to(level_axe, LEFT)

    return VGroup(level_axe, level_label)

def init_dx(ug):
    dx = MathTex(r"\Delta x = 2^{-level}", font_size=20)
    dx.move_to([2, 1, 0])

    level_info = VGroup()
    for i in range(len(ug)):
        itext = MathTex(f"\\frac{{1}}{{ {1<<i} }}" if i != 0 else "1", font_size=16)
        itext.move_to(ug[i].get_center() - [0, .3, 0])
        l_text = MathTex(f"level \quad {i}", font_size=16)
        l_text.move_to([-1, i, 0])
        level_info.add(VGroup(l_text, itext))

    return dx, level_info

class example_1d(MovingCameraScene):
    def construct(self):
        self.next_section()
        mesh, colors, levels = init_mesh()
        self.camera.frame.move_to(mesh)
        self.play(Create(mesh))#, run_time=2)

        self.play(*[m.animate.set_color(c) for m, c in zip(mesh, colors)], run_time=2)
        self.wait()

        self.next_section()

        axe = init_axe()
        all = VGroup(axe, mesh)

        self.play(*[m.animate.shift(2*l*UP) for m, l in zip(mesh, levels)],
                  Create(axe),
                  self.camera.frame.animate.move_to(all)
        )
        self.wait(1)

        self.next_section()

        unit_mesh = [Interval(l, range=[0, 1], stroke_width=stroke_width*0.5, color=c, tick_size=0.05) for m, l, c in zip(mesh, levels, colors)]
        [m.shift(l*UP) for m, l in zip(unit_mesh, levels)]

        ug = VGroup(*unit_mesh)
        mesh_copy = mesh.copy()
        g = VGroup(*mesh)
        self.camera.frame.save_state()
        g.save_state()
        self.play(FadeOut(axe),
                  Transform(g, ug),
                  self.camera.frame.animate.scale(0.5).move_to(ug)
        )

        self.wait()

        dx, level_info = init_dx(ug)
        self.play(Create(dx))

        self.wait()

        self.play(Create(level_info[0]))
        self.play(Create(level_info[1]))
        self.play(Create(level_info[2]))

        self.wait()

        self.next_section()

        ug = VGroup(*mesh_copy)
        self.play(FadeOut(dx),
                  FadeOut(level_info),
                  FadeIn(axe),
                  Transform(g, ug),
                  self.camera.frame.animate.restore()
        )
        self.wait()

        self.play(*[m.animate.add_cell_numbers(font_size=24, color=default_color) for me in mesh_copy for m in me])

        center = MathTex("c_i^l = \\left(i +\\frac{1}{2}\\right)\Delta x", font_size=36)
        center.move_to([7, 5, 0])

        self.play(Create(center))
        self.wait(1)

        r0 = Rectangle(width=1, height=0.5, stroke_color=default_color)
        r0.move_to([0.5, 0.3, 0])
        r1 = Rectangle(width=.5, height=0.5, stroke_color=default_color)
        r1.move_to([2.25, 2.3, 0])

        c0 = MathTex("0.5", font_size=20, color=default_color)
        c0.move_to([0.5, 0.3, 0])
        c1 = MathTex("1.5", font_size=20, color=default_color)
        c1.move_to([1.5, 0.3, 0])
        c2 = MathTex("5.5", font_size=20, color=default_color)
        c2.move_to([5.5, 0.3, 0])
        c3 = MathTex("2.25", font_size=20, color=default_color)
        c3.move_to([2.25, 2.3, 0])
        c4 = MathTex("0.25", font_size=20, color=default_color)
        c4.move_to([0.25, 2.3, 0])

        self.play(Create(r0), Create(c0))
        self.wait()

        self.play(r0.animate.move_to([1.5, 0.3, 0]), Transform(c0, c1))
        self.wait()

        self.play(r0.animate.move_to([5.5, 0.3, 0]), Transform(c0, c2))
        self.wait()
        self.play(Transform(r0, r1), Transform(c0, c3))
        self.wait()

        i = Interval(1, range=[0, 4],  stroke_width=stroke_width, color=intervals_color[1])
        i.set_opacity(0.5)
        i.add_cell_numbers(font_size=24, color=default_color)
        i.shift(2*UP)
        self.play(Create(i))
        self.wait(1)

        self.play(r0.animate.move_to([0.25, 2.3, 0]), Transform(c0, c4))
        self.wait()

class interval_scene(MovingCameraScene):
    def construct(self):
        self.next_section()
        mesh, colors, levels = init_mesh()
        [m.set_color(c) for m, c in zip(mesh, colors)]

        axe = init_axe()
        all = VGroup(axe, mesh)

        [m.shift(2*l*UP) for m, l in zip(mesh, levels)]

        level_tex = VGroup()
        for i in range(3):
            level_tex.add(MathTex(f"Level \\, {i} \\rightarrow ", font_size=36, color=intervals_color[i]))
            level_tex[-1].move_to([0, -1.5 - 0.75*i, 0])

        all_all = VGroup(all, level_tex)
        self.camera.frame.move_to(all_all)
        self.camera.frame.scale(1.1)

        [m.add_cell_numbers(font_size=24, color=default_color) for me in mesh for m in me]

        self.add(all_all)
        self.wait(1)

        r = []
        offsets = VGroup()
        index = 0
        for k, v in intervals.items():
            dx = 1./(1<<k)
            for i, e in enumerate(v):
                r.append(Rectangle(width=dx*(e[1] - e[0]), height=0.3, color=intervals_color[k], fill_opacity=0.5))
                r[-1].move_to([0.5*dx*(e[1] + e[0]), 2*k, 0])
                self.play(Create(r[-1]))
                self.wait(1)

                t = MathTex(f"[{e[0]}, {e[1]}[", font_size=36, color=intervals_color[k])
                if i==0:
                    t.next_to(level_tex[k], RIGHT)
                else:
                    t.next_to(offsets[-1], RIGHT)

                o = MathTex(f"@{index - e[0]}", font_size=36, color=offset_color)
                o.next_to(t, RIGHT)
                offsets += o

                self.play(Transform(r[-1], t))
                index += e[1] - e[0]
                self.wait()

        mesh_field = mesh.copy()
        index = 0
        for k, v in intervals.items():
            for ie, e in enumerate(v):
                for i in range(e[1]-e[0]):
                    center = mesh_field[k][ie].cell_numbers[i] .get_center()
                    mesh_field[k][ie].cell_numbers[i] = MathTex(index, font_size=24)
                    mesh_field[k][ie].cell_numbers[i].move_to(center)
                    index += 1

        self.next_section()
        self.play(Transform(mesh, mesh_field))
        self.wait()
        self.play(Create(offsets))
        self.wait()

class tree(Scene):
    def construct(self):
        self.camera.frame_width = 6
        self.camera.resize_frame_shape()

        r0 = Square(1, color=background_color, fill_color=intervals_color[0], fill_opacity=1)
        s = Square(0.1, color=background_color, fill_color=intervals_color[0], fill_opacity=1)
        r0.next_to(s, UP)

        divide = VGroup(VGroup(r0))
        s_length=1
        level = 2
        for l in range(level):
            divide += VGroup()
            for root in  divide[l]:
                for j in range(2):
                    for i in range(2):
                        r = Square(s_length/2, color=background_color, fill_color=intervals_color[0], fill_opacity=1)
                        r.move_to(root.get_center() + [-s_length/4 + s_length/2*i, -s_length/4 + s_length/2*j, 0])
                        divide[-1] += r
            s_length /= 2

        divide += VGroup()
        for root in [divide[2][3], divide[2][6], divide[2][9], divide[2][12]]:
            for j in range(2):
                for i in range(2):
                    r = Square(s_length/2, color=background_color, fill_color=intervals_color[2], fill_opacity=1)
                    r.move_to(root.get_center() + [-s_length/4 + s_length/2*i, -s_length/4 + s_length/2*j, 0])
                    divide[-1] += r

        nodes = VGroup(VGroup(s))
        lines = VGroup()
        length = 6
        for i in range(level):
            nodes += VGroup()
            lines += VGroup()
            for j in nodes[i]:
                for k in range(4):
                    n = Square(0.1, color=background_color, fill_color=intervals_color[0], fill_opacity=1)
                    n.move_to(j.get_center() + [-length/2 + length/8 + k*length/4, -0.3, 0])
                    nodes[-1] += n
                    lines[-1] += Line(j.get_center(), n.get_center(), stroke_width=1, color=default_color)
            length /= 4

        nodes += VGroup()
        lines += VGroup()

        i = 3
        for j in [nodes[2][3], nodes[2][6], nodes[2][9], nodes[2][12]]:
            for k in range(4):
                n = Square(0.1, color=background_color, fill_color=intervals_color[2], fill_opacity=1)
                n.move_to(j.get_center() + [-length/2 + length/8 + k*length/4, -0.3, 0])
                nodes[-1] += n
                lines[-1] += Line(j.get_center(), n.get_center(), stroke_width=1, color=default_color)

        self.add(r0, s)
        self.wait(2)
        for i in range(3):
            self.add(lines[i], divide[i+1])
            self.add(nodes[i], nodes[i+1])
            self.wait(2)

class mesh_2d_scene(MovingCameraScene):
    def construct(self):
        self.next_section()
        level0 = VGroup()
        for j in range(4):
            for i in range(4):
                r = Square(1, color=default_color, fill_color=intervals_color[0], fill_opacity=1)
                r.move_to([i, j, 0])
                level0.add(r)

        level1 = VGroup()
        for j in range(4):
            for i in range(4):
                r = Square(0.5, color=default_color, fill_color=intervals_color[2], fill_opacity=1)
                r.move_to([0.75+i*0.5, 0.75+j*0.5, 0])
                level1.add(r)

        self.add(level0, level1)

        self.camera.frame.move_to(level0)
        self.camera.frame.scale(0.8)

        path = [
            level0[0],
            level0[1],
            level0[4],
            level1[0],
            level1[1],
            level1[4],
            level1[5],
            level0[2],
            level0[3],
            level1[2],
            level1[3],
            level1[6],
            level1[7],
            level0[7],
            level0[8],
            level1[8],
            level1[9],
            level1[12],
            level1[13],
            level0[12],
            level0[13],
            level1[10],
            level1[11],
            level1[14],
            level1[15],
            level0[11],
            level0[14],
            level0[15],
        ]

        lines = VGroup()
        for i in range(len(path)-1):
            lines.add(Line(path[i].get_center(), path[i+1].get_center(), stroke_width=stroke_width))

        self.play(Create(lines), run_time=4)
        self.wait(1)

        stencil = VGroup()
        for i in range(5):
            r = Square(.5, color=default_color, fill_color=intervals_color[1], fill_opacity=1)
            r.move_to([1.75, 1.75, 0])
            stencil += r

        field = VGroup()
        length = 0.4
        for i in range(len(path)):
            r = Square(length, color=default_color, fill_color=path[i].color, fill_opacity=1)
            r.move_to([1.5 + length/2 - len(path)/2 * length + i*length, -1.25, 0])
            field += r

        self.play(FadeIn(field))
        self.wait(1)
        self.next_section()
        self.play(FadeOut(lines))
        self.wait(1)

        s_i = [12, 17, 22, 23, 24]
        field_s = VGroup()
        for i in s_i:
            r = Square(length, color=default_color, fill_color=intervals_color[1], fill_opacity=1)
            r.move_to([1.5 + length/2 - len(path)/2 * length + i*length, -1.25, 0])
            field_s += r

        self.play(FadeIn(field_s[2]), FadeIn(stencil))
        field_s -= field_s[2]
        self.wait(1)

        self.play(stencil[0].animate.shift([0.5, 0, 0]),
                  stencil[1].animate.shift([-0.5, 0, 0]),
                  stencil[2].animate.shift([0, 0.5, 0]),
                  stencil[3].animate.shift([0, -0.5, 0]),
                  FadeIn(field_s),
        )
        self.wait(1)


class projection(MovingCameraScene):
    def construct(self):
        i0 = Interval(0, [0, 4], color=intervals_color[0])
        i1 = Interval(1, [2, 6], color=intervals_color[1])
        i1.shift(UP)
        intervals = VGroup(i0, i1)

        eq = MathTex("u(l, i) = \\frac{u(l+1, 2i) + u(l+1, 2i+1)}{2}", tex_template=TexTemplate())
        eq.next_to(i0, 2*DOWN)

        arrows = VGroup()
        arrows.add(Arrow([1.25, 1, 0], [1.5, 0, 0], color=default_color))
        arrows.add(Arrow([1.75, 1, 0], [1.5, 0, 0], color=default_color))
        arrows.add(Arrow([2.25, 1, 0], [2.5, 0, 0], color=default_color))
        arrows.add(Arrow([2.75, 1, 0], [2.5, 0, 0], color=default_color))

        self.camera.frame.move_to(intervals)
        self.camera.frame.scale(0.7)
        self.play(Create(intervals))

        self.wait()

        self.play(FadeIn(eq))
        self.play(FadeIn(arrows))
        self.wait()

class mesh(Scene):
    def construct(self):
        self.camera.frame_width = 8.2
        self.camera.resize_frame_shape()
        self.camera.frame_center = [3, 2, 0]
        mesh, colors, levels = init_mesh()
        [m.set_color(c) for m, c in zip(mesh, colors)]
        [m.shift(2*l*UP) for m, l in zip(mesh, levels)]
        self.add(mesh)

class mesh_ghost(Scene):
    def construct(self):
        self.camera.frame_width = 8.2
        self.camera.resize_frame_shape()
        self.camera.frame_center = [3, 2, 0]

        ghost = VGroup()
        for level, v in intervals_ghost.items():
            dx = 1./(1<<level)
            for i in v:
                r = Rectangle(width=dx, height=0.1, color=intervals_color[level], fill_opacity=0.5)
                r.move_to([dx*(i[0]+0.5), 2*level, 0])
                ghost += r

        self.add(ghost)

class mesh_all(Scene):
    def construct(self):
        self.camera.frame_width = 8.2
        self.camera.resize_frame_shape()
        self.camera.frame_center = [3, 2, 0]
        mesh, colors, levels = init_mesh()
        [m.set_color(c) for m, c in zip(mesh, colors)]
        [m.shift(2*l*UP) for m, l in zip(mesh, levels)]

        ghost = VGroup()
        for level, v in intervals_ghost.items():
            dx = 1./(1<<level)
            for i in v:
                r = Rectangle(width=dx, height=0.1, color=intervals_color[level], fill_opacity=0.5)
                r.move_to([dx*(i[0]+0.5), 2*level, 0])
                ghost += r

        self.add(mesh, ghost)

class identify(MovingCameraScene):
    def construct(self):
        mesh, colors, levels = init_mesh()
        [m.set_color(c) for m, c in zip(mesh, colors)]
        [m.shift(2*l*UP) for m, l in zip(mesh, levels)]

        mesh_ghost = VGroup()
        for level, v in intervals_ghost.items():
            dx = 1./(1<<level)
            for i in v:
                r = Rectangle(width=dx, height=0.1, color=intervals_color[level], fill_opacity=0.5)
                r.move_to([dx*(i[0]+0.5), 2*level, 0])
                mesh_ghost += r

        level_tex = VGroup()
        for i in range(3):
            level_tex.add(MathTex(f"Level \\, {i} \\rightarrow ", font_size=36, color=intervals_color[i]))
            level_tex[-1].move_to([0, -1.5 - 0.75*i, 0])

        all_all = VGroup(mesh_ghost, level_tex)
        self.camera.frame.move_to(all_all)
        self.camera.frame.scale(1.1)

        t = VGroup()
        for k, v in intervals.items():
            dx = 1./(1<<k)
            for i, e in enumerate(v):
                t += MathTex(f"[{e[0]}, {e[1]}[", font_size=36, color=intervals_color[k])
                if i==0:
                    t[-1].next_to(level_tex[k], RIGHT)
                else:
                    t[-1].next_to(t[-2], RIGHT)

        t_g = VGroup()
        for k, v in intervals_ghost.items():
            dx = 1./(1<<k)
            for i, e in enumerate(v):
                t_g += MathTex(f"[{e[0]}, {e[1]}[", font_size=36, color=intervals_color[k])
                if i==0:
                    t_g[-1].next_to(level_tex[k], RIGHT)
                else:
                    t_g[-1].next_to(t_g[-2], RIGHT)

        t_a = VGroup()
        for k, v in intervals_all.items():
            dx = 1./(1<<k)
            for i, e in enumerate(v):
                t_a += MathTex(f"[{e[0]}, {e[1]}[", font_size=36, color=intervals_color[k])
                if i==0:
                    t_a[-1].next_to(level_tex[k], RIGHT)
                else:
                    t_a[-1].next_to(t_a[-2], RIGHT)

        self.next_section()
        self.add(mesh, level_tex, t)
        self.wait()

        self.next_section()

        self.play(FadeIn(mesh_ghost), FadeOut(t))
        self.wait()

        self.next_section()

        self.play(FadeOut(mesh), FadeIn(t_g))
        self.wait()

        self.next_section()

        self.play(FadeIn(mesh), Transform(t_g, t_a))
        self.wait()
