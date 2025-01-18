__all__ = ["MyNumberLine"]

from typing import Sequence

from manim import NumberLine
from manim.mobject.text.tex_mobject import MathTex
from manim.mobject.types.vectorized_mobject import VMobject

class MyNumberLine(NumberLine):
    def get_number_mobject(
        self,
        x: float,
        direction: Sequence[float] | None = None,
        buff: float | None = None,
        font_size: float | None = None,
        label_constructor: VMobject | None = None,
        **number_config,
    ) -> VMobject:
        if direction is None:
            direction = self.label_direction
        if buff is None:
            buff = self.line_to_number_buff
        if font_size is None:
            font_size = self.font_size
        if label_constructor is None:
            label_constructor = self.label_constructor

        num_mob = MathTex(
            int(x), font_size=font_size, **number_config
        )

        num_mob.next_to(self.number_to_point(x), direction=direction, buff=buff)
        if x < 0 and self.label_direction[0] == 0:
            # Align without the minus sign
            num_mob.shift(num_mob[0].get_width() * LEFT / 2)
        return num_mob

