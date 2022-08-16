import os

from django import template
from django.utils.safestring import mark_safe

register = template.Library()
cell_color = "2, 155, 194"
cell_color = "255, 210, 112"


@register.filter(name="color_code")
def coverage_col(coverage_value):

    if coverage_value is None or coverage_value == "":
        coverage_value = 0
    ncol = f"background-color: rgba({cell_color}, {int(coverage_value)}%);"
    return ncol


@register.filter(name="round")
def round_num(value):
    if value is None or value == "":
        value = 0
    ncol = round(value, 2)
    return ncol


@register.filter(name="round_to_int")
def round_to_int(value):
    if value is None or value == "":
        value = 0

    value = round(value, 5)
    return value


@register.filter(name="success_code")
def map_success_col(success_count):
    ncol = f"background-color: rgba({cell_color}, {50 * success_count}%);"
    return ncol


@register.simple_tag
def depth_color(depth_value, max_value):

    if depth_value and max_value:
        ncol = float(depth_value) * 100 / float(max_value)
    else:
        ncol = 0

    ncol = f"background-color: rgba({cell_color}, {int(ncol)}%);"
    return ncol


@register.simple_tag
def success_count_color(success):
    counts_dict = {
        "none": 0,
        "reads": 1,
        "contigs": 1,
        "reads and contigs": 2,
    }

    ncol = f"background-color: rgba({cell_color}, {counts_dict[success] * 50}%);"
    return ncol
