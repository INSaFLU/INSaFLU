import os

from django import template

register = template.Library()


@register.filter
def basename(value):
    """Extract the filename from a full path."""
    return os.path.basename(value)
