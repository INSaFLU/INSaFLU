import codecs
import os

from django import template
from django.utils.safestring import mark_safe

register = template.Library()


@register.simple_tag
def igv_app_script(directory, bamfile, indexBam):
    source_data = {}
    with codecs.open(os.path.join(directory, bamfile), "rb") as f:
        source_data = f.read()

    return mark_safe(
        f"""
    <script>
    var source_data = {source_data}
    </script>
    """
    )


@register.filter
def link_ncbi(accid):
    if accid.count("_") > 1:
        accid = accid.split("_")[:-1]
        accid = "_".join(accid)
    return "https://www.ncbi.nlm.nih.gov/nuccore/{}".format(accid)


@register.simple_tag
def dict_select(dict_name, key):
    return dict_name[key]


@register.filter
def read_html_file(html_path):

    if os.path.exists(html_path):
        html_string = codecs.open(html_path, "r").read()
        return mark_safe(html_string)
    return None


@register.filter
def strip_ext(string):

    return string

    string = string.split(".")
    if len(string) > 1:
        return string[1]
    else:
        return string[0]


@register.filter
def get_row_click(row_id):

    return f"showHideRow('{row_id}');"


@register.filter
def get_row_class(row_id):
    return f"{row_id}"


@register.simple_tag
def success_count(covplot_exists, refa_dotplot_exists):
    counts = "none"

    if covplot_exists and refa_dotplot_exists:
        counts = "reads and contigs"
    if covplot_exists and not refa_dotplot_exists:
        counts = "reads"
    if refa_dotplot_exists and not covplot_exists:
        counts = "contigs"

    return counts


@register.simple_tag
def windows_safe(windows_covered):
    if windows_covered:
        return windows_covered
    else:
        return "not calculated"
