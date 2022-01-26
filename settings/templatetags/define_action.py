from django import template
register = template.Library()

"""
{% load define_action %}
{% if item %}
   {% define "Edit" as action %}
{% else %}
   {% define "Create" as action %}
{% endif %}

Would you like to {{action}} this item?

OR
{% define counter|add:1 as counter %
"""
@register.simple_tag
def define(val=None):
    return val

"""
{% load your_template_tag_file %}
{% to_list 1 2 3 4 5 "yes" as my_list %}
{% for i in my_list %}
    {{ i }}
{% endfor %}
"""
@register.simple_tag
def to_list(*args):
    return args
