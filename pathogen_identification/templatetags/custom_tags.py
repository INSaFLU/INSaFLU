from django import template

register = template.Library()


@register.simple_tag
def is_active(request, url):
    if request.path == url:
        return "background-color: #007bff; color: white;"
    return ""
