from django import template

register = template.Library()


@register.simple_tag
def is_active(request, url):
    if request.path == url:
        return "background-color: #007bff; color: white;"
    return ""


@register.simple_tag
def column_name_check(column_name, numbers_dict):

    print("##########################3sdvedfv")

    if str(column_name) in numbers_dict.keys():

        return f"{column_name}   #{numbers_dict[str(column_name)]}"

    return column_name
