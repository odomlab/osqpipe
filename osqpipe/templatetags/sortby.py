from django import template

register = template.Library()

# Source: http://djangosnippets.org/snippets/1614/

@register.filter
def sortby(sequence, attribute):
  """
  Variation on dictsort using attribute access
  Nested attributes can be used, like, "obj.attr.attr_attr"
  """
  def deep_attr(obj, attr_list):
    if len(attr_list) == 1:
      return getattr(obj, attr_list[0])
    return deep_attr(getattr(obj, attr_list[0]), attr_list[1:])
  
  lst = list(sequence)
  lst.sort(key=lambda obj: deep_attr(obj, attribute.split('.')))
  return lst
