
class Registry:
    def __init__(self):
        self._registry = dict()

    def register(self, elem):
        self._registry[elem.__name__] = elem
        return elem

    def get_by_name(self, name):
        return self._registry[name]

    def is_name_valid(self, name):
        return name in self._registry.keys()

    def get_names(self):
        return self._registry.keys()
