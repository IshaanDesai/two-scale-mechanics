class Registry:
    """
    A registry for storing and retrieving objects by name.

    This class provides a simple registry pattern for registering
    objects (typically classes or functions) and retrieving them by name.
    """

    def __init__(self):
        self._registry = dict()

    def register(self, elem):
        """
        Register an element in the registry using its __name__ attribute.

        Parameters
        ----------
        elem : object
            The element to register. Must have a __name__ attribute.

        Returns
        -------
        object
            The registered element (allows use as a decorator).
        """
        self._registry[elem.__name__] = elem
        return elem

    def get_by_name(self, name):
        """
        Retrieve an element from the registry by name.

        Parameters
        ----------
        name : str
            The name of the element to retrieve.

        Returns
        -------
        object
            The registered element.

        Raises
        ------
        KeyError
            If the name is not found in the registry.
        """
        return self._registry[name]

    def is_name_valid(self, name):
        """
        Check if a name exists in the registry.

        Parameters
        ----------
        name : str
            The name to check.

        Returns
        -------
        bool
            True if the name is registered, False otherwise.
        """
        return name in self._registry.keys()

    def get_names(self):
        """
        Get all registered names.

        Returns
        -------
        dict_keys
            A view of all registered names.
        """
        return self._registry.keys()
