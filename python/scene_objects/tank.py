from typing import Union, Optional
class Tank:

    def __init__(self, dimensions: Union[float, float, float]=(None, None, None), 
                 position : Union[float, float, float]=(None, None, None), 
                 radius : float=None, type : str=None) -> None: 
        self._dimensions = dimensions
        self._position = position
        self._radius = radius
        self._type = type

    # ------------------------------ PROPERTY GETTER METHODS -------------------------------
    @property
    def dimensions(self):
        return self._dimensions

    @property
    def position(self):
        return self._position

    @property
    def radius(self):
        return self._radius
    
    @property
    def type(self):
        return self._type
    
    # ------------------------- PROPERTY SETTER METHODS --------------------------

    @dimensions.setter
    def dimensions(self, value):
        self._dimensions = value

    @position.setter
    def position(self, value):
        self._position = value

    @radius.setter
    def radius(self, value):
        self._radius = value

    @type.setter
    def type(self, value):
        self._type = value
    

     
