import os
import Tkinter

_RESOURCE_PATH = os.path.dirname(__file__)

_cache={}

def getImage(name, master):
    if name in _cache:
        return _cache[name]
    
    path = os.path.join(_RESOURCE_PATH, name)
    image = Tkinter.PhotoImage(file=path)
    _cache[name] = image
    return image
