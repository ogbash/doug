import os
import images

_RESOURCE_PATH = os.path.dirname(__file__)

def _load_resources():

    # Find all .png files in the resource directory. This construct is called
    # "list comprehension" and will be covered in detail in episode #3.
    # This returns a list of the names of all files in the resource directory
    # ending with ".png".
    resources = [ f for f in os.listdir(_RESOURCE_PATH)
                  if f.endswith(".tmpl") ]
                    
    # load resources into module namespace
    for name in resources:
          
        # this is the full path to the resource file
        path = os.path.join(_RESOURCE_PATH, name)
          
        # Now we can load the resource into the module namespace.
        # globals() gives us the dictionary of the module-global variables,
        # and we can easily extend it by assigning new keys.
        f = open(path)
        contents = f.read()
        f.close()
        globals()[name.replace('.', '_')] = contents

# load the resources when importing the module
_load_resources()
