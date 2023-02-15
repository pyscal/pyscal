import yaml
import os

doc_fields = ["doc", "url"]

def read_yaml(filename):
    with open(filename, 'r') as fin:
        data = yaml.safe_load(fin)
        return data

def _get_doc_from_key(keydict):
    url =  keydict["url"] if "url" in keydict.keys() else None

    doc = f"""
    {keydict["doc"]}
    url: {url}
    """
    return doc

class MyList(list):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class AttrSetter:
    """
    Class which enables tab completed, contextual attributes
    """    
    def __init__(self):
        self._map_dict = {}
        self.head = None

    def __dir__(self):
        """
        Gives tab completion
        """
        return list(self._map_dict.keys())

    def __repr__(self):
        """
        String representation
        """
        return ", ".join(list(self._map_dict.keys()))

    def __getattr__(self, key):
        """
        Calls attributes from the class when general lookup fails. There are four main types of
        calls that can occur.

        Normal attribute: this function does not get called
        key which is present in `_map_dict`: filter and call it
        key which add ons: for example `position_selected` or `position_not_selected`
        key does not exist: raise `AttributeError`
        """
        if key in self._map_dict.keys():
            if self.head is None:
                self.head = self
            return self._check_filter(key)

        #now if it doesnt exists
        stripped_key, search_key, which_type, reverse = self._has_addons(key)
        #print(stripped_key, search_key, which_type, reverse)
        if (stripped_key == key):
            #this didnt work
            raise AttributeError(f"Attribute {stripped_key} not found")

        else:
            return self._get_key_on_condition(stripped_key, search_key, which_type, reverse)


    def _check_filter(self, key):
        """
        This checks for non-filtered keywords
        """
        if isinstance(self._map_dict[key], str):
            val_raw = self._map_dict[key].split("/")
            if val_raw[-1] == "nofilter":
                return self.head[val_raw[0]]
            else:
                return self.head[self._map_dict[key]][:self.head.nreal]
        else:
            return self._map_dict[key]

    def _has_addons(self, key):
        raw = key.split("_")
        which_type = None
        reverse = False
        search_key = None

        if "for" in raw:
            index = raw.index("for")
            which_type = raw[index+1] if raw[index+1] in ["real", "ghost", "all"] else None
            #now remove these keys
            del raw[index+1]
            del raw[index]

        if "with" in raw:
            index = raw.index("with")
            with_raw = raw[index+1:]

            if with_raw[0] in ["un", "non", "de", "not", "no"]:
                reverse = True
            if with_raw[-1] in ["selection", "select", "selected", "condition"]:
                search_key = "condition"
            elif with_raw[-1] in ["ghost"]:
                search_key = "ghost"
            elif with_raw[-1] in ["real"]:
                search_key = "real"
            elif with_raw[-1] in ["all"]:
                search_key = "all"
            elif with_raw[-1] in ["maskprimary"]:
                search_key = "mask_1"
            elif with_raw[-1] in ["masksecondary"]:
                search_key = "mask_2"
            
            if len(with_raw) == 2:
                del raw[index+2]
            del raw[index+1]
            del raw[index]

        stripped_key = "_".join(raw)
        return stripped_key, search_key, which_type, reverse

    def _get_key_on_condition(self, stripped_key, search_key, which_type, reverse):
        #print(stripped_key, search_key, which_type, reverse)
        search_val = None
        if stripped_key in self._map_dict.keys():
            if which_type is None:
                which_type = "all"
            if which_type == "all":
                stripped_val = self.head[self._map_dict[stripped_key]]
                if search_key is not None:
                    search_val = self.head[self._map_dict[search_key]]
            elif which_type == "real":
                stripped_val = self.head[self._map_dict[stripped_key]][:self.head.nreal]
                if search_key is not None:
                    search_val = self.head[self._map_dict[search_key]][:self.head.nreal]
            elif which_type == "ghost":
                stripped_val = self.head[self._map_dict[stripped_key]][self.head.nreal:]
                if search_key is not None:
                    search_val = self.head[self._map_dict[search_key]][self.head.nreal:]

            #we have gathered stripped_val and search_val
            if search_val is not None:
                if reverse:
                    stripped_val = [x for count, x in enumerate(stripped_val) if not search_val[count]]
                else:
                    stripped_val = [x for count, x in enumerate(stripped_val) if search_val[count]]

            return stripped_val
        
        else:
            raise AttributeError("Attribute not found")

    def _add_attribute(self, indict, head=None):
        if head is None:
            head = self
        self.head = head

        for key, val in indict.items():
            if isinstance(val, dict):
                if key not in self._map_dict.keys():
                    self._map_dict[key] = AttrSetter()

                self._map_dict[key]._add_attribute(val, head=head)
            else:
                self._map_dict[key] = val

def generate_input_template():
    indict = {"attribute_name":{
        "value": None,
        "doc": None,
        "unit": None,
        "url": None,
    }}
    return indict

def _is_sub_input(indict):
    if "value" in indict.keys():
        return False
    return True

def _create_getter(key):
    def _getter(self):
        return getattr(self, f"_{key}")
    return _getter

def _create_setter(key):
    def _setter(self, value):
        setattr(self, f"_{key}", value)
    return _setter

def _get_doc_from_key(keydict):
    url =  keydict["url"] if "url" in keydict.keys() else None
    unit =  keydict["unit"] if "unit" in keydict.keys() else None
    value =  keydict["value"] if "value" in keydict.keys() else None
    doc =  keydict["doc"] if "doc" in keydict.keys() else None
    typeinfo =  keydict["type"] if "type" in keydict.keys() else None
    optional = False if value is None else True    

    doc = f"""
    {doc}  
    Unit: {unit}
    Default: {value}
    url: {url}
    Optional: {optional}
    type: {typeinfo}
    """
    return doc

class PropertySetter(type):
    def __new__(cls, name, bases, dct, input_params):
        for key, val in input_params.items():
            if not _is_sub_input(val):
                dct[key] = property(fget=_create_getter(key), fset=_create_setter(key), doc=_get_doc_from_key(val))
        x = super().__new__(cls, name, bases, dct)
        return x
    
    def __init__(self, name, bases, dct, input_params):
        for key, val in input_params.items():
            if not _is_sub_input(val): 
                setattr(self, f"_{key}", val["value"])
        
class DocumentedKeywords:
    def _create_input_tree(self, indict, name=None, key=None):
        if name is None:
            name = "input"
            
        #create the main input
        _cmain = PropertySetter("CMain", (), {}, input_params=indict)
        
        if key is None:
            setattr(self, name, _cmain())
        else:
            head = getattr(self, name)
            setattr(head, key, _cmain())
        
        for key, val in indict.items():
            if _is_sub_input(val):
                self._create_input_tree(val, name=name, key=key)
                name = key

                

