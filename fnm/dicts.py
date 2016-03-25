# -*- coding: utf-8 -*-

"""Adding dot access for dictionaries

The module exports the classes dotdict, optdict and limdict


"""

__all__ = ['dotdict', 'limdict']

import sys
import inspect
import traceback

class dotdict(dict):
    """
    Dictionary, where keys can be adressed like attributes, do not modify
    """
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        # optional processing here
        super(dotdict, self).__setitem__(key, value)

    def __getitem__(self, key):
        # optional processing here
        return super(dotdict, self).__getitem__(key)

    def __getstate__(self):
        return dict(self.items())

    def __setstate__(self,state):
        self.update(**state)

    def update(self, *args, **kwargs):
        """
        Update the dictionary
        """
        if len(args) > 1:
            raise TypeError("update expected at most 1 arguments, got %d" % len(args))
        other = dict(*args, **kwargs)
        for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def updateset(self, *args, **kwargs):
        """
        Update existing key-value pairs of the dictionary
        """
        valid_keys = self.keys()
        self.update(*args,**kwargs)

        # Remove unwanted keys
        unwanted = set(self.keys()) - set(valid_keys)
        removed = [self.pop(k, None) for k in unwanted]

        if len(removed) > 0:
            s = '';
            for j in unwanted:
                s = s + j + ', '
            s = s[:-2]
            raise Exception("Unwanted keys: '%s'" % s)

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]

    def __getattr__(self, key):
        """
        Expose all keys as attributes
        """
        return self.__getitem__(key)

    def __setattr__(self, key, value):
        """
        Expose all keys as attributes
        """
        self.__setitem__(key,value)

    def __dir__(self):
        res = dir(type(self)) + list(self.keys())
        return res

class _dotdict(dotdict):
    """
    Internal class which is identical to dotdict, but
    prints out a shorter version of functions
    """
    def __repr__(self):
        retval = "{"
        for k,v in self.iteritems():
            if not(inspect.ismethod(self[k])):
                retval = retval + "'" + k + "': " + str(self[k]) + ", "
            else:
                retval = retval + "'" + k + "': " + str(self[k].im_func) + ", "
        if (len(retval) > 2):
            retval = retval[:-2]
        retval = retval + '}'
        return retval

class limdict(dict):
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        # We assume this order of arguments for
        # tuple values: (value,min,max,default)
        if not(isinstance(value,dict)):
            if not(isinstance(value,tuple)):
                values = (value,)
            else:
                values = value
            # We are changing a value
            if self.has_key(key):
                old = self[key]
                values = (value,old.min,old.max,old.default,old.fun)
        else:
            # In case a dictionary or dotdict is used for values
            if not(value.has_key('val')):
                raise Exception('value must be set')
            values = (value['val'],)
            if value.has_key('min'):
                values = values + (value['min'],)
                if values.has_key('max'):
                    values = values + (value['max'],)
                    if values.has_key('default'):
                        values = values + (value['default'],)
                        if (values.has_key('fun')):
                            values = values + (value['fun'],)
        nvalues = len(values)

        if any((type(values[0])==float, type(values[0])==int)):
            if (nvalues < 5):
                if (type(values[0])==float):
                    if (values[0] > 0):
                        dmax = {True: 1.0, False: values[0]}[values[0] <= 1.0]
                        dmin = 0
                    else:
                        dmax = abs(values[0])
                        dmin = float(int(values[0]-0.5))
                    padding = (values[0], dmin, dmax, values[0])
                elif (type(values[0])==int):
                    padding = (values[0],-sys.maxint-1, sys.maxint,values[0])
                else:
                    raise Exception('Unsupported type')
                padding = padding + (self.validate,)
                values = values + padding[nvalues:]

            # Create dotdict and assign value
            value = _dotdict({})
            value.val       = values[0]
            value.min       = values[1]
            value.max       = values[2]
            value.default   = values[3]
            value.fun       = values[4]
            if not(value.fun(key,value)):
                raise Exception("Value outside range")
        else:
            value = values[0]

        super(limdict, self).__setitem__(key, value)

    def __getitem__(self, key):
        return super(limdict, self).__getitem__(key)

    def update(self, *args, **kwargs):
        if len(args) > 1:
            raise TypeError("update expected at most 1 arguments, got %d" % len(args))
        other = dict(*args, **kwargs)
        for key in other:
            self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]

    def validate(self, key=None, value=None):
        """
        Validate an entry during __setitem__ (key!=None, value!=None),
        Validate an existing key-value pair  (key!=None, value==None),
        or
        Validate the entire dictionary (key==None)
        """
        retval = True
        if (key==None):
            # Validate entire dictionary
            retval = all([v.fun(k) for k,v in self.iteritems()])
        else:
            if (value==None):
                # Validate existing entry
                value = self[key]
            retval = (value.min <= value.val) & (value.val <= value.max)
        return retval

    def __getattr__(self, key):
        """
        Expose all keys as attributes
        """
        return self.__getitem__(key)

    def __setattr__(self, key, value):
        """
        Expose all keys as attributes
        """
        self.__setitem__(key,value)

    def __dir__(self):
        res = dir(type(self)) + list(self.keys())
        return res

class optdict(dotdict):
    def __init__(self, *args, **kwargs):
        super(optdict, self).__init__()
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        # We assume order of tuple arguments (value, {options}, default)
        if not(isinstance(value,dict)):
            if not(isinstance(value,tuple)):
                if  isinstance(value,set):
                    # Only one value argument, no values specified
                    values = value
                    default = values.pop()
                    values.add(default)
                    values = (default, values, default)
                else:
                    values = (value, {value},value)
            else:
                # Value plus set of values
                if len(value) > 1:
                    if isinstance(value[1],set):
                        values = value
                        if len(value) > 2:
                            values = values + (value[2],)
                        else:
                            values = values + (value[0],)
                    elif isinstance(value[1],tuple):
                        options = set()
                        [options.add(each) for each in value[1]]
                        values = (value[0],) + (options,)
                        if len(value) > 2:
                            values = values + (value[2],)
                        else:
                            values = values + (value[0],)
                    else:
                        options = set()
                        [options.add(each) for each in value]
                        values = (value[0],) + (options,) + (value[0],)
                else:
                    values = (value[0], {value[0]}, value[0])
            # We are changing a value
            if self.has_key(key):
                old = self[key]
                values = (value,old.options,old.default)
        else:
            # In case a dictionary or dotdict is used for values
            if not(value.has_key('val')):
                raise Exception('value must be set')
            values = (value['val'],)

            if value.has_key('options'):
                values = values + (value['options'],)
            else:
                values = values + ({values},)

            if value.has_key('default'):
                values = values + (value['default'],)
            else:
                values = values + (value['val'],)

        value         = _optdict({})
        value.options = values[1]
        value.val     = values[0]
        value.default = values[2]
        super(optdict, self).__setitem__(key, value)

    def __repr__(self):
        retval = "{"
        for k,v in self.iteritems():
            retval = retval + "'" + k + "': " + str(self[k].__repr__custom__()) + ", "
        if (len(retval) > 2):
            retval = retval[:-2]
        retval = retval + '}'
        return retval

class _optdict(dotdict):
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
    def __setitem__(self, key, value):
        if ((key == 'val') or (key == 'default')):
            if (self.has_key('options')):
                if not(value in self['options']):
                    raise Exception("Invalid option: %s Valid options are: %s" % (str(value),self['options']))
        super(_optdict, self).__setitem__(key, value)

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        else:
            if (key == 'val'):
                if self.has_key('default'):
                    self[key] = self['default']
        return self[key]

    def __repr__(self):
        # Short form showing only value
        for k in ['val']:
            retval = str(self[k])
        return retval

    def __repr__custom__(self):
        # Long form used by optdict
        retval = "{"
        for k,v in self.iteritems():
            if not(inspect.ismethod(self[k])):
                retval = retval + "'" + k + "': " + str(self[k]) + ", "
            else:
                retval = retval + "'" + k + "': " + str(self[k].im_func) + ", "
        if (len(retval) > 2):
            retval = retval[:-2]
        retval = retval + '}'
        return retval

